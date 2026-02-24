# ============================================================================
# FUNCTION: Automated Cell Type Annotation with SingleR
# ============================================================================
# Purpose: Use SingleR to automatically annotate clusters by comparing
#          the full transcriptional profile against curated reference datasets.
#          More robust than marker-based methods for ambiguous populations
#          (e.g. NK vs cytotoxic CD8 T cells).
#          Use alongside annotation_scType() for cross-validation.
#
# Available references (celldex package):
#
#   Human:
#     "MonacoImmune"              - Best for PBMCs/immune: CD56bright/dim NK,
#                                   MAIT, gamma-delta T, all monocyte subtypes
#     "BlueprintEncode"           - Haematopoietic + epithelial, good for PBMCs
#     "HumanPrimaryCellAtlas"     - Broad human cell types, general purpose
#     "DatabaseImmuneCellExpression" - DICE database, immune cells only
#     "NovershternHematopoietic"  - Haematopoietic focus
#
#   Mouse:
#     "ImmGenData"                - Mouse immune cells, very detailed
#     "MouseRNAseqData"           - General mouse cell types
#
# Usage:
#   source("singleR_annotate.R")
#   result <- annotation_singleR(seurat_obj, reference = "MonacoImmune")
#   seurat_obj <- result$seurat_obj
# ============================================================================

annotation_singleR <- function(seurat_obj,
                             reference             = "MonacoImmune",
                             annotation_level      = "both",    # "cell", "cluster", or "both"
                             cluster_col           = NULL,      # NULL = use Idents(); or "seurat_clusters"
                             apply_labels          = TRUE,      # Rename Idents with SingleR cluster labels
                             low_conf_label        = "Unknown", # Label for pruned/low-confidence clusters
                             compare_with_manual   = FALSE,     # TRUE if cell_type column already exists
                             manual_annotation_col = "cell_type",
                             outdir                = "marker_results_singleR",
                             save_prefix           = "singleR") {

  # --------------------------------------------------------------------------
  # Libraries
  # --------------------------------------------------------------------------
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(BiocParallel)

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  message("\n========================================")
  message("AUTOMATED ANNOTATION WITH SingleR")
  message("========================================")
  message("Reference        : ", reference)
  message("Annotation level : ", annotation_level)
  message("Output dir       : ", outdir)

  # ==========================================================================
  # STEP 1: Load reference dataset
  # ==========================================================================
  message("\n📚 STEP 1: Loading reference dataset...")

  ref <- switch(reference,
    "MonacoImmune"                 = celldex::MonacoImmuneData(),
    "BlueprintEncode"              = celldex::BlueprintEncodeData(),
    "HumanPrimaryCellAtlas"        = celldex::HumanPrimaryCellAtlasData(),
    "DatabaseImmuneCellExpression" = celldex::DatabaseImmuneCellExpressionData(),
    "NovershternHematopoietic"     = celldex::NovershternHematopoieticData(),
    "ImmGenData"                   = celldex::ImmGenData(),
    "MouseRNAseqData"              = celldex::MouseRNAseqData(),
    stop("❌ Unknown reference: '", reference, "'\n",
         "   Choose from: MonacoImmune, BlueprintEncode, HumanPrimaryCellAtlas,\n",
         "   DatabaseImmuneCellExpression, NovershternHematopoietic,\n",
         "   ImmGenData, MouseRNAseqData")
  )

  message("  ✓ Reference loaded: ", ncol(ref), " samples, ",
          length(unique(ref$label.main)), " main cell types")

  if (!is.null(ref$label.fine)) {
    message("  ✓ Fine labels available: ", length(unique(ref$label.fine)), " subtypes")
  }

  # ==========================================================================
  # STEP 2: Extract log-normalised expression matrix
  # ==========================================================================
  message("\n🔬 STEP 2: Extracting expression data...")

  # SingleR needs log-normalised counts (the 'data' layer)
  # For SCTransform pipelines the corrected counts live in SCT,
  # but SingleR works better with RNA log-normalised data if available,
  # as SCT values are on a different scale. Try RNA first, fall back to SCT.
  assays_to_try <- c("RNA", "SCT")
  expr_matrix   <- NULL
  assay_used    <- NULL

  for (a in assays_to_try) {
    if (!a %in% names(seurat_obj@assays)) next
    candidate <- tryCatch(
      GetAssayData(seurat_obj, assay = a, layer = "data"),
      error = function(e) NULL
    )
    if (!is.null(candidate) && nrow(candidate) > 0) {
      expr_matrix <- candidate
      assay_used  <- a
      break
    }
  }

  if (is.null(expr_matrix) || nrow(expr_matrix) == 0) {
    stop("❌ No normalised data found in RNA or SCT assay.\n",
         "   Ensure NormalizeData() or SCTransform() has been run.")
  }

  message("  ✓ Assay used    : ", assay_used)
  message("  ✓ Expr matrix   : ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " cells")

  # Get cluster labels -- always character to avoid factor/integer mismatch
  if (is.null(cluster_col)) {
    cluster_labels <- as.character(Idents(seurat_obj))
  } else {
    if (!cluster_col %in% colnames(seurat_obj@meta.data))
      stop("❌ Column '", cluster_col, "' not found in meta.data.\n",
           "   Available: ", paste(colnames(seurat_obj@meta.data), collapse = ", "))
    cluster_labels <- as.character(seurat_obj@meta.data[[cluster_col]])
  }

  unique_clusters <- sort(unique(cluster_labels))
  message("  ✓ Clusters      : ", length(unique_clusters))

  # ==========================================================================
  # STEP 3: Run SingleR
  # ==========================================================================
  pred_cells    <- NULL
  pred_clusters <- NULL
  cluster_summary <- NULL

  # --- Per-cell annotation ---
  if (annotation_level %in% c("cell", "both")) {
    message("\n🔍 STEP 3a: Running per-cell annotation...")
    message("  (This may take a few minutes for large datasets...)")

    pred_cells <- SingleR(
      test     = expr_matrix,
      ref      = ref,
      labels   = ref$label.main,
      BPPARAM  = MulticoreParam(workers = 1)
    )

    n_high <- sum(!is.na(pred_cells$pruned.labels))
    n_low  <- sum(is.na(pred_cells$pruned.labels))
    message("  ✓ Per-cell annotation complete")
    message("    High confidence : ", n_high, " cells")
    message("    Low confidence  : ", n_low,  " cells (pruned -> labelled '", low_conf_label, "')")

    # Add to Seurat object -- unname() prevents "no cell overlap" error
    # (named vector uses cell type names not barcodes as names)
    pruned <- pred_cells$pruned.labels
    pruned[is.na(pruned)] <- low_conf_label

    seurat_obj$singleR.cell.label  <- unname(pred_cells$labels)
    seurat_obj$singleR.cell.pruned <- unname(pruned)
    seurat_obj$singleR.cell.score  <- unname(apply(pred_cells$scores, 1, max, na.rm = TRUE))

    # Save per-cell results
    cell_df <- data.frame(
      cell_barcode  = colnames(seurat_obj),
      singleR_label = pred_cells$labels,
      pruned_label  = pruned,
      max_score     = apply(pred_cells$scores, 1, max, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    write.csv(cell_df,
              file.path(outdir, paste0(save_prefix, "_per_cell_labels.csv")),
              row.names = FALSE)
    message("  ✓ Per-cell results saved")
  }

  # --- Per-cluster annotation ---
  if (annotation_level %in% c("cluster", "both")) {
    message("\n🔍 STEP 3b: Running per-cluster annotation...")

    # Manually aggregate expression into pseudo-bulk profiles per cluster.
    # This avoids passing clusters= to SingleR, which requires the 'scrapper'
    # package that may not be available. Aggregating first is also the
    # recommended approach for cluster-level annotation.
    message("  Aggregating expression per cluster (pseudo-bulk)...")

    agg_matrix <- sapply(unique_clusters, function(cl) {
      cells <- which(cluster_labels == cl)
      if (length(cells) == 1) {
        as.numeric(expr_matrix[, cells])
      } else {
        Matrix::rowMeans(expr_matrix[, cells, drop = FALSE])
      }
    })
    rownames(agg_matrix) <- rownames(expr_matrix)
    colnames(agg_matrix) <- unique_clusters

    message("  ✓ Aggregated matrix: ", nrow(agg_matrix), " genes x ",
            ncol(agg_matrix), " clusters")

    pred_clusters <- SingleR(
      test    = agg_matrix,
      ref     = ref,
      labels  = ref$label.main,
      BPPARAM = MulticoreParam(workers = 1)
    )

    # Build summary table -- rownames of pred_clusters are cluster IDs
    cluster_summary <- data.frame(
      cluster       = as.character(rownames(pred_clusters)),
      singleR_label = pred_clusters$labels,
      pruned_label  = pred_clusters$pruned.labels,
      max_score     = apply(pred_clusters$scores, 1, max, na.rm = TRUE),
      confidence    = ifelse(is.na(pred_clusters$pruned.labels), "Low", "High"),
      stringsAsFactors = FALSE
    )
    cluster_summary$final_label <- ifelse(
      cluster_summary$confidence == "Low",
      low_conf_label,
      cluster_summary$singleR_label
    )

    # Print results
    message("\n  📋 Per-cluster annotation results:\n")
    for (i in seq_len(nrow(cluster_summary))) {
      flag <- if (cluster_summary$confidence[i] == "High") "✓" else "⚠️  LOW CONF"
      message(sprintf("  %s  Cluster %-4s -> %-35s (score: %.3f)",
                      flag,
                      cluster_summary$cluster[i],
                      cluster_summary$final_label[i],
                      cluster_summary$max_score[i]))
    }

    n_low_cl  <- sum(cluster_summary$confidence == "Low")
    n_high_cl <- nrow(cluster_summary) - n_low_cl
    message(sprintf("\n  Summary: %d high-confidence, %d low-confidence clusters",
                    n_high_cl, n_low_cl))

    write.csv(cluster_summary,
              file.path(outdir, paste0(save_prefix, "_per_cluster_labels.csv")),
              row.names = FALSE)
    message("  ✓ Per-cluster results saved")

    # Add per-cluster label to each cell -- unname() critical here
    label_map <- setNames(cluster_summary$final_label, cluster_summary$cluster)
    seurat_obj$singleR.cluster.label <- unname(label_map[as.character(cluster_labels)])
  }

  # ==========================================================================
  # STEP 4: Apply cluster labels to Idents (optional)
  # ==========================================================================
  if (apply_labels && annotation_level %in% c("cluster", "both") &&
      !is.null(cluster_summary)) {

    message("\n🏷️  STEP 4: Applying SingleR labels to Idents...")

    current_clusters        <- as.character(levels(Idents(seurat_obj)))
    cluster_summary$cluster <- as.character(cluster_summary$cluster)

    message("  Ident levels     : ", paste(head(current_clusters, 10), collapse = ", "))
    message("  Summary clusters : ", paste(head(cluster_summary$cluster, 10), collapse = ", "))

    matched_idx <- match(current_clusters, cluster_summary$cluster)

    if (all(is.na(matched_idx))) {
      stop("❌ Cluster labels do not match between Idents and cluster_summary.\n",
           "   Ident levels     : ", paste(head(current_clusters, 5), collapse = ", "), "\n",
           "   Summary clusters : ", paste(head(cluster_summary$cluster, 5), collapse = ", "), "\n",
           "   Try: cluster_col = 'seurat_clusters'")
    }

    final_labels              <- cluster_summary$final_label[matched_idx]
    final_labels[is.na(final_labels)] <- low_conf_label
    rename_vec                <- setNames(final_labels, current_clusters)

    seurat_obj <- RenameIdents(seurat_obj, rename_vec)
    seurat_obj$cell_type <- as.character(Idents(seurat_obj))

    message("  ✓ Idents updated with SingleR labels")

    message("\n  📊 Final cell type counts:")
    cell_counts <- sort(table(seurat_obj$cell_type), decreasing = TRUE)
    for (ct in names(cell_counts)) {
      pct <- round(cell_counts[ct] / ncol(seurat_obj) * 100, 1)
      message(sprintf("    %-35s %5d cells  (%s%%)", ct, cell_counts[ct], pct))
    }
  }

  # ==========================================================================
  # STEP 5: Visualisations
  # ==========================================================================
  message("\n🎨 STEP 5: Creating visualisations...")

  # --- Per-cell UMAP ---
  if (annotation_level %in% c("cell", "both") && !is.null(pred_cells)) {

    p_cell <- DimPlot(seurat_obj,
                      group.by   = "singleR.cell.label",
                      reduction  = "umap",
                      label      = TRUE,
                      label.size = 3.5,
                      repel      = TRUE,
                      pt.size    = 0.5) +
      ggtitle(paste0("SingleR Per-Cell - ", reference)) +
      theme(plot.title  = element_text(hjust = 0.5, size = 14, face = "bold"),
            legend.text = element_text(size = 9)) +
      guides(colour = guide_legend(override.aes = list(size = 4), ncol = 1))

    ggsave(file.path(outdir, paste0(save_prefix, "_umap_per_cell.pdf")),
           p_cell, width = 12, height = 8)

    # Confidence score feature plot
    p_score <- FeaturePlot(seurat_obj,
                           features  = "singleR.cell.score",
                           reduction = "umap",
                           pt.size   = 0.5) +
      scale_colour_viridis_c(option = "plasma") +
      ggtitle("SingleR Confidence Score") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

    ggsave(file.path(outdir, paste0(save_prefix, "_confidence_scores.pdf")),
           p_score, width = 8, height = 7)

    message("  ✓ Per-cell UMAP saved")
  }

  # --- Per-cluster UMAP ---
  if (annotation_level %in% c("cluster", "both") && !is.null(pred_clusters)) {

    p_cluster <- DimPlot(seurat_obj,
                         group.by   = "singleR.cluster.label",
                         reduction  = "umap",
                         label      = TRUE,
                         label.size = 4,
                         repel      = TRUE,
                         pt.size    = 0.5) +
      ggtitle(paste0("SingleR Per-Cluster - ", reference)) +
      theme(plot.title  = element_text(hjust = 0.5, size = 14, face = "bold"),
            legend.text = element_text(size = 9)) +
      guides(colour = guide_legend(override.aes = list(size = 4), ncol = 1))

    ggsave(file.path(outdir, paste0(save_prefix, "_umap_per_cluster.pdf")),
           p_cluster, width = 12, height = 8)

    # Publication quality
    p_pub <- DimPlot(seurat_obj,
                     group.by   = "singleR.cluster.label",
                     reduction  = "umap",
                     label      = TRUE,
                     label.size = 4.5,
                     repel      = TRUE,
                     pt.size    = 0.5) +
      NoLegend() +
      xlab("UMAP 1") + ylab("UMAP 2") +
      ggtitle(paste0("SingleR Annotation - ", reference)) +
      theme(axis.title = element_text(size = 18),
            plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))

    ggsave(file.path(outdir, paste0(save_prefix, "_umap_publication.pdf")),
           p_pub, width = 10, height = 8)

    message("  ✓ Per-cluster UMAPs saved")
  }

  # --- Side-by-side: original clusters vs SingleR ---
  annotation_col <- if (annotation_level %in% c("cluster", "both")) {
    "singleR.cluster.label"
  } else {
    "singleR.cell.label"
  }

  # Use seurat_clusters for left panel to always show numeric IDs
  left_col <- if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    "seurat_clusters"
  } else NULL

  p_left <- if (!is.null(left_col)) {
    DimPlot(seurat_obj, group.by = left_col, reduction = "umap",
            label = TRUE, label.size = 3.5, repel = TRUE, pt.size = 0.5) +
      ggtitle("Original Clusters") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  } else {
    DimPlot(seurat_obj, reduction = "umap",
            label = TRUE, label.size = 3.5, repel = TRUE, pt.size = 0.5) +
      ggtitle("Original Clusters") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  }

  p_right <- DimPlot(seurat_obj,
                     group.by   = annotation_col,
                     reduction  = "umap",
                     label      = TRUE,
                     label.size = 3.5,
                     repel      = TRUE,
                     pt.size    = 0.5) +
    ggtitle(paste0("SingleR - ", reference)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p_side <- p_left | p_right

  ggsave(file.path(outdir, paste0(save_prefix, "_cluster_vs_singleR.pdf")),
         p_side, width = 20, height = 8)
  message("  ✓ Side-by-side comparison saved")

  # --- SingleR score heatmap (built-in SingleR plot) ---
  if (!is.null(pred_clusters)) {
    p_heatmap_singler <- plotScoreHeatmap(pred_clusters,
                                          show.pruned = TRUE,
                                          main        = paste0("SingleR Score Heatmap - ", reference))
    ggsave(file.path(outdir, paste0(save_prefix, "_score_heatmap.pdf")),
           plot   = p_heatmap_singler,
           width  = max(10, ncol(pred_clusters$scores) * 0.5),
           height = max(6,  nrow(pred_clusters$scores) * 0.4))
    message("  ✓ Score heatmap saved")
  }

  # --- Delta distribution (confidence check, per-cell only) ---
  # Shows the margin between the best and second-best score per cell.
  # A wide delta = confident call; narrow delta = ambiguous cell.
  if (!is.null(pred_cells)) {
    p_delta <- plotDeltaDistribution(pred_cells, ncol = 6)
    ggsave(file.path(outdir, paste0(save_prefix, "_delta_distribution.pdf")),
           plot   = p_delta,
           width  = 14,
           height = 6)
    message("  ✓ Delta distribution saved")
  }

  # --- Cell type distribution bar chart ---
  dist_col <- if (!is.null(seurat_obj$singleR.cluster.label)) {
    "singleR.cluster.label"
  } else {
    "singleR.cell.label"
  }

  cell_dist           <- as.data.frame(table(seurat_obj@meta.data[[dist_col]]),
                                        stringsAsFactors = FALSE)
  colnames(cell_dist) <- c("CellType", "Count")
  cell_dist$Percentage <- round(cell_dist$Count / ncol(seurat_obj) * 100, 1)
  cell_dist           <- cell_dist[order(-cell_dist$Count), ]

  p_bar <- ggplot(cell_dist,
                  aes(x = reorder(CellType, Count), y = Count, fill = CellType)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(Count, "  (", Percentage, "%)")),
              hjust = -0.05, size = 3.5) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(title = paste0("Cell Type Distribution - ", reference),
         x = "Cell Type", y = "Number of Cells") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title  = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.text   = element_text(size = 11),
          axis.title  = element_text(size = 13))

  ggsave(file.path(outdir, paste0(save_prefix, "_celltype_distribution.pdf")),
         p_bar,
         width  = 10,
         height = max(5, nrow(cell_dist) * 0.5))

  message("  ✓ Distribution plot saved")

  # ==========================================================================
  # STEP 6: Optional comparison with existing manual annotations
  # ==========================================================================
  if (compare_with_manual &&
      manual_annotation_col %in% colnames(seurat_obj@meta.data) &&
      !is.null(cluster_summary)) {

    message("\n🔄 STEP 6: Comparing with manual annotations...")

    comp_df <- data.frame(
      manual  = as.character(seurat_obj@meta.data[[manual_annotation_col]]),
      singleR = as.character(seurat_obj$singleR.cluster.label),
      stringsAsFactors = FALSE
    )

    conf_table <- table(Manual = comp_df$manual, SingleR = comp_df$singleR)
    write.csv(as.data.frame.matrix(conf_table),
              file.path(outdir, paste0(save_prefix, "_manual_vs_singleR.csv")))

    agreement <- comp_df %>%
      mutate(agree = manual == singleR) %>%
      group_by(manual) %>%
      summarise(agreement_pct = round(mean(agree) * 100, 1),
                n_cells        = n(),
                .groups        = "drop")

    write.csv(agreement,
              file.path(outdir, paste0(save_prefix, "_agreement_stats.csv")),
              row.names = FALSE)

    overall <- round(mean(comp_df$manual == comp_df$singleR) * 100, 1)
    message("  ✓ Overall agreement with manual: ", overall, "%")
  }

  # ==========================================================================
  # STEP 7: Write annotation summary and template
  # ==========================================================================
  message("\n💾 STEP 7: Saving summary and template...")

  if (!is.null(cluster_summary)) {
    # Summary text file
    sink(file.path(outdir, paste0(save_prefix, "_ANNOTATION_SUMMARY.txt")))
    cat("=======================================================\n")
    cat("         SingleR ANNOTATION SUMMARY\n")
    cat("=======================================================\n\n")
    cat("Reference    :", reference, "\n")
    cat("Date run     :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Total cells  :", ncol(seurat_obj), "\n")
    cat("N clusters   :", length(unique_clusters), "\n\n")
    cat("=======================================================\n")
    cat("         CLUSTER -> CELL TYPE MAPPING\n")
    cat("=======================================================\n\n")
    for (i in seq_len(nrow(cluster_summary))) {
      conf_str <- if (cluster_summary$confidence[i] == "Low") "[LOW CONFIDENCE]" else ""
      cat(sprintf("Cluster %-4s -> %-35s (score: %.3f) %s\n",
                  cluster_summary$cluster[i],
                  cluster_summary$final_label[i],
                  cluster_summary$max_score[i],
                  conf_str))
    }
    sink()
    message("  ✓ Summary saved")

    # Annotation template
    template_lines <- c(
      "# ============================================================",
      "# SingleR Annotation Template -- REVIEW AND EDIT BEFORE RUNNING",
      paste0("# Generated  : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0("# Reference  : ", reference),
      "# Clusters marked LOW CONFIDENCE should be cross-checked",
      "# against your find_markers() output before applying.",
      "# ============================================================",
      "",
      "library(Seurat)",
      "",
      "# Load your object if not already in memory",
      "# seurat_obj <- readRDS('your_seurat_object.rds')",
      "",
      "# SingleR cluster annotations -- edit any labels you disagree with",
      "new_cluster_names <- c("
    )

    for (i in seq_len(nrow(cluster_summary))) {
      conf_comment <- if (cluster_summary$confidence[i] == "Low") {
        paste0("  # LOW CONFIDENCE (score: ",
               round(cluster_summary$max_score[i], 3), ") -- check markers!")
      } else {
        paste0("  # score: ", round(cluster_summary$max_score[i], 3))
      }
      template_lines <- c(template_lines,
                          paste0('  "', cluster_summary$cluster[i], '" = "',
                                 cluster_summary$final_label[i], '",', conf_comment))
    }

    template_lines <- c(template_lines,
      ")",
      "",
      "# Apply annotations",
      "names(new_cluster_names) <- levels(seurat_obj)",
      "seurat_obj <- RenameIdents(seurat_obj, new_cluster_names)",
      "seurat_obj$cell_type <- Idents(seurat_obj)",
      "",
      "# Visualise",
      "DimPlot(seurat_obj, reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 0.5) +",
      "  NoLegend() + ggtitle('Annotated Cell Types')",
      "ggsave('final_annotated_umap.pdf', width = 10, height = 8)",
      "",
      "# Save annotated object",
      "# saveRDS(seurat_obj, 'seurat_singleR_annotated.rds')"
    )

    writeLines(template_lines,
               file.path(outdir, paste0(save_prefix, "_annotation_template.R")))
    message("  ✓ Annotation template saved")
  }

  # ==========================================================================
  # Console summary
  # ==========================================================================
  message("\n========================================")
  message("✅ SingleR ANNOTATION COMPLETE!")
  message("========================================")
  message("\n📁 Output directory : ", file.path(getwd(), outdir))
  message("\n📄 KEY FILES:")
  if (annotation_level %in% c("cluster", "both")) {
    message("   -> ", save_prefix, "_per_cluster_labels.csv    <- START HERE")
    message("   -> ", save_prefix, "_annotation_template.R     <- review & run")
    message("   -> ", save_prefix, "_score_heatmap.pdf         <- check confidence")
    message("   -> ", save_prefix, "_umap_per_cluster.pdf")
    message("   -> ", save_prefix, "_umap_publication.pdf")
  }
  if (annotation_level %in% c("cell", "both")) {
    message("   -> ", save_prefix, "_per_cell_labels.csv")
    message("   -> ", save_prefix, "_umap_per_cell.pdf")
    message("   -> ", save_prefix, "_confidence_scores.pdf")
    message("   -> ", save_prefix, "_delta_distribution.pdf    <- confidence check")
  }
  message("   -> ", save_prefix, "_cluster_vs_singleR.pdf")
  message("   -> ", save_prefix, "_celltype_distribution.pdf")
  message("   -> ", save_prefix, "_ANNOTATION_SUMMARY.txt")
  message("\n📚 NEXT STEPS:")
  message("   1. Check score heatmap and delta distribution for confidence")
  message("   2. Cross-check any low-confidence clusters against find_markers()")
  message("   3. Edit and source: ", save_prefix, "_annotation_template.R")
  message("   4. saveRDS(result$seurat_obj, 'seurat_singleR_annotated.rds')")
  message("========================================\n")

  return(list(
    seurat_obj      = seurat_obj,
    pred_cells      = pred_cells,
    pred_clusters   = pred_clusters,
    cluster_summary = cluster_summary
  ))
}