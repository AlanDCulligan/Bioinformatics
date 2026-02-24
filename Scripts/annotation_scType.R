# ============================================================================
# FUNCTION: Automated Cell Type Annotation with scType
# ============================================================================
# Purpose: Score and annotate clusters using the scType marker database,
#          which covers a wide range of tissues beyond PBMC.
#          Fits into pipeline AFTER find_markers() as an alternative or
#          complement to annotate_clusters().
#
# Supported tissue types (pass exactly as shown):
#   "Immune system"  - T cells, B cells, NK, Monocytes, DC, Neutrophils etc.
#   "Pancreas"       - Alpha, Beta, Delta, Acinar, Ductal cells etc.
#   "Liver"          - Hepatocytes, Kupffer cells, HSC, Cholangiocytes etc.
#   "Eye"            - Retinal ganglion, Muller glia, Photoreceptors etc.
#   "Kidney"         - Podocytes, Proximal tubule, Loop of Henle etc.
#   "Brain"          - Neurons, Astrocytes, Oligodendrocytes, Microglia etc.
#   "Lung"           - AT1, AT2, Club cells, Ciliated, Endothelial etc.
#   "Adrenal"        - Chromaffin, Cortical cells etc.
#   "Heart"          - Cardiomyocytes, Fibroblasts, Endothelial etc.
#   "Intestine"      - Enterocytes, Goblet, Paneth, Enteroendocrine etc.
#   "Muscle"         - Myoblasts, Satellite cells, Smooth muscle etc.
#   "Placenta"       - Trophoblasts, Decidual cells, EVT etc.
#   "Stomach"        - Chief, Parietal, Mucous neck cells etc.
#   "Thymus"         - See Immune system (overlapping)
#
# Usage:
#   source("annotate_scType.R")
#   result <- annotation_scType(seurat_obj, tissue_type = "Immune system")
#   seurat_obj <- result$seurat_obj
# ============================================================================

annotation_scType <- function(seurat_obj,
                              tissue_type        = "Immune system",
                              assay              = "SCT",
                              cluster_col        = NULL,       # NULL = use Idents(); or e.g. "seurat_clusters"
                              apply_labels       = TRUE,       # Rename Idents with scType labels
                              low_conf_label     = "Unknown",  # Label for low-confidence clusters
                              low_conf_threshold = 0,          # Norm score at or below this = "Unknown"
                              custom_db_path     = NULL,       # Optional path to your own scType DB xlsx
                              outdir             = "marker_analysis_scType",
                              save_prefix        = "scType") {

  # --------------------------------------------------------------------------
  # Libraries
  # --------------------------------------------------------------------------
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(HGNChelper)
  library(openxlsx)

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  message("\n========================================")
  message("AUTOMATED ANNOTATION WITH scType")
  message("========================================")
  message("Tissue type : ", tissue_type)
  message("Assay       : ", assay)
  message("Output dir  : ", outdir)

  # ==========================================================================
  # STEP 1: Load scType source functions
  # ==========================================================================
  message("\n📦 STEP 1: Loading scType functions...")

  sctype_funcs_loaded <- FALSE

  tryCatch({
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
    sctype_funcs_loaded <- TRUE
    message("  ✓ Loaded scType functions from GitHub")
  }, error = function(e) {
    message("  ⚠️  Could not load from GitHub: ", conditionMessage(e))
  })

  # Fallback: bundled versions if GitHub is unavailable
  if (!sctype_funcs_loaded) {
    message("  ↩️  Using bundled fallback scType functions...")

    gene_sets_prepare <- function(path_to_db_file, cell_type) {
      cell_markers <- openxlsx::read.xlsx(path_to_db_file)
      cell_markers <- cell_markers[cell_markers$tissueType == cell_type, ]
      cell_markers$geneSymbolmore1 <- gsub(" ", "", cell_markers$geneSymbolmore1)
      cell_markers$geneSymbolmore2 <- gsub(" ", "", cell_markers$geneSymbolmore2)

      correct_genes <- function(gene_str) {
        if (is.na(gene_str) || nchar(trimws(gene_str)) == 0) return(character(0))
        genes     <- unlist(strsplit(gene_str, ","))
        genes     <- genes[nchar(trimws(genes)) > 0]
        corrected <- HGNChelper::checkGeneSymbols(genes, unmapped.as.na = FALSE)
        corrected$Suggested.Symbol[!is.na(corrected$Suggested.Symbol)]
      }

      gs_positive <- lapply(seq_len(nrow(cell_markers)), function(i)
        correct_genes(cell_markers$geneSymbolmore1[i]))
      names(gs_positive) <- cell_markers$cellName

      gs_negative <- lapply(seq_len(nrow(cell_markers)), function(i)
        correct_genes(cell_markers$geneSymbolmore2[i]))
      names(gs_negative) <- cell_markers$cellName

      list(gs_positive = gs_positive, gs_negative = gs_negative)
    }

    sctype_score <- function(scRNAseqData, scaled = TRUE, gs, gs2 = NULL, ...) {
      if (!scaled) scRNAseqData <- scale(scRNAseqData)

      score_matrix <- sapply(names(gs), function(ct) {
        pos <- intersect(gs[[ct]],  rownames(scRNAseqData))
        neg <- if (!is.null(gs2)) intersect(gs2[[ct]], rownames(scRNAseqData)) else character(0)

        pos_score <- if (length(pos) > 0) colSums(scRNAseqData[pos, , drop = FALSE]) else
          rep(0, ncol(scRNAseqData))
        neg_score <- if (length(neg) > 0) colSums(scRNAseqData[neg, , drop = FALSE]) else
          rep(0, ncol(scRNAseqData))

        pos_score - neg_score
      })

      t(score_matrix)
    }

    message("  ✓ Fallback functions ready")
  }

  # ==========================================================================
  # STEP 2: Load scType marker database
  # ==========================================================================
  message("\n📚 STEP 2: Loading scType marker database...")

  if (!is.null(custom_db_path)) {
    if (!file.exists(custom_db_path))
      stop("❌ Custom DB not found: ", custom_db_path)
    db_path <- custom_db_path
    message("  ✓ Using custom database: ", custom_db_path)

  } else {
    db_url   <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    local_db <- file.path(outdir, "ScTypeDB_full.xlsx")

    if (!file.exists(local_db)) {
      message("  Downloading scType database...")
      tryCatch({
        download.file(db_url, destfile = local_db, mode = "wb", quiet = TRUE)
        message("  ✓ Database downloaded to: ", local_db)
      }, error = function(e) {
        stop("❌ Could not download scType database.\n",
             "   Check internet or supply custom_db_path.\n",
             "   Manual download: https://github.com/IanevskiAleksandr/sc-type/raw/master/ScTypeDB_full.xlsx")
      })
    } else {
      message("  ✓ Using cached database: ", local_db)
    }
    db_path <- local_db
  }

  # Validate tissue type
  db_preview        <- openxlsx::read.xlsx(db_path)
  available_tissues <- unique(db_preview$tissueType)

  if (!tissue_type %in% available_tissues) {
    stop("❌ Tissue type '", tissue_type, "' not found.\n",
         "   Available options:\n",
         paste0("     - ", available_tissues, collapse = "\n"))
  }

  gs_list      <- gene_sets_prepare(db_path, tissue_type)
  n_cell_types <- length(gs_list$gs_positive)
  n_pos_genes  <- sum(sapply(gs_list$gs_positive, length))
  n_neg_genes  <- sum(sapply(gs_list$gs_negative, length))

  message("  ✓ Tissue           : ", tissue_type)
  message("  ✓ Cell types in DB : ", n_cell_types)
  message("  ✓ Positive markers : ", n_pos_genes)
  message("  ✓ Negative markers : ", n_neg_genes)

  # Save DB used (useful for methods section)
  db_used <- db_preview[db_preview$tissueType == tissue_type, ]
  write.csv(db_used,
            file.path(outdir, paste0(save_prefix, "_marker_database_used.csv")),
            row.names = FALSE)

  # ==========================================================================
  # STEP 3: Extract scaled expression matrix
  # ==========================================================================
  message("\n🔬 STEP 3: Extracting scaled expression data...")

  # SCTransform stores scaled data in the SCT assay, not RNA -- detect automatically
  active_assay  <- DefaultAssay(seurat_obj)
  assays_to_try <- unique(c(assay, active_assay, "SCT", "RNA"))

  scaled_data <- NULL
  assay_used  <- NULL

  for (a in assays_to_try) {
    if (!a %in% names(seurat_obj@assays)) next
    candidate <- tryCatch(
      GetAssayData(seurat_obj, assay = a, layer = "scale.data"),
      error = function(e) NULL
    )
    if (!is.null(candidate) && nrow(candidate) > 0) {
      scaled_data <- candidate
      assay_used  <- a
      break
    }
  }

  if (is.null(scaled_data) || nrow(scaled_data) == 0) {
    stop("❌ No scaled data found (tried: ", paste(assays_to_try, collapse = ", "), ").\n",
         "   SCTransform users: ensure SCTransform() has been run.\n",
         "   Standard pipeline: run ScaleData(seurat_obj) first.")
  }

  if (assay_used != assay)
    message("  ℹ️  Scaled data found in '", assay_used,
            "' (not '", assay, "') -- expected for SCTransform pipelines")

  message("  ✓ Assay used    : ", assay_used)
  message("  ✓ Scaled matrix : ", nrow(scaled_data), " genes x ", ncol(scaled_data), " cells")

  # Get cluster labels -- always as character to avoid factor/integer mismatch later
  if (is.null(cluster_col)) {
    clusters <- as.character(Idents(seurat_obj))
  } else {
    if (!cluster_col %in% colnames(seurat_obj@meta.data))
      stop("❌ Column '", cluster_col, "' not found in meta.data.\n",
           "   Available: ", paste(colnames(seurat_obj@meta.data), collapse = ", "))
    clusters <- as.character(seurat_obj@meta.data[[cluster_col]])
  }

  # Sort numerically if all cluster names are numbers, else alphabetically
  unique_clusters <- sort(unique(clusters))

  message("  ✓ Clusters to annotate: ", length(unique_clusters))

  # ==========================================================================
  # STEP 4: Run scType scoring
  # ==========================================================================
  message("\n🧮 STEP 4: Running scType scoring...")
  message("  (", ncol(scaled_data), " cells x ", n_cell_types, " cell types)")

  es_max <- sctype_score(
    scRNAseqData = scaled_data,
    scaled       = TRUE,
    gs           = gs_list$gs_positive,
    gs2          = gs_list$gs_negative
  )

  # sctype_score returns cell_types x cells; verify and correct if transposed
  if (nrow(es_max) != n_cell_types && ncol(es_max) == n_cell_types) {
    message("  ℹ️  Transposing score matrix to cell_types x cells")
    es_max <- t(es_max)
  }

  message("  ✓ Score matrix: ", nrow(es_max), " cell types x ", ncol(es_max), " cells")

  # ==========================================================================
  # STEP 5: Aggregate scores per cluster and assign labels
  # ==========================================================================
  message("\n📋 STEP 5: Assigning cell types per cluster...")

  cluster_results <- lapply(unique_clusters, function(cl) {
    cell_idx <- which(clusters == cl)
    n_cells  <- length(cell_idx)

    # Sum scores across all cells in cluster, normalise by cluster size
    cluster_scores <- sort(rowSums(es_max[, cell_idx, drop = FALSE]),
                           decreasing = TRUE)

    top_score      <- cluster_scores[1]
    norm_score     <- top_score / n_cells
    best_cell_type <- names(cluster_scores)[1]
    second_best    <- if (length(cluster_scores) > 1) names(cluster_scores)[2] else NA_character_

    score_gap <- if (!is.na(second_best)) {
      (cluster_scores[1] - cluster_scores[2]) / n_cells
    } else Inf

    is_low_conf <- norm_score <= low_conf_threshold
    final_label <- if (is_low_conf) low_conf_label else best_cell_type

    data.frame(
      cluster      = cl,
      scType_label = best_cell_type,
      final_label  = final_label,
      raw_score    = round(top_score, 2),
      norm_score   = round(norm_score, 4),
      score_gap    = round(score_gap, 4),
      second_best  = second_best,
      n_cells      = n_cells,
      low_conf     = is_low_conf,
      stringsAsFactors = FALSE
    )
  })

  cluster_summary           <- do.call(rbind, cluster_results)
  rownames(cluster_summary) <- NULL
  cluster_summary$cluster   <- as.character(cluster_summary$cluster)

  # Print results
  message("\n  📋 scType annotation results:\n")
  for (i in seq_len(nrow(cluster_summary))) {
    flag <- if (cluster_summary$low_conf[i]) "⚠️  LOW CONF" else "✓"
    message(sprintf("  %s  Cluster %-4s -> %-35s (norm score: %6.4f, n=%d)",
                    flag,
                    cluster_summary$cluster[i],
                    cluster_summary$final_label[i],
                    cluster_summary$norm_score[i],
                    cluster_summary$n_cells[i]))
  }

  n_low  <- sum(cluster_summary$low_conf)
  n_high <- nrow(cluster_summary) - n_low
  message(sprintf("\n  Summary: %d high-confidence, %d low-confidence clusters", n_high, n_low))

  write.csv(cluster_summary,
            file.path(outdir, paste0(save_prefix, "_cluster_annotations.csv")),
            row.names = FALSE)
  message("  ✓ Cluster annotations saved")

  # --------------------------------------------------------------------------
  # Add per-cell metadata to Seurat object
  # IMPORTANT: unname() strips cluster-ID names from the vector.
  # Seurat adds metadata by position (one value per cell in order),
  # so a named vector causes "no cell overlap" because names are cluster IDs
  # not cell barcodes.
  # --------------------------------------------------------------------------
  clusters <- as.character(clusters)  # ensure character for map lookup

  label_map  <- setNames(cluster_summary$final_label,  cluster_summary$cluster)
  sctype_map <- setNames(cluster_summary$scType_label, cluster_summary$cluster)
  score_map  <- setNames(cluster_summary$norm_score,   cluster_summary$cluster)

  seurat_obj$scType.label      <- unname(label_map[clusters])
  seurat_obj$scType.raw.label  <- unname(sctype_map[clusters])
  seurat_obj$scType.norm.score <- unname(score_map[clusters])

  # ==========================================================================
  # STEP 6: Apply labels to Idents (optional)
  # ==========================================================================
  if (apply_labels) {
    message("\n🏷️  STEP 6: Applying scType labels to Idents...")

    # Force to character on both sides to avoid factor/integer mismatch
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

    message("  ✓ Idents updated with scType labels")

    message("\n  📊 Final cell type counts:")
    cell_counts <- sort(table(seurat_obj$cell_type), decreasing = TRUE)
    for (ct in names(cell_counts)) {
      pct <- round(cell_counts[ct] / ncol(seurat_obj) * 100, 1)
      message(sprintf("    %-35s %5d cells  (%s%%)", ct, cell_counts[ct], pct))
    }
  }

  # ==========================================================================
  # STEP 7: Visualisations
  # ==========================================================================
  message("\n🎨 STEP 7: Creating visualisations...")

  # --- UMAP: scType labels (with legend) ---
  p_sctype <- DimPlot(seurat_obj,
                      group.by   = "scType.label",
                      reduction  = "umap",
                      label      = TRUE,
                      label.size = 4,
                      repel      = TRUE,
                      pt.size    = 0.5) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    ggtitle(paste0("scType Annotation - ", tissue_type)) +
    theme(axis.title  = element_text(size = 14),
          legend.text = element_text(size = 10),
          plot.title  = element_text(hjust = 0.5, size = 16, face = "bold")) +
    guides(colour = guide_legend(override.aes = list(size = 5), ncol = 1))

  ggsave(file.path(outdir, paste0(save_prefix, "_umap_annotated.pdf")),
         p_sctype, width = 12, height = 8)

  # --- Publication quality (no legend) ---
  p_pub <- DimPlot(seurat_obj,
                   group.by   = "scType.label",
                   reduction  = "umap",
                   label      = TRUE,
                   label.size = 4.5,
                   repel      = TRUE,
                   pt.size    = 0.5) +
    NoLegend() +
    xlab("UMAP 1") + ylab("UMAP 2") +
    ggtitle(paste0("scType Annotation - ", tissue_type)) +
    theme(axis.title = element_text(size = 18),
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))

  ggsave(file.path(outdir, paste0(save_prefix, "_umap_publication.pdf")),
         p_pub, width = 10, height = 8)

  # --- Side-by-side: original cluster numbers vs scType labels ---
  # Use seurat_clusters column directly to always show numeric cluster IDs
  # regardless of whether Idents have been renamed
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    p_clusters <- DimPlot(seurat_obj,
                          group.by   = "seurat_clusters",
                          reduction  = "umap",
                          label      = TRUE,
                          label.size = 3.5,
                          repel      = TRUE,
                          pt.size    = 0.5) +
      ggtitle("Original Clusters") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  } else {
    # Fallback if seurat_clusters column not present
    p_clusters <- DimPlot(seurat_obj,
                          reduction  = "umap",
                          label      = TRUE,
                          label.size = 3.5,
                          repel      = TRUE,
                          pt.size    = 0.5) +
      ggtitle("Original Clusters") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  }

  p_side <- p_clusters | p_sctype

  ggsave(file.path(outdir, paste0(save_prefix, "_cluster_vs_scType.pdf")),
         p_side, width = 20, height = 8)

  message("  ✓ UMAP plots saved")

  # --- Score heatmap (cell types x clusters) ---
  message("  Creating score heatmap...")

  score_matrix <- sapply(unique_clusters, function(cl) {
    cell_idx <- which(clusters == cl)
    rowSums(es_max[, cell_idx, drop = FALSE]) / length(cell_idx)
  })
  colnames(score_matrix) <- unique_clusters

  # Keep top 25 cell types by max score for readability
  top_n_types  <- min(25, nrow(score_matrix))
  top_type_idx <- order(apply(score_matrix, 1, max), decreasing = TRUE)[1:top_n_types]
  score_top    <- score_matrix[top_type_idx, , drop = FALSE]

  score_df        <- as.data.frame(score_top)
  score_df$cell_type <- rownames(score_df)
  score_long      <- tidyr::pivot_longer(score_df,
                                         cols      = -cell_type,
                                         names_to  = "cluster",
                                         values_to = "norm_score")

  # Fix factor ordering so clusters appear in numeric order on x axis
  score_long$cluster   <- factor(score_long$cluster,   levels = unique_clusters)
  score_long$cell_type <- factor(score_long$cell_type,
                                 levels = rownames(score_top)[
                                   order(apply(score_top, 1, which.max))])

  p_heatmap <- ggplot(score_long,
                      aes(x = cluster, y = cell_type, fill = norm_score)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(low      = "white",
                         mid      = "#fee8c8",
                         high     = "#e34a33",
                         midpoint = max(score_long$norm_score) / 2,
                         name     = "Norm.\nScore") +
    geom_text(aes(label = ifelse(norm_score > max(norm_score) * 0.25,
                                 round(norm_score, 2), "")),
              size = 2.5, colour = "black") +
    labs(title = paste0("scType Score Heatmap - ", tissue_type),
         x = "Cluster", y = "Cell Type") +
    theme_minimal() +
    theme(plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x  = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y  = element_text(size = 9),
          axis.title   = element_text(size = 12),
          legend.title = element_text(size = 10))

  ggsave(file.path(outdir, paste0(save_prefix, "_score_heatmap.pdf")),
         p_heatmap,
         width  = max(8, length(unique_clusters) * 0.8),
         height = max(6, top_n_types * 0.35))

  message("  ✓ Score heatmap saved")

  # --- Cell type distribution bar chart ---
  message("  Creating cell type distribution plot...")

  cell_dist           <- as.data.frame(table(seurat_obj$scType.label),
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
    labs(title = "Cell Type Distribution",
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
  # STEP 8: Save annotation summary
  # ==========================================================================
  message("\n💾 STEP 8: Saving annotation summary...")

  summary_path <- file.path(outdir, paste0(save_prefix, "_ANNOTATION_SUMMARY.txt"))
  sink(summary_path)
  cat("=======================================================\n")
  cat("         scType ANNOTATION SUMMARY\n")
  cat("=======================================================\n\n")
  cat("Tissue type  :", tissue_type, "\n")
  cat("Date run     :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Total cells  :", ncol(seurat_obj), "\n")
  cat("N clusters   :", length(unique_clusters), "\n")
  cat("N cell types :", length(unique(cluster_summary$final_label)), "\n\n")
  cat("=======================================================\n")
  cat("         CLUSTER -> CELL TYPE MAPPING\n")
  cat("=======================================================\n\n")
  for (i in seq_len(nrow(cluster_summary))) {
    conf_str <- if (cluster_summary$low_conf[i]) "[LOW CONFIDENCE]" else ""
    cat(sprintf("Cluster %-4s -> %-35s (norm score: %.4f, n=%d) %s\n",
                cluster_summary$cluster[i],
                cluster_summary$final_label[i],
                cluster_summary$norm_score[i],
                cluster_summary$n_cells[i],
                conf_str))
  }
  cat("\n=======================================================\n")
  cat("         CELL TYPE DISTRIBUTION\n")
  cat("=======================================================\n\n")
  print(cell_dist, row.names = FALSE)
  cat("\n=======================================================\n")
  cat("         NEXT STEPS\n")
  cat("=======================================================\n\n")
  cat("1. Review score heatmap to check confidence of assignments\n")
  cat("2. Inspect low-confidence clusters against find_markers() output\n")
  cat("3. Edit the annotation template if any labels need changing\n")
  sink()

  message("  ✓ Summary saved")

  # ==========================================================================
  # STEP 9: Write editable annotation template
  # ==========================================================================
  message("\n📄 STEP 9: Writing annotation template...")

  template_lines <- c(
    "# ============================================================",
    "# scType Annotation Template -- REVIEW AND EDIT BEFORE RUNNING",
    paste0("# Generated : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("# Tissue    : ", tissue_type),
    "# Clusters marked LOW CONF should be checked against",
    "# your marker heatmap from find_markers()",
    "# ============================================================",
    "",
    "library(Seurat)",
    "",
    "# Load your object if not already in memory",
    "# seurat_obj <- readRDS('your_seurat_object.rds')",
    "",
    "# Cluster annotations from scType -- edit any labels you disagree with",
    "new_cluster_names <- c("
  )

  for (i in seq_len(nrow(cluster_summary))) {
    conf_comment <- if (cluster_summary$low_conf[i]) {
      paste0("  # LOW CONF (score: ", round(cluster_summary$norm_score[i], 4), ") -- check markers!")
    } else {
      paste0("  # score: ", round(cluster_summary$norm_score[i], 4))
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
    "# saveRDS(seurat_obj, 'seurat_scType_annotated.rds')"
  )

  writeLines(template_lines,
             file.path(outdir, paste0(save_prefix, "_annotation_template.R")))

  message("  ✓ Template saved")

  # ==========================================================================
  # Console summary
  # ==========================================================================
  message("\n========================================")
  message("✅ scType ANNOTATION COMPLETE!")
  message("========================================")
  message("\n📁 Output directory : ", file.path(getwd(), outdir))
  message("\n📄 KEY FILES:")
  message("   -> ", save_prefix, "_cluster_annotations.csv      <- START HERE")
  message("   -> ", save_prefix, "_annotation_template.R        <- review & run")
  message("   -> ", save_prefix, "_score_heatmap.pdf            <- check confidence")
  message("   -> ", save_prefix, "_umap_annotated.pdf")
  message("   -> ", save_prefix, "_umap_publication.pdf")
  message("   -> ", save_prefix, "_cluster_vs_scType.pdf")
  message("   -> ", save_prefix, "_celltype_distribution.pdf")
  message("   -> ", save_prefix, "_ANNOTATION_SUMMARY.txt")
  message("   -> ", save_prefix, "_marker_database_used.csv     <- for methods section")
  message("\n📚 NEXT STEPS:")
  message("   1. Check score heatmap for ambiguous clusters")
  message("   2. Review low-confidence clusters against find_markers() output")
  message("   3. Edit and source: ", save_prefix, "_annotation_template.R")
  message("   4. saveRDS(result$seurat_obj, 'seurat_scType_annotated.rds')")
  message("========================================\n")

  return(list(
    seurat_obj      = seurat_obj,
    cluster_summary = cluster_summary,
    score_matrix    = score_matrix,
    gs_list         = gs_list
  ))
}