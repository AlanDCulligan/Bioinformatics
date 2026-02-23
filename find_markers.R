# ============================================================================
# FUNCTION: Find Cluster Markers and Assist with Annotation (with Presto)
# ============================================================================
# Purpose: Discover marker genes for all clusters and create annotation-ready visualizations
# Uses presto for fast differential expression analysis
# ============================================================================

find_markers <- function(seurat_obj,
                                      test_use = "wilcox", # Fast with presto
                                      min_pct = 0.25, # Only test genes expressed in at least 25% of cells in either cluster
                                      logfc_threshold = 0.25, # Only test genes with at least 0.25 log fold change
                                      only_pos = TRUE, # Only return positive markers (upregulated in cluster)
                                      top_n_genes = 10, # Number of top markers to return per cluster
                                      outdir = "marker_analysis") {
  
  library(Seurat) 
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  
  # Try to load presto for speed
  if (!requireNamespace("presto", quietly = TRUE)) {
    message("⚠️  Installing presto for faster analysis...")
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github("immunogenomics/presto")
  }
  library(presto)
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  message("\n========================================")
  message("FINDING CLUSTER MARKERS WITH PRESTO")
  message("========================================")
  
  n_clusters <- length(unique(Idents(seurat_obj)))
  message("\nAnalyzing ", n_clusters, " clusters...")
  message("Total cells: ", ncol(seurat_obj))
  
  # ========================================
  # STEP 1: Find All Markers (FAST with presto!)
  # ========================================
  message("\n📊 STEP 1: Finding markers for all clusters...")
  start_time <- Sys.time()
  
  all_markers <- FindAllMarkers(
    seurat_obj,
    only.pos = only_pos,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    test.use = test_use,
    verbose = TRUE
  )
  
  end_time <- Sys.time()
  time_taken <- round(difftime(end_time, start_time, units = "secs"), 2)
  message("✓ Found ", nrow(all_markers), " markers in ", time_taken, " seconds")
  
  # ========================================
  # STEP 2: Get Top Markers per Cluster
  # ========================================
  message("\n📋 STEP 2: Extracting top markers...")
  
  top_markers <- all_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = top_n_genes, wt = avg_log2FC) %>%
    dplyr::arrange(cluster, desc(avg_log2FC))
  
  # Get top 5 for annotation help
  top5_markers <- all_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = 5, wt = avg_log2FC) %>%
    dplyr::arrange(cluster, desc(avg_log2FC))
  
  # ========================================
  # STEP 3: Save Marker Tables
  # ========================================
  message("\n💾 STEP 3: Saving marker tables...")
  
  write.csv(all_markers, 
            file.path(outdir, "all_cluster_markers.csv"),
            row.names = FALSE)
  
  write.csv(top_markers, 
            file.path(outdir, paste0("top", top_n_genes, "_markers_per_cluster.csv")),
            row.names = FALSE)
  
  # ========================================
  # STEP 4: Annotation Guide
  # ========================================
  message("\n🏷️  STEP 4: Creating annotation guide...")
  
  sink(file.path(outdir, "ANNOTATION_GUIDE.txt"))
  cat("=======================================================\n")
  cat("         CLUSTER ANNOTATION GUIDE                      \n")
  cat("=======================================================\n\n")
  
  cat("Use these top markers to identify cell types:\n")
  cat("Look up marker genes in databases like:\n")
  cat("  - CellMarker: http://biocc.hrbmu.edu.cn/CellMarker/\n")
  cat("  - PanglaoDB: https://panglaodb.se/\n")
  cat("  - Literature for your tissue type\n\n")
  
  cat("=======================================================\n")
  cat("         TOP 5 MARKERS PER CLUSTER                     \n")
  cat("=======================================================\n\n")
  
  for (clust in sort(unique(top5_markers$cluster))) {
    cat("\n┌─ CLUSTER", clust, "─────────────────────────────────┐\n")
    clust_markers <- top5_markers %>% dplyr::filter(cluster == clust)  # ← FIXED
    cat("│ Top 5 Marker Genes:                        │\n")
    for (i in 1:min(5, nrow(clust_markers))) {
      cat(sprintf("│  %d. %-12s (Log2FC: %5.2f, p: %.2e) │\n",
                  i,
                  clust_markers$gene[i],
                  clust_markers$avg_log2FC[i],
                  clust_markers$p_val_adj[i]))
    }
    cat("│                                            │\n")
    cat("│ Suggested cell type: _____________________ │\n")
    cat("└────────────────────────────────────────────┘\n")
  }
  
  cat("\n=======================================================\n")
  cat("         ANNOTATION TEMPLATE                           \n")
  cat("=======================================================\n\n")
  
  cat("# Copy this template to annotate your clusters:\n\n")
  cat("new_cluster_names <- c(\n")
  for (clust in sort(unique(top5_markers$cluster))) {
    top_gene <- (top5_markers %>% dplyr::filter(cluster == clust))[1, "gene"]  # ← FIXED
    cat(sprintf('  "%s" = "Cell_Type_%s",  # Top marker: %s\n', 
                clust, clust, top_gene))
  }
  cat(")\n\n")
  cat("names(new_cluster_names) <- levels(seurat_obj)\n")
  cat("seurat_obj <- RenameIdents(seurat_obj, new_cluster_names)\n\n")
  
  sink()
  
  # ========================================
  # STEP 5: Visualization - Heatmap
  # ========================================
  message("\n🎨 STEP 5: Creating visualizations...")
  message("  Creating marker heatmap...")
  
  top_genes_for_heatmap <- top_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::pull(gene) %>%
    unique()
  
  p_heatmap <- DoHeatmap(seurat_obj, 
                         features = top_genes_for_heatmap,
                         size = 3,
                         angle = 90) +
    ggtitle("Top 10 Markers per Cluster") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(file.path(outdir, "01_marker_heatmap.pdf"),
         p_heatmap,
         width = 14, height = 16)
  
  # ========================================
  # STEP 6: Dot Plot
  # ========================================
  message("  Creating dot plot...")
  
  top_genes_for_dot <- top5_markers %>%
    dplyr::pull(gene) %>%
    unique()
  
  p_dotplot <- DotPlot(seurat_obj, 
                       features = top_genes_for_dot) +
    RotatedAxis() +
    ggtitle("Top 5 Markers per Cluster") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(file.path(outdir, "02_marker_dotplot.pdf"),
         p_dotplot,
         width = max(12, length(top_genes_for_dot) * 0.3),
         height = 8)
  
  # ========================================
  # STEP 7: Violin Plots (batched)
  # ========================================
  message("  Creating violin plots...")
  
  # Get top 2 per cluster for violin plots
  top2_genes <- top_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_head(n = 2) %>%
    dplyr::pull(gene) %>%
    unique()
  
  # Create violin plots in batches of 9
  n_batches <- ceiling(length(top2_genes) / 9)
  
  for (batch in 1:n_batches) {
    start_idx <- (batch - 1) * 9 + 1
    end_idx <- min(batch * 9, length(top2_genes))
    genes_batch <- top2_genes[start_idx:end_idx]
    
    p_vln <- VlnPlot(seurat_obj, 
                     features = genes_batch,
                     ncol = 3,
                     pt.size = 0.1)
    
    ggsave(file.path(outdir, paste0("03_violin_batch", batch, ".pdf")),
           p_vln,
           width = 14, height = 10)
  }
  
  # ========================================
  # STEP 8: Feature Plots (top marker per cluster)
  # ========================================
  message("  Creating feature plots...")
  
  top1_per_cluster <- top_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::pull(gene)
  
  p_feature <- FeaturePlot(seurat_obj,
                           features = top1_per_cluster,
                           ncol = 3)
  
  ggsave(file.path(outdir, "04_feature_plots_top_markers.pdf"),
         p_feature,
         width = 15,
         height = ceiling(length(top1_per_cluster) / 3) * 4)
  
  # ========================================
  # STEP 9: Known Cell Type Markers (if PBMC)
  # ========================================
  message("  Creating known marker plots...")
  
  # Common PBMC markers
  known_markers <- list(
    "T_cells" = c("CD3D", "CD3E"),
    "CD4_T" = c("IL7R", "CCR7"),
    "CD8_T" = c("CD8A", "CD8B"),
    "B_cells" = c("MS4A1", "CD79A"),
    "Monocytes" = c("CD14", "LYZ"),
    "NK_cells" = c("GNLY", "NKG7"),
    "DC" = c("FCER1A", "CST3"),
    "Platelets" = c("PPBP")
  )
  
  # Check which markers are present
  all_genes <- rownames(seurat_obj)
  available_known <- lapply(known_markers, function(x) x[x %in% all_genes])
  available_known <- available_known[sapply(available_known, length) > 0]
  
  if (length(available_known) > 0) {
    all_known_genes <- unlist(available_known)
    
    # Dot plot for known markers
    p_known_dot <- DotPlot(seurat_obj, 
                           features = all_known_genes) +
      RotatedAxis() +
      ggtitle("Known Cell Type Markers") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    ggsave(file.path(outdir, "05_known_markers_dotplot.pdf"),
           p_known_dot,
           width = max(10, length(all_known_genes) * 0.4),
           height = 8)
    
    # Feature plots for known markers
    p_known_feature <- FeaturePlot(seurat_obj,
                                   features = all_known_genes,
                                   ncol = 3)
    
    ggsave(file.path(outdir, "05_known_markers_featureplot.pdf"),
           p_known_feature,
           width = 15,
           height = ceiling(length(all_known_genes) / 3) * 4)
  }
  
  # ========================================
  # STEP 10: Summary Statistics
  # ========================================
  message("  Writing summary statistics...")
  
  sink(file.path(outdir, "SUMMARY_STATISTICS.txt"))
  cat("=======================================================\n")
  cat("         MARKER ANALYSIS SUMMARY                       \n")
  cat("=======================================================\n\n")
  
  cat("PARAMETERS:\n")
  cat("  Test used:", test_use, "\n")
  cat("  Min pct:", min_pct, "\n")
  cat("  Log FC threshold:", logfc_threshold, "\n")
  cat("  Only positive:", only_pos, "\n")
  cat("  Time taken:", time_taken, "seconds\n\n")
  
  cat("RESULTS:\n")
  cat("  Total markers found:", nrow(all_markers), "\n")
  cat("  Number of clusters:", n_clusters, "\n")
  cat("  Average markers per cluster:", 
      round(nrow(all_markers) / n_clusters, 1), "\n\n")
  
  cat("MARKERS PER CLUSTER:\n")
  marker_counts <- all_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(n_markers = n(), .groups = "drop")
  
  for (i in 1:nrow(marker_counts)) {
    cat(sprintf("  Cluster %s: %d markers\n",
                marker_counts$cluster[i],
                marker_counts$n_markers[i]))
  }
  
  cat("\n=======================================================\n")
  cat("         TOP MARKER STATISTICS                         \n")
  cat("=======================================================\n\n")
  
  top_stats <- all_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = 1) %>%
    dplyr::select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)
  
  cat("Strongest marker per cluster:\n\n")
  for (i in 1:nrow(top_stats)) {
    cat(sprintf("Cluster %s: %s (Log2FC: %.2f, Pct1: %.1f%%, Pct2: %.1f%%)\n",
                top_stats$cluster[i],
                top_stats$gene[i],
                top_stats$avg_log2FC[i],
                top_stats$pct.1[i] * 100,
                top_stats$pct.2[i] * 100))
  }
  
  cat("\n=======================================================\n")
  cat("         NEXT STEPS                                    \n")
  cat("=======================================================\n\n")
  cat("1. Review ANNOTATION_GUIDE.txt for top markers\n")
  cat("2. Look up marker genes in cell type databases\n")
  cat("3. Use the annotation template to rename clusters\n")
  cat("4. Run: source('annotate_clusters.R')\n\n")
  
  sink()
  
  # ========================================
  # STEP 11: Create Interactive Annotation Script
  # ========================================
  message("  Creating annotation helper script...")
  
  annotation_script <- paste0(
    "# ============================================================================\n",
    "# CLUSTER ANNOTATION HELPER SCRIPT\n",
    "# ============================================================================\n",
    "# Generated: ", Sys.time(), "\n",
    "# Edit the cell type names below based on marker genes\n",
    "# ============================================================================\n\n",
    "library(Seurat)\n\n",
    "# Load your object\n",
    "# seurat_obj <- readRDS('your_object.rds')\n\n",
    "# Define new cluster names based on your annotation\n",
    "new_cluster_names <- c(\n"
  )
  
  for (clust in sort(unique(top5_markers$cluster))) {
    top_genes <- (top5_markers %>% dplyr::filter(cluster == clust) %>% dplyr::slice_head(n = 3))$gene  # ← FIXED
    annotation_script <- paste0(annotation_script,
                                sprintf('  "%s" = "CellType_%s",  # Markers: %s\n',
                                        clust, clust, paste(top_genes, collapse = ", ")))
  }
  
  annotation_script <- paste0(annotation_script,
                              ")\n\n",
                              "# Apply annotations\n",
                              "names(new_cluster_names) <- levels(seurat_obj)\n",
                              "seurat_obj <- RenameIdents(seurat_obj, new_cluster_names)\n\n",
                              "# Visualize annotated clusters\n",
                              "DimPlot(seurat_obj, reduction = 'umap', label = TRUE, pt.size = 0.5) + NoLegend()\n",
                              "ggsave('annotated_umap.pdf', width = 10, height = 8)\n\n",
                              "# Save annotated object\n",
                              "# saveRDS(seurat_obj, 'seurat_annotated.rds')\n"
  )
  
  writeLines(annotation_script, file.path(outdir, "annotate_clusters.R"))
  
  # ========================================
  # Console Summary
  # ========================================
  message("\n========================================")
  message("✅ MARKER ANALYSIS COMPLETE!")
  message("========================================")
  message("\n📊 RESULTS:")
  message("   • Found ", nrow(all_markers), " total markers")
  message("   • Across ", n_clusters, " clusters")
  message("   • Analysis took ", time_taken, " seconds")
  message("\n📁 Output directory: ", file.path(getwd(), outdir))
  message("\n📄 KEY FILES CREATED:")
  message("   → ANNOTATION_GUIDE.txt (START HERE!)")
  message("   → all_cluster_markers.csv")
  message("   → top", top_n_genes, "_markers_per_cluster.csv")
  message("   → 01_marker_heatmap.pdf")
  message("   → 02_marker_dotplot.pdf")
  message("   → 05_known_markers_*.pdf (if PBMC)")
  message("   → annotate_clusters.R (ready to edit)")
  message("\n📚 NEXT STEPS:")
  message("   1. Open ANNOTATION_GUIDE.txt")
  message("   2. Look up marker genes in databases")
  message("   3. Edit annotate_clusters.R with cell type names")
  message("   4. source('", file.path(outdir, "annotate_clusters.R"), "')")
  message("========================================\n")
  
  # Return list of results
  return(list(
    all_markers = all_markers,
    top_markers = top_markers,
    time_taken = time_taken
  ))
  
  return(seurat_obj)
}
