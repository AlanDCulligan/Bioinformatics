# ============================================================================
# FUNCTION: Clustering → UMAP → tSNE → Visualization
# ============================================================================
# Purpose: Cluster cells and visualize with chosen number of PCs
# RUN THIS AFTER reviewing elbow plot from normalise_and_pcs()
# ============================================================================
# DoubletCheck Usage:
#   1st run:  run_clustering(seurat_obj)                          # DoubletCheck off by default
#   2nd run:  run_clustering(seurat_obj, DoubletCheck = TRUE,     # Run on specific clusters
#                            doublet_clusters = c(8))
# ============================================================================

run_clustering <- function(seurat_obj,
                           n_dims = 10,          # Choose based on elbow plot and amount of dims you want to use for clustering and dimensionality reduction (UMAP/tSNE)
                           run_tsne = TRUE,
                           resolution = 0.5,
                           mt_col = "percent.mt",
                           outdir = "qc_post",
                           DoubletCheck = FALSE,    # Set to TRUE to run doublet detection on specified clusters
                           doublet_clusters = NULL, # Vector of cluster IDs suspected to be doublets (e.g., c(3, 7))
                           doublet_pc_n = 10) {
  
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  
  # Ensure output directory exists
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  message("\n========================================")
  message("RUNNING DOWNSTREAM ANALYSIS")
  message(paste0("Using ", n_dims, " PCs"))
  message("========================================")
  
  # ========================================
  # Check if PCA exists
  # ========================================
  if (!"pca" %in% names(seurat_obj@reductions)) {
    stop("❌ PCA not found! Please run normalise_and_pcs() first.")
  }
  
  # ========================================
  # STEP 1: Cluster the cells
  # ========================================
  message("\nSTEP 1: Clustering cells...")
  
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_dims, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, verbose = FALSE)
  
  n_clusters <- length(unique(seurat_obj$seurat_clusters))
  message("  ✓ Found ", n_clusters, " clusters")
  message("  Cluster IDs of first 5 cells:")
  print(head(Idents(seurat_obj), 5))
  
  # ========================================
  # STEP 2: Run UMAP
  # ========================================
  message("\nSTEP 2: Running UMAP...")
  
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_dims, verbose = FALSE)
  message("  ✓ UMAP complete")
  
  # ========================================
  # STEP 3: Run tSNE (optional)
  # ========================================
  if (run_tsne) {
    message("\nSTEP 3: Running tSNE...")
    seurat_obj <- RunTSNE(seurat_obj, dims = 1:n_dims, verbose = FALSE)
    message("  ✓ tSNE complete")
  }
  
  # ========================================
  # STEP 4: Generate Plots and Save PDFs
  # ========================================
  message("\nSTEP 4: Generating visualization plots...")
  
  # PCA plot
  message("  Creating PCA plot...")
  p_pca <- DimPlot(seurat_obj, 
                   reduction = "pca", 
                   dims = c(1, 2),
                   group.by = "seurat_clusters", 
                   label = TRUE, 
                   pt.size = 0.5) +
    ggtitle(paste0("PCA: Clusters (", n_dims, " PCs)")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(file.path(outdir, "07_pca_with_clusters.pdf"), 
         p_pca, width = 10, height = 8)
  
  # UMAP plots
  message("  Creating UMAP plots...")
  p_umap <- DimPlot(seurat_obj, 
                    reduction = "umap", 
                    label = TRUE, 
                    pt.size = 0.5) +
    ggtitle(paste0("UMAP: ", n_clusters, " Clusters (", n_dims, " PCs)")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(file.path(outdir, "08_umap_clusters.pdf"), 
         p_umap, width = 10, height = 8)
  
  p_umap_nolabel <- DimPlot(seurat_obj, 
                            reduction = "umap", 
                            label = FALSE, 
                            pt.size = 0.5) +
    ggtitle("UMAP: Cell Clusters")
  
  ggsave(file.path(outdir, "08_umap_clusters_no_labels.pdf"), 
         p_umap_nolabel, width = 10, height = 8)
  
  # tSNE plot
  if (run_tsne) {
    message("  Creating tSNE plots...")
    p_tsne <- DimPlot(seurat_obj, 
                      reduction = "tsne", 
                      label = TRUE, 
                      pt.size = 0.5) +
      ggtitle(paste0("tSNE: ", n_clusters, " Clusters (", n_dims, " PCs)")) +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    ggsave(file.path(outdir, "09_tsne_clusters.pdf"), 
           p_tsne, width = 10, height = 8)
    
    # Side-by-side comparison
    p_comparison <- p_umap + p_tsne + 
      plot_annotation(title = "UMAP vs tSNE Comparison")
    
    ggsave(file.path(outdir, "10_umap_vs_tsne_comparison.pdf"), 
           p_comparison, width = 18, height = 8)
  }
  
  # QC overlays
  message("  Creating QC metric overlays...")
  p_qc1 <- FeaturePlot(seurat_obj, 
                       features = "nFeature_RNA", 
                       reduction = "umap", 
                       pt.size = 0.5) +
    ggtitle("Genes per Cell") + 
    scale_color_viridis_c(option = "plasma") +
    theme(plot.title = element_text(hjust = 0.5))
  
  p_qc2 <- FeaturePlot(seurat_obj, 
                       features = "nCount_RNA", 
                       reduction = "umap", 
                       pt.size = 0.5) +
    ggtitle("UMI Count per Cell") + 
    scale_color_viridis_c(option = "plasma") +
    theme(plot.title = element_text(hjust = 0.5))
  
  p_qc3 <- FeaturePlot(seurat_obj, 
                       features = mt_col, 
                       reduction = "umap", 
                       pt.size = 0.5) +
    ggtitle("Mitochondrial %") + 
    scale_color_viridis_c(option = "plasma") +
    theme(plot.title = element_text(hjust = 0.5))
  
  p_qc_comb <- p_qc1 + p_qc2 + p_qc3 +
    plot_annotation(
      title = "QC Metrics on UMAP - Check for Technical Artifacts",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
  
  ggsave(file.path(outdir, "11_umap_qc_metrics_overlay.pdf"), 
         p_qc_comb, width = 18, height = 6)
  
  # Violin plots
  message("  Creating violin plots...")
  p_vln <- VlnPlot(seurat_obj, 
                   features = c("nFeature_RNA", "nCount_RNA", mt_col), 
                   ncol = 3, 
                   pt.size = 0.01) +
    plot_annotation(
      title = "QC Metrics Distribution by Cluster",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
  
  ggsave(file.path(outdir, "12_violin_by_cluster.pdf"), 
         p_vln, width = 16, height = 6)
  
  # Top variable genes heatmap
  message("  Creating variable genes heatmap...")
  top_genes <- head(VariableFeatures(seurat_obj), 50)
  
  p_heat <- DoHeatmap(seurat_obj, 
                      features = top_genes, 
                      size = 3, 
                      angle = 90) +
    ggtitle("Top 50 Variable Genes by Cluster") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.y = element_text(size = 6))
  
  ggsave(file.path(outdir, "13_top_variable_genes_heatmap.pdf"), 
         p_heat, width = 12, height = 14)
  
  # Cells per cluster barplot
  message("  Creating cells per cluster plot...")
  cluster_counts <- table(seurat_obj$seurat_clusters)
  cluster_df <- data.frame(
    Cluster = names(cluster_counts),
    Count = as.numeric(cluster_counts),
    Percentage = round(as.numeric(cluster_counts) / ncol(seurat_obj) * 100, 1)
  )
  
  p_bar <- ggplot(cluster_df, aes(x = Cluster, y = Count, fill = Cluster)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), 
              vjust = -0.3, size = 4) +
    ggtitle("Number of Cells per Cluster") + 
    labs(y = "Number of Cells", x = "Cluster") +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  
  ggsave(file.path(outdir, "14_cells_per_cluster.pdf"), 
         p_bar, width = 10, height = 6)
  
  # ========================================
  # STEP 5: scDblFinder Doublet Detection (optional)
  # ========================================
  doublet_results <- NULL
  
  if (DoubletCheck) {
    
    message("\nSTEP 5: Running scDblFinder doublet detection...")
    
    # Validate inputs
    if (is.null(doublet_clusters)) {
      stop("❌ Please specify which clusters to check: doublet_clusters = c(8, 9)")
    }
    
    # Check scDblFinder is installed
    if (!requireNamespace("scDblFinder", quietly = TRUE)) {
      stop("❌ scDblFinder not installed. Run: BiocManager::install('scDblFinder')")
    }
    
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("❌ SingleCellExperiment not installed. Run: BiocManager::install('SingleCellExperiment')")
    }
    
    library(scDblFinder)
    library(SingleCellExperiment)
    
    message("  Checking clusters: ", paste(doublet_clusters, collapse = ", "))
    
    # Subset to clusters of interest
    cells_to_check <- WhichCells(seurat_obj, 
                                 idents = as.character(doublet_clusters))
    seurat_sub <- subset(seurat_obj, cells = cells_to_check)
    
    message("  Subsetting to ", ncol(seurat_sub), " cells from specified clusters")
    
    # Convert to SingleCellExperiment (scDblFinder requires this)
    # Use RNA assay raw counts
    sce <- as.SingleCellExperiment(seurat_sub, assay = "RNA")
    
    # Run scDblFinder
    message("  Running scDblFinder...")
    set.seed(42) # for reproducibility
    sce <- scDblFinder(sce)
    
    message("  ✓ scDblFinder complete")
    
    # ----------------------------------------
    # Extract results
    # ----------------------------------------
    doublet_results <- data.frame(
      Cell_Barcode  = colnames(sce),
      Cluster       = as.character(seurat_sub$seurat_clusters),
      Doublet_Score = sce$scDblFinder.score,
      Classification = sce$scDblFinder.class  # "singlet" or "doublet"
    )
    
    n_doublets  <- sum(doublet_results$Classification == "doublet")
    n_singlets  <- sum(doublet_results$Classification == "singlet")
    doublet_pct <- round(n_doublets / nrow(doublet_results) * 100, 1)
    
    doublet_table <- table(doublet_results$Classification,
                           doublet_results$Cluster)
    
    # Add classifications back to seurat subset for plotting
    seurat_sub$scDblFinder_Classification <- doublet_results$Classification
    seurat_sub$scDblFinder_Score          <- doublet_results$Doublet_Score
    
    # ----------------------------------------
    # Plot 1: UMAP coloured by doublet classification
    # ----------------------------------------
    message("  Creating doublet UMAP plot...")
    
    p_doublet_umap <- DimPlot(seurat_sub,
                              group.by = "scDblFinder_Classification",
                              reduction = "umap",
                              cols = c("singlet" = "grey80", "doublet" = "red"),
                              pt.size = 0.8) +
      ggtitle(paste0("scDblFinder: Clusters ",
                     paste(doublet_clusters, collapse = ", "),
                     "\n", n_doublets, " Doublets (", doublet_pct, "%)")) +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
    
    ggsave(file.path(outdir, "16_scDblFinder_umap.pdf"),
           p_doublet_umap, width = 10, height = 8)
    
    # ----------------------------------------
    # Plot 2: Doublet score feature plot
    # ----------------------------------------
    message("  Creating doublet score plot...")
    
    p_score <- FeaturePlot(seurat_sub,
                           features = "scDblFinder_Score",
                           reduction = "umap",
                           pt.size = 0.8) +
      scale_color_viridis_c(option = "inferno") +
      ggtitle("scDblFinder: Doublet Score") +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
    
    ggsave(file.path(outdir, "17_scDblFinder_score.pdf"),
           p_score, width = 10, height = 8)
    
    # ----------------------------------------
    # Plot 3: Doublets per cluster bar chart
    # ----------------------------------------
    message("  Creating doublets per cluster bar chart...")
    
    doublet_per_cluster <- doublet_results %>%
      dplyr::group_by(Cluster, Classification) %>%
      dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(Cluster) %>%
      dplyr::mutate(Pct = round(Count / sum(Count) * 100, 1))
    
    p_doublet_bar <- ggplot(doublet_per_cluster,
                            aes(x = Cluster, y = Count, fill = Classification)) +
      geom_bar(stat = "identity", position = "stack") +
      geom_text(data = subset(doublet_per_cluster, Classification == "doublet"),
                aes(label = paste0(Pct, "%")),
                position = position_stack(vjust = 1.05),
                size = 4) +
      scale_fill_manual(values = c("singlet" = "grey70", "doublet" = "red")) +
      ggtitle("Doublets per Cluster") +
      labs(y = "Number of Cells", x = "Cluster") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    ggsave(file.path(outdir, "18_doublets_per_cluster_bar.pdf"),
           p_doublet_bar, width = 10, height = 6)
    
    # ----------------------------------------
    # Save doublet results as CSV
    # ----------------------------------------
    write.csv(doublet_results,
              file.path(outdir, "19_scDblFinder_results.csv"),
              row.names = FALSE)
    
    # ----------------------------------------
    # Save doublet summary TXT
    # ----------------------------------------
    sink(file.path(outdir, "20_scDblFinder_summary.txt"))
    cat("=======================================================\n")
    cat("         scDblFinder DOUBLET SUMMARY                   \n")
    cat("=======================================================\n\n")
    
    cat("PARAMETERS:\n")
    cat("  Clusters checked:    ", paste(doublet_clusters, collapse = ", "), "\n\n")
    
    cat("RESULTS:\n")
    cat("  Total cells checked: ", nrow(doublet_results), "\n")
    cat("  Singlets:            ", n_singlets, "\n")
    cat("  Doublets:            ", n_doublets, " (", doublet_pct, "%)\n\n")
    
    cat("DOUBLETS PER CLUSTER:\n")
    print(doublet_table)
    
    cat("\n\nRECOMMENDATION:\n")
    if (doublet_pct > 15) {
      cat("  ⚠️  High doublet rate detected (>15%). Consider removing these clusters\n")
      cat("     or filtering doublet-classified cells before downstream analysis.\n")
    } else if (doublet_pct > 7.5) {
      cat("  ⚠️  Elevated doublet rate (>7.5%). Review these clusters carefully.\n")
    } else {
      cat("  ✓  Doublet rate within expected range for 10x Genomics data.\n")
    }
    
    cat("\nNEXT STEPS:\n")
    cat("  To remove doublets from full object:\n")
    cat("  doublet_barcodes <- doublet_results$Cell_Barcode[doublet_results$Classification == 'doublet']\n")
    cat("  seurat_clean <- subset(seurat_obj, cells = setdiff(colnames(seurat_obj), doublet_barcodes))\n")
    cat("=======================================================\n")
    sink()
    
    message("  ✓ scDblFinder results saved")
    message("  → Doublets found: ", n_doublets, " / ", nrow(doublet_results),
            " cells (", doublet_pct, "%)")
    
  } else {
    message("\nSTEP 5: Doublet detection skipped (DoubletCheck = FALSE)")
    message("  💡 Tip: If you see suspicious clusters (high nFeature/nCount),")
    message("     re-run with DoubletCheck = TRUE, doublet_clusters = c(X, Y)")
  }
  
  # ========================================
  # Summary Statistics
  # ========================================
  message("  Writing summary statistics...")
  
  sink(file.path(outdir, "15_clustering_summary.txt"))
  cat("=======================================================\n")
  cat("         CLUSTERING SUMMARY                            \n")
  cat("=======================================================\n\n")
  
  cat("PARAMETERS:\n")
  cat("  Number of PCs used:", n_dims, "\n")
  cat("  Resolution:", resolution, "\n")
  cat("  tSNE performed:", run_tsne, "\n")
  cat("  DoubletCheck run:", DoubletCheck, "\n\n")
  
  cat("RESULTS:\n")
  cat("  Total cells:", ncol(seurat_obj), "\n")
  cat("  Number of clusters:", n_clusters, "\n\n")
  
  cat("CELLS PER CLUSTER:\n")
  for (i in 1:nrow(cluster_df)) {
    cat(sprintf("  Cluster %s: %d cells (%.1f%%)\n",
                cluster_df$Cluster[i],
                cluster_df$Count[i],
                cluster_df$Percentage[i]))
  }
  
  cat("\n=======================================================\n")
  cat("NEXT STEPS:\n")
  cat("  1. Find marker genes: FindAllMarkers(seurat_obj)\n")
  cat("  2. Annotate cell types based on markers\n")
  if (DoubletCheck) {
    cat("  3. Review doublet results in 20_scDblFinder_summary.txt\n")
    cat("  4. Save object: saveRDS(seurat_obj, 'final.rds')\n")
  } else {
    cat("  3. Save object: saveRDS(seurat_obj, 'final.rds')\n")
  }
  cat("=======================================================\n")
  sink()
  
  # ========================================
  # Console Summary
  # ========================================
  message("\n========================================")
  message("✅ ANALYSIS COMPLETE!")
  message("========================================")
  message("\n📊 RESULTS:")
  message("   • Used ", n_dims, " PCs")
  message("   • Found ", n_clusters, " clusters")
  message("   • Analyzed ", ncol(seurat_obj), " cells")
  message("\n📁 Output directory: ", file.path(getwd(), outdir))
  message("\n📄 FILES CREATED:")
  message("   → 07_pca_with_clusters.pdf")
  message("   → 08_umap_clusters.pdf")
  if (run_tsne) message("   → 09_tsne_clusters.pdf")
  message("   → 11_umap_qc_metrics_overlay.pdf")
  message("   → 12_violin_by_cluster.pdf")
  message("   → 13_top_variable_genes_heatmap.pdf")
  message("   → 14_cells_per_cluster.pdf")
  message("   → 15_clustering_summary.txt")
  if (DoubletCheck) {
    message("   → 16_scDblFinder_umap.pdf")
    message("   → 17_scDblFinder_score.pdf")
    message("   → 18_doublets_per_cluster_bar.pdf")
    message("   → 19_scDblFinder_results.csv")
    message("   → 20_scDblFinder_summary.txt")
  }
  message("\n🔬 Cluster distribution:")
  for (i in 1:nrow(cluster_df)) {
    message(sprintf("   Cluster %s: %d cells (%.1f%%)",
                    cluster_df$Cluster[i],
                    cluster_df$Count[i],
                    cluster_df$Percentage[i]))
  }
  message("\n📚 NEXT: Find marker genes with FindAllMarkers()")
  message("========================================\n")
  
  return(seurat_obj)
}
