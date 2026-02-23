# ============================================================================
# FUNCTION 1: SCTransform Normalization → PCA → Elbow Plot
# ============================================================================
# Purpose: Normalize data and generate diagnostic plots to choose optimal PCs
# YOU MUST look at the elbow plot before proceeding to Function 2
# ============================================================================

normalise_and_pcs<- function(seurat_obj,
                                     mt_col = "percent.mt",
                                     vars_to_regress = c("percent.mt"), # Regress out mitochondrial percentage by default, but can add more (e.g., cell cycle scores)
                                     n_variable_features = 2000, # Choose how many variable features to use for PCA
                                     n_pcs_compute = 50,
                                     outdir = "qc_post") {
  
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  
  # Create output directory
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  message("========================================")
  message("STEP 1: SCTransform Normalization")
  message("========================================")
  
  # ========================================
  # SCTransform Normalization
  # ========================================
  seurat_obj <- SCTransform(
    seurat_obj,
    vars.to.regress = vars_to_regress,
    variable.features.n = n_variable_features,
    verbose = TRUE,
    return.only.var.genes = FALSE
  )
  
  message("\n========================================")
  message("STEP 2: Running PCA")
  message("========================================")
  
  # ========================================
  # Run PCA
  # ========================================
  seurat_obj <- RunPCA(
    seurat_obj,features = VariableFeatures(object = seurat_obj),
    verbose = FALSE
  )
  
  message("\n========================================")
  message("STEP 3: Generating Diagnostic Plots")
  message("========================================")
  
  # ========================================
  # Calculate variance metrics
  # ========================================
  stdev <- seurat_obj@reductions$pca@stdev
  var_explained <- stdev^2 / sum(stdev^2) * 100
  cumvar <- cumsum(var_explained)
  
  # Auto-detect suggestion (but user should verify!)
  elbow_point <- which(diff(stdev) < 0.1)[1]
  if (is.na(elbow_point)) elbow_point <- 15
  var_90 <- which(cumvar > 90)[1]
  suggested_pcs <- min(elbow_point + 5, 30)
  
  # ========================================
  # PLOT 1: Standard Elbow Plot
  # ========================================
  message("  Creating elbow plot...")
  
  p_elbow <- ElbowPlot(seurat_obj, ndims = n_pcs_compute) +
    ggtitle("Elbow Plot: Standard Deviation by PC") +
    labs(subtitle = paste0("Suggested: ~", suggested_pcs, " PCs (verify visually!)")) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "red")
    ) +
    geom_vline(xintercept = suggested_pcs, 
               linetype = "dashed", color = "red", size = 1) +
    annotate("text", 
             x = suggested_pcs + 3, 
             y = max(stdev[1:10]),
             label = paste("Suggested:", suggested_pcs, "PCs"),
             color = "red", size = 5)
  
  ggsave(filename = file.path(outdir, "01_ELBOW_PLOT_CHECK_THIS_FIRST.pdf"),
         plot = p_elbow, 
         width = 10, height = 6)
  
  # ========================================
  # PLOT 2: Cumulative Variance
  # ========================================
  message("  Creating cumulative variance plot...")
  
  var_df <- data.frame(
    PC = 1:length(var_explained),
    Variance_Pct = var_explained,
    Cumulative = cumvar
  )
  
  p_cumvar <- ggplot(var_df[1:n_pcs_compute,], aes(x = PC, y = Cumulative)) +
    geom_line(size = 1.5, color = "steelblue") +
    geom_point(size = 3, color = "steelblue") +
    geom_hline(yintercept = 90, linetype = "dashed", color = "red", size = 1) +
    annotate("text", x = n_pcs_compute * 0.7, y = 92, 
             label = "90% variance threshold", color = "red", size = 5) +
    labs(
      title = "Cumulative Variance Explained by PCs",
      subtitle = paste0("90% variance reached at PC ", var_90),
      x = "Principal Component",
      y = "Cumulative Variance Explained (%)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  
  ggsave(filename = file.path(outdir, "02_cumulative_variance.pdf"),
         plot = p_cumvar, 
         width = 10, height = 6)
  
  # ========================================
  # PLOT 3: Combined Elbow + Cumulative
  # ========================================
  message("  Creating combined plot...")
  
  p_combined <- p_elbow + p_cumvar
  
  ggsave(filename = file.path(outdir, "00_CHOOSE_NUMBER_OF_PCS.pdf"),
         plot = p_combined, 
         width = 18, height = 6)
  
  # ========================================
  # PLOT 4: Variable Features
  # ========================================
  message("  Creating variable features plot...")
  
  top10 <- head(VariableFeatures(seurat_obj), 10)
  p_varfeat <- VariableFeaturePlot(seurat_obj) +
    ggtitle("Top Variable Features") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  p_varfeat <- LabelPoints(plot = p_varfeat,
                           points = top10,
                           repel = TRUE)
  
  ggsave(filename = file.path(outdir, "03_variable_features.pdf"),
         plot = p_varfeat, 
         width = 10, height = 6)
  
  # ========================================
  # PLOT 5: PCA Visualization
  # ========================================
  message("  Creating PCA scatter plot...")
  
  # FIXED VERSION - use DimPlot instead
  p_pca <- DimPlot(seurat_obj, 
                   reduction = "pca",
                   dims = c(1, 2)) +
    ggtitle("PCA: PC1 vs PC2") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(filename = file.path(outdir, "04_pca_scatter.pdf"),
         plot = p_pca, 
         width = 10, height = 8)
  
  # ========================================
  # PLOT 6: DimHeatmap - First 9 PCs
  # ========================================
  # ========================================
  # PLOT 6: DimHeatmap - First 9 PCs
  # ========================================
  message("  Creating PCA heatmaps...")
  
  # Reduce complexity for better PDF rendering
  p_heatmap_1_9 <- DimHeatmap(seurat_obj, 
                              dims = 1:9, 
                              cells = 300,        # Reduced from 500
                              balanced = TRUE,
                              nfeatures = 20,     # Reduced from 30
                              fast = FALSE)       # Add this
  
  ggsave(filename = file.path(outdir, "05_pca_heatmap_PC1-9.pdf"),
         plot = p_heatmap_1_9, 
         width = 16,              # Reduced from 18
         height = 18,             # Reduced from 20
         limitsize = FALSE)       # Add this
  
  # PLOT 7: DimHeatmap - Next 9 PCs
  if (n_pcs_compute >= 18) {
    p_heatmap_10_18 <- DimHeatmap(seurat_obj, 
                                  dims = 10:18, 
                                  cells = 300,        # Reduced
                                  balanced = TRUE,
                                  nfeatures = 20,     # Reduced
                                  fast = FALSE)
    
    ggsave(filename = file.path(outdir, "05_pca_heatmap_PC10-18.pdf"),
           plot = p_heatmap_10_18, 
           width = 16,              # Reduced
           height = 16,             # Reduced
           limitsize = FALSE)
  }
  
  # ========================================
  # PLOT 8: PCA Loadings
  # ========================================
  message("  Creating PCA loadings plots...")
  
  p_load_1_2 <- VizDimLoadings(seurat_obj, 
                               dims = 1:2, 
                               reduction = "pca",
                               ncol = 2)
  
  ggsave(filename = file.path(outdir, "06_pca_loadings_PC1-2.pdf"),
         plot = p_load_1_2, 
         width = 12, height = 8)
  
  p_load_3_4 <- VizDimLoadings(seurat_obj, 
                               dims = 3:4, 
                               reduction = "pca",
                               ncol = 2)
  
  ggsave(filename = file.path(outdir, "06_pca_loadings_PC3-4.pdf"),
         plot = p_load_3_4, 
         width = 12, height = 8)
  
  # ========================================
  # GUIDE FILE: How to Choose PCs
  # ========================================
  message("  Writing PC selection guide...")
  
  sink(file.path(outdir, "PC_SELECTION_GUIDE.txt"))
  cat("=======================================================\n")
  cat("         HOW TO CHOOSE NUMBER OF PCs                   \n")
  cat("=======================================================\n\n")
  
  cat("📊 LOOK AT THESE FILES:\n")
  cat("   1. 00_CHOOSE_NUMBER_OF_PCS.pdf (main decision plot)\n")
  cat("   2. 01_ELBOW_PLOT_CHECK_THIS_FIRST.pdf\n")
  cat("   3. 05_pca_heatmap_PC*.pdf (when patterns become random)\n\n")
  
  cat("🔍 WHAT TO LOOK FOR:\n")
  cat("   1. ELBOW: Where the curve flattens/plateaus\n")
  cat("   2. CUMULATIVE VARIANCE: Typically want 80-90%\n")
  cat("   3. HEATMAPS: When gene patterns look less structured\n\n")
  
  cat("📈 AUTO-DETECTED METRICS:\n")
  cat("   Elbow detected at PC:", elbow_point, "\n")
  cat("   90% variance at PC:", var_90, "\n")
  cat("   SUGGESTED: Use", suggested_pcs, "PCs\n\n")
  
  cat("⚠️  IMPORTANT: These are suggestions only!\n")
  cat("   You should visually verify the elbow plot.\n\n")
  
  cat("📋 TYPICAL RANGES:\n")
  cat("   • Small/simple datasets: 10-15 PCs\n")
  cat("   • Medium datasets: 15-25 PCs\n")
  cat("   • Large/complex datasets: 25-40 PCs\n\n")
  
  cat("✅ NEXT STEP:\n")
  cat("   After reviewing plots, run:\n")
  cat("   \n")
  cat("   seurat_obj <- run_umap_clustering(\n")
  cat("     seurat_obj = seurat_obj,\n")
  cat("     n_dims = YOUR_CHOICE,  # e.g., 20\n")
  cat("     outdir = 'qc_post'\n")
  cat("   )\n\n")
  
  cat("=======================================================\n")
  cat("         VARIANCE BREAKDOWN (First 30 PCs)            \n")
  cat("=======================================================\n\n")
  
  for (i in 1:min(30, n_pcs_compute)) {
    cat(sprintf("PC%-2d: %5.2f%% | Cumulative: %5.2f%%\n", 
                i, var_explained[i], cumvar[i]))
  }
  
  sink()
  
  # ========================================
  # Print Summary to Console
  # ========================================
  message("\n========================================")
  message("✅ NORMALIZATION & PCA COMPLETE!")
  message("========================================")
  message("\n📁 Output directory: ", file.path(getwd(), outdir))
  message("\n📊 KEY FILES TO CHECK:")
  message("   → 00_CHOOSE_NUMBER_OF_PCS.pdf")
  message("   → PC_SELECTION_GUIDE.txt")
  message("\n💡 AUTO-SUGGESTION: Use ", suggested_pcs, " PCs")
  message("   (Based on elbow at PC", elbow_point, ")")
  message("\n⚠️  NEXT STEP:")
  message("   1. Open the PDF files in '", outdir, "'")
  message("   2. Look at the elbow plot")
  message("   3. Choose number of PCs (typically ", suggested_pcs, ")")
  message("   4. Run: run_umap_clustering(seurat_obj, n_dims = YOUR_CHOICE)")
  message("\n========================================\n")
  
  return(seurat_obj)
}