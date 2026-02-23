qc_diagnostic <- function(seurat_obj,
                               mt_col = "percent.mt",
                               outdir = "qc_pre") {
  
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  
  dir.create(outdir, showWarnings = FALSE)
  
  meta <- seurat_obj@meta.data
  
  # -----------------------------
  # 1️⃣ Violin plots
  # -----------------------------
  vln <- VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", mt_col),
    ncol = 3,
    pt.size = 0.1
  )
  
  ggsave(file.path(outdir, "01_violin.pdf"),
         vln, width = 12, height = 5)
  
  # -----------------------------
  # 2️⃣ Histograms + density
  # -----------------------------
  h_feat <- ggplot(meta, aes(nFeature_RNA)) +
    geom_histogram(aes(y = ..density..), bins = 100, fill = "grey80") +
    geom_density(color = "red") +
    ggtitle("nFeature_RNA Distribution")
  
  h_count <- ggplot(meta, aes(nCount_RNA)) +
    geom_histogram(aes(y = ..density..), bins = 100, fill = "grey80") +
    geom_density(color = "red") +
    scale_x_log10() +
    ggtitle("nCount_RNA Distribution (log scale)")
  
  h_mt <- ggplot(meta, aes_string(mt_col)) +
    geom_histogram(aes(y = ..density..), bins = 100, fill = "grey80") +
    geom_density(color = "red") +
    ggtitle(paste(mt_col, "Distribution"))
  
  ggsave(file.path(outdir, "02_histograms_density.pdf"),
         h_feat + h_count + h_mt,
         width = 14, height = 5)
  
  # -----------------------------
  # 3️⃣ Scatter plots
  # -----------------------------
  sc1 <- FeatureScatter(seurat_obj,
                        feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA") +
    ggtitle("nCount vs nFeature (Doublet detection)")
  
  sc2 <- FeatureScatter(seurat_obj,
                        feature1 = "nCount_RNA",
                        feature2 = mt_col) +
    ggtitle("nCount vs Mitochondrial %")
  
  ggsave(file.path(outdir, "03_scatter_plots.pdf"),
         sc1 + sc2,
         width = 10, height = 5)
  
  # -----------------------------
  # 4️⃣ Log-scale feature plot
  # -----------------------------
  log_feat <- ggplot(meta,
                     aes(x = log10(nCount_RNA),
                         y = log10(nFeature_RNA))) +
    geom_point(alpha = 0.3) +
    ggtitle("log10(nCount) vs log10(nFeature)")
  
  ggsave(file.path(outdir, "04_log_feature_relationship.pdf"),
         log_feat,
         width = 6, height = 5)
  
  message("QC diagnostic plots saved in: ", outdir)
  
  return(seurat_obj)
}
