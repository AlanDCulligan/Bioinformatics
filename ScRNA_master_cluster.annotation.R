# ScRNA Seuquencing Worflow [For single / multi sample datasets]
all_times <- list()  # store the time for each chunk
knitr::knit_hooks$set(time_it = local({
  now <- NULL
  function(before, options) {
    if (before) {
      now <<- Sys.time()
    } else {
      res <- difftime(Sys.time(), now, units = "secs")
      all_times[[options$label]] <<- res
    }
  }
}))
knitr::opts_chunk$set(
  tidy = TRUE,
  tidy.opts = list(width.cutoff = 95),
  message = FALSE,
  warning = FALSE,
  time_it = TRUE,
  error = TRUE
)

# Set working directory to project folder. Change this to your project folder.
setwd("~/Desktop/Genetic_Analysis/1_ScRNA/1_Seurat_Guide_Clust_Tut/filtered_gene_bc_matrices/hg19")

# Load Packages
packages <- c("dplyr","Seurat","patchwork","sctransform","glmGamPoi","tidyr",
  "AnnotationHub","ensembldb", "multtest","tibble","openxlsx","HGNChelper",
  "tuple", "ggplot2","ggridges","SeuratWrappers", "scDblFinder", "harmony", "presto"
)
for (pkg in packages) {
  library(pkg, character.only = TRUE)
}

# Load in data 
pbmc.data = Read10X(data.dir = "~/Desktop/Genetic_Analysis/1_ScRNA/1_Seurat_Guide_Clust_Tut/filtered_gene_bc_matrices/hg19") # Load whichever dataset you want to work with. This is the 3k PBMC dataset from 10X Genomics. You can also load in your own data here.
pbmc = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3,
                          min.features = 200 )
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc

# Create backup
pbmc.b = pbmc

# Check class and structure of Seurat object
class(pbmc)
str(pbmc, max.level = 1)
is(pbmc, "Seurat")


# Quality Check 
source("qc_diagnostic.R")   
pbmc_postqc = qc_diagnostic(seurat_obj = pbmc)

# Subset out low quality cells (change settings below as needed for your dataset). This is informed by the QC diagnostic plots.

# pbmc_postqc <- subset(pbmc_postqc, subset = nFeature_RNA > 200 & 
#                       nFeature_RNA < 2500 & 
#                       percent.mt < 5)

# Normalise and PCA. This normalises using SCTransform (look at function to choose)
source("normalise_and_pcs.R")
pbmc_norm = normalise_and_pcs(seurat_obj = pbmc_postqc)

# Clustering and UMAP / tSNE with PC dimensions obtained form previous step.
# Check params in the function script!!!
source("run_clustering.R")
pbmc_postcluster <- run_clustering(
  seurat_obj = pbmc_norm,
  n_dims = 10,
  resolution = 0.5,
  outdir = "qc_post"
)

# [This integration function is a WORK IN PROGRESS]
# --- If we have multiple samples in our dataset we can correct for batch with Harmony or seurat integration. This is not necessary for this dataset but the code is provided in case you have multiple samples. --- #
# If this is the case this script will follow the harmony or Seurat integration workflow and will ignore following marker and annotation steps.

# source("batch_correction.R")
# pbmc_batch_corrected <- batch_correction(
#   seurat_obj = pbmc_postcluster,
#   method = "harmony", # choose between "harmony" and "seurat_integration"
#   batch_var = "orig.ident" # change this to the name of the metadata variable that contains the batch information (e.g. sample name, batch number, etc.)NU

# Now establish markers and annotate clusters. This outputs top 10 markers per cluster and visualisations to assist with annotation.
# Change step 9 settings for tissue type cannonical markers if using such markers.
source("find_markers.R")
marker_results <- find_markers(
  seurat_obj = pbmc_postcluster,
  test_use = "wilcox",      # Fast with presto
  min_pct = 0.25,
  logfc_threshold = 0.25,
  only_pos = TRUE,
  top_n_genes = 10,
  outdir = "marker_analysis"
)

# Annotate clusters (change settings for tissue type cannonical markers if using such markers). 
source("annotate_clusters.R")
pbmc_annotated = annotate_clusters(
  seurat_obj = pbmc_postcluster,
)

