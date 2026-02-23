# ============================================================================
# FUNCTION: Annotate Clusters with Marker-Based Suggestions
# ============================================================================
# Purpose: Suggest cell types based on top markers, then rename and save
# ============================================================================

annotate_clusters <- function(seurat_obj,
                              cluster_names = NULL,
                              marker_file = "marker_analysis/top10_markers_per_cluster.csv",
                              auto_suggest = TRUE,
                              tissue_type = "PBMC",
                              species = "human",
                              outdir = "annotation_results",
                              save_prefix = "annotated") {
  
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  message("\n========================================")
  message("ANNOTATING CLUSTERS")
  message("========================================")
  
  n_clusters <- length(unique(Idents(seurat_obj)))
  message("\nTotal clusters: ", n_clusters)
  
  # ========================================
  # STEP 1: Load Top 10 Markers and Suggest Cell Types
  # ========================================
  
  if (auto_suggest && is.null(cluster_names)) {
    message("\n🤖 Auto-suggesting cell types based on top 10 markers...")
    
    # Load top 10 markers file
    if (file.exists(marker_file)) {
      message("  ✓ Loading: ", marker_file)
      top10_markers <- read.csv(marker_file)
      
      # Validate columns
      if (!("cluster" %in% names(top10_markers)) || !("gene" %in% names(top10_markers))) {
        # Try column indices if names don't match
        if (ncol(top10_markers) >= 7) {
          colnames(top10_markers)[6] <- "cluster"
          colnames(top10_markers)[7] <- "gene"
        } else {
          stop("❌ Cannot find cluster (col 6) and gene (col 7) columns in marker file")
        }
      }
      
      message("  ✓ Loaded ", nrow(top10_markers), " markers")
      
    } else {
      stop("❌ Marker file not found: ", marker_file, 
           "\n   Please run find_markers() first or provide correct path")
    }
    
    # Get top 3 markers per cluster for display
    top3_per_cluster <- top10_markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_head(n = 3) %>%
      dplyr::summarize(
        top_marker = gene[1],
        top3_markers = paste(gene, collapse = ", "),
        .groups = "drop"
      )
    
    # Define marker-to-cell-type mapping (PBMC specific)
    if (tissue_type == "PBMC") {
      marker_database <- list(
        # T cells
        "CCR7" = "Naive CD4 T",
        "TCF7" = "Naive CD4 T",
        "LEF1" = "Naive CD4 T",
        "IL7R" = "Naive CD4 T",
        "S100A4" = "Memory CD4 T",
        "LTB" = "Memory CD4 T",
        "GZMK" = "Memory CD4 T",
        "CD8A" = "CD8 T",
        "CD8B" = "CD8 T",
        "GZMH" = "CD8 T",
        "GZMB" = "Cytotoxic CD8 T",
        "PRF1" = "Cytotoxic CD8 T",
        "NKG7" = "Cytotoxic T/NK",
        
        # B cells
        "MS4A1" = "B",
        "CD79A" = "B",
        "CD79B" = "B",
        "VPREB3" = "B",
        "CD19" = "B",
        "IGLL5" = "B",
        
        # Monocytes
        "CD14" = "CD14+ Mono",
        "LYZ" = "CD14+ Mono",
        "S100A9" = "CD14+ Mono",
        "S100A8" = "CD14+ Mono",
        "AQP3" = "CD14+ Mono",
        "FCGR3A" = "FCGR3A+ Mono",
        "MS4A7" = "FCGR3A+ Mono",
        "CDKN1C" = "FCGR3A+ Mono",
        "CX3CR1" = "FCGR3A+ Mono",
        "FOLR3" = "CD16+ Mono",
        
        # NK cells
        "GNLY" = "NK",
        "GZMB" = "NK/Cytotoxic T",
        "NCAM1" = "NK",
        
        # Dendritic cells
        "FCER1A" = "DC",
        "CST3" = "DC",
        "CD1C" = "DC",
        
        # Platelets
        "PPBP" = "Platelet",
        "PF4" = "Platelet",
        "GNG11" = "Platelet",
        "HIST1H2BJ" = "Platelet",
        "HIST1H2AC" = "Platelet"
      )
      
      # Score each cluster based on top 10 markers
      suggested_names <- c()
      
      message("\n📋 Suggested annotations based on top markers:\n")
      
      for (i in 1:nrow(top3_per_cluster)) {
        clust <- top3_per_cluster$cluster[i]
        top_marker <- top3_per_cluster$top_marker[i]
        top3_list <- top3_per_cluster$top3_markers[i]
        
        # Get all top 10 markers for this cluster
        cluster_top10 <- top10_markers %>%
          dplyr::filter(cluster == clust) %>%
          dplyr::pull(gene)
        
        # Score against marker database
        scores <- list()
        for (cell_type in unique(marker_database)) {
          # Get all markers for this cell type
          type_markers <- names(marker_database)[marker_database == cell_type]
          # Calculate overlap
          overlap <- length(intersect(cluster_top10, type_markers))
          if (overlap > 0) {
            scores[[cell_type]] <- overlap
          }
        }
        
        # Get best match
        if (length(scores) > 0) {
          best_match <- names(which.max(scores))
          confidence <- max(unlist(scores))
        } else {
          # Fallback to single top marker
          if (top_marker %in% names(marker_database)) {
            best_match <- marker_database[[top_marker]]
            confidence <- 1
          } else {
            best_match <- paste0("Unknown_", clust)
            confidence <- 0
          }
        }
        
        suggested_names <- c(suggested_names, best_match)
        
        message(sprintf("  Cluster %s: %-20s (confidence: %d/10, top3: %s)",
                       clust, best_match, confidence, top3_list))
      }
      
      names(suggested_names) <- sort(unique(Idents(seurat_obj)))
      cluster_names <- suggested_names
      
      # Save suggestions to file
      suggestion_df <- data.frame(
        Cluster = names(suggested_names),
        SuggestedType = suggested_names,
        TopMarkers = top3_per_cluster$top3_markers
      )
      write.csv(suggestion_df,
                file.path(outdir, "suggested_annotations.csv"),
                row.names = FALSE)
      message("\n  ✓ Suggestions saved to: ", file.path(outdir, "suggested_annotations.csv"))
      
    } else {
      message("  ⚠️  Auto-suggestion only available for PBMC tissue type")
      message("  Please provide manual cluster_names")
      return(seurat_obj)
    }
  }
  
  # ========================================
  # STEP 2: Validate and Apply Annotations
  # ========================================
  
  if (is.null(cluster_names)) {
    stop("❌ No cluster names provided. Set cluster_names or auto_suggest = TRUE")
  }
  
  current_clusters <- levels(Idents(seurat_obj))
  
  if (length(cluster_names) != length(current_clusters)) {
    stop("❌ Number of names (", length(cluster_names), 
         ") doesn't match number of clusters (", length(current_clusters), ")")
  }
  
  message("\n🏷️  Applying annotations...")
  
  # Apply annotations
  names(cluster_names) <- current_clusters
  seurat_obj <- RenameIdents(seurat_obj, cluster_names)
  seurat_obj$cell_type <- Idents(seurat_obj)
  
  message("  ✓ Annotations applied")
  
  # Print mapping
  message("\n📊 Final Cluster → Cell Type mapping:")
  for (i in 1:length(cluster_names)) {
    n_cells <- sum(seurat_obj$cell_type == cluster_names[i])
    message(sprintf("  Cluster %s → %-20s (%d cells)", 
                   names(cluster_names)[i], 
                   cluster_names[i],
                   n_cells))
  }
  
  # ========================================
  # STEP 3: Create Visualizations
  # ========================================
  message("\n🎨 Creating publication-quality plots...")
  
  # Simple UMAP with labels (like tutorial)
  p_simple <- DimPlot(seurat_obj, 
                     reduction = 'umap', 
                     label = TRUE, 
                     pt.size = 0.5) + 
    NoLegend()
  
  ggsave(file.path(outdir, paste0(save_prefix, "_umap_simple.pdf")),
         p_simple, 
         width = 8, height = 7)
  
  # Publication-quality plot (like tutorial save.img)
  p_publication <- DimPlot(seurat_obj, 
                          reduction = "umap", 
                          label = TRUE, 
                          label.size = 4.5) + 
    xlab("UMAP 1") + 
    ylab("UMAP 2") + 
    ggtitle("Annotated Cell Types") +
    theme(axis.title = element_text(size = 18), 
          legend.text = element_text(size = 18),
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
    guides(colour = guide_legend(override.aes = list(size = 10)))
  
  # Save as high-quality image
  ggsave(
    filename = file.path(outdir, paste0(save_prefix, "_umap_publication.jpg")),
    plot = p_publication,
    height = 7, 
    width = 12, 
    quality = 50
  )
  
  # Also save as PDF
  ggsave(
    filename = file.path(outdir, paste0(save_prefix, "_umap_publication.pdf")),
    plot = p_publication,
    height = 7, 
    width = 12
  )
  
  # UMAP without labels (show legend)
  p_legend <- DimPlot(seurat_obj,
                     reduction = "umap",
                     label = FALSE,
                     pt.size = 0.5) +
    xlab("UMAP 1") + 
    ylab("UMAP 2") +
    ggtitle("Annotated Cell Types") +
    theme(axis.title = element_text(size = 18),
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  ggsave(file.path(outdir, paste0(save_prefix, "_umap_with_legend.pdf")),
         p_legend,
         width = 10, height = 7)
  
  # ========================================
  # STEP 4: Validate with Canonical Markers
  # ========================================
  message("\n✓ Creating marker validation plots...")
  
  if (tissue_type == "PBMC") {
    # Canonical PBMC markers
    canonical_markers <- c(
      "CD3D", "CD3E",      # T cells
      "CD4",               # CD4 T
      "CD8A",              # CD8 T
      "CCR7", "IL7R",      # Naive
      "GZMK",              # Memory
      "CD14", "LYZ",       # Monocytes
      "FCGR3A",            # Non-classical mono
      "MS4A1", "CD79A",    # B cells
      "GNLY", "NKG7",      # NK
      "FCER1A", "CST3",    # DC
      "PPBP"               # Platelets
    )
    
    # Keep only markers present in dataset
    canonical_markers <- canonical_markers[canonical_markers %in% rownames(seurat_obj)]
    
    if (length(canonical_markers) > 0) {
      # Dot plot
      p_canonical_dot <- DotPlot(seurat_obj, features = canonical_markers) +
        RotatedAxis() +
        ggtitle("Canonical Marker Expression") +
        theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      
      ggsave(file.path(outdir, paste0(save_prefix, "_canonical_markers_dotplot.pdf")),
             p_canonical_dot,
             width = max(12, length(canonical_markers) * 0.4),
             height = 8)
      
      # Feature plots
      p_canonical_feature <- FeaturePlot(seurat_obj,
                                         features = canonical_markers,
                                         ncol = 3)
      
      ggsave(file.path(outdir, paste0(save_prefix, "_canonical_markers_featureplot.pdf")),
             p_canonical_feature,
             width = 15,
             height = ceiling(length(canonical_markers) / 3) * 4)
    }
  }
  
  # ========================================
  # STEP 5: Cell Type Distribution
  # ========================================
  message("\n📊 Creating summary statistics...")
  
  # Count cells per type
  cell_counts <- table(seurat_obj$cell_type)
  summary_df <- data.frame(
    CellType = names(cell_counts),
    Count = as.numeric(cell_counts),
    Percentage = round(as.numeric(cell_counts) / ncol(seurat_obj) * 100, 2)
  )
  summary_df <- summary_df[order(-summary_df$Count), ]
  
  # Bar plot
  p_bar <- ggplot(summary_df, 
                 aes(x = reorder(CellType, Count), 
                     y = Count, 
                     fill = CellType)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")),
              hjust = -0.1, size = 3.5) +
    coord_flip() +
    ggtitle("Cell Type Distribution") +
    labs(x = "Cell Type", y = "Number of Cells") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))
  
  ggsave(file.path(outdir, paste0(save_prefix, "_celltype_distribution.pdf")),
         p_bar,
         width = 10, height = max(6, nrow(summary_df) * 0.5))
  
  # ========================================
  # STEP 6: Save Summary Files
  # ========================================
 message("\n💾 Saving annotation summary...")

sink(file.path(outdir, "ANNOTATION_SUMMARY.txt"))

cat("=======================================================\n")
cat("         CELL TYPE ANNOTATION SUMMARY                  \n")
cat("=======================================================\n\n")

print(summary_df)

cat("\nTotal cells:", ncol(seurat_obj), "\n")
cat("Total clusters:", length(unique(seurat_obj$cell_type)), "\n")

sink()  # ← CLOSE CONNECTION

message("  ✓ Summary saved")

return(seurat_obj)

}  # ← CLOSE FUNCTION




