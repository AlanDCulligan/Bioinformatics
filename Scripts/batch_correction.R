# --- This script is used to correct for batch effects in ScRNA Sequencing of multiple samples like (control vs treatment) --- #
batch_correction = function(seurat_obj, 
									method = NULL, # The method can be using Harmony (harmony) or Seurat Integration (seurat_integration)
									batch_var = NULL, # Change this to the metadata object with the batch information (e.g. sample name, batch number, etc)	
									) {
									
	# Load libraries for analysis
	library(harmony)

	# Selecting the method to use
	if (method == "harmony") {

	# Run harmomny 
	seurat_obj.harmony = seurat_obj %>%
		RunHarmony(group.by.vars = 'NULL', plot_convergence = FALSE)
	
	seurat_obj.harmony@reductions # See the new dimensionality reduction with Harmony embeddings for downstream analysis (e.g. clustering, UMAP, etc.)
	seurat_obj.harmony.embed = Embeddings(seurat_obj.harmony, 'harmony') # Extract the Harmony embeddings for downstream analysis (e.g. clustering, UMAP, etc.) 

	# Now run UMAP and clustering on the Harmony embeddings.
	seurat_obj.harmony %>%
		RunUMAP(reduction = 'harmony', dims = 1:20) %>%
		FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
		FindClusters(resolution = 0.5) 


	else if (method == "seurat_integration") {

		# Add Seurat Integration workflow here

	}
	}
	

	# Run Seurat Integration






	return(seurat_obj)	
}
