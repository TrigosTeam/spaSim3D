add_spe_metadata3D <- function(spe, metadata, plot_image = TRUE) {
  
  # Ignore the 'background' element in metadata
  metadata[['background']] <- NULL
  
  for (i in seq(length(metadata))) {
    metadata_cluster <- metadata[[i]]
    
    if (metadata_cluster$cluster_type == "regular") {
      spe <- simulate_clusters3D(spe, list(metadata_cluster), plot_image = plot_image)
    }
    else if (metadata_cluster$cluster_type == "ring") {
      spe <- simulate_rings3D(spe, list(metadata_cluster), plot_image = plot_image)
    }
    else if (metadata_cluster$cluster_type == "double ring") {
      spe <- simulate_double_rings3D(spe, list(metadata_cluster), plot_image = plot_image)
    }
  }
  
  return(spe)
}