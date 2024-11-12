add_spe_metadata3D <- function(spe, metadata, plot_image = TRUE) {
  
  # Ignore the 'background' element in metadata
  metadata[['background']] <- NULL
  
  for (i in seq(length(metadata))) {
    metadata_cluster <- metadata[[i]]
    
    if (!(is.character(metadata_cluster$cluster_type) && length(metadata_cluster$cluster_type) == 1)) {
      stop(paste("cluster_type parameter found in the metadata cluster list", i,"is not a character."))
    }
    
    if (metadata_cluster$cluster_type == "regular") {
      spe <- simulate_clusters3D(spe, list(metadata_cluster), plot_image = F)
    }
    else if (metadata_cluster$cluster_type == "ring") {
      spe <- simulate_rings3D(spe, list(metadata_cluster), plot_image = F)
    }
    else if (metadata_cluster$cluster_type == "double ring") {
      spe <- simulate_double_rings3D(spe, list(metadata_cluster), plot_image = F)
    }
    else {
      stop("cluster_type parameter must be either 'regular', 'ring' or 'double ring'.")
    }
  }
  
  if (plot_image) {
    fig <- plot_cells3D(spe)
    methods::show(fig)
  }
  
  return(spe)
}
