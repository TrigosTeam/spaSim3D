simulate_spe_metadata3D <- function(spe_metadata, plot_image = TRUE) {
  
  # First element should contain background metadata
  bg_metadata <- spe_metadata[[1]]
  if (!(is.character(bg_metadata$background_type) && length(bg_metadata$background_type) == 1)) {
    stop("background_type parameter found in the metadata background list is not a character.")
  }
  
  if (bg_metadata$background_type == "random") {
    spe <- simulate_random_background_cells3D(bg_metadata$n_cells,
                                              bg_metadata$length,
                                              bg_metadata$width,
                                              bg_metadata$height,
                                              bg_metadata$minimum_distance_between_cells,
                                              plot_image = F)    
  }
  else if (bg_metadata$background_type == "ordered") {
    spe <- simulate_ordered_background_cells3D(bg_metadata$n_cells,
                                              bg_metadata$length,
                                              bg_metadata$width,
                                              bg_metadata$height,
                                              bg_metadata$jitter_proportion,
                                              plot_image = F) 
  }
  else {
    stop("background_type parameter found in the first list must be 'random' or 'ordered'.")
  }
  # Apply background mixing
  spe <- simulate_mixing3D(spe,
                           bg_metadata$cell_types,
                           bg_metadata$cell_proportions,
                           plot_image = F)
  
  ### If there is only background metadata, we are done
  if (length(spe_metadata) == 1) {
    
    # Plot
    if (plot_image) {
      fig <- plot_cells3D(spe)
      methods::show(fig)
    }
    
    return(spe)
  }
  
  ### All other elements should help to simulate clusters 
  for (i in 2:length(spe_metadata)) {
    cluster_metadata <- spe_metadata[[i]]
    
    if (!(is.character(cluster_metadata$cluster_type) && length(cluster_metadata$cluster_type) == 1)) {
      stop(paste("cluster_type parameter found in the metadata cluster list", i,"is not a character."))
    }
    
    if (cluster_metadata$cluster_type == "regular") {
      spe <- simulate_clusters3D(spe, list(cluster_metadata), plot_image = F)
    }
    else if (cluster_metadata$cluster_type == "ring") {
      spe <- simulate_rings3D(spe, list(cluster_metadata), plot_image = F)      
    }
    else if (cluster_metadata$cluster_type == "double ring") {
      spe <- simulate_double_rings3D(spe, list(cluster_metadata), plot_image = F)
    }
    else {
      stop("cluster_type parameter must be either 'regular', 'ring' or 'double ring'.")
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe)
    methods::show(fig)
  }
  
  return(spe)
}
