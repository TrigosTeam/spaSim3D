simulate_spe_metadata3D <- function(spe_metadata) {
  
  # First element should contain background metadata
  bg_metadata <- spe_metadata[[1]]
  if (bg_metadata$background_type == "random") {
    spe <- simulate_random_background_cells3D(bg_metadata$n_cells,
                                              bg_metadata$length,
                                              bg_metadata$width,
                                              bg_metadata$height,
                                              bg_metadata$minimum_distance_between_cells)    
  }
  else if (bg_metadata$background_type == "normal") {
    spe <- simulate_normal_background_cells3D(bg_metadata$n_cells,
                                              bg_metadata$length,
                                              bg_metadata$width,
                                              bg_metadata$height,
                                              bg_metadata$jitter_proportion) 
  }
  else {
    stop("background_type parameter found in the first list must be 'random' or 'normal'.")
  }
  # Apply background mixing
  spe <- simulate_mixing3D(spe,
                           bg_metadata$cell_types,
                           bg_metadata$cell_proportions)
  
  ### If there is only background metadata, we are done
  if (length(spe_metadata) == 1) return(spe)
  
  
  ### All other elements should help to simulate clusters 
  for (i in 2:length(spe_metadata)) {
    cluster_metadata <- spe_metadata[[i]]
    if (cluster_metadata$cluster_type == "regular") {
      spe <- simulate_clusters3D(spe, list(cluster_metadata))
    }
    else if (cluster_metadata$cluster_type == "ring") {
      spe <- simulate_rings3D(spe, list(cluster_metadata))      
    }
    else if (cluster_metadata$cluster_type == "double ring") {
      spe <- simulate_double_rings3D(spe, list(cluster_metadata))
    }
  }
  
  return(spe)
}