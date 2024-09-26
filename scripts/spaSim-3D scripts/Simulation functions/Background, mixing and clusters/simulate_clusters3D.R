simulate_clusters3D <- function(bg_spe,
                                cluster_properties,
                                plot_image = TRUE,
                                plot_cell_types = NULL,
                                plot_colours = NULL) {
  
  
  for (k in seq(length(cluster_properties))) { 
    
    # For each cluster, get the shape
    shape <- cluster_properties[[k]]$shape
    
    ### Sphere shape
    if (shape == "sphere") {
      bg_spe <- simulate_sphere_cluster(bg_spe, cluster_properties[[k]])
    } 
    
    ### Ellipsoid shape
    if (shape == "ellipsoid") {
      bg_spe <- simulate_ellipsoid_cluster(bg_spe, cluster_properties[[k]])
    }
    
    ### Cylinder shape
    if (shape == "cylinder") {
      bg_spe <- simulate_cylinder_cluster(bg_spe, cluster_properties[[k]])
    }
    
    ### Network shape
    if (shape == "network") {
      bg_spe <- simulate_network_cluster(bg_spe, cluster_properties[[k]])
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(bg_spe, 
                        plot_cell_types,
                        plot_colours)
    print(fig)
  }
  
  return(bg_spe)
}
