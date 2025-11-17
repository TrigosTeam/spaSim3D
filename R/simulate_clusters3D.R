simulate_clusters3D <- function(spe,
                                cluster_properties_list,
                                plot_image = TRUE,
                                plot_cell_types = NULL,
                                plot_colours = NULL) {
  
  # Check shape variable of cluster_properties
  shapes <- sapply(cluster_properties_list, function(x) {return(x[["shape"]])})
  n_invalid_shapes <- sum(!(shapes %in% c("sphere", "ellipsoid", "cylinder", "network")))
  if (n_invalid_shapes > 0) {
    stop("`cluster_properties_list` contains invalid shape parameters or no shape parameters.")
  }
  
  for (i in seq(length(cluster_properties_list))) { 
    
    shape <- shapes[[i]]
    
    ### Sphere shape
    if (shape == "sphere") {
      spe <- simulate_sphere_cluster(spe, cluster_properties_list[[i]])
    } 
    
    ### Ellipsoid shape
    if (shape == "ellipsoid") {
      spe <- simulate_ellipsoid_cluster(spe, cluster_properties_list[[i]])
    }
    
    ### Cylinder shape
    if (shape == "cylinder") {
      spe <- simulate_cylinder_cluster(spe, cluster_properties_list[[i]])
    }
    
    ### Network shape
    if (shape == "network") {
      spe <- simulate_network_cluster(spe, cluster_properties_list[[i]])
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe, 
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(spe)
}
