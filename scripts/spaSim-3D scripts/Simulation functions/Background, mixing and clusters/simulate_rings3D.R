simulate_rings3D <- function(spe,
                             ring_properties_list,
                             plot_image = TRUE,
                             plot_cell_types = NULL,
                             plot_colours = NULL) {
  
  # Check shape variable of ring_properties_list
  shapes <- sapply(ring_properties_list, function(x) {return(x[["shape"]])})
  n_invalid_shapes <- sum(!(shapes %in% c("sphere", "ellipsoid", "cylinder", "network")))
  if (n_invalid_shapes > 0) {
    stop("`ring_properties_list` contains invalid shape parameters or no shape parameters.")
  }
  
  for (i in seq(length(ring_properties_list))) { 
    
    shape <- shapes[[i]]
    
    ### Sphere shape +  ring
    if (shape == "sphere") {
      spe <- simulate_sphere_ring(spe, ring_properties_list[[i]])
    } 
    
    ### Ellipsoid shape + ring
    else if (shape == "ellipsoid") {
      spe <- simulate_ellipsoid_ring(spe, ring_properties_list[[i]])
    }
    
    ### Cylinder shape + ring
    else if (shape == "cylinder") {
      spe <- simulate_cylinder_ring(spe, ring_properties_list[[i]])
    }
    
    ### Network shape + ring
    else if (shape == "network") {
      spe <- simulate_network_ring(spe, ring_properties_list[[i]])
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
