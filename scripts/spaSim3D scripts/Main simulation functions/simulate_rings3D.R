simulate_rings3D <- function(bg_spe,
                             ring_properties,
                             plot_image = TRUE,
                             plot_cell_types = NULL,
                             plot_colours = NULL) {
  
  for (k in seq(length(ring_properties))) { 
    
    # For each cluster, get the shape
    shape <- ring_properties[[k]]$shape
    
    ### Sphere shape +  ring
    if (shape == "Sphere") {
      bg_spe <- simulate_sphere_ring(bg_spe, ring_properties[[k]])
    } 
    
    ### Ellipsoid shape + ring
    else if (shape == "Ellipsoid") {
      bg_spe <- simulate_ellipsoid_ring(bg_spe, ring_properties[[k]])
    }
    
    ### Cylinder shape + ring
    else if (shape == "Cylinder") {
      bg_spe <- simulate_cylinder_ring(bg_spe, ring_properties[[k]])
    }
    
    ### Network shape + ring
    else if (shape == "Network") {
      bg_spe <- simulate_network_ring(bg_spe, ring_properties[[k]])
    }
    
    else {
      stop("Invalid shape")
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
