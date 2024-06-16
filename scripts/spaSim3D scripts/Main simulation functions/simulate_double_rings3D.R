simulate_double_rings3D <- function(bg_spe,
                                    dr_properties,
                                    plot_image = TRUE,
                                    plot_cell_types = NULL,
                                    plot_colours = NULL) {
  
  for (k in seq(length(dr_properties))) { 
    
    # For each cluster, get the shape
    shape <- dr_properties[[k]]$shape
    
    ### Sphere double ring shape
    if (shape == "Sphere") {
      bg_spe <- simulate_sphere_dr(bg_spe, dr_properties[[k]])
    } 
    
    ### Ellipsoid double ring shape
    if (shape == "Ellipsoid") {
      bg_spe <- simulate_ellipsoid_dr(bg_spe, dr_properties[[k]])
    }
    
    ### Cylinder double ring shape
    if (shape == "Cylinder") {
      bg_spe <- simulate_cylinder_dr(bg_spe, dr_properties[[k]])
    }
    
    ### Network double ring shape
    if (shape == "Network") {
      bg_spe <- simulate_network_dr(bg_spe, dr_properties[[k]])
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
