simulate_double_rings3D <- function(spe,
                                    dr_properties_list,
                                    plot_image = TRUE,
                                    plot_cell_types = NULL,
                                    plot_colours = NULL) {
  
  # Check shape variable of dr_properties_list
  shapes <- sapply(dr_properties_list, function(x) {return(x[["shape"]])})
  n_invalid_shapes <- sum(!(shapes %in% c("sphere", "ellipsoid", "cylinder", "network")))
  if (n_invalid_shapes > 0) {
    stop("`dr_properties_list` contains invalid shape parameters or no shape parameters.")
  }
  
  for (i in seq(length(dr_properties_list))) { 
    
    shape <- shapes[[i]]
    
    ### Sphere shape with double ring
    if (shape == "sphere") {
      spe <- simulate_sphere_dr(spe, dr_properties_list[[i]])
    } 
    
    ### Ellipsoid shape with double ring
    if (shape == "ellipsoid") {
      spe <- simulate_ellipsoid_dr(spe, dr_properties_list[[i]])
    }
    
    ### Cylinder shape with double ring
    if (shape == "cylinder") {
      spe <- simulate_cylinder_dr(spe, dr_properties_list[[i]])
    }
    
    ### Network shape with double ring
    if (shape == "network") {
      spe <- simulate_network_dr(spe, dr_properties_list[[i]])
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
