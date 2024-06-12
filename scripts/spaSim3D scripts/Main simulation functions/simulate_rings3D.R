simulate_rings3D <- function(bg_sample,
                             n_ring = 3,
                             ring_properties = list(
                               R1 = list(
                                 shape = "Sphere",
                                 cluster_cell_types = c("Tumour", "Others"),
                                 cluster_cell_proportions = c(0.95, 0.05),
                                 radius = 20,
                                 centre_loc = c(40, 40, 40),
                                 ring_cell_types = c("Immune", "Others"),
                                 ring_cell_proportions = c(0.85, 0.15),
                                 ring_width = 5
                               ),
                               R2 = list(
                                 shape = "Cylinder",
                                 cluster_cell_types = c("Void"),
                                 cluster_cell_proportions = c(1),
                                 radius = 8,
                                 start_loc = c(0, 0, 0),
                                 end_loc   = c(20, 20 , 100),
                                 ring_cell_types = c("Endothelial", "Others"),
                                 ring_cell_proportions = c(0.85, 0.15),
                                 ring_width = 5
                               ),
                               R3 = list(
                                 shape = "Ellipsoid",
                                 cluster_cell_types = c("Tumour", "Others"),
                                 cluster_cell_proportions = c(0.95, 0.05),
                                 x_radius = 10,
                                 y_radius = 15,
                                 z_radius = 20,
                                 centre_loc = c(70, 70, 70),
                                 x_y_rotation = 0,
                                 x_z_rotation = 0,
                                 y_z_rotation = 0,
                                 ring_cell_types = c("Immune", "Others"),
                                 ring_cell_proportions = c(0.85, 0.15),
                                 ring_width = 5
                               )
                             ),
                             plot_image = TRUE,
                             plot_categories = c("Others", "Tumour", "Immune", "Endothelial"),
                             plot_colours = c("lightgray", "orange", "skyblue", "#FF7F7F")) {
  
  for (k in seq_len(n_ring)) { 
    
    # For each cluster, get the shape
    shape <- ring_properties[[k]]$shape
    
    ### Sphere shape +  ring
    if (shape == "Sphere") {
      bg_sample <- simulate_sphere_ring(bg_sample = bg_sample, ring_properties = ring_properties[[k]])
    } 
    
    ### Ellipsoid shape + ring
    else if (shape == "Ellipsoid") {
      bg_sample <- simulate_ellipsoid_ring(bg_sample = bg_sample, ring_properties = ring_properties[[k]])
    }
    
    ### Cylinder shape + ring
    else if (shape == "Cylinder") {
      bg_sample <- simulate_cylinder_ring(bg_sample = bg_sample, ring_properties = ring_properties[[k]])
    }
    
    ### Network shape + ring
    else if (shape == "Network") {
      bg_sample <- simulate_network_ring(bg_sample = bg_sample, ring_properties = ring_properties[[k]])
    }
    
    else {
      stop("Invalid shape")
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cell_categories3D(bg_sample,
                                  cell_types_of_interest = plot_categories,
                                  colour_vector = plot_colours,
                                  size = 2,
                                  include_cell_types_of_no_interest = FALSE,
                                  feature_colname = "Cell.Type")
    print(fig)
  }
  
  return (bg_sample)
}
