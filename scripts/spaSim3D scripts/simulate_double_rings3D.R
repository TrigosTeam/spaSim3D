simulate_double_rings3D <- function(bg_sample,
                             bg_type = "Others",
                             n_dr = 1,
                             dr_properties = list(
                               D1 = list(
                                 name_of_cluster_cell = "Tumour",
                                 infiltration_types = c("Immune1", "Others"),
                                 infiltration_proportions = c(0.1, 0.05),
                                 shape = "Sphere",
                                 radius = 35,
                                 centre_loc = c(50, 50, 50),
                                 name_of_inner_ring_cell = "Immune1",
                                 inner_ring_width = 6,
                                 inner_ring_infiltration_types = c("Others"),
                                 inner_ring_infiltration_proportions = c(0.15),
                                 name_of_outer_ring_cell = "Immune2",
                                 outer_ring_width = 3,
                                 outer_ring_infiltration_types = c("Others"),
                                 outer_ring_infiltration_proportions = c(0.15)
                               )
                             ),
                             plot_image = TRUE,
                             plot_categories = c("Others", "Tumour", "Immune1", "Immune2"),
                             plot_colours = c("lightgray", "orange", "blue", "green")) {
  
  for (k in seq_len(n_dr)) { 
    
    # for each cluster, get the shape
    shape <- dr_properties[[k]]$shape
    
    
    ### Sphere shape
    if (shape == "Sphere") {
      bg_sample <- simulate_sphere_dr(bg_sample = bg_sample, dr_properties = dr_properties[[k]])
    } 
    
    ### Ellipsoid shape
    if (shape == "Ellipsoid") {
      bg_sample <- simulate_ellipsoid_dr(bg_sample = bg_sample, dr_properties = dr_properties[[k]])
    }
  }
  
  if (plot_image) {
    plot <- plot_cell_categories3D(bg_sample,
                                   cell_types_of_interest = plot_categories,
                                   colour_vector = plot_colours)
    print(plot)
    
  }
  
  return(bg_sample)
}
