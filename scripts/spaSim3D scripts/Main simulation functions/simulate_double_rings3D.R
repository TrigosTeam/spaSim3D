simulate_double_rings3D <- function(bg_sample,
                                    n_dr = 1,
                                    dr_properties = list(
                                      D1 = list(
                                        shape = "Sphere",
                                        cluster_cell_types = c("Tumour", "Immune1", "Others"),
                                        cluster_cell_proportions = c(0.1, 0.05),
                                        radius = 35,
                                        centre_loc = c(50, 50, 50),
                                        inner_ring_cell_types = c("Immune1", "Others"),
                                        inner_ring_cell_proportions = c(0.85, 0.15),
                                        inner_ring_width = 6,
                                        outer_ring_cell_types = c("Immune2", "Others"),
                                        outer_ring_cell_proportions = c(0.85, 0.15),
                                        outer_ring_width = 3
                                      )
                                    ),
                                    plot_image = TRUE,
                                    plot_categories = c("Others", "Tumour", "Immune1", "Immune2"),
                                    plot_colours = c("lightgray", "orange", "blue", "green")) {
  
  for (k in seq_len(n_dr)) { 
    
    # For each cluster, get the shape
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
  
  return(bg_sample)
}
