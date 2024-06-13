simulate_double_rings3D <- function(bg_sample,
                                    n_dr = 3,
                                    dr_properties = list(
                                      D1 = list(
                                        shape = "Sphere",
                                        cluster_cell_types = c("Tumour", "Others"),
                                        cluster_cell_proportions = c(0.95, 0.05),
                                        radius = 20,
                                        centre_loc = c(40, 40, 40),
                                        inner_ring_cell_types = c("Immune1", "Others"),
                                        inner_ring_cell_proportions = c(0.85, 0.15),
                                        inner_ring_width = 5,
                                        outer_ring_cell_types = c("Immune2"),
                                        outer_ring_cell_proportions = c(1),
                                        outer_ring_width = 3
                                      ),
                                      D2 = list(
                                        shape = "Cylinder",
                                        cluster_cell_types = c("Void"),
                                        cluster_cell_proportions = c(1),
                                        radius = 8,
                                        start_loc = c(0, 0, 0),
                                        end_loc   = c(20, 20 , 100),
                                        inner_ring_cell_types = c("Endothelial", "Others"),
                                        inner_ring_cell_proportions = c(0.85, 0.15),
                                        inner_ring_width = 5,
                                        outer_ring_cell_types = c("Immune2"),
                                        outer_ring_cell_proportions = c(1),
                                        outer_ring_width = 3
                                      ),
                                      D3 = list(
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
                                        inner_ring_cell_types = c("Immune1", "Others"),
                                        inner_ring_cell_proportions = c(0.85, 0.15),
                                        inner_ring_width = 5,
                                        outer_ring_cell_types = c("Immune2"),
                                        outer_ring_cell_proportions = c(1),
                                        outer_ring_width = 3
                                      )
                                    ),
                                    plot_image = TRUE,
                                    plot_categories = c("Others", "Tumour", "Immune1", "Immune2", "Endothelial"),
                                    plot_colours = c("lightgray", "orange", "skyblue", "blue", "#FF7F7F")) {
  
  for (k in seq_len(n_dr)) { 
    
    # For each cluster, get the shape
    shape <- dr_properties[[k]]$shape
    
    ### Sphere double ring shape
    if (shape == "Sphere") {
      bg_sample <- simulate_sphere_dr(bg_sample = bg_sample, dr_properties = dr_properties[[k]])
    } 
    
    ### Ellipsoid double ring shape
    if (shape == "Ellipsoid") {
      bg_sample <- simulate_ellipsoid_dr(bg_sample = bg_sample, dr_properties = dr_properties[[k]])
    }
    
    ### Cylinder double ring shape
    if (shape == "Cylinder") {
      bg_sample <- simulate_cylinder_dr(bg_sample = bg_sample, dr_properties = dr_properties[[k]])
    }
    
    ### Network double ring shape
    if (shape == "Network") {
      bg_sample <- simulate_network_dr(bg_sample = bg_sample, dr_properties = dr_properties[[k]])
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
