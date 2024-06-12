simulate_clusters3D <- function(bg_sample,
                                n_clusters = 3,
                                cluster_properties = list(
                                  C1 = list(
                                    shape = "Sphere",
                                    cluster_cell_types = c("Tumour", "Immune", "Others"),
                                    cluster_cell_proportions = c(0.55, 0.4, 0.05),
                                    radius = 25,
                                    centre_loc = c(40, 40, 40)
                                  ),
                                  C2 = list(
                                    shape = "Cylinder",
                                    cluster_cell_types = c("Endothelial", "Others"),
                                    cluster_cell_proportions = c(0.95, 0.05),
                                    radius = 10,
                                    start_loc = c(0, 0, 0),
                                    end_loc   = c(20, 20 , 100)
                                  ),
                                  C3 = list(
                                    shape = "Ellipsoid",
                                    cluster_cell_types = c("Tumour", "Immune", "Others"),
                                    cluster_cell_proportions = c(0.65, 0.3, 0.05),
                                    x_radius = 15,
                                    y_radius = 20,
                                    z_radius = 25,
                                    centre_loc = c(70, 70, 70),
                                    x_y_rotation = 0,
                                    x_z_rotation = 0,
                                    y_z_rotation = 0
                                  )
                                ),
                                plot_image = TRUE,
                                plot_categories = c("Others", "Immune", "Endothelial", "Tumour"),
                                plot_colours = c("lightgray", "skyblue", "#FF7F7F", "orange")) {
  
  
  for (k in seq_len(n_clusters)) { 
    
    # For each cluster, get the shape
    shape <- cluster_properties[[k]]$shape
    
    ### Sphere shape
    if (shape == "Sphere") {
      bg_sample <- simulate_sphere_cluster(bg_sample = bg_sample, cluster_properties = cluster_properties[[k]])
    } 
    
    ### Ellipsoid shape
    if (shape == "Ellipsoid") {
      bg_sample <- simulate_ellipsoid_cluster(bg_sample = bg_sample, cluster_properties = cluster_properties[[k]])
    }
    
    ### Cylinder shape
    if (shape == "Cylinder") {
      bg_sample <- simulate_cylinder_cluster(bg_sample = bg_sample, cluster_properties = cluster_properties[[k]])
    }
    
    ### Network shape
    if (shape == "Network") {
      bg_sample <- simulate_network_cluster(bg_sample = bg_sample, cluster_properties = cluster_properties[[k]])
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
