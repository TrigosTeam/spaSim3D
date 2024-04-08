simulate_clusters3D <- function(bg_sample,
                                bg_type = "Others",
                                n_clusters = 3,
                                cluster_properties = list(
                                  C1 = list(
                                    name_of_cluster_cell = "Tumour",
                                    infiltration_types = c("Immune", "Others"),
                                    infiltration_proportions = c(0.4, 0.05),
                                    shape = "Sphere",
                                    radius = 25,
                                    centre_loc = c(40, 40, 40)
                                  ),
                                  C2 = list(
                                    name_of_cluster_cell = "Endothelial",
                                    infiltration_types = c("Others"),
                                    infiltration_proportions = c(0.05),
                                    shape = "Cylinder",
                                    radius = 10,
                                    start_loc = c(0, 0, 0),
                                    end_loc   = c(20, 20 , 100)
                                  ),
                                  C3 = list(
                                    name_of_cluster_cell = "Tumour",
                                    infiltration_types = c("Immune", "Others"),
                                    infiltration_proportions = c(0.3, 0.05),
                                    shape = "Ellipsoid",
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
    
    # for each cluster, get the shape
    shape <- cluster_properties[[k]]$shape
    
    # # generate a location as the centre of the cluster
    # if (is.null(centre_loc)){
    #   seed_point <- spatstat.random::runifpoint(1, win=win)}
    # else seed_point <- centre_loc
    # a <- seed_point$x
    # b <- seed_point$y
    
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
  }
  
  if (plot_image) {
    plot <- plot_cell_categories3D(bg_sample,
                                   cell_types_of_interest = plot_categories,
                                   colour_vector = plot_colours)
    print(plot)
    
  }
  
  return(bg_sample)
}
