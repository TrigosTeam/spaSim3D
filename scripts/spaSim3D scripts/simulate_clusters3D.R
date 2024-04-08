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
    
    if (is.null(plot_categories)) {
      plot_categories <- unique(bg_sample$Cell.Type)
    }
    
    if (is.null(plot_colours)) {
      plot_colours <- c("red", "orange", "green", "blue", "skyblue", "pink", "purple", "lightgray")[1:length(plot_categories)]
    }
    
    bg_sample <- bg_sample[bg_sample[["Cell.Type"]] %in% plot_categories, ]
    
    ## Factor for feature column
    bg_sample[, "Cell.Type"] <- factor(bg_sample[, "Cell.Type"],
                                  levels = plot_categories)
    
    ## Plot
    fig <- plot_ly(bg_sample,
                   type = "scatter3d",
                   mode = 'markers',
                   x = ~Cell.X.Position,
                   y = ~Cell.Y.Position,
                   z = ~Cell.Z.Position,
                   color = ~Cell.Type,
                   colors = plot_colours,
                   marker = list(size = 2))
    
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                       yaxis = list(title = 'y'),
                                       zaxis = list(title = 'z')))
    print(fig)
  }
  
  return(bg_sample)
}
