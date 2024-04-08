simulate_rings3D <- function(bg_sample,
                             bg_type = "Others",
                             n_ring = 3,
                             ring_properties = list(
                               R1 = list(
                                 name_of_cluster_cell = "Tumour",
                                 infiltration_types = c("Others"),
                                 infiltration_proportions = c(0.05),
                                 shape = "Sphere",
                                 radius = 20,
                                 centre_loc = c(40, 40, 40),
                                 name_of_ring_cell = "Immune",
                                 ring_width = 5,
                                 ring_infiltration_types = c("Others"),
                                 ring_infiltration_proportions = c(0.15)
                               ),
                               R2 = list(
                                 name_of_cluster_cell = "Void",
                                 infiltration_types = NULL,
                                 infiltration_proportions = NULL,
                                 shape = "Cylinder",
                                 radius = 8,
                                 start_loc = c(0, 0, 0),
                                 end_loc   = c(20, 20 , 100),
                                 name_of_ring_cell = "Endothelial",
                                 ring_width = 5,
                                 ring_infiltration_types = c("Others"),
                                 ring_infiltration_proportions = c(0.15)
                               ),
                               R3 = list(
                                 name_of_cluster_cell = "Tumour",
                                 infiltration_types = c("Others"),
                                 infiltration_proportions = c(0.05),
                                 shape = "Ellipsoid",
                                 x_radius = 10,
                                 y_radius = 15,
                                 z_radius = 20,
                                 centre_loc = c(70, 70, 70),
                                 x_y_rotation = 0,
                                 x_z_rotation = 0,
                                 y_z_rotation = 0,
                                 name_of_ring_cell = "Immune",
                                 ring_width = 5,
                                 ring_infiltration_types = c("Others"),
                                 ring_infiltration_proportions = c(0.15)
                               )
                             ),
                             plot_image = TRUE,
                             plot_categories = c("Others", "Tumour", "Immune", "Endothelial"),
                             plot_colours = c("lightgray", "orange", "skyblue", "#FF7F7F")) {
  
  for (k in seq_len(n_ring)) { 
    
    # for each cluster, get the shape
    shape <- ring_properties[[k]]$shape
    
    
    ### Sphere shape + immune ring
    if (shape == "Sphere") {
      bg_sample <- simulate_sphere_ring(bg_sample = bg_sample, ring_properties = ring_properties[[k]])
    } 
    
    ### Ellipsoid shape + ring
    if (shape == "Ellipsoid") {
      bg_sample <- simulate_ellipsoid_ring(bg_sample = bg_sample, ring_properties = ring_properties[[k]])
    }
    
    ### Cylinder shape + ring
    if (shape == "Cylinder") {
      bg_sample <- simulate_cylinder_ring(bg_sample = bg_sample, ring_properties = ring_properties[[k]])
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
  
  return (bg_sample)
}
