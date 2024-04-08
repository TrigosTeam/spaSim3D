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
