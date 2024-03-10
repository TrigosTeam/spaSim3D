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
                             plot_colours = NULL) {
  
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
  
  if (plot_image){
    if(is.null(plot_categories)) plot_categories <- unique(bg_sample$Cell.Type)
    if (is.null(plot_colours)){
      plot_colours <- c("gray","darkgreen", "red", "darkblue", "brown", "purple", "lightblue",
                        "lightgreen", "yellow", "black", "pink")
    }
    phenos <- plot_categories
    
    colors <- c()
    for (i in 1:nrow(bg_sample)) {
      for (j in 1:length(phenos)) {
        if (bg_sample$Cell.Type[i] == phenos[j]) {
          colors <- append(colors, plot_colours[j])
          break
        }
      }
    }
    
    plot3d(bg_sample$Cell.X.Position,
           bg_sample$Cell.Y.Position,
           bg_sample$Cell.Z.Position,
           xlab = "x",
           ylab = "y",
           zlab = "z",
           col = colors,
           size = 4)
    
    # add legend
    legend3d("topright", legend = phenos, pch = 16, col = plot_colours[seq_len(length(phenos))], inset = c(0.02))
  }
  
  return(bg_sample)
}
