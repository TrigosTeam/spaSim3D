simulate_rings3D <- function(bg_sample,
                             bg_type = "Others",
                             n_ring = 1,
                             ring_properties = list(
                               R1 = list(
                               name_of_cluster_cell = "Tumour",
                               infiltration_types = c("Immune1", "Others"),
                               infiltration_proportions = c(0.1, 0.05),
                               shape = "Sphere",
                               radius = 35,
                               centre_loc = c(50, 50, 50),
                               name_of_ring_cell = "Immune1",
                               ring_width = 5,
                               ring_infiltration_types = c("Others"),
                               ring_infiltration_proportions = c(0.15))
                             ),
                             plot_image = TRUE,
                             plot_categories = c("Others", "Tumour", "Immune1"),
                             plot_colours = NULL) {
  
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
