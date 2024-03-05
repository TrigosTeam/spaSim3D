simulate_clusters3D <- function(bg_sample,
                                n_clusters = 2,
                                bg_type = "Others",
                                cluster_properties = list(
                                  C1 = list(
                                    name_of_cluster_cell = "Tumour",
                                    infiltration_types = c("Immune1", "Others"),
                                    infiltration_proportions = c(0.1, 0.05),
                                    shape = "Sphere",
                                    radius = 50,
                                    centre_loc = c(50, 50, 50)),
                                  C2 = list(
                                    name_of_cluster_cell = "Endo",
                                    infiltration_types = c("Immune1", "Others"),
                                    infiltration_proportions = c(0.1, 0.05),
                                    shape = "Cylinder",
                                    radius = 10,
                                    start_loc = c(0, 0, 0),
                                    end_loc   = c(20, 40 ,60)
                                  )
                                ),
                                plot_image = TRUE,
                                plot_categories = NULL,
                                plot_colours = NULL) {
  
  
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
      
    }
    
    ### Cylinder shape
    if (shape == "Cylinder") {
      bg_sample <- simulate_cylinder_cluster(bg_sample = bg_sample, cluster_properties = cluster_properties[[k]])
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
