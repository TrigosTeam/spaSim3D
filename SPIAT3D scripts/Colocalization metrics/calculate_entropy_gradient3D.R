calculate_entropy_gradient3D <- function(data,
                                         radii,
                                         reference_cell_type,
                                         target_cell_types,
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  entropy_gradient <- list()
  entropy_mean <- c()
  
  for (radius in seq(1, radii, 0.5)) {
    entropy_data <- calculate_entropy3D(data,
                                        radius,
                                        reference_cell_type,
                                        target_cell_types,
                                        feature_colname)
    
    entropy_gradient[[paste(radius)]] <- entropy_data
    entropy_mean <- append(entropy_mean, mean(entropy_data$Entropy))
    
  }
  
  if (plot_image) {
    plot(seq(1, radii, 0.5), entropy_mean, type = "l", xlab = "Radius", ylab = "Entropy Mean")
  }
  
  return(entropy_gradient)
  
}
