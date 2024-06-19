calculate_cross_K_gradient3D <- function(spe, 
                                         reference_cell_type, 
                                         target_cell_type, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = 2))
  colnames(result) <- c("observed_cross_K", 
                        "expected_cross_K")
  
  for (radius in seq(radii)) {
    cross_K_data <- calculate_cross_K3D(spe,
                                        reference_cell_type,
                                        target_cell_type,
                                        radius,
                                        feature_colname)
    
    result[radius, ] <- cross_K_data
  }
  
  # Add a radius column to the result
  result$radius <- seq(radii)
  
  if (plot_image) {
    plot(result$radius, result$observed_cross_K, type = "l", col = "red", 
         xlim = c(0, radius), ylim = c(0, max(result)),
         xlab = "Radius", ylab = "Cross K-function value")
    lines(result$radius, result$expected_cross_K, type = "l", col = "blue", lty = 2)
    legend(0, max(result), legend = c("Observed cross K", "Expected cross K"), col = c("red", "blue"), lty = c(1, 2))
  }
  
  return(result)
}