determine_entropy_prevalence3D <- function(entropy_grid_data,
                                           threshold,
                                           above = TRUE) {
  
  if (above) {
    p <- sum(entropy_grid_data$Entropy >= threshold) / nrow(entropy_grid_data) * 100
  }
  else {
    p <- sum(entropy_grid_data$Entropy < threshold) / nrow(entropy_grid_data) * 100    
  }
  
  return (p)
}
