calculate_entropy_gradient3D <- function(spe,
                                         reference_cell_type,
                                         target_cell_types,
                                         radii,
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = 1))
  colnames(result) <- "entropy"
  
  for (radius in seq(radii)) {
    entropy_df <- calculate_entropy3D(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radius,
                                      feature_colname,
                                      FALSE)
    
    result[radius, "entropy"] <- mean(entropy_df$entropy)
  }
  
  # Add a radius column to the result
  result$radius <- seq(radii)
  
  if (plot_image) {
    expected_entropy <- calculate_entropy_background3D(spe, target_cell_types, feature_colname)
    plot_entropy_gradient3D(result, expected_entropy)
  }
  
  return(result)
}
