calculate_cross_K_gradient3D <- function(spe, 
                                         reference_cell_type, 
                                         target_cell_type, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = 3))
  colnames(result) <- c("observed_cross_K", 
                        "expected_cross_K",
                        "cross_K_ratio")
  
  for (radius in seq(radii)) {
    cross_K_df <- calculate_cross_K3D(spe,
                                      reference_cell_type,
                                      target_cell_type,
                                      radius,
                                      feature_colname)
    
    result[radius, ] <- cross_K_df
  }
  
  # Add a radius column to the result
  result$radius <- seq(radii)
  
  if (plot_image) plot_cross_K_gradient3D(result, reference_cell_type, target_cell_type)
  
  return(result)
}