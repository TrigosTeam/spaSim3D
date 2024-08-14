calculate_cells_in_neighbourhood_gradient3D <- function(spe, 
                                                        reference_cell_type, 
                                                        target_cell_types, 
                                                        radii, 
                                                        feature_colname = "Cell.Type",
                                                        plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (radius in seq(radii)) {
    cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                   reference_cell_type,
                                                                   target_cell_types,
                                                                   radius,
                                                                   feature_colname,
                                                                   FALSE,
                                                                   FALSE)
    
    cells_in_neighbourhood_df$ref_cell_id <- NULL
    result[radius, ] <- apply(cells_in_neighbourhood_df, 2, mean)
  }
  # Add a radius column to the result
  result$radius <- seq(radii)
  
  if (plot_image) plot_cells_in_neighbourhood_gradient3D(result, reference_cell_type)
  
  return(result)
}
