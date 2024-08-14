calculate_cells_in_neighbourhood_proportions_gradient3D <- function(spe, 
                                                                    reference_cell_type, 
                                                                    target_cell_types, 
                                                                    radii, 
                                                                    feature_colname = "Cell.Type",
                                                                    plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (radius in seq(radii)) {
    cell_proportions_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood_proportions3D(spe,
                                                                                                    reference_cell_type,
                                                                                                    target_cell_types,
                                                                                                    radius,
                                                                                                    feature_colname)
    
    result[radius, ] <- apply(cell_proportions_neighbourhood_proportions_df[ , paste(target_cell_types, "_prop", sep = "")], 2, mean)
  }
  
  # Add a radius column to the result
  result$radius <- seq(radii)
  
  # Plot
  if (plot_image) plot_cells_in_neighbourhood_proportions_gradient3D(result, reference_cell_type)
  
  return(result)
}
