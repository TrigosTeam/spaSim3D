calculate_entropy_background3D <- function(spe,
                                           cell_types_of_interest = NULL, 
                                           feature_colname = "Cell.Type") {
  
  
  cell_proportions_data <- calculate_cell_proportions3D(spe, cell_types_of_interest, feature_colname, FALSE, FALSE)
  
  # Calculate entropy of the entire image
  cell_proportions <- cell_proportions_data$proportion
  entropy <- -1 * sum(cell_proportions * log(cell_proportions, length(cell_proportions)))

  return(entropy) 
}