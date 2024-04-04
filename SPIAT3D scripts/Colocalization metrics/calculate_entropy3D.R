calculate_entropy3D <- function(data,
                                radius = NULL,
                                reference_cell_type = NULL,
                                target_cell_types,
                                feature_colname = "Cell.Type") {
 
  entropy <- 0
  
  # Calculate entropy of the entire image
  if (is.null(radius)) {
    
    ## Get data for chosen target cells
    data <- data[which(data[, feature_colname] == target_cell_types), ]
    
    n_all_cell_types <- nrow(data)
    
    ## No cells found, return 0
    if (n_all_cell_types == 0) {
      return (0)
    }
    
    for (target_cell_type in target_cell_types) {
      
      ## Get data for current target cell 
      target_cell_type_data <- data[which(data[, feature_colname] == target_cell_type), ]
      n_target_cell_type <- nrow(target_cell_type_data)
      
      ## No cells found for current target cell, move on
      if (n_target_cell_type == 0) {
        next
      }
      
      target_cell_proportion <- n_target_cell_type / n_all_cell_types
      entropy <- entropy + (-1 * (target_cell_proportion) * log2(target_cell_proportion))
    }
    
    return (entropy)
  }
  
  ## Radius has been specified, calculate entropy for chosen reference cell
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  
  cells_in_neighborhood_data <- calculate_cells_in_neighborhood3D(data,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radius,
                                                                  feature_colname)[[1]]
  
  ## Get total number of target cells for each row
  cells_in_neighborhood_data$Total <- apply(cells_in_neighborhood_data[target_cell_types], 1, sum)
  
  ## Get entropy for each row
  cells_in_neighborhood_data$Entropy <- 0
  
  for (target_cell_type in target_cell_types) {
    
    target_cell_type_proportions <- (cells_in_neighborhood_data[target_cell_type] / cells_in_neighborhood_data$Total)[[1]]

    ## If an element in target_cell_type_proportion is 0, just add 0.    
    cells_in_neighborhood_data$Entropy <- cells_in_neighborhood_data$Entropy +
                                          (-1 * target_cell_type_proportions * ifelse(target_cell_type_proportions == 0,
                                                                                      0, log2(target_cell_type_proportions)))
    
  }
  
  ## Case when row has 0 target cells
  cells_in_neighborhood_data[cells_in_neighborhood_data$Total == 0, "Entropy"] <- 0
  
  return (cells_in_neighborhood_data)
}
