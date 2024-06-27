calculate_cell_proportions_in_neighbourhood3D <- function(spe, 
                                                          reference_cell_type, 
                                                          target_cell_types, 
                                                          radius, 
                                                          feature_colname = "Cell.Type") {
  
  ## Get 'count' neighbourhood data
  cells_in_neighbourhood_data <- calculate_cells_in_neighbourhood3D(spe,
                                                                    reference_cell_type,
                                                                    target_cell_types,
                                                                    radius,
                                                                    feature_colname,
                                                                    FALSE,
                                                                    FALSE)

    
  
  result <- data.frame(matrix(nrow = length(target_cell_types), ncol = 4))
  colnames(result) <- c("target_cell_type", "frequency", "proportion", "percentage")
  
  result$target_cell_type <- target_cell_types
  
  ## Get frequency of each target cell type
  result$frequency <- apply(cells_in_neighbourhood_data[ , target_cell_types], 2, sum)
  
  ## Use frequency to get proportion and percentage of each cell type
  total <- sum(result$frequency)
  if (total != 0) {
    result$proportion <- result$frequency / total
    result$percentage <- result$proportion * 100  
  }
  else {
    result$proportion <- NA
    result$percentage <- NA
  }
  
  
  return(result)
}
