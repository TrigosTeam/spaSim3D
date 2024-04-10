calculate_cells_in_neighborhood3D <- function(data, 
                                              reference_cell_types, 
                                              target_cell_types, 
                                              radius = 20, 
                                              feature_colname = "Cell.Type") {
  
  # If the columns are not correct, give error
  required_colnames <- c("Cell.X.Position", 
                         "Cell.Y.Position", 
                         "Cell.Z.Position", 
                         feature_colname,
                         "Cell.ID")
  
  missing_colnames <- setdiff(required_colnames,
                              colnames(data))
  
  if (length(missing_colnames) > 0) {
    stop(paste(paste(missing_colnames, collapse = ', '),
               "are missing as column names in your data")) 
  }
  
  
  # Check if reference_cell_types has cells not found in the data
  incorrect_cell_types <- setdiff(reference_cell_types, unique(data[[feature_colname]]))
  if (length(incorrect_cell_types) > 0) {
    stop(paste(paste(incorrect_cell_types, collapse = ', '),
               "in reference_cell_types don't existin data."))
  }
  
  # Check if target_cell_types has cells not found in the data
  incorrect_cell_types <- setdiff(target_cell_types, unique(data[[feature_colname]]))
  if (length(incorrect_cell_types) > 0) {
    stop(paste(paste(incorrect_cell_types, collapse = ', '),
               "in target_cell_types don't exist in data."))
  }
  
  # Check if radius is numeric
  if (!is.numeric(radius)) {
    stop(paste(radius, " is not of type 'numeric'"))
  }
  
  
  
  result <- list()
  
  for (reference_cell_type in reference_cell_types) {
    ## Get data for reference cells
    reference_cells <- data[which(data[, feature_colname] == reference_cell_type), ]
    reference_cell_coords <- reference_cells[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
    
    ## Set up data frame for current reference cell type
    reference_cell_df <- reference_cells[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
    rownames(reference_cell_df) <- reference_cells$Cell.ID
    
    for (target_cell_type in target_cell_types) {
      ## Get data for target cells
      target_cells <- data[which(data[, feature_colname] == target_cell_type), ]
      target_cell_coords <- target_cells[, c("Cell.X.Position","Cell.Y.Position", "Cell.Z.Position")]
      
      ## Determine number of target cells specified distance for each reference cell
      reference_target_result <- dbscan::frNN(target_cell_coords, 
                                              eps = radius,
                                              query = reference_cell_coords, 
                                              sort = FALSE)
      n_targets <- rapply(reference_target_result$id, length)
      
      ## Add to results data frame
      reference_cell_df[, target_cell_type] <- n_targets
      
    }
    result[[reference_cell_type]] <- reference_cell_df
  }
  
  return (result)
}
