calculate_cells_in_neighborhood3D <- function(data, 
                                              reference_cell_types, 
                                              target_cell_types, 
                                              radius = 20, 
                                              feature_colname = "Cell.Type") {
  
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
