calculate_mixing_scores3D <- function(data, 
                                      reference_cell_types, 
                                      target_cell_types, 
                                      radius = 20, 
                                      feature_colname = "Cell.Type") {
  
  # If the columns are not correct, give error
  required_colnames <- c("Cell.X.Position", 
                         "Cell.Y.Position", 
                         "Cell.Z.Position", 
                         feature_colname)
  
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
  
  
  df <- data.frame(matrix(ncol=8, nrow=0))
  
  for (reference_cell_type in reference_cell_types) {
    
    # Get all info for cells of reference cell_type
    reference_cells <- data[data[, feature_colname] == reference_cell_type, ]
    
    for (target_cell_type in target_cell_types) {
      
      # Get all info for cells of target cell_type      
      target_cells <- data[data[, feature_colname] == target_cell_type, ]
      
      # No point getting mixing scores if comparing the same cell type
      if (reference_cell_type == target_cell_type) {
        next
      }
      
      # Can't get mixing scores if there are no reference cells
      if (nrow(reference_cells) == 0) {
        methods::show(paste("There are no unique reference cells of specified celltype", reference_cell_type, "for target cell", target_cell_type))
        df <-  rbind(df, 
                     c(reference_cell_type, 
                       target_cell_type, 
                       0, 
                       nrow(target_cells), 
                       0, 
                       0, 
                       NA, 
                       NA))
      }
      
      # Can't get mixing scores if there are no target cells
      else if (nrow(target_cells) == 0) {
        methods::show(paste("There are no unique target cells of specified cell type", target_cell_type, "for reference cell", reference_cell_type))
        
        reference_cell_coords <- reference_cells[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
        reference_reference_result <- dbscan::frNN(reference_cell_coords, 
                                                   eps = radius, 
                                                   query = NULL,
                                                   sort = FALSE)
        
        # halve it to avoid counting each ref-ref interaction twice
        reference_reference_interactions <- 0.5 * sum(rapply(reference_reference_result$id, length)) 
        
        df <-  rbind(df[ , df.cols], 
                     c(reference_cell_type, 
                       target_cell_type, 
                       nrow(reference_cells), 
                       0, 
                       0, 
                       reference_reference_interactions, 
                       NA, 
                       NA))
      }
      
      # Generic case: We have reference cells and target cells
      else {
        # Get x,y,z coords for reference cells and target cells
        reference_cell_coords <- reference_cells[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
        target_cell_coords <- target_cells[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
        
        # For each reference cell, find all target cells within the chosen radius
        reference_target_result <- dbscan::frNN(target_cell_coords, 
                                                eps = radius, 
                                                query = reference_cell_coords, 
                                                sort = FALSE)
        
        # Find the total sum of how many target cells were close enough to reference cells
        reference_target_interactions <- sum(rapply(reference_target_result$id, length))
        
        # For each reference cell, find all other reference cells within the chosen radius
        reference_reference_result <- dbscan::frNN(reference_cell_coords, 
                                                   eps = radius,
                                                   query = NULL,
                                                   sort = FALSE)
        
        # Find the the total sum of how many other reference cells were close enough to reference cells
        # Halve it to avoid counting each ref-ref interaction twice
        reference_reference_interactions <- 0.5 * sum(rapply(reference_reference_result$id, length)) 
        
        
        if (reference_reference_interactions != 0) {
          mixing_score <- reference_target_interactions / reference_reference_interactions
          normalised_mixing_score <- 0.5 * mixing_score * (nrow(reference_cells) - 1) / nrow(target_cells)
        }
        else {
          mixing_score <- 0
          normalised_mixing_score <- 0
          methods::show(paste("There are no reference to reference interactions for", target_cell_type, "in the specified radius, cannot calculate mixing score"))
        }
        
        df <-  rbind(df, 
                     c(reference_cell_type, 
                       target_cell_type, 
                       nrow(reference_cells), 
                       nrow(target_cells), 
                       reference_target_interactions, 
                       reference_reference_interactions, 
                       mixing_score, 
                       normalised_mixing_score))
      }
    }
  }
  
  # Required column names of our output data frame
  df.cols <- c("Reference", 
               "Target", 
               "Number_of_reference_cells",
               "Number_of_target_cells", 
               "Reference_target_interaction",
               "Reference_reference_interaction", 
               "Mixing_score", 
               "Normalised_mixing_score")
  colnames(df) <- df.cols
  
  return(df)
}
