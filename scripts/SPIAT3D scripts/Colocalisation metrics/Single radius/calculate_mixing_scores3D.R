calculate_mixing_scores3D <- function(spe, 
                                      reference_cell_types, 
                                      target_cell_types, 
                                      radius, 
                                      feature_colname = "Cell.Type") {
  
  if (is.null(spe[[feature_colname]])) stop(paste("No column called", feature_colname, "found in spe object"))
  
  
  ## For reference_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(reference_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in reference_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  ## For target_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  # Check if radius is numeric
  if (!is.numeric(radius)) stop(paste(radius, " is not of type 'numeric'"))
  
  # Get spe coords
  spe_coords <- spatialCoords(spe)
  
  # Define result
  result <- data.frame()
  
  for (reference_cell_type in reference_cell_types) {
    
    # Get coords for reference_cell_type
    reference_cell_type_coords <- spe_coords[spe[[feature_colname]] == reference_cell_type, ]
    
    for (target_cell_type in target_cell_types) {
      
      # Get coords for target_cell_type
      target_cell_type_coords <- spe_coords[spe[[feature_colname]] == target_cell_type, ]
      
      # No point getting mixing scores if comparing the same cell type
      if (reference_cell_type == target_cell_type) {
        next
      }
      
      # Can't get mixing scores if there are no reference cells
      if (nrow(reference_cell_type_coords) == 0) {
        methods::show(paste("There are no unique reference cells of specified cell type ", reference_cell_type, "for target cell", target_cell_type))
        result <-  rbind(result, 
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
      else if (nrow(target_cell_type_coords) == 0) {
        methods::show(paste("There are no unique target cells of specified cell type", target_cell_type, "for reference cell", reference_cell_type))
        
        ref_ref_result <- dbscan::frNN(reference_cell_type_coords, 
                                       eps = radius, 
                                       query = NULL,
                                       sort = FALSE)
        
        # halve it to avoid counting each ref-ref interaction twice
        n_ref_ref_interactions <- 0.5 * sum(rapply(ref_ref_result$id, length)) 
        
        result <-  rbind(result, 
                         c(reference_cell_type, 
                           target_cell_type, 
                           nrow(reference_cells), 
                           0, 
                           0, 
                           n_ref_ref_interactions, 
                           NA, 
                           NA))
      }
      
      # Generic case: We have reference cells and target cells
      else {
        
        # For each reference cell, find all target cells within the chosen radius
        ref_tar_result <- dbscan::frNN(target_cell_type_coords, 
                                       eps = radius, 
                                       query = reference_cell_type_coords, 
                                       sort = FALSE)
        
        # Find the total sum of how many target cells were close enough to reference cells
        n_ref_tar_interactions <- sum(rapply(ref_tar_result$id, length))
        
        # For each reference cell, find all other reference cells within the chosen radius
        ref_ref_result <- dbscan::frNN(reference_cell_type_coords, 
                                       eps = radius,
                                       query = NULL,
                                       sort = FALSE)
        
        # Find the the total sum of how many other reference cells were close enough to reference cells
        # Halve it to avoid counting each ref-ref interaction twice
        n_ref_ref_interactions <- 0.5 * sum(rapply(ref_ref_result$id, length)) 
        
        
        if (n_ref_ref_interactions != 0) {
          mixing_score <- n_ref_tar_interactions / n_ref_ref_interactions
          normalised_mixing_score <- 0.5 * mixing_score * (nrow(reference_cell_type_coords) - 1) / nrow(target_cell_type_coords)
        }
        else {
          mixing_score <- 0
          normalised_mixing_score <- 0
          methods::show(paste("There are no reference to reference interactions for", target_cell_type, "in the specified radius, cannot calculate mixing score"))
        }
        
        result <-  rbind(result, 
                         c(reference_cell_type, 
                           target_cell_type, 
                           nrow(reference_cell_type_coords), 
                           nrow(target_cell_type_coords), 
                           n_ref_tar_interactions, 
                           n_ref_ref_interactions, 
                           mixing_score, 
                           normalised_mixing_score))
      }
    }
  }
  
  # Required column names of our output data frame
  colnames(result) <- c("ref_cell_type", 
                        "tar_cell_type", 
                        "n_ref_cells",
                        "n_tar_cells", 
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions", 
                        "mixing_score", 
                        "normalised_mixing_score")
  
  # Turn numeric data into numeric type
  result[ , 3:8] <- apply(result[ , 3:8], 2, as.numeric)
  
  return(result)
}
