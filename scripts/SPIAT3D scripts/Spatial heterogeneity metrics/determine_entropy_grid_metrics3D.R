determine_entropy_grid_metrics3D <- function(spe, 
                                             n_splits,
                                             cell_types_of_interest,
                                             feature_colname = "Cell.Type",
                                             plot_image = TRUE) {
  
  
  # Check if n_splits is numeric
  if (!is.numeric(n_splits)) {
    stop(paste(n_splits, " n_splits is not of type 'numeric'"))
  }
  
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, unique(spe[[feature_colname]]))
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  
  spe_coords <- data.frame(spatialCoords(spe))
  
  ## Get dimensions of the window
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  d_lay <- height / n_splits
  
  ## Figure out which 'grid prism number' each cell is inside
  spe$Prism.Num <- floor(spe_coords$Cell.X.Position / d_row) +
    floor(spe_coords$Cell.Y.Position / d_col) * n_splits + 
    floor(spe_coords$Cell.Z.Position / d_lay) * n_splits^2 + 1
  
  ## Get number of grid prisms
  n_grid_prisms <- n_splits^3
  
  ## Define data frame which contains all results
  result <- data.frame(matrix(nrow = n_grid_prisms, ncol = (length(cell_types_of_interest) + 2)))
  colnames(result) <- c(cell_types_of_interest, "total", "entropy")
  
  ## Calculate entropy for each grid prism
  for (grid_prism_num in seq(n_grid_prisms)) {
    
    ## Get spe object for current grid_prism
    spe_temp <- spe[ , spe$Prism.Num == grid_prism_num, ]
    
    ## Get cell_types_of_interest found in the sub-spe object
    temp_cell_types_of_interest <- intersect(cell_types_of_interest, unique(spe_temp[[feature_colname]]))
    
    grid_prism_entropy <- calculate_entropy_background3D(spe_temp,
                                                         temp_cell_types_of_interest)
    result[grid_prism_num, "entropy"] <- grid_prism_entropy
    
    ## Get number of cells of each cell_types_of_interest in each grid prism
    for (cell_type_of_interest in cell_types_of_interest) {
      result[grid_prism_num, cell_type_of_interest] <- sum(spe_temp[[feature_colname]] == cell_type_of_interest)
    }
  }
  
  ## Add column for total cell count for each grid prism
  result$total <- apply(result[ , colnames(result) %in% cell_types_of_interest], 1, sum)
  
  ## Add x, y and z coords of each grid prism to result
  result$prism_x_coord <- ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row
  result$prism_y_coord <- (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col
  result$prism_z_coord <- (floor((seq(n_grid_prisms) - 1) / (n_splits^2)) + 0.5) * d_lay
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous3D(result, "entropy")
    methods::show(fig)
  }
  
  return(result)
}
