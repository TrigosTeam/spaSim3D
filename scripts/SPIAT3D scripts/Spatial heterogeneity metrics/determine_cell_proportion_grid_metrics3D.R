determine_cell_proportion_grid_metrics3D <- function(spe, 
                                                     n_splits,
                                                     reference_cell_types,
                                                     target_cell_types,
                                                     feature_colname = "Cell.Type",
                                                     plot_image = TRUE) {
  

  
  # Check if n_splits is numeric
  if (!is.numeric(n_splits)) {
    stop(paste(n_splits, " n_splits is not of type 'numeric'"))
  }
  
  ## Check reference_cell_types are found in the spe object
  unknown_cell_types <- setdiff(reference_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in reference_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  ## Check target_cell_types are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  # Check if there is intersection between reference_cell_types and target_cell_types
  if (length(intersect(reference_cell_types, target_cell_types)) > 0) {
    stop("Cannot have same cells in both reference_cell_types and target_cell_types")
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
  
  ## Calculate cell_proportions for each grid prism
  n_grid_prisms <- n_splits^3
  n_reference_cells_vec <- c()
  n_target_cells_vec <- c()
  grid_prism_cell_proportions <- c()
  
  ## Define data frame which contains all results
  result <- data.frame(matrix(nrow = n_grid_prisms, ncol = 4))
  colnames(result) <- c("reference", "target", "total", "proportion")
  
  for (grid_prism_num in seq(n_grid_prisms)) {
    
    ## Get spe object for current grid_prism
    spe_temp <- spe[ , spe$Prism.Num == grid_prism_num, ]
    
    ## Get cell_proportion: n_target_cells / (n_target_cells + n_reference_cells)
    n_target_cells <- sum(spe_temp[[feature_colname]] %in% target_cell_types)
    n_reference_cells <- sum(spe_temp[[feature_colname]] %in% reference_cell_types)
    
    ## Case when there are no target or reference cell, result is NA
    if (n_target_cells == 0 && n_reference_cells == 0) {
      grid_prism_cell_proportion <- NA  
    }
    else {
      grid_prism_cell_proportion <- n_target_cells / (n_target_cells + n_reference_cells)
    }
    
    result[grid_prism_num, ] <- c(n_reference_cells, 
                                  n_target_cells, 
                                  n_reference_cells + n_target_cells, 
                                  grid_prism_cell_proportion)
  }
  
  ## Add x, y and z coords of each grid prism to result
  result$prism_x_coord <- ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row
  result$prism_y_coord <- (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col
  result$prism_z_coord <- (floor((seq(n_grid_prisms) - 1) / (n_splits^2)) + 0.5) * d_lay
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous3D(result, "proportion")
    methods::show(fig)
  }
  
  return(result)
}
