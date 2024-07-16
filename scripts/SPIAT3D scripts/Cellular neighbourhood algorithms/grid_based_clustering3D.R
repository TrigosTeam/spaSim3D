grid_based_clustering3D <- function(spe,
                                    cell_types_of_interest,
                                    n_splits,
                                    feature_colname = "Cell.Type",
                                    plot_image = TRUE) {
  
  # Check if n_splits is numeric
  if (!is.numeric(n_splits)) {
    stop(paste(n_splits, " n_splits is not of type 'numeric'"))
  }
  
  ## Check cell_types_of_interest are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
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

  ## Calculate proportions for each grid prism
  n_grid_prisms <- n_splits^3
  grid_prism_cell_proportions <- c()
  for (grid_prism_num in seq(n_grid_prisms)) {
    
    ## Get spe object for current grid_prism
    spe_temp <- spe[ , spe$Prism.Num == grid_prism_num]
    
    ## Get total number of cells and number of cell_types_of_interest in current grid_prism
    n_total <- ncol(spe_temp)
    n_cell_types_of_interest <- sum(spe_temp[[feature_colname]] %in% cell_types_of_interest)
    
    if (n_total == 0) {
      grid_prism_cell_proportion <- 0
    }
    else {
      grid_prism_cell_proportion <- n_cell_types_of_interest / n_total
    }
    
    grid_prism_cell_proportions <- c(grid_prism_cell_proportions, grid_prism_cell_proportion)
  }
  names(grid_prism_cell_proportions) <- seq(n_grid_prisms)
  
  
  ## Create template for final result
  result <- list()
  n_clusters <- 1
  
  
  ### CLUSTER DETECTION RECURSIVE ALGORITHM LOOP ###
  
  # First, remove all 0s from grid_prism_cell_proportions
  grid_prism_cell_proportions <- grid_prism_cell_proportions[grid_prism_cell_proportions != 0]
  
  while (length(grid_prism_cell_proportions) != 0) {
    # Get the maximum cell proportion and its corresponding grid prism number
    maximum_cell_proportion <- max(grid_prism_cell_proportions)
    maximum_cell_proportion_prism_number <- as.numeric(names(which.max(grid_prism_cell_proportions)))
    
    # Break out the loop if maximum cell proportion is less than 0.5
    if (maximum_cell_proportion < 0.5) break 
    
    # Else, find all the grid prisms adjacent to the maximum cell proportion grid prism. 
    # These are potentially apart of the cluster
    # Adjacent grid prisms must have cell proportion > 0.25 * max cell proportion
    grid_prisms_in_cluster <- determine_grid_prism_numbers_in_cluster3D(maximum_cell_proportion_prism_number,
                                                                        grid_prism_cell_proportions,
                                                                        maximum_cell_proportion,
                                                                        n_splits,
                                                                        c())
  
    # Perform the recursive algorithm on each grid prism potentially apart of the cluster to get a more precise shape of each cluster
    result[[n_clusters]] <- data.frame()
    for (grid_prism in as.numeric(grid_prisms_in_cluster)) {
      
      spe_prism <- spe[ , spe$Prism.Num == grid_prism]
      
      grid_prism_x <- ((grid_prism - 1) %% n_splits) * d_row
      grid_prism_y <- (floor(((grid_prism - 1) %% n_splits^2) / n_splits)) * d_col
      grid_prism_z <- (floor((grid_prism - 1) / n_splits^2)) * d_lay
      
      result[[n_clusters]] <- rbind(result[[n_clusters]], 
                                    grid_based_cluster_recursion3D(spe_prism, 
                                                                   cell_types_of_interest, 
                                                                   0.75 * maximum_cell_proportion,
                                                                   grid_prism_x, grid_prism_y, grid_prism_z,
                                                                   d_row, d_col, d_lay,
                                                                   "Cell.Type",
                                                                   data.frame()))
    }
    
    colnames(result[[n_clusters]]) <- c("x", "y", "z", "l", "w", "h")
    n_clusters <- n_clusters + 1
    
    # Remove grid prisms which have just been examined
    grid_prism_cell_proportions <- grid_prism_cell_proportions[setdiff((names(grid_prism_cell_proportions)), 
                                                                       grid_prisms_in_cluster)]
    
  }
  
  ## Add all the information to the spe
  spe@metadata[["grid_prisms"]] <- result
  spe$Prism.Num <- NULL
  spe$grid_based_cluster <- 0
  cluster_number <- 1
  
  for (cluster_info in result) {
    for (i in seq(nrow(cluster_info))) {
      x <- cluster_info$x[i]
      y <- cluster_info$y[i]
      z <- cluster_info$z[i]
      l <- cluster_info$l[i]
      w <- cluster_info$w[i]
      h <- cluster_info$h[i]
      
      spe$grid_based_cluster <- ifelse(spe_coords$Cell.X.Position >= x &
                                         spe_coords$Cell.X.Position < (x + l) &
                                         spe_coords$Cell.Y.Position >= y &
                                         spe_coords$Cell.Y.Position < (y + w) &
                                         spe_coords$Cell.Z.Position >= z &
                                         spe_coords$Cell.Z.Position < (z + h), 
                                       cluster_number, 
                                       spe$grid_based_cluster)
    }
    cluster_number <- cluster_number + 1
  }
  
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_based_clusters3D(spe, feature_colname = feature_colname)
    methods::show(fig)
  }
  
  return(spe)
}
