grid_based_clustering3D <- function(spe,
                                    cell_types_of_interest,
                                    n_splits,
                                    feature_colname = "Cell.Type",
                                    plot_image = TRUE) {
  
  if (is.null(spe[[feature_colname]])) stop(paste("No column called", feature_colname, "found in spe object"))
  
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
  
  # Add grid metrics to spe
  spe <- get_spe_grid_metrics3D(spe, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spe
  grid_prism_cell_matrix <- spe@metadata$grid_metrics$grid_prism_cell_matrix
  
  ## Calculate proportions for each grid prism
  if (length(cell_types_of_interest) == 1) {
    grid_prism_cell_proportions <- grid_prism_cell_matrix[ , cell_types_of_interest]
  }
  else {
    grid_prism_cell_proportions <- rowSums(grid_prism_cell_matrix[ , cell_types_of_interest])
  }
  grid_prism_cell_proportions <- grid_prism_cell_proportions / rowSums(grid_prism_cell_matrix[ , unique(spe[[feature_colname]])])
  n_grid_prisms <- n_splits^3
  names(grid_prism_cell_proportions) <- seq(n_grid_prisms)
  
  
  ## Create template for final result
  result <- list()
  n_clusters <- 1
  
  ## Get dimensions of the window
  spe_coords <- data.frame(spatialCoords(spe))
  
  min_x <- min(spe_coords$Cell.X.Position)
  min_y <- min(spe_coords$Cell.Y.Position)
  min_z <- min(spe_coords$Cell.Z.Position)
  
  max_x <- max(spe_coords$Cell.X.Position)
  max_y <- max(spe_coords$Cell.Y.Position)
  max_z <- max(spe_coords$Cell.Z.Position)
  
  length <- round(max_x - min_x)
  width  <- round(max_y - min_y)
  height <- round(max_z - min_z)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  d_lay <- height / n_splits
  
  
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
    grid_prisms_in_cluster <- calculate_grid_prism_numbers_in_cluster3D(maximum_cell_proportion_prism_number,
                                                                        grid_prism_cell_proportions,
                                                                        maximum_cell_proportion,
                                                                        n_splits,
                                                                        c())
    
    # Perform the recursive algorithm on each grid prism potentially apart of the cluster to get a more precise shape of each cluster
    # Create data frame with spatial coords and cell types as columns. Use this as input
    df <- spe_coords
    df[[feature_colname]] <- spe[[feature_colname]] 
    curr_result <- sapply(as.numeric(grid_prisms_in_cluster), function(x) grid_based_cluster_recursion3D(df,
                                                                                                         cell_types_of_interest,
                                                                                                         0.75 * maximum_cell_proportion,
                                                                                                         ((x - 1) %% n_splits) * d_row + round(min_x),
                                                                                                         (floor(((x - 1) %% n_splits^2) / n_splits)) * d_col + round(min_y),
                                                                                                         (floor((x - 1) / n_splits^2)) * d_lay + round(min_z),
                                                                                                         d_row, d_col, d_lay,
                                                                                                         feature_colname,
                                                                                                         data.frame()))
    
    
    
    if (is.array(curr_result))  {
      curr_result <- data.frame(t(unlist(curr_result)))
      colnames(curr_result) <- c("x", "y", "z", "l", "w", "h")
      result[[n_clusters]] <- curr_result
    }
    else {
      result[[n_clusters]] <- rbindlist(curr_result)
    }
    
    n_clusters <- n_clusters + 1
    
    # Remove grid prisms which have just been examined
    grid_prism_cell_proportions <- grid_prism_cell_proportions[setdiff((names(grid_prism_cell_proportions)), 
                                                                       grid_prisms_in_cluster)]
    
  }
  
  ## Add all the information to the spe
  spe@metadata[["grid_prisms"]] <- result
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
                                         spe_coords$Cell.Z.Position < (z + h) &
                                         spe[[feature_colname]] %in% cell_types_of_interest, 
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
