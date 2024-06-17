simulate_random_background_cells3D <- function(n_cells, 
                                               length, 
                                               width, 
                                               height, 
                                               minimum_distance_between_cells,
                                               background_cell_type = "Others", 
                                               plot_image = TRUE) {
  
  # Check
  if (!is.numeric(n_cells) | !is.numeric(length) | !is.numeric(width) | 
      !is.numeric(height)) {
    stop("One or more of `n_cells`, `length`, width`, `height` is not numeric!")
  }
  if (!is.character(background_cell_type)) {
    stop("`background_cell_type` should be of character type!")
  }
  if(!is.numeric(minimum_distance_between_cells)) {
    stop("`minimum_distance_between_cells` is not numeric!")
  }
  
  # Need to over-sample as cells which are too close will be removed later
  n_cells_inflated <- n_cells * 1.2
  
  # Use poisson distribution to sample points
  pois_df <- poisson_distribution3D(n_cells = n_cells_inflated, 
                                    length = length, 
                                    width = width, 
                                    height = height)
  
  # Give cells a unique ID
  rownames(pois_df) <- paste("Cell_", seq(nrow(pois_df)), sep = "")    
  
  ### Check if all other cells are to close to the current cell 
  # Use frNN function: for each point, get all points within min_d of it
  pois_df_distances <- dbscan::frNN(pois_df, 
                                    eps = minimum_distance_between_cells,
                                    query = NULL, 
                                    sort = FALSE)
  
  # For each cell, get all other cells which were within 'minimum_distance_between_cells'
  pois_df_distances_ids <- pois_df_distances$id
  
  n_cells <- nrow(pois_df)
  i <- 1
  
  while (i < n_cells) {
    cells_too_close <- paste("Cell_", pois_df_distances_ids[[i]], sep = "")
    
    for (cell in cells_too_close) {
      
      ## Remove cell that is too close
      if (!is.null(pois_df_distances_ids[[cell]])) {
        pois_df_distances_ids[cell] <- NULL
        n_cells <- n_cells - 1
      }
    }
    i <- i + 1
  }
  
  # Left over cells are the cells we choose
  cells_chosen <- names(pois_df_distances_ids)
  
  pois_df <- pois_df[cells_chosen, ]
  
  x <- pois_df$x
  y <- pois_df$y
  z <- pois_df$z
  
  # Put data into data frame
  df <- data.frame("Cell.X.Position" = x,
                   "Cell.Y.Position" = y,
                   "Cell.Z.Position" = z,
                   "Cell.Type" = background_cell_type)
  df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")
  
  # Get meta data
  background_metadata <- list("background_type" = "random",
                              "n_cells" = n_cells,
                              "length" = length,
                              "width" = width,
                              "height" = height,
                              "minimum_distance_between_cells" = minimum_distance_between_cells,
                              "cell_types" = background_cell_type,
                              "cell_proportions" = 1)
  
  ## Convert data frame to spe object
  spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = list(background = background_metadata))
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        background_cell_type,
                        "lightgray")
    print(fig)
  }
  
  return (spe)
}
