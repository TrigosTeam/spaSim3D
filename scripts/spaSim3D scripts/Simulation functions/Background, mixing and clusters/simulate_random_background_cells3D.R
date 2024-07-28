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
  
  # Add integer rownames to data frame - each cell is labelled by an integer
  rownames(pois_df) <- seq(nrow(pois_df)) 
  
  ### Check if all other cells are to close to the current cell 
  # Use frNN function: for each point, get all points within min_d of it
  pois_df_distances <- dbscan::frNN(pois_df, 
                                    eps = minimum_distance_between_cells,
                                    query = NULL, 
                                    sort = FALSE)
  
  # For each cell, get all other cells which were within 'minimum_distance_between_cells'
  pois_df_distances_ids <- pois_df_distances$id
  
  # Filter out zero length cells
  pois_df_distances_ids <- Filter(function(x) length(x) != 0, pois_df_distances_ids)
  
  # Get integer labels for the remaining cells
  pois_df_distances_ids_names <- as.integer(names(pois_df_distances_ids))
  
  # Determine which cells should be chosen from pois_df
  cells_chosen <- rep(T, nrow(pois_df))
  for (i in seq_len(length(pois_df_distances_ids))) {
    cells_too_close <- pois_df_distances_ids[[i]]
    
    if (cells_chosen[pois_df_distances_ids_names[i]]) cells_chosen[cells_too_close] <- F
  }
  
  pois_df <- pois_df[cells_chosen, ]
  
  # If number of cells remaining is still higher than n_cells, randomly subset n_cells cells
  if (nrow(pois_df) > n_cells) pois_df <- dplyr::sample_n(pois_df, n_cells)

  # Add Cell.Type and Cell.ID
  pois_df$Cell.Type <- background_cell_type
  pois_df$Cell.ID <- paste("Cell", seq(nrow(pois_df)), sep = "_")
  
  # Get meta data
  background_metadata <- list("background_type" = "random",
                              "n_cells" = n_cells,
                              "length" = length,
                              "width" = width,
                              "height" = height,
                              "minimum_distance_between_cells" = minimum_distance_between_cells,
                              "cell_types" = background_cell_type,
                              "cell_proportions" = 1)
  simulation_metadata <- list(background = background_metadata)
  
  ## Convert data frame to spe object
  spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(pois_df), ncol = nrow(pois_df)),
    colData = pois_df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = list(simulation = simulation_metadata))
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        background_cell_type,
                        "lightgray")
    methods::show(fig)
  }
  
  return(spe)
}

