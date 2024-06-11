simulate_random_background_cells3D <- function(n_cells, 
                                               length, 
                                               width, 
                                               height, 
                                               minimum_distance_between_cells = 2, 
                                               oversampling_rate = 1.2, 
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
  
  
  # Check
  if(!is.numeric(minimum_distance_between_cells) | !is.numeric(oversampling_rate)){
    stop("One or more of `minimum_distance_between_cells`, `oversampling_rate` is not numeric!")
  }
  
  # Need to oversample as we will be removing cells too close later
  n_cells_inflated <- n_cells * oversampling_rate
  
  # Use poisson distribution to sample points
  pois_df <- poisson_distribution3D(n_cells = n_cells_inflated, 
                                    length = length, 
                                    width = width, 
                                    height = height)
  
  rownames(pois_df) <- paste("Cell_", seq(nrow(pois_df)), sep = "")    
  
  
  ### Check if all other cells are to close to the current cell 
  # Use frNN function: for each point, get all points within min_d of it
  pois_df_distances <- dbscan::frNN(pois_df, 
                                    eps = minimum_distance_between_cells,
                                    query = NULL, 
                                    sort = FALSE)
  
  # Check only the cells, don't care about the exact distance (definitely smaller than min_d)
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
  # Left over cells are the cells we keep
  chosen_cells <- names(pois_df_distances_ids)
  
  pois_df <- pois_df[chosen_cells, ]
  
  x <- pois_df$x
  y <- pois_df$y
  z <- pois_df$z
  
  
  # Put data into data frame
  df <- data.frame("Cell.X.Position" = x,
                   "Cell.Y.Position" = y,
                   "Cell.Z.Position" = z,
                   "Cell.Type" = background_cell_type)
  
  # Plot
  if (plot_image) {
    
    ## Plot
    fig <- plot_ly(df,
                   type = "scatter3d",
                   mode = 'markers',
                   x = ~Cell.X.Position,
                   y = ~Cell.Y.Position,
                   z = ~Cell.Z.Position,
                   color = ~Cell.Type,
                   colors = "lightgray",
                   marker = list(size = 2))
    
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                       yaxis = list(title = 'y'),
                                       zaxis = list(title = 'z')))
    
    print(fig)
  }
  
  return (df)
}
