simulate_normal_background_cells3D <- function(n_cells, 
                                               length, 
                                               width, 
                                               height,
                                               jitter_proportion = 0.25,
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
  if (!is.numeric(jitter_proportion)) {
    stop("`jitter` should be numeric!")
  }
  
  # Obtain distance between each point using MAGIC formula
  s <- ((sqrt(2) * length * width * height)/n_cells)^(1/3)
  
  # Get value for x_cells (points in 1 row),
  #               y_cells (points in 1 column) and 
  #               z_cells (points in 1 vertical thing), rounded
  x_cells <- round(length/s)
  y_cells <- round((2 * width)/(sqrt(3) * s))
  z_cells <- round((3 * height)/(sqrt(6) * s))
  
  # Phase 0. Assume points are on a 3D rectangular grid
  x <- rep(1:x_cells, y_cells * z_cells) * s
  y <- rep(rep(1:y_cells, each = x_cells), z_cells) * ((sqrt(3)*s)/2)
  z <- rep(1:z_cells, each = x_cells * y_cells) * ((sqrt(6)*s)/3)
  
  # Phase 1. For every odd sheet, every even row shifts by s/2 right
  if (y_cells %% 2 == 0) {
    shift <- rep(c(rep(0, x_cells), rep(s/2, x_cells)), y_cells/2)
  } else {
    shift <- c(rep(c(rep(0, x_cells), rep(s/2, x_cells)), y_cells/2), rep(0, x_cells))
  }
  x <- x + c(shift, rep(0, x_cells * y_cells)) # Shift each even row by s/2 right
  
  # Phase 2. For every even sheet, odd rows shift s/2 right, all rows shift s/(2*sqrt(3)) up
  if (y_cells %% 2 == 0) {
    shift <- rep(c(rep(s/2, x_cells), rep(0, x_cells)), y_cells/2)
  } else {
    shift <- c(rep(c(rep(s/2, x_cells), rep(0, x_cells)), y_cells/2), rep(s/2, x_cells))
  }
  x <- x + c(rep(0, x_cells * y_cells), shift) # Shift each odd row by s/2 right
  y <- y + rep(c(0, s/(2*sqrt(3))), each = x_cells*y_cells) # Shift all rows by s/(2*sqrt(3)) up
  
  # Get total number of cells (should be roughly equal to n_cells)
  n_total <- x_cells * y_cells * z_cells
  
  # Add randomness to the location of the cells
  jitter <- jitter_proportion * s # Jitter is proportional to distance between points in hexagonal grid
  jitter_x <- runif(n_total, -jitter, jitter)
  jitter_y <- runif(n_total, -jitter, jitter)
  jitter_z <- runif(n_total, -jitter, jitter)
  
  x <- x + jitter_x
  y <- y + jitter_y
  z <- z + jitter_z

  # Put data into data frame
  df <- data.frame("Cell.X.Position" = x,
                   "Cell.Y.Position" = y,
                   "Cell.Z.Position" = z,
                   "Cell.Type" = background_cell_type)
  
  # Plot
  if (plot_image) {
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

