#' Simulate background cells in 3D
#'
#' @description Simulate cell locations. The 3D locations of the cells are
#'   simulated and plotted using plot3D function in RGL. Users can specify the 
#'   window size, cell number and the minimum distance between two cells. 
#'   All cells have the same cell type, specified by the "cell_type" param.
#'
#' @details There are two options for the background cell distribution model. 
#'   Method = 1
#'   Use for tumour tissues. Assumes cells have a random location inside the
#'   window. Cells within the distance of `min_d` of another cell are deleted.
#'   This results fewer cells than specified by the user. Hence, we introduce
#'   paramter `oversampling_rate` to generate more cells than specified.
#'   
#'   Method = 2
#'   Normal tissues use an evenly-spaced model where the cells are distributed 
#'   approximately according to the vertices of a hexagon which are stacked in 
#'   the z-plane. The function accomplishes thisby generating cells on a 
#'   hexagonal grid in 3D and individually applying a random jitter which is 
#'   proportional to the original distance between the points. In our algorithm, 
#'   `jitter_prop` is the parameter to define the uniform distribution of the 
#'   jitter of the cells from the hexagon vertices. We recommend a `jitter_prop`
#'   near 0.25 to produce a sensible outcome.
#'   
#
#' @param n_cells Numeric. Number of cells to simulate in the background.
#' @param length,width,height Numeric. The length, width and height of the window.
#' @param method String. The distribution model for the background cells.
#'   Options are "tumour" for tumour tissues and "normal" for normal tissues 
#' @param min_d (OPTIONAL) Numeric. Use when `method` is "tumour". The minimum
#'   distance between two cells.
#' @param oversampling_rate (OPTIONAL) Numeric. Use when `method` is "tumour".
#'   The multiplier for oversampling. Without oversampling, the simulation
#'   deletes cells that are within `min_d` from each other, resulting in a less
#'   total number of cells than `n_cells`. Default is 1.2 (this should be set
#'   based on `n_cells`, `min_d`, `length`, `width` and `height`; should always 
#'   be larger than 1).
#' @param jitter_prop (OPTIONAL) Numeric. Use when `method` is "normal". The 
#'   uniform distribution parameter to generate the jitter distance for each 
#'   cell from the vertices of the hexagon. Default is 0.25.
#' @param cell_type (OPTIONAL) String. The name of the background cell type.
#'   Default is "Others" since there shouldn't be any identity of the background
#'   cells.
#' @param plot_image (OPTIONAL) Boolean. Default is TRUE.
#'
#' @family simulate pattern functions
#' @seealso \code{\link{simulate_mixing}} for mixed background simulation,
#'   \code{\link{simulate_clusters}} for cluster simulation,
#'   \code{\link{simulate_immune_rings}}/\code{\link{simulate_double_rings}} for
#'   immune ring simulation, and \code{\link{simulate_stripes}} for vessel
#'   simulation.
#'
#' @return A data.frame of the simulated background image
#' @export
#'
#' @examples
#' set.seed(610) # set seed for this background image simulation for reproducibility
#' background_image <- simulate_background_cells(n_cells = 10000, 
#'                                               length = 100
#'                                               width = 100,
#'                                               height = 100, 
#'                                               method = "tumour",
#'                                               min_d = 2,
#'                                               oversampling_rate = 1.2,
#'                                               jitter_prop = 0.25,
#'                                               cell_type = "Others",
#'                                               plot_image = TRUE)


simulate_background_cells3D <- function(n_cells, 
                                        length, 
                                        width, 
                                        height, 
                                        method, 
                                        min_d = 2, 
                                        oversampling_rate = 1.2, 
                                        jitter_prop = 0.25,
                                        cell_type = "Others", 
                                        plot_image = TRUE) {
  
  # Check
  if (!is.numeric(n_cells) | !is.numeric(length) | !is.numeric(width) | 
      !is.numeric(height)) {
    stop("One or more of `n_cells`, `length`, width`, `height` is not numeric!")
  }
  if (!is.character(cell_type)) {
    stop("`cell_type` should be of character type!")
  }
  
  
  ## Method = 1 (Tumour Tissue: Minimum distance constraint) ----------------------------------
  if (method == "tumour") {
    
    # Check
    if(!is.numeric(min_d) | !is.numeric(oversampling_rate)){
      stop("One or more of `min_d`, `oversampling_rate` is not numeric!")
    }
    
    # Need to oversample
    n_cells_inflated <- n_cells*oversampling_rate
    
    # Use poisson distribution to sample points
    pois_df <- poisson_distribution3D(n_cells = n_cells_inflated, 
                                      length = length, 
                                      width = width, 
                                      height = height)
    
    rownames(pois_df) <- paste("Cell_", seq(nrow(pois_df)), sep = "")    
    
    
    ### Check if all other cells are to close to the current cell 
    # Use frNN function: for each point, get all points within min_d of it
    pois_df_distances <- dbscan::frNN(pois_df, 
                                      eps = min_d,
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
    
  } 
  
  
  ## Method = 2 (Normal Tissue: Hexagonal grid pattern) ---------------------------------------
  else if (method == "normal") {
    
    # Check
    if (!is.numeric(jitter_prop)) {
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
    
    # First, assume points are on a 3D rectangular grid
    x <- rep(1:x_cells, y_cells * z_cells) * s
    y <- rep(rep(1:y_cells, each = x_cells), z_cells) * ((sqrt(3)*s)/2)
    z <- rep(1:z_cells, each = x_cells * y_cells) * ((sqrt(6)*s)/3)
    
    # Next:
    # Phase 1. For every odd sheet, every even row shifts by s/2 right
    # Phase 2. For every even sheet, odd rows shift s/2 right,
    #                    all rows shift s/(2*sqrt(3)) up
    
    # Phase 1. Every odd sheet
    if (y_cells %% 2 == 0) {
      shift <- rep(c(rep(0, x_cells), rep(s/2, x_cells)), y_cells/2)
    } else {
      shift <- c(rep(c(rep(0, x_cells), rep(s/2, x_cells)), y_cells/2), rep(0, x_cells))
    }
    
    x <- x + c(shift, rep(0, x_cells * y_cells)) # Shift each even row by s/2 right
    
    
    # Phase 2. Every even sheet
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
    jitter <- jitter_prop * s # Jitter is proportional to distance between points in hexagonal grid
    jitter_x <- runif(n_total, -jitter, jitter)
    jitter_y <- runif(n_total, -jitter, jitter)
    jitter_z <- runif(n_total, -jitter, jitter)
    
    x <- x + jitter_x
    y <- y + jitter_y
    z <- z + jitter_z
  } 
  
  else {
    stop("`method` should be 'tumour' or 'normal'")
  }
  
  
  
  df <- data.frame("Cell.X.Position" = x,
                   "Cell.Y.Position" = y,
                   "Cell.Z.Position" = z,
                   "Cell.Type" = cell_type)
  
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
