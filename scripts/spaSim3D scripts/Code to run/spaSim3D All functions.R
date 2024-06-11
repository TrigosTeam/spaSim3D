library(plotly)


poisson_distribution3D <- function(n_cells, length, width, height)  {
  
  # Choose lambda
  lambda <- 5
  
  # Set number of rows, columns and layers
  nRows <- nCols <- nLays <- round((n_cells/lambda)^(1/3))
  
  # Get number of cubes in grid
  nCubes <- nRows * nCols * nLays
  
  # Get pois vector
  pois <- rpois(nCubes, lambda)
  
  
  # Get points for each square region
  x <- c()
  y <- c()
  z <- c()
  
  for (row in seq(nRows)) {
    
    for (col in seq(nCols)) {
      
      for (lay in seq(nLays)) {
        current_cube_index <- nRows^2 * (row - 1) + nCols * (col - 1) + lay
        
        x <- append(x, runif(pois[current_cube_index], row - 1, row))
        y <- append(y, runif(pois[current_cube_index], col - 1, col))
        z <- append(z, runif(pois[current_cube_index], lay - 1, lay))
        
      }
    }
  }
  x <- x * length / nRows
  y <- y * width / nCols
  z <- z * height / nLays
  
  df <- data.frame(x = x,
                   y = y,
                   z = z)
  return(df)
}


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
  
  # Check
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








simulate_clusters3D <- function(bg_sample,
                                n_clusters = 3,
                                cluster_properties = list(
                                  C1 = list(
                                    name_of_cluster_cell = "Tumour",
                                    infiltration_types = c("Immune", "Others"),
                                    infiltration_proportions = c(0.4, 0.05),
                                    shape = "Sphere",
                                    radius = 25,
                                    centre_loc = c(40, 40, 40)
                                  ),
                                  C2 = list(
                                    name_of_cluster_cell = "Endothelial",
                                    infiltration_types = c("Others"),
                                    infiltration_proportions = c(0.05),
                                    shape = "Cylinder",
                                    radius = 10,
                                    start_loc = c(0, 0, 0),
                                    end_loc   = c(20, 20 , 100)
                                  ),
                                  C3 = list(
                                    name_of_cluster_cell = "Tumour",
                                    infiltration_types = c("Immune", "Others"),
                                    infiltration_proportions = c(0.3, 0.05),
                                    shape = "Ellipsoid",
                                    x_radius = 15,
                                    y_radius = 20,
                                    z_radius = 25,
                                    centre_loc = c(70, 70, 70),
                                    x_y_rotation = 0,
                                    x_z_rotation = 0,
                                    y_z_rotation = 0
                                  )
                                ),
                                plot_image = TRUE,
                                plot_categories = c("Others", "Immune", "Endothelial", "Tumour"),
                                plot_colours = c("lightgray", "skyblue", "#FF7F7F", "orange")) {
  
  
  for (k in seq_len(n_clusters)) { 
    
    # for each cluster, get the shape
    shape <- cluster_properties[[k]]$shape
    
    # # generate a location as the centre of the cluster
    # if (is.null(centre_loc)){
    #   seed_point <- spatstat.random::runifpoint(1, win=win)}
    # else seed_point <- centre_loc
    # a <- seed_point$x
    # b <- seed_point$y
    
    ### Sphere shape
    if (shape == "Sphere") {
      bg_sample <- simulate_sphere_cluster(bg_sample = bg_sample, cluster_properties = cluster_properties[[k]])
    } 
    
    ### Ellipsoid shape
    if (shape == "Ellipsoid") {
      bg_sample <- simulate_ellipsoid_cluster(bg_sample = bg_sample, cluster_properties = cluster_properties[[k]])
    }
    
    ### Cylinder shape
    if (shape == "Cylinder") {
      bg_sample <- simulate_cylinder_cluster(bg_sample = bg_sample, cluster_properties = cluster_properties[[k]])
    }
    
    ### Network shape
    if (shape == "Network") {
      bg_sample <- simulate_network_cluster(bg_sample = bg_sample, cluster_properties = cluster_properties[[k]])
    }
    
  }
  
  if (plot_image) {
    
    if (is.null(plot_categories)) {
      plot_categories <- unique(bg_sample$Cell.Type)
    }
    
    if (is.null(plot_colours)) {
      plot_colours <- c("red", "orange", "green", "blue", "skyblue", "pink", "purple", "lightgray")[1:length(plot_categories)]
    }
    
    bg_sample <- bg_sample[bg_sample[["Cell.Type"]] %in% plot_categories, ]
    
    ## Factor for feature column
    bg_sample_factored <- bg_sample
    
    bg_sample_factored[, "Cell.Type"] <- factor(bg_sample_factored[, "Cell.Type"],
                                                levels = plot_categories)
    
    ## Plot
    fig <- plot_ly(bg_sample_factored,
                   type = "scatter3d",
                   mode = 'markers',
                   x = ~Cell.X.Position,
                   y = ~Cell.Y.Position,
                   z = ~Cell.Z.Position,
                   color = ~Cell.Type,
                   colors = plot_colours,
                   marker = list(size = 2))
    
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                       yaxis = list(title = 'y'),
                                       zaxis = list(title = 'z')))
    print(fig)
  }
  
  return (bg_sample)
}





simulate_cylinder_cluster <- function(bg_sample, cluster_properties) {
  
  # Get cylinder properties
  cell_type <- cluster_properties$name_of_cluster_cell
  infiltration_types <- cluster_properties$infiltration_types
  infiltration_proportions <- cluster_properties$infiltration_proportions
  radius <- cluster_properties$radius
  start_loc <- cluster_properties$start_loc
  end_loc <- cluster_properties$end_loc
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  for (i in seq_len(n_cells)) {
    # Get x, y, z and phenotype of ith cell
    x <- bg_sample[i, "Cell.X.Position"]
    y <- bg_sample[i, "Cell.Y.Position"]
    z <- bg_sample[i, "Cell.Z.Position"]
    pheno <- bg_sample[i, "Cell.Type"]
    
    # Ignore points outside of these planes
    if (sum(v1 *  c(x, y, z)) < d1 || sum(v1 * c(x, y, z)) > d2) {
      next
    }
    
    # Get vector between from point to start_loc
    v2 <- c(x, y, z) - start_loc
    
    # Get perpendicular distance squared between point and line
    D <- ((v1[2]*v2[3] - v1[3]*v2[2])^2 + 
            (v1[1]*v2[3] - v1[3]*v2[1])^2 + 
            (v1[1]*v2[2] - v1[2]*v2[1])^2)/
      (v1[1]^2 + v1[2]^2 + v1[3]^2)
    
    # Get maximum distance squared
    R <- radius^2
    
    if (D < R){ 
      # in the region of cluster, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_infiltration_types <- length(infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this cluster
      pheno <- cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_infiltration_types){
        current_p <- current_p + infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- infiltration_types[n]
          break
        }
        n <- n+1
      }
    }
    bg_sample[i, "Cell.Type"] <- pheno
  }
  return(bg_sample)
}




simulate_cylinder_ring <- function(bg_sample, ring_properties) {
  
  # Get sphere immune ring properties
  cell_type <- ring_properties$name_of_cluster_cell
  infiltration_types <- ring_properties$infiltration_types
  infiltration_proportions <- ring_properties$infiltration_proportions
  radius <- ring_properties$radius
  start_loc <- ring_properties$start_loc
  end_loc <- ring_properties$end_loc
  
  ring_cell_type <- ring_properties$name_of_ring_cell
  ring_width <- ring_properties$ring_width
  ring_infiltration_types <- ring_properties$ring_infiltration_types
  ring_infiltration_proportions <- ring_properties$ring_infiltration_proportions
  
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  i <- 1
  
  while (i <= n_cells) {
    # Get x, y, z and phenotype of ith cell
    x <- bg_sample[i, "Cell.X.Position"]
    y <- bg_sample[i, "Cell.Y.Position"]
    z <- bg_sample[i, "Cell.Z.Position"]
    pheno <- bg_sample[i, "Cell.Type"]
    
    # Ignore points outside of these planes
    if (sum(v1 *  c(x, y, z)) < d1 || sum(v1 * c(x, y, z)) > d2) {
      i <- i + 1
      next
    }
    
    # Get vector between from point to start_loc
    v2 <- c(x, y, z) - start_loc
    
    # Get perpendicular distance squared between point and line
    D <- ((v1[2]*v2[3] - v1[3]*v2[2])^2 + 
            (v1[1]*v2[3] - v1[3]*v2[1])^2 + 
            (v1[1]*v2[2] - v1[2]*v2[1])^2)/
      (v1[1]^2 + v1[2]^2 + v1[3]^2)
    
    # Get maximum distance without and with ring squared
    R1 <- radius^2
    R2 <- (radius + ring_width)^2
    
    if (D < R1) { 
      # in the region of cluster, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_infiltration_types <- length(infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this cluster
      pheno <- cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_infiltration_types){
        current_p <- current_p + infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- infiltration_types[n]
          break
        }
        n <- n + 1
      }
    }
    else if (D < R2) {
      # in the region of ring, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_ring_infiltration_types <- length(ring_infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this ring
      pheno <- ring_cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_ring_infiltration_types){
        current_p <- current_p + ring_infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- ring_infiltration_types[n]
          break
        }
        n <- n + 1
      }
    }
    
    if (pheno == "Void") {
      bg_sample <- bg_sample[-c(i), ]
      n_cells <- n_cells - 1
      
    } else {
      bg_sample[i, "Cell.Type"] <- pheno  
      i <- i + 1
    }
    
  }
  return(bg_sample)
}



simulate_double_rings3D <- function(bg_sample,
                                    n_dr = 1,
                                    dr_properties = list(
                                      D1 = list(
                                        name_of_cluster_cell = "Tumour",
                                        infiltration_types = c("Immune1", "Others"),
                                        infiltration_proportions = c(0.1, 0.05),
                                        shape = "Sphere",
                                        radius = 35,
                                        centre_loc = c(50, 50, 50),
                                        name_of_inner_ring_cell = "Immune1",
                                        inner_ring_width = 6,
                                        inner_ring_infiltration_types = c("Others"),
                                        inner_ring_infiltration_proportions = c(0.15),
                                        name_of_outer_ring_cell = "Immune2",
                                        outer_ring_width = 3,
                                        outer_ring_infiltration_types = c("Others"),
                                        outer_ring_infiltration_proportions = c(0.15)
                                      )
                                    ),
                                    plot_image = TRUE,
                                    plot_categories = c("Others", "Tumour", "Immune1", "Immune2"),
                                    plot_colours = c("lightgray", "orange", "blue", "green")) {
  
  for (k in seq_len(n_dr)) { 
    
    # for each cluster, get the shape
    shape <- dr_properties[[k]]$shape
    
    
    ### Sphere shape
    if (shape == "Sphere") {
      bg_sample <- simulate_sphere_dr(bg_sample = bg_sample, dr_properties = dr_properties[[k]])
    } 
    
    ### Ellipsoid shape
    if (shape == "Ellipsoid") {
      bg_sample <- simulate_ellipsoid_dr(bg_sample = bg_sample, dr_properties = dr_properties[[k]])
    }
  }
  
  if (plot_image) {
    
    if (is.null(plot_categories)) {
      plot_categories <- unique(bg_sample$Cell.Type)
    }
    
    if (is.null(plot_colours)) {
      plot_colours <- c("red", "orange", "green", "blue", "skyblue", "pink", "purple", "lightgray")[1:length(plot_categories)]
    }
    
    bg_sample <- bg_sample[bg_sample[["Cell.Type"]] %in% plot_categories, ]
    
    ## Factor for feature column
    bg_sample_factored <- bg_sample
    
    bg_sample_factored[, "Cell.Type"] <- factor(bg_sample_factored[, "Cell.Type"],
                                                levels = plot_categories)
    
    ## Plot
    fig <- plot_ly(bg_sample_factored,
                   type = "scatter3d",
                   mode = 'markers',
                   x = ~Cell.X.Position,
                   y = ~Cell.Y.Position,
                   z = ~Cell.Z.Position,
                   color = ~Cell.Type,
                   colors = plot_colours,
                   marker = list(size = 2))
    
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                       yaxis = list(title = 'y'),
                                       zaxis = list(title = 'z')))
    print(fig)
  }
  
  return(bg_sample)
}



simulate_ellipsoid_cluster <- function(bg_sample, cluster_properties) {
  
  # Get ellipsoid properties
  cell_type <- cluster_properties$name_of_cluster_cell
  infiltration_types <- cluster_properties$infiltration_types
  infiltration_proportions <- cluster_properties$infiltration_proportions
  x_radius <- cluster_properties$x_radius
  y_radius <- cluster_properties$y_radius
  z_radius <- cluster_properties$z_radius
  centre_loc <- cluster_properties$centre_loc
  
  theta <- cluster_properties$y_z_rotation # in x-axis
  alpha <- cluster_properties$x_z_rotation # in y-axis
  beta  <- cluster_properties$x_y_rotation # in z-axis
  
  a <- cos(alpha) * cos(beta)
  b <- cos(alpha) * sin(beta)
  c <- sin(alpha)
  d <- -sin(theta) * sin(alpha) * cos(beta) - cos(theta) * sin(beta)
  e <- -sin(theta) * sin(alpha) * sin(beta) + cos(theta) * cos(beta)
  f <- sin(theta) * cos(alpha)
  g <- -cos(theta) * sin(alpha) * cos(beta) + sin(theta) * sin(beta)
  h <- -cos(theta) * sin(alpha) * sin(beta) - sin(theta) * cos(beta)
  i <- cos(theta) * cos(alpha)
  
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  for (index in seq_len(n_cells)) {
    # Get x, y, z and phenotype of ith cell
    x <- bg_sample[index, "Cell.X.Position"] - centre_loc[1]
    y <- bg_sample[index, "Cell.Y.Position"] - centre_loc[2]
    z <- bg_sample[index, "Cell.Z.Position"] - centre_loc[3]
    pheno <- bg_sample[index, "Cell.Type"]
    
    x_new <- a * x + b * y + c * z
    y_new <- d * x + e * y + f * z
    z_new <- g * x + h * y + i * z
    
    D <- (x_new/x_radius)^2 + 
      (y_new/y_radius)^2 + 
      (z_new/z_radius)^2
    
    if (D <= 1) { 
      # in the region of cluster, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_infiltration_types <- length(infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this cluster
      pheno <- cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_infiltration_types) {
        current_p <- current_p + infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- infiltration_types[n]
          break
        }
        n <- n+1
      }
    }
    bg_sample[index, "Cell.Type"] <- pheno
  }
  return(bg_sample)
}



simulate_ellipsoid_dr <- function(bg_sample, dr_properties) {
  
  # Get ellipsoid properties
  cell_type <- dr_properties$name_of_cluster_cell
  infiltration_types <- dr_properties$infiltration_types
  infiltration_proportions <- dr_properties$infiltration_proportions
  x_radius <- dr_properties$x_radius
  y_radius <- dr_properties$y_radius
  z_radius <- dr_properties$z_radius
  centre_loc <- dr_properties$centre_loc
  
  inner_ring_cell_type <- dr_properties$name_of_inner_ring_cell
  inner_ring_width <- dr_properties$inner_ring_width
  inner_ring_infiltration_types <- dr_properties$inner_ring_infiltration_types
  inner_ring_infiltration_proportions <- dr_properties$inner_ring_infiltration_proportions
  
  outer_ring_cell_type <- dr_properties$name_of_outer_ring_cell
  outer_ring_width <- dr_properties$outer_ring_width
  outer_ring_infiltration_types <- dr_properties$outer_ring_infiltration_types
  outer_ring_infiltration_proportions <- dr_properties$outer_ring_infiltration_proportions
  
  theta <- dr_properties$y_z_rotation # in x-axis
  alpha <- dr_properties$x_z_rotation # in y-axis
  beta  <- dr_properties$x_y_rotation # in z-axis
  
  a <- cos(alpha) * cos(beta)
  b <- cos(alpha) * sin(beta)
  c <- sin(alpha)
  d <- -sin(theta) * sin(alpha) * cos(beta) - cos(theta) * sin(beta)
  e <- -sin(theta) * sin(alpha) * sin(beta) + cos(theta) * cos(beta)
  f <- sin(theta) * cos(alpha)
  g <- -cos(theta) * sin(alpha) * cos(beta) + sin(theta) * sin(beta)
  h <- -cos(theta) * sin(alpha) * sin(beta) - sin(theta) * cos(beta)
  i <- cos(theta) * cos(alpha)
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  for (index in seq_len(n_cells)) {
    # Get x, y, z and phenotype of ith cell
    x <- bg_sample[index, "Cell.X.Position"] - centre_loc[1]
    y <- bg_sample[index, "Cell.Y.Position"] - centre_loc[2]
    z <- bg_sample[index, "Cell.Z.Position"] - centre_loc[3]
    pheno <- bg_sample[index, "Cell.Type"]
    
    x_new <- a * x + b * y + c * z
    y_new <- d * x + e * y + f * z
    z_new <- g * x + h * y + i * z
    
    D1 <- (x_new/x_radius)^2 + 
      (y_new/y_radius)^2 + 
      (z_new/z_radius)^2
    
    D2 <- (x_new/(x_radius + inner_ring_width))^2 + 
      (y_new/(y_radius + inner_ring_width))^2 + 
      (z_new/(z_radius + inner_ring_width))^2
    
    D3 <- (x_new/(x_radius + inner_ring_width + outer_ring_width))^2 + 
      (y_new/(y_radius + inner_ring_width + outer_ring_width))^2 + 
      (z_new/(z_radius + inner_ring_width + outer_ring_width))^2
    
    if (D1 <= 1) { 
      # in the region of cluster, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_infiltration_types <- length(infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this cluster
      pheno <- cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_infiltration_types) {
        current_p <- current_p + infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- infiltration_types[n]
          break
        }
        n <- n+1
      }
    }
    else if (D2 <= 1) {
      # in the region of ring, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_inner_ring_infiltration_types <- length(inner_ring_infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this ring
      pheno <- inner_ring_cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_inner_ring_infiltration_types) {
        current_p <- current_p + inner_ring_infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- inner_ring_infiltration_types[n]
          break
        }
        n <- n+1
      }
    }
    else if (D3 <= 1) {
      # in the region of ring, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_outer_ring_infiltration_types <- length(outer_ring_infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this ring
      pheno <- outer_ring_cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_outer_ring_infiltration_types) {
        current_p <- current_p + outer_ring_infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- outer_ring_infiltration_types[n]
          break
        }
        n <- n+1
      }
    }
    bg_sample[index, "Cell.Type"] <- pheno
    
  }
  return(bg_sample)
}



simulate_ellipsoid_ring <- function(bg_sample, ring_properties) {
  
  # Get ellipsoid properties
  cell_type <- ring_properties$name_of_cluster_cell
  infiltration_types <- ring_properties$infiltration_types
  infiltration_proportions <- ring_properties$infiltration_proportions
  x_radius <- ring_properties$x_radius
  y_radius <- ring_properties$y_radius
  z_radius <- ring_properties$z_radius
  centre_loc <- ring_properties$centre_loc
  
  ring_cell_type <- ring_properties$name_of_ring_cell
  ring_width <- ring_properties$ring_width
  ring_infiltration_types <- ring_properties$ring_infiltration_types
  ring_infiltration_proportions <- ring_properties$ring_infiltration_proportions
  
  theta <- ring_properties$y_z_rotation # in x-axis
  alpha <- ring_properties$x_z_rotation # in y-axis
  beta  <- ring_properties$x_y_rotation # in z-axis
  
  a <- cos(alpha) * cos(beta)
  b <- cos(alpha) * sin(beta)
  c <- sin(alpha)
  d <- -sin(theta) * sin(alpha) * cos(beta) - cos(theta) * sin(beta)
  e <- -sin(theta) * sin(alpha) * sin(beta) + cos(theta) * cos(beta)
  f <- sin(theta) * cos(alpha)
  g <- -cos(theta) * sin(alpha) * cos(beta) + sin(theta) * sin(beta)
  h <- -cos(theta) * sin(alpha) * sin(beta) - sin(theta) * cos(beta)
  i <- cos(theta) * cos(alpha)
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  for (index in seq_len(n_cells)) {
    # Get x, y, z and phenotype of ith cell
    x <- bg_sample[index, "Cell.X.Position"] - centre_loc[1]
    y <- bg_sample[index, "Cell.Y.Position"] - centre_loc[2]
    z <- bg_sample[index, "Cell.Z.Position"] - centre_loc[3]
    pheno <- bg_sample[index, "Cell.Type"]
    
    x_new <- a * x + b * y + c * z
    y_new <- d * x + e * y + f * z
    z_new <- g * x + h * y + i * z
    
    D1 <- (x_new/x_radius)^2 + 
      (y_new/y_radius)^2 + 
      (z_new/z_radius)^2
    
    D2 <- (x_new/(x_radius + ring_width))^2 + 
      (y_new/(y_radius + ring_width))^2 + 
      (z_new/(z_radius + ring_width))^2
    
    if (D1 <= 1) { 
      # in the region of cluster, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_infiltration_types <- length(infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this cluster
      pheno <- cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_infiltration_types) {
        current_p <- current_p + infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- infiltration_types[n]
          break
        }
        n <- n+1
      }
    }
    else if (D2 <= 1) {
      # in the region of ring, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_ring_infiltration_types <- length(ring_infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this ring
      pheno <- ring_cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_ring_infiltration_types) {
        current_p <- current_p + ring_infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- ring_infiltration_types[n]
          break
        }
        n <- n+1
      }
    }
    bg_sample[index, "Cell.Type"] <- pheno
    
  }
  return(bg_sample)
}



simulate_mixing3D <- function(bg_sample,
                              cell_types = c("Others", "Immune", "Tumour"),
                              props = c(0.5, 0.2, 0.3),
                              plot_image = TRUE,
                              plot_categories = c("Others", "Immune", "Tumour"),
                              plot_colours = c("lightgray", "skyblue", "orange")) {
  
  
  n_types <- length(cell_types)
  
  for (i in 1:nrow(bg_sample)) {
    x <- bg_sample$Cell.X.Position[i]
    y <- bg_sample$Cell.Y.Position[i]
    z <- bg_sample$Cell.Z.Position[i]
    
    random <- runif(1)
    
    # if the random number falls in the range of an infiltration proportion,
    # pheno will be the corresponding infiltraiton type
    n <- 1 # start from the first proportion
    current_p <- 0
    while (n <= n_types){
      current_p <- current_p + props[n]
      if (random <= current_p) {
        pheno <- cell_types[n]
        break
      }
      n <- n+1
    }
    bg_sample[i, "Cell.Type"] <- pheno
  }
  
  if (plot_image) {
    
    if (is.null(plot_categories)) {
      plot_categories <- unique(bg_sample$Cell.Type)
    }
    
    if (is.null(plot_colours)) {
      plot_colours <- c("red", "orange", "green", "blue", "skyblue", "pink", "purple", "lightgray")[1:length(plot_categories)]
    }
    
    bg_sample <- bg_sample[bg_sample[["Cell.Type"]] %in% plot_categories, ]
    
    ## Factor for feature column
    bg_sample_factored <- bg_sample
    
    bg_sample_factored[, "Cell.Type"] <- factor(bg_sample_factored[, "Cell.Type"],
                                                levels = plot_categories)
    
    ## Plot
    fig <- plot_ly(bg_sample_factored,
                   type = "scatter3d",
                   mode = 'markers',
                   x = ~Cell.X.Position,
                   y = ~Cell.Y.Position,
                   z = ~Cell.Z.Position,
                   color = ~Cell.Type,
                   colors = plot_colours,
                   marker = list(size = 2))
    
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                       yaxis = list(title = 'y'),
                                       zaxis = list(title = 'z')))
    print(fig)
  }
  
  return(bg_sample)
}



simulate_rings3D <- function(bg_sample,
                             n_ring = 3,
                             ring_properties = list(
                               R1 = list(
                                 name_of_cluster_cell = "Tumour",
                                 infiltration_types = c("Others"),
                                 infiltration_proportions = c(0.05),
                                 shape = "Sphere",
                                 radius = 20,
                                 centre_loc = c(40, 40, 40),
                                 name_of_ring_cell = "Immune",
                                 ring_width = 5,
                                 ring_infiltration_types = c("Others"),
                                 ring_infiltration_proportions = c(0.15)
                               ),
                               R2 = list(
                                 name_of_cluster_cell = "Void",
                                 infiltration_types = NULL,
                                 infiltration_proportions = NULL,
                                 shape = "Cylinder",
                                 radius = 8,
                                 start_loc = c(0, 0, 0),
                                 end_loc   = c(20, 20 , 100),
                                 name_of_ring_cell = "Endothelial",
                                 ring_width = 5,
                                 ring_infiltration_types = c("Others"),
                                 ring_infiltration_proportions = c(0.15)
                               ),
                               R3 = list(
                                 name_of_cluster_cell = "Tumour",
                                 infiltration_types = c("Others"),
                                 infiltration_proportions = c(0.05),
                                 shape = "Ellipsoid",
                                 x_radius = 10,
                                 y_radius = 15,
                                 z_radius = 20,
                                 centre_loc = c(70, 70, 70),
                                 x_y_rotation = 0,
                                 x_z_rotation = 0,
                                 y_z_rotation = 0,
                                 name_of_ring_cell = "Immune",
                                 ring_width = 5,
                                 ring_infiltration_types = c("Others"),
                                 ring_infiltration_proportions = c(0.15)
                               )
                             ),
                             plot_image = TRUE,
                             plot_categories = c("Others", "Tumour", "Immune", "Endothelial"),
                             plot_colours = c("lightgray", "orange", "skyblue", "#FF7F7F")) {
  
  for (k in seq_len(n_ring)) { 
    
    # for each cluster, get the shape
    shape <- ring_properties[[k]]$shape
    
    
    ### Sphere shape + immune ring
    if (shape == "Sphere") {
      bg_sample <- simulate_sphere_ring(bg_sample = bg_sample, ring_properties = ring_properties[[k]])
    } 
    
    ### Ellipsoid shape + ring
    if (shape == "Ellipsoid") {
      bg_sample <- simulate_ellipsoid_ring(bg_sample = bg_sample, ring_properties = ring_properties[[k]])
    }
    
    ### Cylinder shape + ring
    if (shape == "Cylinder") {
      bg_sample <- simulate_cylinder_ring(bg_sample = bg_sample, ring_properties = ring_properties[[k]])
    }
    
    ### Network shape + ring
    if (shape == "Network") {
      bg_sample <- simulate_network_ring(bg_sample = bg_sample, ring_properties = ring_properties[[k]])
    }
  }
  
  if (plot_image) {
    
    if (is.null(plot_categories)) {
      plot_categories <- unique(bg_sample$Cell.Type)
    }
    
    if (is.null(plot_colours)) {
      plot_colours <- c("red", "orange", "green", "blue", "skyblue", "pink", "purple", "lightgray")[1:length(plot_categories)]
    }
    
    bg_sample <- bg_sample[bg_sample[["Cell.Type"]] %in% plot_categories, ]
    
    ## Factor for feature column
    bg_sample_factored <- bg_sample
    
    bg_sample_factored[, "Cell.Type"] <- factor(bg_sample_factored[, "Cell.Type"],
                                                levels = plot_categories)
    
    ## Plot
    fig <- plot_ly(bg_sample_factored,
                   type = "scatter3d",
                   mode = 'markers',
                   x = ~Cell.X.Position,
                   y = ~Cell.Y.Position,
                   z = ~Cell.Z.Position,
                   color = ~Cell.Type,
                   colors = plot_colours,
                   marker = list(size = 2))
    
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                       yaxis = list(title = 'y'),
                                       zaxis = list(title = 'z')))
    print(fig)
  }
  
  return (bg_sample)
}


simulate_sphere_cluster <- function(bg_sample, cluster_properties) {
  
  
  # Get sphere properties
  cell_type <- cluster_properties$name_of_cluster_cell
  infiltration_types <- cluster_properties$infiltration_types
  infiltration_proportions <- cluster_properties$infiltration_proportions
  radius <- cluster_properties$radius
  centre_loc <- cluster_properties$centre_loc
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  for (i in seq_len(n_cells)) {
    # Get x, y, z and phenotype of ith cell
    x <- bg_sample[i, "Cell.X.Position"]
    y <- bg_sample[i, "Cell.Y.Position"]
    z <- bg_sample[i, "Cell.Z.Position"]
    pheno <- bg_sample[i, "Cell.Type"]
    
    R <- (radius * runif(1, min = 0.7, max = 1.3))^2
    
    D <- (x - centre_loc[1])^2 + (y - centre_loc[2])^2 + (z - centre_loc[3])^2
    
    if (D < R){ 
      # in the region of cluster, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_infiltration_types <- length(infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this cluster
      pheno <- cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_infiltration_types){
        current_p <- current_p + infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- infiltration_types[n]
          break
        }
        n <- n+1
      }
    }
    bg_sample[i, "Cell.Type"] <- pheno
  }
  return(bg_sample)
}



simulate_sphere_dr <- function(bg_sample, dr_properties) {
  
  # Get sphere double ring properties
  cell_type <- dr_properties$name_of_cluster_cell
  infiltration_types <- dr_properties$infiltration_types
  infiltration_proportions <- dr_properties$infiltration_proportions
  radius <- dr_properties$radius
  centre_loc <- dr_properties$centre_loc
  
  inner_ring_cell_type <- dr_properties$name_of_inner_ring_cell
  inner_ring_width <- dr_properties$inner_ring_width
  inner_ring_infiltration_types <- dr_properties$inner_ring_infiltration_types
  inner_ring_infiltration_proportions <- dr_properties$inner_ring_infiltration_proportions
  
  outer_ring_cell_type <- dr_properties$name_of_outer_ring_cell
  outer_ring_width <- dr_properties$outer_ring_width
  outer_ring_infiltration_types <- dr_properties$outer_ring_infiltration_types
  outer_ring_infiltration_proportions <- dr_properties$outer_ring_infiltration_proportions
  
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  for (i in seq_len(n_cells)) {
    # Get x, y, z and phenotype of ith cell
    x <- bg_sample[i, "Cell.X.Position"]
    y <- bg_sample[i, "Cell.Y.Position"]
    z <- bg_sample[i, "Cell.Z.Position"]
    pheno <- bg_sample[i, "Cell.Type"]
    
    R1 <- radius^2
    R2 <- (radius + inner_ring_width)^2
    R3 <- (radius + inner_ring_width + outer_ring_width)^2
    
    D <- (x - centre_loc[1])^2 + (y - centre_loc[2])^2 + (z - centre_loc[3])^2
    
    if (D < R1) { 
      # in the region of cluster, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_infiltration_types <- length(infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this cluster
      pheno <- cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_infiltration_types){
        current_p <- current_p + infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- infiltration_types[n]
          break
        }
        n <- n + 1
      }
    }
    else if (D < R2) {
      # in the region of inner ring, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_inner_ring_infiltration_types <- length(inner_ring_infiltration_types)
      
      # default `Cell.Type` is cell type of interest of inner ring
      pheno <- inner_ring_cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_inner_ring_infiltration_types){
        current_p <- current_p + inner_ring_infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- inner_ring_infiltration_types[n]
          break
        }
        n <- n + 1
      }
    }
    else if (D < R3) {
      # in the region of outer ring, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_outer_ring_infiltration_types <- length(outer_ring_infiltration_types)
      
      # default `Cell.Type` is cell type of interest of outer ring
      pheno <- outer_ring_cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_outer_ring_infiltration_types){
        current_p <- current_p + outer_ring_infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- outer_ring_infiltration_types[n]
          break
        }
        n <- n + 1
      }
    }
    bg_sample[i, "Cell.Type"] <- pheno
    
  }
  return(bg_sample)
}



simulate_sphere_ring <- function(bg_sample, ring_properties) {
  
  # Get sphere immune ring properties
  cell_type <- ring_properties$name_of_cluster_cell
  infiltration_types <- ring_properties$infiltration_types
  infiltration_proportions <- ring_properties$infiltration_proportions
  radius <- ring_properties$radius
  centre_loc <- ring_properties$centre_loc
  
  ring_cell_type <- ring_properties$name_of_ring_cell
  ring_width <- ring_properties$ring_width
  ring_infiltration_types <- ring_properties$ring_infiltration_types
  ring_infiltration_proportions <- ring_properties$ring_infiltration_proportions
  
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  for (i in seq_len(n_cells)) {
    # Get x, y, z and phenotype of ith cell
    x <- bg_sample[i, "Cell.X.Position"]
    y <- bg_sample[i, "Cell.Y.Position"]
    z <- bg_sample[i, "Cell.Z.Position"]
    pheno <- bg_sample[i, "Cell.Type"]
    
    R1 <- radius^2
    R2 <- (radius + ring_width)^2
    
    D <- (x - centre_loc[1])^2 + (y - centre_loc[2])^2 + (z - centre_loc[3])^2
    
    if (D < R1) { 
      # in the region of cluster, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_infiltration_types <- length(infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this cluster
      pheno <- cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_infiltration_types){
        current_p <- current_p + infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- infiltration_types[n]
          break
        }
        n <- n + 1
      }
    }
    else if (D < R2) {
      # in the region of ring, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_ring_infiltration_types <- length(ring_infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this ring
      pheno <- ring_cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_ring_infiltration_types){
        current_p <- current_p + ring_infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- ring_infiltration_types[n]
          break
        }
        n <- n + 1
      }
    }
    bg_sample[i, "Cell.Type"] <- pheno
    
  }
  return(bg_sample)
}


plot_cell_categories3D <- function(data,
                                   cell_types_of_interest = NULL,
                                   colour_vector = NULL,
                                   size = 2,
                                   include_cell_types_of_no_interest = FALSE,
                                   feature_colname = "Cell.Type") {
  
  if (is.null(cell_types_of_interest)) {
    cell_types_of_interest <- unique(data$Cell.Type)
  }
  
  if (is.null(colour_vector)) {
    colour_vector <- hcl.colors(length(cell_types_of_interest), "Batlow")
  }
  
  if (length(cell_types_of_interest) != length(colour_vector)) {
    stop("Length of cell_types_of_interest is not equal to length of colour_vector")
  }
  
  ## Including non-interest cell types
  ## Define cell.id of non-interest cell types as "No Interest"
  cell_types_of_non_interest <- c()
  if (include_cell_types_of_no_interest) {
    cell_types_of_non_interest <- setdiff(unique(data[[feature_colname]]), cell_types_of_interest)
    
    data[data[[feature_colname]] %in% cell_types_of_non_interest, feature_colname] <- "No Interest"
    
    ## Add "No Interest" as a cell type of interest
    cell_types_of_interest <- c(cell_types_of_interest, "No Interest")
    
    ## Use lightgray for "No Interest" cell types
    colour_vector <- c(colour_vector, "#F0F0F0")
  }
  ## Excluding non-interest cell types
  ## Subset data to only include cell types of interest
  else {
    data <- data[data[[feature_colname]] %in% cell_types_of_interest, ]
  }
  
  ## Factor for feature column
  data[, feature_colname] <- factor(data[, feature_colname],
                                    levels = cell_types_of_interest)
  
  ## Plot
  fig <- plot_ly(data,
                 type = "scatter3d",
                 mode = 'markers',
                 x = ~Cell.X.Position,
                 y = ~Cell.Y.Position,
                 z = ~Cell.Z.Position,
                 color = ~Cell.Type,
                 colors = colour_vector,
                 marker = list(size = size))
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                     yaxis = list(title = 'y'),
                                     zaxis = list(title = 'z')))
  
  return (fig)
  
}




## Prim's algorithm function
# Input is the adjacency matrix of the graph (i.e. output from -1 * apcluster::negDistMat(df of coords))
prim <- function(graph) {
  
  # Number of vertices is number of points
  num_vertices <- nrow(graph)
  
  # Start with no vertices selected except first
  selected <- rep(FALSE, num_vertices)
  selected[1] <- TRUE
  
  # Create tree_edge matrix. Currently zero, each row represents the two vertices the edge joins
  tree_edges <- matrix(0, 
                       nrow = num_vertices - 1,
                       ncol = 2)
  
  
  # Iterate until we select enough edges (one less than the number of vertices for a MST)
  num_edges <- 0
  while (num_edges < num_vertices - 1) {
    min_weight <- Inf
    min_vertex <- -1
    
    # Iterate through each currently selected vertex
    for (i in seq(num_vertices)) {
      
      # Found a currently selected vertex
      if (selected[i] == TRUE) {
        
        # Iterate through each unselected vertex and find the nearest one
        for (j in seq(num_vertices)) {
          if (!selected[j] && graph[i, j] < min_weight) {
            min_weight <- graph[i, j]
            min_vertex <- j
            curr_vertex <- i
          }
        }
      }
    }
    
    # Current edge connects the min_vertex and ???
    tree_edges[num_edges + 1, ] <- c(min_vertex, curr_vertex)
    selected[min_vertex] <- TRUE
    num_edges <- num_edges + 1
  }
  return (tree_edges)
}




simulate_network_cluster <- function(bg_sample, cluster_properties) {  
  
  # Get network properties
  n_edges <- cluster_properties$n_edges
  cell_type <- cluster_properties$name_of_cluster_cell
  infiltration_types <- cluster_properties$infiltration_types
  infiltration_proportions <- cluster_properties$infiltration_proportions
  width <- cluster_properties$width
  centre_loc <- cluster_properties$centre_loc
  radius <- cluster_properties$radius
  
  
  # number of vertices is always one more than the number of edges for the MST will we make
  n_vertices <- n_edges + 1 
  
  ## Subset coordinate within the radius of the centre_loc
  R <- radius^2
  
  D <- (bg_sample$Cell.X.Position - centre_loc[1])^2 +
    (bg_sample$Cell.Y.Position - centre_loc[2])^2 +
    (bg_sample$Cell.Z.Position - centre_loc[3])^2
  
  chosen_cells <- bg_sample[D <= R, ]
  
  
  ## Subset further and pick n_vertices number of cells to represent the vertices
  chosen_rows <- sample(seq(nrow(chosen_cells)), n_vertices)
  
  chosen_cells <- chosen_cells[chosen_rows,
                               c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  
  ## Get adjacency matrix from points (pairwise distance between points)
  # Assume all points have an edge between each other
  # Assume weight of each edge is equal to the distance between points
  adj_mat <- -1 * apcluster::negDistMat(chosen_cells)
  
  ## Use prim's algorithm to get edges (i.e. the cells connected by each edge)
  tree_edges <- prim(adj_mat)
  
  
  
  ### Determine width of cylinders so that cylinders further away are thinner
  tree_edges <- data.frame(tree_edges)
  colnames(tree_edges) <- c("Cell1", "Cell2")
  tree_edges$Depth <- NA # If tree_edge is not NA, we have already accounted for it
  
  # Get cells on the 'outskirts' of MST (i.e. leaf_vertices)
  tree_vertices <- c(tree_edges[ , 1], tree_edges[ , 2])
  
  leaf_vertices <- names(table(tree_vertices))[table(tree_vertices) == 1]
  leaf_vertices <- as.numeric(leaf_vertices)
  
  # Start with leaf_vertices
  curr_vertices <- leaf_vertices
  curr_depth <- 1
  
  while (NA %in% tree_edges$Depth) {
    
    # New vertices will be those adjacent to the current vertices
    new_vertices <- c()
    
    # Check each current vertex
    for (vertex in curr_vertices) {
      # Start with Cell1
      curr_edges <- which(tree_edges$Cell1 == vertex)
      
      tree_edges[curr_edges, "Depth"][is.na(tree_edges[curr_edges, "Depth"])] <- curr_depth
      
      new_vertices <- c(new_vertices, tree_edges[curr_edges, "Cell2"])
      
      # Then Cell2
      curr_edges <- which(tree_edges$Cell2 == vertex)
      
      tree_edges[curr_edges, "Depth"][is.na(tree_edges[curr_edges, "Depth"])] <- curr_depth
      
      new_vertices <- c(new_vertices, tree_edges[curr_edges, "Cell1"])
      
      # Only keep unique vertices
      new_vertices <- unique(new_vertices)
    }
    
    curr_depth <- curr_depth + 1
    curr_vertices <- new_vertices
  }
  
  
  
  
  ## Get cluster properties using edge data
  cluster_properties <- list()
  max_depth <- max(tree_edges[["Depth"]])
  
  for (i in seq(n_vertices - 1)) {
    start_loc <- as.numeric(chosen_cells[tree_edges[i, "Cell1"], ])
    end_loc <- as.numeric(chosen_cells[tree_edges[i, "Cell2"], ])
    curr_width <- (1 - 0.10 * (max_depth - tree_edges[i, "Depth"])) * width # 10% decrease with each depth
    
    # Very unlikely case when width is negative, just ignore these cylinders
    if (width < 0) {
      width <- 0
    }
    
    cluster_properties[[i]] <- list(name_of_cluster_cell = cell_type,
                                    infiltration_types = infiltration_types,
                                    infiltration_proportions = infiltration_proportions,
                                    shape = "Cylinder",
                                    radius = curr_width,
                                    start_loc = start_loc,
                                    end_loc = end_loc)
  }
  
  network_bg <- simulate_clusters3D(bg_sample,
                                    n_clusters = n_edges,
                                    cluster_properties = cluster_properties,
                                    plot_image = F)
  
  return (network_bg)
  
}



simulate_network_ring <- function(bg_sample, ring_properties) {  
  
  # Get network ring properties
  n_edges <- ring_properties$n_edges
  cell_type <- ring_properties$name_of_cluster_cell
  infiltration_types <- ring_properties$infiltration_types
  infiltration_proportions <- ring_properties$infiltration_proportions
  width <- ring_properties$width
  centre_loc <- ring_properties$centre_loc
  radius <- ring_properties$radius
  
  ring_cell_type <- ring_properties$name_of_ring_cell
  ring_width <- ring_properties$ring_width
  ring_infiltration_types <- ring_properties$ring_infiltration_types
  ring_infiltration_proportions <- ring_properties$ring_infiltration_proportions
  
  
  # number of vertices is always one more than the number of edges for the MST will we make
  n_vertices <- n_edges + 1 
  
  
  ## Subset coordinate within the radius of the centre_loc
  R <- radius^2
  
  D <- (bg_sample$Cell.X.Position - centre_loc[1])^2 +
    (bg_sample$Cell.Y.Position - centre_loc[2])^2 +
    (bg_sample$Cell.Z.Position - centre_loc[3])^2
  
  chosen_cells <- bg_sample[D <= R, ]
  
  ## Subset further and pick n_vertices number of cells to represent the vertices
  chosen_rows <- sample(seq(nrow(chosen_cells)), n_vertices)
  chosen_cells <- chosen_cells[chosen_rows,
                               c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  
  ## Get adjacency matrix from points (pairwise distance between points)
  # Assume all points have an edge between each other
  # Assume weight of each edge is equal to the distance between points
  adj_mat <- -1 * apcluster::negDistMat(chosen_cells)
  
  ## Use prim's algorithm to get edges (i.e. the cells connected by each edge)
  tree_edges <- prim(adj_mat)
  
  
  
  ### Determine width of cylinders so that cylinders further away are thinner
  tree_edges <- data.frame(tree_edges)
  colnames(tree_edges) <- c("Cell1", "Cell2")
  tree_edges$Depth <- NA # If tree_edge is not NA, we have already accounted for it
  
  # Get cells on the 'outskirts' of MST (i.e. leaf_vertices)
  tree_vertices <- c(tree_edges[ , 1], tree_edges[ , 2])
  
  leaf_vertices <- names(table(tree_vertices))[table(tree_vertices) == 1]
  leaf_vertices <- as.numeric(leaf_vertices)
  
  # Start with leaf_vertices
  curr_vertices <- leaf_vertices
  curr_depth <- 1
  
  while (NA %in% tree_edges$Depth) {
    
    # New vertices will be those adjacent to the current vertices
    new_vertices <- c()
    
    # Check each current vertex
    for (vertex in curr_vertices) {
      # Start with Cell1
      curr_edges <- which(tree_edges$Cell1 == vertex)
      
      tree_edges[curr_edges, "Depth"][is.na(tree_edges[curr_edges, "Depth"])] <- curr_depth
      
      new_vertices <- c(new_vertices, tree_edges[curr_edges, "Cell2"])
      
      # Then Cell2
      curr_edges <- which(tree_edges$Cell2 == vertex)
      
      tree_edges[curr_edges, "Depth"][is.na(tree_edges[curr_edges, "Depth"])] <- curr_depth
      
      new_vertices <- c(new_vertices, tree_edges[curr_edges, "Cell1"])
      
      # Only keep unique vertices
      new_vertices <- unique(new_vertices)
    }
    
    curr_depth <- curr_depth + 1
    curr_vertices <- new_vertices
  }
  
  
  
  
  ## Get cluster properties using edge data
  ring_properties <- list()
  max_depth <- max(tree_edges[["Depth"]])
  
  for (i in seq(n_edges)) {
    start_loc <- as.numeric(chosen_cells[tree_edges[i, "Cell1"], ])
    end_loc <- as.numeric(chosen_cells[tree_edges[i, "Cell2"], ])
    curr_width <- (1 - 0.10 * (max_depth - tree_edges[i, "Depth"])) * width # 10% decrease with each depth
    
    # Very unlikely case when width is negative, just ignore these cylinders
    if (width < 0) {
      width <- 0
    }
    
    ring_properties[[i]] <- list(name_of_cluster_cell = cell_type,
                                 infiltration_types = infiltration_types,
                                 infiltration_proportions = infiltration_proportions,
                                 shape = "Cylinder",
                                 radius = curr_width,
                                 start_loc = start_loc,
                                 end_loc = end_loc,
                                 name_of_ring_cell = ring_cell_type,
                                 ring_width = ring_width,
                                 ring_infiltration_types = ring_infiltration_types,
                                 ring_infiltration_proportions = ring_infiltration_proportions)
  }
  
  network_bg <- simulate_rings3D(bg_sample,
                                 n_ring = n_edges,
                                 ring_properties = ring_properties,
                                 plot_image = F)
  
  return (network_bg)
  
}
