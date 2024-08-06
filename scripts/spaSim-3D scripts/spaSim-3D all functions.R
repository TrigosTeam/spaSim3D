library(plotly)
library(SpatialExperiment)
library(dplyr)

### Utility -------------------------------------------------------------------
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
  
  df <- data.frame("Cell.X.Position" = x, 
                   "Cell.Y.Position" = y, 
                   "Cell.Z.Position" = z)
  
  return(df)
}


prims_algorithm <- function(graph) {
  
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
  return(tree_edges)
}

plot_cells3D <- function(spe,
                         plot_cell_types = NULL,
                         plot_colours = NULL,
                         feature_colname = "Cell.Type") {
  
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), "Cell.Type" = spe[[feature_colname]])
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) {
    plot_cell_types <- unique(df[["Cell.Type"]])
  }
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(plot_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following plot_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  ## If no colours inputted, use rainbow palette
  if (is.null(plot_colours)) {
    plot_colours <- rainbow(length(plot_cell_types))
  }
  
  ## User inputs mismatching cell types and colours
  if (length(plot_cell_types) != length(plot_colours)) {
    stop("Length of plot_cell_types is not equal to length of plot_colours")
  }
  
  ## Factor for feature column
  df[, "Cell.Type"] <- factor(df[, "Cell.Type"],
                              levels = plot_cell_types)
  
  ## Plot
  fig <- plot_ly(df,
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
  
  return (fig)
}



### Main simulation functions -------------------------------------------------
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


simulate_ordered_background_cells3D <- function(n_cells, 
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
  df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")
  
  # Get meta data
  background_metadata <- list("background_type" = "ordered",
                              "n_cells" = n_cells,
                              "length" = length,
                              "width" = width,
                              "height" = height,
                              "jitter_proportion" = jitter_proportion,
                              "cell_types" = background_cell_type,
                              "cell_proportions" = 1)
  simulation_metadata <- list(background = background_metadata)
  
  ## Convert data frame to spe object
  spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = list(simulation = simulation_metadata))
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        background_cell_type,
                        "lightgray")
    print(fig)
  }
  
  return(spe)
}



simulate_mixing3D <- function(bg_spe,
                              cell_types,
                              cell_proportions,
                              plot_image = TRUE,
                              plot_cell_types = NULL,
                              plot_colours = NULL) {
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cell_types) != length(cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cell_proportions < 0 | cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (!all.equal(sum(cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  
  bg_spe[["Cell.Type"]] <- sample(cell_types, size = ncol(bg_spe), replace = TRUE, prob = cell_proportions)
  
  bg_spe@metadata[["simulation"]][["background"]][["cell_types"]] <- cell_types
  bg_spe@metadata[["simulation"]][["background"]][["cell_proportions"]] <- cell_proportions
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(bg_spe,
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(bg_spe)
}

simulate_clusters3D <- function(bg_spe,
                                cluster_properties,
                                plot_image = TRUE,
                                plot_cell_types = NULL,
                                plot_colours = NULL) {
  
  
  for (k in seq(length(cluster_properties))) { 
    
    # For each cluster, get the shape
    shape <- cluster_properties[[k]]$shape
    
    ### Sphere shape
    if (shape == "Sphere") {
      bg_spe <- simulate_sphere_cluster(bg_spe, cluster_properties[[k]])
    } 
    
    ### Ellipsoid shape
    if (shape == "Ellipsoid") {
      bg_spe <- simulate_ellipsoid_cluster(bg_spe, cluster_properties[[k]])
    }
    
    ### Cylinder shape
    if (shape == "Cylinder") {
      bg_spe <- simulate_cylinder_cluster(bg_spe, cluster_properties[[k]])
    }
    
    ### Network shape
    if (shape == "Network") {
      bg_spe <- simulate_network_cluster(bg_spe, cluster_properties[[k]])
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(bg_spe, 
                        plot_cell_types,
                        plot_colours)
    print(fig)
  }
  
  return(bg_spe)
}

simulate_rings3D <- function(bg_spe,
                             ring_properties,
                             plot_image = TRUE,
                             plot_cell_types = NULL,
                             plot_colours = NULL) {
  
  for (k in seq(length(ring_properties))) { 
    
    # For each cluster, get the shape
    shape <- ring_properties[[k]]$shape
    
    ### Sphere shape +  ring
    if (shape == "Sphere") {
      bg_spe <- simulate_sphere_ring(bg_spe, ring_properties[[k]])
    } 
    
    ### Ellipsoid shape + ring
    else if (shape == "Ellipsoid") {
      bg_spe <- simulate_ellipsoid_ring(bg_spe, ring_properties[[k]])
    }
    
    ### Cylinder shape + ring
    else if (shape == "Cylinder") {
      bg_spe <- simulate_cylinder_ring(bg_spe, ring_properties[[k]])
    }
    
    ### Network shape + ring
    else if (shape == "Network") {
      bg_spe <- simulate_network_ring(bg_spe, ring_properties[[k]])
    }
    
    else {
      stop("Invalid shape")
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(bg_spe, 
                        plot_cell_types,
                        plot_colours)
    print(fig)
  }
  
  return(bg_spe)
}

simulate_double_rings3D <- function(bg_spe,
                                    dr_properties,
                                    plot_image = TRUE,
                                    plot_cell_types = NULL,
                                    plot_colours = NULL) {
  
  for (k in seq(length(dr_properties))) { 
    
    # For each cluster, get the shape
    shape <- dr_properties[[k]]$shape
    
    ### Sphere double ring shape
    if (shape == "Sphere") {
      bg_spe <- simulate_sphere_dr(bg_spe, dr_properties[[k]])
    } 
    
    ### Ellipsoid double ring shape
    if (shape == "Ellipsoid") {
      bg_spe <- simulate_ellipsoid_dr(bg_spe, dr_properties[[k]])
    }
    
    ### Cylinder double ring shape
    if (shape == "Cylinder") {
      bg_spe <- simulate_cylinder_dr(bg_spe, dr_properties[[k]])
    }
    
    ### Network double ring shape
    if (shape == "Network") {
      bg_spe <- simulate_network_dr(bg_spe, dr_properties[[k]])
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(bg_spe, 
                        plot_cell_types,
                        plot_colours)
    print(fig)
  }
  
  return(bg_spe)
}




### Sphere --------------------------------------------------------------------
simulate_sphere_cluster <- function(bg_spe, cluster_properties) {
  
  # Get sphere properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  radius <- cluster_properties$radius
  centre_loc <- cluster_properties$centre_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  ## Change cell types in the sphere cluster
  spe_coords <- data.frame(spatialCoords(bg_spe))
  
  bg_spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                    (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                    (spe_coords$Cell.Z.Position - centre_loc[3])^2 < radius^2,
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- cluster_properties
  
  return(bg_spe)
}

simulate_sphere_ring <- function(bg_spe, ring_properties) {
  
  # Get sphere ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  radius <- ring_properties$radius
  centre_loc <- ring_properties$centre_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  ## Check number of ring cell types matches the number of cell proportions
  if (length(ring_cell_types) != length(ring_cell_proportions)) stop("Number of ring cell types doesn't match number of ring cell proportion.")
  
  ## Check ring cell proportions are not negative or greater than 1
  if (sum(ring_cell_proportions < 0 | ring_cell_proportions > 1) != 0) stop("Ring cell proportions cannot be negative or greater than 1")
  
  ## Check ring cell proportions add up to 1
  if (!all.equal(sum(ring_cell_proportions), 1)) stop("Sum of ring cell proportions is NOT 1")
  
  ## Change cell types in the sphere ringed cluster
  spe_coords <- data.frame(spatialCoords(bg_spe))
  
  # Start with cells in ring  
  bg_spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                    (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                    (spe_coords$Cell.Z.Position - centre_loc[3])^2 < (radius + ring_width)^2,
                                  sample(ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Then do cells in the cluster 
  bg_spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                    (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                    (spe_coords$Cell.Z.Position - centre_loc[3])^2 < radius^2,
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- ring_properties
  
  return(bg_spe)
}

simulate_sphere_dr <- function(bg_spe, dr_properties) {
  
  # Get sphere double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  radius <- dr_properties$radius
  centre_loc <- dr_properties$centre_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  
  ## Check number of inner ring cell types matches the number of inner ring cell proportions
  if (length(inner_ring_cell_types) != length(inner_ring_cell_proportions)) stop("Number of inner ring cell types doesn't match number of inner ring cell proportion.")
  
  ## Check inner ring cell proportions are not negative or greater than 1
  if (sum(inner_ring_cell_proportions < 0 | inner_ring_cell_proportions > 1) != 0) stop("Inner ring cell proportions cannot be negative or greater than 1")
  
  ## Check inner ring cell proportions add up to 1
  if (!all.equal(sum(inner_ring_cell_proportions), 1)) stop("Sum of inner ring cell proportions is NOT 1")
  
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  ## Check number of outer ring cell types matches the number of outer ring cell proportions
  if (length(outer_ring_cell_types) != length(outer_ring_cell_proportions)) stop("Number of outer ring cell types doesn't match number of outer ring cell proportion.")
  
  ## Check outer ring cell proportions are not negative or greater than 1
  if (sum(outer_ring_cell_proportions < 0 | outer_ring_cell_proportions > 1) != 0) stop("Outer ring cell proportions cannot be negative or greater than 1")
  
  ## Check outer ring cell proportions add up to 1
  if (!all.equal(sum(outer_ring_cell_proportions), 1)) stop("Sum of outer ring cell proportions is NOT 1")
  
  ## Change cell types in the sphere ringed cluster
  spe_coords <- data.frame(spatialCoords(bg_spe))
  
  # Start with cells in outer ring  
  bg_spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                    (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                    (spe_coords$Cell.Z.Position - centre_loc[3])^2 < (radius + inner_ring_width + outer_ring_width)^2,
                                  sample(outer_ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = outer_ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Then do cells in inner ring  
  bg_spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                    (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                    (spe_coords$Cell.Z.Position - centre_loc[3])^2 < (radius + inner_ring_width)^2,
                                  sample(inner_ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = inner_ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Then do cells in the cluster 
  bg_spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                    (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                    (spe_coords$Cell.Z.Position - centre_loc[3])^2 < radius^2,
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  return(bg_spe)
}



### Ellipsoid -----------------------------------------------------------------
simulate_ellipsoid_cluster <- function(bg_spe, cluster_properties) {
  
  # Get ellipsoid properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  x_radius <- cluster_properties$x_radius
  y_radius <- cluster_properties$y_radius
  z_radius <- cluster_properties$z_radius
  centre_loc <- cluster_properties$centre_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  
  # Rotation angles
  theta <- cluster_properties$y_z_rotation * (pi/180) # rotation in x-axis
  alpha <- cluster_properties$x_z_rotation * (pi/180) # rotation in y-axis
  beta  <- cluster_properties$x_y_rotation * (pi/180) # rotation in z-axis
  
  # 3x3 Transformation matrix (T_M) using rotation angles
  T_M <- matrix(data = c(cos(alpha) * cos(beta), 
                         cos(alpha) * sin(beta), 
                         sin(alpha),
                         -sin(theta) * sin(alpha) * cos(beta) - cos(theta) * sin(beta),
                         -sin(theta) * sin(alpha) * sin(beta) + cos(theta) * cos(beta),
                         sin(theta) * cos(alpha),
                         -cos(theta) * sin(alpha) * cos(beta) + sin(theta) * sin(beta),
                         -cos(theta) * sin(alpha) * sin(beta) - sin(theta) * cos(beta),
                         cos(theta) * cos(alpha)), 
                nrow = 3, 
                ncol = 3, 
                byrow = TRUE)
  
  ## Change cell types in the ellipsoid cluster
  spe_coords <- data.frame(spatialCoords(bg_spe))
  
  # Adjust x, y and z coordinates relative to the ellipsoid centre
  x <- spe_coords$Cell.X.Position - centre_loc[1]
  y <- spe_coords$Cell.Y.Position - centre_loc[2]
  z <- spe_coords$Cell.Z.Position - centre_loc[3]
  
  # Transform  x, y and z coordinates using rotation transformation matrix
  x_new <- T_M[1, 1] * x + T_M[1, 2] * y + T_M[1, 3] * z
  y_new <- T_M[2, 1] * x + T_M[2, 2] * y + T_M[2, 3] * z
  z_new <- T_M[3, 1] * x + T_M[3, 2] * y + T_M[3, 3] * z
  
  bg_spe[["Cell.Type"]] <- ifelse((x_new / x_radius)^2 +
                                    (y_new / y_radius)^2 +
                                    (z_new / z_radius)^2 <= 1,
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  
  # Update current meta data
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- cluster_properties
  
  return(bg_spe)
}

simulate_ellipsoid_ring <- function(bg_spe, ring_properties) {
  
  # Get ellipsoid ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  x_radius <- ring_properties$x_radius
  y_radius <- ring_properties$y_radius
  z_radius <- ring_properties$z_radius
  centre_loc <- ring_properties$centre_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  ## Check number of ring cell types matches the number of cell proportions
  if (length(ring_cell_types) != length(ring_cell_proportions)) stop("Number of ring cell types doesn't match number of ring cell proportion.")
  
  ## Check ring cell proportions are not negative or greater than 1
  if (sum(ring_cell_proportions < 0 | ring_cell_proportions > 1) != 0) stop("Ring cell proportions cannot be negative or greater than 1")
  
  ## Check ring cell proportions add up to 1
  if (!all.equal(sum(ring_cell_proportions), 1)) stop("Sum of ring cell proportions is NOT 1")
  
  # Rotation angles
  theta <- ring_properties$y_z_rotation * (pi/180) # rotation in x-axis
  alpha <- ring_properties$x_z_rotation * (pi/180) # rotation in y-axis
  beta  <- ring_properties$x_y_rotation * (pi/180) # rotation in z-axis
  
  # 3x3 Transformation matrix using rotation angles
  T_M <- matrix(data = c(cos(alpha) * cos(beta), 
                         cos(alpha) * sin(beta), 
                         sin(alpha),
                         -sin(theta) * sin(alpha) * cos(beta) - cos(theta) * sin(beta),
                         -sin(theta) * sin(alpha) * sin(beta) + cos(theta) * cos(beta),
                         sin(theta) * cos(alpha),
                         -cos(theta) * sin(alpha) * cos(beta) + sin(theta) * sin(beta),
                         -cos(theta) * sin(alpha) * sin(beta) - sin(theta) * cos(beta),
                         cos(theta) * cos(alpha)), 
                nrow = 3, 
                ncol = 3, 
                byrow = TRUE)
  
  
  
  ## Change cell types in the ellipsoid cluster
  spe_coords <- data.frame(spatialCoords(bg_spe))
  
  # Adjust x, y and z coordinates relative to the ellipsoid centre
  x <- spe_coords$Cell.X.Position - centre_loc[1]
  y <- spe_coords$Cell.Y.Position - centre_loc[2]
  z <- spe_coords$Cell.Z.Position - centre_loc[3]
  
  # Transform  x, y and z coordinates using rotation transformation matrix
  x_new <- T_M[1, 1] * x + T_M[1, 2] * y + T_M[1, 3] * z
  y_new <- T_M[2, 1] * x + T_M[2, 2] * y + T_M[2, 3] * z
  z_new <- T_M[3, 1] * x + T_M[3, 2] * y + T_M[3, 3] * z
  
  # Start with cells in ring  
  bg_spe[["Cell.Type"]] <- ifelse((x_new / (x_radius + ring_width))^2 +
                                    (y_new / (y_radius + ring_width))^2 +
                                    (z_new / (z_radius + ring_width))^2 <= 1,
                                  sample(ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  
  # Then do cells in the cluster  
  bg_spe[["Cell.Type"]] <- ifelse((x_new / x_radius)^2 +
                                    (y_new / y_radius)^2 +
                                    (z_new / z_radius)^2 <= 1,
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  
  # Update current meta data
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- ring_properties
  
  return(bg_spe)
}

simulate_ellipsoid_dr <- function(bg_spe, dr_properties) {
  
  # Get ellipsoid double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  x_radius <- dr_properties$x_radius
  y_radius <- dr_properties$y_radius
  z_radius <- dr_properties$z_radius
  centre_loc <- dr_properties$centre_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  
  ## Check number of inner ring cell types matches the number of inner ring cell proportions
  if (length(inner_ring_cell_types) != length(inner_ring_cell_proportions)) stop("Number of inner ring cell types doesn't match number of inner ring cell proportion.")
  
  ## Check inner ring cell proportions are not negative or greater than 1
  if (sum(inner_ring_cell_proportions < 0 | inner_ring_cell_proportions > 1) != 0) stop("Inner ring cell proportions cannot be negative or greater than 1")
  
  ## Check inner ring cell proportions add up to 1
  if (!all.equal(sum(inner_ring_cell_proportions), 1)) stop("Sum of inner ring cell proportions is NOT 1")
  
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  ## Check number of outer ring cell types matches the number of outer ring cell proportions
  if (length(outer_ring_cell_types) != length(outer_ring_cell_proportions)) stop("Number of outer ring cell types doesn't match number of outer ring cell proportion.")
  
  ## Check outer ring cell proportions are not negative or greater than 1
  if (sum(outer_ring_cell_proportions < 0 | outer_ring_cell_proportions > 1) != 0) stop("Outer ring cell proportions cannot be negative or greater than 1")
  
  ## Check outer ring cell proportions add up to 1
  if (!all.equal(sum(outer_ring_cell_proportions), 1)) stop("Sum of outer ring cell proportions is NOT 1")
  
  # Rotation angles
  theta <- dr_properties$y_z_rotation * (pi/180) # rotation in x-axis
  alpha <- dr_properties$x_z_rotation * (pi/180) # rotation in y-axis
  beta  <- dr_properties$x_y_rotation * (pi/180) # rotation in z-axis
  
  # 3x3 Transformation matrix using rotation angles
  T_M <- matrix(data = c(cos(alpha) * cos(beta), 
                         cos(alpha) * sin(beta), 
                         sin(alpha),
                         -sin(theta) * sin(alpha) * cos(beta) - cos(theta) * sin(beta),
                         -sin(theta) * sin(alpha) * sin(beta) + cos(theta) * cos(beta),
                         sin(theta) * cos(alpha),
                         -cos(theta) * sin(alpha) * cos(beta) + sin(theta) * sin(beta),
                         -cos(theta) * sin(alpha) * sin(beta) - sin(theta) * cos(beta),
                         cos(theta) * cos(alpha)), 
                nrow = 3, 
                ncol = 3, 
                byrow = TRUE)
  
  ## Change cell types in the ellipsoid cluster
  spe_coords <- data.frame(spatialCoords(bg_spe))
  
  # Adjust x, y and z coordinates relative to the ellipsoid centre
  x_new <- spe_coords$Cell.X.Position - centre_loc[1]
  y_new <- spe_coords$Cell.Y.Position - centre_loc[2]
  z_new <- spe_coords$Cell.Z.Position - centre_loc[3]
  
  # Transform  x, y and z coordinates using rotation transformation matrix
  x <- T_M[1, 1] * x + T_M[1, 2] * y + T_M[1, 3] * z
  y <- T_M[2, 1] * x + T_M[2, 2] * y + T_M[2, 3] * z
  z <- T_M[3, 1] * x + T_M[3, 2] * y + T_M[3, 3] * z
  
  
  # Start with cells in outer ring  
  bg_spe[["Cell.Type"]] <- ifelse((x_new / (x_radius + inner_ring_width + outer_ring_width))^2 +
                                    (y_new / (y_radius + inner_ring_width + outer_ring_width))^2 +
                                    (z_new / (z_radius + inner_ring_width + outer_ring_width))^2 <= 1,
                                  sample(outer_ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = outer_ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Then do cells in inner ring  
  bg_spe[["Cell.Type"]] <- ifelse((x_new / (x_radius + inner_ring_width))^2 +
                                    (y_new / (y_radius + inner_ring_width))^2 +
                                    (z_new / (z_radius + inner_ring_width))^2 <= 1,
                                  sample(inner_ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = inner_ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  
  # Then do cells in the cluster  
  bg_spe[["Cell.Type"]] <- ifelse((x_new / x_radius)^2 +
                                    (y_new / y_radius)^2 +
                                    (z_new / z_radius)^2 <= 1,
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  return(bg_spe)
}



### Cylinder ------------------------------------------------------------------
simulate_cylinder_cluster <- function(bg_spe, cluster_properties) {
  
  # Get cylinder properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  radius <- cluster_properties$radius
  start_loc <- cluster_properties$start_loc
  end_loc <- cluster_properties$end_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1.")
  
  ## Check cell proportions add up to 1
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  ## Check if start and end coordinates of the cylinder are the same
  if (identical(start_loc, end_loc)) warning("Start and end coordinates of the cylinder are the same.")
  
  ## Change cell types in the cylinder cluster
  spe_coords <- spatialCoords(bg_spe)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  # Get vector between from each cell to start_loc
  v2 <- sweep(spe_coords, 2, end_loc, "-")
  
  bg_spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                     rowSums(sweep(spe_coords, 2, v1, "*")) > d1 & rowSums(sweep(spe_coords, 2, v1, "*")) < d2) & # Cell must be between the planes
                                    (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) < radius^2), # Cell must be close enough to the cylinder line
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- cluster_properties
  
  return(bg_spe)
}

simulate_cylinder_ring <- function(bg_spe, ring_properties) {
  
  # Get cylinder ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  radius <- ring_properties$radius
  start_loc <- ring_properties$start_loc
  end_loc <- ring_properties$end_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  ## Check number of ring cell types matches the number of cell proportions
  if (length(ring_cell_types) != length(ring_cell_proportions)) stop("Number of ring cell types doesn't match number of ring cell proportion.")
  
  ## Check ring cell proportions are not negative or greater than 1
  if (sum(ring_cell_proportions < 0 | ring_cell_proportions > 1) != 0) stop("Ring cell proportions cannot be negative or greater than 1")
  
  ## Check ring cell proportions add up to 1
  if (!all.equal(sum(ring_cell_proportions), 1)) stop("Sum of ring cell proportions is NOT 1")
  
  ## Check if start and end coordinates of the cylinder are the same
  if (identical(start_loc, end_loc)) warning("Start and end coordinates of the cylinder are the same.")
  
  ## Change cell types in the cylinder cluster
  spe_coords <- spatialCoords(bg_spe)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  # Get vector between from each cell to start_loc
  v2 <- sweep(spe_coords, 2, end_loc, "-")
  
  # Start with cells in ring
  bg_spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                     rowSums(sweep(spe_coords, 2, v1, "*")) > d1 & rowSums(sweep(spe_coords, 2, v1, "*")) < d2) & # Cell must be between the planes
                                    (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) < (radius + ring_width)^2), # Cell must be close enough to the cylinder line
                                  sample(ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Then do cells in the cluster 
  bg_spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                     rowSums(sweep(spe_coords, 2, v1, "*")) > d1 & rowSums(sweep(spe_coords, 2, v1, "*")) < d2) & # Cell must be between the planes
                                    (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) < radius^2), # Cell must be close enough to the cylinder line
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  
  # Update current meta data
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- ring_properties
  
  return(bg_spe)
}

simulate_cylinder_dr <- function(bg_spe, dr_properties) {
  
  # Get cylinder double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  radius <- dr_properties$radius
  start_loc <- dr_properties$start_loc
  end_loc <- dr_properties$end_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  
  ## Check number of inner ring cell types matches the number of inner ring cell proportions
  if (length(inner_ring_cell_types) != length(inner_ring_cell_proportions)) stop("Number of inner ring cell types doesn't match number of inner ring cell proportion.")
  
  ## Check inner ring cell proportions are not negative or greater than 1
  if (sum(inner_ring_cell_proportions < 0 | inner_ring_cell_proportions > 1) != 0) stop("Inner ring cell proportions cannot be negative or greater than 1")
  
  ## Check inner ring cell proportions add up to 1
  if (!all.equal(sum(inner_ring_cell_proportions), 1)) stop("Sum of inner ring cell proportions is NOT 1")
  
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  ## Check number of outer ring cell types matches the number of outer ring cell proportions
  if (length(outer_ring_cell_types) != length(outer_ring_cell_proportions)) stop("Number of outer ring cell types doesn't match number of outer ring cell proportion.")
  
  ## Check outer ring cell proportions are not negative or greater than 1
  if (sum(outer_ring_cell_proportions < 0 | outer_ring_cell_proportions > 1) != 0) stop("Outer ring cell proportions cannot be negative or greater than 1")
  
  ## Check outer ring cell proportions add up to 1
  if (!all.equal(sum(outer_ring_cell_proportions), 1)) stop("Sum of outer ring cell proportions is NOT 1")
  
  ## Check if start and end coordinates of the cylinder are the same
  if (identical(start_loc, end_loc)) warning("Start and end coordinates of the cylinder are the same.")
  
  ## Change cell types in the cylinder cluster
  spe_coords <- spatialCoords(bg_spe)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  # Get vector between from each cell to start_loc
  v2 <- sweep(spe_coords, 2, end_loc, "-")
  
  # Start with cells in outer ring
  bg_spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                     rowSums(sweep(spe_coords, 2, v1, "*")) > d1 & rowSums(sweep(spe_coords, 2, v1, "*")) < d2) & # Cell must be between the planes
                                    (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) < (radius + inner_ring_width + outer_ring_width)^2), # Cell must be close enough to the cylinder line
                                  sample(outer_ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = outer_ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Start with cells in inner ring
  bg_spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                     rowSums(sweep(spe_coords, 2, v1, "*")) > d1 & rowSums(sweep(spe_coords, 2, v1, "*")) < d2) & # Cell must be between the planes
                                    (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) < (radius + inner_ring_width)^2), # Cell must be close enough to the cylinder line
                                  sample(inner_ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = inner_ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Then do cells in the cluster 
  bg_spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                     rowSums(sweep(spe_coords, 2, v1, "*")) > d1 & rowSums(sweep(spe_coords, 2, v1, "*")) < d2) & # Cell must be between the planes
                                    (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) < radius^2), # Cell must be close enough to the cylinder line
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  return(bg_spe)
}




### Network -------------------------------------------------------------------
simulate_network_cluster <- function(bg_spe, cluster_properties) {  
  
  # Get network properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  n_edges <- cluster_properties$n_edges
  width <- cluster_properties$width
  centre_loc <- cluster_properties$centre_loc
  radius <- cluster_properties$radius
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(bg_spe), 
                   "Cell.Type" = bg_spe[["Cell.Type"]],
                   "Cell.ID" = bg_spe[["Cell.ID"]])
  
  # Number of vertices is always one more than the number of edges for the MST will we make
  n_vertices <- n_edges + 1 
  
  ## Choose cells within the radius of the centre_loc
  R <- radius^2
  
  D <- (df$Cell.X.Position - centre_loc[1])^2 +
    (df$Cell.Y.Position - centre_loc[2])^2 +
    (df$Cell.Z.Position - centre_loc[3])^2
  
  cells_chosen <- df[D <= R, ]
  
  ## Subset further and pick 'n_vertices' cells to represent the vertices
  cells_chosen <- sample_n(cells_chosen, n_vertices)
  
  ## Get coordinates of cells chosen for vertices
  cells_chosen <- cells_chosen[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  ## Get adjacency matrix from points (pairwise distance between points)
  # Assume all points have an edge between each other
  # Assume weight of each edge is equal to the distance between points
  adj_mat <- -1 * apcluster::negDistMat(cells_chosen)
  
  ## Use prim's algorithm to get edges (i.e. the cells connected by each edge)
  tree_edges <- prims_algorithm(adj_mat)
  
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
  network_cluster_properties <- list()
  max_depth <- max(tree_edges[["Depth"]])
  
  for (i in seq(n_vertices - 1)) {
    start_loc <- as.numeric(cells_chosen[tree_edges[i, "Cell1"], ])
    end_loc <- as.numeric(cells_chosen[tree_edges[i, "Cell2"], ])
    curr_width <- (1 - 0.10 * (max_depth - tree_edges[i, "Depth"])) * width # 10% decrease with each depth
    
    # Very unlikely case when width is negative, just ignore these cylinders
    if (width < 0) {
      width <- 0
    }
    
    network_cluster_properties[[i]] <- list(shape = "Cylinder",
                                            cluster_cell_types = cluster_cell_types,
                                            cluster_cell_proportions = cluster_cell_proportions,
                                            radius = curr_width,
                                            start_loc = start_loc,
                                            end_loc = end_loc)
  }
  
  network_spe <- simulate_clusters3D(bg_spe,
                                     cluster_properties = network_cluster_properties,
                                     plot_image = F)
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(network_spe), "Cell.Type" = network_spe[["Cell.Type"]], "Cell.ID" = network_spe[["Cell.ID"]])
  
  # Update current meta data
  metadata <- bg_spe@metadata
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  metadata[["simulation"]][[paste("cluster", length(metadata[["simulation"]]), sep="_")]] <- cluster_properties
  
  # Convert data frame to spe object
  cluster_spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = metadata)
  
  return(cluster_spe)
}


simulate_network_ring <- function(bg_spe, ring_properties) {  
  
  # Get network ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  n_edges <- ring_properties$n_edges
  width <- ring_properties$width
  centre_loc <- ring_properties$centre_loc
  radius <- ring_properties$radius
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  ## Check number of ring cell types matches the number of cell proportions
  if (length(ring_cell_types) != length(ring_cell_proportions)) stop("Number of ring cell types doesn't match number of ring cell proportion.")
  
  ## Check ring cell proportions are not negative or greater than 1
  if (sum(ring_cell_proportions < 0 | ring_cell_proportions > 1) != 0) stop("Ring cell proportions cannot be negative or greater than 1")
  
  ## Check ring cell proportions add up to 1
  if (!all.equal(sum(ring_cell_proportions), 1)) stop("Sum of ring cell proportions is NOT 1")
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(bg_spe), 
                   "Cell.Type" = bg_spe[["Cell.Type"]],
                   "Cell.ID" = bg_spe[["Cell.ID"]])
  
  # number of vertices is always one more than the number of edges for the MST will we make
  n_vertices <- n_edges + 1 
  
  ## Subset coordinate within the radius of the centre_loc
  R <- radius^2
  
  D <- (df$Cell.X.Position - centre_loc[1])^2 +
    (df$Cell.Y.Position - centre_loc[2])^2 +
    (df$Cell.Z.Position - centre_loc[3])^2
  
  cells_chosen <- df[D <= R, ]
  
  ## Subset further and pick 'n_vertices' cells to represent the vertices
  cells_chosen <- sample_n(cells_chosen, n_vertices)
  
  ## Get coordinates of cells chosen for vertices
  cells_chosen <- cells_chosen[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  ## Get adjacency matrix from points (pairwise distance between points)
  # Assume all points have an edge between each other
  # Assume weight of each edge is equal to the distance between points
  adj_mat <- -1 * apcluster::negDistMat(cells_chosen)
  
  ## Use prim's algorithm to get edges (i.e. the cells connected by each edge)
  tree_edges <- prims_algorithm(adj_mat)
  
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
  network_ring_properties <- list()
  max_depth <- max(tree_edges[["Depth"]])
  
  for (i in seq(n_edges)) {
    start_loc <- as.numeric(cells_chosen[tree_edges[i, "Cell1"], ])
    end_loc <- as.numeric(cells_chosen[tree_edges[i, "Cell2"], ])
    curr_width <- (1 - 0.10 * (max_depth - tree_edges[i, "Depth"])) * width # 10% decrease with each depth
    
    # Very unlikely case when width is negative, just ignore these cylinders
    if (width < 0) {
      width <- 0
    }
    
    network_ring_properties[[i]] <- list(shape = "Cylinder",
                                         cluster_cell_types = cluster_cell_types,
                                         cluster_cell_proportions = cluster_cell_proportions,
                                         radius = curr_width,
                                         start_loc = start_loc,
                                         end_loc = end_loc,
                                         ring_cell_types = ring_cell_types,
                                         ring_cell_proportions = ring_cell_proportions,
                                         ring_width = ring_width)
  }
  
  network_spe <- simulate_rings3D(bg_spe,
                                  ring_properties = network_ring_properties,
                                  plot_image = F)
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(network_spe), "Cell.Type" = network_spe[["Cell.Type"]], "Cell.ID" = network_spe[["Cell.ID"]])
  
  # Update current meta data
  metadata <- bg_spe@metadata
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  metadata[["simulation"]][[paste("cluster", length(metadata[["simulation"]]), sep="_")]] <- ring_properties
  
  # Convert data frame to spe object
  cluster_spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = metadata)
  
  return(cluster_spe)
  
}


simulate_network_dr <- function(bg_spe, dr_properties) {  
  
  # Get network double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  n_edges <- dr_properties$n_edges
  width <- dr_properties$width
  centre_loc <- dr_properties$centre_loc
  radius <- dr_properties$radius
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  
  ## Check number of inner ring cell types matches the number of inner ring cell proportions
  if (length(inner_ring_cell_types) != length(inner_ring_cell_proportions)) stop("Number of inner ring cell types doesn't match number of inner ring cell proportion.")
  
  ## Check inner ring cell proportions are not negative or greater than 1
  if (sum(inner_ring_cell_proportions < 0 | inner_ring_cell_proportions > 1) != 0) stop("Inner ring cell proportions cannot be negative or greater than 1")
  
  ## Check inner ring cell proportions add up to 1
  if (!all.equal(sum(inner_ring_cell_proportions), 1)) stop("Sum of inner ring cell proportions is NOT 1")
  
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  ## Check number of outer ring cell types matches the number of outer ring cell proportions
  if (length(outer_ring_cell_types) != length(outer_ring_cell_proportions)) stop("Number of outer ring cell types doesn't match number of outer ring cell proportion.")
  
  ## Check outer ring cell proportions are not negative or greater than 1
  if (sum(outer_ring_cell_proportions < 0 | outer_ring_cell_proportions > 1) != 0) stop("Outer ring cell proportions cannot be negative or greater than 1")
  
  ## Check outer ring cell proportions add up to 1
  if (!all.equal(sum(outer_ring_cell_proportions), 1)) stop("Sum of outer ring cell proportions is NOT 1")
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(bg_spe), 
                   "Cell.Type" = bg_spe[["Cell.Type"]],
                   "Cell.ID" = bg_spe[["Cell.ID"]])
  
  # number of vertices is always one more than the number of edges for the MST will we make
  n_vertices <- n_edges + 1 
  
  ## Subset coordinate within the radius of the centre_loc
  R <- radius^2
  
  D <- (df$Cell.X.Position - centre_loc[1])^2 +
    (df$Cell.Y.Position - centre_loc[2])^2 +
    (df$Cell.Z.Position - centre_loc[3])^2
  
  cells_chosen <- df[D <= R, ]
  
  ## Subset further and pick 'n_vertices' cells to represent the vertices
  cells_chosen <- sample_n(cells_chosen, n_vertices)
  
  ## Get coordinates of cells chosen for vertices
  cells_chosen <- cells_chosen[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  ## Get adjacency matrix from points (pairwise distance between points)
  # Assume all points have an edge between each other
  # Assume weight of each edge is equal to the distance between points
  adj_mat <- -1 * apcluster::negDistMat(cells_chosen)
  
  ## Use prim's algorithm to get edges (i.e. the cells connected by each edge)
  tree_edges <- prims_algorithm(adj_mat)
  
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
  network_dr_properties <- list()
  max_depth <- max(tree_edges[["Depth"]])
  
  for (i in seq(n_edges)) {
    start_loc <- as.numeric(cells_chosen[tree_edges[i, "Cell1"], ])
    end_loc <- as.numeric(cells_chosen[tree_edges[i, "Cell2"], ])
    curr_width <- (1 - 0.10 * (max_depth - tree_edges[i, "Depth"])) * width # 10% decrease with each depth
    
    # Very unlikely case when width is negative, just ignore these cylinders
    if (width < 0) {
      width <- 0
    }
    
    network_dr_properties[[i]] <- list(shape = "Cylinder",
                                       cluster_cell_types = cluster_cell_types,
                                       cluster_cell_proportions = cluster_cell_proportions,
                                       radius = curr_width,
                                       start_loc = start_loc,
                                       end_loc = end_loc,
                                       inner_ring_cell_types = inner_ring_cell_types,
                                       inner_ring_cell_proportions = inner_ring_cell_proportions,
                                       inner_ring_width = inner_ring_width,
                                       outer_ring_cell_types = outer_ring_cell_types,
                                       outer_ring_cell_proportions = outer_ring_cell_proportions,
                                       outer_ring_width = outer_ring_width)
  }
  
  network_spe <- simulate_double_rings3D(bg_spe,
                                         dr_properties = network_dr_properties,
                                         plot_image = F)
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(network_spe), "Cell.Type" = network_spe[["Cell.Type"]], "Cell.ID" = network_spe[["Cell.ID"]])
  
  # Update current meta data
  metadata <- bg_spe@metadata
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  metadata[["simulation"]][[paste("cluster", length(metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  # Convert data frame to spe object
  cluster_spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = metadata)
  
  return(cluster_spe)
  
}




### Message background integrator strings and function ------------------------------------------------------------
message_background <- paste("Hello spaSim-3D user, how do you want your background cells to look like?\n
          1. Random pattern\n
          2. Ordered pattern\n\n",
                            "In a random pattern, cells are placed randomly...\n",
                            "In a ordered pattern, cells follow a regularly spaced in a hexagonal grid\n",
                            "To choose, please enter 1 or 2.\n", sep = "")

message_background_random <- paste("We will need a few parameters before we can obtain the simulation\n",
                                   "    Window size - length, width and height (e.g. 100 x 100 x 100)\n",
                                   "    Number of cells (e.g. 10000 cells)\n",
                                   "    Minimum distance between cells (e.g. minimum distance of 2)\n",
                                   "If you want to change your inputs, you'll be able to at the end.\n", sep = "")

message_background_ordered <- paste("We will need a few parameters before we can obtain the simulation\n",
                                   "    Window size - length, width and height (e.g. 100 x 100 x 100)\n",
                                   "    Number of cells (e.g. 10000 cells)\n",
                                   "    Amount of jitter (choose to give a bit or a lot of randomness)\n",
                                   "If you want to change your inputs, you'll be able to at the end.\n", sep = "")

message_mixing <- paste("Would you like to MIX the background cells with chosen cell types randomly?\n")


spaSim3D_background_integrator <- function() {
  
  # Ask if user wants a 'random' or 'ordered' patterned background
  message(message_background)
  user_input_background <- get_integer_input_from_options(c(1, 2))
  
  ### Simulate random pattern
  if (user_input_background == 1) {
    
    # Get required parameters for a random background from user
    message(message_background_random)
    parameter_values <- list("length" = get_positive_numeric_input("length"),
                             "width" = get_positive_numeric_input("width"),
                             "height" = get_positive_numeric_input("height"),
                             "number of cells" = get_positive_numeric_input("number of cells"),
                             "minimum distance between cells" = get_non_negative_numeric_input("minimum distance between cells"))
    display_parameters(parameter_values)
    
    # Generate random background simulation using these parameters
    message("Generating simulation...")
    simulated_spe <- simulate_random_background_cells3D(parameter_values[["number of cells"]],
                                                        parameter_values[["length"]],
                                                        parameter_values[["width"]],
                                                        parameter_values[["height"]],
                                                        parameter_values[["minimum distance between cells"]])
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["length"]] <- get_positive_numeric_input("length")
      if (user_input_parameter_choice == 2) parameter_values[["width"]] <- get_positive_numeric_input("width")
      if (user_input_parameter_choice == 3) parameter_values[["height"]] <- get_positive_numeric_input("height")
      if (user_input_parameter_choice == 4) parameter_values[["number of cells"]] <- get_positive_numeric_input("number of cells")
      if (user_input_parameter_choice == 5) parameter_values[["minimum distance between cells"]] <- get_non_negative_numeric_input("minimum distance between cells")
      
      # Generate random background simulation using updated parameters
      display_parameters(parameter_values)
      message("Generating simulation...")
      simulated_spe <- simulate_random_background_cells3D(parameter_values[["number of cells"]],
                                                          parameter_values[["length"]],
                                                          parameter_values[["width"]],
                                                          parameter_values[["height"]],
                                                          parameter_values[["minimum distance between cells"]])
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Simulate ordered pattern
  else if (user_input_background == 2) {
    
    # Get required parameters for a ordered background from user
    message(message_background_ordered)
    parameter_values <- list("length" = get_positive_numeric_input("length"),
                             "width" = get_positive_numeric_input("width"),
                             "height" = get_positive_numeric_input("height"),
                             "number of cells" = get_positive_numeric_input("number of cells"),
                             "amount of jitter" = get_numeric_between_input("amount of jitter", 0, 1))
    display_parameters(parameter_values)
    
    # Generate ordered background simulation using these parameters
    message("Generating simulation...")
    simulated_spe <- simulate_ordered_background_cells3D(parameter_values[["number of cells"]],
                                                        parameter_values[["length"]],
                                                        parameter_values[["width"]],
                                                        parameter_values[["height"]],
                                                        parameter_values[["amount of jitter"]])
    
    # Allow user the option to change their input parameters
    message("Would you like to change your inputs?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values))) # 5 different parameters
      
      if (user_input_parameter_choice == 1) parameter_values[["length"]] <- get_positive_numeric_input("length")
      if (user_input_parameter_choice == 2) parameter_values[["width"]] <- get_positive_numeric_input("width")
      if (user_input_parameter_choice == 3) parameter_values[["height"]] <- get_positive_numeric_input("height")
      if (user_input_parameter_choice == 4) parameter_values[["number of cells"]] <- get_positive_numeric_input("number of cells")
      if (user_input_parameter_choice == 5) parameter_values[["amount of jitter"]] <- get_numeric_between_input("amount of jitter", 0, 1)
      
      # Generate ordered background simulation using updated parameters
      display_parameters(parameter_values)
      message("Generating simulation...")
      simulated_spe <- simulate_ordered_background_cells3D(parameter_values[["number of cells"]],
                                                          parameter_values[["length"]],
                                                          parameter_values[["width"]],
                                                          parameter_values[["height"]],
                                                          parameter_values[["amount of jitter"]])
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  
  ### Simulate mixing
  message(message_mixing)
  choose_cell_types_y_or_n <- get_y_or_n_input()
  if (choose_cell_types_y_or_n == "y") {
    simulated_spe <- get_cell_types_and_proportions_for_mixing(simulated_spe) 
  }
  message("All done!")
  
  return(simulated_spe)
}




### Message cluster integrator strings and function ------------------------------------------------------------
message_no_simulated_spe <- paste("Hello spaSim-3D user. Please input your simulated spe object into this function.\n",
                                  "If you don't have any, you can use the spaSim3D_background_integrator function")

message_shape_choice <- paste("Hello spaSim-3D user, hopefully you can see a plot of your current spe object. What type of shape do you want your cluster to be?\n
          1. Sphere\n
          2. Ellipsoid\n
          3. Cylinder\n
          4. Network\n\n",
                              "To choose, please enter 1, 2, 3 or 4.\n", sep = "")

message_sphere_cluster <- paste("We will need a few parameters to generate a sphere cluster\n",
                                "    Radius\n",
                                "    Coordinates of sphere centre: x, y and z\n",
                                "If you want to change your inputs, you'll be able to at the end.\n", sep = "")

message_ellipsoid_cluster <- paste("We will need a few parameters to generate a ellipsoid cluster\n",
                                   "    Radii: x, y and z\n",
                                   "    Coordinates of ellipsoid centre: x, y and z\n",
                                   "    Angle of rotation in the x-axis, y-axis and z-axis\n",
                                   "If you want to change your inputs, you'll be able to at the end.\n", sep = "")

message_cylinder_cluster <- paste("We will need a few parameters to generate a cylinder cluster\n",
                                  "    Radius\n",
                                  "    Coordinates of the cylinder start point: x, y and z\n",
                                  "    Coordinates of the cylinder end point: x, y and z\n",
                                  "If you want to change your inputs, you'll be able to at the end.\n", sep = "")

message_network_cluster <- paste("We will need a few parameters to generate a network cluster\n",
                                 "    Number of branches\n",
                                 "    Width of each branch: x, y and z\n",
                                 "    Radius spanned by the whole network: x, y and z\n",
                                 "    Coordinates of network centre: x, y and z\n",
                                 "If you want to change your inputs, you'll be able to at the end.\n", sep = "")


message_cluster_choice <- paste("You can customise your cluster further if you'd like:\n
          1. Add a ring\n
          2. Add a double ring\n
          3. Continue\n\n",
                                "To choose, please enter 1, 2 or 3.\n", sep = "")


spaSim3D_cluster_integrator <- function(simulated_spe = NULL) {
  
  ## Start with checking if the user has inputted spe object
  if (class(simulated_spe) != "SpatialExperiment") {
    stop(message_no_simulated_spe)
  }
  
  ## Plot the user's data so they can see what they already have
  fig <- plot_cells3D(simulated_spe)
  print(fig)
  
  ## Get user's choice for shape type (sphere, ellipsoid, cylinder or network)
  message(message_shape_choice)
  user_input_shape <- get_integer_input_from_options(1:4)
  
  ### Sphere
  if (user_input_shape == 1) {
    # Get required parameters for a sphere cluster from user
    message(message_sphere_cluster)
    parameter_values <- list("radius" = get_positive_numeric_input("radius"),
                             "centre x coordinate" = get_non_negative_numeric_input("centre x coordinate"),
                             "centre y coordinate" = get_non_negative_numeric_input("centre y coordinate"),
                             "centre z coordinate" = get_non_negative_numeric_input("centre z coordinate"))
    display_parameters(parameter_values)
    
    # Generate sphere cluster simulation using these parameters
    cluster_properties <- list(list(shape = "Sphere",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    radius = parameter_values[["radius"]],
                                    centre_loc = c(parameter_values[["centre x coordinate"]],
                                                   parameter_values[["centre y coordinate"]],
                                                   parameter_values[["centre z coordinate"]])))
    message("Generating simulation...")
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
                                             plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["radius"]] <- get_positive_numeric_input("radius")
      if (user_input_parameter_choice == 2) parameter_values[["centre x coordinate"]] <- get_non_negative_numeric_input("centre x coordinate")
      if (user_input_parameter_choice == 3) parameter_values[["centre y coordinate"]] <- get_non_negative_numeric_input("centre y coordinate")
      if (user_input_parameter_choice == 4) parameter_values[["centre z coordinate"]] <- get_non_negative_numeric_input("centre z coordinate")
      
      display_parameters(parameter_values)
      
      # Generate sphere cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "Sphere",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      radius = parameter_values[["radius"]],
                                      centre_loc = c(parameter_values[["centre x coordinate"]],
                                                     parameter_values[["centre y coordinate"]],
                                                     parameter_values[["centre z coordinate"]])))
      message("Generating simulation...")
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
                                               plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Ellipsoid
  else if (user_input_shape == 2) {
    # Get required parameters for an ellipsoid cluster from user
    message(message_ellipsoid_cluster)
    parameter_values <- list("x radius" = get_positive_numeric_input("x radius"),
                             "y radius" = get_positive_numeric_input("y radius"),
                             "z radius" = get_positive_numeric_input("z radius"),
                             "centre x coordinate" = get_non_negative_numeric_input("centre x coordinate"),
                             "centre y coordinate" = get_non_negative_numeric_input("centre y coordinate"),
                             "centre z coordinate" = get_non_negative_numeric_input("centre z coordinate"),
                             "x-axis rotation angle" = get_non_negative_numeric_input("x-axis rotation angle"),
                             "y-axis rotation angle" = get_non_negative_numeric_input("y-axis rotation angle"),
                             "z-axis rotation angle" = get_non_negative_numeric_input("z-axis rotation angle"))
    display_parameters(parameter_values)
    
    # Generate ellipsoid cluster simulation using these parameters
    cluster_properties <- list(list(shape = "Ellipsoid",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    x_radius = parameter_values[["x radius"]],
                                    y_radius = parameter_values[["y radius"]],
                                    z_radius = parameter_values[["z radius"]],
                                    centre_loc = c(parameter_values[["centre x coordinate"]],
                                                   parameter_values[["centre y coordinate"]],
                                                   parameter_values[["centre z coordinate"]]),
                                    y_z_rotation = parameter_values[["x-axis rotation angle"]],
                                    x_z_rotation = parameter_values[["y-axis rotation angle"]],
                                    x_y_rotation = parameter_values[["z-axis rotation angle"]]))
    
    message("Generating simulation...")
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
                                             plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["x radius"]] <- get_positive_numeric_input("x radius")
      if (user_input_parameter_choice == 2) parameter_values[["y radius"]] <- get_positive_numeric_input("y radius")
      if (user_input_parameter_choice == 3) parameter_values[["z radius"]] <- get_positive_numeric_input("z radius")
      if (user_input_parameter_choice == 4) parameter_values[["centre x coordinate"]] <- get_non_negative_numeric_input("centre x coordinate")
      if (user_input_parameter_choice == 5) parameter_values[["centre y coordinate"]] <- get_non_negative_numeric_input("centre y coordinate")
      if (user_input_parameter_choice == 6) parameter_values[["centre z coordinate"]] <- get_non_negative_numeric_input("centre z coordinate")
      if (user_input_parameter_choice == 7) parameter_values[["x-axis rotation angle"]] <- get_non_negative_numeric_input("x-axis rotation angle")
      if (user_input_parameter_choice == 8) parameter_values[["y-axis rotation angle"]] <- get_non_negative_numeric_input("y-axis rotation angle")
      if (user_input_parameter_choice == 9) parameter_values[["z-axis rotation angle"]] <- get_non_negative_numeric_input("z-axis rotation angle")
      
      display_parameters(parameter_values)
      
      # Generate ellipsoid cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "Ellipsoid",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      x_radius = parameter_values[["x radius"]],
                                      y_radius = parameter_values[["y radius"]],
                                      z_radius = parameter_values[["z radius"]],
                                      centre_loc = c(parameter_values[["centre x coordinate"]],
                                                     parameter_values[["centre y coordinate"]],
                                                     parameter_values[["centre z coordinate"]]),
                                      y_z_rotation = parameter_values[["x-axis rotation angle"]],
                                      x_z_rotation = parameter_values[["y-axis rotation angle"]],
                                      x_y_rotation = parameter_values[["z-axis rotation angle"]]))
      
      message("Generating simulation...")
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
                                               plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Cylinder
  else if (user_input_shape == 3) {
    # Get required parameters for a cylinder cluster from user
    message(message_cylinder_cluster)
    parameter_values <- list("radius" = get_positive_numeric_input("radius"),
                             "start x coordinate" = get_non_negative_numeric_input("start x coordinate"),
                             "start y coordinate" = get_non_negative_numeric_input("start y coordinate"),
                             "start z coordinate" = get_non_negative_numeric_input("start z coordinate"),
                             "end x coordinate" = get_non_negative_numeric_input("end x coordinate"),
                             "end y coordinate" = get_non_negative_numeric_input("end y coordinate"),
                             "end z coordinate" = get_non_negative_numeric_input("end z coordinate"))
    display_parameters(parameter_values)
    
    # Generate cylinder cluster simulation using these parameters
    cluster_properties <- list(list(shape = "Cylinder",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    radius = parameter_values[["radius"]],
                                    start_loc = c(parameter_values[["start x coordinate"]],
                                                  parameter_values[["start y coordinate"]],
                                                  parameter_values[["start z coordinate"]]),
                                    end_loc = c(parameter_values[["end x coordinate"]],
                                                parameter_values[["end y coordinate"]],
                                                parameter_values[["end z coordinate"]])))
    message("Generating simulation...")
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
                                             plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["radius"]] <- get_positive_numeric_input("radius")
      if (user_input_parameter_choice == 2) parameter_values[["start x coordinate"]] <- get_non_negative_numeric_input("start x coordinate")
      if (user_input_parameter_choice == 3) parameter_values[["start y coordinate"]] <- get_non_negative_numeric_input("start y coordinate")
      if (user_input_parameter_choice == 4) parameter_values[["start z coordinate"]] <- get_non_negative_numeric_input("start z coordinate")
      if (user_input_parameter_choice == 5) parameter_values[["end x coordinate"]] <- get_non_negative_numeric_input("end x coordinate")
      if (user_input_parameter_choice == 6) parameter_values[["end y coordinate"]] <- get_non_negative_numeric_input("end y coordinate")
      if (user_input_parameter_choice == 7) parameter_values[["end z coordinate"]] <- get_non_negative_numeric_input("end z coordinate")
      
      display_parameters(parameter_values)
      
      # Generate cylinder cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "Cylinder",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      radius = parameter_values[["radius"]],
                                      start_loc = c(parameter_values[["start x coordinate"]],
                                                    parameter_values[["start y coordinate"]],
                                                    parameter_values[["start z coordinate"]]),
                                      end_loc = c(parameter_values[["end x coordinate"]],
                                                  parameter_values[["end y coordinate"]],
                                                  parameter_values[["end z coordinate"]])))
      message("Generating simulation...")
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
                                               plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Network
  else if (user_input_shape == 4) {
    # Get required parameters for a network cluster from user
    message(message_network_cluster)
    parameter_values <- list("number of branches" = get_integer_greater_than_or_equal_input("number of branches", 2),
                             "width of branch" = get_positive_numeric_input("width of branch"),
                             "radius spanned by network" = get_positive_numeric_input("radius spanned by network"),
                             "centre x coordinate" = get_non_negative_numeric_input("centre x coordinate"),
                             "centre y coordinate" = get_non_negative_numeric_input("centre y coordinate"),
                             "centre z coordinate" = get_non_negative_numeric_input("centre z coordinate"))
    display_parameters(parameter_values)
    
    # Generate network cluster simulation using these parameters
    cluster_properties <- list(list(shape = "Network",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    n_edges = parameter_values[["number of branches"]],
                                    width = parameter_values[["width of branch"]],
                                    radius = parameter_values[["radius spanned by network"]],
                                    centre_loc = c(parameter_values[["centre x coordinate"]],
                                                   parameter_values[["centre y coordinate"]],
                                                   parameter_values[["centre z coordinate"]])))
    message("Generating simulation...")
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
                                             plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["number of branches"]] <- get_integer_greater_than_or_equal_input("number of branches", 2)
      if (user_input_parameter_choice == 2) parameter_values[["width of branch"]] <- get_positive_numeric_input("width of branch")
      if (user_input_parameter_choice == 3) parameter_values[["radius spanned by network"]] <- get_positive_numeric_input("radius spanned by network")
      if (user_input_parameter_choice == 4) parameter_values[["centre x coordinate"]] <- get_non_negative_numeric_input("centre x coordinate")
      if (user_input_parameter_choice == 5) parameter_values[["centre y coordinate"]] <- get_non_negative_numeric_input("centre y coordinate")
      if (user_input_parameter_choice == 6) parameter_values[["centre z coordinate"]] <- get_non_negative_numeric_input("centre z coordinate")
      
      display_parameters(parameter_values)
      
      # Generate sphere cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "Network",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      n_edges = parameter_values[["number of branches"]],
                                      width = parameter_values[["width of branch"]],
                                      radius = parameter_values[["radius spanned by network"]],
                                      centre_loc = c(parameter_values[["centre x coordinate"]],
                                                     parameter_values[["centre y coordinate"]],
                                                     parameter_values[["centre z coordinate"]])))
      
      message("Generating simulation...")
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
                                               plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  
  # Allow user to change the cell composition of the cluster
  message("Let's change the cell composition of this cluster")
  simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                  simulate_clusters3D,
                                                                                  cluster_properties,
                                                                                  "cluster_cell_types",
                                                                                  "cluster_cell_proportions",
                                                                                  "Cluster")
  simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
  cluster_properties <- simulated_spe_new_and_properties[["properties"]]
  
  
  
  ## Get user's choice for cluster type (ringed, double ringed or continue)
  message(message_cluster_choice)
  user_input_cluster <- get_integer_input_from_options(1:3)
  
  ### Ring
  if (user_input_cluster == 1) {
    # Get width of ring from user
    message("For a single ring, we needs its width.\n")
    
    # Generate cluster with ring simulation using this width
    cluster_properties[[1]][["ring_width"]] <- get_positive_numeric_input("ring width")
    cluster_properties[[1]][["ring_cell_types"]] <- c("Ring")
    cluster_properties[[1]][["ring_cell_proportions"]] <- 1
    
    message("Generating simulation...")
    simulated_spe_new <- simulate_rings3D(simulated_spe,
                                          cluster_properties,
                                          plot_image = TRUE,
                                          plot_cell_types = NULL,
                                          plot_colours = NULL)
    
    # Allow user the option to change the ring width
    message("Would you like to change the ring width?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      cluster_properties[[1]][["ring_width"]] <- get_positive_numeric_input("ring width")
      
      message("Generating simulation...")
      simulated_spe_new <- simulate_rings3D(simulated_spe,
                                            cluster_properties,
                                            plot_image = TRUE,
                                            plot_cell_types = NULL,
                                            plot_colours = NULL)
      
      message("Would you like to the ring width?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
    
    # Allow user to change the cell composition of the ring
    message("Let's change the cell composition of the ring")
    
    simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                    simulate_rings3D,
                                                                                    cluster_properties,
                                                                                    "ring_cell_types",
                                                                                    "ring_cell_proportions",
                                                                                    "Ring")
    
    simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
  }
  ### Double ring
  else if (user_input_cluster == 2) {
    # Get width of inner and outer ring from user
    message("For a double ring, we needs the width of the inner and outer ring.\n")
    
    # Generate cluster with double ring simulation using both widths
    parameter_values <- list("inner ring width" = get_positive_numeric_input("inner ring width"),
                             "outer ring width" = get_positive_numeric_input("outer ring width"))
    display_parameters(parameter_values)
    
    cluster_properties[[1]][["inner_ring_width"]] <- parameter_values[["inner ring width"]]
    cluster_properties[[1]][["outer_ring_width"]] <- parameter_values[["outer ring width"]]
    cluster_properties[[1]][["inner_ring_cell_types"]] <- c("Inner ring")
    cluster_properties[[1]][["inner_ring_cell_proportions"]] <- 1
    cluster_properties[[1]][["outer_ring_cell_types"]] <- c("Outer ring")
    cluster_properties[[1]][["outer_ring_cell_proportions"]] <- 1
    
    message("Generating simulation...")
    simulated_spe_new <- simulate_double_rings3D(simulated_spe,
                                                 cluster_properties,
                                                 plot_image = TRUE,
                                                 plot_cell_types = NULL,
                                                 plot_colours = NULL)
    
    # Allow user the option to change the widths of the inner or outer ring
    message("Would you like to change the widths of the inner or outer ring?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      if (user_input_parameter_choice == 1) parameter_values[["inner ring width"]] <- get_positive_numeric_input("inner ring width")
      if (user_input_parameter_choice == 2) parameter_values[["outer ring width"]] <- get_positive_numeric_input("outer ring width")
      
      cluster_properties[[1]][["inner_ring_width"]] <- parameter_values[["inner ring width"]]
      cluster_properties[[1]][["outer_ring_width"]] <- parameter_values[["outer ring width"]]
      
      display_parameters(parameter_values)
      
      message("Generating simulation...")
      simulated_spe_new <- simulate_double_rings3D(simulated_spe,
                                                   cluster_properties,
                                                   plot_image = TRUE,
                                                   plot_cell_types = NULL,
                                                   plot_colours = NULL)
      
      message("Would you like to change the widths of the inner or outer ring?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
    
    # Allow user to change the cell composition of the inner ring
    message("Let's change the cell composition of the inner ring")
    simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                    simulate_double_rings3D,
                                                                                    cluster_properties,
                                                                                    "inner_ring_cell_types",
                                                                                    "inner_ring_cell_proportions",
                                                                                    "Inner ring")
    
    simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
    cluster_properties <- simulated_spe_new_and_properties[["properties"]]
    
    # Allow user to change the cell composition of the outer ring
    message("Let's change the cell composition of the outer ring")
    simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                    simulate_double_rings3D,
                                                                                    cluster_properties,
                                                                                    "outer_ring_cell_types",
                                                                                    "outer_ring_cell_proportions",
                                                                                    "Outer ring")
    
    simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
  }
  ### Continue
  else if (user_input_cluster == 3) {
    
  }
  
  message("All done!")
  return(simulated_spe_new) 
}


### Utility functions -----------------------------------------------------------------

get_y_or_n_input <- function() {
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = "Enter either y or n: ")
    
    if (user_input %in% c("y", "n")) {
      valid_input <- TRUE
      message("Valid input received!")
    }
    else {
      message("Invalid input. Please enter either y or n.")
    }
  }
  
  return(user_input)
}


get_random_cell_type <- function(cell_types, cell_proportions) {
  
  random <- runif(n = 1, min = 0, max = 1)
  i <- 1
  current_proportion <- 0
  
  while (i <= length(cell_types)){
    current_proportion <- current_proportion + cell_proportions[i]
    if (random <= current_proportion) {
      return(cell_types[i])
    }
    i <- i + 1
  }
}


get_positive_numeric_input <- function(parameter) {
  
  prompt <- paste("Enter a positive numeric value for the ", parameter, ": ", sep = "")
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = prompt)
    # Try converting to numeric
    positive_numeric_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
    
    # Non-numeric input
    if (is.na(positive_numeric_value)) {
      message("Invalid input. Please enter a numeric value.")
    }
    # Non-positive input
    else if (positive_numeric_value <= 0) {
      message("Non-positive input. Please enter a positive number") 
    }
    # Should be correct input
    else {
      valid_input <- TRUE
      message("Valid input received!") 
    }
  }
  
  return(positive_numeric_value)
}


get_numeric_between_input <- function(parameter, lower, upper) {
  
  prompt <- paste("Enter a numeric value between ", lower, " and ", upper, " for the ", parameter, ": ", sep = "")
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = prompt)
    # Try converting to numeric
    numeric_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
    
    # Non-numeric input
    if (is.na(numeric_value)) {
      message("Invalid input. Please enter a numeric value.")
    }
    # Out of bounds input
    else if (numeric_value < lower || numeric_value > upper) {
      message("Out of bounds input. Please a number between ", lower, " and ", upper, ".", sep = "")
    }
    # Should be correct
    else {
      valid_input <- TRUE
      message("Valid input received!")
    }
  }
  
  return(numeric_value)
}


get_non_negative_numeric_input <- function(parameter) {
  
  prompt <- paste("Enter a non-negative numeric value for the ", parameter, ": ", sep = "")
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = prompt)
    # Try converting to numeric
    non_negative_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
    
    # Non-numeric input
    if (is.na(non_negative_value)) {
      message("Invalid input. Please enter a numeric value.")
    }
    # Negative input
    else if (non_negative_value < 0) {
      message("Negative input. Please enter a non-negative number") 
    }
    # Should be correct input
    else {
      valid_input <- TRUE
      message("Valid input received!") 
    }
  }
  
  return(non_negative_value)
}



get_integer_input_from_options <- function(integer_options) {
  
  first_integers <- integer_options[1:(length(integer_options) - 1)]
  last_integer <- integer_options[length(integer_options)]
  integers_string <- paste(paste(first_integers, collapse = ", "), "or", last_integer)
  
  prompt <- paste("Enter either ", integers_string, ": ", sep = "")
  invalid_input_message <- paste("Invalid input. Please enter only", integers_string)
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = prompt)
    # Try converting to integer
    int_value <- tryCatch({as.integer(user_input)}, error = function(e) NA)
    
    # Check if conversion was successful and value is in integer_options
    if (!is.na(int_value) && int_value %in% integer_options) {
      valid_input <- TRUE
      message("Valid input received!")
    } 
    else {
      message(invalid_input_message)
    }
  }
  
  return(int_value)
}



get_integer_greater_than_or_equal_input <- function(parameter, lower) {
  
  prompt <- paste("Enter an integer in value greater than or equal to ", lower, " for the ", parameter, ": ", sep = "")
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = prompt)
    # Try converting to numeric
    integer_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
    
    # Non-numeric input
    if (is.na(integer_value)) {
      message("Invalid input. Please enter a numeric integer value.")
    }
    # Numeric but not integer
    else if (integer_value%%1 != 0) {
      message("Non-integer input. Please enter an integer value.")
    }
    # Integer but below the lower bound
    else if (integer_value < lower) {
      message("Out of bounds input. Please a number greater than or equal to ", lower)
    }
    else {
      valid_input <- TRUE
      message("Valid input received!")
    }
  }
  
  return(integer_value)
}



message_get_cell_types <- "Keep entering the name of cell types you would like (e.g. Tumour, Immune, etc.).\n    enter 'stop' to move on."



get_cell_types_and_proportions_for_mixing <- function(simulated_spe) {
  
  ## Get cell types from user
  cell_types <- c()
  user_input <- ""
  message(message_get_cell_types)
  while (user_input != "stop") {
    
    user_input <- readline(prompt = "Enter a cell type, or enter 'stop': ")
    
    ## Ignore if user enters a blank string
    if (user_input == "") {
      
    }
    ## Add inputted cell type to cell_types vector
    else if (user_input != "stop") {
      cell_types <- c(cell_types, user_input)
      message(paste("Cell type added:", user_input))
    }
    ## User wants to stop but hasn't entered any cell types
    else if (user_input == "stop" && length(cell_types) == 0) {
      message("You have not entered any cell types. Try again\n")
      user_input <- ""
    }
    ## User wants to stop
    else {
      message(paste("Your cell types chosen are:", paste(cell_types, collapse = ", ")))
      
      ## Allow user to re-choose cell types
      message("Would like to re-choose these cell types?\n")
      user_input_y_or_n <- get_y_or_n_input()
      if (user_input_y_or_n == "y") {
        cell_types <- c()
        message(message_get_cell_types)
        user_input <- ""
      }
    }
  }
  
  ## Get cell proportions from user
  cell_proportions <- c()
  max_proportion <- 1
  i <- 1
  message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
  while (i <= length(cell_types)) {
    
    ## For the last cell type, we can figure out what the cell proportion must be
    if (i == length(cell_types)) {
      cell_proportions <- c(cell_proportions, max_proportion)
      message("Cell proportion for ", cell_types[i], " must be ", round(max_proportion, 5))
    }
    ## Add inputted cell proportion to cell_proportions vector
    else {
      cell_proportion <- get_numeric_between_input(paste("cell proportion of", cell_types[i], "cells"), 0, max_proportion)
      cell_proportions <- c(cell_proportions, cell_proportion)
      max_proportion <- 1 - sum(cell_proportions)
      message("Cell proportion for ", cell_types[i], " is ", cell_proportion)
    }
    i <- i + 1
    
    if (i > length(cell_types)) {
      ## Generate simulation
      message("Generating simulation...")
      simulated_spe <- simulate_mixing3D(simulated_spe,
                                         cell_types,
                                         cell_proportions,
                                         plot_image = F)
      
      fig <- plot_cells3D(simulated_spe)
      print(fig)
      
      if (length(cell_types) == 1) break # If there is only one cell type, proportion is always 1
      
      ## Allow user to re-choose cell proportions  
      message("Would like to re-choose these cell proportions?\n")
      user_input_y_or_n <- get_y_or_n_input()
      if (user_input_y_or_n == "y") {
        cell_proportions <- c()
        max_proportion <- 1
        i <- 1
        message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
      }
    }
  }
  
  return(simulated_spe)
}



get_cell_types_and_proportions_for_clusters <- function(simulated_spe, simulate_function, properties, cell_type_option, cell_proportion_option, temp_cell_type) {
  
  ## Display the cell types currently found in simulated_spe to the user
  current_cell_types <- setdiff(unique(simulated_spe[["Cell.Type"]]), temp_cell_type)
  message("Your data currently has the following cell types:\n", paste(current_cell_types, collapse = ", "), "\n")
  
  ## Get cell types from user
  cell_types <- c()
  user_input <- ""
  message(message_get_cell_types)
  while (user_input != "stop") {
    
    user_input <- readline(prompt = "Enter a cell type, or enter 'stop': ")
    
    ## Ignore if user enters a blank string
    if (user_input == "") {
      
    }
    ## Add inputted cell type to cell_types vector
    else if (user_input != "stop") {
      cell_types <- c(cell_types, user_input)
      message(paste("Cell type added:", user_input))
    }
    ## User wants to stop but hasn't entered any cell types
    else if (user_input == "stop" && length(cell_types) == 0) {
      message("You have not entered any cell types. Try again\n")
      user_input <- ""
    }
    ## User wants to stop
    else {
      message(paste("Your cell types chosen are:", paste(cell_types, collapse = ", ")))
      
      ## Allow user to re-choose cell types
      message("Would like to re-choose these cell types?\n")
      user_input_y_or_n <- get_y_or_n_input()
      if (user_input_y_or_n == "y") {
        message("Your data currently has the following cell types:\n", paste(current_cell_types, collapse = ", "), "\n")
        cell_types <- c()
        message(message_get_cell_types)
        user_input <- ""
      }
    }
  }
  properties[[1]][[cell_type_option]] <- cell_types
  
  ## Get cell proportions from user
  cell_proportions <- c()
  max_proportion <- 1
  i <- 1
  message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
  while (i <= length(cell_types)) {
    
    ## For the last cell type, we can figure out what the cell proportion must be
    if (i == length(cell_types)) {
      cell_proportions <- c(cell_proportions, max_proportion)
      message("Cell proportion for ", cell_types[i], " must be ", round(max_proportion, 5))
    }
    ## Add inputted cell proportion to cell_proportions vector
    else {
      cell_proportion <- get_numeric_between_input(paste("cell proportion of", cell_types[i], "cells"), 0, max_proportion)
      cell_proportions <- c(cell_proportions, cell_proportion)
      max_proportion <- 1 - sum(cell_proportions)
      message("Cell proportion for ", cell_types[i], " is ", cell_proportion)
    }
    i <- i + 1
    
    if (i > length(cell_types)) {
      properties[[1]][[cell_proportion_option]] <- cell_proportions
      
      ## Convert spe object to data frame
      df <- data.frame(spatialCoords(simulated_spe), "Cell.Type" = simulated_spe[["Cell.Type"]])
      
      ## Just change the cell type of the temp_cell_type, no need to actually re-simulate
      for (i in seq(nrow(df))) {
        if (df[i, "Cell.Type"] == temp_cell_type) {
          df[i, "Cell.Type"] <- get_random_cell_type(cell_types, cell_proportions)
        }
      }
      
      # Add Cell.ID column to data frame
      df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")
      
      # Update current meta data
      metadata <- simulated_spe@metadata
      metadata[["simulation"]][[length(metadata[["simulation"]])]][[cell_type_option]] <- cell_types
      metadata[["simulation"]][[length(metadata[["simulation"]])]][[cell_proportion_option]] <- cell_proportions
      
      # Convert data frame to spe object
      simulated_spe_new <- SpatialExperiment(
        assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
        colData = df,
        spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
        metadata = metadata)
      
      ## Generate simulation
      message("Generating simulation...")
      fig <- plot_cells3D(simulated_spe_new)
      print(fig)
      
      if (length(cell_types) == 1) break # If there is only one cell type, proportion is always 1
      
      ## Allow user to re-choose cell proportions  
      message("Would like to re-choose these cell proportions?\n")
      user_input_y_or_n <- get_y_or_n_input()
      if (user_input_y_or_n == "y") {
        cell_proportions <- c()
        max_proportion <- 1
        i <- 1
        message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
      }
    }
  }
  
  return(list(data = simulated_spe_new, properties = properties))
}



display_parameters <- function(parameter_values) {
  
  message("Your current inputs are:\n")
  
  display_message <- ""
  
  for (i in seq(length(parameter_values))) {
    display_message <- paste(display_message, "    ", i, ". ", names(parameter_values)[i], ": ", parameter_values[[i]], '\n', sep = "")
  }
  message(display_message)
}


### spe_metadata functions -----------------------------------------------------------------

spe_metadata_background_template <- function(background_type) {
  
  if (background_type == "random") {
    background_metadata <- list(background = list(background_type = "random",
                                                  n_cells = 10000,
                                                  length = 100,
                                                  width = 100,
                                                  height = 100,
                                                  minimum_distance_between_cells = 2,
                                                  cell_types = c("Tumour", "Others"),
                                                  cell_proportions = c(0.05, 0.95)))
  }
  else if (background_type == "ordered") {
    background_metadata <- list(background = list(background_type = "ordered",
                                                  n_cells = 10000,
                                                  length = 100,
                                                  width = 100,
                                                  height = 100,
                                                  jitter_proportion = 0.25,
                                                  cell_types = c("Immune", "Others"),
                                                  cell_proportions = c(0.05, 0.95)))
  }
  else {
    stop("background_type parameter must be 'random' or 'ordered'.")
  }
  
  return(background_metadata)
}



spe_metadata_cluster_template <- function(background_metadata, cluster_type, shape) {
  
  
  ### Get template for different shapes
  if (shape == "Sphere") {
    cluster_metadata <- list(shape = "Sphere",
                             cluster_cell_types = c("Tumour", "Immune", "Others"),
                             cluster_cell_proportions = c(0.8, 0.15, 0.05),
                             radius = 25,
                             centre_loc = c(40, 40, 40))
  }
  else if (shape == "Ellipsoid") {
    cluster_metadata <- list(shape = "Ellipsoid",
                             cluster_cell_types = c("Tumour", "Immune", "Others"),
                             cluster_cell_proportions = c(0.8, 0.15, 0.05),
                             x_radius = 15,
                             y_radius = 20,
                             z_radius = 25,
                             centre_loc = c(70, 70, 70),
                             x_y_rotation = 0,
                             x_z_rotation = 45,
                             y_z_rotation = 0)
  }
  else if (shape == "Cylinder") {
    cluster_metadata <- list(shape = "Cylinder",
                             cluster_cell_types = c("Endothelial", "Others"),
                             cluster_cell_proportions = c(0.95, 0.05),
                             radius = 10,
                             start_loc = c(0, 0, 0),
                             end_loc   = c(20, 20, 100)) 
  }
  else if (shape == "Network") {
    cluster_metadata <- list(shape = "Network",
                             cluster_cell_types = c("Tumour", "Others"),
                             cluster_cell_proportions = c(0.95, 0.05),
                             n_edges = 18,
                             width = 9,
                             centre_loc = c(50, 50, 50),
                             radius = 50)
  }
  else {
    stop("shape parameter must be 'Sphere', 'Ellipsoid', 'Cylinder' or 'Network'")
  }
  
  ### Add extra metadata for different cluster types
  if (cluster_type == "regular") {
    cluster_metadata <- append(list(cluster_type = "regular"), cluster_metadata)    
  }
  else if (cluster_type == "ring") {
    cluster_metadata <- append(list(cluster_type = "ring"), cluster_metadata)
    cluster_metadata$ring_cell_types <- c("Immune", "Others")
    cluster_metadata$ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$ring_width <- 5
  }
  else if (cluster_type == "double ring") {
    cluster_metadata <- append(list(cluster_type = "double ring"), cluster_metadata)
    cluster_metadata$inner_ring_cell_types <- c("Immune1", "Others")
    cluster_metadata$inner_ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$inner_ring_width <- 3
    cluster_metadata$outer_ring_cell_types <- c("Immune2", "Others")
    cluster_metadata$outer_ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$outer_ring_width <- 2
  }
  else {
    stop("cluster_type parameter must be 'regular', 'ring' or 'double ring'")
  }
  
  background_metadata[[paste("cluster", length(background_metadata), sep="_")]] <- cluster_metadata
  
  return(background_metadata)
}


simulate_spe_metadata3D <- function(spe_metadata, plot_image = TRUE) {
  
  # First element should contain background metadata
  bg_metadata <- spe_metadata[[1]]
  if (bg_metadata$background_type == "random") {
    spe <- simulate_random_background_cells3D(bg_metadata$n_cells,
                                              bg_metadata$length,
                                              bg_metadata$width,
                                              bg_metadata$height,
                                              bg_metadata$minimum_distance_between_cells,
                                              plot_image = plot_image)    
  }
  else if (bg_metadata$background_type == "ordered") {
    spe <- simulate_ordered_background_cells3D(bg_metadata$n_cells,
                                              bg_metadata$length,
                                              bg_metadata$width,
                                              bg_metadata$height,
                                              bg_metadata$jitter_proportion,
                                              plot_image = plot_image) 
  }
  else {
    stop("background_type parameter found in the first list must be 'random' or 'ordered'.")
  }
  # Apply background mixing
  spe <- simulate_mixing3D(spe,
                           bg_metadata$cell_types,
                           bg_metadata$cell_proportions,
                           plot_image = plot_image)
  
  ### If there is only background metadata, we are done
  if (length(spe_metadata) == 1) return(spe)
  
  
  ### All other elements should help to simulate clusters 
  for (i in 2:length(spe_metadata)) {
    cluster_metadata <- spe_metadata[[i]]
    if (cluster_metadata$cluster_type == "regular") {
      spe <- simulate_clusters3D(spe, list(cluster_metadata), plot_image = plot_image)
    }
    else if (cluster_metadata$cluster_type == "ring") {
      spe <- simulate_rings3D(spe, list(cluster_metadata), plot_image = plot_image)      
    }
    else if (cluster_metadata$cluster_type == "double ring") {
      spe <- simulate_double_rings3D(spe, list(cluster_metadata), plot_image = plot_image)
    }
  }
  
  return(spe)
}



add_spe_metadata3D <- function(spe, metadata, plot_image = TRUE) {
  
  # Ignore the 'background' element in metadata
  metadata[['background']] <- NULL
  
  for (i in seq(length(metadata))) {
    metadata_cluster <- metadata[[i]]
    
    if (metadata_cluster$cluster_type == "regular") {
      spe <- simulate_clusters3D(spe, list(metadata_cluster), plot_image = plot_image)
    }
    else if (metadata_cluster$cluster_type == "ring") {
      spe <- simulate_rings3D(spe, list(metadata_cluster), plot_image = plot_image)
    }
    else if (metadata_cluster$cluster_type == "double ring") {
      spe <- simulate_double_rings3D(spe, list(metadata_cluster), plot_image = plot_image)
    }
  }
  
  return(spe)
}










