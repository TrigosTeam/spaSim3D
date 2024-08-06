library(parallel)

cl <- makeCluster(detectCores())

radii <- 1:50

result <- parLapply(cl, radii, function(radius) {
  ### Functions -------------------------------------------------------------------
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(MatrixGenerics)
  library(matrixStats)
  library(plotly)
  library(dplyr)
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
    if (sum(cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
    
    
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
    if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
    
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
    if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
    
    ring_cell_types <- ring_properties$ring_cell_types
    ring_cell_proportions <- ring_properties$ring_cell_proportions
    ring_width <- ring_properties$ring_width
    
    ## Check number of ring cell types matches the number of cell proportions
    if (length(ring_cell_types) != length(ring_cell_proportions)) stop("Number of ring cell types doesn't match number of ring cell proportion.")
    
    ## Check ring cell proportions are not negative or greater than 1
    if (sum(ring_cell_proportions < 0 | ring_cell_proportions > 1) != 0) stop("Ring cell proportions cannot be negative or greater than 1")
    
    ## Check ring cell proportions add up to 1
    if (sum(ring_cell_proportions) != 1) stop("Sum of ring cell proportions is NOT 1")
    
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
    if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
    
    inner_ring_cell_types <- dr_properties$inner_ring_cell_types
    inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
    inner_ring_width <- dr_properties$inner_ring_width
    
    ## Check number of inner ring cell types matches the number of inner ring cell proportions
    if (length(inner_ring_cell_types) != length(inner_ring_cell_proportions)) stop("Number of inner ring cell types doesn't match number of inner ring cell proportion.")
    
    ## Check inner ring cell proportions are not negative or greater than 1
    if (sum(inner_ring_cell_proportions < 0 | inner_ring_cell_proportions > 1) != 0) stop("Inner ring cell proportions cannot be negative or greater than 1")
    
    ## Check inner ring cell proportions add up to 1
    if (sum(inner_ring_cell_proportions) != 1) stop("Sum of inner ring cell proportions is NOT 1")
    
    outer_ring_cell_types <- dr_properties$outer_ring_cell_types
    outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
    outer_ring_width <- dr_properties$outer_ring_width
    
    ## Check number of outer ring cell types matches the number of outer ring cell proportions
    if (length(outer_ring_cell_types) != length(outer_ring_cell_proportions)) stop("Number of outer ring cell types doesn't match number of outer ring cell proportion.")
    
    ## Check outer ring cell proportions are not negative or greater than 1
    if (sum(outer_ring_cell_proportions < 0 | outer_ring_cell_proportions > 1) != 0) stop("Outer ring cell proportions cannot be negative or greater than 1")
    
    ## Check outer ring cell proportions add up to 1
    if (sum(outer_ring_cell_proportions) != 1) stop("Sum of outer ring cell proportions is NOT 1")
    
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
    if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
    
    
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
    x <- T_M[1, 1] * x + T_M[1, 2] * y + T_M[1, 3] * z
    y <- T_M[2, 1] * x + T_M[2, 2] * y + T_M[2, 3] * z
    z <- T_M[3, 1] * x + T_M[3, 2] * y + T_M[3, 3] * z
    
    bg_spe[["Cell.Type"]] <- ifelse((x / x_radius)^2 +
                                      (y / y_radius)^2 +
                                      (z / z_radius)^2 <= 1,
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
    if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
    
    ring_cell_types <- ring_properties$ring_cell_types
    ring_cell_proportions <- ring_properties$ring_cell_proportions
    ring_width <- ring_properties$ring_width
    
    ## Check number of ring cell types matches the number of cell proportions
    if (length(ring_cell_types) != length(ring_cell_proportions)) stop("Number of ring cell types doesn't match number of ring cell proportion.")
    
    ## Check ring cell proportions are not negative or greater than 1
    if (sum(ring_cell_proportions < 0 | ring_cell_proportions > 1) != 0) stop("Ring cell proportions cannot be negative or greater than 1")
    
    ## Check ring cell proportions add up to 1
    if (sum(ring_cell_proportions) != 1) stop("Sum of ring cell proportions is NOT 1")
    
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
    x <- T_M[1, 1] * x + T_M[1, 2] * y + T_M[1, 3] * z
    y <- T_M[2, 1] * x + T_M[2, 2] * y + T_M[2, 3] * z
    z <- T_M[3, 1] * x + T_M[3, 2] * y + T_M[3, 3] * z
    
    # Start with cells in ring  
    bg_spe[["Cell.Type"]] <- ifelse((x / (x_radius + ring_width))^2 +
                                      (y / (y_radius + ring_width))^2 +
                                      (z / (z_radius + ring_width))^2 <= 1,
                                    sample(ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = ring_cell_proportions),
                                    bg_spe[["Cell.Type"]])
    
    
    # Then do cells in the cluster  
    bg_spe[["Cell.Type"]] <- ifelse((x / x_radius)^2 +
                                      (y / y_radius)^2 +
                                      (z / z_radius)^2 <= 1,
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
    if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
    
    inner_ring_cell_types <- dr_properties$inner_ring_cell_types
    inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
    inner_ring_width <- dr_properties$inner_ring_width
    
    ## Check number of inner ring cell types matches the number of inner ring cell proportions
    if (length(inner_ring_cell_types) != length(inner_ring_cell_proportions)) stop("Number of inner ring cell types doesn't match number of inner ring cell proportion.")
    
    ## Check inner ring cell proportions are not negative or greater than 1
    if (sum(inner_ring_cell_proportions < 0 | inner_ring_cell_proportions > 1) != 0) stop("Inner ring cell proportions cannot be negative or greater than 1")
    
    ## Check inner ring cell proportions add up to 1
    if (sum(inner_ring_cell_proportions) != 1) stop("Sum of inner ring cell proportions is NOT 1")
    
    outer_ring_cell_types <- dr_properties$outer_ring_cell_types
    outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
    outer_ring_width <- dr_properties$outer_ring_width
    
    ## Check number of outer ring cell types matches the number of outer ring cell proportions
    if (length(outer_ring_cell_types) != length(outer_ring_cell_proportions)) stop("Number of outer ring cell types doesn't match number of outer ring cell proportion.")
    
    ## Check outer ring cell proportions are not negative or greater than 1
    if (sum(outer_ring_cell_proportions < 0 | outer_ring_cell_proportions > 1) != 0) stop("Outer ring cell proportions cannot be negative or greater than 1")
    
    ## Check outer ring cell proportions add up to 1
    if (sum(outer_ring_cell_proportions) != 1) stop("Sum of outer ring cell proportions is NOT 1")
    
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
    x <- spe_coords$Cell.X.Position - centre_loc[1]
    y <- spe_coords$Cell.Y.Position - centre_loc[2]
    z <- spe_coords$Cell.Z.Position - centre_loc[3]
    
    # Transform  x, y and z coordinates using rotation transformation matrix
    x <- T_M[1, 1] * x + T_M[1, 2] * y + T_M[1, 3] * z
    y <- T_M[2, 1] * x + T_M[2, 2] * y + T_M[2, 3] * z
    z <- T_M[3, 1] * x + T_M[3, 2] * y + T_M[3, 3] * z
    
    
    # Start with cells in outer ring  
    bg_spe[["Cell.Type"]] <- ifelse((x / (x_radius + inner_ring_width + outer_ring_width))^2 +
                                      (y / (y_radius + inner_ring_width + outer_ring_width))^2 +
                                      (z / (z_radius + inner_ring_width + outer_ring_width))^2 <= 1,
                                    sample(outer_ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = outer_ring_cell_proportions),
                                    bg_spe[["Cell.Type"]])
    
    # Then do cells in inner ring  
    bg_spe[["Cell.Type"]] <- ifelse((x / (x_radius + inner_ring_width))^2 +
                                      (y / (y_radius + inner_ring_width))^2 +
                                      (z / (z_radius + inner_ring_width))^2 <= 1,
                                    sample(inner_ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = inner_ring_cell_proportions),
                                    bg_spe[["Cell.Type"]])
    
    
    # Then do cells in the cluster  
    bg_spe[["Cell.Type"]] <- ifelse((x / x_radius)^2 +
                                      (y / y_radius)^2 +
                                      (z / z_radius)^2 <= 1,
                                    sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                    bg_spe[["Cell.Type"]])
    
    # Update current meta data
    if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
    bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
    
    return(bg_spe)
  }
  
  
  
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
    if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1.")
    
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
    if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
    
    ring_cell_types <- ring_properties$ring_cell_types
    ring_cell_proportions <- ring_properties$ring_cell_proportions
    ring_width <- ring_properties$ring_width
    
    ## Check number of ring cell types matches the number of cell proportions
    if (length(ring_cell_types) != length(ring_cell_proportions)) stop("Number of ring cell types doesn't match number of ring cell proportion.")
    
    ## Check ring cell proportions are not negative or greater than 1
    if (sum(ring_cell_proportions < 0 | ring_cell_proportions > 1) != 0) stop("Ring cell proportions cannot be negative or greater than 1")
    
    ## Check ring cell proportions add up to 1
    if (sum(ring_cell_proportions) != 1) stop("Sum of ring cell proportions is NOT 1")
    
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
    if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
    
    inner_ring_cell_types <- dr_properties$inner_ring_cell_types
    inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
    inner_ring_width <- dr_properties$inner_ring_width
    
    ## Check number of inner ring cell types matches the number of inner ring cell proportions
    if (length(inner_ring_cell_types) != length(inner_ring_cell_proportions)) stop("Number of inner ring cell types doesn't match number of inner ring cell proportion.")
    
    ## Check inner ring cell proportions are not negative or greater than 1
    if (sum(inner_ring_cell_proportions < 0 | inner_ring_cell_proportions > 1) != 0) stop("Inner ring cell proportions cannot be negative or greater than 1")
    
    ## Check inner ring cell proportions add up to 1
    if (sum(inner_ring_cell_proportions) != 1) stop("Sum of inner ring cell proportions is NOT 1")
    
    outer_ring_cell_types <- dr_properties$outer_ring_cell_types
    outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
    outer_ring_width <- dr_properties$outer_ring_width
    
    ## Check number of outer ring cell types matches the number of outer ring cell proportions
    if (length(outer_ring_cell_types) != length(outer_ring_cell_proportions)) stop("Number of outer ring cell types doesn't match number of outer ring cell proportion.")
    
    ## Check outer ring cell proportions are not negative or greater than 1
    if (sum(outer_ring_cell_proportions < 0 | outer_ring_cell_proportions > 1) != 0) stop("Outer ring cell proportions cannot be negative or greater than 1")
    
    ## Check outer ring cell proportions add up to 1
    if (sum(outer_ring_cell_proportions) != 1) stop("Sum of outer ring cell proportions is NOT 1")
    
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
    if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
    
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
    if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
    
    ring_cell_types <- ring_properties$ring_cell_types
    ring_cell_proportions <- ring_properties$ring_cell_proportions
    ring_width <- ring_properties$ring_width
    
    ## Check number of ring cell types matches the number of cell proportions
    if (length(ring_cell_types) != length(ring_cell_proportions)) stop("Number of ring cell types doesn't match number of ring cell proportion.")
    
    ## Check ring cell proportions are not negative or greater than 1
    if (sum(ring_cell_proportions < 0 | ring_cell_proportions > 1) != 0) stop("Ring cell proportions cannot be negative or greater than 1")
    
    ## Check ring cell proportions add up to 1
    if (sum(ring_cell_proportions) != 1) stop("Sum of ring cell proportions is NOT 1")
    
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
    if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
    
    inner_ring_cell_types <- dr_properties$inner_ring_cell_types
    inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
    inner_ring_width <- dr_properties$inner_ring_width
    
    ## Check number of inner ring cell types matches the number of inner ring cell proportions
    if (length(inner_ring_cell_types) != length(inner_ring_cell_proportions)) stop("Number of inner ring cell types doesn't match number of inner ring cell proportion.")
    
    ## Check inner ring cell proportions are not negative or greater than 1
    if (sum(inner_ring_cell_proportions < 0 | inner_ring_cell_proportions > 1) != 0) stop("Inner ring cell proportions cannot be negative or greater than 1")
    
    ## Check inner ring cell proportions add up to 1
    if (sum(inner_ring_cell_proportions) != 1) stop("Sum of inner ring cell proportions is NOT 1")
    
    outer_ring_cell_types <- dr_properties$outer_ring_cell_types
    outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
    outer_ring_width <- dr_properties$outer_ring_width
    
    ## Check number of outer ring cell types matches the number of outer ring cell proportions
    if (length(outer_ring_cell_types) != length(outer_ring_cell_proportions)) stop("Number of outer ring cell types doesn't match number of outer ring cell proportion.")
    
    ## Check outer ring cell proportions are not negative or greater than 1
    if (sum(outer_ring_cell_proportions < 0 | outer_ring_cell_proportions > 1) != 0) stop("Outer ring cell proportions cannot be negative or greater than 1")
    
    ## Check outer ring cell proportions add up to 1
    if (sum(outer_ring_cell_proportions) != 1) stop("Sum of outer ring cell proportions is NOT 1")
    
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
  
  
  ### Analysis ---------------------------------------------------------------

  spe_metadata <- spe_metadata_background_template("random")
  spe_metadata$background$n_cells <- 50000
  spe_metadata$background$minimum_distance_between_cells <- 0
  spe_metadata <- spe_metadata_cluster_template(spe_metadata, "ring", "Network")
  spe_clusters <- simulate_spe_metadata3D(spe_metadata)
  
  calculate_mixing_scores3D <- function(spe, 
                                        reference_cell_types, 
                                        target_cell_types, 
                                        radius, 
                                        feature_colname = "Cell.Type") {
    
    if (is.null(spe[[feature_colname]])) stop(paste("No column called", feature_colname, "found in spe object"))
    
    
    ## For reference_cell_types, check they are found in the spe object
    unknown_cell_types <- setdiff(reference_cell_types, spe[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in reference_cell_types are not found in the spe object:\n   ",
                 paste(unknown_cell_types, collapse = ", ")))
    }
    
    ## For target_cell_types, check they are found in the spe object
    unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
                 paste(unknown_cell_types, collapse = ", ")))
    }
    
    # Check if radius is numeric
    if (!is.numeric(radius)) stop(paste(radius, " is not of type 'numeric'"))
    
    # Get spe coords
    spe_coords <- spatialCoords(spe)
    
    # Define result
    result <- data.frame()
    
    for (reference_cell_type in reference_cell_types) {
      
      # Get coords for reference_cell_type
      reference_cell_type_coords <- spe_coords[spe[[feature_colname]] == reference_cell_type, ]
      
      for (target_cell_type in target_cell_types) {
        
        # Get coords for target_cell_type
        target_cell_type_coords <- spe_coords[spe[[feature_colname]] == target_cell_type, ]
        
        # No point getting mixing scores if comparing the same cell type
        if (reference_cell_type == target_cell_type) {
          next
        }
        
        # Can't get mixing scores if there are no reference cells
        if (nrow(reference_cell_type_coords) == 0) {
          methods::show(paste("There are no unique reference cells of specified cell type ", reference_cell_type, "for target cell", target_cell_type))
          result <-  rbind(result, 
                           c(reference_cell_type, 
                             target_cell_type, 
                             0, 
                             nrow(target_cells), 
                             0, 
                             0, 
                             NA, 
                             NA))
        }
        
        # Can't get mixing scores if there are no target cells
        else if (nrow(target_cell_type_coords) == 0) {
          methods::show(paste("There are no unique target cells of specified cell type", target_cell_type, "for reference cell", reference_cell_type))
          
          ref_ref_result <- dbscan::frNN(reference_cell_type_coords, 
                                         eps = radius, 
                                         query = NULL,
                                         sort = FALSE)
          
          # halve it to avoid counting each ref-ref interaction twice
          n_ref_ref_interactions <- 0.5 * sum(rapply(ref_ref_result$id, length)) 
          
          result <-  rbind(result, 
                           c(reference_cell_type, 
                             target_cell_type, 
                             nrow(reference_cells), 
                             0, 
                             0, 
                             n_ref_ref_interactions, 
                             NA, 
                             NA))
        }
        
        # Generic case: We have reference cells and target cells
        else {
          
          # For each reference cell, find all target cells within the chosen radius
          ref_tar_result <- dbscan::frNN(target_cell_type_coords, 
                                         eps = radius, 
                                         query = reference_cell_type_coords, 
                                         sort = FALSE)
          
          # Find the total sum of how many target cells were close enough to reference cells
          n_ref_tar_interactions <- sum(rapply(ref_tar_result$id, length))
          
          # For each reference cell, find all other reference cells within the chosen radius
          ref_ref_result <- dbscan::frNN(reference_cell_type_coords, 
                                         eps = radius,
                                         query = NULL,
                                         sort = FALSE)
          
          # Find the the total sum of how many other reference cells were close enough to reference cells
          # Halve it to avoid counting each ref-ref interaction twice
          n_ref_ref_interactions <- 0.5 * sum(rapply(ref_ref_result$id, length)) 
          
          
          if (n_ref_ref_interactions != 0) {
            mixing_score <- n_ref_tar_interactions / n_ref_ref_interactions
            normalised_mixing_score <- 0.5 * mixing_score * (nrow(reference_cell_type_coords) - 1) / nrow(target_cell_type_coords)
          }
          else {
            mixing_score <- 0
            normalised_mixing_score <- 0
            methods::show(paste("There are no reference to reference interactions for", target_cell_type, "in the specified radius, cannot calculate mixing score"))
          }
          
          result <-  rbind(result, 
                           c(reference_cell_type, 
                             target_cell_type, 
                             nrow(reference_cell_type_coords), 
                             nrow(target_cell_type_coords), 
                             n_ref_tar_interactions, 
                             n_ref_ref_interactions, 
                             mixing_score, 
                             normalised_mixing_score))
        }
      }
    }
    
    # Required column names of our output data frame
    colnames(result) <- c("ref_cell_type", 
                          "tar_cell_type", 
                          "n_ref_cells",
                          "n_tar_cells", 
                          "n_ref_tar_interactions",
                          "n_ref_ref_interactions", 
                          "mixing_score", 
                          "normalised_mixing_score")
    
    # Turn numeric data into numeric type
    result[ , 3:8] <- apply(result[ , 3:8], 2, as.numeric)
    
    return(result)
  }
  calculate_mixing_scores3D(spe_clusters, "Tumour", "Immune", radius)
})

result <- do.call(rbind, result)

stopCluster(cl)
