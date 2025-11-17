simulate_network_ring <- function(spe, 
                                  ring_properties) {  
  
  # Check input parameters
  input_parameters <- ring_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get network ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  n_edges <- ring_properties$n_edges
  width <- ring_properties$width
  centre_loc <- ring_properties$centre_loc
  radius <- ring_properties$radius
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  # Number of vertices is always one more than the number of edges for the MST will we make
  n_vertices <- n_edges + 1 
  
  ## Generate n_vertices random points with coords inside a sphere with given radius and centre loc. 
  # Starting with 1000 points inside a cube should be a good enough buffer, unless the user wants more than 1000 edges...
  # Lets stop them from inputting more than 99
  max_edges <- 99
  if (n_edges > max_edges) stop("Only networks with less than 100 edges can be simulated")
  random_coords <- data.frame(x = runif(1000, centre_loc[1] - radius, centre_loc[1] + radius),
                              y = runif(1000, centre_loc[2] - radius, centre_loc[2] + radius),
                              z = runif(1000, centre_loc[3] - radius, centre_loc[3] + radius))
  
  # Then subset points which are inside the sphere
  random_coords <- random_coords[(random_coords$x - centre_loc[1])^2 +
                                   (random_coords$y - centre_loc[2])^2 +
                                   (random_coords$z- centre_loc[3])^2 <= radius^2, ]
  
  ## Subset further and pick 'n_vertices' coords to represent the vertices
  random_coords <- sample_n(random_coords, n_vertices)
  
  ## Get adjacency matrix from points (pairwise distance between points)
  # Assume all points have an edge between each other
  # Assume weight of each edge is equal to the distance between points
  adj_mat <- -1 * apcluster::negDistMat(random_coords)
  
  ## Use prim's algorithm to get edges (i.e. the cells connected by each edge)
  # Input is the adjacency matrix of the graph (i.e. output from -1 * apcluster::negDistMat(df of coords))
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
      # Set initial temp values for weight and vertex
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
      
      # Current edge connects the min_vertex and curr_vertex
      tree_edges[num_edges + 1, ] <- c(min_vertex, curr_vertex)
      selected[min_vertex] <- TRUE
      num_edges <- num_edges + 1
    }
    return(tree_edges)
  }
  
  tree_edges <- prims_algorithm(adj_mat)
  
  ### Determine width of cylinders so that cylinders further away are thinner
  get_tree_depth <- function(tree_edges) {
    
    tree_edges <- data.frame(tree_edges)
    colnames(tree_edges) <- c("vertex1", "vertex2")
    
    # Set the initial depth of each tree_edge to be NA.
    tree_edges$depth <- NA
    
    # Get vertices on the 'outskirts' of MST (leaf_vertices which have a depth of 1)
    tree_vertices <- c(tree_edges[ , 1], tree_edges[ , 2])
    leaf_vertices <- as.numeric(names(table(tree_vertices))[table(tree_vertices) == 1])
    
    # Start with leaf_vertices
    curr_vertices <- leaf_vertices
    curr_depth <- 1
    
    while (NA %in% tree_edges$depth) {
      
      # New vertices will be those adjacent to the current vertices
      new_vertices <- c()
      
      # Check each current vertex
      for (vertex in curr_vertices) {
        # Start with vertex1
        curr_edges <- which(tree_edges$vertex1 == vertex)
        tree_edges[curr_edges, "depth"][is.na(tree_edges[curr_edges, "depth"])] <- curr_depth
        new_vertices <- c(new_vertices, tree_edges[curr_edges, "vertex2"])
        
        # Then vertex2
        curr_edges <- which(tree_edges$vertex2 == vertex)
        tree_edges[curr_edges, "depth"][is.na(tree_edges[curr_edges, "depth"])] <- curr_depth
        new_vertices <- c(new_vertices, tree_edges[curr_edges, "vertex1"])
        
        # Only keep unique vertices
        new_vertices <- unique(new_vertices)
      }
      
      curr_depth <- curr_depth + 1
      curr_vertices <- new_vertices
    }
    
    return(tree_edges)
  }
  
  tree_edges <- get_tree_depth(tree_edges)
  
  ## Get cluster properties using edge data
  network_ring_properties <- list()
  max_depth <- max(tree_edges[["depth"]])
  
  for (i in seq(n_edges)) {
    start_loc <- as.numeric(random_coords[tree_edges[i, "vertex1"], ])
    end_loc <- as.numeric(random_coords[tree_edges[i, "vertex2"], ])
    curr_width <- (1 - 0.10 * (max_depth - tree_edges[i, "depth"])) * width # 10% decrease with each depth
    
    # Very unlikely case when width is negative, just ignore these cylinders
    if (width < 0) {
      width <- 0
    }
    
    network_ring_properties[[i]] <- list(shape = "cylinder",
                                         cluster_cell_types = cluster_cell_types,
                                         cluster_cell_proportions = cluster_cell_proportions,
                                         radius = curr_width,
                                         start_loc = start_loc,
                                         end_loc = end_loc,
                                         ring_cell_types = ring_cell_types,
                                         ring_cell_proportions = ring_cell_proportions,
                                         ring_width = ring_width)
  }
  
  network_spe <- simulate_rings3D(spe,
                                  ring_properties = network_ring_properties,
                                  plot_image = F)
  
  # Update current meta data
  metadata <- spe@metadata
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  ring_properties[["cylinders"]] <- network_ring_properties # Include metadata of cylinders used to make up network
  metadata[["simulation"]][[paste("cluster", length(metadata[["simulation"]]), sep = "_")]] <- ring_properties
  
  network_spe@metadata <- metadata
  
  return(network_spe)
}
