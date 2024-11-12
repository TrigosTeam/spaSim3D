simulate_network_dr <- function(spe, dr_properties) {  
  
  # Check input parameters
  input_parameters <- dr_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get network double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  n_edges <- dr_properties$n_edges
  width <- dr_properties$width
  centre_loc <- dr_properties$centre_loc
  radius <- dr_properties$radius
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  # Number of vertices is always one more than the number of edges for the MST will we make
  n_vertices <- n_edges + 1 
  
  ## Generate n_vertices random points with coords inside a sphere with given radius and centre loc. 
  # Starting with 1000 points inside a cube should be a good enough buffer, unless the user wants more than 1000 edges...
  # Let's stop them from inputting more than 99
  max_edges <- 99
  if (n_edges > max_edges) stop("Only networks with less than 100 edges can be simulated.")
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
  tree_edges <- prims_algorithm(adj_mat)
  
  ### Determine width of cylinders so that cylinders further away are thinner
  tree_edges <- get_tree_depth(tree_edges)
  
  ## Get cluster properties using edge data
  network_dr_properties <- list()
  max_depth <- max(tree_edges[["depth"]])
  
  for (i in seq(n_edges)) {
    start_loc <- as.numeric(random_coords[tree_edges[i, "vertex1"], ])
    end_loc <- as.numeric(random_coords[tree_edges[i, "vertex2"], ])
    curr_width <- (1 - 0.10 * (max_depth - tree_edges[i, "depth"])) * width # 10% decrease with each depth
    
    # Very unlikely case when width is negative, just ignore these cylinders
    if (width < 0) {
      width <- 0
    }
    
    network_dr_properties[[i]] <- list(shape = "cylinder",
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
  
  network_spe <- simulate_double_rings3D(spe,
                                         dr_properties = network_dr_properties,
                                         plot_image = F)
  
  # Update current meta data
  metadata <- spe@metadata
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  dr_properties[["cylinders"]] <- network_dr_properties # Include metadata of cylinders used to make up network
  metadata[["simulation"]][[paste("cluster", length(metadata[["simulation"]]), sep = "_")]] <- dr_properties
  
  network_spe@metadata <- metadata
  
  return(network_spe)
}
