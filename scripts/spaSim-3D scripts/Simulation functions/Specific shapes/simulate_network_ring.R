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
  tree_edges <- prims_algorithm(adj_mat)
  
  ### Determine width of cylinders so that cylinders further away are thinner
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
  
  network_spe <- simulate_rings3D(bg_spe,
                                  ring_properties = network_ring_properties,
                                  plot_image = F)
  
  # Update current meta data
  metadata <- bg_spe@metadata
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  ring_properties[["cylinders"]] <- network_ring_properties # Include metadata of cylinders used to make up network
  metadata[["simulation"]][[paste("cluster", length(metadata[["simulation"]]), sep = "_")]] <- ring_properties
  
  network_spe@metadata <- metadata
  
  return(network_spe)
  
}
