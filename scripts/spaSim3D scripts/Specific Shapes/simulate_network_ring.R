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
