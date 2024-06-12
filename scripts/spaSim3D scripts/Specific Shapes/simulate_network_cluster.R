simulate_network_cluster <- function(bg_sample, cluster_properties) {  
  
  # Get network properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  n_edges <- cluster_properties$n_edges
  width <- cluster_properties$width
  centre_loc <- cluster_properties$centre_loc
  radius <- cluster_properties$radius
  
  # Number of vertices is always one more than the number of edges for the MST will we make
  n_vertices <- n_edges + 1 

  ## Choose cells within the radius of the centre_loc
  R <- radius^2
  
  D <- (bg_sample$Cell.X.Position - centre_loc[1])^2 +
       (bg_sample$Cell.Y.Position - centre_loc[2])^2 +
       (bg_sample$Cell.Z.Position - centre_loc[3])^2
  
  cells_chosen <- bg_sample[D <= R, ]
  
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
  cluster_properties <- list()
  max_depth <- max(tree_edges[["Depth"]])
  
  for (i in seq(n_vertices - 1)) {
    start_loc <- as.numeric(cells_chosen[tree_edges[i, "Cell1"], ])
    end_loc <- as.numeric(cells_chosen[tree_edges[i, "Cell2"], ])
    curr_width <- (1 - 0.10 * (max_depth - tree_edges[i, "Depth"])) * width # 10% decrease with each depth
    
    # Very unlikely case when width is negative, just ignore these cylinders
    if (width < 0) {
      width <- 0
    }
    
    cluster_properties[[i]] <- list(shape = "Cylinder",
                                    cluster_cell_types = cluster_cell_types,
                                    cluster_cell_proportions = cluster_cell_proportions,
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
