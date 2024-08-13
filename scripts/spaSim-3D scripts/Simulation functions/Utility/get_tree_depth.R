get_tree_depth <- function(tree_edges) {
  
  tree_edges <- data.frame(tree_edges)
  colnames(tree_edges) <- c("vertex1", "vertex2")
  tree_edges$depth <- NA # If tree_edge is not NA, we have already accounted for it
  
  # Get cells on the 'outskirts' of MST (i.e. leaf_vertices)
  tree_vertices <- c(tree_edges[ , 1], tree_edges[ , 2])
  
  leaf_vertices <- names(table(tree_vertices))[table(tree_vertices) == 1]
  leaf_vertices <- as.numeric(leaf_vertices)
  
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