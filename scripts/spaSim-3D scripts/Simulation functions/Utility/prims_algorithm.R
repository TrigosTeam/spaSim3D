## Prim's algorithm function
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
