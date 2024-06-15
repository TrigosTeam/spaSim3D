simulate_cylinder_cluster <- function(bg_sample, cluster_properties) {
  
  # Get cylinder properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  radius <- cluster_properties$radius
  start_loc <- cluster_properties$start_loc
  end_loc <- cluster_properties$end_loc
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  # Get number of unique cell types
  n_cluster_cell_types <- length(cluster_cell_types)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  i <- 1
  
  while (i <= n_cells) {
    # Get x, y, z coordinate of current cell
    x <- bg_sample[i, "Cell.X.Position"]
    y <- bg_sample[i, "Cell.Y.Position"]
    z <- bg_sample[i, "Cell.Z.Position"]
    
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
          (v1[1]*v2[2] - v1[2]*v2[1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2)
    
    # Dumb case where the start and end loc is the same
    if (is.nan(D)) D <- Inf
    
    # Get maximum distance squared
    R <- radius^2
    
    if (D < R) { 
      # Random number will determine the cluster_cell_type of the cell
      random <- stats::runif(1)
      
      # Start with the first cell
      n <- 1
      current_proportion <- 0
      
      while (n <= n_cluster_cell_types){
        current_proportion <- current_proportion + cluster_cell_proportions[n]
        if (random <= current_proportion) {
          bg_sample[i, "Cell.Type"] <- cluster_cell_types[n]
          break
        }
        n <- n + 1
      }
    }
    
    if (bg_sample[i, "Cell.Type"] == "Void") { 
      bg_sample <- bg_sample[-c(i), ]
      n_cells <- n_cells - 1
      
    } else {
      i <- i + 1
    }
  }
  
  return (bg_sample)
}
