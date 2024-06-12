simulate_sphere_ring <- function(bg_sample, ring_properties) {
  
  # Get sphere ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  radius <- ring_properties$radius
  centre_loc <- ring_properties$centre_loc
  
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  # Get number of unique cluster cell types
  n_cluster_cell_types <- length(cluster_cell_types)
  
  # Get number of unique ring cell types
  n_ring_cell_types <- length(ring_cell_types)
  
  for (i in seq_len(n_cells)) {
    # Get x, y, z coordinate of current cell
    x <- bg_sample[i, "Cell.X.Position"]
    y <- bg_sample[i, "Cell.Y.Position"]
    z <- bg_sample[i, "Cell.Z.Position"]
    
    # Using radius of sphere
    R1 <- radius^2
    
    # Using radius of sphere with ring
    R2 <- (radius + ring_width)^2
    
    # Calculate distance of current cell from sphere centre
    D <- (x - centre_loc[1])^2 + (y - centre_loc[2])^2 + (z - centre_loc[3])^2
    
    if (D < R1) { 
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
    else if (D < R2) {
      # Random number will determine the ring_cell_type of the cell
      random <- stats::runif(1)
      
      # Start with the first cell
      n <- 1
      current_proportion <- 0
      
      while (n <= n_ring_cell_types){
        current_proportion <- current_proportion + ring_cell_proportions[n]
        if (random <= current_proportion) {
          bg_sample[i, "Cell.Type"] <- ring_cell_types[n]
          break
        }
        n <- n + 1
      }
    }
  }
  
  return(bg_sample)
}
