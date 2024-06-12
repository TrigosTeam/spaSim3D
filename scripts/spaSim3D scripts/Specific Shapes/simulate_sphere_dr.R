simulate_sphere_dr <- function(bg_sample, dr_properties) {
  
  # Get sphere double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  radius <- dr_properties$radius
  centre_loc <- dr_properties$centre_loc
  
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  # Get number of unique cluster cell types
  n_cluster_cell_types <- length(cluster_cell_types)
  
  # Get number of unique inner ring cell types
  n_inner_ring_cell_types <- length(inner_ring_cell_types)
  
  # Get number of unique outer ring cell types
  n_outer_ring_cell_types <- length(outer_ring_cell_types)
  
  for (i in seq_len(n_cells)) {
    # Get x, y, z coordinate of current cell
    x <- bg_sample[i, "Cell.X.Position"]
    y <- bg_sample[i, "Cell.Y.Position"]
    z <- bg_sample[i, "Cell.Z.Position"]
    
    # Using radius of sphere
    R1 <- radius^2
    
    # Using radius of sphere with inner ring
    R2 <- (radius + inner_ring_width)^2
    
    # Using radius of sphere with inner and outer ring
    R3 <- (radius + inner_ring_width + outer_ring_width)^2
    
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
      # Random number will determine the inner_ring_cell_type of the cell
      random <- stats::runif(1)
      
      # Start with the first cell
      n <- 1
      current_proportion <- 0
      
      while (n <= n_inner_ring_cell_types){
        current_proportion <- current_proportion + inner_ring_cell_proportions[n]
        if (random <= current_proportion) {
          bg_sample[i, "Cell.Type"] <- inner_ring_cell_types[n]
          break
        }
        n <- n + 1
      }
    }
    else if (D < R3) {
      # Random number will determine the outer_ring_cell_type of the cell
      random <- stats::runif(1)
      
      # Start with the first cell
      n <- 1
      current_proportion <- 0

      while (n <= n_outer_ring_cell_types){
        current_proportion <- current_proportion + outer_ring_cell_proportions[n]
        if (random <= current_proportion) {
          bg_sample[i, "Cell.Type"] <- outer_ring_cell_types[n]
          break
        }
        n <- n + 1
      }
    }
  }
  
  return (bg_sample)
}
