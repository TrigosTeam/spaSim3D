simulate_ellipsoid_dr <- function(bg_sample, dr_properties) {
  
  # Get ellipsoid double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  x_radius <- dr_properties$x_radius
  y_radius <- dr_properties$y_radius
  z_radius <- dr_properties$z_radius
  centre_loc <- dr_properties$centre_loc
  
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  # Rotation angles
  theta <- dr_properties$y_z_rotation * (pi/180) # rotation in x-axis
  alpha <- dr_properties$x_z_rotation * (pi/180) # rotation in y-axis
  beta  <- dr_properties$x_y_rotation * (pi/180) # rotation in z-axis
  
  # 3x3 Transformation matrix using rotation angles
  T_M <- matrix(data = c(cos(alpha) * cos(beta), 
                         cos(alpha) * sin(beta), 
                         sin(alpha),
                         -sin(theta) * sin(alpha) * cos(beta) - cos(theta) * sin(beta),
                         -sin(theta) * sin(alpha) * sin(beta) + cos(theta) * cos(beta),
                         sin(theta) * cos(alpha),
                         -cos(theta) * sin(alpha) * cos(beta) + sin(theta) * sin(beta),
                         -cos(theta) * sin(alpha) * sin(beta) - sin(theta) * cos(beta),
                         cos(theta) * cos(alpha)), 
                nrow = 3, 
                ncol = 3, 
                byrow = TRUE)
  
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
    x <- bg_sample[i, "Cell.X.Position"] - centre_loc[1]
    y <- bg_sample[i, "Cell.Y.Position"] - centre_loc[2]
    z <- bg_sample[i, "Cell.Z.Position"] - centre_loc[3]
    
    x_new <- T_M[1, 1] * x + T_M[1, 2] * y + T_M[1, 3] * z
    y_new <- T_M[2, 1] * x + T_M[2, 2] * y + T_M[2, 3] * z
    z_new <- T_M[3, 1] * x + T_M[3, 2] * y + T_M[3, 3] * z
    
    # Using radius of ellipsoid
    D1 <- (x_new/x_radius)^2 + 
          (y_new/y_radius)^2 + 
          (z_new/z_radius)^2
    
    # Using radius of ellipsoid with inner ring
    D2 <- (x_new/(x_radius + inner_ring_width))^2 + 
          (y_new/(y_radius + inner_ring_width))^2 + 
          (z_new/(z_radius + inner_ring_width))^2
    
    # Using radius of ellipsoid with inner and outer ring
    D3 <- (x_new/(x_radius + inner_ring_width + outer_ring_width))^2 + 
          (y_new/(y_radius + inner_ring_width + outer_ring_width))^2 + 
          (z_new/(z_radius + inner_ring_width + outer_ring_width))^2
    
    if (D1 <= 1) { 
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
    else if (D2 <= 1) {
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
    else if (D3 <= 1) {
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
