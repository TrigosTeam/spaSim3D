simulate_ellipsoid_dr <- function(bg_sample, dr_properties) {
  
  # Get ellipsoid properties
  cell_type <- dr_properties$name_of_cluster_cell
  infiltration_types <- dr_properties$infiltration_types
  infiltration_proportions <- dr_properties$infiltration_proportions
  x_radius <- dr_properties$x_radius
  y_radius <- dr_properties$y_radius
  z_radius <- dr_properties$z_radius
  centre_loc <- dr_properties$centre_loc
  
  inner_ring_cell_type <- dr_properties$name_of_inner_ring_cell
  inner_ring_width <- dr_properties$inner_ring_width
  inner_ring_infiltration_types <- dr_properties$inner_ring_infiltration_types
  inner_ring_infiltration_proportions <- dr_properties$inner_ring_infiltration_proportions
  
  outer_ring_cell_type <- dr_properties$name_of_outer_ring_cell
  outer_ring_width <- dr_properties$outer_ring_width
  outer_ring_infiltration_types <- dr_properties$outer_ring_infiltration_types
  outer_ring_infiltration_proportions <- dr_properties$outer_ring_infiltration_proportions
  
  theta <- dr_properties$y_z_rotation # in x-axis
  alpha <- dr_properties$x_z_rotation # in y-axis
  beta  <- dr_properties$x_y_rotation # in z-axis
  
  a <- cos(alpha) * cos(beta)
  b <- cos(alpha) * sin(beta)
  c <- sin(alpha)
  d <- -sin(theta) * sin(alpha) * cos(beta) - cos(theta) * sin(beta)
  e <- -sin(theta) * sin(alpha) * sin(beta) + cos(theta) * cos(beta)
  f <- sin(theta) * cos(alpha)
  g <- -cos(theta) * sin(alpha) * cos(beta) + sin(theta) * sin(beta)
  h <- -cos(theta) * sin(alpha) * sin(beta) - sin(theta) * cos(beta)
  i <- cos(theta) * cos(alpha)
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  for (index in seq_len(n_cells)) {
    # Get x, y, z and phenotype of ith cell
    x <- bg_sample[index, "Cell.X.Position"] - centre_loc[1]
    y <- bg_sample[index, "Cell.Y.Position"] - centre_loc[2]
    z <- bg_sample[index, "Cell.Z.Position"] - centre_loc[3]
    pheno <- bg_sample[index, "Cell.Type"]
    
    x_new <- a * x + b * y + c * z
    y_new <- d * x + e * y + f * z
    z_new <- g * x + h * y + i * z
    
    D1 <- (x_new/x_radius)^2 + 
      (y_new/y_radius)^2 + 
      (z_new/z_radius)^2
    
    D2 <- (x_new/(x_radius + inner_ring_width))^2 + 
      (y_new/(y_radius + inner_ring_width))^2 + 
      (z_new/(z_radius + inner_ring_width))^2
    
    D3 <- (x_new/(x_radius + inner_ring_width + outer_ring_width))^2 + 
      (y_new/(y_radius + inner_ring_width + outer_ring_width))^2 + 
      (z_new/(z_radius + inner_ring_width + outer_ring_width))^2
    
    if (D1 <= 1) { 
      # in the region of cluster, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_infiltration_types <- length(infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this cluster
      pheno <- cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_infiltration_types) {
        current_p <- current_p + infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- infiltration_types[n]
          break
        }
        n <- n+1
      }
    }
    else if (D2 <= 1) {
      # in the region of ring, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_inner_ring_infiltration_types <- length(inner_ring_infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this ring
      pheno <- inner_ring_cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_inner_ring_infiltration_types) {
        current_p <- current_p + inner_ring_infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- inner_ring_infiltration_types[n]
          break
        }
        n <- n+1
      }
    }
    else if (D3 <= 1) {
      # in the region of ring, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_outer_ring_infiltration_types <- length(outer_ring_infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this ring
      pheno <- outer_ring_cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_outer_ring_infiltration_types) {
        current_p <- current_p + outer_ring_infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- outer_ring_infiltration_types[n]
          break
        }
        n <- n+1
      }
    }
    bg_sample[index, "Cell.Type"] <- pheno
    
  }
  return(bg_sample)
}
