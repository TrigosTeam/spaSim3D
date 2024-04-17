simulate_ellipsoid_ring <- function(bg_sample, ring_properties) {
  
  # Get ellipsoid properties
  cell_type <- ring_properties$name_of_cluster_cell
  infiltration_types <- ring_properties$infiltration_types
  infiltration_proportions <- ring_properties$infiltration_proportions
  x_radius <- ring_properties$x_radius
  y_radius <- ring_properties$y_radius
  z_radius <- ring_properties$z_radius
  centre_loc <- ring_properties$centre_loc
  
  ring_cell_type <- ring_properties$name_of_ring_cell
  ring_width <- ring_properties$ring_width
  ring_infiltration_types <- ring_properties$ring_infiltration_types
  ring_infiltration_proportions <- ring_properties$ring_infiltration_proportions
  
  theta <- ring_properties$y_z_rotation # in x-axis
  alpha <- ring_properties$x_z_rotation # in y-axis
  beta  <- ring_properties$x_y_rotation # in z-axis
  
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
    
    D2 <- (x_new/(x_radius + ring_width))^2 + 
          (y_new/(y_radius + ring_width))^2 + 
          (z_new/(z_radius + ring_width))^2
    
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
      
      n_ring_infiltration_types <- length(ring_infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this ring
      pheno <- ring_cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_ring_infiltration_types) {
        current_p <- current_p + ring_infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- ring_infiltration_types[n]
          break
        }
        n <- n+1
      }
    }
    bg_sample[index, "Cell.Type"] <- pheno
    
  }
  return(bg_sample)
}
