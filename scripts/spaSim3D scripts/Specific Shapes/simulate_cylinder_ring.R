simulate_cylinder_ring <- function(bg_sample, ring_properties) {
  
  # Get cylinder ring properties
  cell_type <- ring_properties$name_of_cluster_cell
  infiltration_types <- ring_properties$infiltration_types
  infiltration_proportions <- ring_properties$infiltration_proportions
  radius <- ring_properties$radius
  start_loc <- ring_properties$start_loc
  end_loc <- ring_properties$end_loc
  
  ring_cell_type <- ring_properties$name_of_ring_cell
  ring_width <- ring_properties$ring_width
  ring_infiltration_types <- ring_properties$ring_infiltration_types
  ring_infiltration_proportions <- ring_properties$ring_infiltration_proportions
  
  
  # Get number of cells
  n_cells <- nrow(bg_sample)

  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  i <- 1
  
  while (i <= n_cells) {
    # Get x, y, z and phenotype of ith cell
    x <- bg_sample[i, "Cell.X.Position"]
    y <- bg_sample[i, "Cell.Y.Position"]
    z <- bg_sample[i, "Cell.Z.Position"]
    pheno <- bg_sample[i, "Cell.Type"]
    
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
            (v1[1]*v2[2] - v1[2]*v2[1])^2)/
      (v1[1]^2 + v1[2]^2 + v1[3]^2)
    
    # Get maximum distance without and with ring squared
    R1 <- radius^2
    R2 <- (radius + ring_width)^2
    
    if (D < R1) { 
      # in the region of cluster, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_infiltration_types <- length(infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this cluster
      pheno <- cell_type 
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_infiltration_types){
        current_p <- current_p + infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- infiltration_types[n]
          break
        }
        n <- n + 1
      }
    }
    else if (D < R2) {
      # in the region of ring, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_ring_infiltration_types <- length(ring_infiltration_types)
      
      # default `Cell.Type` is cell type of interest of this ring
      pheno <- ring_cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_ring_infiltration_types){
        current_p <- current_p + ring_infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- ring_infiltration_types[n]
          break
        }
        n <- n + 1
      }
    }

    if (pheno == "Void") { 
      bg_sample <- bg_sample[-c(i), ]
      n_cells <- n_cells - 1
        
    } else {
      bg_sample[i, "Cell.Type"] <- pheno  
      i <- i + 1
    }
  
  }
  return(bg_sample)
}
