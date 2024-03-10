simulate_sphere_ring <- function(bg_sample, ring_properties) {
  
  # Get sphere immune ring properties
  cell_type <- ring_properties$name_of_cluster_cell
  infiltration_types <- ring_properties$infiltration_types
  infiltration_proportions <- ring_properties$infiltration_proportions
  radius <- ring_properties$radius
  centre_loc <- ring_properties$centre_loc
  
  ring_cell_type <- ring_properties$name_of_ring_cell
  ring_width <- ring_properties$ring_width
  ring_infiltration_types <- ring_properties$ring_infiltration_types
  ring_infiltration_proportions <- ring_properties$ring_infiltration_proportions
  
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  for (i in seq_len(n_cells)) {
    # Get x, y, z and phenotype of ith cell
    x <- bg_sample[i, "Cell.X.Position"]
    y <- bg_sample[i, "Cell.Y.Position"]
    z <- bg_sample[i, "Cell.Z.Position"]
    pheno <- bg_sample[i, "Cell.Type"]
    
    R1 <- radius^2
    R2 <- (radius + ring_width)^2
    
    D <- (x - centre_loc[1])^2 + (y - centre_loc[2])^2 + (z - centre_loc[3])^2
    
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
    bg_sample[i, "Cell.Type"] <- pheno
    
  }
  return(bg_sample)
}
