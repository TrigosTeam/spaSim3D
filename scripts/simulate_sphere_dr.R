simulate_sphere_dr <- function(bg_sample, dr_properties) {
  
  # Get sphere double ring properties
  cell_type <- dr_properties$name_of_cluster_cell
  infiltration_types <- dr_properties$infiltration_types
  infiltration_proportions <- dr_properties$infiltration_proportions
  radius <- dr_properties$radius
  centre_loc <- dr_properties$centre_loc
  
  inner_ring_cell_type <- dr_properties$name_of_inner_ring_cell
  inner_ring_width <- dr_properties$inner_ring_width
  inner_ring_infiltration_types <- dr_properties$inner_ring_infiltration_types
  inner_ring_infiltration_proportions <- dr_properties$inner_ring_infiltration_proportions
  
  outer_ring_cell_type <- dr_properties$name_of_outer_ring_cell
  outer_ring_width <- dr_properties$outer_ring_width
  outer_ring_infiltration_types <- dr_properties$outer_ring_infiltration_types
  outer_ring_infiltration_proportions <- dr_properties$outer_ring_infiltration_proportions
  
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  for (i in seq_len(n_cells)) {
    # Get x, y, z and phenotype of ith cell
    x <- bg_sample[i, "Cell.X.Position"]
    y <- bg_sample[i, "Cell.Y.Position"]
    z <- bg_sample[i, "Cell.Z.Position"]
    pheno <- bg_sample[i, "Cell.Type"]
    
    R1 <- radius^2
    R2 <- (radius + inner_ring_width)^2
    R3 <- (radius + inner_ring_width + outer_ring_width)^2
    
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
      # in the region of inner ring, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_inner_ring_infiltration_types <- length(inner_ring_infiltration_types)
      
      # default `Cell.Type` is cell type of interest of inner ring
      pheno <- inner_ring_cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_inner_ring_infiltration_types){
        current_p <- current_p + inner_ring_infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- inner_ring_infiltration_types[n]
          break
        }
        n <- n + 1
      }
    }
    else if (D < R3) {
      # in the region of outer ring, generate random number to decide the `Cell.Type`
      random <- stats::runif(1)
      
      n_outer_ring_infiltration_types <- length(outer_ring_infiltration_types)
      
      # default `Cell.Type` is cell type of interest of outer ring
      pheno <- outer_ring_cell_type
      # if the random number falls in the range of an infiltration proportion,
      # pheno will be the corresponding infiltraiton type
      n <- 1 # start from the first proportion
      current_p <- 0
      while (n <= n_outer_ring_infiltration_types){
        current_p <- current_p + outer_ring_infiltration_proportions[n]
        if (random <= current_p) {
          pheno <- outer_ring_infiltration_types[n]
          break
        }
        n <- n + 1
      }
    }
    bg_sample[i, "Cell.Type"] <- pheno
    
  }
  return(bg_sample)
}
