
simulate_cylinder_cluster <- function(bg_sample, cluster_properties) {
  
  # Get cylinder properties
  cell_type <- cluster_properties$name_of_cluster_cell
  infiltration_types <- cluster_properties$infiltration_types
  infiltration_proportions <- cluster_properties$infiltration_proportions
  radius <- cluster_properties$radius
  start_loc <- cluster_properties$start_loc
  end_loc <- cluster_properties$end_loc
  
  # Get number of cells
  n_cells <- nrow(bg_sample)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  for (i in seq_len(n_cells)) {
    # Get x, y, z and phenotype of ith cell
    x <- bg_sample[i, "Cell.X.Position"]
    y <- bg_sample[i, "Cell.Y.Position"]
    z <- bg_sample[i, "Cell.Z.Position"]
    pheno <- bg_sample[i, "Cell.Type"]
    
    # Get vector between from point to start_loc
    v2 <- c(x, y, z) - start_loc
    
    # Get perpendicular distance squared between point and line
    D <- ((v1[2]*v2[3] - v1[3]*v2[2])^2 + 
            (v1[1]*v2[3] - v1[3]*v2[1])^2 + 
            (v1[1]*v2[2] - v1[2]*v2[1])^2)/
      (v1[1]^2 + v1[2]^2 + v1[3]^2)
    
    # Get maximum distance squared
    R <- radius^2
    
    if (D < R){ 
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
        n <- n+1
      }
    }
    bg_sample[i, "Cell.Type"] <- pheno
  }
  return(bg_sample)
}
