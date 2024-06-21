simulate_ellipsoid_ring <- function(bg_spe, ring_properties) {
  
  # Get ellipsoid ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  x_radius <- ring_properties$x_radius
  y_radius <- ring_properties$y_radius
  z_radius <- ring_properties$z_radius
  centre_loc <- ring_properties$centre_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
  
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  ## Check number of ring cell types matches the number of cell proportions
  if (length(ring_cell_types) != length(ring_cell_proportions)) stop("Number of ring cell types doesn't match number of ring cell proportion.")
  
  ## Check ring cell proportions are not negative or greater than 1
  if (sum(ring_cell_proportions < 0 | ring_cell_proportions > 1) != 0) stop("Ring cell proportions cannot be negative or greater than 1")
  
  ## Check ring cell proportions add up to 1
  if (sum(ring_cell_proportions) != 1) stop("Sum of ring cell proportions is NOT 1")
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(bg_spe), 
                   "Cell.Type" = bg_spe[["Cell.Type"]],
                   "Cell.ID" = bg_spe[["Cell.ID"]])
  
  # Rotation angles
  theta <- ring_properties$y_z_rotation * (pi/180) # rotation in x-axis
  alpha <- ring_properties$x_z_rotation * (pi/180) # rotation in y-axis
  beta  <- ring_properties$x_y_rotation * (pi/180) # rotation in z-axis
  
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
  n_cells <- nrow(df)
  
  # Get number of unique cluster cell types
  n_cluster_cell_types <- length(cluster_cell_types)
  
  # Get number of unique ring cell types
  n_ring_cell_types <- length(ring_cell_types)
  
  for (i in seq_len(n_cells)) {
    # Get x, y, z coordinate of current cell
    x <- df[i, "Cell.X.Position"] - centre_loc[1]
    y <- df[i, "Cell.Y.Position"] - centre_loc[2]
    z <- df[i, "Cell.Z.Position"] - centre_loc[3]
    
    x_new <- T_M[1, 1] * x + T_M[1, 2] * y + T_M[1, 3] * z
    y_new <- T_M[2, 1] * x + T_M[2, 2] * y + T_M[2, 3] * z
    z_new <- T_M[3, 1] * x + T_M[3, 2] * y + T_M[3, 3] * z
    
    # Using radii of ellipsoid
    D1 <- (x_new / x_radius)^2 + 
          (y_new / y_radius)^2 + 
          (z_new / z_radius)^2
    
    # Using radii of ellipsoid with ring
    D2 <- (x_new / (x_radius + ring_width))^2 + 
          (y_new / (y_radius + ring_width))^2 + 
          (z_new / (z_radius + ring_width))^2
    
    if (D1 <= 1) { 
      # Random number will determine the cluster_cell_type of the cell
      random <- stats::runif(1)
      
      # Start with the first cell
      n <- 1
      current_proportion <- 0
      
      while (n <= n_cluster_cell_types){
        current_proportion <- current_proportion + cluster_cell_proportions[n]
        if (random <= current_proportion) {
          df[i, "Cell.Type"] <- cluster_cell_types[n]
          break
        }
        n <- n + 1
      }
    }
    else if (D2 <= 1) {
      # Random number will determine the ring_cell_type of the cell
      random <- stats::runif(1)
      
      # Start with the first cell
      n <- 1
      current_proportion <- 0
      
      while (n <= n_ring_cell_types){
        current_proportion <- current_proportion + ring_cell_proportions[n]
        if (random <= current_proportion) {
          df[i, "Cell.Type"] <- ring_cell_types[n]
          break
        }
        n <- n + 1
      }
    }
  }
  
  # Update current meta data
  metadata <- bg_spe@metadata
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  metadata[[paste("cluster", length(metadata), sep="_")]] <- ring_properties
  
  # Convert data frame to spe object
  cluster_spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = metadata)
  
  return(cluster_spe)
}
