simulate_ellipsoid_dr <- function(bg_spe, dr_properties) {
  
  # Get ellipsoid double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  x_radius <- dr_properties$x_radius
  y_radius <- dr_properties$y_radius
  z_radius <- dr_properties$z_radius
  centre_loc <- dr_properties$centre_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
  
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  
  ## Check number of inner ring cell types matches the number of inner ring cell proportions
  if (length(inner_ring_cell_types) != length(inner_ring_cell_proportions)) stop("Number of inner ring cell types doesn't match number of inner ring cell proportion.")
  
  ## Check inner ring cell proportions are not negative or greater than 1
  if (sum(inner_ring_cell_proportions < 0 | inner_ring_cell_proportions > 1) != 0) stop("Inner ring cell proportions cannot be negative or greater than 1")
  
  ## Check inner ring cell proportions add up to 1
  if (sum(inner_ring_cell_proportions) != 1) stop("Sum of inner ring cell proportions is NOT 1")
  
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  ## Check number of outer ring cell types matches the number of outer ring cell proportions
  if (length(outer_ring_cell_types) != length(outer_ring_cell_proportions)) stop("Number of outer ring cell types doesn't match number of outer ring cell proportion.")
  
  ## Check outer ring cell proportions are not negative or greater than 1
  if (sum(outer_ring_cell_proportions < 0 | outer_ring_cell_proportions > 1) != 0) stop("Outer ring cell proportions cannot be negative or greater than 1")
  
  ## Check outer ring cell proportions add up to 1
  if (sum(outer_ring_cell_proportions) != 1) stop("Sum of outer ring cell proportions is NOT 1")
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(bg_spe), 
                   "Cell.Type" = bg_spe[["Cell.Type"]],
                   "Cell.ID" = bg_spe[["Cell.ID"]])
  
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
  n_cells <- nrow(df)
  
  # Get number of unique cluster cell types
  n_cluster_cell_types <- length(cluster_cell_types)
  
  # Get number of unique inner ring cell types
  n_inner_ring_cell_types <- length(inner_ring_cell_types)
  
  # Get number of unique outer ring cell types
  n_outer_ring_cell_types <- length(outer_ring_cell_types)
  
  for (i in seq_len(n_cells)) {
    # Get x, y, z coordinate of current cell
    x <- df[i, "Cell.X.Position"] - centre_loc[1]
    y <- df[i, "Cell.Y.Position"] - centre_loc[2]
    z <- df[i, "Cell.Z.Position"] - centre_loc[3]
    
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
          df[i, "Cell.Type"] <- cluster_cell_types[n]
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
          df[i, "Cell.Type"] <- inner_ring_cell_types[n]
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
          df[i, "Cell.Type"] <- outer_ring_cell_types[n]
          break
        }
        n <- n + 1
      }
    }
  }
  
  # Update current meta data
  metadata <- bg_spe@metadata
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  metadata[[paste("cluster", length(metadata), sep="_")]] <- dr_properties
  
  # Convert data frame to spe object
  cluster_spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = metadata)
  
  return(cluster_spe)
}
