simulate_ellipsoid_cluster <- function(bg_spe, cluster_properties) {
  
  # Get ellipsoid properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  x_radius <- cluster_properties$x_radius
  y_radius <- cluster_properties$y_radius
  z_radius <- cluster_properties$z_radius
  centre_loc <- cluster_properties$centre_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
  
  
  # Rotation angles
  theta <- cluster_properties$y_z_rotation * (pi/180) # rotation in x-axis
  alpha <- cluster_properties$x_z_rotation * (pi/180) # rotation in y-axis
  beta  <- cluster_properties$x_y_rotation * (pi/180) # rotation in z-axis
  
  # 3x3 Transformation matrix (T_M) using rotation angles
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
  
  ## Change cell types in the ellipsoid cluster
  spe_coords <- data.frame(spatialCoords(bg_spe))
  
  # Adjust x, y and z coordinates relative to the ellipsoid centre
  x <- spe_coords$Cell.X.Position - centre_loc[1]
  y <- spe_coords$Cell.Y.Position - centre_loc[2]
  z <- spe_coords$Cell.Z.Position - centre_loc[3]
  
  # Transform  x, y and z coordinates using rotation transformation matrix
  x <- T_M[1, 1] * x + T_M[1, 2] * y + T_M[1, 3] * z
  y <- T_M[2, 1] * x + T_M[2, 2] * y + T_M[2, 3] * z
  z <- T_M[3, 1] * x + T_M[3, 2] * y + T_M[3, 3] * z
  
  bg_spe[["Cell.Type"]] <- ifelse((x / x_radius)^2 +
                                    (y / y_radius)^2 +
                                    (z / z_radius)^2 <= 1,
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  
  # Update current meta data
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- cluster_properties
  
  return(bg_spe)
}