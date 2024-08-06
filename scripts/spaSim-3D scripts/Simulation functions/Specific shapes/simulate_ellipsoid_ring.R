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
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  ## Check number of ring cell types matches the number of cell proportions
  if (length(ring_cell_types) != length(ring_cell_proportions)) stop("Number of ring cell types doesn't match number of ring cell proportion.")
  
  ## Check ring cell proportions are not negative or greater than 1
  if (sum(ring_cell_proportions < 0 | ring_cell_proportions > 1) != 0) stop("Ring cell proportions cannot be negative or greater than 1")
  
  ## Check ring cell proportions add up to 1
  if (!all.equal(sum(ring_cell_proportions), 1)) stop("Sum of ring cell proportions is NOT 1")
  
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
  
  
  
  ## Change cell types in the ellipsoid cluster
  spe_coords <- data.frame(spatialCoords(bg_spe))
  
  # Adjust x, y and z coordinates relative to the ellipsoid centre
  x <- spe_coords$Cell.X.Position - centre_loc[1]
  y <- spe_coords$Cell.Y.Position - centre_loc[2]
  z <- spe_coords$Cell.Z.Position - centre_loc[3]
  
  # Transform  x, y and z coordinates using rotation transformation matrix
  x_new <- T_M[1, 1] * x + T_M[1, 2] * y + T_M[1, 3] * z
  y_new <- T_M[2, 1] * x + T_M[2, 2] * y + T_M[2, 3] * z
  z_new <- T_M[3, 1] * x + T_M[3, 2] * y + T_M[3, 3] * z
  
  # Start with cells in ring  
  bg_spe[["Cell.Type"]] <- ifelse((x_new / (x_radius + ring_width))^2 +
                                    (y_new / (y_radius + ring_width))^2 +
                                    (z_new / (z_radius + ring_width))^2 <= 1,
                                  sample(ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  
  # Then do cells in the cluster  
  bg_spe[["Cell.Type"]] <- ifelse((x_new / x_radius)^2 +
                                    (y_new / y_radius)^2 +
                                    (z_new / z_radius)^2 <= 1,
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  
  # Update current meta data
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- ring_properties
  
  return(bg_spe)
}