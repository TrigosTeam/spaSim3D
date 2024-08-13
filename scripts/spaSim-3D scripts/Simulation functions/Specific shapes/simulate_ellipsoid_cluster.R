library(matlib)

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
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  # Rotation angles
  theta <- cluster_properties$y_z_rotation * (pi/180) # rotation in x-axis
  alpha <- cluster_properties$x_z_rotation * (pi/180) # rotation in y-axis
  beta  <- cluster_properties$x_y_rotation * (pi/180) # rotation in z-axis

  # Get rotation matrices for rotation in the y-z plane (T2), x-z plane (T3) and x-y plane (T4)
  T1 <- matrix(data = c(1, 0, 0,
                        0, cos(theta), -sin(theta),
                        0, sin(theta), cos(theta)), nrow = 3, ncol = 3, byrow = TRUE)
  T2 <- matrix(data = c(cos(alpha), 0, -sin(alpha),
                        0, 1, 0,
                        sin(alpha), 0, cos(alpha)), nrow = 3, ncol = 3, byrow = TRUE)
  T3 <- matrix(data = c(cos(beta), -sin(beta), 0,
                        sin(beta), cos(beta), 0,
                        0, 0, 1), nrow = 3, ncol = 3, byrow = TRUE)
  
  # Get translation matrix from ellipsoid centre (same as centre...)
  T4 <- centre_loc

  ## Change cell types in the ellipsoid cluster
  # Get spatial coords from spe (rows are x, y, z, columns are each cell)
  spe_coords <- t(spatialCoords(bg_spe))
  
  # Apply transformations to spe_coords'
  spe_coords <- inv(T1) %*% inv(T2) %*% inv(T3) %*% (spe_coords - T4)
  x <- spe_coords[1, ]
  y <- spe_coords[2, ]
  z <- spe_coords[3, ]
  
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
