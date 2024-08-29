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
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  
  ## Check number of inner ring cell types matches the number of inner ring cell proportions
  if (length(inner_ring_cell_types) != length(inner_ring_cell_proportions)) stop("Number of inner ring cell types doesn't match number of inner ring cell proportion.")
  
  ## Check inner ring cell proportions are not negative or greater than 1
  if (sum(inner_ring_cell_proportions < 0 | inner_ring_cell_proportions > 1) != 0) stop("Inner ring cell proportions cannot be negative or greater than 1")
  
  ## Check inner ring cell proportions add up to 1
  if (!all.equal(sum(inner_ring_cell_proportions), 1)) stop("Sum of inner ring cell proportions is NOT 1")
  
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  ## Check number of outer ring cell types matches the number of outer ring cell proportions
  if (length(outer_ring_cell_types) != length(outer_ring_cell_proportions)) stop("Number of outer ring cell types doesn't match number of outer ring cell proportion.")
  
  ## Check outer ring cell proportions are not negative or greater than 1
  if (sum(outer_ring_cell_proportions < 0 | outer_ring_cell_proportions > 1) != 0) stop("Outer ring cell proportions cannot be negative or greater than 1")
  
  ## Check outer ring cell proportions add up to 1
  if (!all.equal(sum(outer_ring_cell_proportions), 1)) stop("Sum of outer ring cell proportions is NOT 1")
  
  # Rotation angles
  theta <- dr_properties$y_z_rotation * (pi/180) # rotation in x-axis
  alpha <- dr_properties$x_z_rotation * (pi/180) # rotation in y-axis
  beta  <- dr_properties$x_y_rotation * (pi/180) # rotation in z-axis
  
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
  spe_coords <- solve(T1) %*% solve(T2) %*% solve(T3) %*% (spe_coords - T4)
  x <- spe_coords[1, ]
  y <- spe_coords[2, ]
  z <- spe_coords[3, ]
  
  
  # Start with cells in outer ring  
  bg_spe[["Cell.Type"]] <- ifelse((x / (x_radius + inner_ring_width + outer_ring_width))^2 +
                                    (y / (y_radius + inner_ring_width + outer_ring_width))^2 +
                                    (z / (z_radius + inner_ring_width + outer_ring_width))^2 <= 1,
                                  sample(outer_ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = outer_ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
    
  # Then do cells in inner ring  
  bg_spe[["Cell.Type"]] <- ifelse((x / (x_radius + inner_ring_width))^2 +
                                    (y / (y_radius + inner_ring_width))^2 +
                                    (z / (z_radius + inner_ring_width))^2 <= 1,
                                  sample(inner_ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = inner_ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  
  # Then do cells in the cluster  
  bg_spe[["Cell.Type"]] <- ifelse((x / x_radius)^2 +
                                    (y / y_radius)^2 +
                                    (z / z_radius)^2 <= 1,
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  return(bg_spe)
}