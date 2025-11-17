simulate_ellipsoid_dr <- function(spe, 
                                  dr_properties) {
  
  # Check input parameters
  input_parameters <- dr_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get ellipsoid double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  x_radius <- dr_properties$radii[1]
  y_radius <- dr_properties$radii[2]
  z_radius <- dr_properties$radii[3]
  centre_loc <- dr_properties$centre_loc
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  theta <- dr_properties$axes_rotation[1] * (pi/180) # rotation in x-axis
  alpha <- dr_properties$axes_rotation[2] * (pi/180) # rotation in y-axis
  beta  <- dr_properties$axes_rotation[3] * (pi/180) # rotation in z-axis
  
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
  spe_coords <- t(spatialCoords(spe))
  
  # Apply transformations to spe_coords'
  spe_coords <- solve(T1) %*% solve(T2) %*% solve(T3) %*% (spe_coords - T4)
  x <- spe_coords[1, ]
  y <- spe_coords[2, ]
  z <- spe_coords[3, ]
  
  
  # Start with cells in outer ring  
  spe[["Cell.Type"]] <- ifelse((x / (x_radius + inner_ring_width + outer_ring_width))^2 +
                                 (y / (y_radius + inner_ring_width + outer_ring_width))^2 +
                                 (z / (z_radius + inner_ring_width + outer_ring_width))^2 <= 1,
                               sample(outer_ring_cell_types, size = ncol(spe), replace = TRUE, prob = outer_ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Then do cells in inner ring  
  spe[["Cell.Type"]] <- ifelse((x / (x_radius + inner_ring_width))^2 +
                                 (y / (y_radius + inner_ring_width))^2 +
                                 (z / (z_radius + inner_ring_width))^2 <= 1,
                               sample(inner_ring_cell_types, size = ncol(spe), replace = TRUE, prob = inner_ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  
  # Then do cells in the cluster  
  spe[["Cell.Type"]] <- ifelse((x / x_radius)^2 +
                                 (y / y_radius)^2 +
                                 (z / z_radius)^2 <= 1,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  return(spe)
}
