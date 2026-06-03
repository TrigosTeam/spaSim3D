#' @title Simulate an ellipsoid cluster with a ring in spaSim3D.
#'
#' @description This function simulates an ellipsoid cluster with a ring onto an
#'     existing SpatialExperiment object. The parameters of the ellipsoid and
#'     its ring are fully customisable by the user.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for
#'     the cells. Should be generated using the output of one of the following
#'     functions: simulate_random_background_cells3D,
#'     simulate_ordered_background_cells3D,
#'     simulate_mixing3D or any of the other simulate_* functions. This is
#'     because the metadata of the SpatialExperiment object needs to already
#'     contain spaSim3D specific data relating to the background of the
#'     SpatialExperiment object, and any clusters.
#' @param ring_properties A list containing the properties of the ellipsoid
#'     cluster and ring desired. The list should contain the following elements:
#'
#'     "shape": Must be equal to the character "ellipsoid".
#'
#'     "cluster_cell_types": A character vector representing the cell types that
#'         make up the cluster. E.g. c("Tumour", "Immune").
#'
#'     "cluster_cell_proportions": A numeric vector representing the proportion
#'         of each cell type in the cluster. Its elements must each be
#'         greater than 0, sum to 1 and the vector must be the same length as
#'         "cluster_cell_types". E.g. c(0.6, 0.4) corresponds to a cluster made
#'         up of 60% Tumour and 40% Immune.
#'
#'     "radii": A numerical vector of length 3 containing only positive numbers,
#'         representing the radii of the ellipsoid in the x,y,z direction
#'         respectively. E.g. c(15, 20, 25).
#'
#'     "centre_loc": A numerical vector of length 3 representing the centre
#'         x,y,z coordinate of the ellipsoid E.g. (0, 0, 0).
#'
#'     "axes_rotation": A numerical vector of length 3 representing axes
#'         rotation of the ellipsoid in the y-z plane, x-z plane and x-y plane
#'         respectively. Values should be in degrees unit. E.g. (30, 45, 0).
#'
#'     "ring_cell_types": A character vector representing the cell types
#'         that make up the ring. E.g. c("Immune1", "Immune2", "Immune3").
#'
#'     "ring_cell_proportions": A numeric vector representing the
#'         proportion of each cell type in the ring. Its elements must each be
#'         greater than 0, sum to 1 and the vector must be the same length as
#'         "ring_cell_types". E.g. c(0.3, 0.4, 0.3) corresponds to an
#'         ring made up of 30% Immune1, 40% Immune2 and 30% Immune3.
#'
#'     "ring_width": A positive number representing the width of the ring.
#'
#' @return The same 3D SpatialExperiment object used as input for spe, updated
#'     with the new ellipsoid with ring cluster and the corresponding metadata.
#'
#' @examples
#' # Simulate background
#' bg_r <- simulate_random_background_cells3D(n_cells = 10000,
#'                                            length = 100,
#'                                            width = 100,
#'                                            height = 100,
#'                                            minimum_distance_between_cells = 0.5,
#'                                            background_cell_type = "Others",
#'                                            plot_image = FALSE)
#'
#' # Simulate cluster
#' ellipsoid_ring_cluster <- simulate_ellipsoid_ring(bg_r,
#'                                                   ring_properties = list(
#'                                                     shape = "ellipsoid",
#'                                                     cluster_cell_types = c("Tumour", "Others"),
#'                                                     cluster_cell_proportions = c(0.95, 0.05),
#'                                                     radii = c(20, 30, 40),
#'                                                     centre_loc = c(70, 70, 70),
#'                                                     axes_rotation = c(0, 0, 0),
#'                                                     ring_cell_types = c("Immune", "Others"),
#'                                                     ring_cell_proportions = c(0.85, 0.15),
#'                                                     ring_width = 5
#'                                                   ))
#' # Plot
#' plot_cells3D(ellipsoid_ring_cluster)
#'
#' @export

simulate_ellipsoid_ring <- function(spe,
                                    ring_properties) {

  # Check input parameters
  input_parameters <- ring_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))

  # Get ellipsoid ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  x_radius <- ring_properties$radii[1]
  y_radius <- ring_properties$radii[2]
  z_radius <- ring_properties$radii[3]
  centre_loc <- ring_properties$centre_loc
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  theta <- ring_properties$axes_rotation[1] * (pi/180) # rotation in x-axis
  alpha <- ring_properties$axes_rotation[2] * (pi/180) # rotation in y-axis
  beta  <- ring_properties$axes_rotation[3] * (pi/180) # rotation in z-axis

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
  spe_coords <- t(SpatialExperiment::spatialCoords(spe))

  # Apply transformations to spe_coords'
  spe_coords <- solve(T1) %*% solve(T2) %*% solve(T3) %*% (spe_coords - T4)
  x <- spe_coords[1, ]
  y <- spe_coords[2, ]
  z <- spe_coords[3, ]

  # Start with cells in ring
  spe[["Cell.Type"]] <- ifelse((x / (x_radius + ring_width))^2 +
                                 (y / (y_radius + ring_width))^2 +
                                 (z / (z_radius + ring_width))^2 <= 1,
                               sample(ring_cell_types, size = ncol(spe), replace = TRUE, prob = ring_cell_proportions),
                               spe[["Cell.Type"]])


  # Then do cells in the cluster
  spe[["Cell.Type"]] <- ifelse((x / x_radius)^2 +
                                 (y / y_radius)^2 +
                                 (z / z_radius)^2 <= 1,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])


  # Update current meta data
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- ring_properties

  return(spe)
}
