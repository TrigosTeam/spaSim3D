#' @title Simulate an ellipsoid cluster with a double ring in spaSim3D.
#'
#' @description This function simulates an ellipsoid cluster with a double ring
#'     onto an existing SpatialExperiment object. The parameters of the
#'     ellipsoid and its double ring are fully customisable by the user.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for
#'     the cells. Should be generated using the output of one of the following
#'     functions: simulate_random_background_cells3D,
#'     simulate_ordered_background_cells3D,
#'     simulate_mixing3D or any of the other simulate_* functions. This is
#'     because the metadata of the SpatialExperiment object needs to already
#'     contain spaSim3D specific data relating to the background of the
#'     SpatialExperiment object, and any clusters.
#' @param dr_properties A list containing the properties of the ellipsoid
#'     cluster and double ring desired. The list should contain the following
#'     elements:
#'     "shape": Must be equal to the character "ellipsoid".
#'     "cluster_cell_types": A character vector representing the cell types that
#'         make up the cluster. E.g. c("Tumour", "Immune").
#'     "cluster_cell_proportions": A numeric vector representing the proportion
#'         of each cell type in the cluster. Its elements must each be
#'         greater than 0, sum to 1 and the vector must be the same length as
#'         "cluster_cell_types". E.g. c(0.6, 0.4) corresponds to a cluster made
#'         up of 60% Tumour and 40% Immune.
#'     "radii": A numerical vector of length 3 containing only positive numbers,
#'         representing the radii of the ellipsoid in the x,y,z direction
#'         respectively. E.g. c(15, 20, 25).
#'     "centre_loc": A numerical vector of length 3 representing the centre
#'         x,y,z coordinate of the ellipsoid E.g. (0, 0, 0).
#'     "axes_rotation": A numerical vector of length 3 representing axes
#'         rotation of the ellipsoid in the y-z plane, x-z plane and x-y plane
#'         respectively. Values should be in degrees unit. E.g. (30, 45, 0).
#'     "inner_ring_cell_types": A character vector representing the cell types
#'         that make up the inner ring. E.g. c("Immune1", "Immune2", "Immune3").
#'     "inner_ring_cell_proportions": A numeric vector representing the
#'         proportion of each cell type in the inner ring. Its elements must
#'         each be greater than 0, sum to 1 and the vector must be the same
#'         length as "inner_ring_cell_types". E.g. c(0.3, 0.4, 0.3) corresponds
#'         to an inner ring made up of 30% Immune1, 40% Immune2 and 30% Immune3.
#'     "inner_ring_width": A positive number representing the width of the inner
#'         ring.
#'     "outer_ring_cell_types": A character vector representing the cell types
#'         that make up the outer ring. E.g. c("T cell", "B cell").
#'     "outer_ring_cell_proportions": A numeric vector representing the
#'         proportion of each cell type in the outer ring. Its elements must
#'         each be greater than 0, sum to 1 and the vector must be the same
#'         length as "outer_ring_cell_types". E.g. c(0.5, 0.5) corresponds to an
#'         outer ring made up of 50% T cell and 50% B cell.
#'     "outer_ring_width": A positive number representing the width of the outer
#'         ring.
#'
#' @return The same 3D SpatialExperiment object used as input for spe, updated
#'     with the new ellipsoid with double ring cluster and the corresponding
#'     metadata.
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
#' ellipsoid_dr_cluster <- simulate_ellipsoid_dr(bg_r,
#'                                               dr_properties = list(
#'                                                 shape = "ellipsoid",
#'                                                 cluster_cell_types = c("Tumour", "Others"),
#'                                                 cluster_cell_proportions = c(0.95, 0.05),
#'                                                 radii = c(10, 15, 20),
#'                                                 centre_loc = c(70, 70, 70),
#'                                                 axes_rotation = c(0, 0, 0),
#'                                                 inner_ring_cell_types = c("Immune1", "Others"),
#'                                                 inner_ring_cell_proportions = c(0.85, 0.15),
#'                                                 inner_ring_width = 5,
#'                                                 outer_ring_cell_types = c("Immune2"),
#'                                                 outer_ring_cell_proportions = c(1),
#'                                                 outer_ring_width = 3
#'                                               ))
#'
#' # Plot
#' plots_cells3D(ellipsoid_dr_cluster)
#'
#' @export

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
  spe_coords <- t(SpatialExperiment::spatialCoords(spe))

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
