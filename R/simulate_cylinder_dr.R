#' @title Simulate a cylinder cluster with a double ring in spaSim3D.
#'
#' @description This function simulates a cylinder cluster with a double ring
#'     onto an existing SpatialExperiment object. The parameters of the cylinder
#'     and its double ring are fully customisable by the user.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for
#'     the cells. Should be generated using the output of one of the following
#'     functions: simulate_random_background_cells3D,
#'     simulate_ordered_background_cells3D,
#'     simulate_mixing3D or any of the other simulate_* functions. This is
#'     because the metadata of the SpatialExperiment object needs to already
#'     contain spaSim3D specific data relating to the background of the
#'     SpatialExperiment object, and any clusters.
#' @param dr_properties A list containing the properties of the cylinder cluster
#'     and double ring desired. The list should contain the following elements:
#'
#'     "shape": Must be equal to the character "cylinder".
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
#'     "radius": A positive number representing the radius of the cylinder.
#'
#'     "start_loc": A numerical vector of length 3 representing the starting
#'         x,y,z coordinate of the cylinder. E.g. (0, 0, 0).
#'
#'     "end_loc": A numerical vector of length 3 representing the ending x,y,z
#'         coordinate of the cylinder. E.g. (20, 20, 20).
#'
#'     "inner_ring_cell_types": A character vector representing the cell types
#'         that make up the inner ring. E.g. c("Immune1", "Immune2", "Immune3").
#'
#'     "inner_ring_cell_proportions": A numeric vector representing the
#'         proportion of each cell type in the inner ring. Its elements must
#'         each be greater than 0, sum to 1 and the vector must be the same
#'         length as "inner_ring_cell_types". E.g. c(0.3, 0.4, 0.3) corresponds
#'         to an inner ring made up of 30% Immune1, 40% Immune2 and 30% Immune3.
#'
#'     "inner_ring_width": A positive number representing the width of the inner
#'         ring.
#'
#'     "outer_ring_cell_types": A character vector representing the cell types
#'         that make up the outer ring. E.g. c("T cell", "B cell").
#'
#'     "outer_ring_cell_proportions": A numeric vector representing the
#'         proportion of each cell type in the outer ring. Its elements must
#'         each be greater than 0, sum to 1 and the vector must be the same
#'         length as "outer_ring_cell_types". E.g. c(0.5, 0.5) corresponds to an
#'         outer ring made up of 50% T cell and 50% B cell.
#'
#'     "outer_ring_width": A positive number representing the width of the outer
#'         ring.
#'
#' @return The same 3D SpatialExperiment object used as input for spe, updated
#'     with the new cylinder with double ring cluster and the corresponding
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
#' cylinder_dr_cluster <- simulate_cylinder_dr(bg_r,
#'                                             dr_properties = list(
#'                                               shape = "cylinder",
#'                                               cluster_cell_types = c("Endothelial", "Others"),
#'                                               cluster_cell_proportions = c(0.95, 0.05),
#'                                               radius = 10,
#'                                               start_loc = c(0, 0, 0),
#'                                               end_loc   = c(20, 20, 100),
#'                                               inner_ring_cell_types = c("Immune1", "Others"),
#'                                               inner_ring_cell_proportions = c(0.85, 0.15),
#'                                               inner_ring_width = 5,
#'                                               outer_ring_cell_types = c("Immune2"),
#'                                               outer_ring_cell_proportions = c(1),
#'                                               outer_ring_width = 3
#'                                             ))
#' # Plot
#' plot_cells3D(cylinder_dr_cluster)
#'
#' @export

simulate_cylinder_dr <- function(spe,
                                 dr_properties) {

  # Check input parameters
  input_parameters <- dr_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))

  # Get cylinder double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  radius <- dr_properties$radius
  start_loc <- dr_properties$start_loc
  end_loc <- dr_properties$end_loc
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width

  ## Check if start and end coordinates of the cylinder are the same
  if (identical(start_loc, end_loc)) warning("Start and end coordinates of the cylinder are the same.")

  ## Change cell types in the cylinder cluster
  spe_coords <- SpatialExperiment::spatialCoords(spe)

  # Get directional vector
  v1 <- end_loc - start_loc

  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)

  # Get vector between from each cell to start_loc
  v2 <- sweep(spe_coords, 2, end_loc, "-")

  # Start with cells in outer ring
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= (radius + inner_ring_width + outer_ring_width)^2), # Cell must be close enough to the cylinder line
                               sample(outer_ring_cell_types, size = ncol(spe), replace = TRUE, prob = outer_ring_cell_proportions),
                               spe[["Cell.Type"]])

  # Start with cells in inner ring
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= (radius + inner_ring_width)^2), # Cell must be close enough to the cylinder line
                               sample(inner_ring_cell_types, size = ncol(spe), replace = TRUE, prob = inner_ring_cell_proportions),
                               spe[["Cell.Type"]])

  # Then do cells in the cluster
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= radius^2), # Cell must be close enough to the cylinder line
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])

  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- dr_properties

  return(spe)
}
