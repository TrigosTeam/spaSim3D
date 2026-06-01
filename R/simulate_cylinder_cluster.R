#' @title Simulate a cylinder cluster in spaSim3D.
#'
#' @description This function simulates a cylinder cluster onto an existing
#'     SpatialExperiment object. The parameters of the cylinder are fully
#'     customisable by the user.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for
#'     the cells. Should be generated using the output of one of the following
#'     functions: simulate_random_background_cells3D,
#'     simulate_ordered_background_cells3D,
#'     simulate_mixing3D or any of the other simulate_* functions. This is
#'     because the metadata of the SpatialExperiment object needs to already
#'     contain spaSim3D specific data relating to the background of the
#'     SpatialExperiment object, and any clusters.
#' @param cluster_properties A list containing the properties of the cylinder
#'     cluster desired. The list should contain the following elements:
#'     "shape": Must be equal to the character "cylinder".
#'     "cluster_cell_types": A character vector representing the cell types that
#'         make up the cluster. E.g. c("Tumour", "Immune").
#'     "cluster_cell_proportions": A numeric vector representing the proportion
#'         of each cell type in the cluster. Its elements must each be
#'         greater than 0, sum to 1 and the vector must be the same length as
#'         "cluster_cell_types". E.g. c(0.6, 0.4) corresponds to a cluster made
#'         up of 60% Tumour and 40% Immune.
#'     "radius": A positive number representing the radius of the cylinder.
#'     "start_loc": A numerical vector of length 3 representing the starting
#'         x,y,z coordinate of the cylinder. E.g. (0, 0, 0).
#'     "end_loc": A numerical vector of length 3 representing the ending x,y,z
#'         coordinate of the cylinder. E.g. (20, 20, 20).
#'
#' @return The same 3D SpatialExperiment object used as input for spe, updated
#'     with the new cylinder cluster and the corresponding metadata.
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
#' cylinder_cluster <- simulate_cylinder_cluster(bg_r,
#'                                               cluster_properties = list(
#'                                                 shape = "cylinder",
#'                                                 cluster_cell_types = c("Endothelial", "Others"),
#'                                                 cluster_cell_proportions = c(0.95, 0.05),
#'                                                 radius = 10,
#'                                                 start_loc = c(0, 0, 0),
#'                                                 end_loc   = c(20, 20, 100)
#'                                               ))
#' # Plot
#' plots_cells3D(cylinder_cluster)
#'
#' @export

simulate_cylinder_cluster <- function(spe,
                                      cluster_properties) {

  # Check input parameters
  input_parameters <- cluster_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))

  # Get cylinder properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  radius <- cluster_properties$radius
  start_loc <- cluster_properties$start_loc
  end_loc <- cluster_properties$end_loc

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

  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                  rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                 (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= radius^2), # Cell must be close enough to the cylinder line
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])

  # Update current meta data
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- cluster_properties

  return(spe)
}
