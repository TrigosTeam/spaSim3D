#' @title Simulate a sphere cluster in spaSim3D.
#'
#' @description This functions simulates a sphere cluster onto an existing 
#'     SpatialExperiment object. The parameters of the sphere are fully
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
#' @param cluster_properties A list containing the properties of the sphere 
#'     cluster desired. The list should contain the following elements:
#'     "shape": Must be equal to the character "sphere".
#'     "cluster_cell_types": A character vector representing the cell types that 
#'         make up the cluster. E.g. c("Tumour", "Immune").
#'     "cluster_cell_proportions": A numeric vector representing the proportion 
#'         of each cell type in the cluster. Its elements must each be 
#'         greater than 0, sum to 1 and the vector must be the same length as
#'         "cluster_cell_types". E.g. c(0.6, 0.4) corresponds to a cluster made
#'         up of 60% Tumour and 40% Immune.
#'     "radius": A positive number representing the radius of the sphere.
#'     "centre_loc": A numerical vector of length 3 representing the centre
#'         x,y,z coordinate of the sphere E.g. (40, 50, 60).
#'
#' @return The same 3D SpatialExperiment object used as input for spe, updated
#'     with the new sphere cluster and the corresponding metadata.
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
#' sphere_cluster <- simulate_sphere_cluster(bg_r,
#'                                           cluster_properties = list(
#'                                             shape = "sphere",
#'                                             cluster_cell_types = c("Tumour", "Immune", "Others"),
#'                                             cluster_cell_proportions = c(0.55, 0.4, 0.05),
#'                                             radius = 25,
#'                                             centre_loc = c(40, 40, 40)
#'                                           ))
#' # Plot
#' plots_cells3D(sphere_cluster)
#'                                            
#' @export

simulate_sphere_cluster <- function(spe, 
                                    cluster_properties) {
  
  # Check input parameters
  input_parameters <- cluster_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get sphere properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  radius <- cluster_properties$radius
  centre_loc <- cluster_properties$centre_loc
  
  # Change cell types in the sphere cluster
  spe_coords <- data.frame(spatialCoords(spe))
  
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= radius^2,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- cluster_properties
  
  return(spe)
}
