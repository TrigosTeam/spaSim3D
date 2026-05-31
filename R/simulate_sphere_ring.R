#' @title Simulate a sphere cluster with a ring in spaSim3D.
#'
#' @description This function simulates a sphere cluster with a ring onto an 
#'     existing SpatialExperiment object. The parameters of the sphere and its 
#'     ring are fully customisable by the user.
#' 
#' @param spe A SpatialExperiment object containing 3D spatial information for 
#'     the cells. Should be generated using the output of one of the following
#'     functions: simulate_random_background_cells3D, 
#'     simulate_ordered_background_cells3D,
#'     simulate_mixing3D or any of the other simulate_* functions. This is 
#'     because the metadata of the SpatialExperiment object needs to already 
#'     contain spaSim3D specific data relating to the background of the 
#'     SpatialExperiment object, and any clusters.
#' @param ring_properties A list containing the properties of the sphere 
#'     cluster and ring desired. The list should contain the following elements:
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
#'     "ring_cell_types": A character vector representing the cell types 
#'         that make up the ring. E.g. c("Immune1", "Immune2", "Immune3").
#'     "ring_cell_proportions": A numeric vector representing the 
#'         proportion of each cell type in the ring. Its elements must each be 
#'         greater than 0, sum to 1 and the vector must be the same length as
#'         "ring_cell_types". E.g. c(0.3, 0.4, 0.3) corresponds to an
#'         ring made up of 30% Immune1, 40% Immune2 and 30% Immune3.
#'     "ring_width": A positive number representing the width of the ring.
#'
#' @return The same 3D SpatialExperiment object used as input for spe, updated
#'     with the new sphere with ring cluster and the corresponding metadata.
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
#' sphere_ring_cluster <- simulate_cylinder_ring(bg_r,
#'                                               ring_properties = list(
#'                                                 shape = "sphere",
#'                                                 cluster_cell_types = c("Tumour", "Others"),
#'                                                 cluster_cell_proportions = c(0.95, 0.05),
#'                                                 radius = 20,
#'                                                 centre_loc = c(40, 40, 40),
#'                                                 ring_cell_types = c("Immune", "Others"),
#'                                                 ring_cell_proportions = c(0.85, 0.15),
#'                                                 ring_width = 5
#'                                               ))
#' # Plot
#' plots_cells3D(sphere_ring_cluster)
#'                                             
#' @export

simulate_sphere_ring <- function(spe, 
                                 ring_properties) {
  
  # Check input parameters
  input_parameters <- ring_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get sphere ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  radius <- ring_properties$radius
  centre_loc <- ring_properties$centre_loc
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  ## Change cell types in the sphere ringed cluster
  spe_coords <- data.frame(spatialCoords(spe))
  
  # Start with cells in ring  
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= (radius + ring_width)^2,
                               sample(ring_cell_types, size = ncol(spe), replace = TRUE, prob = ring_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Then do cells in the cluster 
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                 (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                 (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= radius^2,
                               sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                               spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- ring_properties
  
  return(spe)
}
