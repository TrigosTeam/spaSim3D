#' @title Create spaSim3D cluster metadata template.
#'
#' @description This function creates a spaSim3D cluster metadata template with
#'     the predefined parameters required. Users can change the parameters after
#'     creating the template. Can be used as input for the following functions:
#'     spe_metadata_background_template, spe_metadata_cluster_template,
#'     simulate_spe_metadata3D, add_spe_metadata3D.
#'
#' @param cluster_type Either "regular", "ring" or "double ring". These are the
#'     three types of cluster types available in spaSim3D, each with their own
#'     set of parameters.
#' @param shape Either "sphere", "ellipsoid", "cylinder" or "network". These are
#'     the four types of cluster shapes available in spaSim3D, each with their
#'     own set of parameters.
#' @param original_spe_metadata A list containing the metadata for a spaSim3D
#'     SpatialExperiment object. This is generated from the output of the
#'     spaSim3D function 'spe_metadata_background_template',
#'     'spe_metadata_cluster_template'. Defaults to NULL. If original metadata
#'     is inputted, the cluster metadata template will be added to it. If NULL,
#'     the cluster metadata template will be returned, without a background
#'     element. Hence, it is permissible to generate metadata that does not
#'     include a background element. This can be inputted into the
#'     spe_metadata_background_template to add a background element, however,
#'     cannot be inputed into simulate_spe_metadata3D function as this function
#'     requires a background element.
#'
#' @return A list containing the cluster metadata template or the original
#'     metadata updated with the cluster metadata template.
#'
#' @examples
#' # Get background metadata
#' bg_metadata <- spe_metadata_background_template("random")
#'
#' # Change background metadata
#' bg_metadata$background$n_cells <- 25000
#'
#' # Cluster metadata (using background metadata as background)
#' cluster_metadata <- spe_metadata_cluster_template("regular", "sphere", bg_metadata)
#' cluster_metadata <- spe_metadata_cluster_template("ring", "ellipsoid", cluster_metadata)
#' cluster_metadata <- spe_metadata_cluster_template("double ring", "cylinder", cluster_metadata)
#'
#' # Change cluster metadata
#' cluster_metadata$cluster_1$radius <- 120
#'
#' # Get spe from updated metadata
#' spe_clusters <- simulate_spe_metadata3D(cluster_metadata)
#' plot_cells3D(spe_clusters,
#'              plot_cell_types = c("Others", "Tumour", "Immune", "Immune1", "Endothelial"),
#'              plot_colours = c("lightgray", "orange", "skyblue", "lightgreen", "tomato"))
#'
#' @export

spe_metadata_cluster_template <- function(cluster_type,
                                          shape,
                                          original_spe_metadata = NULL) {

  ### Get template for different shapes
  if (shape == "sphere") {
    cluster_metadata <- list(shape = "sphere",
                             cluster_cell_types = c("Tumour", "Immune", "Others"),
                             cluster_cell_proportions = c(0.8, 0.15, 0.05),
                             radius = 100,
                             centre_loc = c(200, 150, 200))
  }
  else if (shape == "ellipsoid") {
    cluster_metadata <- list(shape = "ellipsoid",
                             cluster_cell_types = c("Tumour", "Immune", "Others"),
                             cluster_cell_proportions = c(0.8, 0.15, 0.05),
                             radii = c(75, 100, 125),
                             centre_loc = c(450, 300, 100),
                             axes_rotation = c(0, 45, 0))
  }
  else if (shape == "cylinder") {
    cluster_metadata <- list(shape = "cylinder",
                             cluster_cell_types = c("Endothelial", "Others"),
                             cluster_cell_proportions = c(0.95, 0.05),
                             radius = 40,
                             start_loc = c(400, 0, 0),
                             end_loc   = c(600, 400, 200))
  }
  else if (shape == "network") {
    cluster_metadata <- list(shape = "network",
                             cluster_cell_types = c("Immune", "Others"),
                             cluster_cell_proportions = c(0.95, 0.05),
                             n_edges = 20,
                             width = 30,
                             centre_loc = c(200, 400, 150),
                             radius = 200)
  }
  else {
    stop("shape parameter must be 'sphere', 'ellipsoid', 'cylinder' or 'network'")
  }

  ### Add extra metadata for different cluster types
  if (cluster_type == "regular") {
    cluster_metadata <- append(list(cluster_type = "regular"), cluster_metadata)
  }
  else if (cluster_type == "ring") {
    cluster_metadata <- append(list(cluster_type = "ring"), cluster_metadata)
    cluster_metadata$ring_cell_types <- c("Immune1", "Others")
    cluster_metadata$ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$ring_width <- 12
  }
  else if (cluster_type == "double ring") {
    cluster_metadata <- append(list(cluster_type = "double ring"), cluster_metadata)
    cluster_metadata$inner_ring_cell_types <- c("Immune1", "Others")
    cluster_metadata$inner_ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$inner_ring_width <- 10
    cluster_metadata$outer_ring_cell_types <- c("Immune2", "Others")
    cluster_metadata$outer_ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$outer_ring_width <- 10
  }
  else {
    stop("cluster_type parameter must be 'regular', 'ring' or 'double ring'")
  }

  # If original_spe_metadata input is not null, add new cluster_metadata to it
  if (!is.null(original_spe_metadata) && !is.null(original_spe_metadata[["background"]])) {
    original_spe_metadata[[paste("cluster", length(original_spe_metadata), sep="_")]] <- cluster_metadata
    return(original_spe_metadata)
  }
  else if (!is.null(original_spe_metadata) && is.null(original_spe_metadata[["background"]])) {
    original_spe_metadata[[paste("cluster", length(original_spe_metadata) + 1, sep="_")]] <- cluster_metadata
    return(original_spe_metadata)
  }

  # Else, just return the new cluster_metadata
  return(list("cluster_1" = cluster_metadata))
}
