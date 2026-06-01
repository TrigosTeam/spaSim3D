#' @title Simulate clusters in spaSim3D.
#'
#' @description This function simulates clusters specified by the user onto an
#'     existing SpatialExperiment object. This includes customisable spheres,
#'     ellipsoids, cylinders and network clusters.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for
#'     the cells. Should be generated using the output of one of the following
#'     functions: simulate_random_background_cells3D,
#'     simulate_ordered_background_cells3D,
#'     simulate_mixing3D or any of the other simulate_* functions. This is
#'     because the metadata of the SpatialExperiment object needs to already
#'     contain spaSim3D specific data relating to the background of the
#'     SpatialExperiment object, and any clusters.
#' @param cluster_properties_list A list of lists containing the cluster
#'     properties of each cluster desired. See the example to know what
#'     properties are required for each cluster shape. If you want more detail
#'     on each property, I recommend viewing the other simulate_* functions for
#'     each particular shape.
#' @param plot_image A logical indicating whether to plot 3D spatial data with
#'     the added metadata. Defaults to TRUE.
#' @param plot_cell_types A string vector specifying the cell types to plot. If
#'     NULL, all cell types in the `feature_colname` column will be considered.
#'     Defaults to NULL.
#' @param plot_colours A string vector specifying the colours of the cell types
#'     when plotting. Must match the number of cell types specified in
#'     `plot_cell_types`. If NULL, the viridis color pallete will be used.
#'     Defaults to NULL.
#'
#' @return The same 3D SpatialExperiment object used as input for spe, updated
#'     with the new clusters and corresponding metadata.
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
#' # Simulate clusters
#' clusters <- simulate_clusters3D(bg_r,
#'                                 cluster_properties = list(
#'                                   # Sphere example
#'                                   C1 = list(
#'                                     shape = "sphere",
#'                                     cluster_cell_types = c("Tumour", "Immune", "Others"),
#'                                     cluster_cell_proportions = c(0.55, 0.4, 0.05),
#'                                     radius = 25,
#'                                     centre_loc = c(40, 40, 40)
#'                                   ),
#'                                   # Cylinder example
#'                                   C2 = list(
#'                                     shape = "cylinder",
#'                                     cluster_cell_types = c("Endothelial", "Others"),
#'                                     cluster_cell_proportions = c(0.95, 0.05),
#'                                     radius = 10,
#'                                     start_loc = c(0, 0, 0),
#'                                     end_loc   = c(20, 20 , 100)
#'                                   ),
#'                                   # Ellipsoid example
#'                                   C3 = list(
#'                                     shape = "ellipsoid",
#'                                     cluster_cell_types = c("Tumour", "Immune", "Others"),
#'                                     cluster_cell_proportions = c(0.65, 0.3, 0.05),
#'                                     radii = c(15, 20, 25),
#'                                     centre_loc = c(70, 70, 70),
#'                                     axes_rotation = c(0, 0, 0)
#'                                   ),
#'                                   # Network example
#'                                   C4 = list(
#'                                     shape = "network",
#'                                     cluster_cell_types = c("Immune"),
#'                                     cluster_cell_proportions = c(1),
#'                                     n_edges = 15,
#'                                     width = 8,
#'                                     centre_loc = c(75, 75, 25), # Rough centre of network cluster
#'                                     radius = 40 # Rough radius spanned by the network cluster
#'                                   )
#'                                 ),
#'                                 plot_image = TRUE,
#'                                 plot_cell_types = c("Others", "Immune", "Endothelial", "Tumour"),
#'                                 plot_colours = c("lightgray", "skyblue", "#FF7F7F", "orange"))
#'
#' @export

simulate_clusters3D <- function(spe,
                                cluster_properties_list,
                                plot_image = TRUE,
                                plot_cell_types = NULL,
                                plot_colours = NULL) {

  # Check shape variable of cluster_properties
  shapes <- sapply(cluster_properties_list, function(x) {return(x[["shape"]])})
  n_invalid_shapes <- sum(!(shapes %in% c("sphere", "ellipsoid", "cylinder", "network")))
  if (n_invalid_shapes > 0) {
    stop("`cluster_properties_list` contains invalid shape parameters or no shape parameters.")
  }

  for (i in seq(length(cluster_properties_list))) {

    shape <- shapes[[i]]

    ### Sphere shape
    if (shape == "sphere") {
      spe <- simulate_sphere_cluster(spe, cluster_properties_list[[i]])
    }

    ### Ellipsoid shape
    if (shape == "ellipsoid") {
      spe <- simulate_ellipsoid_cluster(spe, cluster_properties_list[[i]])
    }

    ### Cylinder shape
    if (shape == "cylinder") {
      spe <- simulate_cylinder_cluster(spe, cluster_properties_list[[i]])
    }

    ### Network shape
    if (shape == "network") {
      spe <- simulate_network_cluster(spe, cluster_properties_list[[i]])
    }
  }

  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }

  return(spe)
}
