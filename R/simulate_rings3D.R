#' @title Simulate clusters with rings in spaSim3D.
#'
#' @description This function simulates clusters specified by the user onto an
#'     existing SpatialExperiment object. This includes customisable spheres,
#'     ellipsoids, cylinders and network clusters, all with rings.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for
#'     the cells. Should be generated using the output of one of the following
#'     functions: simulate_random_background_cells3D,
#'     simulate_ordered_background_cells3D,
#'     simulate_mixing3D or any of the other simulate_* functions. This is
#'     because the metadata of the SpatialExperiment object needs to already
#'     contain spaSim3D specific data relating to the background of the
#'     SpatialExperiment object, and any clusters.
#' @param ring_properties_list A list of lists containing the cluster and ring
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
#'     with the new cluster with rings and corresponding metadata.
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
#' # Simulate clusters with rings
#' clusters_with_rings <- simulate_rings3D(bg_r,
#'                                         ring_properties = list(
#'                                           R1 = list(
#'                                             shape = "sphere",
#'                                             cluster_cell_types = c("Tumour", "Others"),
#'                                             cluster_cell_proportions = c(0.95, 0.05),
#'                                             radius = 20,
#'                                             centre_loc = c(40, 40, 40),
#'                                             ring_cell_types = c("Immune", "Others"),
#'                                             ring_cell_proportions = c(0.85, 0.15),
#'                                             ring_width = 5
#'                                           ),
#'                                           R2 = list(
#'                                             shape = "cylinder",
#'                                             cluster_cell_types = c("Void"),
#'                                             cluster_cell_proportions = c(1),
#'                                             radius = 8,
#'                                             start_loc = c(0, 0, 0),
#'                                             end_loc = c(20, 20 , 100),
#'                                             ring_cell_types = c("Endothelial", "Others"),
#'                                             ring_cell_proportions = c(0.85, 0.15),
#'                                             ring_width = 5
#'                                           ),
#'                                           R3 = list(
#'                                             shape = "ellipsoid",
#'                                             cluster_cell_types = c("Tumour", "Others"),
#'                                             cluster_cell_proportions = c(0.95, 0.05),
#'                                             radii = c(20, 30, 40),
#'                                             centre_loc = c(70, 70, 70),
#'                                             axes_rotation = c(0, 0, 0),
#'                                             ring_cell_types = c("Immune", "Others"),
#'                                             ring_cell_proportions = c(0.85, 0.15),
#'                                             ring_width = 5
#'                                           ),
#'                                           R4 = list(
#'                                             shape = "network",
#'                                             cluster_cell_types = c("Tumour"),
#'                                             cluster_cell_proportions = c(1),
#'                                             n_edges = 15,
#'                                             width = 5,
#'                                             centre_loc = c(75, 75, 25),
#'                                             radius = 50,
#'                                             ring_cell_types = c("Endothelial"),
#'                                             ring_cell_proportions = c(1),
#'                                             ring_width = 2
#'                                           )
#'                                         ),
#'                                         plot_image = TRUE,
#'                                         plot_cell_types = c("Others", "Tumour", "Immune", "Endothelial"),
#'                                         plot_colours = c("lightgray", "orange", "skyblue", "tomato"))
#'
#' # Re-plotting so it shows in GitHub pages
#' plot_cells3D(clusters_with_rings,
#'              plot_cell_types = c("Others", "Tumour", "Immune", "Endothelial"),
#'              plot_colours = c("lightgray", "orange", "skyblue", "tomato"))
#'
#' @export

simulate_rings3D <- function(spe,
                             ring_properties_list,
                             plot_image = TRUE,
                             plot_cell_types = NULL,
                             plot_colours = NULL) {

  # Check shape variable of ring_properties_list
  shapes <- sapply(ring_properties_list, function(x) {return(x[["shape"]])})
  n_invalid_shapes <- sum(!(shapes %in% c("sphere", "ellipsoid", "cylinder", "network")))
  if (n_invalid_shapes > 0) {
    stop("`ring_properties_list` contains invalid shape parameters or no shape parameters.")
  }

  for (i in seq(length(ring_properties_list))) {

    shape <- shapes[[i]]

    ### Sphere shape with ring
    if (shape == "sphere") {
      spe <- simulate_sphere_ring(spe, ring_properties_list[[i]])
    }

    ### Ellipsoid shape with ring
    else if (shape == "ellipsoid") {
      spe <- simulate_ellipsoid_ring(spe, ring_properties_list[[i]])
    }

    ### Cylinder shape with ring
    else if (shape == "cylinder") {
      spe <- simulate_cylinder_ring(spe, ring_properties_list[[i]])
    }

    ### Network shape with ring
    else if (shape == "network") {
      spe <- simulate_network_ring(spe, ring_properties_list[[i]])
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
