#' @title Simulate random background cells in spaSim3D.
#'
#' @description This function simulates background cells in a 'random' fashion,
#'     so that cells follows a Poisson distribution. The parameters of the
#'     background are completely customisable by the user.
#'
#' @param n_cells A positive number representing the number of cells in the
#'     background.
#' @param length A positive number representing the length of the 3D window of
#'     the background. E.g. If you want a 100 by 150 by 200 unit window, set
#'     length to 100.
#' @param width A positive number representing the width of the 3D window of
#'     the background. E.g. If you want a 100 by 150 by 200 unit window, set
#'     width to 150.
#' @param height A positive number representing the height of the 3D window of
#'     the background. E.g. If you want a 100 by 150 by 200 unit window, set
#'     height to 200.
#' @param minimum_distance_between_cells A positive number representing the
#'     minimum distance between cells.
#' @param background_cell_type A character representing the cell type label of
#'     all the cells in the background. Defaults to "Others".
#' @param plot_image A logical indicating whether to plot 3D spatial data.
#'     Defaults to TRUE.
#'
#' @return A 3D SpatialExperiment object with the background cells and
#'     corresponding metadata.
#'
#' @examples
#' # Simulate background
#' bg_r <- simulate_random_background_cells3D(n_cells = 10000,
#'                                            length = 100,
#'                                            width = 100,
#'                                            height = 100,
#'                                            minimum_distance_between_cells = 0.5,
#'                                            background_cell_type = "Others",
#'                                            plot_image = TRUE)
#'
#' # Re-plotting so it shows in GitHub pages
#' plot_cells3D(bg_r,
#'              plot_cell_types = "Others",
#'              plot_colours = "lightgray")
#'
#' @export

simulate_random_background_cells3D <- function(n_cells,
                                               length,
                                               width,
                                               height,
                                               minimum_distance_between_cells,
                                               background_cell_type = "Others",
                                               plot_image = TRUE) {

  # Check input parameters
  input_parameters <- list("n_cells" = n_cells,
                           "length" = length,
                           "width" = width,
                           "height" = height,
                           "minimum_distance_between_cells" = minimum_distance_between_cells,
                           "background_cell_type" = background_cell_type,
                           "plot_image" = plot_image)
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))

  # Need to over-sample as cells which are too close will be removed later
  n_cells_inflated <- n_cells * 2

  spatial_df <- data.frame("Cell.X.Position" = runif(n_cells_inflated, 0, length),
                           "Cell.Y.Position" = runif(n_cells_inflated, 0, width),
                           "Cell.Z.Position" = runif(n_cells_inflated, 0, height))

  # Add integer rownames to data frame - each cell is labelled by an integer
  rownames(spatial_df) <- seq(nrow(spatial_df))

  ### Check if all other cells are to close to the current cell
  # Use frNN function: for each point, get all points within min_d of it
  spatial_df_distances <- dbscan::frNN(spatial_df,
                                       eps = minimum_distance_between_cells,
                                       query = NULL,
                                       sort = FALSE)

  # For each cell, get all other cells which were within 'minimum_distance_between_cells'
  spatial_df_distances_ids <- spatial_df_distances$id

  # Filter out zero length cells
  spatial_df_distances_ids <- Filter(function(x) length(x) != 0, spatial_df_distances_ids)

  # Get integer labels for the remaining cells
  spatial_df_distances_ids_names <- as.integer(names(spatial_df_distances_ids))

  # Determine which cells should be chosen from spatial_df
  cells_chosen <- rep(T, nrow(spatial_df))
  for (i in seq_len(length(spatial_df_distances_ids))) {
    cells_too_close <- spatial_df_distances_ids[[i]]

    if (cells_chosen[spatial_df_distances_ids_names[i]]) cells_chosen[cells_too_close] <- F
  }

  spatial_df <- spatial_df[cells_chosen, ]

  # If number of cells remaining is still higher than n_cells, randomly subset n_cells
  if (nrow(spatial_df) > n_cells) {
    spatial_df <- dplyr::sample_n(spatial_df, n_cells)
  }
  else {
    warning("The number of cells in your object could not reach the desired number of cells, consider lowering the number of cells, the minimum distance between cells, or increasing length, width or height.")
  }

  # Add Cell.Type and Cell.ID
  spatial_df$Cell.Type <- background_cell_type
  spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")

  # Get metadata
  background_metadata <- list("background_type" = "random",
                              "n_cells" = n_cells,
                              "length" = length,
                              "width" = width,
                              "height" = height,
                              "minimum_distance_between_cells" = minimum_distance_between_cells,
                              "cell_types" = background_cell_type,
                              "cell_proportions" = 1)
  simulation_metadata <- list(background = background_metadata)

  ## Convert data frame to spe object
  spe <- SpatialExperiment::SpatialExperiment(
    assay = matrix(data = NA, nrow = 0, ncol = nrow(spatial_df)),
    colData = spatial_df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = list(simulation = simulation_metadata))

  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        background_cell_type,
                        "lightgray")
    methods::show(fig)
  }

  return(spe)
}
