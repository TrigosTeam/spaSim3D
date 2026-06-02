#' @title Simulate ordered background cells in spaSim3D.
#'
#' @description This function simulates background cells in an 'ordered'
#'     fashion, so that cells are geometrically equidistant apart.
#'     The parameters of the background are completely customisable by the user.
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
#' @param jitter_proportion A number between 0 and 1 representing the amount of
#'     'jitter' to apply to all cells. jitter_proportion = 0.5 and the original
#'     distance between each cell is 5 units, then each cell will move UP TO
#'     0.5 * 5 = 2.5 units away from its original position in the x, y and z
#'     directions. The 'original distance between each cell' is calculated by
#'     this function and will depend on the previous parameters, not something
#'     to worry about :)
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
#' bg_n <- simulate_normal_background_cells3D(n_cells = 10000,
#'                                            length = 100,
#'                                            width = 100,
#'                                            height = 100,
#'                                            jitter_proportion = 0,
#'                                            background_cell_type = "Others",
#'                                            plot_image = TRUE)
#'
#' @export

simulate_ordered_background_cells3D <- function(n_cells,
                                                length,
                                                width,
                                                height,
                                                jitter_proportion = 0.25,
                                                background_cell_type = "Others",
                                                plot_image = TRUE) {

  # Check input parameters
  input_parameters <- list("n_cells" = n_cells,
                           "length" = length,
                           "width" = width,
                           "height" = height,
                           "jitter_proportion" = jitter_proportion,
                           "background_cell_type" = background_cell_type,
                           "plot_image" = plot_image)
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))

  # Obtain distance between each point using MAGIC formula
  d_cells <- ((sqrt(2) * length * width * height)/n_cells)^(1/3)

  # Get distance between rows, columns and layers using 'd_cells'
  d_rows <- d_cells
  d_cols <- (sqrt(3) / 2) * d_cells
  d_lays <- (sqrt(6) / 3) * d_cells

  # Get number of rows, columns and layers
  n_rows <- round(length / d_rows)
  n_cols <- round(width / d_cols)
  n_lays <- round(height / d_lays)

  # Step 0. Assume points are on a 3D rectangular grid
  rows <- rep(seq(n_rows), n_cols * n_lays) * d_rows
  cols <- rep(rep(seq(n_cols), each = n_rows), n_lays) * d_cols
  lays <- rep(seq(n_lays), each = n_rows * n_cols) * d_lays

  # Step 1. For every odd sheet, every even row shifts by d_cells/2 right
  if (n_cols %% 2 == 0) {
    shift <- rep(c(rep(0, n_rows), rep(d_cells/2, n_rows)), n_cols/2)
  }
  else {
    shift <- c(rep(c(rep(0, n_rows), rep(d_cells/2, n_rows)), n_cols/2), rep(0, n_rows))
  }
  rows <- rows + c(shift, rep(0, n_rows * n_cols)) # Shift each even row by d_cells/2 right

  # Step 2. For every even sheet, odd rows shift d_cells/2 right, all rows shift d_cells/(2*sqrt(3)) up
  if (n_cols %% 2 == 0) {
    shift <- rep(c(rep(d_cells/2, n_rows), rep(0, n_rows)), n_cols/2)
  }
  else {
    shift <- c(rep(c(rep(d_cells/2, n_rows), rep(0, n_rows)), n_cols/2), rep(d_cells/2, n_rows))
  }
  rows <- rows + c(rep(0, n_rows * n_cols), shift) # Shift each odd row by d_cells/2 right
  cols <- cols + rep(c(0, d_cells/(2 * sqrt(3))), each = n_rows * n_cols) # Shift all rows by d_cells/(2*sqrt(3)) up

  # Get total number of cells (should be roughly equal to n_cells)
  n_total <- n_rows * n_cols * n_lays

  # Add randomness to the location of the cells
  jitter <- jitter_proportion * d_cells # Jitter is proportional to distance between points in hexagonal grid
  jitter_row <- runif(n_total, -jitter, jitter)
  jitter_col <- runif(n_total, -jitter, jitter)
  jitter_lay <- runif(n_total, -jitter, jitter)

  rows <- rows + jitter_row
  cols <- cols + jitter_col
  lays <- lays + jitter_lay

  # Put data into data frame
  df <- data.frame("Cell.X.Position" = rows,
                   "Cell.Y.Position" = cols,
                   "Cell.Z.Position" = lays,
                   "Cell.Type" = background_cell_type)
  df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")

  # Get metadata
  background_metadata <- list("background_type" = "normal",
                              "n_cells" = n_cells,
                              "length" = length,
                              "width" = width,
                              "height" = height,
                              "jitter_proportion" = jitter_proportion,
                              "cell_types" = background_cell_type,
                              "cell_proportions" = 1)
  simulation_metadata <- list(background = background_metadata)

  ## Convert data frame to spe object
  spe <- SpatialExperiment::SpatialExperiment(
    assay = matrix(data = NA, nrow = 0, ncol = nrow(df)),
    colData = df,
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
