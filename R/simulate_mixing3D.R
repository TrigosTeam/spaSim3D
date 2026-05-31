#' @title Simulate mixing in spaSim3D.
#'
#' @description This function simulates mixing of an existing SpatialExperiment 
#'     object. This means completely updating the cell types that make up the
#'     background of the SpatialExperiment object.
#' 
#' @param spe A SpatialExperiment object containing 3D spatial information for 
#'     the cells. Should be generated using the output of one of the following
#'     functions: simulate_random_background_cells3D, 
#'     simulate_ordered_background_cells3D,
#'     simulate_mixing3D or any of the other simulate_* functions. This is 
#'     because the metadata of the SpatialExperiment object needs to already 
#'     contain spaSim3D specific data relating to the background of the 
#'     SpatialExperiment object, and any clusters.
#' @param cell_types A character vector representing the cell types that will
#'         make up the background of the SpatialExperiment object. E.g. 
#'         c("Tumour", "Immune").
#' @param cell_proportions  A numeric vector representing the proportion 
#'         of each cell type in the background of the SpatialExperiment object. 
#'         Its elements must each be greater than 0, sum to 1 and the vector 
#'         must be the same length as "cell_types". E.g. c(0.6, 0.4) corresponds 
#'         to a background made up of 60% Tumour and 40% Immune.
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
#'     with the new mixing of the background and corresponding metadata.
#'
#' @examples
#' # Simulate background
#' bg_n <- simulate_normal_background_cells3D(n_cells = 10000,
#'                                            length = 100,
#'                                            width = 100,
#'                                            height = 100,
#'                                            jitter_proportion = 0,
#'                                            background_cell_type = "Others",
#'                                            plot_image = FALSE)
#'                                            
#' # Simulate mixing
#' bg_mix <- simulate_mixing3D(bg_n,
#'                             cell_types = c("Others", "Immune", "Tumour"),
#'                             cell_proportions = c(0.5, 0.25, 0.25),
#'                             plot_image = TRUE,
#'                             plot_cell_types = c("Others", "Immune", "Tumour"),
#'                             plot_colours = c("lightgray", "skyblue", "orange"))
#'                                                      
#' @export

simulate_mixing3D <- function(spe,
                              cell_types,
                              cell_proportions,
                              plot_image = TRUE,
                              plot_cell_types = NULL,
                              plot_colours = NULL) {
  
  # Check input parameters
  input_parameters <- list("spe" = spe,
                           "cell_types" = cell_types,
                           "cell_proportions" = cell_proportions,
                           "plot_image" = plot_image)
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Apply mixing
  spe[["Cell.Type"]] <- sample(cell_types, size = ncol(spe), replace = TRUE, prob = cell_proportions)
  
  spe@metadata[["simulation"]][["background"]][["cell_types"]] <- cell_types
  spe@metadata[["simulation"]][["background"]][["cell_proportions"]] <- cell_proportions
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(spe)
}
