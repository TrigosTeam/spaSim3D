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
