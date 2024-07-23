simulate_mixing3D <- function(bg_spe,
                              cell_types,
                              cell_proportions,
                              plot_image = TRUE,
                              plot_cell_types = NULL,
                              plot_colours = NULL) {
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cell_types) != length(cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cell_proportions < 0 | cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (sum(cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
  
  
  bg_spe[["Cell.Type"]] <- sample(cell_types, size = ncol(bg_spe), replace = TRUE, prob = cell_proportions)
  
  bg_spe@metadata[["background"]][["cell_types"]] <- cell_types
  bg_spe@metadata[["background"]][["cell_proportions"]] <- cell_proportions
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(bg_spe,
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
    
  return(bg_spe)
}