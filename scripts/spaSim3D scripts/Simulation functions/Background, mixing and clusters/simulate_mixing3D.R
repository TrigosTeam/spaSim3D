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
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(bg_spe), 
                   "Cell.Type" = bg_spe[["Cell.Type"]],
                   "Cell.ID" = bg_spe[["Cell.ID"]])
  
  n_cell_types <- length(cell_types)
  
  for (i in 1:nrow(df)) {
    x <- df$Cell.X.Position[i]
    y <- df$Cell.Y.Position[i]
    z <- df$Cell.Z.Position[i]
    
    # Random number will determine the cell_type of the cell
    random <- runif(n = 1, min = 0, max = 1)
    
    # Start with the first cell
    n <- 1 
    current_proportion <- 0
    
    while (n <= n_cell_types){
      current_proportion <- current_proportion + cell_proportions[n]
      if (random <= current_proportion) {
        chosen_cell_type <- cell_types[n]
        break
      }
      n <- n + 1
    }
    df[i, "Cell.Type"] <- chosen_cell_type
  }
  
  # Get meta data
  metadata <- bg_spe@metadata
  metadata[["background"]][["cell_types"]] <- cell_types
  metadata[["background"]][["cell_proportions"]] <- cell_proportions
  
  # Convert data frame to spe object
  mixed_spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = metadata)
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(mixed_spe,
                        plot_cell_types,
                        plot_colours)
    print(fig)
  }
    
  return(mixed_spe)
}
