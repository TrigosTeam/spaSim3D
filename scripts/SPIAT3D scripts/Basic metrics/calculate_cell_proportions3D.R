calculate_cell_proportions3D <- function(spe,
                                         cell_types_chosen = NULL, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), "Cell.Type" = spe[[feature_colname]])
  
  # Check
  if (nrow(df) == 0) stop("No cells found for calculating cell proportions")
  
  # Creates frequency/bar plot of all cell types in the entire image
  cell_proportions <- data.frame(table(data[, feature_colname]))
  names(cell_proportions) <- c("Cell.Type", 'Frequency')
  
  # Only include cell types the user has chosen
  if (!is.null(cell_types_chosen)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_chosen, cell_proportions$Cell.Type)
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in cell_types_chosen are not found in the spe object:\n   ",
                 paste(unknown_cell_types, collapse = ", ")))
    }
    
    # Subset for 
    cell_proportions <- cell_proportions[(cell_proportions$Cell.Type %in% cell_types_chosen), ]
    
    # Check if the user has excluded all cell types
    if (nrow(cell_proportions) == 0) {
      stop("All cells have been excluded")
    }
  }
  
  # Get frequency total for all cells
  cell_type_frequency_total <- sum(cell_proportions$Frequency)
  
  # Get proportions and percentages
  cell_proportions$Proportion <- cell_proportions$Frequency / cell_type_frequency_total
  cell_proportions$Percentage <- cell_proportions$Proportion * 100

  # Order the cell types by proportion (highest cell proportion is first)
  cell_proportions <- cell_proportions[rev(order(cell_proportions$Proportion)), ]
  
  if (plot_image) {
    
    labels <- paste(round(cell_proportions$Percentage, 1), "%", sep = "")
    
    g <- ggplot(cell_proportions, aes(x = factor(Cell.Type, Cell.Type), y = Percentage, fill = Cell.Type)) +
      geom_bar(stat='identity') + 
      theme_bw() +
      labs(title="Cell proportions", x = "Cell type", y = "Percentage") +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "none") +
      geom_text(aes(label = labels), vjust = 0)
    
    methods::show(g)
  }
  
  return(cell_proportions)
}
