calculate_entropy3D <- function(spe,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname = "Cell.Type",
                                plot_image = TRUE) {
  
  # Check if radius is numeric
  if (!is.numeric(radius)) {
    stop(paste(radius, " is not of type 'numeric'"))
  }
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighborhood_data <- calculate_cells_in_neighborhood3D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  FALSE,
                                                                  FALSE)
  
  ## Get total number of target cells for each row
  cells_in_neighborhood_data$total <- apply(cells_in_neighborhood_data[ , c(-1)], 1, sum)
  
  ## Get entropy for each row
  cells_in_neighborhood_data$entropy <- 0
  
  for (target_cell_type in target_cell_types) {
    
    target_cell_type_proportions <- (cells_in_neighborhood_data[[target_cell_type]] / cells_in_neighborhood_data$total)

    ## If an element in target_cell_type_proportion is 0, just add 0.    
    target_cell_entropy <- ifelse(target_cell_type_proportions == 0,
                                  0,
                                  -1 * target_cell_type_proportions * log(target_cell_type_proportions, length(target_cell_types)))
    
    cells_in_neighborhood_data$entropy <- cells_in_neighborhood_data$entropy + target_cell_entropy
    
  }
  
  ## Case when row has 0 target cells
  cells_in_neighborhood_data[cells_in_neighborhood_data$Total == 0, "entropy"] <- 0
  
  if (plot_image) {
    fig <- plot_entropy_violin3D(cells_in_neighborhood_data)
    methods::show(fig)
  }
  
  return(cells_in_neighborhood_data)
}
