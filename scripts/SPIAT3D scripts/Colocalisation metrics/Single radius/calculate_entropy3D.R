calculate_entropy3D <- function(spe,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname = "Cell.Type",
                                plot_image = TRUE) {
  
  # Check
  if (length(target_cell_types) < 2) stop("Need at least two target cell types")
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighborhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                 reference_cell_type,
                                                                 target_cell_types,
                                                                 radius,
                                                                 feature_colname,
                                                                 FALSE,
                                                                 FALSE)
  
  ## Get total number of target cells for each row (first column is the reference cell id column, so we exclude it)
  cells_in_neighborhood_df$total <- apply(cells_in_neighborhood_df[ , c(-1)], 1, sum)
  
  ## Get entropy for each row
  cells_in_neighborhood_df$entropy <- 0
  
  for (target_cell_type in target_cell_types) {
    
    target_cell_type_proportions <- (cells_in_neighborhood_df[[target_cell_type]] / cells_in_neighborhood_df$total)
    
    ## If an element in target_cell_type_proportion is 0, just add 0.    
    target_cell_entropy <- ifelse(target_cell_type_proportions == 0,
                                  0,
                                  -1 * target_cell_type_proportions * log(target_cell_type_proportions, length(target_cell_types)))
    
    cells_in_neighborhood_df$entropy <- cells_in_neighborhood_df$entropy + target_cell_entropy
    
  }
  
  ## Plot
  if (plot_image) {
    fig <- plot_entropy_violin3D(cells_in_neighborhood_df)
    methods::show(fig)
  }
  
  return(cells_in_neighborhood_df)
}
