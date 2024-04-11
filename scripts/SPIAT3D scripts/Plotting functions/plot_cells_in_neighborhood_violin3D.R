## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_cells_in_neighborhood_violin3D <- function(cells_in_neighborhood, scales = "free_x") {
  
  df <- data.frame(matrix(nrow = 0, ncol = 3))
  colnames(df) <- c("Target", "Count", "Reference")
  
  for (i in seq(length(cells_in_neighborhood))) {
    
    # Get reference cell type for current index
    reference_cell_type <- names(cells_in_neighborhood)[i]
    
    # Get data for current index
    cells_in_neighborhood_df <- cells_in_neighborhood[[i]]
    
    # Get columns which contain cell count data (4th column onwards)
    cells_in_neighborhood_df <- cells_in_neighborhood_df[ , 4:ncol(cells_in_neighborhood_df)]
    
    # Melt
    cells_in_neighborhood_df <- reshape2::melt(cells_in_neighborhood_df, id.vars = 0)
    colnames(cells_in_neighborhood_df) <- c("Target", "Count")
    
    # Add reference cell type column
    cells_in_neighborhood_df$Reference <- reference_cell_type
    
    # Add result to main df
    df <- rbind(df, cells_in_neighborhood_df)
    
  }
  
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Reference <- Count <- Target <- NULL
  
  ggplot(df, aes(x = Reference, y = Count, fill = Target)) + geom_violin() +
    facet_wrap(~Reference, scales=scales) +
    theme_bw() +
    theme(axis.text.x=element_blank())
  
}
