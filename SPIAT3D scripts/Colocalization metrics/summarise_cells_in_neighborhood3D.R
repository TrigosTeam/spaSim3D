summarise_cells_in_neighborhood3D <- function(cells_in_neighborhood_data) {
 
  summarised_results <- list()
  
  for (i in seq(length(cells_in_neighborhood_data))) {
    df <- cells_in_neighborhood_data[[i]]
    
    ## Target cell types will be the fourth column onwards
    target_cell_types <- colnames(df)[4:ncol(df)]
    
    ## Set up data frame for summarised_results list
    df_results <- data.frame(row.names = c("Mean", "Min", "Max", "Median", "St.Dev"))
    
    for (target_cell_type in target_cell_types) {
      
      ## Get statistical measures for each target cell type
      target_cell_type_values <- df[[target_cell_type]]
      df_results[[target_cell_type]] <- c(mean(target_cell_type_values),
                                          min(target_cell_type_values),
                                          max(target_cell_type_values),
                                          median(target_cell_type_values),
                                          sd(target_cell_type_values))
      
      
    }
    
    ## Add data frame result to summarised_results
    summarised_results[[names(cells_in_neighborhood_data)[i]]] <- df_results
  }
  
  return (summarised_results)
}



