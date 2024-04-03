summarise_cells_in_neighborhood <- function(cells_in_neighborhood_df) {
 
  summarised_results <- data.frame(matrix(nrow = 0, ncol = 8))
  
  for (pair in unique(cells_in_neighborhood_df$Pair)) {
    reference_cell_type <- strsplit(pair, "/")[[1]][1]
    target_cell_type <- strsplit(pair, "/")[[1]][2]
    
    pair_data <- cells_in_neighborhood_df[cells_in_neighborhood_df$Pair == pair, "nTarget.Cell"]
    
    summarised_results <- rbind(summarised_results,
                                c(reference_cell_type,
                                  target_cell_type,
                                  pair,
                                  mean(pair_data),
                                  min(pair_data),
                                  max(pair_data),
                                  median(pair_data),
                                  sd(pair_data)))
    
  }
   
  colnames(summarised_results) <- c("Reference.Cell",
                                    "Target.Cell",
                                    "Pair", 
                                    "Mean", 
                                    "Min", 
                                    "Max", 
                                    "Median", 
                                    "Std.Dev")
  
  
  return (summarised_results)
}



