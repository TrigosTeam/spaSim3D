summarise_distances_between_cell_types3D <- function(df) {
  
  Pair <- Distance <- NULL
  
  # summarise the results
  summarised_dists <- df %>% 
    dplyr::group_by(Pair) %>%
    dplyr::summarise(mean(Distance, na.rm = TRUE), 
                     min(Distance, na.rm = TRUE), 
                     max(Distance, na.rm = TRUE),
                     stats::median(Distance, na.rm = TRUE), 
                     stats::sd(Distance, na.rm = TRUE))
  
  summarised_dists <- data.frame(summarised_dists)
  
  colnames(summarised_dists) <- c("Pair", 
                                  "Mean", 
                                  "Min", 
                                  "Max", 
                                  "Median", 
                                  "Std.Dev")
  
  for (i in seq(nrow(summarised_dists))) {
    # Get cell_types for each pair
    cell_types <- strsplit(summarised_dists[i,"Pair"], "/")[[1]]
    
    summarised_dists[i, "Reference"] <- cell_types[1]
    summarised_dists[i, "Target"] <- cell_types[2]
  }
  
  return(summarised_dists)
}
