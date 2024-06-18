summarise_distances_between_cell_types3D <- function(df) {
  
  pair <- distance <- NULL
  
  # summarise the results
  summarised_dists <- df %>% 
    dplyr::group_by(pair) %>%
    dplyr::summarise(mean(distance, na.rm = TRUE), 
                     min(distance, na.rm = TRUE), 
                     max(distance, na.rm = TRUE),
                     stats::median(distance, na.rm = TRUE), 
                     stats::sd(distance, na.rm = TRUE))
  
  summarised_dists <- data.frame(summarised_dists)
  
  colnames(summarised_dists) <- c("pair", 
                                  "mean", 
                                  "min", 
                                  "max", 
                                  "median", 
                                  "std_dev")
  
  for (i in seq(nrow(summarised_dists))) {
    # Get cell_types for each pair
    cell_types <- strsplit(summarised_dists[i,"pair"], "/")[[1]]
    
    summarised_dists[i, "reference"] <- cell_types[1]
    summarised_dists[i, "target"] <- cell_types[2]
  }
  
  return(summarised_dists)
}
