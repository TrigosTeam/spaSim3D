determine_prevalence3D <- function(grid_data,
                                   metric_colname,
                                   threshold,
                                   above = TRUE) {
  
  ## Exclude rows with NA values
  grid_data <- grid_data[!is.na(grid_data[[metric_colname]]), ]
  
  if (above) {
    p <- sum(grid_data[[metric_colname]] >= threshold) / nrow(grid_data) * 100
  }
  else {
    p <- sum(grid_data[[metric_colname]] < threshold) / nrow(grid_data) * 100    
  }
  
  return(p)
}
