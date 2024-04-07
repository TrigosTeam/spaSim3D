average_nearest_neighbor_index3D <- function(data,
                                             cell_types_of_interest,
                                             feature_colname = "Cell.Type",
                                             n_simulations = 100) {
  
  
  
  ## Subset for cell types of interest
  data <- data[which(data[, feature_colname] %in% cell_types_of_interest), ]
  
  ## Assume all cell types of interest are the same
  ## Aiming to find the nearest neighbor for each cell, disregarding cell type
  data[feature_colname] <- "temp"
  
  ## Calculate nearest neighbor and minimum distance for each cell
  nearest_neighbor_data <- calculate_minimum_distances_between_cell_types3D(data,
                                                                            "temp",
                                                                            feature_colname)
  
  ## Calculate observed average distance
  observed_average_distance <- mean(nearest_neighbor_data$Distance)
  
  ## Calculate expected average distance
  ## Use Monte Carlo simulations: Simulate n_simulations samples using poisson distribution and get average
  expected_average_distances <- c()
  
  ## Get number of cells and dimensions of the window
  n_cells <- nrow(data)
  length <- round(max(data$Cell.X.Position) - min(data$Cell.X.Position))
  width  <- round(max(data$Cell.Y.Position) - min(data$Cell.Y.Position))
  height <- round(max(data$Cell.Z.Position) - min(data$Cell.Z.Position))
  
  # Need to think about the min_d parameter
  for (i in 1:n_simulations) {
    simulation <- simulate_background_cells3D(n_cells = nrow(data),
                                              length = length,
                                              width  = width,
                                              height = height,
                                              method = "tumour",
                                              min_d = 2,
                                              oversampling_rate = 1,
                                              jitter_prop = 0,
                                              cell_type = "Others",
                                              plot_image = F)
    
    ## Adding Cell.ID column
    simulation$Cell.ID <- (paste("Cell_", seq(nrow(simulation)), sep="")) 
    
    simulation_nearest_neighbor_data <- calculate_minimum_distances_between_cell_types3D(simulation,
                                                                                         "Others",
                                                                                         "Cell.Type")
    
    expected_average_distances <- c(expected_average_distances, 
                                    mean(simulation_nearest_neighbor_data$Distance))
  }
  
  ## Determine position of observed_average_distance in expected_average_distances vector
  expected_average_distances <- expected_average_distances[order(expected_average_distances)]
  position <- which(expected_average_distances >= observed_average_distance)[1]
  
  ## Calculate p_value
  p_value <- ifelse(position < n_simulations/2, 
                    position/n_simulations, (n_simulations-position)/n_simulations)
  
  
  ## Calculate average nearest neighbor index for each expected_average_distance simulated
  ANNIs <- observed_average_distance / expected_average_distances
  
  ## Get confidence interval for average nearest neighbor index
  ANNI_mean <- mean(ANNIs)

  ANNI_CI_95 <- quantile(ANNIs, c(0.025, 0.975))
  
  return (list(simulated_distances = expected_average_distances,
               observed_distance = observed_average_distance,
               p_value = p_value, 
               ANNI_mean = ANNI_mean,
               ANNI_95confidence_interval = ANNI_CI_95))
}
 
