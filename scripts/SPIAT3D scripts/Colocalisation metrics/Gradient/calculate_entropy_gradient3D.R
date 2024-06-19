calculate_entropy_gradient3D <- function(spe,
                                         reference_cell_type,
                                         target_cell_types,
                                         radii,
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (radius in seq(radii)) {
    cells_in_neighborhood_data <- calculate_cells_in_neighborhood3D(spe,
                                                                    reference_cell_type,
                                                                    target_cell_types,
                                                                    radius,
                                                                    feature_colname,
                                                                    FALSE,
                                                                    FALSE)
    
    cells_in_neighborhood_data$ref_cell_id <- NULL
    result[radius, ] <- apply(cells_in_neighborhood_data, 2, sum)
  }
  
  ## Get total number of target cells for each row
  result$total <- apply(result, 1, sum)
  
  ## Set intial entropy to 0
  result$entropy <- 0
  
  for (target_cell_type in target_cell_types) {
    
    target_cell_type_proportions <- (result[[target_cell_type]] / result$total)
    
    ## If an element in target_cell_type_proportion is 0, just add 0.    
    target_cell_entropy <- ifelse(target_cell_type_proportions == 0,
                                  0,
                                  -1 * target_cell_type_proportions * log(target_cell_type_proportions, length(target_cell_types)))
    
    result$entropy <- result$entropy + target_cell_entropy
    
  }
  
  # Add a radius column to the result
  result$radius <- seq(radii)
  
  if (plot_image) {
    plot(result$radius, result$entropy, type = "l", col = "red", 
         xlim = c(0, radius), ylim = c(0, max(result$entropy)),
         xlab = "Radius", ylab = "Entropy")
    # lines(result$radius, result$expected_cross_K, type = "l", col = "blue", lty = 2)
    # legend(0, max(result), legend = c("Observed cross K", "Expected cross K"), col = c("red", "blue"), lty = c(1, 2))
  }
  
  return(result)
}
