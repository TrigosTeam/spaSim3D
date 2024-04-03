calculate_k_cross3D <- function(data, 
                                reference_cell_type,
                                target_cell_type,
                                distance,
                                feature_colname = "Cell.Type", 
                                plot_results = TRUE) {
  
  
  ## Get x, y, z coords for cells of reference cell type and target cell type
  reference_cell_data <- data[data[[feature_colname]] == reference_cell_type, ]
  target_cell_data <- data[data[[feature_colname]] == target_cell_type, ]
  
  ## Get rough dimensions of the window the points are in
  length <- round(max(data$Cell.X.Position) - min(data$Cell.X.Position))
  width  <- round(max(data$Cell.Y.Position) - min(data$Cell.Y.Position))
  height <- round(max(data$Cell.Z.Position) - min(data$Cell.Z.Position))
  
  ## Get distance of each reference cell from the border of the window
  ## Defined as the minimum of the x, y and z distance
  reference_cell_data$Cell.X.Position.Inverse <- length - reference_cell_data$Cell.X.Position
  reference_cell_data$Cell.Y.Position.Inverse <- width  - reference_cell_data$Cell.Y.Position
  reference_cell_data$Cell.Z.Position.Inverse <- height - reference_cell_data$Cell.Z.Position
  reference_cell_data$Border.Distance <- as.numeric(apply(reference_cell_data, 1, min))
  
  ## Only use reference cells which a far enough away from the border
  chosen_rows <- rownames(reference_cell_data[reference_cell_data$Border.Distance > distance, ])
  reference_cell_data <- reference_cell_data[chosen_rows, c("Cell.X.Position",
                                                            "Cell.Y.Position",
                                                            "Cell.Z.Position",
                                                            "Cell.Type")]
  
  ## Get number of reference and target cells
  n_reference_cells <- nrow(reference_cell_data)
  n_target_cells <- nrow(target_cell_data)
  
  ## Combine together
  combined_cell_data <- rbind(reference_cell_data, target_cell_data)
  
  ## Get distances between chosen cell types
  reference_target_distances <- -1 * apcluster::negDistMat(combined_cell_data[, c("Cell.X.Position",
                                                                                  "Cell.Y.Position",
                                                                                  "Cell.Z.Position")])
  
  # Only need distances between reference cells and target cells
  # Ignore ref-ref or target-target distances (i.e. top right part of matrix)
  reference_target_distances <- reference_target_distances[1:n_reference_cells, 
                                                           (n_reference_cells + 1):ncol(reference_target_distances)]
  
  
  # Calculate observed cross-k value for a sequence of distances
  # i.e. the number of ref-target distances less than the chosen distance
  distances <- seq(1, distance, 0.25)
  observed_k <- unlist(lapply(distances, function(x) sum(reference_target_distances < x)))
  
  # Get volume of the window the cells are in
  volume <- length * width * height
  
  # Calculate expected cross-k value for a sequence of distances (using the formula?)
  expected_k <- n_reference_cells * n_target_cells * ((4/3) * pi * distances^3) / volume
  
  result <- data.frame(Distance = distances,
                       Observed = observed_k,
                       Expected = expected_k)
  
  if (plot_results) {
    plot(result$Distance, result$Observed, type = "l", col = "red", 
         xlim = c(0, distance), ylim = c(0, max(result)),
         xlab = "Distance", ylab = "Cross K-function Value")
    lines(result$Distance, result$Expected, type = "l", col = "blue", lty = 2)
    legend(0, max(result), legend = c("Observed K", "Expected K"), col = c("red", "blue"), lty = c(1, 2))
  }
  
  return (result)
}
