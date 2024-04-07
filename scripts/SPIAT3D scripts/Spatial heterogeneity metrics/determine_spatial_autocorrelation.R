determine_spatial_autocorrelation <- function(data,
                                              entropy_grid_data, 
                                              n_split) {
  
  ## Get dimensions of the window
  length <- round(max(data$Cell.X.Position) - min(data$Cell.X.Position))
  width  <- round(max(data$Cell.Y.Position) - min(data$Cell.Y.Position))
  height <- round(max(data$Cell.Z.Position) - min(data$Cell.Z.Position))
  
  ## Get distance of row, col and lay
  d_row <- length / n_split
  d_col <- width / n_split
  d_lay <- height / n_split
  
  ## Get number of grid prisms
  n_grid_prisms <- nrow(entropy_grid_data)
  
  ## Find the coordinates of each grid prism
  x <- ((seq(n_grid_prisms) - 1) %% n_split) * d_row
  y <- (floor(((seq(n_grid_prisms) - 1) %% (n_split)^2) / n_split)) * d_col
  z <- (floor((seq(n_grid_prisms) - 1) / (n_split^2))) * d_lay
  grid_prism_coords <- data.frame(x = x, y = y, z = z)
  
  ## Use the inverse distance between two points as the weight 
  weight_matrix <- -1 * apcluster::negDistMat(grid_prism_coords)
  weight_matrix <- 1 / weight_matrix
  ## Points along the diagonal are comparing the same point so its weight is zero
  diag(weight_matrix) <- 0
  
  entropy_mean <- mean(entropy_grid_data$Entropy)
  
  numerator <- 0
  denominator <- 0
  
  for (i in seq(n_grid_prisms)) {
    
    for (j in seq(n_grid_prisms)) {
      
      numerator <- numerator + weight_matrix[i, j] * 
                              (entropy_grid_data[i, "Entropy"] - entropy_mean) * 
                              (entropy_grid_data[j, "Entropy"] - entropy_mean)
      
    }
    denominator <- denominator + (entropy_grid_data[i, "Entropy"] - entropy_mean)^2
  }
  
  
  I <- (n_grid_prisms * numerator) / (sum(weight_matrix) * denominator)
  
  return (I)
  
}
