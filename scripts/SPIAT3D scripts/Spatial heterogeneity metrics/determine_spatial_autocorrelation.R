determine_spatial_autocorrelation <- function(grid_data,
                                              metric_colname = "Entropy",
                                              weight_method = "IDW") {
  
  
  ## Get number of grid prisms
  n_grid_prisms <- nrow(grid_data)
  
  ## Get splitting number (should be the cube root of n_grid_prisms)
  n_split <- (n_grid_prisms)^(1/3)
  
  ## Find the coordinates of each grid prism
  x <- ((seq(n_grid_prisms) - 1) %% n_split)
  y <- (floor(((seq(n_grid_prisms) - 1) %% (n_split)^2) / n_split))
  z <- (floor((seq(n_grid_prisms) - 1) / (n_split^2)))
  grid_prism_coords <- data.frame(x = x, y = y, z = z)
  
  
  weight_matrix <- -1 * apcluster::negDistMat(grid_prism_coords)
  ## Use the inverse distance between two points as the weight (IDW is 'inverse distance weighting')
  if (weight_method == "IDW") {
    weight_matrix <- 1 / weight_matrix
  }
  ## Use binary method: adjacent points get a weight of 1, otherwise, weight of 0
  ## Adjacent points are within sqrt(3) units apart. e.g. (0, 0, 0) vs (1, 1, 1)
  else if (weight_method == "Binary") {
    weight_matrix <- ifelse(weight_matrix > sqrt(3), 0, 1)  
  }
  
  ## Points along the diagonal are comparing the same point so its weight is zero
  diag(weight_matrix) <- 0
  
  data_mean <- mean(grid_data[!is.na(grid_data[[metric_colname]]), metric_colname])
  
  numerator <- 0
  denominator <- 0
  
  for (i in seq(n_grid_prisms)) {
    
    if (is.na(grid_data[i, metric_colname])) {
      next
    }
    
    for (j in seq(n_grid_prisms)) {
      
      if (is.na(grid_data[j, metric_colname])) {
        next
      }
      
      numerator <- numerator + weight_matrix[i, j] * 
                              (grid_data[i, metric_colname] - data_mean) * 
                              (grid_data[j, metric_colname] - data_mean)
      
    }
    denominator <- denominator + (grid_data[i, metric_colname] - data_mean)^2
  }
  
  
  I <- (n_grid_prisms * numerator) / (sum(weight_matrix) * denominator)
  
  return (I)
  
}
