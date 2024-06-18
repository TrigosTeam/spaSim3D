simulate_cylinder_cluster <- function(bg_spe, cluster_properties) {
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(bg_spe), 
                   "Cell.Type" = bg_spe[["Cell.Type"]],
                   "Cell.ID" = bg_spe[["Cell.ID"]])
  
  # Get cylinder properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  radius <- cluster_properties$radius
  start_loc <- cluster_properties$start_loc
  end_loc <- cluster_properties$end_loc
  
  # Get number of cells
  n_cells <- nrow(df)
  
  # Get number of unique cell types
  n_cluster_cell_types <- length(cluster_cell_types)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  i <- 1
  
  while (i <= n_cells) {
    # Get x, y, z coordinate of current cell
    x <- df[i, "Cell.X.Position"]
    y <- df[i, "Cell.Y.Position"]
    z <- df[i, "Cell.Z.Position"]
    
    # Ignore points outside of these planes
    if (sum(v1 *  c(x, y, z)) < d1 || sum(v1 * c(x, y, z)) > d2) {
      i <- i + 1
      next
    }
    
    # Get vector between from point to start_loc
    v2 <- c(x, y, z) - start_loc
    
    # Get perpendicular distance squared between point and line
    D <- ((v1[2]*v2[3] - v1[3]*v2[2])^2 + 
          (v1[1]*v2[3] - v1[3]*v2[1])^2 + 
          (v1[1]*v2[2] - v1[2]*v2[1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2)
    
    # Dumb case where the start and end loc is the same
    if (is.nan(D)) D <- Inf
    
    # Get maximum distance squared
    R <- radius^2
    
    if (D < R) { 
      # Random number will determine the cluster_cell_type of the cell
      random <- stats::runif(1)
      
      # Start with the first cell
      n <- 1
      current_proportion <- 0
      
      while (n <= n_cluster_cell_types){
        current_proportion <- current_proportion + cluster_cell_proportions[n]
        if (random <= current_proportion) {
          df[i, "Cell.Type"] <- cluster_cell_types[n]
          break
        }
        n <- n + 1
      }
    }
    
    if (df[i, "Cell.Type"] == "Void") { 
      df <- df[-c(i), ]
      n_cells <- n_cells - 1
      
    } else {
      i <- i + 1
    }
  }
  
  # Update current meta data
  metadata <- bg_spe@metadata
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  metadata[[paste("cluster", length(metadata), sep="_")]] <- cluster_properties
  
  # Convert data frame to spe object
  cluster_spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = metadata)
  
  return(cluster_spe)
}
