calculate_border_of_clusters3D <- function(spe, 
                                           radius,
                                           cluster_colname, 
                                           feature_colname = "Cell.Type", 
                                           plot_image = T) {
  
  ## Get spatial coords of spe
  spe_coords <- data.frame(spatialCoords(spe))
  
  ## Get coords of non-cluster cells
  non_cluster_coords <- spe_coords[spe[[cluster_colname]] == 0, ]
  
  # New column for spe object: 'cluster_border'. Default is 'outside'
  spe$cluster_border <- "outside"
  
  # Label cells part of a cluster (e.g. 'cluster1')
  spe$cluster_border[spe[[cluster_colname]] != 0] <- paste("inside_C", spe[[cluster_colname]][spe[[cluster_colname]] != 0], sep = "")
  
  ## Iterate for each cluster
  n_clusters <- max(spe[[cluster_colname]])
  
  for (i in seq_len(n_clusters)) {
    
    ## Subset for cells in the current cluster of interest
    cluster_coords <- spe_coords[spe[[cluster_colname]] == i, ]
    
    # For each cell in the current cluster, check how many other cells in the cluster are in its radius
    cluster_to_cluster_interactions <- dbscan::frNN(cluster_coords, radius)
    
    # Determine the median minimum number of cluster cells found in the radius of cluster cell. Use this as the threshold for non-cluster cells.
    non_cluster_threshold <- quantile(unlist(lapply(cluster_to_cluster_interactions$dist, length)), 0.5)
    
    # For each non-cluster cell, check how many cluster cells are in its radius.
    non_cluster_to_cluster_interactions <- dbscan::frNN(cluster_coords, radius, non_cluster_coords)
    
    # If number of cluster cells found in the radius of non-cluster cells is greater than threshold, non-cluster cell has probably infiltrated cluster too
    n_cluster_cells_in_non_cluster_cell_radius <- unlist(lapply(non_cluster_to_cluster_interactions$id, length))
    
    spe$cluster_border[as.numeric(names(non_cluster_to_cluster_interactions$id)[n_cluster_cells_in_non_cluster_cell_radius > non_cluster_threshold])] <- paste("infiltrated_C", i, sep = "")
    
    # If number of cluster cells found in the radius of non-cluster cells is less than threshold, but greater than 0, non-cluster cell is probably on the border
    spe$cluster_border[as.numeric(names(non_cluster_to_cluster_interactions$id)[n_cluster_cells_in_non_cluster_cell_radius > 0 & n_cluster_cells_in_non_cluster_cell_radius < non_cluster_threshold])] <- paste("border_C", i, sep = "")
  }
  
  ## Plot
  if (plot_image) {
    fig <- plot_cells3D(spe, feature_colname = "cluster_border")
    methods::show(fig)
  }
      
  return(spe)
}