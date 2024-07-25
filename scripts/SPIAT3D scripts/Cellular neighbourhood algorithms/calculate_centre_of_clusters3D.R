### Assume that clusters have uniform density and that the centre of each cluster is defined by its centre of mass
### Centre of mass can be estimated by taking the average of the x, y, and z coordinates of cells in the cluster

calculate_center_of_clusters3D <- function(spe, cluster_colname) {
  
  # Get number of clusters
  n_clusters <- max(spe[[cluster_colname]])
  
  ## For each cluster, determine the number of cells in each cluster of each cluster
  result <- data.frame(matrix(nrow = n_clusters, ncol = 4))
  colnames(result) <- c("cluster_number", "Centre.X.Position", "Centre.Y.Position", "Centre.Z.Position")
  
  result$cluster_number <- as.character(seq(n_clusters))
  for (i in seq(n_clusters)) {
    spe_cluster <- spe[ , spe[[cluster_colname]] == i]
    result[i, c("Centre.X.Position", "Centre.Y.Position", "Centre.Z.Position")] <- 
      apply(spatialCoords(spe_cluster), 2, mean)
  }
  
  return(result)
}

