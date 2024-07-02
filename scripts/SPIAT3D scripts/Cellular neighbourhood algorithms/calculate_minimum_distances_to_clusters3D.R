calculate_minimum_distances_to_clusters3D <- function(spe, cell_types_inside_cluster, cell_types_outside_cluster, cluster_colname, feature_colname = "Cell.Type", plot_image = T) {
  
  ## Get cluster numbers (ignoring 0)
  cluster_numbers <- spe[[cluster_colname]][spe[[cluster_colname]] != 0]
  
  ## Get number of clusters
  n_clusters <- length(unique(cluster_numbers))

  ## For each cell type outside clusters, get their set of coords. These exclude cell types in clusters
  spe_coords <- spatialCoords(spe)
  cluster_rows <- rep(FALSE, nrow(spe_coords))
  for (i in seq(n_clusters)) {
    cluster_rows <- cluster_rows | (spe[[cluster_colname]] == i)
  }

  spe_outside_cluster <- spe[ , !cluster_rows]
  cell_types_outside_cluster_coords <- list()
  for (cell_type in cell_types_outside_cluster) {
    cell_types_outside_cluster_coords[[cell_type]] <- spatialCoords(spe_outside_cluster)[spe_outside_cluster[[feature_colname]] == cell_type, ]
  }
  
  ## For each cluster, determine the minimum distance of each cell_type_of_interest  
  result <- vector()
  
  for (i in seq(n_clusters)) {
    cluster_coords <- spe_coords[spe[[cluster_colname]] == i & spe[[feature_colname]] %in% cell_types_inside_cluster, ]
    
    for (cell_type in cell_types_outside_cluster) {
      curr_cell_type_coords <- cell_types_outside_cluster_coords[[cell_type]]
      
      all_closest <- RANN::nn2(data = cluster_coords, 
                               query = curr_cell_type_coords, 
                               k = 1)  
      
      local_dist_mins <- data.frame(
        cluster_number = i,
        cell_type_of_interest = cell_type,
        distance = all_closest$nn.dists
      )
      ## Remove any 0 distance rows
      local_dist_mins <- local_dist_mins[local_dist_mins$distance != 0, ]
      result <- rbind(result, local_dist_mins)
    }
    

    ## Plot
    if (plot_image) {

      cluster_number_labs <- paste("cluster_", seq(n_clusters), sep = "")
      names(cluster_number_labs) <- seq(n_clusters)
      
      fig <- ggplot(result, aes(x = cell_type_of_interest, y = distance, fill = cell_type_of_interest)) + 
        geom_violin() +
        facet_grid(cluster_number~., scales="free_x", labeller = labeller(cluster_number = cluster_number_labs)) +
        theme_bw() +
        theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none") +
        labs(title="Minimum cell distances to clusters", x = "Cell type", y = "Distance") +
        stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
      
      methods::show(fig)
    }
    
  }
  return(result)
}
