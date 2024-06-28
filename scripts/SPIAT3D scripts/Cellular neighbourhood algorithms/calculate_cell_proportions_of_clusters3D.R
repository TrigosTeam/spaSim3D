calculate_cell_proportions_of_clusters3D <- function(spe, cluster_colname, feature_colname = "Cell.Type", plot_image = T) {
  
  ## Get cluster numbers (ignoring 0)
  cluster_numbers <- spe[[cluster_colname]][spe[[cluster_colname]] != 0]
  
  ## Get number of clusters
  n_clusters <- length(unique(cluster_numbers))
  
  ## Get different cell types found in the clusters (alphabetical for consistency)
  cell_types <- unique(spe[[feature_colname]][spe[[cluster_colname]] != 0])
  cell_types <- cell_types[order(cell_types)]
  
  ## For each cluster, determine the size and cell proportion of each cluster
  result <- data.frame(matrix(nrow = n_clusters, ncol = 2 + length(cell_types)))
  colnames(result) <- c("cluster_number", "n_cells", cell_types)
  
  for (i in seq(n_clusters)) {
    cells_in_clusters <- spe[[feature_colname]][spe[[cluster_colname]] == i]
    result[i, "n_cells"] <- length(cells_in_clusters)
    
    for (cell_type in cell_types) {
      result[i, cell_type] <- sum(cells_in_clusters == cell_type) / result[i, "n_cells"]
    }
  }
  
  result <- result[order(result$n_cells), ]
  rownames(result) <- seq(n_clusters)
  result$cluster_number <- as.character(seq(n_clusters))
  
  ## Plot
  if (plot_image) {
    plot_result <- reshape2::melt(result, id.vars = c("cluster_number", "n_cells"))
    fig <- ggplot(plot_result, aes(cluster_number, value, fill = variable)) +
      geom_bar(stat = "identity") +
      labs(title = "Cell proportions of each cluster", x = "", y = "Cell proportion") +
      scale_x_discrete(labels = paste("cluster_", result$cluster_number, ", n = ", result$n_cells, sep = "")) +
      guides(fill = guide_legend(title="Cell type")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    
    methods::show(fig)
  }
  
  return(result)
}

