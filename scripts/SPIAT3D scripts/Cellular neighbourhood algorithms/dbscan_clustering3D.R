library(dbscan)

dbscan_clustering3D <- function(spe,
                                cell_types_of_interest,
                                radius,
                                minimum_cells_in_radius,
                                feature_colname = "Cell.Type",
                                plot_image = T) {
  
  spe_subset <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  spe_subset_coords <- spatialCoords(spe_subset)
  
  db <- dbscan::dbscan(spe_subset_coords, eps = radius, minPts = minimum_cells_in_radius, borderPoints = F)
  
  
  
  # Convert spe object to data frame
  df <- data.frame(spatialCoords(spe),
                   "Cell.Type" = spe[[feature_colname]],
                   "Cell.ID" = spe[["Cell.ID"]])

  df_cell_types_of_interest <- df[df$Cell.Type %in% cell_types_of_interest, ]
  df_other_cell_types <- df[!(df$Cell.Type %in% cell_types_of_interest), ]
  df_cell_types_of_interest$dbscan_cluster <- db$cluster
  df_other_cell_types$dbscan_cluster <- 0
  
  ## Convert data frame to spe object
  df <- rbind(df_cell_types_of_interest, df_other_cell_types)
  
  spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = spe@metadata)
  
  ## Plot
  if (plot_image) {
    df$dbscan_cluster <- ifelse(df$dbscan_cluster == 0, "non_cluster", paste("cluster_", df$dbscan_cluster, sep = ""))
    
    fig <- plot_ly(df,
                   type = "scatter3d",
                   mode = 'markers',
                   x = ~Cell.X.Position,
                   y = ~Cell.Y.Position,
                   z = ~Cell.Z.Position,
                   color = ~dbscan_cluster,
                   colors = rainbow(length(unique(df$dbscan_cluster))),
                   marker = list(size = 2)) %>% 
      layout(scene = list(xaxis = list(title = 'x'),
                          yaxis = list(title = 'y'),
                          zaxis = list(title = 'z')))
    
    methods::show(fig)
  }
  
  return(spe)
}
