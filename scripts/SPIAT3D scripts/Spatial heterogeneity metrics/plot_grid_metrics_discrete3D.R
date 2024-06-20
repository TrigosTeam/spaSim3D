plot_grid_metrics_discrete3D <- function(spe, grid_metrics, metric_colname) {
  
  spe_coords <- data.frame(spatialCoords(spe))
  
  ## Get dimensions of the window
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  
  
  ## Get distance of row, col and lay
  n_grid_prisms <- nrow(grid_metrics)
  n_split <- n_grid_prisms^(1/3)
  d_row <- length / n_split
  d_col <- width / n_split
  d_lay <- height / n_split
  
  ## Add x, y and z coords of each grid prism to data
  grid_metrics$x <- ((seq(n_grid_prisms) - 1) %% n_split + 0.5) * d_row
  grid_metrics$y <- (floor(((seq(n_grid_prisms) - 1) %% (n_split)^2) / n_split) + 0.5) * d_col
  grid_metrics$z <- (floor((seq(n_grid_prisms) - 1) / (n_split^2)) + 0.5) * d_lay
  
  
  ## Define low, medium and high categories
  # Low: between 0 and 1/3
  # Medium: between 1/3 and 2/3
  # High: between 2/3 and 1
  
  grid_metrics$rank <- ifelse(is.na(grid_metrics[[metric_colname]]), "na",
                              ifelse(grid_metrics[[metric_colname]] < 1/3, "low",
                                     ifelse(grid_metrics[[metric_colname]] < 2/3, "medium", "high")))
  grid_metrics$rank <- factor(grid_metrics$rank, c("low", "medium", "high", "na"))
  
  fig <- plot_ly(grid_metrics,
                 type = "scatter3d",
                 mode = 'markers',
                 x = ~x,
                 y = ~y,
                 z = ~z,
                 color = ~rank,
                 colors = c("#AEB6E5", "#BC6EB9", "#A93154", "gray"),
                 symbol = 1,
                 symbols = "square",
                 marker = list(size = 4))
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                     yaxis = list(title = 'y'),
                                     zaxis = list(title = 'z')))
  return(fig)
}
