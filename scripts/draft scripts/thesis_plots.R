bg_metadata <- spe_metadata_background_template("ordered")
bg_metadata$background$n_cells <- 20
bg_metadata$background$length <- 10
bg_metadata$background$width <- 10
bg_metadata$background$height <- 8
bg_metadata$background$jitter_proportion <- 0
bg_metadata$background$cell_types <- "Others"
bg_metadata$background$cell_proportions <- 1

bg_spe <- simulate_spe_metadata3D(bg_metadata)

cells_z_coords <- spatialCoords(bg_spe)[ , "Cell.Z.Position"]
layers_z_coords <- unique(cells_z_coords)

bg_spe$Cell.Type <- ifelse(cells_z_coords == layers_z_coords[1], "A",
                           ifelse(cells_z_coords == layers_z_coords[2], "B", "C"))


plot_cells3D_modified <- function(spe,
                                  plot_cell_types = NULL,
                                  plot_colours = NULL,
                                  feature_colname = "Cell.Type") {
  
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), "Cell.Type" = spe[[feature_colname]])
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) {
    plot_cell_types <- unique(df[["Cell.Type"]])
  }
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(plot_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following plot_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  ## If no colours inputted, use rainbow palette
  if (is.null(plot_colours)) {
    plot_colours <- rainbow(length(plot_cell_types))
  }
  
  ## User inputs mismatching cell types and colours
  if (length(plot_cell_types) != length(plot_colours)) {
    stop("Length of plot_cell_types is not equal to length of plot_colours")
  }
  
  ## Factor for feature column
  df[, "Cell.Type"] <- factor(df[, "Cell.Type"],
                              levels = plot_cell_types)
  
  ## Plot
  fig <- plot_ly(df,
                 type = "scatter3d",
                 mode = 'markers',
                 x = ~Cell.X.Position,
                 y = ~Cell.Y.Position,
                 z = ~Cell.Z.Position,
                 color = ~Cell.Type,
                 colors = plot_colours,
                 marker = list(size = 10))
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x', showticklabels=FALSE),
                                     yaxis = list(title = 'y', showticklabels=FALSE),
                                     zaxis = list(title = 'z', showticklabels=FALSE)))
  
  return (fig)
}

plot_cells3D_modified(bg_spe, c("A", "B", "C"), c("lightgreen", "orange", "tomato"))

