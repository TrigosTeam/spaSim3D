library(alphashape3d)

determine_alpha_hull3D <- function(spe, 
                                   cell_types_of_interest, 
                                   alpha, 
                                   feature_colname = "Cell.Type", 
                                   plot_image = T) {
  
  ## Check cell types of interst are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, spe$Cell.Type)
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  ## Subset for the chosen cell_types_of_interest
  spe_subset <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  spe_coords <- spatialCoords(spe_subset)
  
  ## Get the alpha hull
  alpha_hull <- ashape3d(as.matrix(spe_coords), alpha = 15)
  
  ## Get the information of the vertices and faces of the alpha hull (what 3 vertices make up each face triangle?)
  vertices <- alpha_hull$x
  faces <- alpha_hull$triang[alpha_hull$triang[, 9] != 0, c("tr1", "tr2", "tr3")]
  
  
  ## Plot
  if (plot_image) {
    
    ## Convert spe object to data frame
    df <- data.frame(spatialCoords(spe), "Cell.Type" = spe[[feature_colname]])
    
    fig <- plot_ly() %>%
      add_trace(
        data = df,
        type = "scatter3d",
        mode = 'markers',
        x = ~Cell.X.Position,
        y = ~Cell.Y.Position,
        z = ~Cell.Z.Position,
        marker = list(size = 2),
        color = as.formula(paste0('~', feature_colname)),
        colors = rainbow(length(unique(df[[feature_colname]])))
      ) %>%
      add_trace(
        type = 'mesh3d',
        x = vertices[, 1], 
        y = vertices[, 2], 
        z = vertices[, 3],
        i = faces[, 1] - 1, 
        j = faces[, 2] - 1, 
        k = faces[, 3] - 1,
        opacity = 0.05,
        facecolor = rep("red", nrow(faces))
      )
    
    methods::show(fig)
  }
  return(1)
}



determine_alpha_hull3D(spe1, c("Tumour"))
