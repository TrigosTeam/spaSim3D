library(alphashape3d)

determine_alpha_hull3D <- function(spe, 
                                   cell_types_of_interest, 
                                   alpha = NULL, 
                                   minimum_cells_in_alpha_hull,
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
  spe_subset_coords <- spatialCoords(spe_subset)
  
  ## Get alpha value if not specified by user
  if (is.null(alpha)) {
    spe_coords <- spatialCoords(spe)
    window_volume <- 
      (max(spe_coords[, "Cell.X.Position"]) - min(spe_coords[, "Cell.X.Position"])) * 
      (max(spe_coords[, "Cell.Y.Position"]) - min(spe_coords[, "Cell.Y.Position"])) * 
      (max(spe_coords[, "Cell.Z.Position"]) - min(spe_coords[, "Cell.Z.Position"]))
    n_cells <- nrow(spe_coords)
    
    ### Estimated alpha is 10% of the ratio between the window volume and the number of cells
    alpha <- 0.1 * (window_volume / n_cells) 
    print(paste("No alpha inputted. Choosing alpha to be", round(alpha, 2)))
  }
  
  ## Get the alpha hull
  alpha_hull <- ashape3d(as.matrix(spe_subset_coords), alpha = alpha)
  
  ## Determine which alpha hull each cell_type_of_interest belongs to
  alpha_hull_numbers <- components_ashape3d(alpha_hull)
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), 
                   "Cell.Type" = spe[[feature_colname]],
                   "Cell.ID" = spe[["Cell.ID"]])
  
  df_cell_types_of_interest <- df[df$Cell.Type %in% cell_types_of_interest, ]
  df_other_cell_types <- df[!(df$Cell.Type %in% cell_types_of_interest), ]
  df_cell_types_of_interest$alpha_hull_number <- alpha_hull_numbers
  df_other_cell_types$alpha_hull_number <- -1
  
  ## Ignore cell_types_of_interest which belong to an alpha hull cluster with less than minimum_cells_in_alpha_hull
  alpha_hull_numbers_table <- table(alpha_hull_numbers)
  maximium_alpha_hull_number <- Position(function(x) x < minimum_cells_in_alpha_hull, alpha_hull_numbers_table)
  maximium_alpha_hull_number <- as.numeric(names(alpha_hull_numbers_table[maximium_alpha_hull_number]))
  
  if (!is.na(maximium_alpha_hull_number) && maximium_alpha_hull_number != -1) {
    spe_subset_coords <- spe_subset_coords[alpha_hull_numbers >= 1 & alpha_hull_numbers < maximium_alpha_hull_number, ]
    
    df_cell_types_of_interest$alpha_hull_number <- ifelse(alpha_hull_numbers >= 1 & alpha_hull_numbers < maximium_alpha_hull_number, 
                                                           alpha_hull_numbers, -1)
  
    ## Get the alpha hull again...
    alpha_hull <- ashape3d(as.matrix(spe_subset_coords), alpha = alpha)
  }
  
  ## Convert data frame to spe object
  df <- rbind(df_cell_types_of_interest, df_other_cell_types)
  
  spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = spe@metadata)
  
  ## Get the information of the vertices and faces of the alpha hull (what 3 vertices make up each face triangle?)
  vertices <- alpha_hull$x
  faces <- alpha_hull$triang[alpha_hull$triang[, 9] != 0, c("tr1", "tr2", "tr3")]
  spe@metadata$alpha_hull <- list(vertices = vertices, faces = faces)
  
  ## Plot
  if (plot_image) {
    fig <- plot_alpha_hull3D(spe, feature_colname = feature_colname)
    methods::show(fig)
  }
  
  return(spe)
}
