calculate_cross_K3D <- function(spe, 
                                reference_cell_type, 
                                target_cell_type, 
                                radius, 
                                feature_colname = "Cell.Type") {
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), 
                   "Cell.Type" = spe[[feature_colname]], 
                   "Cell.ID" = spe[["Cell.ID"]])
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% df$Cell.Type)) {
    stop(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
  }
  
  ## For target_cell_type, check it is found in the spe object
  if (!(target_cell_type %in% df$Cell.Type)) {
    stop(paste("The target_cell_type", target_cell_type,"is not found in the spe object"))
  }
  
  
  cells_in_neighbourhood_data <- calculate_cells_in_neighborhood3D(spe,
                                                                   reference_cell_type,
                                                                   target_cell_type,
                                                                   radius,
                                                                   feature_colname,
                                                                   show_summary = FALSE,
                                                                   plot_image = FALSE)
  
  n_ref_tar_interactions <- sum(cells_in_neighbourhood_data[[target_cell_type]])
  n_ref_cells <- sum(df$Cell.Type == reference_cell_type)
  n_tar_cells <- sum(df$Cell.Type == target_cell_type)
  
  ## Get rough dimensions of the window the points are in
  length <- round(max(df$Cell.X.Position) - min(df$Cell.X.Position))
  width  <- round(max(df$Cell.Y.Position) - min(df$Cell.Y.Position))
  height <- round(max(df$Cell.Z.Position) - min(df$Cell.Z.Position))
  ## Get volume of the window the cells are in
  volume <- length * width * height
  
  
  ## Get observed cross K-function
  observed_cross_K <- (volume * n_ref_tar_interactions) / (n_ref_cells * n_tar_cells)
  
  ## Get expected cross K-function
  expected_cross_K <- (4/3) * (pi * radius^3)
  
  result <- data.frame(observed_cross_K = observed_cross_K,
                       expected_cross_K = expected_cross_K)
  
  return(result)
}