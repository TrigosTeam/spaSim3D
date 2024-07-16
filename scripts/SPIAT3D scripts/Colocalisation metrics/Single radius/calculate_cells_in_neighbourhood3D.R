calculate_cells_in_neighbourhood3D <- function(spe, 
                                               reference_cell_type, 
                                               target_cell_types, 
                                               radius, 
                                               feature_colname = "Cell.Type",
                                               show_summary = TRUE,
                                               plot_image = TRUE) {
  
  if (is.null(spe[["Cell.ID"]])) stop("No Cell.ID column. Add a Cell.ID columnt to your spe.")
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), 
                   "Cell.Type" = spe[[feature_colname]], 
                   "Cell.ID" = spe[["Cell.ID"]])

  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% df$Cell.Type)) {
    stop(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
  }
  
  ## For target_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, df$Cell.Type)
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  # Check if radius is numeric
  if (!is.numeric(radius)) {
    stop(paste(radius, " is not of type 'numeric'"))
  }
  
  ## Get data for reference cells
  reference_cells <- df[df[[feature_colname]] == reference_cell_type, ]
  reference_cell_coords <- reference_cells[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  rownames(reference_cell_coords) <- reference_cells$Cell.ID
  
  result <- data.frame(matrix(nrow = nrow(reference_cells), ncol = 0))
  
  for (target_cell_type in target_cell_types) {
    ## Get df for target cells
    target_cells <- df[df[[feature_colname]] == target_cell_type, ]
    target_cell_coords <- target_cells[, c("Cell.X.Position","Cell.Y.Position", "Cell.Z.Position")]
    
    ## Determine number of target cells specified distance for each reference cell
    reference_target_result <- dbscan::frNN(target_cell_coords, 
                                            eps = radius,
                                            query = reference_cell_coords, 
                                            sort = FALSE)
    n_targets <- rapply(reference_target_result$id, length)
    
    ## Add to data frame
    result[[target_cell_type]] <- n_targets
  }
  
  result <- data.frame(ref_cell_id = reference_cells$Cell.ID, result)
  
  if (show_summary) {
    ## Show summarised results
    print(summarise_cells_in_neighbourhood3D(result))    
  }

  
  ## Plot
  if (plot_image) {
    fig <- plot_cells_in_neighbourhood_violin3D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}
