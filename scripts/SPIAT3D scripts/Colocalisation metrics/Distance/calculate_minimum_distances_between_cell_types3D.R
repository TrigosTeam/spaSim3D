## Please ensure there is no factoring in any of the columns!!!

calculate_minimum_distances_between_cell_types3D <- function(spe,
                                                             cell_types_of_interest = NULL,
                                                             feature_colname = "Cell.Type",
                                                             show_summary = TRUE,
                                                             plot_image = TRUE) {
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.Id column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), 
                   "Cell.Type" = spe[[feature_colname]], 
                   "Cell.ID" = spe[["Cell.ID"]])
  
  # If there are no cells, give error
  if (nrow(df) == 0) stop("There are no cells in spe")
  
  # Select all rows in data frame which only contains the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, df$Cell.Type)
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                 paste(unknown_cell_types, collapse = ", ")))
    }
    
    df <- df[df[["Cell.Type"]] %in% cell_types_of_interest, ]
  }
  
  # Create a list of the number of cell types with their corresponding cell ID's
  cell_types <- list()
  for (cell_type in unique(df[["Cell.Type"]])) {
    cell_types[[cell_type]] <- as.character(df$Cell.ID[df[["Cell.Type"]] == cell_type])
  }
  
  # Get different possible cell type combinations
  # Each row represents a combination
  # If a row is [1 , 2], then we are comparing cell type 1 and cell type 2
  unique_cells <- unique(df[["Cell.Type"]]) # unique cell types
  permu <- gtools::permutations(length(unique_cells), 2, repeats.allowed = TRUE)

  result <- vector()
  
  for (i in seq(nrow(permu))) {
    name1 <- unique_cells[permu[i, 1]]
    name2 <- unique_cells[permu[i, 2]]
    
    # Get x,y,z coords for all cells of cell_type1 and cell_type2
    all_cell_type1_coord <- df[df[, "Cell.Type"] == name1, 
                               c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
    
    all_cell_type2_coord <- df[df[, "Cell.Type"] == name2, 
                                c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
    
    # Find all of closest points
    # For each cell of cell_type1, find the closest cell of cell_type2
    if (name1 != name2) {
      all_closest <- RANN::nn2(data = all_cell_type2_coord, 
                               query = all_cell_type1_coord, 
                               k = 1)  
    }
    else {
      # If we are comparing the same cell_type, use the second closest neighbour
      all_closest <- RANN::nn2(data = all_cell_type2_coord, 
                               query = all_cell_type1_coord, 
                               k = 2)
      all_closest[['nn.idx']] <- all_closest[['nn.idx']][, 2]
      all_closest[['nn.dists']] <- all_closest[['nn.dists']][, 2]
    }
    
    # Create the data frame containing the chosen cells and their ids, as well as
    # the nearest cell to them and their ids, and the distance between
    cell_type2_cell_IDs <- df[df[ , "Cell.Type"] == name2, "Cell.ID"]
    
    local_dist_mins <- data.frame(
      ref_cell_id = cell_types[[name1]],
      ref_cell_type = name1,
      nearest_cell_id = cell_type2_cell_IDs[as.vector(all_closest$nn.idx)],
      nearest_cell_type = name2,
      distance = all_closest$nn.dists
    )
    result <- rbind(result, local_dist_mins)
  }
  
  result$pair <- paste(result$ref_cell_type, result$nearest_cell_type,sep = "/")
  
  # Plot
  if (plot_image) {
    fig <- plot_cell_distances_violin3D(result)
    methods::show(fig)
  }
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types3D(result))  
  }
  
  return(result)
}
