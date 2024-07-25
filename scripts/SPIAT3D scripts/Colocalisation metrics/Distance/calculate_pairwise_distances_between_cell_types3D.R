
calculate_pairwise_distances_between_cell_types3D <- function(spe,
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
  
  # If there are less than two cells, give error
  if (nrow(df) <= 1) stop("There must be at least two cells in spe")
  
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
  
  # Create a list of the number of cell types with their
  # corresponding cell ID's
  cell_types <- list()
  for (cell_type in unique(df[["Cell.Type"]])) {
    cell_types[[cell_type]] <- as.character(df$Cell.ID[df[["Cell.Type"]] == cell_type])
  }
  
  # Calculate cell to cell distances
  dist_all <- -1 * apcluster::negDistMat(df[, c("Cell.X.Position",
                                                "Cell.Y.Position",
                                                "Cell.Z.Position")])
  
  cell_id_vector <- df$Cell.ID
  colnames(dist_all) <- cell_id_vector
  rownames(dist_all) <- cell_id_vector
  
  cell_to_cell_dist_all <- vector()

  for (i in seq(length(cell_types))) {
    
    for (j in i:length(cell_types)) {
  
      cell_name1 <- names(cell_types)[i]
      cell_name2 <- names(cell_types)[j]

      cell_ids1 <- cell_types[[cell_name1]]
      cell_ids2 <- cell_types[[cell_name2]]
      
      ## Same cell type, only one cell
      if (cell_name1 == cell_name2 && length(cell_ids1) == 1) next
        
      cell_to_cell <- dist_all[cell_id_vector %in% cell_ids1, 
                               cell_id_vector %in% cell_ids2]
      
      ## Different cell types, each only has one cell
      if (length(cell_ids1) == 1 && length(cell_ids2) == 1) {
        cell_to_cell <- as.matrix(cell_to_cell)
        rownames(cell_to_cell) <- cell_ids1
        colnames(cell_to_cell) <- cell_ids2
      }    
      ## Different cell types, only one cell of cell_type1
      else if (length(cell_ids1) == 1) {
        cell_to_cell <- as.matrix(cell_to_cell)
        colnames(cell_to_cell) <- cell_ids1
      }
      ## Different cell types, only one cell of cell_type2
      else if (length(cell_ids2) == 1) {
        cell_to_cell <- as.matrix(cell_to_cell)
        colnames(cell_to_cell) <- cell_ids2
      }
      
      ## Same cell type, only need part of the matrix
      if (cell_name1 == cell_name2) cell_to_cell[upper.tri(cell_to_cell, diag = TRUE)] <- NA
      
      # Melts dist_all to produce dataframe of target and nearest 
      # cell ID's columns and distance column
      cell_to_cell_dist <- reshape2::melt(cell_to_cell, na.rm = TRUE)
      cell_to_cell_dist$cell_type1 <- cell_name1
      cell_to_cell_dist$cell_type2 <- cell_name2
      cell_to_cell_dist$pair <- paste(cell_name1, cell_name2, sep="/")

      cell_to_cell_dist_all <- rbind(cell_to_cell_dist_all, 
                                     cell_to_cell_dist)
    }
  }
  
  colnames(cell_to_cell_dist_all)[c(1,2,3)] <- c("cell_id1", "cell_id2", "distance")
 
  # Plot
  if (plot_image) {
    fig <- plot_cell_distances_violin3D(cell_to_cell_dist_all)
    methods::show(fig)
  }
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types3D(cell_to_cell_dist_all))  
  }
  
  return(cell_to_cell_dist_all)
}
