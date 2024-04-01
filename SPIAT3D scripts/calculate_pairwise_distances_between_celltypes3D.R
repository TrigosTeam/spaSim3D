## data is a dataframe with colnames:
## "Cell.X.Position" "Cell.Y.Position" "Cell.Z.Position" "Cell.Type" "Cell.ID"     

calculate_pairwise_distances_between_celltypes3D <- function(
    data,
    cell_type_of_interest = NULL,
    feature_colname = "Cell.Type") {
 
  # Select all rows in data which only contains the cells of interest
  if (!is.null(cell_type_of_interest)) {
    data <- data[data[ , feature_colname] %in% cell_type_of_interest, ]
  }
   
  
  # If there are no cells which match cell_type_of_interest, give error:
  if (nrow(data) == 0) {
    stop("There are no cells or no cells of specified cell types")
  }
  
  # Create a list of the number of cell types with their
  # corresponding cell ID's
  cell_types <- list()
  for (eachType in unique(data[ , feature_colname])) {
    cell_types[[eachType]] <- as.character(data$Cell.ID[data[ , feature_colname] == eachType])
  }
  
  # Calculate cell to cell distances
  dist_all <- -1 * apcluster::negDistMat(data[ c("Cell.X.Position",
                                                 "Cell.Y.Position",
                                                 "Cell.Z.Position")])
  
  cell_id_vector <- data$Cell.ID
  colnames(dist_all) <- cell_id_vector
  rownames(dist_all) <- cell_id_vector
  
  cell_to_cell_dist_all <- vector()

  for (i in seq(length(cell_types))) {
    
    for (j in i:length(cell_types)) {
  
      cell_name1 <- names(cell_types)[i]
      cell_name2 <- names(cell_types)[j]

      cell_ids1 <- cell_types[[cell_name1]]
      cell_ids2 <- cell_types[[cell_name2]]
      
      ## Need to investigate this
      if (length(cell_ids1) < 2 & length(cell_ids2) < 2) {
        next
      }
        
      cell_to_cell <- dist_all[cell_id_vector %in% cell_ids1, 
                               cell_id_vector %in% cell_ids2]
      
      if (cell_name1 == cell_name2) {
        cell_to_cell[upper.tri(cell_to_cell, diag = TRUE)] <- NA
      }
      
      # Melts dist_all to produce dataframe of target and nearest 
      # cell ID's columns and distance column
      cell_to_cell_dist <- reshape2::melt(cell_to_cell, na.rm = TRUE)
      cell_to_cell_dist$Type1 <- cell_name1
      cell_to_cell_dist$Type2 <- cell_name2
      cell_to_cell_dist$Pair <- paste(cell_name1, cell_name2, sep="/")
      
      cell_to_cell_dist_all <- rbind(cell_to_cell_dist_all, 
                                     cell_to_cell_dist)
    }
  }
  
  colnames(cell_to_cell_dist_all)[c(1,2,3)] <- c("Cell1", "Cell2", "Distance")
 
  return (cell_to_cell_dist_all)
}
