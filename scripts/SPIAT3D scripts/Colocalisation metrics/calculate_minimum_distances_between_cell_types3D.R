## data is a dataframe with colnames:
## "Cell.X.Position" "Cell.Y.Position" "Cell.Z.Position" "Cell.Type" "Cell.ID"

## Please ensure there is no factoring in any of the columns!!!

calculate_minimum_distances_between_cell_types3D <- function(data,
                                                             cell_types_of_interest = NULL,
                                                             feature_colname = "Cell.Type") {
  
  
  # If the columns are not correct, give error
  required_colnames <- c("Cell.X.Position", 
                         "Cell.Y.Position", 
                         "Cell.Z.Position", 
                         feature_colname, 
                         "Cell.ID")
  
  missing_colnames <- setdiff(required_colnames,
                              colnames(data))
  
  if (length(missing_colnames) > 0) {
    stop(paste(paste(missing_colnames, collapse = ', '),
               "are missing as column names in your data")) 
  }
  
  # If there are no cells, give error
  if (nrow(data) == 0) {
    stop("There are no cells in data")
  }
  
  
  # Select all rows in data which only contains the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    # Check if cell_types_of_interest has cells not found in the data
    incorrect_cell_types <- setdiff(cell_types_of_interest, unique(data[[feature_colname]]))
    if (length(incorrect_cell_types) > 0) {
      stop(paste(paste(incorrect_cell_types, collapse = ', '),
                 "in cell_types_of_interest don't exist."))
    }
    
    data <- data[data[ , feature_colname] %in% cell_types_of_interest, ]
  }
  
  # Create a list of the number of cell types with their
  # corresponding cell ID's
  cell_types <- list()
  for (eachType in unique(data[ , feature_colname])) {
    cell_types[[eachType]] <- as.character(data$Cell.ID[data[ , feature_colname] == eachType])
  }
  
  # Get different possible cell type combinations
  # Each row represents a combination
  # If a row is [1 , 2], then we are comparing cell type 1 and cell type 2
  unique_cells <- unique(data[[feature_colname]]) # unique cell types
  permu <- gtools::permutations(length(unique_cells), 2, repeats.allowed = TRUE)

  result <- vector()
  
  for (i in seq(nrow(permu))) {
    name1 <- unique_cells[permu[i, 1]]
    name2 <- unique_cells[permu[i, 2]]
    
    # Get x,y,z coords for all cells of cell_type1 and cell_type2
    all_cell_type1_coord <- data[data[, feature_colname] == name1, 
                               c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
    
    all_cell_type2_coord <- data[data[, feature_colname] == name2, 
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
    
    # Create the data.frame containing the chosen cells and their ids, as well as
    # the nearest cell to them and their ids, and the distance between
    cell_type2_cell_IDs <- data[data[ , feature_colname] == name2, "Cell.ID"]
    
    local_dist_mins <- data.frame(
      RefCell = cell_types[[name1]],
      RefType = name1,
      NearestCell = cell_type2_cell_IDs[as.vector(all_closest$nn.idx)],
      NearestType = name2,
      Distance = all_closest$nn.dists
    )

    result <- rbind(result, local_dist_mins)
    
  }
  
  result$Pair <- paste(result$RefType, result$NearestType,sep = "/")
  
  return (result)
}
