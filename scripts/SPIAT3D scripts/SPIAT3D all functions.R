calculate_cell_proportions3D <- function(data,
                                         reference_cell_types = NULL, 
                                         cell_types_to_exclude = NULL, 
                                         feature_colname = "Cell.Type",
                                         plot.image = TRUE) {
  
  # If the columns are not correct, give error
  required_colnames <- c("Cell.X.Position", 
                         "Cell.Y.Position", 
                         "Cell.Z.Position", 
                         feature_colname)
  
  missing_colnames <- setdiff(required_colnames,
                              colnames(data))
  
  if (length(missing_colnames) > 0) {
    stop(paste(paste(missing_colnames, collapse = ', '),
               "are missing as column names in your data")) 
  }
  
  # Check
  if (nrow(data) == 0) {
    stop("No cells found for calculating cell proportions")
  }
  
  # Creates frequency/bar plot of all cell types in the entire image
  cell_proportions <- as.data.frame(table(data[, feature_colname]))
  names(cell_proportions) <- c("Cell.Type", 'Frequency')
  
  # Exclude any cell types not wanted
  if (!is.null(cell_types_to_exclude)) {
    
    # Check
    incorrect_cell_types <- setdiff(cell_types_to_exclude, cell_proportions$Cell.Type)
    if (length(incorrect_cell_types) > 0) {
      stop(paste(paste(incorrect_cell_types, collapse = ', '),
                 "in cell_types_to_exclude don't exist."))
    }
    
    cell_proportions <- cell_proportions[!(cell_proportions$Cell.Type %in% cell_types_to_exclude), ]
    
    # Check
    if (nrow(cell_proportions) == 0) {
      stop("All cells have been excluded")
    }
  }
  
  # Find proportion of each cell type against all cells
  if (is.null(reference_cell_types)) {
    
    # Get frequency total for all cells
    cell_type_frequency_total <- sum(cell_proportions$Frequency)
    
    cell_proportions$Proportion <- cell_proportions$Frequency / cell_type_frequency_total
    cell_proportions$Percentage <- cell_proportions$Proportion * 100
    cell_proportions$Proportion_Name <- "/Total"
  }
  # Find proportion of each cell type against the chosen reference cell types
  else {
    
    # Check
    incorrect_cell_types <- setdiff(reference_cell_types, cell_proportions$Cell.Type)
    if (length(incorrect_cell_types) > 0) {
      stop(paste(paste(incorrect_cell_types, collapse = ', '),
                 "in reference_cell_types have been excluded or don't exist."))
    }
    
    # Get frequency total for chosen reference cells
    cell_type_frequency_total <- sum(cell_proportions$Frequency[cell_proportions[['Cell.Type']] %in% reference_cell_types])
    
    cell_proportions$Proportion <- cell_proportions$Frequency/cell_type_frequency_total
    cell_proportions$Percentage <- cell_proportions$Proportion * 100
    cell_proportions$Proportion_Name <- "/Custom"  
    cell_proportions$Reference <- paste(reference_cell_types, collapse=",")
  }
  
  # Order by Reference cell type (reverse to have Total first if present) then by highest proportion
  cell_proportions <- cell_proportions[rev(order(cell_proportions$Proportion)), ]
  
  if (plot.image) {
    g <- ggplot(cell_proportions, aes(x=Cell.Type, y=Percentage)) +
      geom_bar(stat='identity') + theme_bw()
    methods::show(g)
  }
  
  return (cell_proportions)
}


summarise_distances_between_cell_types3D <- function(df) {
  
  Pair <- Distance <- NULL
  
  # summarise the results
  summarised_dists <- df %>% 
    dplyr::group_by(Pair) %>%
    dplyr::summarise(mean(Distance, na.rm = TRUE), 
                     min(Distance, na.rm = TRUE), 
                     max(Distance, na.rm = TRUE),
                     stats::median(Distance, na.rm = TRUE), 
                     stats::sd(Distance, na.rm = TRUE))
  
  summarised_dists <- data.frame(summarised_dists)
  
  colnames(summarised_dists) <- c("Pair", 
                                  "Mean", 
                                  "Min", 
                                  "Max", 
                                  "Median", 
                                  "Std.Dev")
  
  for (i in seq(nrow(summarised_dists))) {
    # Get cell_types for each pair
    cell_types <- strsplit(summarised_dists[i,"Pair"], "/")[[1]]
    
    summarised_dists[i, "Reference"] <- cell_types[1]
    summarised_dists[i, "Target"] <- cell_types[2]
  }
  
  return(summarised_dists)
}


plot_cell_percentages_bar3D <- function(cell_proportions) {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Cell.Type <- Percentage <- Percentage_label <- NULL
  
  cell_proportions$Percentage_label <- round(cell_proportions$Percentage, digits=1)
  
  cell_percentages_full_plot <-
    ggplot(cell_proportions,
           aes(x = stats::reorder(Proportion_Name, Percentage), 
               y = Percentage, fill = Cell.Type)) +
    geom_bar(stat = 'identity') +
    theme_bw() +
    theme() +
    xlab("Cell Type") + 
    ylab("Proportion of cells") +
    geom_text(aes(label = Percentage_label), 
              position = position_stack(vjust = 0.5), size = 3) +
    coord_flip() 
  
  return (cell_percentages_full_plot)
}



## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_cell_distances_violin3D <- function(cell_to_cell_dist, scales = "free_x") {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Pair <- Distance <- NULL
  
  ggplot(cell_to_cell_dist, aes(x = Pair, y = Distance)) + geom_violin() +
    facet_wrap(~Pair, scales=scales) +
    theme_bw() +
    theme(axis.text.x=element_blank())
  
}



plot_cell_distances_summary_heatmap3D <- function(distances_summary_df, 
                                                  metric = "Mean"){
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Reference <- Target <- Min <- Max <- Mean <- Std.Dev <- Median <- NULL
  
  if (metric == "Mean") {
    limit <- range(unlist(distances_summary_df$Mean), na.rm=TRUE)
    g <- ggplot(distances_summary_df, aes(x = Reference, y = Target, fill = Mean))
    
  }
  else if (metric == "Std.Dev") {
    limit <- range(unlist(distances_summary_df$Std.Dev), na.rm=TRUE)
    g <- ggplot(distances_summary_df, aes(x = Reference, y = Target, fill = Std.Dev))
    
  }
  else if (metric == "Median") {
    limit <- range(unlist(distances_summary_df$Median), na.rm=TRUE)
    g <- ggplot(distances_summary_df, aes(x = Reference, y = Target, fill = Median))
  }
  else if (metric == "Min") {
    limit <- range(unlist(distances_summary_df$Min), na.rm=TRUE)
    g <- ggplot(distances_summary_df, aes(x = Reference, y = Target, fill = Min))
  }
  else if (metric == "Max") {
    limit <- range(unlist(distances_summary_df$Max), na.rm=TRUE)
    g <- ggplot(distances_summary_df, aes(x = Reference, y = Target, fill = Max))
  }
  else {
    stop(paste(metric," is not a valid metric"))
  }
  
  g <- g +
    geom_tile() +
    xlab("Reference cell type") +
    ylab("Target cell type") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_viridis_c(limits = limit, direction = -1)
  
  return (g)
}



## data is a dataframe with colnames:
## "Cell.X.Position" "Cell.Y.Position" "Cell.Z.Position" "Cell.Type" "Cell.ID"     

calculate_pairwise_distances_between_cell_types3D <- function(data,
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
    cell_types[[eachType]] <- as.character(data$Cell.ID[data[, feature_colname] == eachType])
  }
  
  # Calculate cell to cell distances
  dist_all <- -1 * apcluster::negDistMat(data[, c("Cell.X.Position",
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


calculate_mixing_scores3D <- function(data, 
                                      reference_cell_types, 
                                      target_cell_types, 
                                      radius = 20, 
                                      feature_colname = "Cell.Type") {
  
  # If the columns are not correct, give error
  required_colnames <- c("Cell.X.Position", 
                         "Cell.Y.Position", 
                         "Cell.Z.Position", 
                         feature_colname)
  
  missing_colnames <- setdiff(required_colnames,
                              colnames(data))
  
  if (length(missing_colnames) > 0) {
    stop(paste(paste(missing_colnames, collapse = ', '),
               "are missing as column names in your data")) 
  }
  
  
  # Check if reference_cell_types has cells not found in the data
  incorrect_cell_types <- setdiff(reference_cell_types, unique(data[[feature_colname]]))
  if (length(incorrect_cell_types) > 0) {
    stop(paste(paste(incorrect_cell_types, collapse = ', '),
               "in reference_cell_types don't existin data."))
  }
  
  # Check if target_cell_types has cells not found in the data
  incorrect_cell_types <- setdiff(target_cell_types, unique(data[[feature_colname]]))
  if (length(incorrect_cell_types) > 0) {
    stop(paste(paste(incorrect_cell_types, collapse = ', '),
               "in target_cell_types don't exist in data."))
  }
  
  # Check if radius is numeric
  if (!is.numeric(radius)) {
    stop(paste(radius, " is not of type 'numeric'"))
  }
  
  
  df <- data.frame(matrix(ncol=8, nrow=0))
  
  for (reference_cell_type in reference_cell_types) {
    
    # Get all info for cells of reference cell_type
    reference_cells <- data[data[, feature_colname] == reference_cell_type, ]
    
    for (target_cell_type in target_cell_types) {
      
      # Get all info for cells of target cell_type      
      target_cells <- data[data[, feature_colname] == target_cell_type, ]
      
      # No point getting mixing scores if comparing the same cell type
      if (reference_cell_type == target_cell_type) {
        next
      }
      
      # Can't get mixing scores if there are no reference cells
      if (nrow(reference_cells) == 0) {
        methods::show(paste("There are no unique reference cells of specified celltype", reference_cell_type, "for target cell", target_cell_type))
        df <-  rbind(df, 
                     c(reference_cell_type, 
                       target_cell_type, 
                       0, 
                       nrow(target_cells), 
                       0, 
                       0, 
                       NA, 
                       NA))
      }
      
      # Can't get mixing scores if there are no target cells
      else if (nrow(target_cells) == 0) {
        methods::show(paste("There are no unique target cells of specified celltype", target_cell_type, "for reference cell", reference_cell_type))
        
        reference_cell_coords <- reference_cells[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
        reference_reference_result <- dbscan::frNN(reference_cell_coords, 
                                                   eps = radius, 
                                                   query = NULL,
                                                   sort = FALSE)
        
        # halve it to avoid counting each ref-ref interaction twice
        reference_reference_interactions <- 0.5 * sum(rapply(reference_reference_result$id, length)) 
        
        df <-  rbind(df[ , df.cols], 
                     c(reference_cell_type, 
                       target_cell_type, 
                       nrow(reference_cells), 
                       0, 
                       0, 
                       reference_reference_interactions, 
                       NA, 
                       NA))
      }
      
      # Generic case: We have reference cells and target cells
      else {
        # Get x,y,z coords for reference cells and target cells
        reference_cell_coords <- reference_cells[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
        target_cell_coords <- target_cells[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
        
        # For each reference cell, find all target cells within the chosen radius
        reference_target_result <- dbscan::frNN(target_cell_coords, 
                                                eps = radius, 
                                                query = reference_cell_coords, 
                                                sort = FALSE)
        
        # Find the total sum of how many target cells were close enough to reference cells
        reference_target_interactions <- sum(rapply(reference_target_result$id, length))
        
        # For each reference cell, find all other reference cells within the chosen radius
        reference_reference_result <- dbscan::frNN(reference_cell_coords, 
                                                   eps = radius,
                                                   query = NULL,
                                                   sort = FALSE)
        
        # Find the the total sum of how many other reference cells were close enough to reference cells
        # Halve it to avoid counting each ref-ref interaction twice
        reference_reference_interactions <- 0.5 * sum(rapply(reference_reference_result$id, length)) 
        
        
        if (reference_reference_interactions != 0) {
          mixing_score <- reference_target_interactions / reference_reference_interactions
          normalised_mixing_score <- 0.5 * mixing_score * (nrow(reference_cells) - 1) / nrow(target_cells)
        }
        else {
          mixing_score <- 0
          normalised_mixing_score <- 0
          methods::show(paste("There are no reference to reference interactions for", target_cell_type, "in the specified radius, cannot calculate mixing score"))
        }
        
        df <-  rbind(df, 
                     c(reference_cell_type, 
                       target_cell_type, 
                       nrow(reference_cells), 
                       nrow(target_cells), 
                       reference_target_interactions, 
                       reference_reference_interactions, 
                       mixing_score, 
                       normalised_mixing_score))
      }
    }
  }
  
  # Required column names of our output data frame
  df.cols <- c("Reference", 
               "Target", 
               "Number_of_reference_cells",
               "Number_of_target_cells", 
               "Reference_target_interaction",
               "Reference_reference_interaction", 
               "Mixing_score", 
               "Normalised_mixing_score")
  colnames(df) <- df.cols
  
  return(df)
}



calculate_mixing_scores_gradient3D <- function(data, 
                                               reference_cell_type, 
                                               target_cell_type, 
                                               radii = 20, 
                                               feature_colname = "Cell.Type",
                                               plot_image = TRUE) {
  
  
  df <- data.frame(matrix(nrow = radii, ncol = 8))
  df.cols <- c("Reference", 
               "Target", 
               "Number_of_reference_cells",
               "Number_of_target_cells", 
               "Reference_target_interaction",
               "Reference_reference_interaction", 
               "Mixing_score", 
               "Normalised_mixing_score")
  colnames(df) <- df.cols
  
  for (radius in seq(radii)) {
    mixing_scores <- calculate_mixing_scores3D(data,
                                               reference_cell_type,
                                               target_cell_type,
                                               radius,
                                               feature_colname)
    
    df[radius, ] <- mixing_scores
  }
  
  
  if (plot_image) {
    plot(seq(radii), df[["Normalised_mixing_score"]], type = "l", xlab = "Radius", ylab = "Normalised Mixing Score")
    abline(a = 1, b = 0, col = "red", lwd = 2, lty = 2)
  }
  
  return (df)
}



calculate_cells_in_neighborhood3D <- function(data, 
                                              reference_cell_types, 
                                              target_cell_types, 
                                              radius = 20, 
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
  
  
  # Check if reference_cell_types has cells not found in the data
  incorrect_cell_types <- setdiff(reference_cell_types, unique(data[[feature_colname]]))
  if (length(incorrect_cell_types) > 0) {
    stop(paste(paste(incorrect_cell_types, collapse = ', '),
               "in reference_cell_types don't existin data."))
  }
  
  # Check if target_cell_types has cells not found in the data
  incorrect_cell_types <- setdiff(target_cell_types, unique(data[[feature_colname]]))
  if (length(incorrect_cell_types) > 0) {
    stop(paste(paste(incorrect_cell_types, collapse = ', '),
               "in target_cell_types don't exist in data."))
  }
  
  # Check if radius is numeric
  if (!is.numeric(radius)) {
    stop(paste(radius, " is not of type 'numeric'"))
  }
  
  
  
  result <- list()
  
  for (reference_cell_type in reference_cell_types) {
    ## Get data for reference cells
    reference_cells <- data[which(data[, feature_colname] == reference_cell_type), ]
    reference_cell_coords <- reference_cells[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
    
    ## Set up data frame for current reference cell type
    reference_cell_df <- reference_cells[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
    rownames(reference_cell_df) <- reference_cells$Cell.ID
    
    for (target_cell_type in target_cell_types) {
      ## Get data for target cells
      target_cells <- data[which(data[, feature_colname] == target_cell_type), ]
      target_cell_coords <- target_cells[, c("Cell.X.Position","Cell.Y.Position", "Cell.Z.Position")]
      
      ## Determine number of target cells specified distance for each reference cell
      reference_target_result <- dbscan::frNN(target_cell_coords, 
                                              eps = radius,
                                              query = reference_cell_coords, 
                                              sort = FALSE)
      n_targets <- rapply(reference_target_result$id, length)
      
      ## Add to results data frame
      reference_cell_df[, target_cell_type] <- n_targets
      
    }
    result[[reference_cell_type]] <- reference_cell_df
  }
  
  return (result)
}



## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_cells_in_neighborhood_violin3D <- function(cells_in_neighborhood, scales = "free_x") {
  
  df <- data.frame(matrix(nrow = 0, ncol = 3))
  colnames(df) <- c("Target", "Count", "Reference")
  
  for (i in seq(length(cells_in_neighborhood))) {
    
    # Get reference cell type for current index
    reference_cell_type <- names(cells_in_neighborhood)[i]
    
    # Get data for current index
    cells_in_neighborhood_df <- cells_in_neighborhood[[i]]
    
    # Get columns which contain cell count data (4th column onwards)
    cells_in_neighborhood_df <- cells_in_neighborhood_df[ , 4:ncol(cells_in_neighborhood_df)]
    
    # Melt
    cells_in_neighborhood_df <- reshape2::melt(cells_in_neighborhood_df, id.vars = 0)
    colnames(cells_in_neighborhood_df) <- c("Target", "Count")
    
    # Add reference cell type column
    cells_in_neighborhood_df$Reference <- reference_cell_type
    
    # Add result to main df
    df <- rbind(df, cells_in_neighborhood_df)
    
  }
  
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Reference <- Count <- Target <- NULL
  
  ggplot(df, aes(x = Reference, y = Count, fill = Target)) + geom_violin() +
    facet_wrap(~Reference, scales=scales) +
    theme_bw() +
    theme(axis.text.x=element_blank())
  
}


summarise_cells_in_neighborhood3D <- function(cells_in_neighborhood_data) {
  
  summarised_results <- list()
  
  for (i in seq(length(cells_in_neighborhood_data))) {
    df <- cells_in_neighborhood_data[[i]]
    
    ## Target cell types will be the fourth column onwards
    target_cell_types <- colnames(df)[4:ncol(df)]
    
    ## Set up data frame for summarised_results list
    df_results <- data.frame(row.names = c("Mean", "Min", "Max", "Median", "St.Dev"))
    
    for (target_cell_type in target_cell_types) {
      
      ## Get statistical measures for each target cell type
      target_cell_type_values <- df[[target_cell_type]]
      df_results[[target_cell_type]] <- c(mean(target_cell_type_values),
                                          min(target_cell_type_values),
                                          max(target_cell_type_values),
                                          median(target_cell_type_values),
                                          sd(target_cell_type_values))
      
      
    }
    
    ## Add data frame result to summarised_results
    summarised_results[[names(cells_in_neighborhood_data)[i]]] <- df_results
  }
  
  return (summarised_results)
}



calculate_Kcross3D <- function(data, 
                               reference_cell_type,
                               target_cell_type,
                               distance,
                               feature_colname = "Cell.Type", 
                               plot_results = TRUE) {
  
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
  
  
  # Check if reference_cell_type is in the data
  if (!reference_cell_type %in% unique(data[[feature_colname]])) {
    stop(paste(reference_cell_type, " reference_cell_type does not exist in data"))
  }
  
  # Check if target_cell_type is in the data
  if (!target_cell_type %in% unique(data[[feature_colname]])) {
    stop(paste(target_cell_type, " target_cell_type does not exist in data"))
  }
  
  
  # Check if distance is numeric
  if (!is.numeric(distance)) {
    stop(paste(distance, " is not of type 'numeric'"))
  }
  
  
  
  ## Get x, y, z coords for cells of reference cell type and target cell type
  reference_cell_data <- data[data[[feature_colname]] == reference_cell_type, ]
  target_cell_data <- data[data[[feature_colname]] == target_cell_type, ]
  
  ## Get rough dimensions of the window the points are in
  length <- round(max(data$Cell.X.Position) - min(data$Cell.X.Position))
  width  <- round(max(data$Cell.Y.Position) - min(data$Cell.Y.Position))
  height <- round(max(data$Cell.Z.Position) - min(data$Cell.Z.Position))
  
  ## Get distance of each reference cell from the border of the window
  ## Defined as the minimum of the x, y and z distance
  reference_cell_data$Cell.X.Position.Inverse <- length - reference_cell_data$Cell.X.Position
  reference_cell_data$Cell.Y.Position.Inverse <- width  - reference_cell_data$Cell.Y.Position
  reference_cell_data$Cell.Z.Position.Inverse <- height - reference_cell_data$Cell.Z.Position
  reference_cell_data$Border.Distance <- as.numeric(apply(reference_cell_data, 1, min))
  
  ## Only use reference cells which a far enough away from the border
  if (sum(reference_cell_data$Border.Distance > distance) == 0) {
    stop(paste(distance, " distance is too large, all reference cells have been deleted"))
  }
  else if (sum(reference_cell_data$Border.Distance > distance) < nrow(reference_cell_data)) {
    warning("Some reference cells will be ignored due to border effects, consider reducing distance parameter")
  }
  
  chosen_rows <- rownames(reference_cell_data[reference_cell_data$Border.Distance > distance, ])
  reference_cell_data <- reference_cell_data[chosen_rows, c("Cell.X.Position",
                                                            "Cell.Y.Position",
                                                            "Cell.Z.Position",
                                                            "Cell.Type",
                                                            "Cell.ID")]
  
  ## Get number of reference and target cells
  n_reference_cells <- nrow(reference_cell_data)
  n_target_cells <- nrow(target_cell_data)
  
  
  ## Combine together
  combined_cell_data <- rbind(reference_cell_data, target_cell_data)
  
  ## Get distances between chosen cell types
  reference_target_distances <- -1 * apcluster::negDistMat(combined_cell_data[, c("Cell.X.Position",
                                                                                  "Cell.Y.Position",
                                                                                  "Cell.Z.Position")])
  
  # Only need distances between reference cells and target cells
  # Ignore ref-ref or target-target distances (i.e. top right part of matrix)
  reference_target_distances <- reference_target_distances[1:n_reference_cells, 
                                                           (n_reference_cells + 1):ncol(reference_target_distances)]
  
  # Calculate observed cross-k value for a sequence of distances
  # i.e. the number of ref-target distances less than the chosen distance
  distances <- seq(1, distance, 0.25)
  observed_k <- unlist(lapply(distances, function(x) sum(reference_target_distances < x)))
  
  # Get volume of the window the cells are in
  volume <- length * width * height
  
  # Calculate expected cross-k value for a sequence of distances (using the formula?)
  expected_k <- n_reference_cells * n_target_cells * ((4/3) * pi * distances^3) / volume
  
  result <- data.frame(Distance = distances,
                       Observed = observed_k,
                       Expected = expected_k)
  
  if (plot_results) {
    plot(result$Distance, result$Observed, type = "l", col = "red", 
         xlim = c(0, distance), ylim = c(0, max(result)),
         xlab = "Distance", ylab = "Cross K-function Value")
    lines(result$Distance, result$Expected, type = "l", col = "blue", lty = 2)
    legend(0, max(result), legend = c("Observed K", "Expected K"), col = c("red", "blue"), lty = c(1, 2))
  }
  
  return (result)
}


calculate_Kcross_intersection3D <- function(Kcross_df) {
  
  Kcross_df$sign <- Kcross_df$Observed - Kcross_df$Expected
  
  ## Determine when sign flips from -ve to +ve OR +ve to -ve
  change_of_sign <- diff(sign(Kcross_df$sign))
  
  ## Determine indices for when observed curve goes above or below the expected curve
  observed_goes_above_indices <- which(change_of_sign == 2)
  observed_goes_below_indices <- which(change_of_sign == -2)
  
  if (length(observed_goes_above_indices) + length(observed_goes_below_indices) == 0) {
    warning("No cross-K intersections occur")
    return (0)
  }
  
  above_distance <- c()
  below_distance <- c()
  
  if (length(observed_goes_above_indices) != 0) {
    above_distance <- Kcross_df$Distance[observed_goes_above_indices]
    print("The observed curve goes ABOVE the expected curve at the following distances:")
    print(paste(above_distance))
  }
  if (length(observed_goes_below_indices) != 0) {
    below_distance <- Kcross_df$Distance[observed_goes_below_indices]
    print("The observed curve goes BELOW the expected curve at the following distances:")
    print(paste(below_distance))
  }
  
  result <- data.frame(Distance = c(above_distance, below_distance),
                       Change = c(rep("Observed goes above Expected", length(above_distance)),
                                  rep("Observed goes below Expected", length(below_distance)))) 
  
  return (result)
}



calculate_AUC_of_Kcross3D <- function(Kcross_df) {
  
  ## Get difference in area between 
  AUC <- pracma::trapz(Kcross_df$Distance, Kcross_df$Observed) - 
    pracma::trapz(Kcross_df$Distance, Kcross_df$Expected)
  
  ## Get the cross-k result image size
  max_distance <- max(Kcross_df$Distance)
  max_Kcross <- max(c(Kcross_df$Observed, Kcross_df$Expected))
  
  ## Calculate normalised AUC
  n_AUC <- AUC / (max_distance * max_Kcross)
  
  return (n_AUC)
}



calculate_entropy3D <- function(data,
                                radius = NULL,
                                reference_cell_type = NULL,
                                target_cell_types,
                                log_base = NULL,
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
  
  
  # Check if target_cell_types has cells not found in the data
  incorrect_cell_types <- setdiff(target_cell_types, unique(data[[feature_colname]]))
  if (length(incorrect_cell_types) > 0) {
    stop(paste(paste(incorrect_cell_types, collapse = ', '),
               "in target_cell_types don't exist in data."))
  }
  
  
  # Assume log_base is the length of target_cell_types
  # This ensures that entropy calculated is between 0 and 1, allowing for comparison
  if (is.null(log_base)) {
    log_base <- length(target_cell_types)
  }
  
  
  # Calculate entropy of the entire image
  if (is.null(radius) && is.null(reference_cell_type)) {
    
    entropy <- 0
    
    ## Get data for chosen target cells
    data <- data[(data[, feature_colname] %in% target_cell_types), ]
    
    n_all_cell_types <- nrow(data)
    
    ## No cells found or only one cell type present, return 0
    if (n_all_cell_types == 0 || length(target_cell_types) == 1) {
      return (0)
    }
    
    for (target_cell_type in target_cell_types) {
      
      ## Get data for current target cell 
      target_cell_type_data <- data[data[, feature_colname] == target_cell_type, ]
      n_target_cell_type <- nrow(target_cell_type_data)
      
      ## No cells found for current target cell, move on
      if (n_target_cell_type == 0) {
        next
      }
      
      target_cell_proportion <- n_target_cell_type / n_all_cell_types
      entropy <- entropy + (-1 * (target_cell_proportion) * log(target_cell_proportion, log_base))
    }
    
    return (entropy)
  }
  
  else if (is.null(radius) || is.null(reference_cell_type)) {
    stop("one of radius and reference_cell_type is NULL. 
         Both must be NULL to calculate entropy of whole image or 
         both must be specified to calculate entropy for each reference cell")
  }
  
  ## Radius has been specified, calculate entropy for chosen reference cell
  
  # Check if reference_cell_type is in the data
  if (!reference_cell_type %in% unique(data[[feature_colname]])) {
    stop(paste(reference_cell_type, " reference_cell_type does not exist in data"))
  }
  
  # Check if radius is numeric
  if (!is.numeric(radius)) {
    stop(paste(radius, " is not of type 'numeric'"))
  }
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighborhood_data <- calculate_cells_in_neighborhood3D(data,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radius,
                                                                  feature_colname)[[1]]
  
  ## Get total number of target cells for each row
  cells_in_neighborhood_data$Total <- apply(cells_in_neighborhood_data[target_cell_types], 1, sum)
  
  ## Get entropy for each row
  cells_in_neighborhood_data$Entropy <- 0
  
  for (target_cell_type in target_cell_types) {
    
    target_cell_type_proportions <- (cells_in_neighborhood_data[target_cell_type] / cells_in_neighborhood_data$Total)[[1]]
    
    ## If an element in target_cell_type_proportion is 0, just add 0.    
    cells_in_neighborhood_data$Entropy <- cells_in_neighborhood_data$Entropy +
      (-1 * target_cell_type_proportions * ifelse(target_cell_type_proportions == 0,
                                                  0, log(target_cell_type_proportions,
                                                         log_base)))
    
  }
  
  ## Case when row has 0 target cells
  cells_in_neighborhood_data[cells_in_neighborhood_data$Total == 0, "Entropy"] <- 0
  
  return (cells_in_neighborhood_data)
}



calculate_entropy_gradient3D <- function(data,
                                         radii,
                                         reference_cell_type,
                                         target_cell_types,
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
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
  
  # Check if radii is numeric
  if (!is.numeric(radii)) {
    stop(paste(radii, " radii is not of type 'numeric'"))
  }
  
  # Check if reference_cell_type is in the data
  if (!reference_cell_type %in% unique(data[[feature_colname]])) {
    stop(paste(reference_cell_type, " reference_cell_type does not exist in data"))
  }
  
  # Check if target_cell_types has cells not found in the data
  incorrect_cell_types <- setdiff(target_cell_types, unique(data[[feature_colname]]))
  if (length(incorrect_cell_types) > 0) {
    stop(paste(paste(incorrect_cell_types, collapse = ', '),
               "in target_cell_types don't exist in data."))
  }
  
  
  entropy_gradient <- list()
  entropy_mean <- c()
  
  for (radius in seq(radii)) {
    entropy_data <- calculate_entropy3D(data,
                                        radius,
                                        reference_cell_type,
                                        target_cell_types,
                                        length(target_cell_types),
                                        feature_colname)
    
    entropy_gradient[[paste(radius)]] <- entropy_data
    entropy_mean <- append(entropy_mean, mean(entropy_data$Entropy))
    
  }
  
  if (plot_image) {
    plot(seq(radii), entropy_mean, type = "l", xlab = "Radius", ylab = "Entropy Mean")
  }
  
  return (entropy_gradient)
  
}



calculate_entropy_gradient_aggregated3D <- function(data,
                                                    radii,
                                                    reference_cell_type,
                                                    target_cell_types,
                                                    feature_colname = "Cell.Type",
                                                    plot_image = TRUE) {
  
  
  ## Get entropy gradient data
  ## Entropy column is useless but the total column is nice.
  entropy_gradient_data <- calculate_entropy_gradient3D(data,
                                                        radii,
                                                        reference_cell_type,
                                                        target_cell_types,
                                                        feature_colname,
                                                        FALSE)
  
  
  result <- data.frame()
  
  for (i in seq(length(entropy_gradient_data))) {
    
    ## Subset each data frame from the entropy_gradient_data list so it only
    ## includes the cells and the total columns
    
    # reference cell type might also be in target cell type, no need to double up
    all_cell_types <- unique(c(reference_cell_type, target_cell_types))
    
    entropy_gradient_data[[i]] <- entropy_gradient_data[[i]][c(all_cell_types, "Total")]
    
    ## Add the summed values of each column to result data frame
    result <- rbind(result, t(apply(entropy_gradient_data[[i]], 2, sum)))
    
  }
  rownames(result) <- names(entropy_gradient_data)
  
  ## Get entropies for each element in the data frame
  result_entropies <- result / result$Total # Make cell count values into cell proportion values
  result_entropies <- -1 * result_entropies  * log(result_entropies, length(target_cell_types))
  result_entropies <- apply(result_entropies, 2, function(x) replace(x, is.nan(x), 0))
  
  ## Calculate total entropy for each row in result data frame
  result$Entropy <- apply(result_entropies, 1, sum)
  
  # Plot
  if (plot_image) {
    plot(rownames(result), result$Entropy, type = "l", xlab = "Radius", ylab = "Entropy Aggregated")
  }
  
  return (result)
  
}






determine_entropy_grid_metrics3D <- function(data, 
                                             n_split,
                                             target_cell_types,
                                             feature_colname = "Cell.Type",
                                             size = NULL,
                                             plot_image = TRUE) {
  
  # If the columns are not correct, give error
  required_colnames <- c("Cell.X.Position", 
                         "Cell.Y.Position", 
                         "Cell.Z.Position", 
                         feature_colname)
  
  missing_colnames <- setdiff(required_colnames,
                              colnames(data))
  
  if (length(missing_colnames) > 0) {
    stop(paste(paste(missing_colnames, collapse = ', '),
               "are missing as column names in your data")) 
  }
  
  # Check if n_split is numeric
  if (!is.numeric(n_split)) {
    stop(paste(n_split, " n_split is not of type 'numeric'"))
  }
  
  # Check if target_cell_types has cells not found in the data
  incorrect_cell_types <- setdiff(target_cell_types, unique(data[[feature_colname]]))
  if (length(incorrect_cell_types) > 0) {
    stop(paste(paste(incorrect_cell_types, collapse = ', '),
               "in target_cell_types don't exist in data."))
  }
  
  
  ## Get dimensions of the window
  length <- round(max(data$Cell.X.Position) - min(data$Cell.X.Position))
  width  <- round(max(data$Cell.Y.Position) - min(data$Cell.Y.Position))
  height <- round(max(data$Cell.Z.Position) - min(data$Cell.Z.Position))
  
  ## Get distance of row, col and lay
  d_row <- length / n_split
  d_col <- width / n_split
  d_lay <- height / n_split
  
  ## Figure out which 'grid prism number' each cell is inside
  data$Prism.Num <- floor(data$Cell.X.Position / d_row) +
    floor(data$Cell.Y.Position / d_col) * n_split + 
    floor(data$Cell.Z.Position / d_lay) * n_split^2 + 1
  
  
  ## Calculate entropy for each grid prism
  n_grid_prisms <- n_split^3
  cell_type_list <- vector(mode = 'list', length = length(target_cell_types))
  grid_prism_entropies <- c()
  
  for (grid_prism_num in seq(n_grid_prisms)) {
    
    ## Get data of cells in the current grid_prism
    data_temp <- data[data$Prism.Num == grid_prism_num, ]
    
    temp_target_cell_types <- intersect(target_cell_types, unique(data_temp[[feature_colname]]))
    
    grid_prism_entropy <- calculate_entropy3D(data_temp,
                                              target_cell_types = temp_target_cell_types,
                                              log_base = length(target_cell_types))
    
    ## Get number of cells of each target cell type in each grid prism
    for (target_cell_type in target_cell_types) {
      cell_type_list[[target_cell_type]] <- append(cell_type_list[[target_cell_type]], 
                                                   sum(data_temp$Cell.Type == target_cell_type))
    }
    
    grid_prism_entropies <- c(grid_prism_entropies, grid_prism_entropy)
    
  }
  
  result <- data.frame(row.names = seq(n_grid_prisms))
  
  ## Add column for each target cell type representing the number of cells in each grid prism
  for (target_cell_type in target_cell_types) {
    result[[target_cell_type]] <- cell_type_list[[target_cell_type]]
  }
  
  ## Add column for total cell count for each grid prism
  result$Total <- apply(result, 1, sum)
  
  ## Add entropy column
  result$Entropy = grid_prism_entropies
  
  ## Plot
  if (plot_image) {
    
    # Check if size is numeric or not
    if (!is.numeric(size)) {
      stop(paste(size, " size is not numeric"))
    }
    
    plot_data <- result
    
    ## Place a dot at the center of each grid prism to represent entropy
    ## Use the grid prism number to figure out their location...
    plot_data$x <- ((seq(n_grid_prisms) - 1) %% n_split + 0.5) * d_row
    plot_data$y <- (floor(((seq(n_grid_prisms) - 1) %% (n_split)^2) / n_split) + 0.5) * d_col
    plot_data$z <- (floor((seq(n_grid_prisms) - 1) / (n_split^2)) + 0.5) * d_lay
    
    ## Color of each dot is related to its entropy
    pal <- colorRampPalette(hcl.colors(n = 5, palette = "Terrain", rev = TRUE))
    
    ## Add size column and for 0 cell proportion values, make the size small
    plot_data$size <- ifelse(plot_data$Entropy == 0, 5, size)
    
    fig <- plot_ly(plot_data,
                   type = "scatter3d",
                   mode = 'markers',
                   x = ~x,
                   y = ~y,
                   z = ~z,
                   color = ~Entropy,
                   colors = pal(nrow(plot_data)),
                   marker = list(size = ~size))
    
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                       yaxis = list(title = 'y'),
                                       zaxis = list(title = 'z')))
    
    print(fig)
    
  }
  
  return (result)
}



determine_cell_proportion_grid_metrics3D <- function(data, 
                                                     n_split,
                                                     reference_cell_types,
                                                     target_cell_types,
                                                     feature_colname = "Cell.Type",
                                                     size = NULL,
                                                     plot_image = TRUE) {
  
  # If the columns are not correct, give error
  required_colnames <- c("Cell.X.Position", 
                         "Cell.Y.Position", 
                         "Cell.Z.Position", 
                         feature_colname)
  
  missing_colnames <- setdiff(required_colnames,
                              colnames(data))
  
  if (length(missing_colnames) > 0) {
    stop(paste(paste(missing_colnames, collapse = ', '),
               "are missing as column names in your data")) 
  }
  
  # Check if n_split is numeric
  if (!is.numeric(n_split)) {
    stop(paste(n_split, " n_split is not of type 'numeric'"))
  }
  
  # Check if reference_cell_types has cells not found in the data
  incorrect_cell_types <- setdiff(reference_cell_types, unique(data[[feature_colname]]))
  if (length(incorrect_cell_types) > 0) {
    stop(paste(paste(incorrect_cell_types, collapse = ', '),
               "in reference_cell_types don't existin data."))
  }
  # Check if target_cell_types has cells not found in the data
  incorrect_cell_types <- setdiff(target_cell_types, unique(data[[feature_colname]]))
  if (length(incorrect_cell_types) > 0) {
    stop(paste(paste(incorrect_cell_types, collapse = ', '),
               "in target_cell_types don't exist in data."))
  }
  # Check if there is intersection between reference_cell_types and target_cell_types
  if (length(intersect(reference_cell_types, target_cell_types)) > 0) {
    stop("Cannot have same cells in both reference_cell_types and target_cell_types")
  }
  
  
  
  
  ## Get dimensions of the window
  length <- round(max(data$Cell.X.Position) - min(data$Cell.X.Position))
  width  <- round(max(data$Cell.Y.Position) - min(data$Cell.Y.Position))
  height <- round(max(data$Cell.Z.Position) - min(data$Cell.Z.Position))
  
  ## Get distance of row, col and lay
  d_row <- length / n_split
  d_col <- width / n_split
  d_lay <- height / n_split
  
  ## Figure out which 'grid prism number' each cell is inside
  data$Prism.Num <- floor(data$Cell.X.Position / d_row) +
    floor(data$Cell.Y.Position / d_col) * n_split + 
    floor(data$Cell.Z.Position / d_lay) * n_split^2 + 1
  
  ## Calculate cell_proportions for each grid prism
  n_grid_prisms <- n_split^3
  n_reference_cells_vec <- c()
  n_target_cells_vec <- c()
  grid_prism_cell_proportions <- c()
  
  for (grid_prism_num in seq(n_grid_prisms)) {
    
    ## Get data of cells in the current grid_prism
    data_temp <- data[data$Prism.Num == grid_prism_num, ]
    
    ## Get cell_proportion: n_target_cells / (n_target_cells + n_reference_cells)
    n_target_cells <- sum(data_temp[[feature_colname]] %in% target_cell_types)
    n_reference_cells <- sum(data_temp[[feature_colname]] %in% reference_cell_types)
    
    ## Case when there are no target or reference cell, result is NA
    if (n_target_cells == 0 && n_reference_cells == 0) {
      grid_prism_cell_proportion <- NA  
    }
    else {
      grid_prism_cell_proportion <- n_target_cells / (n_target_cells + n_reference_cells)
    }
    
    n_reference_cells_vec <- c(n_reference_cells_vec, n_reference_cells)
    n_target_cells_vec <- c(n_target_cells_vec, n_target_cells)
    
    grid_prism_cell_proportions <- c(grid_prism_cell_proportions, grid_prism_cell_proportion)
    
  }
  
  result <- data.frame(row.names = seq(n_grid_prisms))
  
  ## Add column for reference cell type and target cell type representing the number of cells in each grid prism
  result[["Reference"]] <- n_reference_cells_vec
  result[["Target"]] <- n_target_cells_vec
  
  ## Add column for total cell count for each grid prism
  result$Total <- apply(result, 1, sum)
  
  ## Add cell proportion column
  result$Proportion = grid_prism_cell_proportions
  
  ## Plot
  if (plot_image) {
    
    # Check if size is numeric or not
    if (!is.numeric(size)) {
      stop(paste(size, " size is not numeric"))
    }
    
    plot_data <- result
    
    ## Place a dot at the center of each grid prism to represent cell proportion
    ## Use the grid prism number to figure out their location...
    plot_data$x <- ((seq(n_grid_prisms) - 1) %% n_split + 0.5) * d_row
    plot_data$y <- (floor(((seq(n_grid_prisms) - 1) %% (n_split)^2) / n_split) + 0.5) * d_col
    plot_data$z <- (floor((seq(n_grid_prisms) - 1) / (n_split^2)) + 0.5) * d_lay
    
    ## Color of each dot is related to its cell proportion
    pal <- colorRampPalette(hcl.colors(n = 5, palette = "Red-Blue", rev = TRUE))
    
    
    ## Add size column and for NA cell proportion values, make the size small
    plot_data$size <- ifelse(is.na(plot_data$Proportion), 3, size)
    
    fig <- plot_ly(plot_data,
                   type = "scatter3d",
                   mode = 'markers',
                   x = ~x,
                   y = ~y,
                   z = ~z,
                   color = ~Proportion,
                   colors = pal(nrow(plot_data)),
                   marker = list(size = ~size))
    
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                       yaxis = list(title = 'y'),
                                       zaxis = list(title = 'z')))
    
    print(fig)
    
  }
  
  return (result)
}




determine_prevalence3D <- function(grid_data,
                                   metric_colname = "Entropy",
                                   threshold,
                                   above = TRUE) {
  
  ## Exclude rows with NA values
  grid_data <- grid_data[!is.na(grid_data[[metric_colname]]), ]
  
  if (above) {
    p <- sum(grid_data[[metric_colname]] >= threshold) / nrow(grid_data) * 100
  }
  else {
    p <- sum(grid_data[[metric_colname]] < threshold) / nrow(grid_data) * 100    
  }
  
  return (p)
}




determine_spatial_autocorrelation <- function(grid_data,
                                              metric_colname = "Entropy",
                                              weight_method = "IDW") {
  
  
  ## Get number of grid prisms
  n_grid_prisms <- nrow(grid_data)
  
  ## Get splitting number (should be the cube root of n_grid_prisms)
  n_split <- (n_grid_prisms)^(1/3)
  
  ## Find the coordinates of each grid prism
  x <- ((seq(n_grid_prisms) - 1) %% n_split)
  y <- (floor(((seq(n_grid_prisms) - 1) %% (n_split)^2) / n_split))
  z <- (floor((seq(n_grid_prisms) - 1) / (n_split^2)))
  grid_prism_coords <- data.frame(x = x, y = y, z = z)
  
  
  weight_matrix <- -1 * apcluster::negDistMat(grid_prism_coords)
  ## Use the inverse distance between two points as the weight (IDW is 'inverse distance weighting')
  if (weight_method == "IDW") {
    weight_matrix <- 1 / weight_matrix
  }
  ## Use binary method: adjacent points get a weight of 1, otherwise, weight of 0
  ## Adjacent points are within sqrt(3) units apart. e.g. (0, 0, 0) vs (1, 1, 1)
  else if (weight_method == "Binary") {
    weight_matrix <- ifelse(weight_matrix > sqrt(3), 0, 1)  
  }
  else {
    stop(paste(weight_method, " weight_method is not an appropriate method"))
  }
  
  ## Points along the diagonal are comparing the same point so its weight is zero
  diag(weight_matrix) <- 0
  
  data_mean <- mean(grid_data[!is.na(grid_data[[metric_colname]]), metric_colname])
  
  numerator <- 0
  denominator <- 0
  
  for (i in seq(n_grid_prisms)) {
    
    if (is.na(grid_data[i, metric_colname])) {
      next
    }
    
    for (j in seq(n_grid_prisms)) {
      
      if (is.na(grid_data[j, metric_colname])) {
        next
      }
      
      numerator <- numerator + weight_matrix[i, j] * 
        (grid_data[i, metric_colname] - data_mean) * 
        (grid_data[j, metric_colname] - data_mean)
      
    }
    denominator <- denominator + (grid_data[i, metric_colname] - data_mean)^2
  }
  
  
  I <- (n_grid_prisms * numerator) / (sum(weight_matrix) * denominator)
  
  return (I)
  
}
