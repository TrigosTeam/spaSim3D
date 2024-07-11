### Basic metrics -------------------------------------------------------------

calculate_cell_proportions3D <- function(spe,
                                         cell_types_of_interest = NULL, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), "Cell.Type" = spe[[feature_colname]])
  
  # Check
  if (nrow(df) == 0) stop("No cells found for calculating cell proportions")
  
  # Creates frequency/bar plot of all cell types in the entire image
  cell_proportions <- data.frame(table(df[, feature_colname]))
  names(cell_proportions) <- c("cell_type", 'frequency')
  
  # Only include cell types the user has chosen
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, cell_proportions$cell_type)
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                 paste(unknown_cell_types, collapse = ", ")))
    }
    
    # Subset for cell types chosen by user
    cell_proportions <- cell_proportions[(cell_proportions$cell_type %in% cell_types_of_interest), ]
    
    # Check if the user has excluded all cell types
    if (nrow(cell_proportions) == 0) stop("All cells have been excluded")
  }
  
  # Get frequency total for all cells
  cell_type_frequency_total <- sum(cell_proportions$frequency)
  
  # Get proportions and percentages
  cell_proportions$proportion <- cell_proportions$frequency / cell_type_frequency_total
  cell_proportions$percentage <- cell_proportions$proportion * 100
  
  # Order the cell types by proportion (highest cell proportion is first)
  cell_proportions <- cell_proportions[rev(order(cell_proportions$proportion)), ]
  rownames(cell_proportions) <- seq(nrow(cell_proportions))
  
  
  # Plot
  if (plot_image) {
    
    labels <- paste(round(cell_proportions$percentage, 1), "%", sep = "")
    
    fig <- ggplot(cell_proportions, aes(x = factor(cell_type, cell_type), y = percentage, fill = cell_type)) +
      geom_bar(stat='identity') + 
      theme_bw() +
      labs(title="Cell proportions", x = "Cell type", y = "Percentage") +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "none") +
      geom_text(aes(label = labels), vjust = 0)
    
    methods::show(fig)
  }
  
  return(cell_proportions)
}


calculate_entropy_background3D <- function(spe,
                                           cell_types_of_interest, 
                                           feature_colname = "Cell.Type") {
  
  if (length(cell_types_of_interest) == 0) return(NA)
  if (length(cell_types_of_interest) == 1) return(0)
  
  cell_proportions_data <- calculate_cell_proportions3D(spe, cell_types_of_interest, feature_colname, FALSE)
  
  # Calculate entropy of the entire image
  entropy <- -1 * sum(cell_proportions_data$proportion * log(cell_proportions_data$proportion, length(cell_proportions_data$proportion)))
  
  return(entropy) 
}


### Cell colocalisation metrics -----------------------------------------------
calculate_pairwise_distances_between_cell_types3D <- function(spe,
                                                              cell_types_of_interest = NULL,
                                                              feature_colname = "Cell.Type",
                                                              show_summary = TRUE,
                                                              plot_image = TRUE) {
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), 
                   "Cell.Type" = spe[[feature_colname]], 
                   "Cell.ID" = spe[["Cell.ID"]])
  
  # If there are no cells, give error
  if (nrow(df) == 0) {
    stop("There are no cells in data")
  }
  
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
  for (eachType in unique(df[ , feature_colname])) {
    cell_types[[eachType]] <- as.character(df$Cell.ID[df[, feature_colname] == eachType])
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


## Please ensure there is no factoring in any of the columns!!!

calculate_minimum_distances_between_cell_types3D <- function(spe,
                                                             cell_types_of_interest = NULL,
                                                             feature_colname = "Cell.Type",
                                                             show_summary = TRUE,
                                                             plot_image = TRUE) {
  
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
  for (eachType in unique(df[["Cell.Type"]])) {
    cell_types[[eachType]] <- as.character(df$Cell.ID[df[["Cell.Type"]] == eachType])
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


summarise_distances_between_cell_types3D <- function(df) {
  
  pair <- distance <- NULL
  
  # summarise the results
  summarised_dists <- df %>% 
    dplyr::group_by(pair) %>%
    dplyr::summarise(mean(distance, na.rm = TRUE), 
                     min(distance, na.rm = TRUE), 
                     max(distance, na.rm = TRUE),
                     stats::median(distance, na.rm = TRUE), 
                     stats::sd(distance, na.rm = TRUE))
  
  summarised_dists <- data.frame(summarised_dists)
  
  colnames(summarised_dists) <- c("pair", 
                                  "mean", 
                                  "min", 
                                  "max", 
                                  "median", 
                                  "std_dev")
  
  for (i in seq(nrow(summarised_dists))) {
    # Get cell_types for each pair
    cell_types <- strsplit(summarised_dists[i,"pair"], "/")[[1]]
    
    summarised_dists[i, "reference"] <- cell_types[1]
    summarised_dists[i, "target"] <- cell_types[2]
  }
  
  return(summarised_dists)
}


## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_cell_distances_violin3D <- function(cell_to_cell_dist, scales = "free_x") {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  pair <- distance <- NULL
  
  fig <- ggplot(cell_to_cell_dist, aes(x = pair, y = distance)) + 
    geom_violin() +
    facet_wrap(~pair, scales=scales, strip.position="bottom") +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
    labs(title="Cell distances", x = "Reference/Target pair", y = "Distance") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plots show mean ± sd")
  
  return(fig)
}




calculate_mixing_scores3D <- function(spe, 
                                      reference_cell_types, 
                                      target_cell_types, 
                                      radius, 
                                      feature_colname = "Cell.Type") {
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), "Cell.Type" = spe[[feature_colname]])
  
  ## For reference_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(reference_cell_types, df$Cell.Type)
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in reference_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  ## For target_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, df$Cell.Type)
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  # Check if radius is numeric
  if (!is.numeric(radius)) stop(paste(radius, " is not of type 'numeric'"))
  
  
  result <- data.frame(matrix(ncol=8, nrow=0))
  
  for (reference_cell_type in reference_cell_types) {
    
    # Get all info for cells of reference cell_type
    reference_cells <- df[df[[feature_colname]] == reference_cell_type, ]
    
    for (target_cell_type in target_cell_types) {
      
      # Get all info for cells of target cell_type      
      target_cells <- df[df[[feature_colname]] == target_cell_type, ]
      
      # No point getting mixing scores if comparing the same cell type
      if (reference_cell_type == target_cell_type) {
        next
      }
      
      # Can't get mixing scores if there are no reference cells
      if (nrow(reference_cells) == 0) {
        methods::show(paste("There are no unique reference cells of specified cell type ", reference_cell_type, "for target cell", target_cell_type))
        result <-  rbind(result, 
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
        methods::show(paste("There are no unique target cells of specified cell type", target_cell_type, "for reference cell", reference_cell_type))
        
        reference_cell_coords <- reference_cells[, c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
        reference_reference_result <- dbscan::frNN(reference_cell_coords, 
                                                   eps = radius, 
                                                   query = NULL,
                                                   sort = FALSE)
        
        # halve it to avoid counting each ref-ref interaction twice
        reference_reference_interactions <- 0.5 * sum(rapply(reference_reference_result$id, length)) 
        
        result <-  rbind(result, 
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
        
        result <-  rbind(result, 
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
  colnames(result) <- c("ref_cell_type", 
                        "tar_cell_type", 
                        "n_ref_cells",
                        "n_tar_cells", 
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions", 
                        "mixing_score", 
                        "normalised_mixing_score")
  
  # Turn numeric data into numeric type
  result[ , 3:8] <- apply(result[ , 3:8], 2, as.numeric)
  
  return(result)
}


calculate_cells_in_neighbourhood3D <- function(spe, 
                                              reference_cell_type, 
                                              target_cell_types, 
                                              radius, 
                                              feature_colname = "Cell.Type",
                                              show_summary = TRUE,
                                              plot_image = TRUE) {
  
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


calculate_cells_in_neighbourhood_gradient3D <- function(spe, 
                                                        reference_cell_type, 
                                                        target_cell_types, 
                                                        radii, 
                                                        feature_colname = "Cell.Type",
                                                        plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (radius in seq(radii)) {
    cells_in_neighborhood_data <- calculate_cells_in_neighbourhood3D(spe,
                                                                     reference_cell_type,
                                                                     target_cell_types,
                                                                     radius,
                                                                     feature_colname,
                                                                     FALSE,
                                                                     FALSE)
    
    cells_in_neighborhood_data$ref_cell_id <- NULL
    result[radius, ] <- apply(cells_in_neighborhood_data, 2, mean)
  }
  # Add a radius column to the result
  result$radius <- seq(radii)
  
  if (plot_image) {
    plot_result <- reshape2::melt(result, "radius")
    
    fig <- ggplot(plot_result, aes(radius, value, color = variable)) + 
      geom_line() + 
      labs(x = "Radius", y = "Average cells in neighbourhood") + 
      scale_color_discrete(name = "Cell type")
    
    methods::show(fig)
  }
  
  return(result)
}



summarise_cells_in_neighbourhood3D <- function(cells_in_neighbourhood_data) {
  
  df <- cells_in_neighbourhood_data
  
  ## Target cell types will be all the columns except the first column
  target_cell_types <- colnames(cells_in_neighbourhood_data)[c(-1)]
  
  ## Set up data frame for summarised_results list
  df <- data.frame(row.names = c("mean", "min", "max", "median", "st_dev"))
  
  for (target_cell_type in target_cell_types) {
    
    ## Get statistical measures for each target cell type
    target_cell_type_values <- cells_in_neighbourhood_data[[target_cell_type]]
    df[[target_cell_type]] <- c(mean(target_cell_type_values),
                                min(target_cell_type_values),
                                max(target_cell_type_values),
                                median(target_cell_type_values),
                                sd(target_cell_type_values))
    
  }
  
  return (data.frame(t(df)))
}


calculate_cells_in_neighbourhood_proportions3D <- function(spe, 
                                                           reference_cell_type, 
                                                           target_cell_types, 
                                                           radius, 
                                                           feature_colname = "Cell.Type") {
  
  ## Get 'count' neighbourhood data
  cells_in_neighbourhood_data <- calculate_cells_in_neighbourhood3D(spe,
                                                                    reference_cell_type,
                                                                    target_cell_types,
                                                                    radius,
                                                                    feature_colname,
                                                                    FALSE,
                                                                    FALSE)
  
  
  
  result <- data.frame(matrix(nrow = length(target_cell_types), ncol = 4))
  colnames(result) <- c("target_cell_type", "frequency", "proportion", "percentage")
  
  result$target_cell_type <- target_cell_types
  
  ## Get frequency of each target cell type
  result$frequency <- apply(cells_in_neighbourhood_data[ , target_cell_types], 2, sum)
  
  ## Use frequency to get proportion and percentage of each cell type
  total <- sum(result$frequency)
  if (total != 0) {
    result$proportion <- result$frequency / total
    result$percentage <- result$proportion * 100  
  }
  else {
    result$proportion <- NA
    result$percentage <- NA
  }
  
  
  return(result)
}


calculate_cells_in_neighbourhood_proportions_gradient3D <- function(spe, 
                                                                    reference_cell_type, 
                                                                    target_cell_types, 
                                                                    radii, 
                                                                    feature_colname = "Cell.Type",
                                                                    plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (radius in seq(radii)) {
    cell_proportions_neighbourhood_data <- calculate_cells_in_neighbourhood_proportions3D(spe,
                                                                                          reference_cell_type,
                                                                                          target_cell_types,
                                                                                          radius,
                                                                                          feature_colname)
    
    result[radius, ] <- cell_proportions_neighbourhood_data$proportion
  }
  
  # Add a radius column to the result
  result$radius <- seq(radii)
  
  # Plot
  if (plot_image) {
    plot_result <- reshape2::melt(result, id.vars = c("radius"))
    fig <- ggplot(plot_result, aes(radius, value, color = variable)) +
      geom_point() +
      geom_line() +
      labs(title = "Neighbourhood cell proportion gradients", x = "Radius", y = "Cell proportion", color = "Cell type") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(0, 1)
    
    methods::show(fig)
  }
  
  return(result)
}





## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_cells_in_neighbourhood_violin3D <- function(cells_in_neighbourhood_data, reference_cell_type, scales = "free_x") {
  
  ## Target cell types will be all the columns except the first column
  target_cell_types <- colnames(cells_in_neighbourhood_data)[c(-1)]
  
  df <- reshape2::melt(cells_in_neighbourhood_data, measure.vars = target_cell_types)
  colnames(df) <- c("ref_cell_id", "tar_cell_type", "count")
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  tar_cell_type <- count <- NULL
  
  fig <- ggplot(df, aes(x = tar_cell_type, y = count)) + 
    geom_violin() +
    facet_wrap(~tar_cell_type, scales=scales, strip.position="bottom") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(title=paste("Cells in neighbourhood of", reference_cell_type, "cells"), x = "Target cell type", y = "Number of cells") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plots show mean ± sd")
  
  return(fig)
}




calculate_entropy3D <- function(spe,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname = "Cell.Type",
                                plot_image = TRUE) {
  
  # Check if radius is numeric
  if (!is.numeric(radius)) {
    stop(paste(radius, " is not of type 'numeric'"))
  }
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighbourhood_data <- calculate_cells_in_neighbourhood3D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  FALSE,
                                                                  FALSE)
  
  ## Get total number of target cells for each row
  cells_in_neighbourhood_data$total <- apply(cells_in_neighbourhood_data[ , c(-1)], 1, sum)
  
  ## Get entropy for each row
  cells_in_neighbourhood_data$entropy <- 0
  
  for (target_cell_type in target_cell_types) {
    
    target_cell_type_proportions <- (cells_in_neighbourhood_data[[target_cell_type]] / cells_in_neighbourhood_data$total)
    
    ## If an element in target_cell_type_proportion is 0, just add 0.    
    target_cell_entropy <- ifelse(target_cell_type_proportions == 0,
                                  0,
                                  -1 * target_cell_type_proportions * log(target_cell_type_proportions, length(target_cell_types)))
    
    cells_in_neighbourhood_data$entropy <- cells_in_neighbourhood_data$entropy + target_cell_entropy
    
  }
  
  ## Case when row has 0 target cells
  cells_in_neighbourhood_data[cells_in_neighbourhood_data$Total == 0, "entropy"] <- 0
  
  if (plot_image) {
    fig <- plot_entropy_violin3D(cells_in_neighbourhood_data)
    methods::show(fig)
  }
  
  return(cells_in_neighbourhood_data)
}


## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_entropy_violin3D <- function(entropy_data, scales = "free_x") {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  entropy <- NULL
  
  fig <- ggplot(entropy_data, aes(x = "", y = entropy)) +
    geom_violin() +
    theme_bw() +
    labs(x = "", y = "Entropy") +
    theme(axis.ticks.x = element_blank()) +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plot shows mean ± sd")
  
  return(fig)
}







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
  
  
  cells_in_neighbourhood_data <- calculate_cells_in_neighbourhood3D(spe,
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


calculate_mixing_scores_gradient3D <- function(spe, 
                                               reference_cell_type, 
                                               target_cell_type, 
                                               radii, 
                                               feature_colname = "Cell.Type",
                                               plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = 8))
  colnames(result) <- c("ref_cell_type", 
                        "tar_cell_type", 
                        "n_ref_cells",
                        "n_tar_cells", 
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions", 
                        "mixing_score", 
                        "normalised_mixing_score")
  
  for (radius in seq(radii)) {
    mixing_scores <- calculate_mixing_scores3D(spe,
                                               reference_cell_type,
                                               target_cell_type,
                                               radius,
                                               feature_colname)
    
    result[radius, ] <- mixing_scores
  }
  
  # Add a radius column to the result
  result$radius <- seq(radii)
  
  if (plot_image) {
    plot(result[["radius"]], result[["normalised_mixing_score"]], 
         type = "l", 
         xlab = "Radius", 
         ylab = "normalised mixing score",
         ylim = c(0, max(result[["normalised_mixing_score"]], 1)),
         col = "red")
    abline(a = 1, b = 0, col = "blue", lty = 2)
    legend(0, 0.95, legend = c("Observed normalised mixing score", "Expected CSR normalised mixing score"), col = c("red", "blue"), lty = c(1, 2))
  }
  
  return(result)
}





calculate_cross_K_gradient3D <- function(spe, 
                                         reference_cell_type, 
                                         target_cell_type, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = 2))
  colnames(result) <- c("observed_cross_K", 
                        "expected_cross_K")
  
  for (radius in seq(radii)) {
    cross_K_data <- calculate_cross_K3D(spe,
                                        reference_cell_type,
                                        target_cell_type,
                                        radius,
                                        feature_colname)
    
    result[radius, ] <- cross_K_data
  }
  
  # Add a radius column to the result
  result$radius <- seq(radii)
  
  if (plot_image) {
    plot(result$radius, result$observed_cross_K, type = "l", col = "red", 
         xlim = c(0, radius), ylim = c(0, max(result)),
         xlab = "Radius", ylab = "Cross K-function value")
    lines(result$radius, result$expected_cross_K, type = "l", col = "blue", lty = 2)
    legend(0, max(result), legend = c("Observed cross K", "Expected CSR cross K"), col = c("red", "blue"), lty = c(1, 2))
  }
  
  return(result)
}


calculate_entropy_gradient3D <- function(spe,
                                         reference_cell_type,
                                         target_cell_types,
                                         radii,
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (radius in seq(radii)) {
    cells_in_neighbourhood_data <- calculate_cells_in_neighbourhood3D(spe,
                                                                    reference_cell_type,
                                                                    target_cell_types,
                                                                    radius,
                                                                    feature_colname,
                                                                    FALSE,
                                                                    FALSE)
    
    cells_in_neighbourhood_data$ref_cell_id <- NULL
    result[radius, ] <- apply(cells_in_neighbourhood_data, 2, sum)
  }
  
  ## Get total number of target cells for each row
  result$total <- apply(result, 1, sum)
  
  ## Set intial entropy to 0
  result$entropy <- 0
  
  for (target_cell_type in target_cell_types) {
    
    target_cell_type_proportions <- (result[[target_cell_type]] / result$total)
    
    ## If an element in target_cell_type_proportion is 0, just add 0.    
    target_cell_entropy <- ifelse(target_cell_type_proportions == 0,
                                  0,
                                  -1 * target_cell_type_proportions * log(target_cell_type_proportions, length(target_cell_types)))
    
    result$entropy <- result$entropy + target_cell_entropy
    
  }
  
  # Add a radius column to the result
  result$radius <- seq(radii)
  
  if (plot_image) {
    expected_entropy <- calculate_entropy_background3D(spe, target_cell_types, feature_colname)
    
    plot(result$radius, result$entropy, type = "l", col = "red", 
         xlim = c(0, radius), ylim = c(0, max(result$entropy)),
         xlab = "Radius", ylab = "Entropy")
    abline(a = expected_entropy, b = 0, col = "blue", lty = 2)
    legend(0, max(result$entropy, expected_entropy), legend = c("Observed entropy", "Expected CSR entropy"), col = c("red", "blue"), lty = c(1, 2))
  }
  
  return(result)
}


plot_cross_K_gradient_ratio3D <- function(cross_K_gradient_results) {
  
  plot(cross_K_gradient$radius, 
       cross_K_gradient$observed_cross_K / cross_K_gradient$expected_cross_K, 
       type = "l", 
       col = "red", 
       xlim = c(0, max(cross_K_gradient$radius)), ylim = c(0, 1.2 * max((cross_K_gradient$observed_cross_K / cross_K_gradient$expected_cross_K), 1)),
       xlab = "Radius", ylab = "Cross K-function ratio")
  abline(a = 1, b = 0, col = "blue", lty = 2)
  legend(0, 1.2 * max((cross_K_gradient$observed_cross_K / cross_K_gradient$expected_cross_K), 1), 
         legend = c("Observed cross K ratio", "Expected CSR cross K ratio"), col = c("red", "blue"), lty = c(1, 2))
  
}



### Spatial heterogeneity metrics ---------------------------------------------
determine_entropy_grid_metrics3D <- function(spe, 
                                             n_splits,
                                             cell_types_of_interest,
                                             feature_colname = "Cell.Type",
                                             plot_image = TRUE) {
  
  
  # Check if n_splits is numeric
  if (!is.numeric(n_splits)) {
    stop(paste(n_splits, " n_splits is not of type 'numeric'"))
  }
  
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, unique(spe[[feature_colname]]))
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  
  spe_coords <- data.frame(spatialCoords(spe))
  
  ## Get dimensions of the window
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  d_lay <- height / n_splits
  
  ## Figure out which 'grid prism number' each cell is inside
  spe$Prism.Num <- floor(spe_coords$Cell.X.Position / d_row) +
    floor(spe_coords$Cell.Y.Position / d_col) * n_splits + 
    floor(spe_coords$Cell.Z.Position / d_lay) * n_splits^2 + 1
  
  ## Get number of grid prisms
  n_grid_prisms <- n_splits^3
  
  ## Define data frame which contains all results
  result <- data.frame(matrix(nrow = n_grid_prisms, ncol = (length(cell_types_of_interest) + 2)))
  colnames(result) <- c(cell_types_of_interest, "total", "entropy")
  
  ## Calculate entropy for each grid prism
  for (grid_prism_num in seq(n_grid_prisms)) {
    
    ## Get spe object for current grid_prism
    spe_temp <- spe[ , spe$Prism.Num == grid_prism_num]
    
    ## Get cell_types_of_interest found in the sub-spe object
    temp_cell_types_of_interest <- intersect(cell_types_of_interest, unique(spe_temp[[feature_colname]]))
    
    grid_prism_entropy <- calculate_entropy_background3D(spe_temp,
                                                         temp_cell_types_of_interest)
    result[grid_prism_num, "entropy"] <- grid_prism_entropy
    
    ## Get number of cells of each cell_types_of_interest in each grid prism
    for (cell_type_of_interest in cell_types_of_interest) {
      result[grid_prism_num, cell_type_of_interest] <- sum(spe_temp[[feature_colname]] == cell_type_of_interest)
    }
  }
  
  ## Add column for total cell count for each grid prism
  result$total <- apply(result[ , colnames(result) %in% cell_types_of_interest], 1, sum)
  
  ## Add x, y and z coords of each grid prism to result
  result$prism_x_coord <- ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row
  result$prism_y_coord <- (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col
  result$prism_z_coord <- (floor((seq(n_grid_prisms) - 1) / (n_splits^2)) + 0.5) * d_lay
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous3D(result, "entropy")
    methods::show(fig)
  }
  
  return(result)
}




determine_cell_proportion_grid_metrics3D <- function(spe, 
                                                     n_splits,
                                                     reference_cell_types,
                                                     target_cell_types,
                                                     feature_colname = "Cell.Type",
                                                     plot_image = TRUE) {
  
  
  
  # Check if n_splits is numeric
  if (!is.numeric(n_splits)) {
    stop(paste(n_splits, " n_splits is not of type 'numeric'"))
  }
  
  ## Check reference_cell_types are found in the spe object
  unknown_cell_types <- setdiff(reference_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in reference_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  ## Check target_cell_types are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  # Check if there is intersection between reference_cell_types and target_cell_types
  if (length(intersect(reference_cell_types, target_cell_types)) > 0) {
    stop("Cannot have same cells in both reference_cell_types and target_cell_types")
  }
  
  
  spe_coords <- data.frame(spatialCoords(spe))
  
  ## Get dimensions of the window
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  d_lay <- height / n_splits
  
  ## Figure out which 'grid prism number' each cell is inside
  spe$Prism.Num <- floor(spe_coords$Cell.X.Position / d_row) +
    floor(spe_coords$Cell.Y.Position / d_col) * n_splits + 
    floor(spe_coords$Cell.Z.Position / d_lay) * n_splits^2 + 1
  
  ## Calculate cell_proportions for each grid prism
  n_grid_prisms <- n_splits^3
  n_reference_cells_vec <- c()
  n_target_cells_vec <- c()
  grid_prism_cell_proportions <- c()
  
  ## Define data frame which contains all results
  result <- data.frame(matrix(nrow = n_grid_prisms, ncol = 4))
  colnames(result) <- c("reference", "target", "total", "proportion")
  
  for (grid_prism_num in seq(n_grid_prisms)) {
    
    ## Get spe object for current grid_prism
    spe_temp <- spe[ , spe$Prism.Num == grid_prism_num]
    
    ## Get cell_proportion: n_target_cells / (n_target_cells + n_reference_cells)
    n_target_cells <- sum(spe_temp[[feature_colname]] %in% target_cell_types)
    n_reference_cells <- sum(spe_temp[[feature_colname]] %in% reference_cell_types)
    
    ## Case when there are no target or reference cell, result is NA
    if (n_target_cells == 0 && n_reference_cells == 0) {
      grid_prism_cell_proportion <- NA  
    }
    else {
      grid_prism_cell_proportion <- n_target_cells / (n_target_cells + n_reference_cells)
    }
    
    result[grid_prism_num, ] <- c(n_reference_cells, 
                                  n_target_cells, 
                                  n_reference_cells + n_target_cells, 
                                  grid_prism_cell_proportion)
  }
  
  ## Add x, y and z coords of each grid prism to result
  result$prism_x_coord <- ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row
  result$prism_y_coord <- (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col
  result$prism_z_coord <- (floor((seq(n_grid_prisms) - 1) / (n_splits^2)) + 0.5) * d_lay
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous3D(result, "proportion")
    methods::show(fig)
  }
  
  return(result)
}




plot_grid_metrics_continuous3D <- function(grid_metrics, metric_colname) {
  
  ## Color of each dot is related to its entropy
  pal <- colorRampPalette(hcl.colors(n = 5, palette = "Red-Blue", rev = TRUE))
  
  ## Add size column and for NA entropy values, make the size small
  grid_metrics$size <- ifelse(is.na(grid_metrics[[metric_colname]]), 3, 10)
  
  fig <- plot_ly(grid_metrics,
                 type = "scatter3d",
                 mode = 'markers',
                 x = ~prism_x_coord,
                 y = ~prism_y_coord,
                 z = ~prism_z_coord,
                 color = as.formula(paste0('~', metric_colname)),
                 colors = pal(nrow(grid_metrics)),
                 marker = list(size = ~size),
                 symbol = 1,
                 symbols = "square")
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                     yaxis = list(title = 'y'),
                                     zaxis = list(title = 'z')))
  
  return(fig)
}


plot_grid_metrics_discrete3D <- function(grid_metrics, metric_colname) {
  
  ## Define low, medium and high categories
  # Low: between 0 and 1/3
  # Medium: between 1/3 and 2/3
  # High: between 2/3 and 1
  
  grid_metrics$rank <- ifelse(is.na(grid_metrics[[metric_colname]]), "na",
                              ifelse(grid_metrics[[metric_colname]] < 1/3, "low",
                                     ifelse(grid_metrics[[metric_colname]] < 2/3, "medium", "high")))
  grid_metrics$rank <- factor(grid_metrics$rank, c("low", "medium", "high", "na"))
  
  fig <- plot_ly(grid_metrics,
                 type = "scatter3d",
                 mode = 'markers',
                 x = ~prism_x_coord,
                 y = ~prism_y_coord,
                 z = ~prism_z_coord,
                 color = ~rank,
                 colors = c("#AEB6E5", "#BC6EB9", "#A93154", "gray"),
                 symbol = 1,
                 symbols = "square",
                 marker = list(size = 4))
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                     yaxis = list(title = 'y'),
                                     zaxis = list(title = 'z')))
  return(fig)
}



determine_prevalence3D <- function(grid_data,
                                   metric_colname,
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
  
  return(p)
}

determine_prevalence_gradient3D <- function(grid_data,
                                            metric_colname,
                                            plot_image = T) {
  
  # Thresholds range from 0 to 1
  thresholds <- seq(0, 1, 0.01)
  
  # Define result
  result <- data.frame(threshold = thresholds)
  
  prevalences <- c()
  
  for (threshold in thresholds) {
    prevalences <- c(prevalences, determine_prevalence3D(grid_data,
                                                         metric_colname,
                                                         threshold))
  }
  result$prevalence <- prevalences
  
  # Plot
  if (plot_image) {
    fig <- ggplot(result, aes(threshold, prevalence)) +
      geom_line() +
      theme_bw() +
      labs(x = "Threshold",
           y = "Prevalence",
           title = paste("Prevalence vs Threshold (", metric_colname, ")", sep = "")) +
      theme(plot.title = element_text(hjust = 0.5))
    methods::show(fig)
  }
  
  return(result)
}


determine_spatial_autocorrelation3D <- function(grid_data,
                                              metric_colname,
                                              weight_method = "IDW") {
  
  
  ## Get number of grid prisms
  n_grid_prisms <- nrow(grid_data)
  
  ## Get splitting number (should be the cube root of n_grid_prisms)
  n_splits <- (n_grid_prisms)^(1/3)
  
  ## Find the coordinates of each grid prism
  x <- ((seq(n_grid_prisms) - 1) %% n_splits)
  y <- (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits))
  z <- (floor((seq(n_grid_prisms) - 1) / (n_splits^2)))
  grid_prism_coords <- data.frame(x = x, y = y, z = z)
  
  
  weight_matrix <- -1 * apcluster::negDistMat(grid_prism_coords)
  ## Use the inverse distance between two points as the weight (IDW is 'inverse distance weighting')
  if (weight_method == "IDW") {
    weight_matrix <- 1 / weight_matrix
  }
  ## Use binary method: adjacent points get a weight of 1, otherwise, weight of 0
  ## Adjacent points are within sqrt(3) units apart. e.g. (0, 0, 0) vs (1, 1, 1)
  else if (weight_method == "binary") {
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
  
  return(I)
  
}


### Clustering algorithms ----------------------------------------------------
library(alphashape3d)

alpha_hull_clustering3D <- function(spe, 
                                    cell_types_of_interest, 
                                    alpha, 
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
  
  ## Get the alpha hull
  alpha_hull <- ashape3d(as.matrix(spe_subset_coords), alpha = alpha)
  
  if (sum(alpha_hull$triang[, 9]) == 0) stop("alpha value is too small? No alpha hulls identified")
  
  ## Determine which alpha hull cluster each cell_type_of_interest belongs to
  alpha_hull_clusters <- components_ashape3d(alpha_hull)
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), 
                   "Cell.Type" = spe[[feature_colname]],
                   "Cell.ID" = spe[["Cell.ID"]])
  
  df_cell_types_of_interest <- df[df$Cell.Type %in% cell_types_of_interest, ]
  df_other_cell_types <- df[!(df$Cell.Type %in% cell_types_of_interest), ]
  df_cell_types_of_interest$alpha_hull_cluster <- alpha_hull_clusters
  df_other_cell_types$alpha_hull_cluster <- 0
  
  ## Ignore cell_types_of_interest which belong to an alpha hull cluster with less than minimum_cells_in_alpha_hull
  alpha_hull_clusters_table <- table(alpha_hull_clusters)
  maximium_alpha_hull_cluster <- Position(function(x) x < minimum_cells_in_alpha_hull, alpha_hull_clusters_table)
  maximium_alpha_hull_cluster <- as.numeric(names(alpha_hull_clusters_table[maximium_alpha_hull_cluster]))
  
  if (!is.na(maximium_alpha_hull_cluster) && maximium_alpha_hull_cluster != -1) {
    spe_subset_coords <- spe_subset_coords[alpha_hull_clusters >= 1 & alpha_hull_clusters < maximium_alpha_hull_cluster, ]
    
    df_cell_types_of_interest$alpha_hull_cluster <- ifelse(alpha_hull_clusters >= 1 & alpha_hull_clusters < maximium_alpha_hull_cluster, 
                                                           alpha_hull_clusters, 0)
    
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
  spe@metadata$alpha_hull <- list(vertices = vertices, faces = faces, ashape3d_object = alpha_hull)
  
  ## Plot
  if (plot_image) {
    fig <- plot_alpha_hull_clusters3D(spe, feature_colname = feature_colname)
    methods::show(fig)
  }
  
  return(spe)
}



plot_alpha_hull_clusters3D <- function(spe_with_alpha_hull, 
                                       plot_cell_types = NULL,
                                       plot_colours = NULL,
                                       feature_colname = "Cell.Type") {
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe_with_alpha_hull), "Cell.Type" = spe_with_alpha_hull[[feature_colname]])
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) {
    plot_cell_types <- unique(df[["Cell.Type"]])
  }
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(plot_cell_types, spe_with_alpha_hull[[feature_colname]])
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
  
  ## Add points to fig
  fig <- plot_ly() %>%
    add_trace(
      data = df,
      type = "scatter3d",
      mode = 'markers',
      x = ~Cell.X.Position,
      y = ~Cell.Y.Position,
      z = ~Cell.Z.Position,
      marker = list(size = 2),
      color = ~Cell.Type,
      colors = plot_colours
    ) %>% 
    layout(scene = list(xaxis = list(title = 'x'),
                        yaxis = list(title = 'y'),
                        zaxis = list(title = 'z')))
  
  
  ## Get alpha hull numbers (ignoring 0)
  alpha_hull_clusters <- spe_with_alpha_hull$alpha_hull_cluster[spe_with_alpha_hull$alpha_hull_cluster != 0]
  
  # Get number of alpha hulls
  n_alpha_hulls <- length(unique(alpha_hull_clusters))
  
  vertices <- spe_with_alpha_hull@metadata$alpha_hull$vertices
  faces <- data.frame(spe_with_alpha_hull@metadata$alpha_hull$faces)
  alpha_hull_colours <- rainbow(n_alpha_hulls)
  
  ## Add alpha hulls to fig, one by one  
  for (i in seq(n_alpha_hulls)) {
    faces_temp <- faces[faces[ , 1] %in% which(alpha_hull_clusters == i) , ]
    
    ## Ignore the weird cases where some cells represent clusters, but no faces are associated with them??
    if (nrow(faces_temp) == 0) next
    
    # Large alpha hulls should have a lower opacity so they are more visible
    opacity_level <- ifelse(nrow(faces_temp) > 50, 0.05, 0.25)
    
    fig <- fig %>%
      add_trace(
        type = 'mesh3d',
        x = vertices[, 1], 
        y = vertices[, 2], 
        z = vertices[, 3],
        i = faces_temp[, 1] - 1, 
        j = faces_temp[, 2] - 1, 
        k = faces_temp[, 3] - 1,
        opacity = opacity_level,
        facecolor = rep(alpha_hull_colours[i], nrow(faces_temp))
      )
  }
  
  return(fig)
}



library(dbscan)

dbscan_clustering3D <- function(spe,
                                cell_types_of_interest,
                                radius,
                                minimum_cells_in_radius,
                                feature_colname = "Cell.Type",
                                plot_image = T) {
  
  spe_subset <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  spe_subset_coords <- spatialCoords(spe_subset)
  
  db <- dbscan::dbscan(spe_subset_coords, eps = radius, minPts = minimum_cells_in_radius, borderPoints = F)
  
  
  
  # Convert spe object to data frame
  df <- data.frame(spatialCoords(spe),
                   "Cell.Type" = spe[[feature_colname]],
                   "Cell.ID" = spe[["Cell.ID"]])
  
  df_cell_types_of_interest <- df[df$Cell.Type %in% cell_types_of_interest, ]
  df_other_cell_types <- df[!(df$Cell.Type %in% cell_types_of_interest), ]
  df_cell_types_of_interest$dbscan_cluster <- db$cluster
  df_other_cell_types$dbscan_cluster <- 0
  
  ## Convert data frame to spe object
  df <- rbind(df_cell_types_of_interest, df_other_cell_types)
  
  spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = spe@metadata)
  
  ## Plot
  if (plot_image) {
    df$dbscan_cluster <- ifelse(df$dbscan_cluster == 0, "non_cluster", paste("cluster_", df$dbscan_cluster, sep = ""))
    
    fig <- plot_ly(df,
                   type = "scatter3d",
                   mode = 'markers',
                   x = ~Cell.X.Position,
                   y = ~Cell.Y.Position,
                   z = ~Cell.Z.Position,
                   color = ~dbscan_cluster,
                   colors = rainbow(length(unique(df$dbscan_cluster))),
                   marker = list(size = 2)) %>% 
      layout(scene = list(xaxis = list(title = 'x'),
                          yaxis = list(title = 'y'),
                          zaxis = list(title = 'z')))
    
    methods::show(fig)
  }
  
  return(spe)
}



grid_based_clustering3D <- function(spe,
                                    cell_types_of_interest,
                                    n_splits,
                                    feature_colname = "Cell.Type",
                                    plot_image = TRUE) {
  
  # Check if n_splits is numeric
  if (!is.numeric(n_splits)) {
    stop(paste(n_splits, " n_splits is not of type 'numeric'"))
  }
  
  ## Check cell_types_of_interest are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  spe_coords <- data.frame(spatialCoords(spe))
  
  ## Get dimensions of the window
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  d_lay <- height / n_splits
  
  ## Figure out which 'grid prism number' each cell is inside
  spe$Prism.Num <- floor(spe_coords$Cell.X.Position / d_row) +
    floor(spe_coords$Cell.Y.Position / d_col) * n_splits + 
    floor(spe_coords$Cell.Z.Position / d_lay) * n_splits^2 + 1
  
  ## Calculate proportions for each grid prism
  n_grid_prisms <- n_splits^3
  grid_prism_cell_proportions <- c()
  for (grid_prism_num in seq(n_grid_prisms)) {
    
    ## Get spe object for current grid_prism
    spe_temp <- spe[ , spe$Prism.Num == grid_prism_num]
    
    ## Get total number of cells and number of cell_types_of_interest in current grid_prism
    n_total <- ncol(spe_temp)
    n_cell_types_of_interest <- sum(spe_temp[[feature_colname]] %in% cell_types_of_interest)
    
    if (n_total == 0) {
      grid_prism_cell_proportion <- 0
    }
    else {
      grid_prism_cell_proportion <- n_cell_types_of_interest / n_total
    }
    
    grid_prism_cell_proportions <- c(grid_prism_cell_proportions, grid_prism_cell_proportion)
  }
  names(grid_prism_cell_proportions) <- seq(n_grid_prisms)
  
  
  ## Create template for final result
  result <- list()
  n_clusters <- 1
  
  
  ### CLUSTER DETECTION RECURSIVE ALGORITHM LOOP ###
  
  # First, remove all 0s from grid_prism_cell_proportions
  grid_prism_cell_proportions <- grid_prism_cell_proportions[grid_prism_cell_proportions != 0]
  
  while (length(grid_prism_cell_proportions) != 0) {
    # Get the maximum cell proportion and its corresponding grid prism number
    maximum_cell_proportion <- max(grid_prism_cell_proportions)
    maximum_cell_proportion_prism_number <- as.numeric(names(which.max(grid_prism_cell_proportions)))
    
    # Break out the loop if maximum cell proportion is less than 0.25
    if (maximum_cell_proportion < 0.5) break 
    
    # Else, find all the grid prisms adjacent to the maximum cell proportion grid prism. 
    # These are potentially apart of the cluster
    # Adjacent grid prisms must have cell proportion > 0.25 * max cell proportion
    grid_prisms_in_cluster <- determine_grid_prism_numbers_in_cluster3D(maximum_cell_proportion_prism_number,
                                                                        grid_prism_cell_proportions,
                                                                        maximum_cell_proportion,
                                                                        n_splits,
                                                                        c())
    
    # Perform the recursive algorithm on each grid prism potentially apart of the cluster to get a more precise shape of each cluster
    result[[n_clusters]] <- data.frame()
    for (grid_prism in as.numeric(grid_prisms_in_cluster)) {
      
      spe_prism <- spe[ , spe$Prism.Num == grid_prism]
      
      grid_prism_x <- ((grid_prism - 1) %% n_splits) * d_row
      grid_prism_y <- (floor(((grid_prism - 1) %% n_splits^2) / n_splits)) * d_col
      grid_prism_z <- (floor((grid_prism - 1) / n_splits^2)) * d_lay
      
      result[[n_clusters]] <- rbind(result[[n_clusters]], 
                                    grid_based_cluster_recursion3D(spe_prism, 
                                                                   cell_types_of_interest, 
                                                                   0.75 * maximum_cell_proportion,
                                                                   grid_prism_x, grid_prism_y, grid_prism_z,
                                                                   d_row, d_col, d_lay,
                                                                   "Cell.Type",
                                                                   data.frame()))
    }
    
    colnames(result[[n_clusters]]) <- c("x", "y", "z", "l", "w", "h")
    n_clusters <- n_clusters + 1
    
    # Remove grid prisms which have just been examined
    grid_prism_cell_proportions <- grid_prism_cell_proportions[setdiff((names(grid_prism_cell_proportions)), 
                                                                       grid_prisms_in_cluster)]
    
  }
  
  ## Add all the information to the spe
  spe@metadata[["grid_prisms"]] <- result
  spe$Prism.Num <- NULL
  spe$grid_based_cluster <- 0
  cluster_number <- 1
  
  for (cluster_info in result) {
    for (i in seq(nrow(cluster_info))) {
      x <- cluster_info$x[i]
      y <- cluster_info$y[i]
      z <- cluster_info$z[i]
      l <- cluster_info$l[i]
      w <- cluster_info$w[i]
      h <- cluster_info$h[i]
      
      spe$grid_based_cluster <- ifelse(spe_coords$Cell.X.Position >= x &
                                         spe_coords$Cell.X.Position < (x + l) &
                                         spe_coords$Cell.Y.Position >= y &
                                         spe_coords$Cell.Y.Position < (y + w) &
                                         spe_coords$Cell.Z.Position >= z &
                                         spe_coords$Cell.Z.Position < (z + h), 
                                       cluster_number, 
                                       spe$grid_based_cluster)
    }
    cluster_number <- cluster_number + 1
  }
  
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_based_clusters3D(spe, feature_colname = feature_colname)
    methods::show(fig)
  }
  
  return(spe)
}


### Start from the grid_prism with the maximum cell proportion.
## Look left, right, forward, back, up and down and see if that grid_prism has at least 25% of the maximum cell proportion
## If it does, add it to the answer
## Keep doing this until adjacent grid prisms don't have above 25%, or if you hit a boundary, or it has already been removed
## Return a vector containing all the grid prism numbers which COULD be part of the cluster
determine_grid_prism_numbers_in_cluster3D <- function(curr_grid_prism_number, 
                                                      grid_prism_cell_proportions, 
                                                      maximum_cell_proportion,
                                                      n_splits,
                                                      answer) {
  
  ## If answer already has curr_grid_prism_number, go back
  if (as.character(curr_grid_prism_number) %in% answer) return(answer)
  
  grid_prism_numbers <- names(grid_prism_cell_proportions)
  
  ## If curr_grid_prism_number has already been removed from grid_prism_numbers, go back
  if (!(as.character(curr_grid_prism_number) %in% grid_prism_numbers)) return(answer)
  
  
  if (grid_prism_cell_proportions[as.character(curr_grid_prism_number)] > 0.25 * maximum_cell_proportion) {
    
    answer <- c(answer, as.character(curr_grid_prism_number))
    
    ### CHECK RIGHT, LEFT, FORWARD, BACKWARD, UP, DOWN
    ## Need to check if going right, left, forward, backward, up or down is possible
    
    # Right
    if (curr_grid_prism_number%%n_splits != 0) {
      answer <- determine_grid_prism_numbers_in_cluster3D(curr_grid_prism_number + 1,
                                                          grid_prism_cell_proportions,
                                                          maximum_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Left
    if (curr_grid_prism_number%%n_splits != 1) {
      answer <- determine_grid_prism_numbers_in_cluster3D(curr_grid_prism_number - 1,
                                                          grid_prism_cell_proportions,
                                                          maximum_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Forward
    if ((curr_grid_prism_number - 1)%%(n_splits^2) < n_splits^2 - n_splits) {
      answer <- determine_grid_prism_numbers_in_cluster3D(curr_grid_prism_number + n_splits,
                                                          grid_prism_cell_proportions,
                                                          maximum_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Backward
    if (curr_grid_prism_number%%(n_splits^2) > n_splits) {
      answer <- determine_grid_prism_numbers_in_cluster3D(curr_grid_prism_number - n_splits,
                                                          grid_prism_cell_proportions,
                                                          maximum_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Up
    if (curr_grid_prism_number <= n_splits^3 - n_splits^2) {
      answer <- determine_grid_prism_numbers_in_cluster3D(curr_grid_prism_number + n_splits^2,
                                                          grid_prism_cell_proportions,
                                                          maximum_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Down
    if (curr_grid_prism_number > n_splits^2) {
      answer <- determine_grid_prism_numbers_in_cluster3D(curr_grid_prism_number - n_splits^2,
                                                          grid_prism_cell_proportions,
                                                          maximum_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
  }
  
  return(answer)
}



grid_based_cluster_recursion3D <- function(spe, 
                                           cell_types_of_interest,
                                           threshold_cell_proportion,
                                           x, y, z, l, w, h,
                                           feature_colname,
                                           answer) {
  
  # Look at cells only in the current grid prism
  spe_coords <- data.frame(spatialCoords(spe))
  
  spe <- spe[ , spe_coords$Cell.X.Position >= x &
                spe_coords$Cell.X.Position < (x + l) &
                spe_coords$Cell.Y.Position >= y &
                spe_coords$Cell.Y.Position < (y + w) &
                spe_coords$Cell.Z.Position >= z &
                spe_coords$Cell.Z.Position < (z + h)]
  
  # Get cell types from spe grid prism
  spe_cell_types <- spe[[feature_colname]]
  
  # Number of cells in prism is getting too small
  if (length(spe_cell_types) <= 2) return(data.frame())
  
  # Get total cell proportion for chosen cell_types_of_interest
  cell_proportion <- sum(spe_cell_types %in% cell_types_of_interest) / length(spe_cell_types)
  
  # Keep grid prism if cell proportion is above the threshold cell proportion
  if (cell_proportion >= threshold_cell_proportion) {
    return(data.frame(x, y, z, l, w, h))
  }
  
  # some cell_types_of_interest still in the grid prism, check sub-grid prisms (8 to check)
  else if (cell_proportion > 0) {
    # (0, 0, 0)
    answer <- rbind(answer, grid_based_cluster_recursion3D(spe,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x, y, z, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0.5, 0, 0)
    answer <- rbind(answer, grid_based_cluster_recursion3D(spe,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x + l/2, y, z, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0, 0.5, 0)
    answer <- rbind(answer, grid_based_cluster_recursion3D(spe,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x, y + w/2, z, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    # (0.5, 0.5, 0)
    answer <- rbind(answer, grid_based_cluster_recursion3D(spe,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x + l/2, y + w/2, z, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0, 0, 0.5)
    answer <- rbind(answer, grid_based_cluster_recursion3D(spe,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x, y, z + h/2, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0.5, 0, 0.5)
    answer <- rbind(answer, grid_based_cluster_recursion3D(spe,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x + l/2, y, z + h/2, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0, 0.5, 0.5)
    answer <- rbind(answer, grid_based_cluster_recursion3D(spe,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x, y + w/2, z + h/2, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    # (0.5, 0.5, 0.5)
    answer <- rbind(answer, grid_based_cluster_recursion3D(spe,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x + l/2, y + w/2, z + h/2, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    return(answer)
  }
  
  # cell proportion is zero
  else {
    return(data.frame())
  }
}



plot_grid_based_clusters3D <- function(spe_with_grid, 
                                       plot_cell_types = NULL,
                                       plot_colours = NULL,
                                       feature_colname = "Cell.Type") {
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe_with_grid), "Cell.Type" = spe_with_grid[[feature_colname]])
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) {
    plot_cell_types <- unique(df[["Cell.Type"]])
  }
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(plot_cell_types, spe_with_grid[[feature_colname]])
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
  
  ## Add points to fig
  fig <- plot_ly() %>%
    add_trace(
      data = df,
      type = "scatter3d",
      mode = 'markers',
      x = ~Cell.X.Position,
      y = ~Cell.Y.Position,
      z = ~Cell.Z.Position,
      marker = list(size = 2),
      color = ~Cell.Type,
      colors = plot_colours
    ) %>% 
    layout(scene = list(xaxis = list(title = 'x'),
                        yaxis = list(title = 'y'),
                        zaxis = list(title = 'z')))
  
  # Get number of grid-based clusters
  n_grid_based_clusters <- length(spe_with_grid@metadata[["grid_prisms"]])
  
  faces <- data.frame(edge1 = c(1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8),
                      edge2 = c(2, 5, 2, 3, 3, 5, 6, 4 ,7, 6, 7, 4),
                      edge3 = c(6, 6, 4, 4, 7, 7, 2, 2, 5, 5, 3, 3))
  grid_based_colours <- rainbow(n_grid_based_clusters)
  
  ## Add grid-based clusters to fig, one by one  
  for (i in seq(n_grid_based_clusters)) {
    
    grid_based_cluster <- spe_with_grid@metadata[["grid_prisms"]][[i]]
    
    for (j in seq(nrow(grid_based_cluster))) {
      
      x <- grid_based_cluster$x[j]
      y <- grid_based_cluster$y[j]
      z <- grid_based_cluster$z[j]
      l <- grid_based_cluster$l[j]
      w <- grid_based_cluster$w[j]
      h <- grid_based_cluster$h[j]
      vertices <- data.frame(x = c(x, x + l, x, x + l, x, x + l, x, x + l),
                             y = c(y, y, y + w, y + w, y, y, y + w, y + w),
                             z = c(z, z, z, z, z + h, z + h, z + h, z + h))
      
      fig <- fig %>%
        add_trace(
          type = 'mesh3d',
          x = vertices[, 1], 
          y = vertices[, 2], 
          z = vertices[, 3],
          i = faces[, 1] - 1, 
          j = faces[, 2] - 1, 
          k = faces[, 3] - 1,
          opacity = 0.2,
          facecolor = rep(grid_based_colours[i], 12) # Always 12 faces per grid prism
        )      
    }
  }
  
  return(fig)
}








calculate_cell_proportions_of_clusters3D <- function(spe, cluster_colname, feature_colname = "Cell.Type", plot_image = T) {
  
  ## Get cluster numbers (ignoring 0)
  cluster_numbers <- spe[[cluster_colname]][spe[[cluster_colname]] != 0]
  
  ## Get number of clusters
  n_clusters <- length(unique(cluster_numbers))
  
  ## Get different cell types found in the clusters (alphabetical for consistency)
  cell_types <- unique(spe[[feature_colname]][spe[[cluster_colname]] != 0])
  cell_types <- cell_types[order(cell_types)]
  
  ## For each cluster, determine the size and cell proportion of each cluster
  result <- data.frame(matrix(nrow = n_clusters, ncol = 2 + length(cell_types)))
  colnames(result) <- c("cluster_number", "n_cells", cell_types)
  result$cluster_number <- as.character(seq(n_clusters))
  
  for (i in seq(n_clusters)) {
    cells_in_cluster <- spe[[feature_colname]][spe[[cluster_colname]] == i]
    result[i, "n_cells"] <- length(cells_in_cluster)
    
    for (cell_type in cell_types) {
      result[i, cell_type] <- sum(cells_in_cluster == cell_type) / result[i, "n_cells"]
    }
  }
  
  ## Plot
  if (plot_image) {
    plot_result <- reshape2::melt(result, id.vars = c("cluster_number", "n_cells"))
    fig <- ggplot(plot_result, aes(cluster_number, value, fill = variable)) +
      geom_bar(stat = "identity") +
      labs(title = "Cell proportions of each cluster", x = "", y = "Cell proportion") +
      scale_x_discrete(labels = paste("cluster_", result$cluster_number, ", n = ", result$n_cells, sep = "")) +
      guides(fill = guide_legend(title="Cell type")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    
    methods::show(fig)
  }
  
  return(result)
}



calculate_minimum_distances_to_clusters3D <- function(spe, cell_types_inside_cluster, cell_types_outside_cluster, cluster_colname, feature_colname = "Cell.Type", plot_image = T) {
  
  ## Get cluster numbers (ignoring 0)
  cluster_numbers <- spe[[cluster_colname]][spe[[cluster_colname]] != 0]
  
  ## Get number of clusters
  n_clusters <- length(unique(cluster_numbers))
  
  ## For each cell type outside clusters, get their set of coords. These exclude cell types in clusters
  spe_coords <- spatialCoords(spe)
  cluster_rows <- rep(FALSE, nrow(spe_coords))
  for (i in seq(n_clusters)) {
    cluster_rows <- cluster_rows | (spe[[cluster_colname]] == i)
  }
  
  spe_outside_cluster <- spe[ , !cluster_rows]
  cell_types_outside_cluster_coords <- list()
  for (cell_type in cell_types_outside_cluster) {
    cell_types_outside_cluster_coords[[cell_type]] <- spatialCoords(spe_outside_cluster)[spe_outside_cluster[[feature_colname]] == cell_type, ]
  }
  
  ## For each cluster, determine the minimum distance of each cell_type_of_interest  
  result <- vector()
  
  for (i in seq(n_clusters)) {
    cluster_coords <- spe_coords[spe[[cluster_colname]] == i & spe[[feature_colname]] %in% cell_types_inside_cluster, ]
    
    for (cell_type in cell_types_outside_cluster) {
      curr_cell_type_coords <- cell_types_outside_cluster_coords[[cell_type]]
      
      all_closest <- RANN::nn2(data = cluster_coords, 
                               query = curr_cell_type_coords, 
                               k = 1)  
      
      local_dist_mins <- data.frame(
        cluster_number = i,
        cell_type_of_interest = cell_type,
        distance = all_closest$nn.dists
      )
      ## Remove any 0 distance rows
      local_dist_mins <- local_dist_mins[local_dist_mins$distance != 0, ]
      result <- rbind(result, local_dist_mins)
    }
    
    
    ## Plot
    if (plot_image) {
      
      cluster_number_labs <- paste("cluster_", seq(n_clusters), sep = "")
      names(cluster_number_labs) <- seq(n_clusters)
      
      fig <- ggplot(result, aes(x = cell_type_of_interest, y = distance, fill = cell_type_of_interest)) + 
        geom_violin() +
        facet_grid(cluster_number~., scales="free_x", labeller = labeller(cluster_number = cluster_number_labs)) +
        theme_bw() +
        theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none") +
        labs(title="Minimum cell distances to clusters", x = "Cell type", y = "Distance") +
        stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
      
      methods::show(fig)
    }
    
  }
  return(result)
}

calculate_volume_of_clusters3D <- function(spe, cluster_colname, feature_colname = "Cell.Type") {
  
  ## Get cluster numbers (ignoring 0)
  cluster_numbers <- spe[[cluster_colname]][spe[[cluster_colname]] != 0]
  
  ## Get number of clusters
  n_clusters <- length(unique(cluster_numbers))
  
  
  ### 1. Estimate volume of each cluster by density of the window. ------------
  
  ## For each cluster, determine the number of cells in each cluster of each cluster
  result <- data.frame(matrix(nrow = n_clusters, ncol = 2))
  colnames(result) <- c("cluster_number", "n_cells")
  
  for (i in seq(n_clusters)) {
    cells_in_cluster <- spe[[feature_colname]][spe[[cluster_colname]] == i]
    result[i, "n_cells"] <- length(cells_in_cluster)
    
  }
  # result <- result[order(result$n_cells), ]
  # rownames(result) <- seq(n_clusters)
  result$cluster_number <- as.character(seq(n_clusters))
  
  
  ## Assume window is a rectangular prism
  spe_coords <- data.frame(spatialCoords(spe))
  
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  
  window_volume <- length * width * height
  
  result$volume_by_density <- (result$n_cells / ncol(spe)) * window_volume
  
  
  ### 2. If cluster_colname == "alpha_hull_cluster", use the volume method found in the alphashape3d package
  if (cluster_colname == "alpha_hull_cluster") {
    result$volume_by_alpha_hull <- volume_ashape3d(spe@metadata$alpha_hull$ashape3d_object, byComponents = T)
  }
  
  
  ### 3. If cluster_colname == "grid_based_cluster", sum the volume of each grid prism to get volume of each cluster
  if (cluster_colname == "grid_based_cluster") {
    result$volume_by_grid <- 0
    i <- 1
    for (grid_cluster in spe@metadata$grid_prisms) {
      result[i, "volume_by_grid"] <- sum(grid_cluster$l * grid_cluster$w * grid_cluster$h)
      i <- i + 1
    }
  }
  
  return(result)
}



### Assume that clusters have uniform density and that the centre of each cluster is defined by its centre of mass
### Centre of mass can be estimated by taking the average of the x, y, and z coordinates of cells in the cluster

calculate_center_of_clusters3D <- function(spe, cluster_colname) {
  
  ## Get cluster numbers (ignoring 0)
  cluster_numbers <- spe[[cluster_colname]][spe[[cluster_colname]] != 0]
  
  ## Get number of clusters
  n_clusters <- length(unique(cluster_numbers))
  
  ## For each cluster, determine the number of cells in each cluster of each cluster
  result <- data.frame(matrix(nrow = n_clusters, ncol = 4))
  colnames(result) <- c("cluster_number", "Centre.X.Position", "Centre.Y.Position", "Centre.Z.Position")
  
  result$cluster_number <- as.character(seq(n_clusters))
  for (i in seq(n_clusters)) {
    spe_cluster <- spe[ , spe[[cluster_colname]] == i]
    result[i, c("Centre.X.Position", "Centre.Y.Position", "Centre.Z.Position")] <- 
      apply(spatialCoords(spe_cluster), 2, mean)
  }
  
  return(result)
}



plot_cells3D <- function(spe,
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
                 marker = list(size = 2))
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                     yaxis = list(title = 'y'),
                                     zaxis = list(title = 'z')))
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = '', showgrid = T, showaxeslabels = F, showticklabels = F),
                                     yaxis = list(title = '', showgrid = T, showaxeslabels = F, showticklabels = F),
                                     zaxis = list(title = '', showgrid = T, showaxeslabels = F, showticklabels = F)))
  
  
  return (fig)
}

