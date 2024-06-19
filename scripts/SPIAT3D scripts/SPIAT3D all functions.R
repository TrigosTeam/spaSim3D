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
    
    g <- ggplot(cell_proportions, aes(x = factor(cell_type, cell_type), y = percentage, fill = cell_type)) +
      geom_bar(stat='identity') + 
      theme_bw() +
      labs(title="Cell proportions", x = "Cell type", y = "Percentage") +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "none") +
      geom_text(aes(label = labels), vjust = 0)
    
    methods::show(g)
  }
  
  # Print short summary description
  print(cell_proportions[ , c("cell_type", "frequency", "percentage")])
  
  return(cell_proportions)
}

calculate_pairwise_distances_between_cell_types3D <- function(spe,
                                                              cell_types_of_interest = NULL,
                                                              feature_colname = "Cell.Type",
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
  print(summarise_distances_between_cell_types3D(cell_to_cell_dist_all))
  
  return(cell_to_cell_dist_all)
}

## Please ensure there is no factoring in any of the columns!!!

calculate_minimum_distances_between_cell_types3D <- function(spe,
                                                             cell_types_of_interest = NULL,
                                                             feature_colname = "Cell.Type",
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
  print(summarise_distances_between_cell_types3D(result))
  
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
    facet_wrap(~pair, scales=scales) +
    theme_bw() +
    theme(axis.text.x=element_blank(), plot.title = element_text(hjust = 0.5)) +
    labs(title="Cell distances", x = "Reference/Target pair", y = "Distance") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plots show mean ± sd")
  
  return(fig)
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
  
  return (fig)
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

calculate_mixing_scores_gradient3D <- function(spe, 
                                               reference_cell_type, 
                                               target_cell_type, 
                                               radii = 20, 
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
    plot(result[["radius"]], result[["normalised_mixing_score"]], type = "l", xlab = "Radius", ylab = "Normalised Mixing Score")
    abline(a = 1, b = 0, col = "red", lwd = 2, lty = 2)
  }
  
  return(result)
}


calculate_cells_in_neighborhood3D <- function(spe, 
                                              reference_cell_type, 
                                              target_cell_types, 
                                              radius, 
                                              feature_colname = "Cell.Type",
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
  
  ## Show summarised results
  print(summarise_cells_in_neighborhood3D(result))
  
  ## Plot
  if (plot_image) {
    fig <- plot_cells_in_neighborhood_violin3D(result)
    methods::show(fig)
  }
  
  return(result)
}


summarise_cells_in_neighborhood3D <- function(cells_in_neighborhood_data) {
  
  df <- cells_in_neighborhood_data
  
  ## Target cell types will be all the columns except the first column
  target_cell_types <- colnames(cells_in_neighborhood_data)[c(-1)]
  
  ## Set up data frame for summarised_results list
  df <- data.frame(row.names = c("mean", "min", "max", "median", "st_dev"))
  
  for (target_cell_type in target_cell_types) {
    
    ## Get statistical measures for each target cell type
    target_cell_type_values <- cells_in_neighborhood_data[[target_cell_type]]
    df[[target_cell_type]] <- c(mean(target_cell_type_values),
                                min(target_cell_type_values),
                                max(target_cell_type_values),
                                median(target_cell_type_values),
                                sd(target_cell_type_values))
    
  }
  
  return (data.frame(t(df)))
}



## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_cells_in_neighborhood_violin3D <- function(cells_in_neighborhood_data, scales = "free_x") {
  
  ## Target cell types will be all the columns except the first column
  target_cell_types <- colnames(cells_in_neighborhood_data)[c(-1)]
  
  df <- reshape2::melt(cells_in_neighborhood_data, measure.vars = target_cell_types)
  colnames(df) <- c("ref_cell_id", "tar_cell_type", "count")
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  tar_cell_type <- count <- NULL
  
  fig <- ggplot(df, aes(x = tar_cell_type, y = count)) + 
    geom_violin() +
    facet_wrap(~tar_cell_type, scales=scales) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(title="Cells in neighbourhood", x = "Target cell type", y = "Number of cells") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plots show mean ± sd")
  
  return(fig)
}




