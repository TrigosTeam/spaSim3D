### Basic metrics -------------------------------------------------------------

calculate_cell_proportions2D <- function(spe,
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


calculate_entropy_background2D <- function(spe,
                                           cell_types_of_interest, 
                                           feature_colname = "Cell.Type") {
  
  if (length(cell_types_of_interest) == 0) return(NA)
  if (length(cell_types_of_interest) == 1) return(0)
  
  cell_proportions_data <- calculate_cell_proportions2D(spe, cell_types_of_interest, feature_colname, FALSE)
  
  # Calculate entropy of the entire image
  entropy <- -1 * sum(cell_proportions_data$proportion * log(cell_proportions_data$proportion, length(cell_proportions_data$proportion)))
  
  return(entropy) 
}


### Cell colocalisation metrics -----------------------------------------------
calculate_pairwise_distances_between_cell_types2D <- function(spe,
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
  for (cell_type in unique(df[["Cell.Type"]])) {
    cell_types[[cell_type]] <- as.character(df$Cell.ID[df[["Cell.Type"]] == cell_type])
  }
  
  # Calculate cell to cell distances
  dist_all <- -1 * apcluster::negDistMat(df[, c("Cell.X.Position",
                                                "Cell.Y.Position")])
  
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
    fig <- plot_cell_distances_violin2D(cell_to_cell_dist_all)
    methods::show(fig)
  }
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types2D(cell_to_cell_dist_all))  
  }
  
  return(cell_to_cell_dist_all)
}


## Please ensure there is no factoring in any of the columns!!!

calculate_minimum_distances_between_cell_types2D <- function(spe,
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
                               c("Cell.X.Position", "Cell.Y.Position")]
    
    all_cell_type2_coord <- df[df[, "Cell.Type"] == name2, 
                               c("Cell.X.Position", "Cell.Y.Position")]
    
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
    fig <- plot_cell_distances_violin2D(result)
    methods::show(fig)
  }
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types2D(result))  
  }
  
  return(result)
}


summarise_distances_between_cell_types2D <- function(df) {
  
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
plot_cell_distances_violin2D <- function(cell_to_cell_dist, scales = "free_x") {
  
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




calculate_mixing_scores2D <- function(spe, 
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
    reference_cells <- df[df[["Cell.Type"]] == reference_cell_type, ]
    
    for (target_cell_type in target_cell_types) {
      
      # Get all info for cells of target cell_type      
      target_cells <- df[df[["Cell.Type"]] == target_cell_type, ]
      
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
        
        reference_cell_coords <- reference_cells[, c("Cell.X.Position", "Cell.Y.Position")]
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
        reference_cell_coords <- reference_cells[, c("Cell.X.Position", "Cell.Y.Position")]
        target_cell_coords <- target_cells[, c("Cell.X.Position", "Cell.Y.Position")]
        
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


calculate_cells_in_neighbourhood2D <- function(spe, 
                                               reference_cell_type, 
                                               target_cell_types, 
                                               radius, 
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
  reference_cell_coords <- reference_cells[, c("Cell.X.Position", "Cell.Y.Position")]
  rownames(reference_cell_coords) <- reference_cells$Cell.ID
  
  result <- data.frame(matrix(nrow = nrow(reference_cells), ncol = 0))
  
  for (target_cell_type in target_cell_types) {
    ## Get df for target cells
    target_cells <- df[df[[feature_colname]] == target_cell_type, ]
    target_cell_coords <- target_cells[, c("Cell.X.Position","Cell.Y.Position")]
    
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
    print(summarise_cells_in_neighbourhood2D(result))    
  }
  
  
  ## Plot
  if (plot_image) {
    fig <- plot_cells_in_neighbourhood_violin2D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}


calculate_cells_in_neighbourhood_gradient2D <- function(spe, 
                                                        reference_cell_type, 
                                                        target_cell_types, 
                                                        radii, 
                                                        feature_colname = "Cell.Type",
                                                        plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (radius in seq(radii)) {
    cells_in_neighborhood_data <- calculate_cells_in_neighbourhood2D(spe,
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



summarise_cells_in_neighbourhood2D <- function(cells_in_neighbourhood_data) {
  
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


calculate_cells_in_neighbourhood_proportions2D <- function(spe, 
                                                           reference_cell_type, 
                                                           target_cell_types, 
                                                           radius, 
                                                           feature_colname = "Cell.Type") {
  
  ## Get 'count' neighbourhood data
  cells_in_neighbourhood_data <- calculate_cells_in_neighbourhood2D(spe,
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


calculate_cells_in_neighbourhood_proportions_gradient2D <- function(spe, 
                                                                    reference_cell_type, 
                                                                    target_cell_types, 
                                                                    radii, 
                                                                    feature_colname = "Cell.Type",
                                                                    plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (radius in seq(radii)) {
    cell_proportions_neighbourhood_data <- calculate_cells_in_neighbourhood_proportions2D(spe,
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
plot_cells_in_neighbourhood_violin2D <- function(cells_in_neighbourhood_data, reference_cell_type, scales = "free_x") {
  
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




calculate_entropy2D <- function(spe,
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
  cells_in_neighbourhood_data <- calculate_cells_in_neighbourhood2D(spe,
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
    fig <- plot_entropy_violin2D(cells_in_neighbourhood_data)
    methods::show(fig)
  }
  
  return(cells_in_neighbourhood_data)
}


## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_entropy_violin2D <- function(entropy_data, scales = "free_x") {
  
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







calculate_cross_K2D <- function(spe, 
                                reference_cell_type, 
                                target_cell_type, 
                                radius, 
                                feature_colname = "Cell.Type") {
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), 
                   "Cell.Type" = spe[[feature_colname]])
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% df$Cell.Type)) {
    stop(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
  }
  
  ## For target_cell_type, check it is found in the spe object
  if (!(target_cell_type %in% df$Cell.Type)) {
    stop(paste("The target_cell_type", target_cell_type,"is not found in the spe object"))
  }
  
  
  cells_in_neighbourhood_data <- calculate_cells_in_neighbourhood2D(spe,
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
  ## Get Area of the window the cells are in
  area <- length * width
  
  
  ## Get observed cross K-function
  observed_cross_K <- (area * n_ref_tar_interactions) / (n_ref_cells * n_tar_cells)
  
  ## Get expected cross K-function
  expected_cross_K <- pi * radius^2
  
  result <- data.frame(observed_cross_K = observed_cross_K,
                       expected_cross_K = expected_cross_K)
  
  return(result)
}


calculate_mixing_scores_gradient2D <- function(spe,
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
    mixing_scores <- calculate_mixing_scores2D(spe,
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





calculate_cross_K_gradient2D <- function(spe, 
                                         reference_cell_type, 
                                         target_cell_type, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = 2))
  colnames(result) <- c("observed_cross_K", 
                        "expected_cross_K")
  
  for (radius in seq(radii)) {
    cross_K_data <- calculate_cross_K2D(spe,
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


calculate_entropy_gradient2D <- function(spe,
                                         reference_cell_type,
                                         target_cell_types,
                                         radii,
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (radius in seq(radii)) {
    cells_in_neighbourhood_data <- calculate_cells_in_neighbourhood2D(spe,
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
    expected_entropy <- calculate_entropy_background2D(spe, target_cell_types, feature_colname)
    
    plot(result$radius, result$entropy, type = "l", col = "red", 
         xlim = c(0, radius), ylim = c(0, max(result$entropy)),
         xlab = "Radius", ylab = "Entropy")
    abline(a = expected_entropy, b = 0, col = "blue", lty = 2)
    legend(0, max(result$entropy, expected_entropy), legend = c("Observed entropy", "Expected CSR entropy"), col = c("red", "blue"), lty = c(1, 2))
  }
  
  return(result)
}


plot_cross_K_gradient_ratio2D <- function(cross_K_gradient_results) {
  
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
determine_entropy_grid_metrics2D <- function(spe, 
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
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  
  ## Figure out which 'grid prism number' each cell is inside
  spe$Prism.Num <- floor(spe_coords$Cell.X.Position / d_row) +
    floor(spe_coords$Cell.Y.Position / d_col) * n_splits
  
  ## Get number of grid prisms
  n_grid_prisms <- n_splits^2
  
  ## Define data frame which contains all results
  result <- data.frame(matrix(nrow = n_grid_prisms, ncol = (length(cell_types_of_interest) + 2)))
  colnames(result) <- c(cell_types_of_interest, "total", "entropy")
  
  ## Calculate entropy for each grid prism
  for (grid_prism_num in seq(n_grid_prisms)) {
    
    ## Get spe object for current grid_prism
    spe_temp <- spe[ , spe$Prism.Num == grid_prism_num]
    
    ## Get cell_types_of_interest found in the sub-spe object
    temp_cell_types_of_interest <- intersect(cell_types_of_interest, unique(spe_temp[[feature_colname]]))
    
    grid_prism_entropy <- calculate_entropy_background2D(spe_temp,
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
  
  ## Plot
  if (plot_image) {
    
    fig <- plot_grid_metrics_continuous2D(result, "entropy")
    methods::show(fig)
  }
  
  return(result)
}




determine_cell_proportion_grid_metrics2D <- function(spe, 
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
  
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  
  ## Figure out which 'grid prism number' each cell is inside
  spe$Prism.Num <- floor(spe_coords$Cell.X.Position / d_row) +
    floor(spe_coords$Cell.Y.Position / d_col) * n_splits
  
  ## Calculate cell_proportions for each grid prism
  n_grid_prisms <- n_splits^2
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
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous2D(result, "proportion")
    methods::show(fig)
  }
  
  return(result)
}




plot_grid_metrics_continuous2D <- function(grid_metrics, metric_colname) {
  
  ## Color of each dot is related to its entropy
  pal <- colorRampPalette(hcl.colors(n = 5, palette = "Red-Blue", rev = TRUE))
  
  ## Plot
  fig <- ggplot(grid_metrics, aes(prism_x_coord, prism_y_coord, color = !!sym(metric_colname))) +
    geom_point(size = 25, shape = 15) +
    scale_colour_gradientn(colours = pal(100), limits = c(0, 1))
    labs(x = 'x', y = 'y') +
    theme_bw()
  
  return(fig)
}


plot_grid_metrics_discrete2D <- function(grid_metrics, metric_colname) {
  
  ## Define low, medium and high categories
  # Low: between 0 and 1/3
  # Medium: between 1/3 and 2/3
  # High: between 2/3 and 1
  
  grid_metrics$rank <- ifelse(is.na(grid_metrics[[metric_colname]]), "na",
                              ifelse(grid_metrics[[metric_colname]] < 1/3, "low",
                                     ifelse(grid_metrics[[metric_colname]] < 2/3, "medium", "high")))
  grid_metrics$rank <- factor(grid_metrics$rank, c("low", "medium", "high", "na"))
  
  ## Plot
  fig <- ggplot(grid_metrics, aes(prism_x_coord, prism_y_coord, color = rank)) +
    geom_point(size = 25, shape = 15) +
    scale_color_manual(values = c("#AEB6E5", "#BC6EB9", "#A93154", "gray")) +
    labs(x = 'x', y = 'y') +
    theme_bw()

  return(fig)
}



determine_prevalence2D <- function(grid_data,
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

determine_prevalence_gradient2D <- function(grid_data,
                                            metric_colname,
                                            plot_image = T) {
  
  # Thresholds range from 0 to 1
  thresholds <- seq(0.01, 1, 0.01)
  
  # Define result
  result <- data.frame(threshold = thresholds)
  
  prevalences <- c()
  
  for (threshold in thresholds) {
    prevalences <- c(prevalences, determine_prevalence2D(grid_data,
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
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(0, 100)
    methods::show(fig)
  }
  
  return(result)
}


determine_spatial_autocorrelation2D <- function(grid_data,
                                                metric_colname,
                                                weight_method = "IDW") {
  
  ## Get number of grid prisms
  n_grid_prisms <- nrow(grid_data)
  
  ## Get splitting number (should be the square root of n_grid_prisms)
  n_splits <- (n_grid_prisms)^(1/2)
  
  ## Find the coordinates of each grid prism
  x <- ((seq(n_grid_prisms) - 1) %% n_splits)
  y <- (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits))
  grid_prism_coords <- data.frame(x = x, y = y)
  
  ## Subset for non NA rows
  grid_prism_coords <- grid_prism_coords[!is.na(grid_data[[metric_colname]]), ]
  grid_data <- grid_data[!is.na(grid_data[[metric_colname]]), ]
  
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
  
  data_mean <- mean(grid_data[, metric_colname])
  
  numerator <- 0
  denominator <- 0
  
  for (i in seq(nrow(grid_data))) {
    
    for (j in seq(nrow(grid_data))) {
      
      numerator <- numerator + weight_matrix[i, j] * 
        (grid_data[i, metric_colname] - data_mean) * 
        (grid_data[j, metric_colname] - data_mean)
      
    }
    denominator <- denominator + (grid_data[i, metric_colname] - data_mean)^2
  }
  
  
  I <- (nrow(grid_data) * numerator) / (sum(weight_matrix) * denominator)
  
  return(I)
  
}


### Plotting function ---------------------------------------------------------

plot_cells2D <- function(spe,
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
  fig <- ggplot(df, aes(Cell.X.Position, Cell.Y.Position, color = Cell.Type), size = 3) +
    geom_point() +
    scale_color_manual(values = plot_colours) +
    labs(x = 'x', y = 'y') +
    theme_bw()
  
  
  return(fig)
}





