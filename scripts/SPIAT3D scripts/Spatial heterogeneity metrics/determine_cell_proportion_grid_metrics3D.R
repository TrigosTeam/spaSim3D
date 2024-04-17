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
                   marker = list(size = ~size),
                   symbol = 1,
                   symbols = "square")
    
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                       yaxis = list(title = 'y'),
                                       zaxis = list(title = 'z')))
    
    print(fig)
    
  }
  
  return (result)
}
