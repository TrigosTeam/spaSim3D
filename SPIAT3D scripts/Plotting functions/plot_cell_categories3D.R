plot_cell_categories3D <- function(data,
                                   cell_types_of_interest = NULL,
                                   colour_vector = NULL,
                                   size = 4,
                                   include_cell_types_of_no_interest = FALSE,
                                   feature_colname = "Cell.Type") {
  
  if (is.null(cell_types_of_interest)) {
    cell_types_of_interest <- unique(data$Cell.Type)
  }
  
  if (is.null(colour_vector)) {
    colour_vector <- hcl.colors(length(cell_types_of_interest), "Batlow")
  }
  
  if (length(cell_types_of_interest) != length(colour_vector)) {
    stop("Length of cell_types_of_interest is not equal to length of colour_vector")
  }
  
  ## Including non-interest cell types
  ## Define cell.id of non-interest cell types as "No Interest"
  if (include_cell_types_of_no_interest) {
    cell_types_of_non_interest <- setdiff(unique(data[[feature_colname]]), cell_types_of_interest)
    
    data[data[[feature_colname]] %in% cell_types_of_non_interest, feature_colname] <- "No Interest"
    
    ## Add "No Interest" as a cell type of interest
    cell_types_of_interest <- c(cell_types_of_interest, "No Interest")
    
    ## Use lightgray for "No Interest" cell types
    colour_vector <- c(colour_vector, "lightgray")
  }
  ## Excluding non-interest cell types
  ## Subset data to only include cell types of interest
  else {
    data <- data[data[[feature_colname]] %in% cell_types_of_interest, ]
  }
  
  ## Add color column to data
  data$Color <- colour_vector[match(data$Cell.Type, cell_types_of_interest)]
  
  ## Plot
  plot3d(data$Cell.X.Position,
         data$Cell.Y.Position,
         data$Cell.Z.Position, 
         "x",
         "y",
         "z",
         col = data$Color, 
         size = size)

  # Add Legend
  legend3d("topright", 
           legend = cell_types_of_interest, 
           pch = 16, 
           col = colour_vector, 
           cex = 0.8,
           inset = 0.05)
}


