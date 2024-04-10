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
