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
