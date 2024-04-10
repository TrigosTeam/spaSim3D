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
