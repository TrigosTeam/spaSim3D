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
  result_entropies <- result / result$Total
  result_entropies <- -1 * result_entropies  * log(result_entropies, length(target_cell_types))
  
  ## Calculate total entropy for each row in result data frame
  result$Entropy <- apply(result_entropies, 1, sum)
  result$Entropy <- replace(result$Entropy, is.nan(result$Entropy), 0)
  
  # Plot
  if (plot_image) {
    plot(rownames(result), result$Entropy, type = "l", xlab = "Radius", ylab = "Entropy Aggregated")
  }
  
  return (result)
    
}

