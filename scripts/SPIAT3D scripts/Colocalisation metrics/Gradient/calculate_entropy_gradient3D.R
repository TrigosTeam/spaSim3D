calculate_entropy_gradient3D <- function(spe,
                                         reference_cell_type,
                                         target_cell_types,
                                         radii,
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (radius in seq(radii)) {
    cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                    reference_cell_type,
                                                                    target_cell_types,
                                                                    radius,
                                                                    feature_colname,
                                                                    FALSE,
                                                                    FALSE)
    
    cells_in_neighbourhood_df$ref_cell_id <- NULL
    result[radius, ] <- apply(cells_in_neighbourhood_df, 2, sum)
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

    plot_result <- result
    plot_result$expected_entropy <- calculate_entropy_background3D(spe, target_cell_types, feature_colname)
    plot_result <- reshape2::melt(plot_result, "radius", c("entropy", "expected_entropy"))
    
    fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
      geom_line() +
      labs(x = "Radius", y = "Entropy") +
      scale_colour_discrete(name = "", labels = c("Observed entropy", "Expected CSR entropy")) +
      theme_bw()
    
    methods::show(fig)
  }
  
  return(result)
}
