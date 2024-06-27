calculate_cell_proportions_in_neighbourhood_gradient3D <- function(spe, 
                                                                   reference_cell_type, 
                                                                   target_cell_types, 
                                                                   radii, 
                                                                   feature_colname = "Cell.Type",
                                                                   plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (radius in seq(radii)) {
    cell_proportions_neighbourhood_data <- calculate_cell_proportions_in_neighbourhood3D(spe,
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