calculate_entropy_gradient3D <- function(spe,
                                         reference_cell_type,
                                         target_cell_types,
                                         radii,
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = 1))
  colnames(result) <- "entropy"
  
  for (radius in seq(radii)) {
    entropy_df <- calculate_entropy3D(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radius,
                                      feature_colname,
                                      FALSE)
    
    result[radius, "entropy"] <- mean(entropy_df$entropy)
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
