calculate_cross_K_gradient3D <- function(spe, 
                                         reference_cell_type, 
                                         target_cell_type, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = 2))
  colnames(result) <- c("observed_cross_K", 
                        "expected_cross_K")
  
  for (radius in seq(radii)) {
    cross_K_data <- calculate_cross_K3D(spe,
                                        reference_cell_type,
                                        target_cell_type,
                                        radius,
                                        feature_colname)
    
    result[radius, ] <- cross_K_data
  }
  
  # Add a radius column to the result
  result$radius <- seq(radii)
  
  if (plot_image) {
    plot_result <- reshape2::melt(result, "radius", c("observed_cross_K", "expected_cross_K"))
    
    fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
      geom_line() +
      labs(x = "Radius", y = "Cross K-function value") +
      scale_colour_discrete(name = "", labels = c("Observed cross K", "Expected CSR cross K")) +
      theme_bw()
    
    methods::show(fig)

  }
  
  return(result)
}