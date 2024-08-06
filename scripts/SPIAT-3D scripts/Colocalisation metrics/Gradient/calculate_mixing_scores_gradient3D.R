calculate_mixing_scores_gradient3D <- function(spe, 
                                               reference_cell_type, 
                                               target_cell_type, 
                                               radii, 
                                               feature_colname = "Cell.Type",
                                               plot_image = TRUE) {
  
  result <- data.frame(matrix(nrow = radii, ncol = 8))
  colnames(result) <- c("ref_cell_type", 
                        "tar_cell_type", 
                        "n_ref_cells",
                        "n_tar_cells", 
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions", 
                        "mixing_score", 
                        "normalised_mixing_score")
  
  for (radius in seq(radii)) {
    mixing_scores <- calculate_mixing_scores3D(spe,
                                               reference_cell_type,
                                               target_cell_type,
                                               radius,
                                               feature_colname)

    result[radius, ] <- mixing_scores
  }
  
  # Add a radius column to the result
  result$radius <- seq(radii)
  
  if (plot_image) {
    plot_result1 <- result
    plot_result1$expected_normalised_mixing_score <- 1
    plot_result1 <- reshape2::melt(plot_result1, "radius", c("normalised_mixing_score", "expected_normalised_mixing_score"))
    
    fig1 <- ggplot(plot_result1, aes(x = radius, y = value, color = variable)) +
      geom_line() +
      labs(x = "Radius", y = "Normalised mixing score (NMS)") +
      scale_colour_discrete(name = "", labels = c("Observed NMS", "Expected CSR NMS")) +
      theme_bw()
    
    plot_result2 <- reshape2::melt(result, "radius", "mixing_score")
    
    fig2 <- ggplot(plot_result2, aes(x = radius, y = value, color = variable)) +
      geom_line() +
      labs(x = "Radius", y = "Mixing score (MS)") +
      scale_colour_discrete(name = "", labels = c("Observed MS")) +
      theme_bw()
    
    combined_fig <- plot_grid(fig1, fig2, nrow = 2)
    
    methods::show(combined_fig)
  }
  
  return(result)
}
