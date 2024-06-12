calculate_mixing_scores_gradient3D <- function(data, 
                                               reference_cell_type, 
                                               target_cell_type, 
                                               radii = 20, 
                                               feature_colname = "Cell.Type",
                                               plot_image = TRUE) {
  
  
  df <- data.frame(matrix(nrow = radii, ncol = 8))
  df.cols <- c("Reference", 
               "Target", 
               "Number_of_reference_cells",
               "Number_of_target_cells", 
               "Reference_target_interaction",
               "Reference_reference_interaction", 
               "Mixing_score", 
               "Normalised_mixing_score")
  colnames(df) <- df.cols
  
  for (radius in seq(radii)) {
    mixing_scores <- calculate_mixing_scores3D(data,
                                               reference_cell_type,
                                               target_cell_type,
                                               radius,
                                               feature_colname)

    df[radius, ] <- mixing_scores
  }
  
  
  if (plot_image) {
    plot(seq(radii), df[["Normalised_mixing_score"]], type = "l", xlab = "Radius", ylab = "Normalised Mixing Score")
    abline(a = 1, b = 0, col = "red", lwd = 2, lty = 2)
  }
  
  return (df)
}
