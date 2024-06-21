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
    plot(result[["radius"]], result[["normalised_mixing_score"]], 
         type = "l", 
         xlab = "Radius", 
         ylab = "Normalised mixing score",
         ylim = c(0, max(result[["normalised_mixing_score"]], 1)),
         col = "red")
    abline(a = 1, b = 0, col = "blue", lty = 2)
    legend(0, 0.95, legend = c("Observed normalised mixing score", "Expected CSR normalised mixing score"), col = c("red", "blue"), lty = c(1, 2))
  }
  
  return(result)
}
