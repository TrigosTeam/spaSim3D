calculate_alpha_hull_cell_proportions3D <- function(spe_with_alpha_hull, feature_colname = "Cell.Type", plot_image = T) {
  
  ## Get alpha hull numbers (ignoring -1)
  alpha_hull_numbers <- spe_alpha_hull$alpha_hull_number[spe_alpha_hull$alpha_hull_number != -1]
  
  ## Get number of alpha hulls
  n_alpha_hulls <- length(unique(alpha_hull_numbers))
  
  ## Get different cell types found in the alpha hulls
  cell_types <- unique(spe_with_alpha_hull[[feature_colname]][spe_alpha_hull$alpha_hull_number != -1])
  
  ## For each alpha hull, determine the size and cell proportion of each alpha hull
  result <- data.frame(matrix(nrow = n_alpha_hulls, ncol = 2 + length(cell_types)))
  colnames(result) <- c("alpha_hull_number", "n_cells", cell_types)
  
  for (i in seq(n_alpha_hulls)) {
    cells_in_alpha_hull <- spe_with_alpha_hull[[feature_colname]][spe_with_alpha_hull$alpha_hull_number == i]
    result[i, "n_cells"] <- length(cells_in_alpha_hull)
    
    for (cell_type in cell_types) {
      result[i, cell_type] <- sum(cells_in_alpha_hull == cell_type) / result[i, "n_cells"]
    }
  }
  
  result <- result[order(result$n_cells), ]
  rownames(result) <- seq(n_alpha_hulls)
  result$alpha_hull_number <- as.character(seq(n_alpha_hulls))
  
  ## Plot
  if (plot_image) {
    plot_result <- reshape2::melt(result, id.vars = c("alpha_hull_number", "n_cells"))
    fig <- ggplot(plot_result, aes(alpha_hull_number, value, fill = variable)) +
      geom_bar(stat = "identity") +
      labs(title = "Cell proportions of each alpha hull", x = "", y = "Cell proportion") +
      scale_x_discrete(labels = paste("n =", result$n_cells)) +
      guides(fill = guide_legend(title="Cell type")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    
    methods::show(fig)
  }
  
  return(result)
}

