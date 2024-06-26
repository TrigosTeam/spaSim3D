calculate_alpha_hull_cell_proportions3D <- function(spe_with_alpha_hull, feature_colname = "Cell.Type", plot_image = T) {
  
  ## Get alpha hull numbers (ignoring -1)
  alpha_hull_numbers <- spe_alpha_hull$alpha_hull_number[spe_alpha_hull$alpha_hull_number != -1]
  
  ## Get number of alpha hulls
  n_alpha_hulls <- length(unique(alpha_hull_numbers))
  
  ## Get different cell types found in the alpha hulls
  cell_types <- unique(spe_with_alpha_hull[[feature_colname]][spe_alpha_hull$alpha_hull_number != -1])
  
  ## For each alpha hull, determine the size and cell proportion of each alpha hull
  result <- data.frame(matrix(nrow = n_alpha_hulls, ncol = 1 + length(cell_types)))
  colnames(result) <- c("n_cells", cell_types)
  
  for (i in seq(n_alpha_hulls)) {
    cells_in_alpha_hull <- spe_with_alpha_hull[[feature_colname]][spe_with_alpha_hull$alpha_hull_number == i]
    result[i, "n_cells"] <- length(cells_in_alpha_hull)
    
    for (cell_type in cell_types) {
      result[i, cell_type] <- sum(cells_in_alpha_hull == cell_type) / result[i, "n_cells"]
    }
  }
  
  result <- result[order(result$n_cells), ]
  rownames(result) <- seq(n_alpha_hulls)
  
  ## Plot
  if (plot_image) {
    plot_result <- reshape2::melt(result, id.vars = c("n_cells"))
    
    curr_n_cell <- plot_result[1, "n_cells"]
    len <- 0
    for (i in seq(nrow(plot_result) / length(cell_types))) {
      n_cell <- plot_result[i, "n_cells"]
      
      if (curr_n_cell == n_cell) len <- len + 1
      else {
        plot_result[plot_result$n_cells == curr_n_cell, "value"] <- plot_result[plot_result$n_cells == curr_n_cell, "value"] / len
        len <- 1
        curr_n_cell <- n_cell
      }
    }
    
    plot_result$n_cells <- factor(as.character(plot_result$n_cells), 
                                  levels = as.character(unique(plot_result$n_cells)[order(unique(plot_result$n_cells))]))
    
    fig <- ggplot(plot_result, aes(n_cells, value, fill = variable)) +
      geom_bar(stat = "identity") +
      labs(title = "Cell proportions of each alpha hull", x = "Number of cells in alpha hull", y = "Cell proportion") +
      guides(fill = guide_legend(title="Cell type")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    
    methods::show(fig)
  }
  
  return(result)
}

