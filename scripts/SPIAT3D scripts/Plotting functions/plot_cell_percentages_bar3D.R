plot_cell_percentages_bar3D <- function(cell_proportions) {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Cell.Type <- Percentage <- Percentage_label <- NULL
  
  cell_proportions$Percentage_label <- round(cell_proportions$Percentage, digits=1)
  
  cell_percentages_full_plot <-
    ggplot(cell_proportions,
           aes(x = stats::reorder(Proportion_Name, Percentage), 
               y = Percentage, fill = Cell.Type)) +
    geom_bar(stat = 'identity') +
    theme_bw() +
    theme() +
    xlab("Cell Type") + 
    ylab("Proportion of cells") +
    geom_text(aes(label = Percentage_label), 
              position = position_stack(vjust = 0.5), size = 3) +
    coord_flip() 
  
  return (cell_percentages_full_plot)
}
