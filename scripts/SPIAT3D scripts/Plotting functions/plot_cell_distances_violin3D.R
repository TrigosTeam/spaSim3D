## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_cell_distances_violin3D <- function(cell_to_cell_dist, scales = "free_x") {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Pair <- Distance <- NULL
  
  ggplot(cell_to_cell_dist, aes(x = Pair, y = Distance)) + geom_violin() +
    facet_wrap(~Pair, scales=scales) +
    theme_bw() +
    theme(axis.text.x=element_blank())
  
}
