## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_cell_distances_violin3D <- function(cell_to_cell_dist, scales = "free_x") {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  pair <- distance <- NULL
  
  fig <- ggplot(cell_to_cell_dist, aes(x = pair, y = distance)) + 
    geom_violin() +
    facet_wrap(~pair, scales=scales) +
    theme_bw() +
    theme(axis.text.x=element_blank()) +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plots show mean ± sd")
  
  return(fig)
}
