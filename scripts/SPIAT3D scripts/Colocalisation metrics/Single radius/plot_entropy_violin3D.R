## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_entropy_violin3D <- function(entropy_data, scales = "free_x") {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  entropy <- NULL
  
  fig <- ggplot(entropy_result, aes(x = "", y = entropy)) +
    geom_violin() +
    theme_bw() +
    labs(x = "", y = "Entropy") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plot shows mean ± sd")
  
  return(fig)
}



