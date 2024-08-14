plot_cells_in_neighbourhood_gradient3D <- function(cells_in_neighbourhood_gradient_df) {
  
  plot_result <- reshape2::melt(cells_in_neighbourhood_gradient_df, "radius")
  
  fig <- ggplot(plot_result, aes(radius, value, color = variable)) + 
    geom_line() + 
    labs(x = "Radius", y = "Average cells in neighbourhood") + 
    scale_color_discrete(name = "Cell type") +
    theme_bw()
  
  methods::show(fig)
  
}