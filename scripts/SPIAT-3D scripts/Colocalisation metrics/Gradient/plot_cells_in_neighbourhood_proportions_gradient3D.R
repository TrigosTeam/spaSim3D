plot_cells_in_neighbourhood_proportions_gradient3D <- function(cells_in_neighbourhood_proportions_gradient_df) {
  
  plot_result <- reshape2::melt(cells_in_neighbourhood_proportions_gradient_df, id.vars = c("radius"))
  fig <- ggplot(plot_result, aes(radius, value, color = variable)) +
    geom_point() +
    geom_line() +
    labs(title = "Neighbourhood cell proportion gradients", x = "Radius", y = "Cell proportion", color = "Cell type") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(0, 1)
  
  methods::show(fig)
  
}