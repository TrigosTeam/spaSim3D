plot_cross_K_gradient3D <- function(cross_K_gradient_df) {
  
  plot_result <- reshape2::melt(cross_K_gradient_df, "radius", c("observed_cross_K", "expected_cross_K", "cross_K_ratio"))
  plot_result <- plot_result[plot_result$variable != "cross_K_ratio", ]
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(x = "Radius", y = "Cross K-function value") +
    scale_colour_discrete(name = "", labels = c("Observed cross K", "Expected CSR cross K")) +
    theme_bw()
  
  methods::show(fig)
  
}