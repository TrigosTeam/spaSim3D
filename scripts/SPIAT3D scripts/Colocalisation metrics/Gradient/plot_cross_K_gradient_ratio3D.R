plot_cross_K_gradient_ratio3D <- function(cross_K_gradient_results) {
  
  plot_result <- data.frame(radius = cross_K_gradient_results$radius,
                            observed_cross_K_gradient_ratio = cross_K_gradient_results$observed_cross_K / cross_K_gradient_results$expected_cross_K,
                            expected_cross_K_gradient_ratio = 1)
  
  plot_result <- reshape2::melt(plot_result, "radius", c("observed_cross_K_gradient_ratio", "expected_cross_K_gradient_ratio"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(x = "Radius", y = "Cross K-function ratio") +
    scale_colour_discrete(name = "", labels = c("Observed cross K ratio", "Expected CSR cross K ratio")) +
    theme_bw()
  
  methods::show(fig)
  
}