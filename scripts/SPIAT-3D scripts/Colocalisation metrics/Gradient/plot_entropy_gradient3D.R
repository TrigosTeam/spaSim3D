plot_entropy_gradient3D <- function(entropy_gradient_df, expected_entropy = NULL) {
  
  plot_result <- entropy_gradient_df

  if (!is.null(expected_entropy)) {
    plot_result$expected_entropy <- expected_entropy
    plot_result <- reshape2::melt(plot_result, "radius", c("entropy", "expected_entropy"))
    
    fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
      geom_line() +
      labs(x = "Radius", y = "Entropy") +
      scale_colour_discrete(name = "", labels = c("Observed entropy", "Expected CSR entropy")) +
      theme_bw()
    
  }
  else {
    plot_result <- reshape2::melt(plot_result, "radius", c("entropy"))
    
    fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
      geom_line() +
      labs(x = "Radius", y = "Entropy") +
      scale_colour_discrete(name = "", labels = c("Observed entropy")) +
      theme_bw()
  }
    
  methods::show(fig)
  
}