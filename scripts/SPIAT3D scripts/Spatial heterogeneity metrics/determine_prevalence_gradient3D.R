determine_prevalence_gradient3D <- function(grid_data,
                                            metric_colname,
                                            show_AUC = T,
                                            plot_image = T) {
  
  # Thresholds range from 0 to 1
  thresholds <- seq(0.01, 1, 0.01)
  
  # Define result
  result <- data.frame(threshold = thresholds)
  
  prevalences <- c()
  
  for (threshold in thresholds) {
    prevalences <- c(prevalences, determine_prevalence3D(grid_data,
                                                         metric_colname,
                                                         threshold))
  }
  result$prevalence <- prevalences
  
  # Show AUC of prevalence gradient graph
  if (show_AUC) {
    print(paste("AUC:", round(calculate_prevalence_gradient_AUC3D(result), 2)))
  }
  
  # Plot
  if (plot_image) {
    fig <- ggplot(result, aes(threshold, prevalence)) +
      geom_line() +
      theme_bw() +
      labs(x = "Threshold",
           y = "Prevalence",
           title = paste("Prevalence vs Threshold (", metric_colname, ")", sep = "")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(0, 100)
    methods::show(fig)
  }

  return(result)
}
