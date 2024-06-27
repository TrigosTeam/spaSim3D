plot_cross_K_gradient_ratio3D <- function(cross_K_gradient_results) {
  
  plot(cross_K_gradient$radius, 
              cross_K_gradient$observed_cross_K / cross_K_gradient$expected_cross_K, 
              type = "l", 
              col = "red", 
              xlim = c(0, max(cross_K_gradient$radius)), ylim = c(0, 1.2 * max((cross_K_gradient$observed_cross_K / cross_K_gradient$expected_cross_K), 1)),
              xlab = "Radius", ylab = "Cross K-function ratio")
  abline(a = 1, b = 0, col = "blue", lty = 2)
  legend(0, 1.2 * max((cross_K_gradient$observed_cross_K / cross_K_gradient$expected_cross_K), 1), 
         legend = c("Observed cross K ratio", "Expected CSR cross K ratio"), col = c("red", "blue"), lty = c(1, 2))
  
}