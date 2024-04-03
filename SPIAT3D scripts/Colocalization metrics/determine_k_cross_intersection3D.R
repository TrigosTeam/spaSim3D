determine_k_cross_intersection3D <- function(k_cross_df) {
  
  k_cross_df$sign <- k_cross_df$Observed - k_cross_df$Expected
  
  ## Determine when sign flips from -ve to +ve OR +ve to -ve
  change_of_sign <- diff(sign(k_cross_df$sign))
  
  ## Determine indices for when observed curve goes above or below the expected curve
  observed_goes_above_indices <- which(change_of_sign == 2)
  observed_goes_below_indices <- which(change_of_sign == -2)
  
  if (length(indices) == 0) {
    stop("No cross-K intersections occur")
  }

  above_distance <- c()
  below_distance <- c()
  
  if (length(observed_goes_above_indices) != 0) {
    above_distance <- k_cross_df$Distance[observed_goes_above_indices]
    print("The observed curve goes above the expected curve at the following distances:")
    print(paste(above_distance))
  }
  if (length(observed_goes_below_indices) != 0) {
    below_distance <- k_cross_df$Distance[observed_goes_below_indices]
    print("The observed curve goes below the expected curve at the following distances:")
    print(paste(below_distance))
  }
  
  result <- data.frame(Distance = c(above_distance, below_distance),
                       Change = c(rep("Observed goes above Expected", length(above_distance)),
                                  rep("Observed goes below Expected", length(below_distance)))) 
  
  return (result)
}