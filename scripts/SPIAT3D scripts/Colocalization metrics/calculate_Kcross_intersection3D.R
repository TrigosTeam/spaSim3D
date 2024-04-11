calculate_Kcross_intersection3D <- function(Kcross_df) {
  
  Kcross_df$sign <- Kcross_df$Observed - Kcross_df$Expected
  
  ## Determine when sign flips from -ve to +ve OR +ve to -ve
  change_of_sign <- diff(sign(Kcross_df$sign))
  
  ## Determine indices for when observed curve goes above or below the expected curve
  observed_goes_above_indices <- which(change_of_sign == 2)
  observed_goes_below_indices <- which(change_of_sign == -2)
  
  if (length(observed_goes_above_indices) + length(observed_goes_below_indices) == 0) {
    warning("No cross-K intersections occur")
    return (0)
  }

  above_distance <- c()
  below_distance <- c()
  
  if (length(observed_goes_above_indices) != 0) {
    above_distance <- Kcross_df$Distance[observed_goes_above_indices]
    print("The observed curve goes ABOVE the expected curve at the following distances:")
    print(paste(above_distance))
  }
  if (length(observed_goes_below_indices) != 0) {
    below_distance <- Kcross_df$Distance[observed_goes_below_indices]
    print("The observed curve goes BELOW the expected curve at the following distances:")
    print(paste(below_distance))
  }
  
  result <- data.frame(Distance = c(above_distance, below_distance),
                       Change = c(rep("Observed goes above Expected", length(above_distance)),
                                  rep("Observed goes below Expected", length(below_distance)))) 
  
  return (result)
}
