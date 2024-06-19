calculate_AUC_of_Kcross3D <- function(Kcross_df) {
  
  ## Get difference in area between 
  AUC <- pracma::trapz(Kcross_df$Distance, Kcross_df$Observed) - 
         pracma::trapz(Kcross_df$Distance, Kcross_df$Expected)
  
  ## Get the cross-k result image size
  max_distance <- max(Kcross_df$Distance)
  max_Kcross <- max(c(Kcross_df$Observed, Kcross_df$Expected))
  
  ## Calculate normalised AUC
  n_AUC <- AUC / (max_distance * max_Kcross)
  
  return (n_AUC)
}
