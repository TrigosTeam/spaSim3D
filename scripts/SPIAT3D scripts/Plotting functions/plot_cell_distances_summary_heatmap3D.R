plot_cell_distances_summary_heatmap3D <- function(distances_summary_df, 
                                                  metric = "Mean"){
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  Reference <- Target <- Min <- Max <- Mean <- Std.Dev <- Median <- NULL
  
  if (metric == "Mean") {
    limit <- range(unlist(distances_summary_df$Mean), na.rm=TRUE)
    g <- ggplot(distances_summary_df, aes(x = Reference, y = Target, fill = Mean))
    
  }
  else if (metric == "Std.Dev") {
    limit <- range(unlist(distances_summary_df$Std.Dev), na.rm=TRUE)
    g <- ggplot(distances_summary_df, aes(x = Reference, y = Target, fill = Std.Dev))
    
  }
  else if (metric == "Median") {
    limit <- range(unlist(distances_summary_df$Median), na.rm=TRUE)
    g <- ggplot(distances_summary_df, aes(x = Reference, y = Target, fill = Median))
  }
  else if (metric == "Min") {
    limit <- range(unlist(distances_summary_df$Min), na.rm=TRUE)
    g <- ggplot(distances_summary_df, aes(x = Reference, y = Target, fill = Min))
  }
  else if (metric == "Max") {
    limit <- range(unlist(distances_summary_df$Max), na.rm=TRUE)
    g <- ggplot(distances_summary_df, aes(x = Reference, y = Target, fill = Max))
  }
  else {
    stop(paste(metric," is not a valid metric"))
  }
  
  g <- g +
    geom_tile() +
    xlab("Reference cell type") +
    ylab("Target cell type") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_viridis_c(limits = limit, direction = -1)
  
  return (g)
}
