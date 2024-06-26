calculate_minimum_distances_to_alpha_hull3D <- function(spe_with_alpha_hull, cell_types_of_interest, feature_colname = "Cell.Type", plot_image = T) {
  
  ## Get alpha hull numbers (ignoring -1)
  alpha_hull_numbers <- spe_alpha_hull$alpha_hull_number[spe_alpha_hull$alpha_hull_number != -1]
  
  ## Get number of alpha hulls
  n_alpha_hulls <- length(unique(alpha_hull_numbers))
  
  ## For each alpha hull, determine the minimum distance of each cell_type_of_interest
  result <- data.frame()
  
  
  spe_coords <- spatialCoords(spe_with_alpha_hull)
  cell_types_of_interest_coords <- list()
  for (cell_type in cell_types_of_interest) {
    cell_types_of_interest_coords[[cell_type]] <- spe_coords[spe_with_alpha_hull[[feature_colname]] == cell_type, ]
  }
  
    
  
  result <- vector()
  
  for (i in seq(n_alpha_hulls)) {
    alpha_hull_coords <- spe_coords[spe_with_alpha_hull$alpha_hull_number == i, ]
    
    for (cell_type in cell_types_of_interest) {
      curr_cell_type_coords <- cell_types_of_interest_coords[[cell_type]]
      
      all_closest <- RANN::nn2(data = alpha_hull_coords, 
                               query = curr_cell_type_coords, 
                               k = 1)  
      
      local_dist_mins <- data.frame(
        alpha_hull_number = i,
        cell_type_of_interest = cell_type,
        distance = all_closest$nn.dists
      )
      ## Remove any 0 distance rows
      local_dist_mins <- local_dist_mins[local_dist_mins$distance != 0, ]
      result <- rbind(result, local_dist_mins)
    }
    

    ## Plot
    if (plot_image) {
      
      fig <- ggplot(result, aes(x = cell_type_of_interest, y = distance, fill = cell_type_of_interest)) + 
        geom_violin() +
        facet_grid(alpha_hull_number~., scales="free_x") +
        theme_bw() +
        theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none") +
        labs(title="Minimum cell distances to alpha hulls", x = "Cell type", y = "Distance") +
        stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
      
      methods::show(fig)
    }
    
  }
  return(result)
}