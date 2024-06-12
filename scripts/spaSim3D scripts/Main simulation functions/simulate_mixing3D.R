simulate_mixing3D <- function(bg_sample,
                              cell_types = c("Others", "Immune", "Tumour"),
                              cell_proportions = c(0.5, 0.2, 0.3),
                              plot_image = TRUE,
                              plot_categories = c("Others", "Immune", "Tumour"),
                              plot_colours = c("lightgray", "skyblue", "orange")) {
  
  
  n_cell_types <- length(cell_types)
  
  for (i in 1:nrow(bg_sample)) {
    x <- bg_sample$Cell.X.Position[i]
    y <- bg_sample$Cell.Y.Position[i]
    z <- bg_sample$Cell.Z.Position[i]
    
    # Random number will determine the cell_type of the cell
    random <- runif(n = 1, min = 0, max = 1)
    
    # Start with the first cell
    n <- 1 
    current_proportion <- 0
    
    while (n <= n_cell_types){
      current_proportion <- current_proportion + cell_proportions[n]
      if (random <= current_proportion) {
        chosen_cell_type <- cell_types[n]
        break
      }
      n <- n + 1
    }
    bg_sample[i, "Cell.Type"] <- chosen_cell_type
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cell_categories3D(data = bg_sample,
                                  cell_types_of_interest = plot_categories,
                                  colour_vector = plot_colours,
                                  size = 2,
                                  include_cell_types_of_no_interest = FALSE,
                                  feature_colname = "Cell.Type")
    print(fig)
  }
    
  return(bg_sample)
}
