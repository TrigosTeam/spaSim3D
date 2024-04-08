library(rgl)


simulate_mixing3D <- function(bg_sample,
                              idents = c("Others", "Immune", "Tumour"),
                              props = c(0.5, 0.2, 0.3),
                              plot_image = TRUE,
                              plot_colours = NULL) {
  
  
  n_types <- length(idents)
  
  for (i in 1:nrow(bg_sample)) {
    x <- bg_sample$Cell.X.Position[i]
    y <- bg_sample$Cell.Y.Position[i]
    z <- bg_sample$Cell.Z.Position[i]
    
    random <- runif(1)
    
    # if the random number falls in the range of an infiltration proportion,
    # pheno will be the corresponding infiltraiton type
    n <- 1 # start from the first proportion
    current_p <- 0
    while (n <= n_types){
      current_p <- current_p + props[n]
      if (random <= current_p) {
        pheno <- idents[n]
        break
      }
      n <- n+1
    }
    bg_sample[i, "Cell.Type"] <- pheno
  }
  
  if (plot_image) {
    plot <- plot_cell_categories3D(bg_sample,
                                   cell_types_of_interest = plot_categories,
                                   colour_vector = plot_colours)
    print(plot)
    
  }
  
  return(bg_sample)
}
