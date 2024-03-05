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
  
  if (plot_image){
    if (is.null(plot_colours)){
      plot_colours <- c("gray","darkgreen", "red", "darkblue", "brown", "purple", "lightblue",
                        "lightgreen", "yellow", "black", "pink")
    }
    
    colors <- c()
    for (i in 1:nrow(bg_sample)) {
      for (j in 1:length(idents)) {
        if (bg_sample$Cell.Type[i] == idents[j]) {
          colors <- append(colors, plot_colours[j])
          break
        }
      }
    }
    
    plot3d(bg_sample$Cell.X.Position,
           bg_sample$Cell.Y.Position,
           bg_sample$Cell.Z.Position,
           xlab = "x",
           ylab = "y",
           zlab = "z",
           col = colors,
           size = 4)
    
    # add legend
    legend3d("topright", legend = idents, pch = 16, col = plot_colours[seq_len(length(idents))], inset = c(0.02))
    
  }
}
