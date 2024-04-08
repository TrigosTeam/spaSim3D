simulate_mixing3D <- function(bg_sample,
                              cell_types = c("Others", "Immune", "Tumour"),
                              props = c(0.5, 0.2, 0.3),
                              plot_image = TRUE,
                              plot_categories = c("Others", "Immune", "Tumour"),
                              plot_colours = c("lightgray", "skyblue", "orange")) {
  
  
  n_types <- length(cell_types)
  
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
        pheno <- cell_types[n]
        break
      }
      n <- n+1
    }
    bg_sample[i, "Cell.Type"] <- pheno
  }
  
  if (plot_image) {
    
    if (is.null(plot_categories)) {
      plot_categories <- unique(bg_sample$Cell.Type)
    }
    
    if (is.null(plot_colours)) {
      plot_colours <- c("red", "orange", "green", "blue", "skyblue", "pink", "purple", "lightgray")[1:length(plot_categories)]
    }
    
    bg_sample <- bg_sample[bg_sample[["Cell.Type"]] %in% plot_categories, ]
    
    ## Factor for feature column
    bg_sample[, "Cell.Type"] <- factor(bg_sample[, "Cell.Type"],
                                  levels = plot_categories)
    
    
    ## Plot
    fig <- plot_ly(bg_sample,
                   type = "scatter3d",
                   mode = 'markers',
                   x = ~Cell.X.Position,
                   y = ~Cell.Y.Position,
                   z = ~Cell.Z.Position,
                   color = ~Cell.Type,
                   colors = plot_colours,
                   marker = list(size = 2))
    
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                       yaxis = list(title = 'y'),
                                       zaxis = list(title = 'z')))
    print(fig)
  }
  
  return(bg_sample)
}
