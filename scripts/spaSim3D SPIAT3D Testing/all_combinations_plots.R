### Parameters to test

# Position: Separate clusters (3), Ring (3), Mixing (3)
positions <- c(paste("Cluster", 1:3, sep = ""),
               paste("Ring", 1:3, sep = ""),
               paste("Mixing", 1:3, sep = ""))

# Shape: Sphere (S), Ellipsoid (E), Cylinder (C)
separate_cluster_shapes <- c("SS", "EE", "CC", "SE", "ES", "SC", "CS", "EC", "CE")
ring_or_mixing_shapes <- c("S", "E", "C")

# Size: Small (s), Medium (m), Large ("l)
separate_cluster_sizes <- c("ss", "mm", "ll", "sm", "ms", "sl", "ls", "ml", "lm")
ring_or_mixing_sizes <- c("s", "m", "l")

# Cell Type: Tumour, Immune

### Set up data frame for plots_meta_data
df_colnames <- c("Position", "Shape", "Size")
plots_meta_data <- data.frame(matrix(nrow = 297, ncol = length(df_colnames)))
colnames(plots_meta_data) <- df_colnames

i <- 1

for (position in positions) {
  
  # Treat clusters and ring/mixing differently
  if (grepl("Cluster", position)) {
    
    for (shape in separate_cluster_shapes) {
      
      for (size in separate_cluster_sizes) {
        
        plots_meta_data[i, ] <- c(position, shape, size)
        i <- i + 1
        
      }
    }
  }
  
  else if (grepl("Ring", position) || grepl("Mixing", position)) {
    
    for (shape in ring_or_mixing_shapes) {
      
      for (size in ring_or_mixing_sizes) {
        
        plots_meta_data[i, ] <- c(position, shape, size)
        i <- i + 1
        
      }
    }
  }
}


### Notes about position, size and shape of clusters
# Shape:          small,                medium,             large
# Sphere:         r = 15,               r = 20,             r = 25
# Ellipsoid:  rxyz = (10, 15, 20), rxyz = (15, 20, 25), rxyz = (20, 25, 30) 
# Cylinder:     Cr = 15, r = 10,      CR = 20, r = 15,    CR = 25, r = 20     # CR = cluster radius

# Position: Cluster1,   Cluster2,   Cluster3,   Ring1,    Ring2,    Ring3,    Mixing1,    Mixing2,    Mixing3
#           ~20 apart   ~10 apart   ~0 apart    4 width   6 width   8 width   0.2 prop    0.4 prop    0.6 prop





### Window size
n_cells <- 35000
length <- 150
width <- 150
height <- 150

bg <- simulate_background_cells3D(n_cells = n_cells,
                                  length = length,
                                  width  = width,
                                  height = height,
                                  method = "tumour",
                                  min_d = 2,
                                  oversampling_rate = 1.2,
                                  jitter_prop = 0,
                                  cell_type = "Others",
                                  plot_image = T)


