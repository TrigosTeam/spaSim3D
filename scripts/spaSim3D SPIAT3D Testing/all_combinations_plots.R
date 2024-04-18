### Parameters to test

# Cluster type: Separate clusters (3), Ring (3), Mixing (3)
cluster_types <- c(paste("Separate", 1:3, sep = ""),
                   paste("Ring", 1:3, sep = ""),
                   paste("Mixing", 1:3, sep = ""))

# Shape: Sphere (S), Ellipsoid (E), Network (N)
separate_cluster_shapes <- c("SS", "EE", "NN", "SE", "ES", "SN", "NS", "EN", "NE")
ring_or_mixing_shapes <- c("S", "E", "N")

# Size: Small (s), Medium (m), Large ("l)
separate_cluster_sizes <- c("ss", "mm", "ll", "sm", "ms", "sl", "ls", "ml", "lm")
ring_or_mixing_sizes <- c("s", "m", "l")

# Cell Type: Tumour, Immune

### Set up data frame for plots_meta_data
df_colnames <- c("ClusterType", "Shape", "Size")
plots_meta_data <- data.frame(matrix(nrow = 297, ncol = length(df_colnames)))
colnames(plots_meta_data) <- df_colnames

i <- 1

for (cluster_type in cluster_types) {
  
  # Treat clusters and ring/mixing differently
  if (grepl("Separate", cluster_type)) {
    
    for (shape in separate_cluster_shapes) {
      
      for (size in separate_cluster_sizes) {
        
        plots_meta_data[i, ] <- c(cluster_type, shape, size)
        i <- i + 1
        
      }
    }
  }
  
  else if (grepl("Ring", cluster_type) || grepl("Mixing", cluster_type)) {
    
    for (shape in ring_or_mixing_shapes) {
      
      for (size in ring_or_mixing_sizes) {
        
        plots_meta_data[i, ] <- c(cluster_type, shape, size)
        i <- i + 1
        
      }
    }
  }
}


### Notes about position, size and shape of clusters
# Shape:          small,                medium,             large
# Sphere:         r = 15,               r = 20,             r = 25
# Ellipsoid:  rxyz = (10, 15, 20), rxyz = (15, 20, 25), rxyz = (20, 25, 30) 
# Network:     Cr = 35, r = 5,      CR = 40, r = 7.5,    CR = 45, r = 10     # CR = cluster radius

# Cluster type: Separate1,  Separate2,  Separate3,   Ring1,    Ring2,    Ring3,    Mixing1,    Mixing2,    Mixing3
#               ~20 apart   ~10 apart   ~0 apart    2 width   4 width   6 width   0.2 prop    0.4 prop    0.6 prop


shape_size_data <- list(
  S = list(
    s = 15,
    m = 20,
    l = 25
  ),
  E = list(
    s = c(10, 15, 20),
    m = c(15, 20, 25),
    l = c(20, 25, 30)
  ),
  N = list(
    s = c(35, 4),
    m = c(40, 7.5),
    l = c(45, 10)
  )
)

cluster_data <- list(
  Separate1 = list(
    cell1_centre = c(40, 40, 40),
    cell2_centre = c(110, 110, 110)
  ),
  Separate2 = list(
    cell1_centre = c(50, 50, 50),
    cell2_centre = c(100, 100, 100)
  ),
  Separate3 = list(
    cell1_centre = c(60, 60, 60),
    cell2_centre = c(90, 90, 90)
  ),
  Ring1 = 2,
  Ring2 = 4,
  Ring3 = 6,
  Mixing1 = 0.2,
  Mixing2 = 0.4,
  Mixing3 = 0.6
)



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



### Get plot data for each combination of parameters
cell_type1 <- "Tumour"
cell_type2 <- "Immune"
ring_or_mixing_centre <- c(75, 75, 75)

plot_data <- list()

for (i in 1:nrow(plots_meta_data)) {
  cluster_type <- plots_meta_data[i, "ClusterType"]
  shape <- plots_meta_data[i, "Shape"]
  size <- plots_meta_data[i, "Size"]
  
  if (grepl("Separate", cluster_type) == TRUE) {
    ## Get data for cluster1 and cluster2
    shape1 <- substr(shape, 1, 1)
    size1 <- substr(size, 1, 1)
    shape2 <- substr(shape, 2, 2)
    size2 <- substr(size, 2, 2)
    
    # Get cluster1 properties
    size_data <- shape_size_data[[shape1]][[size1]]
    centre <- cluster_data[[cluster_type]][["cell1_centre"]]
    if (shape1 != "N") centre <- adjust_centre(centre, size1, TRUE)
    cluster_properties <- get_cluster_properties(list(), shape1, cell_type1, size_data, centre)
  
    # Get cluster2 properties
    size_data <- shape_size_data[[shape2]][[size2]]
    centre <- cluster_data[[cluster_type]][["cell2_centre"]]
    if (shape2 != "N") centre <- adjust_centre(centre, size2, FALSE)
    cluster_properties <- get_cluster_properties(cluster_properties, shape2, cell_type2, size_data, centre)
    
    plot_data[[i]] <- simulate_clusters3D(bg,
                                          n_clusters = 2,
                                          cluster_properties = cluster_properties,
                                          plot_image = F)
    
  }
  
  else if (grepl("Ring", cluster_type) == TRUE) {
    size_data <- shape_size_data[[shape]][[size]]
    ring_width <- cluster_data[[cluster_type]]
    if (shape == "N") ring_width <- ring_width / 2
    
    ring_properties <- get_ring_properties(shape,
                                           cell_type1,
                                           size_data,
                                           ring_or_mixing_centre,
                                           cell_type2,
                                           ring_width)
    
    plot_data[[i]] <- simulate_rings3D(bg,
                                       n_ring = 1,
                                       ring_properties = ring_properties,
                                       plot_image = F)
  }
  
  else if (grepl("Mixing", cluster_type) == TRUE) {
    size_data <- shape_size_data[[shape]][[size]]
    infiltration_prop <- cluster_data[[cluster_type]]
    
    cluster_properties <- get_cluster_properties(list(), 
                                                 shape, 
                                                 cell_type1, 
                                                 size_data, 
                                                 ring_or_mixing_centre,
                                                 cell_type2,
                                                 infiltration_prop)
    plot_data[[i]] <- simulate_clusters3D(bg,
                                          n_clusters = 1,
                                          cluster_properties = cluster_properties,
                                          plot_image = F)
    
  }
  
}

setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/Code/spaSim 3D/objects")
all_plot_data <- readRDS(file="all_plot_test_data.rda")



get_cluster_properties <- function(cluster_properties, 
                                   shape, 
                                   cell_type, 
                                   size_data, 
                                   centre,
                                   infiltration_types = NULL,
                                   infiltration_proportions = NULL) {
  
  if (shape == "S") {
    cluster_properties[[length(cluster_properties) + 1]] <- list(
      shape = get_shape(shape),
      name_of_cluster_cell = cell_type,
      infiltration_types = infiltration_types,
      infiltration_proportions = infiltration_proportions,
      radius = size_data,
      centre_loc = centre
    )
  }
  else if (shape == "E") {
    cluster_properties[[length(cluster_properties) + 1]] <- list(
      shape = get_shape(shape),
      name_of_cluster_cell = cell_type,
      infiltration_types = infiltration_types,
      infiltration_proportions = infiltration_proportions,
      x_radius = size_data[1],
      y_radius = size_data[2],
      z_radius = size_data[3],
      y_z_rotation = 0,
      x_z_rotation = 0,
      x_y_rotation = 0,
      centre_loc = centre
    )
    
  }
  else if (shape == "N") {
    cluster_properties[[length(cluster_properties) + 1]] <- list(
      shape = get_shape(shape),
      n_edges = 12,
      name_of_cluster_cell = cell_type,
      infiltration_types = infiltration_types,
      infiltration_proportions = infiltration_proportions,
      radius = size_data[1],
      width = size_data[2],
      centre_loc = centre
    )
  }
  
  return (cluster_properties)
}



get_ring_properties <- function(shape, 
                                cluster_cell_type, 
                                size_data, 
                                centre,
                                ring_cell_type,
                                ring_width) {
 
  ring_properties <- list()
  
  if (shape == "S") {
    ring_properties[[1]] <- list(
      shape = get_shape(shape),
      name_of_cluster_cell = cluster_cell_type,
      infiltration_types = NULL,
      infiltration_proportions = NULL,
      radius = size_data,
      centre_loc = centre,
      name_of_ring_cell = ring_cell_type,
      ring_width = ring_width,
      ring_infiltration_types = NULL,
      ring_infiltration_proportions = NULL 
    )
  }
  else if (shape == "E") {
    ring_properties[[1]] <- list(
      shape = get_shape(shape),
      name_of_cluster_cell = cluster_cell_type,
      infiltration_types = NULL,
      infiltration_proportions = NULL,
      x_radius = size_data[1],
      y_radius = size_data[2],
      z_radius = size_data[3],
      y_z_rotation = 0,
      x_z_rotation = 0,
      x_y_rotation = 0,
      centre_loc = centre,
      name_of_ring_cell = ring_cell_type,
      ring_width = ring_width,
      ring_infiltration_types = NULL,
      ring_infiltration_proportions = NULL 
    )
    
  }
  else if (shape == "N") {
    ring_properties[[1]] <- list(
      shape = get_shape(shape),
      n_edges = 12,
      name_of_cluster_cell = cluster_cell_type,
      infiltration_types = NULL,
      infiltration_proportions = NULL,
      radius = size_data[1],
      width = size_data[2],
      centre_loc = centre,
      name_of_ring_cell = ring_cell_type,
      ring_width = ring_width,
      ring_infiltration_types = NULL,
      ring_infiltration_proportions = NULL 
    )
  }
  
  return (ring_properties) 
}



get_shape <- function(shape) {
  
  answer <- switch(
    shape,
    "S" = "Sphere",
    "E" = "Ellipsoid",
    "N" = "Network"
  )
  return (answer)
}


adjust_centre <- function(centre, size, increase) {
  
  # Increase centre
  if (increase == T) {
    if (size == "s") {
      centre <- centre + 8
    }
    else if (size == "m") {
      centre <- centre + 4
    }
    else if (size == "l") {
      # centre remains unchanged
    }
  }
  
  # Decrease centre
  else {
    if (size == "s") {
      centre <- centre - 8
    }
    else if (size == "m") {
      centre <- centre - 4
    }
    else if (size == "l") {
      # centre remains unchanged
    }
  }
  
  return (centre)
}
