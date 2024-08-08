## 1.1. Generate the different background cell arrangements metadata for mixed and ringed arrangements -----------------
# 100 x 100 x 100, 20000 cells

bg_meta <- spe_metadata_background_template("random")
bg_meta$background$n_cells <- 20000
bg_meta$background$length <- 100
bg_meta$background$width <- 100
bg_meta$background$height <- 100
bg_meta$background$minimum_distance_between_cells <- 0

cell_types <- c("O", "A", "B")
cell_props <- data.frame(O = c(1, 0, 0),
                         A = c(0.95, 0.05, 0),
                         B = c(0.95, 0, 0.05),
                         AB = c(0.9, 0.05, 0.05))
bg_types <- c("O", "A", "B", "AB")

bg_spes_metadata <- list()
for (i in seq(ncol(cell_props))) {
  
  bg_meta$background$cell_types <- cell_types
  bg_meta$background$cell_proportions <- cell_props[[bg_types[i]]]
  
  bg_spes_metadata[[bg_types[i]]] <- bg_meta
}

# setwd("~/Objects/supervised/spes_metadata")
# saveRDS(bg_spes_metadata, "bg_spes_metadata_100_100_100.rds")

## 1.2. Generate the different background cell arrangements metadata for separate cluster arrangements -----------------
# 200 x 100 x 100, 40000 cells

bg_meta <- spe_metadata_background_template("random")
bg_meta$background$n_cells <- 30000
bg_meta$background$length <- 150
bg_meta$background$width <- 100
bg_meta$background$height <- 100
bg_meta$background$minimum_distance_between_cells <- 0

cell_types <- c("O", "A", "B")
cell_props <- data.frame(O = c(1, 0, 0),
                         A = c(0.95, 0.05, 0),
                         B = c(0.95, 0, 0.05),
                         AB = c(0.9, 0.05, 0.05))
bg_types <- c("O", "A", "B", "AB")

bg_spes_metadata <- list()
for (i in seq(ncol(cell_props))) {
  
  bg_meta$background$cell_types <- cell_types
  bg_meta$background$cell_proportions <- cell_props[[bg_types[i]]]
  
  bg_spes_metadata[[bg_types[i]]] <- bg_meta
}

# setwd("~/Objects/supervised/spes_metadata")
# saveRDS(bg_spes_metadata, "bg_spes_metadata_150_100_100.rds")

## 2.1. Generate table for mixed clusters ------------------------------------------

# Determine possible values for different parameters
bg_types <- c("O", "A", "B", "AB")
shapes <- c("Sphere", "Ellipsoid", "Network")
sizes <- c("Small", "Medium", "Large")
arrangements <- c("M1", "M2", "M3")


# Set up the mixed_spes_table
df <- data.frame(matrix(nrow = length(bg_types) * length(shapes) * length(sizes) * length(arrangements),
                        ncol = 4))
colnames(df) <- c("bg_type", "shape", "size", "arrangement")


# Fill in the mixed_spes_table
i <- 1
for (bg_type in bg_types) {
  for (shape in shapes) {
    for (size in sizes) {
      for (arrangement in arrangements) {
        
        df[i, "bg_type"] <- bg_type
        df[i, "shape"] <- shape
        df[i, "size"] <- size
        df[i, "arrangement"] <- arrangement
        i <- i + 1
      }
    }
  }
}

# setwd("~/Objects/supervised/spes_table")
# write.table(mixed_spes_table, "mixed_spes_table.csv")

### 2.2. Generate metadata for mixed clusters -----------------------

# Define shape parameters based on the size (Small, Medium, Large)
shape_parameters_list <- list()
shape_parameters_list[["Sphere"]] <- data.frame(Small = c(20), Medium = c(30), Large = c(40))
rownames(shape_parameters_list[["Sphere"]]) <- c("radius")
shape_parameters_list[["Ellipsoid"]] <- data.frame(Small = c(20, 25, 15), Medium = c(30, 35, 25), Large = c(40, 45, 35))
rownames(shape_parameters_list[["Ellipsoid"]]) <- c("x_radius", "y_radius", "z_radius")
shape_parameters_list[["Network"]] <- data.frame(Small = c(5), Medium = c(7.5), Large = c(10))
rownames(shape_parameters_list[["Network"]]) <- c("width")


# Define constant shape parameters (constant regardless of size)
centre_loc <- c(50, 50, 50)
network_radius <- 40
network_n_edges <- 15
cell_types <- c("A", "B")


### Generate metadata for mixed simulations

# Choose mixing cell proportions for each arrangement
mixing_props <- data.frame(A = c(0.9, 0.7, 0.5),
                           B = c(0.1, 0.3, 0.5))
rownames(mixing_props) <- c("M1", "M2", "M3")

# Get background data
setwd("~/Objects/supervised/spes_metadata")
bg_spes_metadata <- readRDS("bg_spes_metadata_100_100_100.rds")

# Get table for mixed simulations
setwd("~/Objects/supervised/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

mixed_spes_metadata <- list()

for (i in seq(nrow(mixed_spes_table))) {
  # Get general simulation parameters
  bg_type <- mixed_spes_table$bg_type[i]
  shape <- mixed_spes_table$shape[i]
  size <- mixed_spes_table$size[i]
  arrangement <- mixed_spes_table$arrangement[i]
  
  
  # Get bg_spe_metadata from bg_type
  bg_spe_metadata <- bg_spes_metadata[[bg_type]]
  
  # Get cluster mixing cell proportions from arrangement
  cell_props <- as.numeric(mixing_props[arrangement, ])
  
  # Get metadata template for current simulation parameters
  curr_metadata <- spe_metadata_cluster_template(bg_spe_metadata, "regular", shape)
  curr_metadata$cluster_1$cluster_cell_types <- cell_types
  curr_metadata$cluster_1$cluster_cell_proportions <- cell_props
  curr_metadata$cluster_1$centre_loc <- centre_loc
  
  # Get specific shape-size parameters
  shape_parameters <- shape_parameters_list[[shape]]
  
  # Specify metadata for each shape and size
  if (shape == "Sphere") {
    curr_metadata$cluster_1$radius <- shape_parameters["radius", size]
  }
  else if (shape == "Ellipsoid") {
    curr_metadata$cluster_1$x_radius <- shape_parameters["x_radius", size]
    curr_metadata$cluster_1$y_radius <- shape_parameters["y_radius", size]
    curr_metadata$cluster_1$z_radius <- shape_parameters["z_radius", size]
    curr_metadata$cluster_1$x_y_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_1$x_z_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_1$y_z_rotation <- runif(1, 0, 180)
    
  }
  else if (shape == "Network") {
    curr_metadata$cluster_1$n_edges <- network_n_edges
    curr_metadata$cluster_1$width <- shape_parameters["width", size]
    curr_metadata$cluster_1$radius <- network_radius
  }
  else {
    stop(paste(shape, "Shape not found"))
  }
  
  mixed_spes_metadata[[i]] <- curr_metadata
}

setwd("~/Objects/supervised/spes_metadata")
# saveRDS(mixed_spes_metadata, "mixed_spes_metadata_supervised.rds")

## 2.3. Generate table for ringed clusters ------------------------------------------

# Determine possible values for different parameters
bg_types <- c("O", "A", "B", "AB")
shapes <- c("Sphere", "Ellipsoid", "Network")
sizes <- c("Small", "Medium", "Large")
arrangements <- c("R1", "R2", "R3")


# Set up the ringed_spes_table
df <- data.frame(matrix(nrow = length(bg_types) * length(shapes) * length(sizes) * length(arrangements),
                        ncol = 4))
colnames(df) <- c("bg_type", "shape", "size", "arrangement")


# Fill in the ringed_spes_table
i <- 1
for (bg_type in bg_types) {
  for (shape in shapes) {
    for (size in sizes) {
      for (arrangement in arrangements) {
        
        df[i, "bg_type"] <- bg_type
        df[i, "shape"] <- shape
        df[i, "size"] <- size
        df[i, "arrangement"] <- arrangement
        i <- i + 1
      }
    }
  }
}

setwd("~/Objects/supervised/spes_table")
# write.table(ringed_spes_table, "ringed_spes_table.csv")


## 2.4. Generate metadata for ringed clusters -------------

# Define shape parameters based on the size (Small, Medium, Large)
shape_parameters_list <- list()
shape_parameters_list[["Sphere"]] <- data.frame(Small = c(20), Medium = c(30), Large = c(40))
rownames(shape_parameters_list[["Sphere"]]) <- c("radius")
shape_parameters_list[["Ellipsoid"]] <- data.frame(Small = c(20, 25, 15), Medium = c(30, 35, 25), Large = c(40, 45, 35))
rownames(shape_parameters_list[["Ellipsoid"]]) <- c("x_radius", "y_radius", "z_radius")
shape_parameters_list[["Network"]] <- data.frame(Small = c(5), Medium = c(7.5), Large = c(10))
rownames(shape_parameters_list[["Network"]]) <- c("width")


# Choose ring width for each arrangement and shape
ring_width_factors <- c("R1" = 0.1, "R2" = 0.15, "R3" = 0.2)

# Define constant shape parameters (constant regardless of size)
centre_loc <- c(50, 50, 50)
network_radius <- 40
network_n_edges <- 15
cell_types <- c("A")
cell_props <- c(1)
ring_cell_types <- c("B")
ring_cell_props <- c(1)

### Generate metadata for ringed simulations

# Get background data
setwd("~/Objects/supervised/spes_metadata")
bg_spes_metadata <- readRDS("bg_spes_metadata_100_100_100.rds")

# Get table for ringed simulations
setwd("~/Objects/supervised/spes_table")
ringed_spes_table <- read.table("ringed_spes_table.csv")

ringed_spes_metadata <- list()

for (i in seq(nrow(ringed_spes_table))) {
  # Get general simulation parameters
  bg_type <- ringed_spes_table$bg_type[i]
  shape <- ringed_spes_table$shape[i]
  size <- ringed_spes_table$size[i]
  arrangement <- ringed_spes_table$arrangement[i]
  
  
  # Get bg_spe_metadata from bg_type
  bg_spe_metadata <- bg_spes_metadata[[bg_type]]
  
  # Get metadata template for current simulation parameters
  curr_metadata <- spe_metadata_cluster_template(bg_spe_metadata, "ring", shape)
  curr_metadata$cluster_1$cluster_cell_types <- cell_types
  curr_metadata$cluster_1$cluster_cell_proportions <- cell_props
  curr_metadata$cluster_1$centre_loc <- centre_loc
  curr_metadata$cluster_1$ring_cell_types <- ring_cell_types
  curr_metadata$cluster_1$ring_cell_proportions <- ring_cell_props
  
  # Get specific shape-size parameters
  shape_parameters <- shape_parameters_list[[shape]]
  
  # Specify metadata for each shape and size
  if (shape == "Sphere") {
    curr_metadata$cluster_1$radius <- shape_parameters["radius", size]
    curr_metadata$cluster_1$ring_width <- ring_width_factors[arrangement] * shape_parameters["radius", size]
  }
  else if (shape == "Ellipsoid") {
    curr_metadata$cluster_1$x_radius <- shape_parameters["x_radius", size]
    curr_metadata$cluster_1$y_radius <- shape_parameters["y_radius", size]
    curr_metadata$cluster_1$z_radius <- shape_parameters["z_radius", size]
    curr_metadata$cluster_1$x_y_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_1$x_z_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_1$y_z_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_1$ring_width <- ring_width_factors[arrangement] * shape_parameters["x_radius", size]
    
  }
  else if (shape == "Network") {
    curr_metadata$cluster_1$n_edges <- network_n_edges
    curr_metadata$cluster_1$width <- shape_parameters["width", size]
    curr_metadata$cluster_1$radius <- network_radius
    curr_metadata$cluster_1$ring_width <- ring_width_factors[arrangement] * shape_parameters["width", size]
  }
  else {
    stop(paste(shape, "Shape not found"))
  }
  
  ringed_spes_metadata[[i]] <- curr_metadata
}


setwd("~/Objects/supervised/spes_metadata")
# saveRDS(ringed_spes_metadata, "ringed_spes_metadata_supervised.rds")


## 2.5. Generate table for separated clusters ----------------------------------
# Determine possible values for different parameters
bg_types <- c("O", "A", "B", "AB")
shapes <- c("Sphere", "Ellipsoid", "Network")
sizes <- c("Small", "Medium", "Large")
arrangements <- c("S1", "S2", "S3")


# Set up the separated_spes_table
df <- data.frame(matrix(nrow = length(bg_types) * length(shapes)^2 * length(sizes)^2 * length(arrangements),
                        ncol = 6))
colnames(df) <- c("bg_type", "shapeA", "sizeA", "shapeB", "sizeB", "arrangement")


# Fill in the separated_spes_table
i <- 1
for (bg_type in bg_types) {
  for (shapeA in shapes) {
    for (sizeA in sizes) {
      for (shapeB in shapes) {
        for (sizeB in sizes) {
          for (arrangement in arrangements) {
            
            df[i, "bg_type"] <- bg_type
            df[i, "shapeA"] <- shapeA
            df[i, "sizeA"] <- sizeA
            df[i, "shapeB"] <- shapeB
            df[i, "sizeB"] <- sizeB
            df[i, "arrangement"] <- arrangement
            i <- i + 1
          }   
        } 
      }
    }
  }
}

setwd("~/Objects/supervised/spes_table")
# write.table(separated_spes_table, "separated_spes_table.csv")


## 2.6. Generate metadata for separated clusters ----------------------------------
# Define shape parameters based on the size (Small, Medium, Large)
shape_parameters_list <- list()
shape_parameters_list[["Sphere"]] <- data.frame(Small = c(20), Medium = c(30), Large = c(40))
rownames(shape_parameters_list[["Sphere"]]) <- c("radius")
shape_parameters_list[["Ellipsoid"]] <- data.frame(Small = c(20, 25, 15), Medium = c(30, 35, 25), Large = c(40, 45, 35))
rownames(shape_parameters_list[["Ellipsoid"]]) <- c("x_radius", "y_radius", "z_radius")
shape_parameters_list[["Network"]] <- data.frame(Small = c(5), Medium = c(7.5), Large = c(10))
rownames(shape_parameters_list[["Network"]]) <- c("width")


# Choose cluster centre location for each arrangement
centre_x_loc <- data.frame(clusterA = c(25, 37.5, 50),
                           clusterB = c(125, 112.5, 100))
rownames(centre_x_loc) <- c("S1", "S2", "S3")

# Define constant shape parameters (constant regardless of size)
network_n_edges <- 15
network_radius <- 40
clusterA_cell_types <- c("A")
clusterA_cell_props <- c(1)
clusterB_cell_types <- c("B")
clusterB_cell_props <- c(1)

### Generate metadata for separated simulations

# Get background data
setwd("~/Objects/supervised/spes_metadata")
bg_spes_metadata <- readRDS("bg_spes_metadata_150_100_100.rds")

# Get table for ringed simulations
setwd("~/Objects/supervised/spes_table")
separated_spes_table <- read.table("separated_spes_table.csv")

separated_spes_metadata <- list()

for (i in seq(nrow(separated_spes_table))) {
  # Get general simulation parameters
  bg_type <- separated_spes_table$bg_type[i]
  shapeA <- separated_spes_table$shapeA[i]
  sizeA <- separated_spes_table$sizeA[i]
  shapeB <- separated_spes_table$shapeB[i]
  sizeB <- separated_spes_table$sizeB[i]
  arrangement <- separated_spes_table$arrangement[i]
  
  # Get bg_spe_metadata from bg_type
  bg_spe_metadata <- bg_spes_metadata[[bg_type]]
  
  # Get centre location of clusters from arrangement
  centre_locA <- c(centre_x_loc[arrangement, "clusterA"], 50, 50)
  centre_locB <- c(centre_x_loc[arrangement, "clusterB"], 50, 50)
  
  # Get metadata template for current simulation parameters for both clusters
  curr_metadata <- spe_metadata_cluster_template(bg_spe_metadata, "regular", shapeA)
  curr_metadata$cluster_1$cluster_cell_types <- clusterA_cell_types
  curr_metadata$cluster_1$cluster_cell_proportions <- clusterA_cell_props
  curr_metadata$cluster_1$centre_loc <- centre_locA
  
  curr_metadata <- spe_metadata_cluster_template(curr_metadata, "regular", shapeB)
  curr_metadata$cluster_2$cluster_cell_types <- clusterB_cell_types
  curr_metadata$cluster_2$cluster_cell_proportions <- clusterB_cell_props
  curr_metadata$cluster_2$centre_loc <- centre_locB
  
  
  # Get specific shape-size parameters
  shape_parametersA <- shape_parameters_list[[shapeA]]
  shape_parametersB <- shape_parameters_list[[shapeB]]
  
  # Specify metadata for each shape and size
  if (shapeA == "Sphere") {
    curr_metadata$cluster_1$radius <- shape_parametersA["radius", sizeA]
  }
  else if (shapeA == "Ellipsoid") {
    curr_metadata$cluster_1$x_radius <- shape_parametersA["x_radius", sizeA]
    curr_metadata$cluster_1$y_radius <- shape_parametersA["y_radius", sizeA]
    curr_metadata$cluster_1$z_radius <- shape_parametersA["z_radius", sizeA]
    curr_metadata$cluster_1$x_y_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_1$x_z_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_1$y_z_rotation <- runif(1, 0, 180)
  }
  else if (shapeA == "Network") {
    curr_metadata$cluster_1$n_edges <- network_n_edges
    curr_metadata$cluster_1$width <- shape_parametersA["width", sizeA]
    curr_metadata$cluster_1$radius <- network_radius
  }
  else {
    stop(paste(shapeA, "ShapeA not found"))
  }
  
  if (shapeB == "Sphere") {
    curr_metadata$cluster_2$radius <- shape_parametersB["radius", sizeB]
  }
  else if (shapeB == "Ellipsoid") {
    curr_metadata$cluster_2$x_radius <- shape_parametersB["x_radius", sizeB]
    curr_metadata$cluster_2$y_radius <- shape_parametersB["y_radius", sizeB]
    curr_metadata$cluster_2$z_radius <- shape_parametersB["z_radius", sizeB]
    curr_metadata$cluster_2$x_y_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_2$x_z_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_2$y_z_rotation <- runif(1, 0, 180)
  }
  else if (shapeB == "Network") {
    curr_metadata$cluster_2$n_edges <- network_n_edges
    curr_metadata$cluster_2$width <- shape_parametersB["width", sizeB]
    curr_metadata$cluster_2$radius <- network_radius
  }
  else {
    stop(paste(shapeB, "ShapeB not found"))
  }
  separated_spes_metadata[[i]] <- curr_metadata
}
setwd("~/Objects/supervised/spes_metadata")
saveRDS(separated_spes_metadata, "separated_spes_metadata_supervised.rds")
