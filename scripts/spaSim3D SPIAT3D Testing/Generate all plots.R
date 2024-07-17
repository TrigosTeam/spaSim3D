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

## 1.2. Generate the different background cell arrangements metadata for separate cluster arrangements -----------------
# 200 x 100 x 100, 40000 cells

bg_meta <- spe_metadata_background_template("random")
bg_meta$background$n_cells <- 40000
bg_meta$background$length <- 200
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

## 2.1. Generate metadata for mixed clusters ------------------------------------------

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

setwd("~/Objects/spes_table")
# write.table(mixed_spes_table, "mixed_spes_table.csv")


# Define shape parameters based on the size (Small, Medium, Large)
shape_parameters_list <- list()
shape_parameters_list[["Sphere"]] <- data.frame(Small = c(20), Medium = c(30), Large = c(40))
rownames(shape_parameters_list[["Sphere"]]) <- c("radius")
shape_parameters_list[["Ellipsoid"]] <- data.frame(Small = c(25, 20, 15), Medium = c(37.5, 30, 22.5), Large = c(50, 40, 30))
rownames(shape_parameters_list[["Ellipsoid"]]) <- c("x_radius", "y_radius", "z_radius")
shape_parameters_list[["Network"]] <- data.frame(Small = c(3, 30), Medium = c(6, 40), Large = c(9, 50))
rownames(shape_parameters_list[["Network"]]) <- c("width", "radius")


# Define constant shape parameters (constant regardless of size)
centre_loc <- c(50, 50, 50)
ellipsoid_x_y_rotation <- 0
ellipsoid_x_z_rotation <- 0
ellipsoid_y_z_rotation <- 0
network_n_edges <- 15
cell_types <- c("A", "B")


### Generate metadata for mixed simulations

# Choose mixing cell proportions for each arrangement
mixing_props <- data.frame(A = c(0.9, 0.7, 0.5),
                           B = c(0.1, 0.3, 0.5))
rownames(mixing_props) <- c("M1", "M2", "M3")

# Get background data
setwd("~/Objects/spes_metadata")
bg_spes_metadata <- readRDS("bg_spes_metadata_100_100_100.rds")

# Get table for mixed simulations
setwd("~/Objects/spes_table")
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
    curr_metadata$cluster_1$x_y_rotation <- ellipsoid_x_y_rotation
    curr_metadata$cluster_1$x_z_rotation <- ellipsoid_x_z_rotation
    curr_metadata$cluster_1$y_z_rotation <- ellipsoid_y_z_rotation
    
  }
  else if (shape == "Network") {
    curr_metadata$cluster_1$n_edges <- network_n_edges
    curr_metadata$cluster_1$width <- shape_parameters["width", size]
    curr_metadata$cluster_1$radius <- shape_parameters["radius", size]
  }
  else {
    stop(paste(shape, "Shape not found"))
  }
  
  mixed_spes_metadata[[i]] <- curr_metadata
}


## 2.2. Generate metadata for ringed clusters ------------------------------------------

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

setwd("~/Objects/spes_table")
# write.table(ringed_spes_table, "ringed_spes_table.csv")



# Define shape parameters based on the size (Small, Medium, Large)
shape_parameters_list <- list()
shape_parameters_list[["Sphere"]] <- data.frame(Small = c(20), Medium = c(30), Large = c(40))
rownames(shape_parameters_list[["Sphere"]]) <- c("radius")
shape_parameters_list[["Ellipsoid"]] <- data.frame(Small = c(25, 20, 15), Medium = c(37.5, 30, 22.5), Large = c(50, 40, 30))
rownames(shape_parameters_list[["Ellipsoid"]]) <- c("x_radius", "y_radius", "z_radius")
shape_parameters_list[["Network"]] <- data.frame(Small = c(3, 30), Medium = c(6, 40), Large = c(9, 50))
rownames(shape_parameters_list[["Network"]]) <- c("width", "radius")


# Choose ring width for each arrangement and shape
ring_widths <- data.frame(Sphere = c(3, 6, 9),
                          Ellipsoid = c(3, 6, 9),
                          Network = c(1, 2, 3))
rownames(ring_widths) <- c("R1", "R2", "R3")

# Define constant shape parameters (constant regardless of size)
centre_loc <- c(50, 50, 50)
ellipsoid_x_y_rotation <- 0
ellipsoid_x_z_rotation <- 0
ellipsoid_y_z_rotation <- 0
network_n_edges <- 15
cell_types <- c("A")
cell_props <- c(1)
ring_cell_types <- c("B")
ring_cell_props <- c(1)

### Generate metadata for ringed simulations

# Get background data
setwd("~/Objects/spes_metadata")
bg_spes_metadata <- readRDS("bg_spes_metadata_100_100_100.rds")

# Get table for ringed simulations
setwd("~/Objects/spes_table")
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
  
  # Get cluster mixing cell proportions from arrangement and shape
  ring_width <- ring_widths[arrangement, shape]
  
  # Get metadata template for current simulation parameters
  curr_metadata <- spe_metadata_cluster_template(bg_spe_metadata, "ring", shape)
  curr_metadata$cluster_1$cluster_cell_types <- cell_types
  curr_metadata$cluster_1$cluster_cell_proportions <- cell_props
  curr_metadata$cluster_1$centre_loc <- centre_loc
  curr_metadata$cluster_1$ring_cell_types <- ring_cell_types
  curr_metadata$cluster_1$ring_cell_proportions <- ring_cell_props
  curr_metadata$cluster_1$ring_width <- ring_width
  
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
    curr_metadata$cluster_1$x_y_rotation <- ellipsoid_x_y_rotation
    curr_metadata$cluster_1$x_z_rotation <- ellipsoid_x_z_rotation
    curr_metadata$cluster_1$y_z_rotation <- ellipsoid_y_z_rotation
    
  }
  else if (shape == "Network") {
    curr_metadata$cluster_1$n_edges <- network_n_edges
    curr_metadata$cluster_1$width <- shape_parameters["width", size]
    curr_metadata$cluster_1$radius <- shape_parameters["radius", size]
  }
  else {
    stop(paste(shape, "Shape not found"))
  }
  
  ringed_spes_metadata[[i]] <- curr_metadata
}



## 2.3. Generate metadata for separated clusters ----------------------------------
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

setwd("~Objects/spes_table")
# write.table(separated_spes_table, "separated_spes_table.csv")


# Define shape parameters based on the size (Small, Medium, Large)
shape_parameters_list <- list()
shape_parameters_list[["Sphere"]] <- data.frame(Small = c(20), Medium = c(30), Large = c(40))
rownames(shape_parameters_list[["Sphere"]]) <- c("radius")
shape_parameters_list[["Ellipsoid"]] <- data.frame(Small = c(25, 20, 15), Medium = c(37.5, 30, 22.5), Large = c(50, 40, 30))
rownames(shape_parameters_list[["Ellipsoid"]]) <- c("x_radius", "y_radius", "z_radius")
shape_parameters_list[["Network"]] <- data.frame(Small = c(3, 30), Medium = c(6, 40), Large = c(9, 50))
rownames(shape_parameters_list[["Network"]]) <- c("width", "radius")


# Choose cluster centre location for each arrangement
centre_x_loc <- data.frame(clusterA = c(40, 50, 60),
                           clusterB = c(160, 150, 140))
rownames(centre_x_loc) <- c("S1", "S2", "S3")

# Define constant shape parameters (constant regardless of size)
ellipsoidA_x_y_rotation <- 45
ellipsoidA_x_z_rotation <- 0
ellipsoidA_y_z_rotation <- 0
ellipsoidB_x_y_rotation <- 0
ellipsoidB_x_z_rotation <- 45
ellipsoidB_y_z_rotation <- 0
network_n_edges <- 15
clusterA_cell_types <- c("A")
clusterA_cell_props <- c(1)
clusterB_cell_types <- c("B")
clusterB_cell_props <- c(1)

### Generate metadata for ringed simulations

# Get background data
setwd("~/Objects/spes_metadata")
bg_spes_metadata <- readRDS("bg_spes_metadata_200_100_100.rds")

# Get table for ringed simulations
setwd("~/Objects/spes_table")
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
  curr_metadata <- spe_metadata_cluster_template(bg_spe_metadata, "regular", shape)
  curr_metadata$cluster_1$cluster_cell_types <- clusterA_cell_types
  curr_metadata$cluster_1$cluster_cell_proportions <- clusterA_cell_props
  curr_metadata$cluster_1$centre_loc <- centre_locA
  
  curr_metadata <- spe_metadata_cluster_template(curr_metadata, "regular", shape)
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
    curr_metadata$cluster_1$x_y_rotation <- ellipsoidA_x_y_rotation
    curr_metadata$cluster_1$x_z_rotation <- ellipsoidA_x_z_rotation
    curr_metadata$cluster_1$y_z_rotation <- ellipsoidA_y_z_rotation
  }
  else if (shapeA == "Network") {
    curr_metadata$cluster_1$n_edges <- network_n_edges
    curr_metadata$cluster_1$width <- shape_parametersA["width", sizeA]
    curr_metadata$cluster_1$radius <- shape_parametersA["radius", sizeA]
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
    curr_metadata$cluster_2$x_y_rotation <- ellipsoidB_x_y_rotation
    curr_metadata$cluster_2$x_z_rotation <- ellipsoidB_x_z_rotation
    curr_metadata$cluster_2$y_z_rotation <- ellipsoidB_y_z_rotation
  }
  else if (shapeB == "Network") {
    curr_metadata$cluster_2$n_edges <- network_n_edges
    curr_metadata$cluster_2$width <- shape_parametersB["width", sizeB]
    curr_metadata$cluster_2$radius <- shape_parametersB["radius", sizeB]
  }
  else {
    stop(paste(shapeB, "ShapeB not found"))
  }
  separated_spes_metadata[[i]] <- curr_metadata
}


## 3.1. Simulate mixed clusters from metadata --------------------------------
setwd("~/Objects/spes_metadata")
mixed_spes_metadata <- readRDS("mixed_spes_metadata.rds")

setwd("~/Objects/mixed_spes")
i <- 1
for (mixed_spe_metadata in mixed_spes_metadata) {
  
  curr_spe <- simulate_spe_metadata3D(mixed_spe_metadata)
  file_name <- paste("mixed_spe_", i, ".rds", sep = "")
  
  saveRDS(curr_spe, file = file_name)
  
  i <- i + 1
}


## 3.2. Simulate ringed clusters from metadata --------------------------------
setwd("~/Objects/spes_metadata")
ringed_spes_metadata <- readRDS("ringed_spes_metadata.rds")

setwd("~/Objects/ringed_spes")
i <- 1
for (ringed_spe_metadata in ringed_spes_metadata) {
  
  curr_spe <- simulate_spe_metadata3D(ringed_spe_metadata)
  file_name <- paste("ringed_spe_", i, ".rds", sep = "")
  
  saveRDS(curr_spe, file = file_name)
  
  i <- i + 1
}
## 3.3. Simulate separated clusters from metadata --------------------------------
setwd("~/Objects/spes_metadata")
separated_spes_metadata <- readRDS("separated_spes_metadata.rds")

setwd("~/Objects/spes_table")
separated_spes_table <- read.table("separated_spes_table.csv")

setwd("~/Objects/separated_spes")
i <- 1
for (separated_spe_metadata in separated_spes_metadata) {

  if (i <= 384 || (
      separated_spes_table[i, "shapeA"] == separated_spes_table[i, "shapeB"] &&
      separated_spes_table[i, "sizeA"] == separated_spes_table[i, "sizeB"] &&
      separated_spes_table[i, "arrangement"] == "S2")) {
    i <- i + 1
    next
  }
  if (separated_spes_table[i, "shapeA"] != separated_spes_table[i, "shapeB"] ||
      separated_spes_table[i, "sizeA"] != separated_spes_table[i, "sizeB"] ||
      separated_spes_table[i, "arrangement"] == "S2") {
    i <- i + 1
    next
  }
  
  print(i)
  
  curr_spe <- simulate_spe_metadata3D(separated_spe_metadata, plot_image = FALSE)
  file_name <- paste("separated_spe_", i, ".rds", sep = "")
  
  saveRDS(curr_spe, file = file_name)
  
  i <- i + 1
}

## Spacer --------------------------------------------------------------------



