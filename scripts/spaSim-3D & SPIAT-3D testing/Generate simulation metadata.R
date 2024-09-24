### 1.0. Set up all parameters for tables-----------------------------------------------

## Background
bg_prop_A_range <- c("min" = 0, "max" = 0.10)   # bg_prop_A ranges from 0 to 0.1
bg_prop_A_equals_zero_prob <- 1/2               # Give higher probability of bg_prop_A equals 0
bg_prop_B_range <- c("min" = 0, "max" = 0.10)   # bg_prop_B ranges from 0 to 0.1
bg_prop_B_equals_zero_prob <- 1/2               # Give higher probability of bg_prop_B equals 0

# Function to get bg props using bg_prop_range and probability prop equals 0
get_bg_props <- function(n_simulations, bg_prop_range, bg_prop_equals_zero_prob) {
  probs <- runif(n_simulations)
  props <- ifelse(probs <= bg_prop_equals_zero_prob, 0, round(runif(n_simulations, min = bg_prop_range["min"], max = bg_prop_range["max"]), 3))
  return(props)
}


## Shapes
shapes <- c("Ellipsoid", "Network")

radius_x_E_range <- c("min" = 75, "max" = 125)
radius_y_E_range <- c("min" = 75, "max" = 125)
radius_z_E_range <- c("min" = 75, "max" = 125)

width_N_range <- c("min" = 30, "max" = 45)


## Mixed clusters
cluster_prop_A_range <- c("min" = 0.5, "max" = 0.9)   # Proportion of A cells in cluster will range from 0.5 to 0.9, proportion of B cells in cluster will be 1 - prop(A)


## Ringed clusters
ring_width_range_factor <- c("min" = 0.1, "max" = 0.2)  # The width of the ring is proportional to the size (radius/width) of the cluster, for an ellipsoid, use the average of the three radii


## Separated clusters
centre_x_coord_A_range <- c("min" = 125, "max" = 175)   # centre_x_coord_B_range = 600 - centre_x_coord_A_range


### 1.1. Generate mixed_spes_table----------------------------------------------
n_mixed_simulations <- 2000
mixed_spes_table_colnames <- c("bg_prop_A", "bg_prop_B", "shape", "radius_x_E", "radius_y_E", "radius_z_E", "width_N", "cluster_prop_A", "cluster_prop_B")
mixed_spes_table <- data.frame(matrix(nrow = n_mixed_simulations, ncol = length(mixed_spes_table_colnames))) 
colnames(mixed_spes_table) <- mixed_spes_table_colnames

mixed_spes_table$bg_prop_A <- get_bg_props(n_mixed_simulations, bg_prop_A_range, bg_prop_A_equals_zero_prob)
mixed_spes_table$bg_prop_B <- get_bg_props(n_mixed_simulations, bg_prop_B_range, bg_prop_B_equals_zero_prob)

mixed_spes_table$shape <- sample(shapes, n_mixed_simulations, replace = T)
mixed_spes_table$radius_x_E <- ifelse(mixed_spes_table$shape == "Ellipsoid", runif(n_mixed_simulations, radius_x_E_range["min"], radius_x_E_range["max"]), NA)
mixed_spes_table$radius_y_E <- ifelse(mixed_spes_table$shape == "Ellipsoid", runif(n_mixed_simulations, radius_y_E_range["min"], radius_y_E_range["max"]), NA)
mixed_spes_table$radius_z_E <- ifelse(mixed_spes_table$shape == "Ellipsoid", runif(n_mixed_simulations, radius_z_E_range["min"], radius_z_E_range["max"]), NA)
mixed_spes_table$width_N <- ifelse(mixed_spes_table$shape == "Network", runif(n_mixed_simulations, width_N_range["min"], width_N_range["max"]), NA)

mixed_spes_table$cluster_prop_A <- round(runif(n_mixed_simulations, cluster_prop_A_range["min"], cluster_prop_A_range["max"]), 3)
mixed_spes_table$cluster_prop_B <- 1 - mixed_spes_table$cluster_prop_A


### 1.2. Generate ringed_spes_table  -------------------------------------
n_ringed_simulations <- 2000
ringed_spes_table_colnames <- c("bg_prop_A", "bg_prop_B", "shape", "radius_x_E", "radius_y_E", "radius_z_E", "width_N", "ring_width")
ringed_spes_table <- data.frame(matrix(nrow = n_ringed_simulations, ncol = length(ringed_spes_table_colnames))) 
colnames(ringed_spes_table) <- ringed_spes_table_colnames

ringed_spes_table$bg_prop_A <- get_bg_props(n_ringed_simulations, bg_prop_A_range, bg_prop_A_equals_zero_prob)
ringed_spes_table$bg_prop_B <- get_bg_props(n_ringed_simulations, bg_prop_B_range, bg_prop_B_equals_zero_prob)

ringed_spes_table$shape <- sample(shapes, n_ringed_simulations, replace = T)
ringed_spes_table$radius_x_E <- ifelse(ringed_spes_table$shape == "Ellipsoid", runif(n_ringed_simulations, radius_x_E_range["min"], radius_x_E_range["max"]), NA)
ringed_spes_table$radius_y_E <- ifelse(ringed_spes_table$shape == "Ellipsoid", runif(n_ringed_simulations, radius_y_E_range["min"], radius_y_E_range["max"]), NA)
ringed_spes_table$radius_z_E <- ifelse(ringed_spes_table$shape == "Ellipsoid", runif(n_ringed_simulations, radius_z_E_range["min"], radius_z_E_range["max"]), NA)
ringed_spes_table$width_N <- ifelse(ringed_spes_table$shape == "Network", runif(n_ringed_simulations, width_N_range["min"], width_N_range["max"]), NA)

radius_width_colnames <- c("radius_x_E", "radius_y_E", "radius_z_E", "width_N") # Ignore radius_N as this refers to the radius spanned by the cluster, not the width of a branch
ringed_spes_table$ring_width_factor <- runif(n_ringed_simulations, ring_width_range_factor["min"], ring_width_range_factor["max"])
ringed_spes_table$ring_width <- ringed_spes_table$ring_width_factor * apply(ringed_spes_table[, radius_width_colnames], 1, mean, na.rm = T)


### 1.3. Generate separated_spes_table --------------------------------------
n_separated_simulations <- 2000
separated_spes_table_colnames <- c("bg_prop_A", "bg_prop_B", 
                                   "shape_A", "radius_x_E_A", "radius_y_E_A", "radius_z_E_A", "width_N_A", "centre_x_coord_A",
                                   "shape_B", "radius_x_E_B", "radius_y_E_B", "radius_z_E_B", "width_N_B", "centre_x_coord_B")
separated_spes_table <- data.frame(matrix(nrow = n_separated_simulations, ncol = length(separated_spes_table_colnames))) 
colnames(separated_spes_table) <- separated_spes_table_colnames

separated_spes_table$bg_prop_A <- get_bg_props(n_separated_simulations, bg_prop_A_range, bg_prop_A_equals_zero_prob)
separated_spes_table$bg_prop_B <- get_bg_props(n_separated_simulations, bg_prop_B_range, bg_prop_B_equals_zero_prob)

separated_spes_table$shape_A <- sample(shapes, n_separated_simulations, replace = T)
separated_spes_table$radius_x_E_A <- ifelse(separated_spes_table$shape_A == "Ellipsoid", runif(n_separated_simulations, radius_x_E_range["min"], radius_x_E_range["max"]), NA)
separated_spes_table$radius_y_E_A <- ifelse(separated_spes_table$shape_A == "Ellipsoid", runif(n_separated_simulations, radius_y_E_range["min"], radius_y_E_range["max"]), NA)
separated_spes_table$radius_z_E_A <- ifelse(separated_spes_table$shape_A == "Ellipsoid", runif(n_separated_simulations, radius_z_E_range["min"], radius_z_E_range["max"]), NA)
separated_spes_table$width_N_A <- ifelse(separated_spes_table$shape_A == "Network", runif(n_separated_simulations, width_N_range["min"], width_N_range["max"]), NA)

separated_spes_table$shape_B <- sample(shapes, n_separated_simulations, replace = T)
separated_spes_table$radius_x_E_B <- ifelse(separated_spes_table$shape_B == "Ellipsoid", runif(n_separated_simulations, radius_x_E_range["min"], radius_x_E_range["max"]), NA)
separated_spes_table$radius_y_E_B <- ifelse(separated_spes_table$shape_B == "Ellipsoid", runif(n_separated_simulations, radius_y_E_range["min"], radius_y_E_range["max"]), NA)
separated_spes_table$radius_z_E_B <- ifelse(separated_spes_table$shape_B == "Ellipsoid", runif(n_separated_simulations, radius_z_E_range["min"], radius_z_E_range["max"]), NA)
separated_spes_table$width_N_B <- ifelse(separated_spes_table$shape_B == "Network", runif(n_separated_simulations, width_N_range["min"], width_N_range["max"]), NA)

separated_spes_table$centre_x_coord_A <- runif(n_separated_simulations, centre_x_coord_A_range["min"], centre_x_coord_A_range["max"])
separated_spes_table$centre_x_coord_B <- 600 - separated_spes_table$centre_x_coord_A # Where 600 is the 'length' of the separated window

### 1.4. Save tables -------------------------------------------------
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/spe_tables")
write.table(mixed_spes_table, "mixed_spes_table.csv")
write.table(ringed_spes_table, "ringed_spes_table.csv")
write.table(separated_spes_table, "separated_spes_table.csv")


### 2.0. Set up all parameters for metadata -----------------------------

# bg metadata for simulations
bg_metadata <- spe_metadata_background_template("random")
bg_metadata$background$n_cells <- 20000
bg_metadata$background$length <- 600
bg_metadata$background$width <- 600
bg_metadata$background$height <- 600
bg_metadata$background$minimum_distance_between_cells <- 10
bg_metadata$background$cell_types <- c("A", "B", "O") # Cell proportions will change later


# Network
radius_N <- 125
n_edges_N <- 20

# Mixed clusters
mixed_cluster_cell_types <- c("A", "B")
mixed_cluster_centre_loc <- c(300, 300, 300)

# Ringed clusters
ringed_cluster_cell_type <- "A"
ringed_cluster_cell_prop <- 1
ringed_ring_cell_type <- "B"
ringed_ring_cell_prop <- 1
ringed_cluster_centre_loc <- c(300, 300, 300)

# Separated clusters
separated_cluster_A_cell_type <- "A"
separated_cluster_A_cell_prop <- 1
separated_cluster_B_cell_type <- "B"
separated_cluster_B_cell_prop <- 1
# Get centre coords for separated clusters later


### 2.1. Generate mixed_spes_metadata -----------------------------------
# Get table for mixed simulations
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/spe_tables")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Set up metadata list
mixed_spes_metadata <- list()

for (i in seq(nrow(mixed_spes_table))) {
  
  # Get metadata template for current simulation parameters
  shape <- mixed_spes_table$shape[i]
  curr_metadata <- spe_metadata_cluster_template("regular", shape, bg_metadata)
  
  curr_metadata$background$cell_proportions <- c(mixed_spes_table$bg_prop_A[i], 
                                                 mixed_spes_table$bg_prop_B[i],
                                                 1 - mixed_spes_table$bg_prop_A[i] - mixed_spes_table$bg_prop_B[i]) # prop(O) = 1 - prop(A) - prop(B)
  
  curr_metadata$cluster_1$cluster_cell_types <- mixed_cluster_cell_types
  curr_metadata$cluster_1$cluster_cell_proportions <- c(mixed_spes_table$cluster_prop_A[i], mixed_spes_table$cluster_prop_B[i])
  curr_metadata$cluster_1$centre_loc <- mixed_cluster_centre_loc
  
  # Specify metadata for each shape and size
  if (shape == "Ellipsoid") {
    curr_metadata$cluster_1$x_radius <- mixed_spes_table$radius_x_E[i]
    curr_metadata$cluster_1$y_radius <- mixed_spes_table$radius_y_E[i]
    curr_metadata$cluster_1$z_radius <- mixed_spes_table$radius_z_E[i]
    curr_metadata$cluster_1$x_y_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_1$x_z_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_1$y_z_rotation <- runif(1, 0, 180)
    
  }
  else if (shape == "Network") {
    curr_metadata$cluster_1$n_edges <- n_edges_N
    curr_metadata$cluster_1$width <- mixed_spes_table$width_N[i]
    curr_metadata$cluster_1$radius <- radius_N
  }
  mixed_spes_metadata[[i]] <- curr_metadata
}






### 2.2. Generate ringed_spes_metadata ----------------------------------
# Get table for ringed simulations
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/spe_tables")
ringed_spes_table <- read.table("ringed_spes_table.csv")

# Set up metadata list
ringed_spes_metadata <- list()

for (i in seq(nrow(ringed_spes_table))) {
  
  # Get metadata template for current simulation parameters
  shape <- ringed_spes_table$shape[i]
  curr_metadata <- spe_metadata_cluster_template("ring", shape, bg_metadata)
  
  curr_metadata$background$cell_proportions <- c(ringed_spes_table$bg_prop_A[i], 
                                                 ringed_spes_table$bg_prop_B[i],
                                                 1 - ringed_spes_table$bg_prop_A[i] - ringed_spes_table$bg_prop_B[i]) # prop(O) = 1 - prop(A) - prop(B)
  
  curr_metadata$cluster_1$cluster_cell_types <- ringed_cluster_cell_type
  curr_metadata$cluster_1$cluster_cell_proportions <- ringed_cluster_cell_prop
  curr_metadata$cluster_1$centre_loc <- ringed_cluster_centre_loc
  curr_metadata$cluster_1$ring_cell_types <- ringed_ring_cell_type
  curr_metadata$cluster_1$ring_cell_proportions <- ringed_ring_cell_prop
  curr_metadata$cluster_1$ring_width <- ringed_spes_table$ring_width[i]
  
  # Specify metadata for each shape and size
  if (shape == "Ellipsoid") {
    curr_metadata$cluster_1$x_radius <- ringed_spes_table$radius_x_E[i]
    curr_metadata$cluster_1$y_radius <- ringed_spes_table$radius_y_E[i]
    curr_metadata$cluster_1$z_radius <- ringed_spes_table$radius_z_E[i]
    curr_metadata$cluster_1$x_y_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_1$x_z_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_1$y_z_rotation <- runif(1, 0, 180)
    
  }
  else if (shape == "Network") {
    curr_metadata$cluster_1$n_edges <- n_edges_N
    curr_metadata$cluster_1$width <- ringed_spes_table$width_N[i]
    curr_metadata$cluster_1$radius <- radius_N
  }
  ringed_spes_metadata[[i]] <- curr_metadata
}


### 2.3. Generate separated_spes_metadata -------------------------------
# Get table for separated simulations
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/spe_tables")
separated_spes_table <- read.table("separated_spes_table.csv")

# Set up metadata list
separated_spes_metadata <- list()

for (i in seq(nrow(separated_spes_table))) {
  
  # Get metadata template for current simulation parameters
  shape_A <- separated_spes_table$shape_A[i]
  curr_metadata <- spe_metadata_cluster_template("regular", shape_A, bg_metadata)
  curr_metadata$cluster_1$cluster_cell_types <- separated_cluster_A_cell_type
  curr_metadata$cluster_1$cluster_cell_proportions <- separated_cluster_A_cell_prop
  curr_metadata$cluster_1$centre_loc <- c(separated_spes_table$centre_x_coord_A[i], 300, 300)
  
  shape_B <- separated_spes_table$shape_B[i]
  curr_metadata <- spe_metadata_cluster_template("regular", shape_B, curr_metadata)
  curr_metadata$cluster_2$cluster_cell_types <- separated_cluster_B_cell_type
  curr_metadata$cluster_2$cluster_cell_proportions <- separated_cluster_B_cell_prop
  curr_metadata$cluster_2$centre_loc <- c(separated_spes_table$centre_x_coord_B[i], 300, 300)
  
  curr_metadata$background$cell_proportions <- c(separated_spes_table$bg_prop_A[i], 
                                                 separated_spes_table$bg_prop_B[i],
                                                 1 - separated_spes_table$bg_prop_A[i] - separated_spes_table$bg_prop_B[i]) # prop(O) = 1 - prop(A) - prop(B)
  
  # Specify metadata for each shape and size for cluster_A
  if (shape_A == "Ellipsoid") {
    curr_metadata$cluster_1$x_radius <- separated_spes_table$radius_x_E_A[i]
    curr_metadata$cluster_1$y_radius <- separated_spes_table$radius_y_E_A[i]
    curr_metadata$cluster_1$z_radius <- separated_spes_table$radius_z_E_A[i]
    curr_metadata$cluster_1$x_y_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_1$x_z_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_1$y_z_rotation <- runif(1, 0, 180)
    
  }
  else if (shape_A == "Network") {
    curr_metadata$cluster_1$n_edges <- n_edges_N
    curr_metadata$cluster_1$width <- separated_spes_table$width_N_A[i]
    curr_metadata$cluster_1$radius <- radius_N
  }
  
  # Specify metadata for each shape and size for cluster_B
  if (shape_B == "Ellipsoid") {
    curr_metadata$cluster_2$x_radius <- separated_spes_table$radius_x_E_B[i]
    curr_metadata$cluster_2$y_radius <- separated_spes_table$radius_y_E_B[i]
    curr_metadata$cluster_2$z_radius <- separated_spes_table$radius_z_E_B[i]
    curr_metadata$cluster_2$x_y_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_2$x_z_rotation <- runif(1, 0, 180)
    curr_metadata$cluster_2$y_z_rotation <- runif(1, 0, 180)
    
  }
  else if (shape_B == "Network") {
    curr_metadata$cluster_2$n_edges <- n_edges_N
    curr_metadata$cluster_2$width <- separated_spes_table$width_N_B[i]
    curr_metadata$cluster_2$radius <- radius_N
  }
  
  separated_spes_metadata[[i]] <- curr_metadata
}

### 2.4. Save metadatas -------------------------------------------------
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/spe_metadata")
saveRDS(mixed_spes_metadata, "mixed_spes_metadata.rds")
saveRDS(ringed_spes_metadata, "ringed_spes_metadata.rds")
saveRDS(separated_spes_metadata, "separated_spes_metadata.rds")
