### 1.0. Set up all parameters for tables-----------------------------------------------

## Background
bg_prop_A_range <- c("min" = 0, "max" = 0.10)   # bg_prop_A ranges from 0 to 0.1
bg_prop_A_equals_zero_prob <- 1/3               # Give higher probability of bg_prop_A equals 0
bg_prop_B_range <- c("min" = 0, "max" = 0.10)   # bg_prop_B ranges from 0 to 0.1
bg_prop_B_equals_zero_prob <- 1/3               # Give higher probability of bg_prop_B equals 0

# Function to get bg props using bg_prop_range and probability prop equals 0
get_bg_props <- function(n_simulations, bg_prop_range, bg_prop_equals_zero_prob) {
  probs <- runif(n_simulations)
  props <- ifelse(probs <= bg_prop_equals_zero_prob, 0, round(runif(n_simulations, min = bg_prop_range["min"], max = bg_prop_range["max"]), 3))
  return(props)
}


## Shapes
shapes <- c("Sphere", "Ellipsoid", "Network")
radius_S_range <- c("min" = 25, "max" = 50)

radius_x_E_range <- c("min" = 25, "max" = 50)
radius_y_E_range <- c("min" = 25, "max" = 50)
radius_z_E_range <- c("min" = 25, "max" = 50)

radius_N_range <- c("min" = 25, "max" = 50)
width_N_range <- c("min" = 5, "max" = 10)


## Mixed clusters
cluster_prop_A_range <- c("min" = 0, "max" = 1)   # Proportion of A cells in cluster will range from 0 to 1, proportion of B cells in cluster will be 1 - prop(A)


## Ringed clusters
width_ring_range_factor <- c("min" = 0.1, "max" = 0.2)  # The width of the ring is proportional to the size (radius/width) of the cluster, for an ellipsoid, use the average of the three radii


## Separated clusters
centre_x_coord_A_range <- c("min" = 25, "max" = 75)   # centre_x_coord_B_range = 200 - centre_x_coord_A_range


### 1.1. Generate mixed_spes_table----------------------------------------------
n_mixed_simulations <- 5000
mixed_spes_table_colnames <- c("bg_prop_A", "bg_prop_B", "shape", "radius_S", "radius_x_E", "radius_y_E", "radius_z_E", "radius_N", "width_N", "cluster_prop_A", "cluster_prop_B")
mixed_spes_table <- data.frame(matrix(nrow = n_mixed_simulations, ncol = length(mixed_spes_table_colnames))) 
colnames(mixed_spes_table) <- mixed_spes_table_colnames

mixed_spes_table$bg_prop_A <- get_bg_props(n_mixed_simulations, bg_prop_A_range, bg_prop_A_equals_zero_prob)
mixed_spes_table$bg_prop_B <- get_bg_props(n_mixed_simulations, bg_prop_B_range, bg_prop_B_equals_zero_prob)

mixed_spes_table$shape <- sample(shapes, n_mixed_simulations, replace = T)
mixed_spes_table$radius_S <- ifelse(mixed_spes_table$shape == "Sphere", runif(n_mixed_simulations, radius_S_range["min"], radius_S_range["max"]), NA)
mixed_spes_table$radius_x_E <- ifelse(mixed_spes_table$shape == "Ellipsoid", runif(n_mixed_simulations, radius_x_E_range["min"], radius_x_E_range["max"]), NA)
mixed_spes_table$radius_y_E <- ifelse(mixed_spes_table$shape == "Ellipsoid", runif(n_mixed_simulations, radius_y_E_range["min"], radius_y_E_range["max"]), NA)
mixed_spes_table$radius_z_E <- ifelse(mixed_spes_table$shape == "Ellipsoid", runif(n_mixed_simulations, radius_z_E_range["min"], radius_z_E_range["max"]), NA)
mixed_spes_table$radius_N <- ifelse(mixed_spes_table$shape == "Network", runif(n_mixed_simulations, radius_N_range["min"], radius_N_range["max"]), NA)
mixed_spes_table$width_N <- ifelse(mixed_spes_table$shape == "Network", runif(n_mixed_simulations, width_N_range["min"], width_N_range["max"]), NA)

mixed_spes_table$cluster_prop_A <- round(runif(n_mixed_simulations, cluster_prop_A_range["min"], cluster_prop_A_range["max"]), 3)
mixed_spes_table$cluster_prop_B <- 1 - mixed_spes_table$cluster_prop_A


### 1.2. Generate ringed_spes_table  -------------------------------------
n_ringed_simulations <- 5000
ringed_spes_table_colnames <- c("bg_prop_A", "bg_prop_B", "shape", "radius_S", "radius_x_E", "radius_y_E", "radius_z_E", "radius_N", "width_N", "width_ring")
ringed_spes_table <- data.frame(matrix(nrow = n_ringed_simulations, ncol = length(ringed_spes_table_colnames))) 
colnames(ringed_spes_table) <- ringed_spes_table_colnames

ringed_spes_table$bg_prop_A <- get_bg_props(n_ringed_simulations, bg_prop_A_range, bg_prop_A_equals_zero_prob)
ringed_spes_table$bg_prop_B <- get_bg_props(n_ringed_simulations, bg_prop_B_range, bg_prop_B_equals_zero_prob)

ringed_spes_table$shape <- sample(shapes, n_ringed_simulations, replace = T)
ringed_spes_table$radius_S <- ifelse(ringed_spes_table$shape == "Sphere", runif(n_ringed_simulations, radius_S_range["min"], radius_S_range["max"]), NA)
ringed_spes_table$radius_x_E <- ifelse(ringed_spes_table$shape == "Ellipsoid", runif(n_ringed_simulations, radius_x_E_range["min"], radius_x_E_range["max"]), NA)
ringed_spes_table$radius_y_E <- ifelse(ringed_spes_table$shape == "Ellipsoid", runif(n_ringed_simulations, radius_y_E_range["min"], radius_y_E_range["max"]), NA)
ringed_spes_table$radius_z_E <- ifelse(ringed_spes_table$shape == "Ellipsoid", runif(n_ringed_simulations, radius_z_E_range["min"], radius_z_E_range["max"]), NA)
ringed_spes_table$radius_N <- ifelse(ringed_spes_table$shape == "Network", runif(n_ringed_simulations, radius_N_range["min"], radius_N_range["max"]), NA)
ringed_spes_table$width_N <- ifelse(ringed_spes_table$shape == "Network", runif(n_ringed_simulations, width_N_range["min"], width_N_range["max"]), NA)

radius_width_colnames <- c( "radius_S", "radius_x_E", "radius_y_E", "radius_z_E", "width_N") # Ignore radius_N as this refers to the radius spanned by the cluster, not the width of a branch 
ringed_spes_table$width_ring <- runif(n_ringed_simulations, 
                                      width_ring_range_factor["min"], width_ring_range_factor["max"]) * apply(ringed_spes_table[, radius_width_colnames], 1, mean, na.rm = T)


### 1.3. Generate separated_spes_table --------------------------------------
n_separated_simulations <- 5000
separated_spes_table_colnames <- c("bg_prop_A", "bg_prop_B", 
                                   "shape_A", "radius_S_A", "radius_x_E_A", "radius_y_E_A", "radius_z_E_A", "radius_N_A", "width_N_A", "centre_x_coord_A",
                                   "shape_B", "radius_S_B", "radius_x_E_B", "radius_y_E_B", "radius_z_E_B", "radius_N_B", "width_N_B", "centre_x_coord_B")
separated_spes_table <- data.frame(matrix(nrow = n_separated_simulations, ncol = length(separated_spes_table_colnames))) 
colnames(separated_spes_table) <- separated_spes_table_colnames

separated_spes_table$bg_prop_A <- get_bg_props(n_separated_simulations, bg_prop_A_range, bg_prop_A_equals_zero_prob)
separated_spes_table$bg_prop_B <- get_bg_props(n_separated_simulations, bg_prop_B_range, bg_prop_B_equals_zero_prob)

separated_spes_table$shape_A <- sample(shapes, n_separated_simulations, replace = T)
separated_spes_table$radius_S_A <- ifelse(separated_spes_table$shape_A == "Sphere", runif(n_separated_simulations, radius_S_range["min"], radius_S_range["max"]), NA)
separated_spes_table$radius_x_E_A <- ifelse(separated_spes_table$shape_A == "Ellipsoid", runif(n_separated_simulations, radius_x_E_range["min"], radius_x_E_range["max"]), NA)
separated_spes_table$radius_y_E_A <- ifelse(separated_spes_table$shape_A == "Ellipsoid", runif(n_separated_simulations, radius_y_E_range["min"], radius_y_E_range["max"]), NA)
separated_spes_table$radius_z_E_A <- ifelse(separated_spes_table$shape_A == "Ellipsoid", runif(n_separated_simulations, radius_z_E_range["min"], radius_z_E_range["max"]), NA)
separated_spes_table$radius_N_A <- ifelse(separated_spes_table$shape_A == "Network", runif(n_separated_simulations, radius_N_range["min"], radius_N_range["max"]), NA)
separated_spes_table$width_N_A <- ifelse(separated_spes_table$shape_A == "Network", runif(n_separated_simulations, width_N_range["min"], width_N_range["max"]), NA)

separated_spes_table$shape_B <- sample(shapes, n_separated_simulations, replace = T)
separated_spes_table$radius_S_B <- ifelse(separated_spes_table$shape_B == "Sphere", runif(n_separated_simulations, radius_S_range["min"], radius_S_range["max"]), NA)
separated_spes_table$radius_x_E_B <- ifelse(separated_spes_table$shape_B == "Ellipsoid", runif(n_separated_simulations, radius_x_E_range["min"], radius_x_E_range["max"]), NA)
separated_spes_table$radius_y_E_B <- ifelse(separated_spes_table$shape_B == "Ellipsoid", runif(n_separated_simulations, radius_y_E_range["min"], radius_y_E_range["max"]), NA)
separated_spes_table$radius_z_E_B <- ifelse(separated_spes_table$shape_B == "Ellipsoid", runif(n_separated_simulations, radius_z_E_range["min"], radius_z_E_range["max"]), NA)
separated_spes_table$radius_N_B <- ifelse(separated_spes_table$shape_B == "Network", runif(n_separated_simulations, radius_N_range["min"], radius_N_range["max"]), NA)
separated_spes_table$width_N_B <- ifelse(separated_spes_table$shape_B == "Network", runif(n_separated_simulations, width_N_range["min"], width_N_range["max"]), NA)

separated_spes_table$centre_x_coord_A <- runif(n_separated_simulations, centre_x_coord_A_range["min"], centre_x_coord_A_range["max"])
separated_spes_table$centre_x_coord_B <- 200 - separated_spes_table$centre_x_coord_A # Where 200 is the 'length' of the separated window

### 2.4. Save tables
setwd("~/Objects/unsupervised/spes_table")
# write.table(mixed_spes_table, "mixed_spes_table_unsupervised.csv")
# write.table(ringed_spes_table, "ringed_spes_table_unsupervised.csv")
# write.table(separated_spes_table, "separated_spes_table_unsupervised.csv")


### 2.0. Set up all parameters for metadata -----------------------------

# bg metadata for mixed or ringed simulations
bg_mixed_ringed_metadata <- spe_metadata_background_template("random")
bg_mixed_ringed_metadata$background$n_cells <- 20000
bg_mixed_ringed_metadata$background$length <- 100
bg_mixed_ringed_metadata$background$width <- 100
bg_mixed_ringed_metadata$background$height <- 100
bg_mixed_ringed_metadata$background$minimum_distance_between_cells <- 0
bg_mixed_ringed_metadata$background$cell_types <- c("A", "B", "O") # Cell proportions will change later

# bg metadata for separated simulations
bg_separated_metadata <- spe_metadata_background_template("random")
bg_separated_metadata$background$n_cells <- 40000
bg_separated_metadata$background$length <- 200
bg_separated_metadata$background$width <- 100
bg_separated_metadata$background$height <- 100
bg_separated_metadata$background$minimum_distance_between_cells <- 0
bg_separated_metadata$background$cell_types <- c("A", "B", "O") # Cell proportions will change later

# Network
n_edges_N <- 15

# Mixed clusters
mixed_cluster_cell_types <- c("A", "B")
mixed_cluster_centre_loc <- c(50, 50, 50)


### 2.1. Generate mixed_spes_metadata -----------------------------------
# Get table for mixed simulations
setwd("~/Objects/unsupervised/spes_table")
mixed_spes_table <- read.table("mixed_spes_table_unsupervised.csv")

# Set up metadata list
mixed_spes_metadata <- list()

for (i in seq(nrow(mixed_spes_table))) {
  
  # Get metadata template for current simulation parameters
  shape <- mixed_spes_table$shape[i]
  curr_metadata <- spe_metadata_cluster_template(bg_mixed_ringed_metadata, "regular", shape)
  
  curr_metadata$background$cell_proportions <- c(mixed_spes_table$bg_prop_A[i], 
                                                 mixed_spes_table$bg_prop_B[i],
                                                 1 - mixed_spes_table$bg_prop_A[i] - mixed_spes_table$bg_prop_B[i]) # prop(O) = 1 - prop(A) - prop(B)
  
  curr_metadata$cluster_1$cluster_cell_types <- mixed_cluster_cell_types
  curr_metadata$cluster_1$cluster_cell_proportions <- c(mixed_spes_table$cluster_prop_A[i], mixed_spes_table$cluster_prop_B[i])
  curr_metadata$cluster_1$centre_loc <- mixed_cluster_centre_loc
  
  # Specify metadata for each shape and size
  if (shape == "Sphere") {
    curr_metadata$cluster_1$radius <- mixed_spes_table$radius_S[i]
  }
  else if (shape == "Ellipsoid") {
    curr_metadata$cluster_1$x_radius <- mixed_spes_table$radius_x_E[i]
    curr_metadata$cluster_1$y_radius <- mixed_spes_table$radius_y_E[i]
    curr_metadata$cluster_1$z_radius <- mixed_spes_table$radius_z_E[i]
    curr_metadata$cluster_1$x_y_rotation <- runif(1, 0, 90)
    curr_metadata$cluster_1$x_z_rotation <- runif(1, 0, 90)
    curr_metadata$cluster_1$y_z_rotation <- runif(1, 0, 90)
    
  }
  else if (shape == "Network") {
    curr_metadata$cluster_1$n_edges <- n_edges_N
    curr_metadata$cluster_1$width <- mixed_spes_table$width_N[i]
    curr_metadata$cluster_1$radius <- mixed_spes_table$radius_N[i]
  }
  mixed_spes_metadata[[i]] <- curr_metadata
}


i <- 1492
md <- mixed_spes_metadata[[i]]
print(md)
simulate_spe_metadata3D(md)  



### 2.2. Generate ringed_spes_metadata ----------------------------------
### 2.3. Generate separated_spes_metadata -------------------------------
### 2.4. Save metadatas -------------------------------------------------