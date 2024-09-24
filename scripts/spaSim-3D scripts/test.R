# Get background metadata
bg_metadata <- spe_metadata_background_template("random")

# Change background metadata
bg_metadata$background$n_cells <- 20000
bg_metadata$background$length <- 600
bg_metadata$background$width <- 600
bg_metadata$background$height <- 600
bg_metadata$background$minimum_distance_between_cells <- 10
bg_metadata$background$cell_types <- "O"
bg_metadata$background$cell_proportions <- 1




# Test mixed cluster ellipsoid looks
mixed_metadata <- spe_metadata_cluster_template("regular", "Ellipsoid", bg_metadata)
mixed_metadata$cluster_1$cluster_cell_types <- c("A", "B")
mixed_metadata$cluster_1$cluster_cell_proportions <- c(0.6, 0.4)
mixed_metadata$cluster_1$x_radius <- 75
mixed_metadata$cluster_1$y_radius <- 75
mixed_metadata$cluster_1$z_radius <- 75
mixed_metadata$cluster_1$centre_loc <- c(300, 300, 300)

mixed_clusters <- simulate_spe_metadata3D(mixed_metadata)
plot_cells3D(mixed_clusters, c("O", "A", "B"), c("lightgray", "orange", "skyblue"))


# Test mixed cluster network looks
mixed_metadata <- spe_metadata_cluster_template("regular", "Network", bg_metadata)
mixed_metadata$cluster_1$cluster_cell_types <- c("A", "B")
mixed_metadata$cluster_1$cluster_cell_proportions <- c(0.6, 0.4)
mixed_metadata$cluster_1$radius <- 125
mixed_metadata$cluster_1$width <- 35
mixed_metadata$cluster_1$n_edges <- 20
mixed_metadata$cluster_1$centre_loc <- c(300, 300, 300)

mixed_clusters <- simulate_spe_metadata3D(mixed_metadata)
plot_cells3D(mixed_clusters, c("O", "A", "B"), c("lightgray", "orange", "skyblue"))




# Test ringed cluster ellipsoid looks
ring_width_factor <- 0.2
e_radius <- 125

ringed_metadata <- spe_metadata_cluster_template("ring", "Ellipsoid", bg_metadata)
ringed_metadata$cluster_1$cluster_cell_types <- c("A")
ringed_metadata$cluster_1$cluster_cell_proportions <- c(1)
ringed_metadata$cluster_1$x_radius <- e_radius
ringed_metadata$cluster_1$y_radius <- e_radius
ringed_metadata$cluster_1$z_radius <- e_radius
ringed_metadata$cluster_1$centre_loc <- c(300, 300, 300)
ringed_metadata$cluster_1$ring_cell_types <- c("B")
ringed_metadata$cluster_1$ring_cell_proportions <- c(1)
ringed_metadata$cluster_1$ring_width <- ring_width_factor * e_radius

ringed_clusters <- simulate_spe_metadata3D(ringed_metadata)
plot_cells3D(ringed_clusters, c("O", "A", "B"), c("lightgray", "orange", "skyblue"))


# Test ringed cluster network looks
ring_width_factor <- 0.15
n_width <- 45

ringed_metadata <- spe_metadata_cluster_template("ring", "Network", bg_metadata)
ringed_metadata$cluster_1$cluster_cell_types <- c("A")
ringed_metadata$cluster_1$cluster_cell_proportions <- c(1)
ringed_metadata$cluster_1$radius <- 125
ringed_metadata$cluster_1$width <- 35
ringed_metadata$cluster_1$n_edges <- 20
ringed_metadata$cluster_1$centre_loc <- c(300, 300, 300)
ringed_metadata$cluster_1$ring_cell_types <- c("B")
ringed_metadata$cluster_1$ring_cell_proportions <- c(1)
ringed_metadata$cluster_1$ring_width <- ring_width_factor * n_width

ringed_clusters <- simulate_spe_metadata3D(ringed_metadata)
plot_cells3D(ringed_clusters, c("O", "A", "B"), c("lightgray", "orange", "skyblue"))



# Test separated cluster E and N
x_coord <- 125
e_radius <- 100
n_width <- 35

separated_metadata <- spe_metadata_cluster_template("regular", "Ellipsoid", bg_metadata)
separated_metadata$cluster_1$cluster_cell_types <- c("A", "B")
separated_metadata$cluster_1$cluster_cell_proportions <- c(0.6, 0.4)
separated_metadata$cluster_1$x_radius <- e_radius
separated_metadata$cluster_1$y_radius <- e_radius
separated_metadata$cluster_1$z_radius <- e_radius
separated_metadata$cluster_1$centre_loc <- c(x_coord, 300, 300)

separated_metadata <- spe_metadata_cluster_template("regular", "Network", separated_metadata)
separated_metadata$cluster_2$cluster_cell_types <- c("A", "B")
separated_metadata$cluster_2$cluster_cell_proportions <- c(0.6, 0.4)
separated_metadata$cluster_2$radius <- 125
separated_metadata$cluster_2$width <- n_width
separated_metadata$cluster_2$n_edges <- 20
separated_metadata$cluster_2$centre_loc <- c(600 - x_coord, 300, 300)

separated_clusters <- simulate_spe_metadata3D(separated_metadata)
plot_cells3D(separated_clusters, c("O", "A", "B"), c("lightgray", "orange", "skyblue"))




# Test separated cluster E and E
x_coord <- 175
e_radius <- 125

separated_metadata <- spe_metadata_cluster_template("regular", "Ellipsoid", bg_metadata)
separated_metadata$cluster_1$cluster_cell_types <- c("A", "B")
separated_metadata$cluster_1$cluster_cell_proportions <- c(0.6, 0.4)
separated_metadata$cluster_1$x_radius <- e_radius
separated_metadata$cluster_1$y_radius <- e_radius
separated_metadata$cluster_1$z_radius <- e_radius
separated_metadata$cluster_1$centre_loc <- c(x_coord, 300, 300)

separated_metadata <- spe_metadata_cluster_template("regular", "Ellipsoid", separated_metadata)
separated_metadata$cluster_2$cluster_cell_types <- c("A", "B")
separated_metadata$cluster_2$cluster_cell_proportions <- c(0.6, 0.4)
separated_metadata$cluster_2$x_radius <- e_radius
separated_metadata$cluster_2$y_radius <- e_radius
separated_metadata$cluster_2$z_radius <- e_radius
separated_metadata$cluster_2$centre_loc <- c(600 - x_coord, 300, 300)

separated_clusters <- simulate_spe_metadata3D(separated_metadata)
plot_cells3D(separated_clusters, c("O", "A", "B"), c("lightgray", "orange", "skyblue"))




# Test separated cluster N and N
x_coord <- 175
n_width <- 35

separated_metadata <- spe_metadata_cluster_template("regular", "Network", bg_metadata)
separated_metadata$cluster_1$cluster_cell_types <- c("A", "B")
separated_metadata$cluster_1$cluster_cell_proportions <- c(0.6, 0.4)
separated_metadata$cluster_1$radius <- 125
separated_metadata$cluster_1$width <- n_width
separated_metadata$cluster_1$n_edges <- 20
separated_metadata$cluster_1$centre_loc <- c(x_coord, 300, 300)

separated_metadata <- spe_metadata_cluster_template("regular", "Network", separated_metadata)
separated_metadata$cluster_2$cluster_cell_types <- c("A", "B")
separated_metadata$cluster_2$cluster_cell_proportions <- c(0.6, 0.4)
separated_metadata$cluster_2$radius <- 125
separated_metadata$cluster_2$width <- n_width
separated_metadata$cluster_2$n_edges <- 20
separated_metadata$cluster_2$centre_loc <- c(600 - x_coord, 300, 300)

separated_clusters <- simulate_spe_metadata3D(separated_metadata)
plot_cells3D(separated_clusters, c("O", "A", "B"), c("lightgray", "orange", "skyblue"))

