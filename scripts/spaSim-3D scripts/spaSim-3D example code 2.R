### Using the integrator functions --------------------------------------------
spe_bg <- spaSim3D_background_integrator()

spe_cluster <- spaSim3D_cluster_integrator(spe_bg)
spe_cluster1 <- spaSim3D_cluster_integrator(spe_cluster)

# Plot with chosen colours
plot_cells3D(spe_bg, 
             plot_cell_types = c("Others", "Tumour", "Immune"),
             plot_colours = c("lightgray", "orange", "skyblue"))

plot_cells3D(spe_cluster,
             plot_cell_types = c("Tumour", "Immune1", "Immune2", "Others"),
             plot_colours = c("orange", "lightgreen", "skyblue", "lightgray"))

plot_cells3D(spe_cluster1,
             plot_cell_types = c("Tumour", "Immune", "Others"),
             plot_colours = c("orange", "skyblue", "lightgray"))

### Generating duplicate simulation using metadata from previous simulation -----
prev_spe <- spe_cluster
prev_metadata <- prev_spe@metadata
spe_cluster_dup <- simulate_spe_metadata3D(prev_metadata)


### Create your own spe metadata  -------------------------------------------
# Background metadata
bg_metadata <- spe_metadata_background_template("random")

# Change background metadata
bg_metadata$background$length <- 600
bg_metadata$background$width <- 600
bg_metadata$background$height <- 250
bg_metadata$background$minimum_distance_between_cells <- 10
bg_metadata$background$n_cells <- 25000

# Cluster metadata (using background metadata as background)
cluster_metadata <- spe_metadata_cluster_template("regular", "sphere", bg_metadata)
cluster_metadata <- spe_metadata_cluster_template("ring", "ellipsoid", cluster_metadata)
cluster_metadata <- spe_metadata_cluster_template("double ring", "cylinder", cluster_metadata)

# Get spe from updated metadata
spe_clusters <- simulate_spe_metadata3D(cluster_metadata)
plot_cells3D(spe_clusters,
             plot_cell_types = c("Others", "Tumour", "Immune", "Immune1", "Endothelial"),
             plot_colours = c("lightgray", "orange", "skyblue", "lightgreen", "tomato"))



# Add metadata to spe
new_metadata <- spe_metadata_cluster_template("regular", "network")

spe_clusters1 <- add_spe_metadata3D(spe_clusters, new_metadata)
plot_cells3D(spe_clusters1,
             plot_cell_types = c("Others", "Tumour", "Immune", "Immune1", "Endothelial"),
             plot_colours = c("lightgray", "orange", "skyblue", "lightgreen", "tomato"))




### Random simuation ---------------------------------------------------------
metadata_bg_r <- spe_metadata_background_template("random")
metadata_bg_clusters <- spe_metadata_cluster_template("double ring", "sphere", metadata_bg_r)
metadata_bg_clusters <- spe_metadata_cluster_template("regular", "network", metadata_bg_clusters)
metadata_bg_clusters <- spe_metadata_cluster_template("regular", "network", metadata_bg_clusters)

metadata_bg_clusters$cluster_1$centre_loc <- c(50, 50, 50)
metadata_bg_clusters$cluster_2$centre_loc <- c(25, 25, 25)
metadata_bg_clusters$cluster_3$centre_loc <- c(75, 75, 75)


spe_clusters <- simulate_spe_metadata3D(metadata_bg_clusters)
plot_cells3D(spe_clusters,
             plot_cell_types = c("Others", "Tumour", "Immune", "Immune1", "Immune2"),
             plot_colours = c("lightgray", "orange", "skyblue", "lightgreen", "purple"))


### Testing every default clustering -----------------------------------------
spe_metadata <- spe_metadata_background_template("random")
spe_metadata <- spe_metadata_cluster_template("regular", "sphere", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("regular", "ellipsoid", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("regular", "cylinder", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("regular", "network", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("ring", "sphere", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("ring", "ellipsoid", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("ring", "cylinder", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("ring", "network", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("double ring", "sphere", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("double ring", "ellipsoid", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("double ring", "cylinder", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("double ring", "network", spe_metadata)
spe_clusters <- simulate_spe_metadata3D(spe_metadata)

