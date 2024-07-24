### Using the scripts found in the "MAIN USER FUNCTIONS" folder


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
bg_metadata$background$minimum_distance_between_cells <- 0
bg_metadata$background$n_cells <- 15000

# Cluster metadata (using background metadata as background)
cluster_metadata <- spe_metadata_cluster_template(bg_metadata, "regular", "Sphere")
cluster_metadata <- spe_metadata_cluster_template(cluster_metadata, "ring", "Ellipsoid")
cluster_metadata <- spe_metadata_cluster_template(cluster_metadata, "double ring", "Cylinder")

# Get spe from updated metadata
spe_clusters <- simulate_spe_metadata3D(cluster_metadata)
plot_cells3D(spe_clusters,
             plot_cell_types = c("Others", "Tumour", "Immune", "Immune1", "Endothelial"),
             plot_colours = c("lightgray", "orange", "skyblue", "lightgreen", "tomato"))



# Add metadata to spe
new_metadata <- spe_metadata_cluster_template(bg_metadata, "regular", "Network")

spe_clusters1 <- add_spe_metadata3D(spe_clusters, new_metadata)
plot_cells3D(spe_clusters1,
             plot_cell_types = c("Others", "Tumour", "Immune", "Immune1", "Endothelial"),
             plot_colours = c("lightgray", "orange", "skyblue", "lightgreen", "tomato"))




### Random simuation ---------------------------------------------------------
metadata_bg_r <- spe_metadata_background_template("random")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_r, "double ring", "Sphere")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "regular", "Network")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "regular", "Network")

metadata_bg_clusters$cluster_1$centre_loc <- c(50, 50, 50)
metadata_bg_clusters$cluster_2$centre_loc <- c(25, 25, 25)
metadata_bg_clusters$cluster_3$centre_loc <- c(75, 75, 75)


spe_clusters <- simulate_spe_metadata3D(metadata_bg_clusters)
plot_cells3D(spe_clusters,
             plot_cell_types = c("Others", "Tumour", "Immune", "Immune1", "Immune2"),
             plot_colours = c("lightgray", "orange", "skyblue", "lightgreen", "purple"))


### Testing every default clustering -----------------------------------------
spe_metadata <- spe_metadata_background_template("random")
spe_metadata <- spe_metadata_cluster_template(spe_metadata, "regular", "Ellipsoid")
spe_metadata <- spe_metadata_cluster_template(spe_metadata, "regular", "Cylinder")
spe_metadata <- spe_metadata_cluster_template(spe_metadata, "regular", "Network")
spe_metadata <- spe_metadata_cluster_template(spe_metadata, "ring", "Sphere")
spe_metadata <- spe_metadata_cluster_template(spe_metadata, "ring", "Ellipsoid")
spe_metadata <- spe_metadata_cluster_template(spe_metadata, "ring", "Cylinder")
spe_metadata <- spe_metadata_cluster_template(spe_metadata, "ring", "Network")
spe_metadata <- spe_metadata_cluster_template(spe_metadata, "double ring", "Sphere")
spe_metadata <- spe_metadata_cluster_template(spe_metadata, "double ring", "Ellipsoid")
spe_metadata <- spe_metadata_cluster_template(spe_metadata, "double ring", "Cylinder")
spe_metadata <- spe_metadata_cluster_template(spe_metadata, "double ring", "Network")
spe_clusters <- simulate_spe_metadata3D(spe_metadata)

