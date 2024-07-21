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

### Using the metadata-based functions ----------------------------------------

# Generating duplicate simulation using metadata from previous simulation
prev_spe <- spe_cluster
prev_metadata <- prev_spe@metadata
spe_cluster_dup <- simulate_spe_metadata3D(prev_metadata)


# Get default background metadata
metadata_bg_r <- spe_metadata_background_template("random")


# Change background metadata
metadata_bg_r$background$minimum_distance_between_cells <- 0
metadata_bg_r$background$n_cells <- 11000

# Get spe from background metadata
# spe_bg1 <- simulate_spe_metadata3D(metadata_bg_r)



# Add to background metadata to get cluster metadata
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_r, "regular", "Sphere")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "ring", "Ellipsoid")
# metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "double ring", "Cylinder")
# metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "regular", "Network")

# Change cluster metadata
# metadata_bg_clusters$cluster_3$outer_ring_cell_types <- c("Immune", "Others")
# metadata_bg_clusters$cluster_4$width <- 10
# metadata_bg_clusters$cluster_4$cluster_cell_types <- "Immune"
# metadata_bg_clusters$cluster_4$cluster_cell_proportions <- 1


# Get spe from updated metadata
spe_clusters <- simulate_spe_metadata3D(metadata_bg_clusters)
plot_cells3D(spe_clusters,
             plot_cell_types = c("Others", "Tumour", "Immune", "Immune1", "Endothelial"),
             plot_colours = c("lightgray", "orange", "skyblue", "lightgreen", "tomato"))



# Add metadata to spe
metadata_new <- spe_metadata_cluster_template(metadata_bg_r, "regular", "Network")
metadata_new$cluster_1$width <- 10
metadata_new$cluster_1$cluster_cell_types <- "Immune"
metadata_new$cluster_1$cluster_cell_proportions <- 1

spe_clusters1 <- add_spe_metadata3D(spe_clusters, metadata_new)
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
metadata_bg_r <- spe_metadata_background_template("random")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_r, "regular", "Sphere")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "regular", "Ellipsoid")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "regular", "Cylinder")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "regular", "Network")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "ring", "Sphere")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "ring", "Ellipsoid")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "ring", "Cylinder")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "ring", "Network")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "double ring", "Sphere")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "double ring", "Ellipsoid")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "double ring", "Cylinder")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "double ring", "Network")
spe_clusters <- simulate_spe_metadata3D(metadata_bg_clusters)


### Error checking
metadata_bg_r <- spe_metadata_background_template("random")
metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_r, "double ring", "Sphere")
spe_clusters <- simulate_spe_metadata3D(metadata_bg_clusters)
