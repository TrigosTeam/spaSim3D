metadata_bg_r <- spe_metadata_background_template("random")
metadata_bg_r[["background"]][["n_cells"]] <- 20000
metadata_bg_r[["background"]][["cell_types"]] <- c("Tumour", "Immune", "Others")
metadata_bg_r[["background"]][["cell_proportions"]] <- c(0.03, 0.03, 0.94)
metadata_bg_r[["background"]][["minimum_distance_between_cells"]] <- 0


metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_r, "regular", "Sphere")
metadata_bg_clusters[["cluster_1"]][["cluster_cell_types"]] <- c("Tumour")
metadata_bg_clusters[["cluster_1"]][["cluster_cell_proportions"]] <- c(1)
metadata_bg_clusters[["cluster_1"]][["centre_loc"]] <- c(30, 30, 30)
metadata_bg_clusters[["cluster_1"]][["radius"]] <- 20

metadata_bg_clusters <- spe_metadata_cluster_template(metadata_bg_clusters, "ring", "Sphere")
metadata_bg_clusters[["cluster_2"]][["cluster_cell_types"]] <- c("Tumour")
metadata_bg_clusters[["cluster_2"]][["cluster_cell_proportions"]] <- c(1)
metadata_bg_clusters[["cluster_2"]][["centre_loc"]] <- c(70, 70, 70)
metadata_bg_clusters[["cluster_2"]][["radius"]] <- 20

spe_clusters <- simulate_spe_metadata3D(metadata_bg_clusters)
plot_cells3D(spe_clusters, c("Tumour", "Immune", "Others"), c("orange", "skyblue", "lightgray"))


spe_alpha_hull <- determine_alpha_hull3D(spe_clusters, c("Tumour"), alpha = 3, minimum_cells_in_alpha_hull = 10)

plot_alpha_hull3D(spe_alpha_hull, c("Tumour", "Immune", "Others"), c("orange", "skyblue", "lightgray"))

alpha_hull_props <- calculate_alpha_hull_cell_proportions3D(spe_alpha_hull)

alpha_hull_min_distances <- calculate_minimum_distances_to_alpha_hull3D(spe_alpha_hull, cell_types_of_interest = c("Tumour", "Immune"))
