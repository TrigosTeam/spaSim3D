### Generating simulations ---------------------------------------------------

metadata_bg_r <- spe_metadata_background_template("random")
metadata_bg_r[["background"]][["cell_types"]] <- "Others"
metadata_bg_r[["background"]][["cell_proportions"]] <- 1


# Mixed sphere
metadata_mixed_sphere <- spe_metadata_cluster_template(metadata_bg_r, "regular", "Sphere")
metadata_mixed_sphere[["cluster_1"]][["cluster_cell_types"]] <- c("Tumour", "Immune")
metadata_mixed_sphere[["cluster_1"]][["cluster_cell_proportions"]] <- c(0.5, 0.5)
metadata_mixed_sphere[["cluster_1"]][["radius"]] <- 40
metadata_mixed_sphere[["cluster_1"]][["centre_loc"]] <- c(50, 50, 50)

spe_mixed_sphere <- simulate_spe_metadata3D(metadata_mixed_sphere)
plot_cells3D(spe_mixed_sphere,
             plot_cell_types = c("Others", "Tumour", "Immune"),
             plot_colours = c("lightgray", "orange", "skyblue"))



# Sphere within sphere
metadata_sphere_ception <- spe_metadata_cluster_template(metadata_bg_r, "regular", "Sphere")
metadata_sphere_ception <- spe_metadata_cluster_template(metadata_sphere_ception, "regular", "Sphere")

metadata_sphere_ception[["cluster_1"]][["cluster_cell_types"]] <- c("Tumour")
metadata_sphere_ception[["cluster_1"]][["cluster_cell_proportions"]] <- c(1)
metadata_sphere_ception[["cluster_1"]][["radius"]] <- 40
metadata_sphere_ception[["cluster_1"]][["centre_loc"]] <- c(50, 50, 50)

metadata_sphere_ception[["cluster_2"]][["cluster_cell_types"]] <- c("Immune")
metadata_sphere_ception[["cluster_2"]][["cluster_cell_proportions"]] <- c(1)
metadata_sphere_ception[["cluster_2"]][["radius"]] <- 25
metadata_sphere_ception[["cluster_2"]][["centre_loc"]] <- c(60, 60, 60)

spe_sphere_ception <- simulate_spe_metadata3D(metadata_sphere_ception)
plot_cells3D(spe_sphere_ception,
             plot_cell_types = c("Others", "Tumour", "Immune"),
             plot_colours = c("lightgray", "orange", "skyblue"))


# Separate spheres
metadata_separate_spheres <- spe_metadata_cluster_template(metadata_bg_r, "regular", "Sphere")
metadata_separate_spheres <- spe_metadata_cluster_template(metadata_separate_spheres, "regular", "Sphere")

metadata_separate_spheres[["cluster_1"]][["cluster_cell_types"]] <- c("Tumour")
metadata_separate_spheres[["cluster_1"]][["cluster_cell_proportions"]] <- c(1)
metadata_separate_spheres[["cluster_1"]][["radius"]] <- 30
metadata_separate_spheres[["cluster_1"]][["centre_loc"]] <- c(30, 30, 30)

metadata_separate_spheres[["cluster_2"]][["cluster_cell_types"]] <- c("Immune")
metadata_separate_spheres[["cluster_2"]][["cluster_cell_proportions"]] <- c(1)
metadata_separate_spheres[["cluster_2"]][["radius"]] <- 30
metadata_separate_spheres[["cluster_2"]][["centre_loc"]] <- c(70, 70, 70)

spe_separate_spheres <- simulate_spe_metadata3D(metadata_separate_spheres)
plot_cells3D(spe_separate_spheres,
             plot_cell_types = c("Others", "Tumour", "Immune"),
             plot_colours = c("lightgray", "orange", "skyblue"))



### Analysis - mixed sphere --------------------------------------------------
entropy_grid_metrics <- determine_entropy_grid_metrics3D(spe_mixed_sphere,
                                                         n_split = 8,
                                                         cell_types_of_interest = c("Tumour", "Immune"),
                                                         plot_image = TRUE)
plot_grid_metrics_discrete3D(spe_mixed_sphere, entropy_grid_metrics, "entropy")



cell_proportion_grid_metrics <- determine_cell_proportion_grid_metrics3D(spe_mixed_sphere,
                                                                         n_split = 8,
                                                                         reference_cell_types = c("Tumour"),
                                                                         target_cell_types = c("Immune"),
                                                                         plot_image = TRUE)
plot_grid_metrics_discrete3D(spe_mixed_sphere, cell_proportion_grid_metrics, "proportion")




### Analysis - sphere within sphere --------------------------------------------------
entropy_grid_metrics <- determine_entropy_grid_metrics3D(spe_sphere_ception,
                                                         n_split = 8,
                                                         cell_types_of_interest = c("Tumour", "Immune"),
                                                         plot_image = TRUE)
plot_grid_metrics_discrete3D(spe_sphere_ception, entropy_grid_metrics, "entropy")



cell_proportion_grid_metrics <- determine_cell_proportion_grid_metrics3D(spe_sphere_ception,
                                                                         n_split = 8,
                                                                         reference_cell_types = c("Tumour"),
                                                                         target_cell_types = c("Immune"),
                                                                         plot_image = TRUE)
plot_grid_metrics_discrete3D(spe_sphere_ception, cell_proportion_grid_metrics, "proportion")




### Analysis - separate spheres  --------------------------------------------------
entropy_grid_metrics <- determine_entropy_grid_metrics3D(spe_separate_spheres,
                                                         n_split = 8,
                                                         cell_types_of_interest = c("Tumour", "Immune"),
                                                         plot_image = TRUE)
plot_grid_metrics_discrete3D(spe_separate_spheres, entropy_grid_metrics, "entropy")



cell_proportion_grid_metrics <- determine_cell_proportion_grid_metrics3D(spe_separate_spheres,
                                                                         n_split = 8,
                                                                         reference_cell_types = c("Tumour"),
                                                                         target_cell_types = c("Immune"),
                                                                         plot_image = TRUE)
plot_grid_metrics_discrete3D(spe_separate_spheres, cell_proportion_grid_metrics, "proportion")
