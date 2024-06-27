### Background ---------------------------------------------------------------
md_bg_r <- spe_metadata_background_template("random")
md_bg_r[["background"]][["n_cells"]] <- 15000
md_bg_r[["background"]][["cell_types"]] <- c("Others", "Tumour", "Immune")
md_bg_r[["background"]][["cell_proportions"]] <- c(0.96, 0.03, 0.01)
md_bg_r[["background"]][["minimum_distance_between_cells"]] <- 0


# Mixed sphere ---------------------------------------------------------------
md_mixed_sphere <- spe_metadata_cluster_template(md_bg_r, "regular", "Sphere")
md_mixed_sphere[["cluster_1"]][["cluster_cell_types"]] <- c("Tumour", "Immune")
md_mixed_sphere[["cluster_1"]][["cluster_cell_proportions"]] <- c(0.5, 0.5)
md_mixed_sphere[["cluster_1"]][["radius"]] <- 40
md_mixed_sphere[["cluster_1"]][["centre_loc"]] <- c(50, 50, 50)

spe_mixed_sphere <- simulate_spe_metadata3D(md_mixed_sphere)
plot_cells3D(spe_mixed_sphere,
             plot_cell_types = c("Others", "Tumour", "Immune"),
             plot_colours = c("lightgray", "orange", "skyblue"))



# Sphere within sphere -------------------------------------------------------
md_sphere_ception <- spe_metadata_cluster_template(md_bg_r, "regular", "Sphere")
md_sphere_ception <- spe_metadata_cluster_template(md_sphere_ception, "regular", "Sphere")

md_sphere_ception[["cluster_1"]][["cluster_cell_types"]] <- c("Tumour")
md_sphere_ception[["cluster_1"]][["cluster_cell_proportions"]] <- c(1)
md_sphere_ception[["cluster_1"]][["radius"]] <- 40
md_sphere_ception[["cluster_1"]][["centre_loc"]] <- c(50, 50, 50)

md_sphere_ception[["cluster_2"]][["cluster_cell_types"]] <- c("Immune")
md_sphere_ception[["cluster_2"]][["cluster_cell_proportions"]] <- c(1)
md_sphere_ception[["cluster_2"]][["radius"]] <- 25
md_sphere_ception[["cluster_2"]][["centre_loc"]] <- c(60, 60, 60)

spe_sphere_ception <- simulate_spe_metadata3D(md_sphere_ception)
plot_cells3D(spe_sphere_ception,
             plot_cell_types = c("Others", "Tumour", "Immune"),
             plot_colours = c("lightgray", "orange", "skyblue"))


# Separate spheres -----------------------------------------------------------
md_separate_spheres <- spe_metadata_cluster_template(md_bg_r, "regular", "Sphere")
md_separate_spheres <- spe_metadata_cluster_template(md_separate_spheres, "regular", "Sphere")

md_separate_spheres[["cluster_1"]][["cluster_cell_types"]] <- c("Tumour")
md_separate_spheres[["cluster_1"]][["cluster_cell_proportions"]] <- c(1)
md_separate_spheres[["cluster_1"]][["radius"]] <- 30
md_separate_spheres[["cluster_1"]][["centre_loc"]] <- c(30, 30, 30)

md_separate_spheres[["cluster_2"]][["cluster_cell_types"]] <- c("Immune")
md_separate_spheres[["cluster_2"]][["cluster_cell_proportions"]] <- c(1)
md_separate_spheres[["cluster_2"]][["radius"]] <- 30
md_separate_spheres[["cluster_2"]][["centre_loc"]] <- c(70, 70, 70)

spe_separate_spheres <- simulate_spe_metadata3D(md_separate_spheres)
plot_cells3D(spe_separate_spheres,
             plot_cell_types = c("Others", "Tumour", "Immune"),
             plot_colours = c("lightgray", "orange", "skyblue"))


# Ringed sphere --------------------------------------------------------------
md_ringed_sphere <- spe_metadata_cluster_template(md_bg_r, "ring", "Sphere")
md_ringed_sphere[["cluster_1"]][["cluster_cell_types"]] <- c("Tumour")
md_ringed_sphere[["cluster_1"]][["cluster_cell_proportions"]] <- c(1)
md_ringed_sphere[["cluster_1"]][["radius"]] <- 40
md_ringed_sphere[["cluster_1"]][["centre_loc"]] <- c(50, 50, 50)
md_ringed_sphere[["cluster_1"]][["ring_cell_types"]] <- "Immune"
md_ringed_sphere[["cluster_1"]][["ring_cell_proportions"]] <- 1

spe_ringed_sphere <- simulate_spe_metadata3D(md_ringed_sphere)
plot_cells3D(spe_ringed_sphere,
             plot_cell_types = c("Others", "Tumour", "Immune"),
             plot_colours = c("lightgray", "orange", "skyblue"))



# Ringed sphere and non-ringed sphere ----------------------------------------
md_ringed_and_non_ringed_spheres <- spe_metadata_cluster_template(md_bg_r, "regular", "Sphere")
md_ringed_and_non_ringed_spheres <- spe_metadata_cluster_template(md_ringed_and_non_ringed_spheres, "ring", "Sphere")

md_ringed_and_non_ringed_spheres[["cluster_1"]][["cluster_cell_types"]] <- c("Tumour", "Immune")
md_ringed_and_non_ringed_spheres[["cluster_1"]][["cluster_cell_proportions"]] <- c(0.5, 0.5)
md_ringed_and_non_ringed_spheres[["cluster_1"]][["radius"]] <- 20
md_ringed_and_non_ringed_spheres[["cluster_1"]][["centre_loc"]] <- c(30, 30, 30)

md_ringed_and_non_ringed_spheres[["cluster_2"]][["cluster_cell_types"]] <- c("Tumour", "Immune1")
md_ringed_and_non_ringed_spheres[["cluster_2"]][["cluster_cell_proportions"]] <- c(0.9, 0.1)
md_ringed_and_non_ringed_spheres[["cluster_2"]][["radius"]] <- 20
md_ringed_and_non_ringed_spheres[["cluster_2"]][["centre_loc"]] <- c(70, 70, 70)
md_ringed_and_non_ringed_spheres[["cluster_2"]][["ring_cell_types"]] <- c("Immune1")
md_ringed_and_non_ringed_spheres[["cluster_2"]][["ring_cell_proportions"]] <- c(1)
md_ringed_and_non_ringed_spheres[["cluster_2"]][["ring_width"]] <- 5

ringed_and_non_ringed_spheres <- simulate_spe_metadata3D(md_ringed_and_non_ringed_spheres)
plot_cells3D(ringed_and_non_ringed_spheres,
             plot_cell_types = c("Others", "Tumour", "Immune", "Immune1"),
             plot_colours = c("lightgray", "orange", "skyblue", "lightgreen"))
