library(alphashape3d)

# Separate ellipsoids
metadata_bg_r <- spe_metadata_background_template("random")
metadata_bg_r[["background"]][["n_cells"]] <- 20000
metadata_bg_r[["background"]][["minimum_distance_between_cells"]] <- 0
metadata_bg_r[["background"]][["cell_types"]] <- "Others"
metadata_bg_r[["background"]][["cell_proportions"]] <- 1

metadata_separate_ellipsoids <- spe_metadata_cluster_template(metadata_bg_r, "regular", "Ellipsoid")
metadata_separate_ellipsoids <- spe_metadata_cluster_template(metadata_separate_ellipsoids, "regular", "Ellipsoid")

metadata_separate_ellipsoids[["cluster_1"]][["cluster_cell_types"]] <- c("Tumour")
metadata_separate_ellipsoids[["cluster_1"]][["cluster_cell_proportions"]] <- c(1)
metadata_separate_ellipsoids[["cluster_1"]][["centre_loc"]] <- c(30, 30, 30)

metadata_separate_ellipsoids[["cluster_2"]][["cluster_cell_types"]] <- c("Immune")
metadata_separate_ellipsoids[["cluster_2"]][["cluster_cell_proportions"]] <- c(1)
metadata_separate_ellipsoids[["cluster_2"]][["centre_loc"]] <- c(70, 70, 70)

spe_separate_ellipsoids <- simulate_spe_metadata3D(metadata_separate_ellipsoids)
plot_cells3D(spe_separate_ellipsoids,
             plot_cell_types = c("Others", "Tumour", "Immune"),
             plot_colours = c("lightgray", "orange", "skyblue"))



spe1 <- spe_separate_ellipsoids
df1 <- data.frame(spatialCoords(spe1), "Cell.Type" = spe1$Cell.Type, "Cell.ID" = spe1$Cell.ID)
df1$Cell.Type <- factor(df1$Cell.Type, c("Others", "Tumour", "Immune"))

df1_tumour_immune_coords <- df1[df1$Cell.Type %in% c("Tumour", "Immune"), c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]

## Get alpha hull
alpha_hull <- ashape3d(as.matrix(df1_tumour_immune_coords), alpha = 15)
vertices <- alpha_hull$x
faces <- alpha_hull$triang[alpha_hull$triang[, 9] %in% c(1, 2, 3), 1:3]


plot_ly() %>%
  add_trace(
    data = df1,
    type = "scatter3d",
    mode = 'markers',
    x = ~Cell.X.Position,
    y = ~Cell.Y.Position,
    z = ~Cell.Z.Position,
    marker = list(size = 2),
    color = ~Cell.Type,
    colors = c("lightgray", "orange", "skyblue")
  ) %>%
  add_trace(
    type = 'mesh3d',
    x = vertices[, 1], 
    y = vertices[, 2], 
    z = vertices[, 3],
    i = faces[, 1] - 1, 
    j = faces[, 2] - 1, 
    k = faces[, 3] - 1,
    opacity = 0.05,
    facecolor = rep("red", nrow(faces))
  )



## Get clusters from alpha hull
alpha_hull_clusters <- components_ashape3d(alpha_hull)
df1_tumour_immune_coords$cluster_number <- paste("Cluster_", alpha_hull_clusters, sep = "")

plot_ly() %>%
  add_trace(
    data = df1_tumour_immune_coords,
    type = "scatter3d",
    mode = 'markers',
    x = ~Cell.X.Position,
    y = ~Cell.Y.Position,
    z = ~Cell.Z.Position,
    marker = list(size = 2),
    color = ~cluster_number,
    colors = c("tomato", "skyblue")
  ) %>%
  add_trace(
    type = 'mesh3d',
    x = vertices[, 1], 
    y = vertices[, 2], 
    z = vertices[, 3],
    i = faces[, 1] - 1, 
    j = faces[, 2] - 1, 
    k = faces[, 3] - 1,
    opacity = 0.05,
    facecolor = rep("yellow", nrow(faces))
  )



## Check if points are inside alpha hull
points_inside_alpha_hull <- inashape3d(alpha_hull, indexAlpha = 1, points = matrix(spatialCoords(spe1)))

df1$points_inside <- points_inside_alpha_hull

plot_ly() %>%
  add_trace(
    data = df1,
    type = "scatter3d",
    mode = 'markers',
    x = ~Cell.X.Position,
    y = ~Cell.Y.Position,
    z = ~Cell.Z.Position,
    marker = list(size = 2),
    color = ~points_inside,
    colors = c("tomato", "skyblue")
  ) %>%
  add_trace(
    type = 'mesh3d',
    x = vertices[, 1], 
    y = vertices[, 2], 
    z = vertices[, 3],
    i = faces[, 1] - 1, 
    j = faces[, 2] - 1, 
    k = faces[, 3] - 1,
    opacity = 0.05,
    facecolor = rep("yellow", nrow(faces))
  )




## Get volume of alpha hull
alpha_hull_volume <- volume_ashape3d(alpha_hull)

expected_volume <- (4/3) * pi * 15 * 20 * 25
