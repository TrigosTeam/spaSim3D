library(plotly)
library(cxhull)

spe_bg <- spaSim3D_background_integrator()
spe_cluster <- spaSim3D_cluster_integrator(spe_bg)
spe1 <- spe_cluster
df1 <- data.frame(spatialCoords(spe1), "Cell.Type" = spe1$Cell.Type, "Cell.ID" = spe1$Cell.ID)
df1$Cell.Type <- factor(df1$Cell.Type, c("Others", "Tumour"))

df1_tumour_coords <- df1[df1$Cell.Type == "Tumour", c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
hull <- cxhull(as.matrix(df1_tumour_coords))
mesh <- hullMesh(hull)
vertices <- mesh$vertices
faces <- mesh$faces


plot_ly() %>%
  add_trace(
    type = "scatter3d",
    mode = 'markers',
    x = df1$Cell.X.Position,
    y = df1$Cell.Y.Position,
    z = df1$Cell.Z.Position,
    marker = list(size = 2),
    color = df1$Cell.Type,
    colors = c("lightgray", "orange")
  ) %>%
  add_trace(
    type = 'mesh3d',
    x = vertices[, 1], 
    y = vertices[, 2], 
    z = vertices[, 3],
    i = faces[, 1] - 1, 
    j = faces[, 2] - 1, 
    k = faces[, 3] - 1,
    opacity = 0.2,
    facecolor = rep("red", nrow(faces))
  )




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
    colors = c("lightgray", "orange")
  ) %>%
  add_trace(
    type = 'mesh3d',
    x = vertices[, 1], 
    y = vertices[, 2], 
    z = vertices[, 3],
    i = faces[, 1] - 1, 
    j = faces[, 2] - 1, 
    k = faces[, 3] - 1,
    opacity = 0.2,
    facecolor = rep("red", nrow(faces))
  )







library(ptinpoly)
point_in_polygon <- pip3d(as.matrix(vertices), as.matrix(faces), as.matrix(spatialCoords(spe1)))

df_temp <- df1
df_temp$in_polygon <- point_in_polygon
df_temp$location <- ifelse(df_temp$in_polygon == -1, "Outside", "Inside")
df_temp$predicted_cell_type <- ifelse(df_temp$in_polygon == -1, "Others", "Tumour")


plot_ly() %>%
  add_trace(
    data = df_temp,
    type = "scatter3d",
    mode = 'markers',
    x = ~Cell.X.Position,
    y = ~Cell.Y.Position,
    z = ~Cell.Z.Position,
    marker = list(size = 2),
    color = ~location,
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
    opacity = 0.2,
    facecolor = rep("yellow", nrow(faces))
  )



df_temp$matched_cell_type <- (df_temp$Cell.Type == df_temp$predicted_cell_type)
plot_ly() %>%
  add_trace(
    data = df_temp,
    type = "scatter3d",
    mode = 'markers',
    x = ~Cell.X.Position,
    y = ~Cell.Y.Position,
    z = ~Cell.Z.Position,
    marker = list(size = 2),
    color = ~matched_cell_type,
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
    opacity = 0.2,
    facecolor = rep("yellow", nrow(faces))
  )
