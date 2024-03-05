
bg <- simulate_background_cells3D(n_cells = 10000,
                                  length = 100,
                                  width  = 100,
                                  height = 100,
                                  method = "tumour",
                                  min_d = 2,
                                  oversampling_rate = 1.2,
                                  jitter_prop = 0,
                                  cell_type = "Others",
                                  plot_image = T)

# color <- ifelse(bg$Cell.Z.Position == bg$Cell.Z.Position[1],
#                 "blue",
#                 ifelse(bg$Cell.Z.Position == bg$Cell.Z.Position[500],
#                 "red",
#                 "lightgray"))
#
# plot3d(bg$Cell.X.Position,
#        bg$Cell.Y.Position,
#        bg$Cell.Z.Position,
#        xlab = "x",
#        ylab = "y",
#        zlab = "z",
#        col = color,
#        size = 4)

bg_mix <- simulate_mixing3D(bg)

bg_sphere <- simulate_clusters3D(bg_sample = bg)


bg_cylinder <- simulate_clusters3D(bg_sample = bg,
                                   n_clusters = 4,
                                   bg_type = "Others",
                                   cluster_properties = list(
                                     C1 = list(
                                       name_of_cluster_cell = "Endo",
                                       infiltration_types = c("Immune1", "Others"),
                                       infiltration_proportions = c(0.1, 0.05),
                                       shape = "Cylinder",
                                       radius = 8,
                                       start_loc = c(0, 0, 0),
                                       end_loc   = c(30, 30 ,30)
                                     ),
                                     C2 = list(
                                       name_of_cluster_cell = "Endo",
                                       infiltration_types = c("Immune1", "Others"),
                                       infiltration_proportions = c(0.1, 0.05),
                                       shape = "Cylinder",
                                       radius = 5,
                                       start_loc = c(30, 30, 30),
                                       end_loc   = c(30, 30 ,70)
                                     ),
                                     C3 = list(
                                       name_of_cluster_cell = "Endo",
                                       infiltration_types = c("Immune1", "Others"),
                                       infiltration_proportions = c(0.1, 0.05),
                                       shape = "Cylinder",
                                       radius = 5,
                                       start_loc = c(30, 30, 30),
                                       end_loc   = c(60, 40 ,30)
                                     ),
                                     C4 = list(
                                       name_of_cluster_cell = "Endo",
                                       infiltration_types = c("Immune1", "Others"),
                                       infiltration_proportions = c(0.1, 0.05),
                                       shape = "Cylinder",
                                       radius = 4,
                                       start_loc = c(60, 40, 30),
                                       end_loc   = c(100, 100 ,100)
                                     )
                                   ),
                                   plot_image = TRUE,
                                   plot_categories = c("Others", "Immune1", "Endo"),
                                   plot_colours = NULL)

bg_cluster <- simulate_clusters3D(bg_sample = bg, plot_categories = c("Others", "Tumour", "Endo", "Immune1"))




start_loc <- c(0, 0, 0)

end_loc <- c(50, 50, 50)

v1 <- end_loc - start_loc
d <- sum(v1 * end_loc)

color <- c()

for (i in 1:1000) {
  point <- c(x[i], y[i], z[i])
  
  if (sum(v1 * point) >= d) {
    color <- append(color, "blue")
  }
  else {
    color <- append(color, "lightgray")
  }
}

plot3d(x, y, z, col = color, size = 4)
