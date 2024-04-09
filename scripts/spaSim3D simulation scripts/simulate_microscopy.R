# In 3D
library(rgl)

bg <- simulate_background_cells3D(n_cells = 20000,
                                  length = 100,
                                  width  = 100,
                                  height = 100,
                                  method = "tumour",
                                  min_d = 0.5,
                                  oversampling_rate = 1.2,
                                  jitter_prop = 0,
                                  cell_type = "Others",
                                  plot_image = T)

bg_mix <- simulate_mixing3D(bg,
                            idents = c("Others", "Immune1"),
                            props = c(0.95, 0.05),
                            plot_image = TRUE,
                            plot_colours = c("#0077B6", "lightgreen"))



bg_cylinder_ring <- simulate_rings3D(bg_sample = bg_mix,
                                     n_ring = 4,
                                     ring_properties = list(
                                       R1 = list(
                                         name_of_cluster_cell = "Void",
                                         infiltration_types = NULL,
                                         infiltration_proportions = NULL,
                                         shape = "Cylinder",
                                         radius = 10,
                                         start_loc = c(30, 20, 0),
                                         end_loc = c(30, 60, 40),
                                         name_of_ring_cell = "Connective",
                                         ring_width = 5,
                                         ring_infiltration_types = c("Others"),
                                         ring_infiltration_proportions = c(0.30)
                                       ),
                                       R2 = list(
                                         name_of_cluster_cell = "Void",
                                         infiltration_types = NULL,
                                         infiltration_proportions = NULL,
                                         shape = "Cylinder",
                                         radius = 10,
                                         start_loc = c(30, 60, 40),
                                         end_loc = c(30, 60, 100),
                                         name_of_ring_cell = "Connective",
                                         ring_width = 5,
                                         ring_infiltration_types = c("Others"),
                                         ring_infiltration_proportions = c(0.30)
                                       ),
                                       R3 = list(
                                         name_of_cluster_cell = "Void",
                                         infiltration_types = NULL,
                                         infiltration_proportions = NULL,
                                         shape = "Cylinder",
                                         radius = 8,
                                         start_loc = c(30, 70, 40),
                                         end_loc = c(30, 70, 100),
                                         name_of_ring_cell = "Connective",
                                         ring_width = 5,
                                         ring_infiltration_types = c("Others"),
                                         ring_infiltration_proportions = c(0.30)
                                       ),
                                       R4 = list(
                                         name_of_cluster_cell = "Void",
                                         infiltration_types = NULL,
                                         infiltration_proportions = NULL,
                                         shape = "Cylinder",
                                         radius = 8,
                                         start_loc = c(30, 70, 40),
                                         end_loc = c(30, 100, 0),
                                         name_of_ring_cell = "Connective",
                                         ring_width = 5,
                                         ring_infiltration_types = c("Others"),
                                         ring_infiltration_proportions = c(0.30)
                                       )
                                     ),
                                     plot_image = TRUE,
                                     plot_categories = c("Others", "Immune1", "Immune2", "Connective", "Tumour1", "Tumour2"),
                                     plot_colours = c("#0077B6", "lightgreen", "green", "blue", "darkblue", "red"))



bg_ring <- simulate_rings3D(bg_sample = bg_cylinder_ring,
                                  n_ring = 1,
                                  ring_properties = list(
                                    R1 = list(
                                      name_of_cluster_cell = "Tumour1",
                                      infiltration_types = c("Immune1", "Others"),
                                      infiltration_proportions = c(0.1, 0.05),
                                      shape = "Ellipsoid",
                                      x_radius = 18,
                                      y_radius = 18,
                                      z_radius = 25,
                                      centre_loc = c(80, 40, 50),
                                      y_z_rotation = -pi/2,
                                      x_z_rotation = 0,
                                      x_y_rotation = -pi/4,
                                      name_of_ring_cell = "Immune2",
                                      ring_width = 2,
                                      ring_infiltration_types = c("Others"),
                                      ring_infiltration_proportions = c(0.15)
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_categories = c("Others", "Immune1", "Immune2", "Connective", "Tumour1", "Tumour2"),
                                  plot_colours = c("#0077B6", "lightgreen", "green", "blue", "darkblue", "red"))


bg_cluster <- simulate_clusters3D(bg_sample = bg_ring,
                                    n_clusters = 1,
                                    cluster_properties = list(
                                      C1 = list(
                                          name_of_cluster_cell = "Tumour2",
                                          infiltration_types = c("Immune1", "Tumour1"),
                                          infiltration_proportions = c(0.20, 0.20),
                                          shape = "Sphere",
                                          radius = 15,
                                          centre_loc = c(85, 45, 50)
                                        )
                                    ),
                                    plot_image = TRUE,
                                    plot_categories = c("Others", "Immune1", "Immune2", "Connective", "Tumour1", "Tumour2"),
                                    plot_colours = c("#0077B6", "lightgreen", "green", "blue", "darkblue", "red"))


## Get 2D slice from 3D data (between z = 47.5 and z = 52.5)
df_2d <- bg_cluster[bg_cluster$Cell.Z.Position >= 45 &
                    bg_cluster$Cell.Z.Position <= 55, ]

library(ggplot2)

ggplot(data = df_2d,
       aes(x = Cell.X.Position, y = Cell.Y.Position, color = Cell.Type)) + 
  geom_point() + 
  scale_color_manual(values=c("blue", "lightgreen", "green", "#0077B6","darkblue", "red"))




## ---------------------------------------------------------------------------
# ## In 2D
# 
# library(spaSim)
# 
# bg <- simulate_background_cells(n_cells = 5000,
#                                 width = 1600,
#                                 height = 1200,
#                                 method = "Hardcore",
#                                 min_d = 10,
#                                 #oversampling_rate,
#                                 #jitter,
#                                 #Cell.Type,
#                                 plot_image = T)
# 
# bgMix <- simulate_mixing(bg_sample = bg,
#                          idents = c("Others", "Immune1"),
#                          props = c(0.9, 0.1),
#                          plot_image = T,
#                          plot_colours = c("skyblue", "green"))
# 
# 
# 
# cluster_prop <- list(
#   C1 = list(
#     name_of_cluster_cell = "Immune2",
#     size = 200,
#     shape = "Oval",
#     centre_loc = data.frame("x" = 1100, "y" = 400),
#     infiltration_types = c("Immune1", "Tumour"),
#     infiltration_proportions = c(0.10, 0)),
#   C2 = list(
#     name_of_cluster_cell = "Tumour",
#     size = 150,
#     shape = "Circle",
#     centre_loc = data.frame("x" = 1000, "y" = 500),
#     infiltration_types = c("Immune1", "Immune2"),
#     infiltration_proportions = c(0.10, 0.35))
# )
# 
# 
# bgCluster <- simulate_clusters(bg_sample = bgMix,
#                                n_clusters = 2,
#                                cluster_properties = cluster_prop,
#                                plot_image = T,
#                                plot_categories = c("Others", "Immune1", "Immune2", "Tumour"),
#                                plot_colours = c("skyblue", "green", "blue", "red"))
# 
# 
# 
# ir_prop <- list(I1 = list(name_of_cluster_cell = "Void", size = 150, shape = "Circle", 
#                           centre_loc = data.frame(x = 500, y = 700), 
#                           infiltration_types = c("Immune1"), 
#                           infiltration_proportions = c(0),
#                           name_of_ring_cell = "Immune2", immune_ring_width = 40,
#                           immune_ring_infiltration_types = c("Others"), 
#                           immune_ring_infiltration_proportions = c(0)))
# 
# bgRings <- simulate_immune_rings(bg_sample = bgCluster,
#                                  n_ir = 1,
#                                  ir_properties = ir_prop,
#                                  plot_image = T,
#                                  plot_categories = c("Others", "Immune1", "Immune2", "Tumour", "Void"),
#                                  plot_colours = c("skyblue", "green", "blue", "red", "black"))
# 
# 
