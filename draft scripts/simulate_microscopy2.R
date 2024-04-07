library(rgl)

bg <- simulate_background_cells3D(n_cells = 20000,
                                  length = 100,
                                  width  = 100,
                                  height = 100,
                                  method = "tumour",
                                  min_d = 0.5,
                                  oversampling_rate = 1.1,
                                  jitter_prop = 0,
                                  cell_type = "Others",
                                  plot_image = T)

bg_mix <- simulate_mixing3D(bg,
                            idents = c("Others", "Immune1"),
                            props = c(0.95, 0.05),
                            plot_image = TRUE,
                            plot_colours = c("#0077B6", "lightgreen"))




bg_ring <- simulate_rings3D(bg_sample = bg_mix,
                            bg_type = "Others",
                            n_ring = 1,
                            ring_properties = list(
                              R1 = list(
                                name_of_cluster_cell = "Tumour",
                                infiltration_types = c("Immune1", "Others"),
                                infiltration_proportions = c(0.1, 0.05),
                                shape = "Ellipsoid",
                                x_radius = 15,
                                y_radius = 10,
                                z_radius = 10,
                                centre_loc = c(40, 30, 50),
                                y_z_rotation = 0,
                                x_z_rotation = 0,
                                x_y_rotation = -pi/4,
                                name_of_ring_cell = "Immune2",
                                ring_width = 2,
                                ring_infiltration_types = c("Others"),
                                ring_infiltration_proportions = c(0.15)
                              )
                            ),
                            plot_image = TRUE,
                            plot_categories = c("Others", "Immune1", "Immune2", "Connective", "Tumour"),
                            plot_colours = c("#0077B6", "lightgreen", "green", "blue", "red"))



bg_cluster <- simulate_clusters3D(bg_sample = bg_ring,
                                  n_clusters = 4,
                                  bg_type = "Others",
                                  cluster_properties = list(
                                    C1 = list(
                                      name_of_cluster_cell = "Tumour",
                                      infiltration_types = c("Immune1", "Others"),
                                      infiltration_proportions = c(0.20, 0.05),
                                      shape = "Sphere",
                                      radius = 10,
                                      centre_loc = c(70, 35, 50)
                                    ),
                                    C2 = list(
                                      name_of_cluster_cell = "Tumour",
                                      infiltration_types = c("Immune1", "Others"),
                                      infiltration_proportions = c(0.20, 0.05),
                                      shape = "Ellipsoid",
                                      x_radius = 15,
                                      y_radius = 10,
                                      z_radius = 10,
                                      centre_loc = c(10, 40, 50),
                                      y_z_rotation = 0,
                                      x_z_rotation = 0,
                                      x_y_rotation = -pi/4
                                    ),
                                    C3 = list(
                                      name_of_cluster_cell = "Tumour",
                                      infiltration_types = c("Immune1", "Others"),
                                      infiltration_proportions = c(0.20, 0.05),
                                      shape = "Sphere",
                                      radius = 8,
                                      centre_loc = c(15, 50, 50)
                                    ),
                                    C4 = list(
                                      name_of_cluster_cell = "Tumour",
                                      infiltration_types = c("Immune1", "Others"),
                                      infiltration_proportions = c(0.20, 0.05),
                                      shape = "Ellipsoid",
                                      x_radius = 20,
                                      y_radius = 15,
                                      z_radius = 15,
                                      centre_loc = c(85, 40, 50),
                                      y_z_rotation = 0,
                                      x_z_rotation = 0,
                                      x_y_rotation = -pi/4
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_categories = c("Others", "Immune1", "Immune2", "Connective", "Tumour"),
                                  plot_colours = c("#0077B6", "lightgreen", "green", "blue", "red"))



bg_cylinder_ring <- simulate_rings3D(bg_sample = bg_cluster,
                                     bg_type = "Others",
                                     n_ring = 4,
                                     ring_properties = list(
                                       R1 = list(
                                         name_of_cluster_cell = "Void",
                                         infiltration_types = NULL,
                                         infiltration_proportions = NULL,
                                         shape = "Cylinder",
                                         radius = 10,
                                         start_loc = c(80, 80, 0),
                                         end_loc = c(80, 80, 100),
                                         name_of_ring_cell = "Connective",
                                         ring_width = 6,
                                         ring_infiltration_types = c("Others"),
                                         ring_infiltration_proportions = c(0.30)
                                       ),
                                       R2 = list(
                                         name_of_cluster_cell = "Void",
                                         infiltration_types = NULL,
                                         infiltration_proportions = NULL,
                                         shape = "Cylinder",
                                         radius = 5,
                                         start_loc = c(90, 90, 0),
                                         end_loc = c(90, 90, 100),
                                         name_of_ring_cell = "Connective",
                                         ring_width = 6,
                                         ring_infiltration_types = c("Others"),
                                         ring_infiltration_proportions = c(0.30)
                                       ),
                                       R3 = list(
                                         name_of_cluster_cell = "Void",
                                         infiltration_types = NULL,
                                         infiltration_proportions = NULL,
                                         shape = "Cylinder",
                                         radius = 6,
                                         start_loc = c(80, 40, 0),
                                         end_loc = c(80, 40, 100),
                                         name_of_ring_cell = "Connective",
                                         ring_width = 3,
                                         ring_infiltration_types = c("Others"),
                                         ring_infiltration_proportions = c(0.25)
                                       ),
                                       R4 = list(
                                         name_of_cluster_cell = "Void",
                                         infiltration_types = NULL,
                                         infiltration_proportions = NULL,
                                         shape = "Cylinder",
                                         radius = 10,
                                         start_loc = c(5, 30, 0),
                                         end_loc = c(5, 30, 100),
                                         name_of_ring_cell = "Connective",
                                         ring_width = 6,
                                         ring_infiltration_types = c("Others"),
                                         ring_infiltration_proportions = c(0.30)
                                       )
                                     ),
                                     plot_image = TRUE,
                                     plot_categories = c("Others", "Immune1", "Immune2", "Connective", "Tumour"),
                                     plot_colours = c("#0077B6", "lightgreen", "green", "blue", "red"))



## Get 2D slice from 3D data (between z = 47.5 and z = 52.5)

data <- bg_cylinder_ring

df_2d <- data[data$Cell.Z.Position >= 45 &
              data$Cell.Z.Position <= 55, ]

library(ggplot2)

ggplot(data = df_2d,
       aes(x = Cell.X.Position, y = Cell.Y.Position, color = Cell.Type)) + 
  geom_point() + 
  scale_color_manual(values=c("blue", "lightgreen", "green", "#0077B6", "red"))
