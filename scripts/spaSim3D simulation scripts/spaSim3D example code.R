###-------------------------------------------------------------------------###
### Background
###-------------------------------------------------------------------------###
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
bg$Cell.ID <- (paste("Cell_", seq(nrow(bg)), sep="")) ## adding Cell.ID column


###-------------------------------------------------------------------------###
### Mixing
###-------------------------------------------------------------------------###
bg_mix <- simulate_mixing3D(bg)



###-------------------------------------------------------------------------###
### Sphere
###-------------------------------------------------------------------------###
bg_sphere <- simulate_clusters3D(bg_sample = bg)
bg_sphere$Cell.ID <- (paste("Cell_", seq(nrow(bg_sphere)), sep="")) ## adding Cell.ID column


###-------------------------------------------------------------------------###
### Two separate spheres
###-------------------------------------------------------------------------###
bg_cluster <- simulate_clusters3D(bg_sample = bg,
                                  n_clusters = 2,
                                  cluster_properties = list(
                                    C1 = list(
                                      name_of_cluster_cell = "Tumour",
                                      infiltration_types = NULL,
                                      infiltration_proportions = NULL,
                                      shape = "Sphere",
                                      radius = 20,
                                      centre_loc = c(30, 30, 30)
                                    ),
                                    C2 = list(
                                      name_of_cluster_cell = "Immune",
                                      infiltration_types = NULL,
                                      infiltration_proportions = NULL,
                                      shape = "Sphere",
                                      radius = 20,
                                      centre_loc = c(70, 70, 70)
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_categories = c("Others", "Tumour", "Immune"),
                                  plot_colours = NULL)
bg_cluster$Cell.ID <- (paste("Cell_", seq(nrow(bg_cluster)), sep="")) ## adding Cell.ID column



###-------------------------------------------------------------------------###
### One sphere with mixing
###-------------------------------------------------------------------------###
bg_cluster <- simulate_clusters3D(bg_sample = bg,
                                  n_clusters = 1,
                                  cluster_properties = list(
                                    C1 = list(
                                      name_of_cluster_cell = "Tumour",
                                      infiltration_types = "Immune",
                                      infiltration_proportions = 0.5,
                                      shape = "Sphere",
                                      radius = 30,
                                      centre_loc = c(50, 50, 50)
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_categories = c("Others", "Tumour", "Immune"),
                                  plot_colours = NULL)
bg_cluster$Cell.ID <- (paste("Cell_", seq(nrow(bg_cluster)), sep="")) ## adding Cell.ID column




###-------------------------------------------------------------------------###
### Cylinder
###-------------------------------------------------------------------------###
bg_cylinder <- simulate_clusters3D(bg_sample = bg,
                                   n_clusters = 4,
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



###-------------------------------------------------------------------------###
### Ellipsoid
###-------------------------------------------------------------------------###
bg_ellipsoid <- simulate_clusters3D(bg_sample = bg,
                                   n_clusters = 2,
                                   cluster_properties = list(
                                     C1 = list(
                                       name_of_cluster_cell = "Tumour",
                                       infiltration_types = c("Immune1", "Others"),
                                       infiltration_proportions = c(0.1, 0.05),
                                       shape = "Ellipsoid",
                                       x_radius = 20,
                                       y_radius = 20,
                                       z_radius = 30,
                                       centre_loc = c(50, 50, 50),
                                       y_z_rotation = pi/4,
                                       x_z_rotation = 0,
                                       x_y_rotation = 0
                                     ), 
                                     C2 = list(
                                       name_of_cluster_cell = "Tumour",
                                       infiltration_types = c("Immune1", "Others"),
                                       infiltration_proportions = c(0.1, 0.05),
                                       shape = "Ellipsoid",
                                       x_radius = 20,
                                       y_radius = 20,
                                       z_radius = 30,
                                       centre_loc = c(50, 66, 50),
                                       y_z_rotation = -pi/4,
                                       x_z_rotation = 0,
                                       x_y_rotation = 0
                                     )
                                   ),
                                   plot_image = TRUE,
                                   plot_categories = c("Others", "Immune1", "Tumour"),
                                   plot_colours = NULL)



###-------------------------------------------------------------------------###
### Heart
###-------------------------------------------------------------------------###
bg_heart <- bg_ellipsoid[bg_ellipsoid$Cell.Type != "Others" ,]
plot3d(bg_heart$Cell.X.Position, 
       bg_heart$Cell.Y.Position, 
       bg_heart$Cell.Z.Position, 
       col = "red", 
       size = 4,
       xlim = c(0, 100),
       ylim = c(0, 100),
       zlim = c(0, 100))



###-------------------------------------------------------------------------###
### Ring
###-------------------------------------------------------------------------###
bg_ring <- simulate_rings3D(bg_sample = bg)



###-------------------------------------------------------------------------###
### Heart with ring
###-------------------------------------------------------------------------###
bg_heart_ring <- simulate_rings3D(bg_sample = bg,
                                     n_ring = 2,
                                     ring_properties = list(
                                       R1 = list(
                                          name_of_cluster_cell = "Tumour",
                                          infiltration_types = c("Immune1", "Others"),
                                          infiltration_proportions = c(0.1, 0.05),
                                          shape = "Ellipsoid",
                                          x_radius = 15,
                                          y_radius = 15,
                                          z_radius = 25,
                                          centre_loc = c(50, 50, 50),
                                          y_z_rotation = pi/4,
                                          x_z_rotation = 0,
                                          x_y_rotation = 0,
                                          name_of_ring_cell = "Immune1",
                                          ring_width = 5,
                                          ring_infiltration_types = c("Others"),
                                          ring_infiltration_proportions = c(0.15)
                                       ), 
                                       R2 = list(
                                          name_of_cluster_cell = "Tumour",
                                          infiltration_types = c("Immune1", "Others"),
                                          infiltration_proportions = c(0.1, 0.05),
                                          shape = "Ellipsoid",
                                          x_radius = 15,
                                          y_radius = 15,
                                          z_radius = 25,
                                          centre_loc = c(50, 66, 50),
                                          y_z_rotation = -pi/4,
                                          x_z_rotation = 0,
                                          x_y_rotation = 0,
                                          name_of_ring_cell = "Immune1",
                                          ring_width = 5,
                                          ring_infiltration_types = c("Others"),
                                          ring_infiltration_proportions = c(0.15)
                                       )
                                     ),
                                     plot_image = TRUE,
                                     plot_categories = c("Others", "Immune1", "Tumour"),
                                     plot_colours = NULL)



###-------------------------------------------------------------------------###
### Hollow cylinder
###-------------------------------------------------------------------------###
bg_cylinder_ring <- simulate_rings3D(bg_sample = bg,
                                  n_ring = 1,
                                  ring_properties = list(
                                    R1 = list(
                                      name_of_cluster_cell = "Void",
                                      infiltration_types = c("Immune1", "Others"),
                                      infiltration_proportions = c(0, 0),
                                      shape = "Cylinder",
                                      radius = 10,
                                      start_loc = c(0, 0, 0),
                                      end_loc = c(100, 80, 60),
                                      name_of_ring_cell = "Immune1",
                                      ring_width = 5,
                                      ring_infiltration_types = c("Others"),
                                      ring_infiltration_proportions = c(0.15)
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_categories = c("Others", "Immune1", "Tumour"),
                                  plot_colours = NULL)



###-------------------------------------------------------------------------###
### Double rings
###-------------------------------------------------------------------------###
bg_dr <- simulate_double_rings3D(bg_sample = bg)




###-------------------------------------------------------------------------###
### Heart with double rings
###-------------------------------------------------------------------------###
bg_heart_dr <- simulate_double_rings3D(bg_sample = bg,
                                       n_dr = 2,
                                       dr_properties = list(
                                       D1 = list(
                                        name_of_cluster_cell = "Tumour",
                                        infiltration_types = c("Immune1", "Others"),
                                        infiltration_proportions = c(0.1, 0.05),
                                        shape = "Ellipsoid",
                                        x_radius = 15,
                                        y_radius = 15,
                                        z_radius = 25,
                                        centre_loc = c(50, 50, 50),
                                        y_z_rotation = pi/4,
                                        x_z_rotation = 0,
                                        x_y_rotation = 0,
                                        name_of_inner_ring_cell = "Immune1",
                                        inner_ring_width = 3,
                                        inner_ring_infiltration_types = c("Others"),
                                        inner_ring_infiltration_proportions = c(0.15),
                                        name_of_outer_ring_cell = "Immune2",
                                        outer_ring_width = 2,
                                        outer_ring_infiltration_types = c("Others"),
                                        outer_ring_infiltration_proportions = c(0.15)
                                       ), 
                                       D2 = list(
                                        name_of_cluster_cell = "Tumour",
                                        infiltration_types = c("Immune1", "Others"),
                                        infiltration_proportions = c(0.1, 0.05),
                                        shape = "Ellipsoid",
                                        x_radius = 15,
                                        y_radius = 15,
                                        z_radius = 25,
                                        centre_loc = c(50, 66, 50),
                                        y_z_rotation = -pi/4,
                                        x_z_rotation = 0,
                                        x_y_rotation = 0,
                                        name_of_inner_ring_cell = "Immune1",
                                        inner_ring_width = 3,
                                        inner_ring_infiltration_types = c("Others"),
                                        inner_ring_infiltration_proportions = c(0.15),
                                        name_of_outer_ring_cell = "Immune2",
                                        outer_ring_width = 2,
                                        outer_ring_infiltration_types = c("Others"),
                                        outer_ring_infiltration_proportions = c(0.15)
                                       )
                                     ),
                                     plot_image = TRUE,
                                     plot_categories = c("Others", "Tumour", "Immune1", "Immune2"),
                                     plot_colours = NULL)






###-------------------------------------------------------------------------###
### Plotting a slice of 3D data in 2D
###-------------------------------------------------------------------------###
bg_chosen <- bg_sphere
bg_slice <- bg_chosen[bg_chosen$Cell.Z.Position > 45 & 
                      bg_chosen$Cell.Z.Position < 55, ]

library(ggplot2)

ggplot(bg_slice,
       aes(Cell.X.Position, Cell.Y.Position, color = Cell.Type)) +
  geom_point() + 
  scale_color_manual(values=c("red", "green", "grey", "blue"))
