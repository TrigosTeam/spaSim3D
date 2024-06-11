### Add Cell.ID
chosen_bg <- bg_sphere
chosen_bg$Cell.ID <- (paste("Cell_", seq(nrow(chosen_bg)), sep="")) ## adding Cell.ID column
data <- chosen_bg


###-------------------------------------------------------------------------###
### Background
###-------------------------------------------------------------------------###
bg <- simulate_random_background_cells3D(n_cells = 10000,
                                         length = 100,
                                         width = 100,
                                         height = 100,
                                         minimum_distance_between_cells = 2,
                                         background_cell_type = "Others",
                                         plot_image = TRUE)

bg_normal <- simulate_normal_background_cells3D(n_cells = 3000,
                                                length = 100,
                                                width = 100,
                                                height = 100,
                                                jitter_proportion = 0,
                                                background_cell_type = "Others",
                                                plot_image = TRUE)

###-------------------------------------------------------------------------###
### Mixing
###-------------------------------------------------------------------------###
bg_mix <- simulate_mixing3D(bg,
                            cell_types = c("Others", "Immune", "Tumour"),
                            props = c(0.5, 0.2, 0.3),
                            plot_image = TRUE,
                            plot_categories = c("Others", "Immune", "Tumour"),
                            plot_colours = c("lightgray", "skyblue", "orange"))

###-------------------------------------------------------------------------###
### Clusters
###-------------------------------------------------------------------------###
bg_cluster <- simulate_clusters3D(bg,
                                  n_clusters = 3,
                                  cluster_properties = list(
                                    C1 = list(
                                      name_of_cluster_cell = "Tumour",
                                      infiltration_types = c("Immune", "Others"),
                                      infiltration_proportions = c(0.4, 0.05),
                                      shape = "Sphere",
                                      radius = 25,
                                      centre_loc = c(40, 40, 40)
                                    ),
                                    C2 = list(
                                      name_of_cluster_cell = "Endothelial",
                                      infiltration_types = c("Others"),
                                      infiltration_proportions = c(0.05),
                                      shape = "Cylinder",
                                      radius = 10,
                                      start_loc = c(0, 0, 0),
                                      end_loc   = c(20, 20 , 100)
                                    ),
                                    C3 = list(
                                      name_of_cluster_cell = "Tumour",
                                      infiltration_types = c("Immune", "Others"),
                                      infiltration_proportions = c(0.3, 0.05),
                                      shape = "Ellipsoid",
                                      x_radius = 15,
                                      y_radius = 20,
                                      z_radius = 25,
                                      centre_loc = c(70, 70, 70),
                                      x_y_rotation = 0,
                                      x_z_rotation = 0,
                                      y_z_rotation = 0
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_categories = c("Others", "Immune", "Endothelial", "Tumour"),
                                  plot_colours = c("lightgray", "skyblue", "#FF7F7F", "orange"))

###-------------------------------------------------------------------------###
### Two separate spheres
###-------------------------------------------------------------------------###
bg_spheres <- simulate_clusters3D(bg_sample = bg,
                                  n_clusters = 2,
                                  cluster_properties = list(
                                    C1 = list(
                                      name_of_cluster_cell = "Tumour",
                                      infiltration_types = NULL,
                                      infiltration_proportions = NULL,
                                      shape = "Sphere",
                                      radius = 30,
                                      centre_loc = c(30, 30, 30)
                                    ),
                                    C2 = list(
                                      name_of_cluster_cell = "Immune",
                                      infiltration_types = NULL,
                                      infiltration_proportions = NULL,
                                      shape = "Sphere",
                                      radius = 25,
                                      centre_loc = c(70, 70, 70)
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_categories = c("Others", "Tumour", "Immune"),
                                  plot_colours = c("lightgray", "orange", "skyblue"))

###-------------------------------------------------------------------------###
### One sphere with mixing
###-------------------------------------------------------------------------###
bg_sphere <-  simulate_clusters3D(bg_sample = bg,
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
                                  plot_colours = c("lightgray", "orange", "skyblue"))


###-------------------------------------------------------------------------###
### Cylinder
###-------------------------------------------------------------------------###
bg_cylinder <- simulate_clusters3D(bg_sample = bg,
                                   n_clusters = 4,
                                   cluster_properties = list(
                                     C1 = list(
                                       name_of_cluster_cell = "Endothelial",
                                       infiltration_types = c("Immune", "Others"),
                                       infiltration_proportions = c(0.1, 0.05),
                                       shape = "Cylinder",
                                       radius = 8,
                                       start_loc = c(0, 0, 0),
                                       end_loc   = c(30, 30 ,30)
                                     ),
                                     C2 = list(
                                       name_of_cluster_cell = "Endothelial",
                                       infiltration_types = c("Immune", "Others"),
                                       infiltration_proportions = c(0.1, 0.05),
                                       shape = "Cylinder",
                                       radius = 5,
                                       start_loc = c(30, 30, 30),
                                       end_loc   = c(30, 30 ,70)
                                     ),
                                     C3 = list(
                                       name_of_cluster_cell = "Endothelial",
                                       infiltration_types = c("Immune", "Others"),
                                       infiltration_proportions = c(0.1, 0.05),
                                       shape = "Cylinder",
                                       radius = 5,
                                       start_loc = c(30, 30, 30),
                                       end_loc   = c(60, 40 ,30)
                                     ),
                                     C4 = list(
                                       name_of_cluster_cell = "Endothelial",
                                       infiltration_types = c("Immune", "Others"),
                                       infiltration_proportions = c(0.1, 0.05),
                                       shape = "Cylinder",
                                       radius = 4,
                                       start_loc = c(60, 40, 30),
                                       end_loc   = c(100, 100 ,100)
                                     )
                                   ),
                                   plot_image = TRUE,
                                   plot_categories = c("Others", "Immune", "Endothelial"),
                                   plot_colours = c("lightgray", "skyblue", "tomato"))

###-------------------------------------------------------------------------###
### Ellipsoid
###-------------------------------------------------------------------------###
bg_ellipsoid <- simulate_clusters3D(bg_sample = bg,
                                   n_clusters = 2,
                                   cluster_properties = list(
                                     C1 = list(
                                       name_of_cluster_cell = "Tumour",
                                       infiltration_types = c("Immune", "Others"),
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
                                       infiltration_types = c("Immune", "Others"),
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
                                   plot_categories = c("Others", "Immune", "Tumour"),
                                   plot_colours = c("lightgray", "skyblue", "red"))




###-------------------------------------------------------------------------###
### Network
###-------------------------------------------------------------------------###
bg_network <- simulate_clusters3D(bg,
                                  n_clusters = 1,
                                  cluster_properties = list(
                                    C1 = list(
                                      shape = "Network",
                                      n_edges = 15,
                                      name_of_cluster_cell = "Immune1",
                                      infiltration_types = c("Immune2", "Immune3"),
                                      infiltration_proportions = c(0.10, 0.05),
                                      width = 8,
                                      centre_loc = c(50, 50, 50), # Rough centre of network cluster
                                      radius = 50 # Rough radius spanned by the network cluster
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_categories = c("Others", "Immune1", "Immune2", "Immune3"),
                                  plot_colours = c("lightgray", "skyblue", "green", "tomato"))




###-------------------------------------------------------------------------###
### Ring
###-------------------------------------------------------------------------###
bg_ring <- simulate_rings3D(bg,
                            n_ring = 3,
                            ring_properties = list(
                              R1 = list(
                                name_of_cluster_cell = "Tumour",
                                infiltration_types = c("Others"),
                                infiltration_proportions = c(0.05),
                                shape = "Sphere",
                                radius = 20,
                                centre_loc = c(40, 40, 40),
                                name_of_ring_cell = "Immune",
                                ring_width = 5,
                                ring_infiltration_types = c("Others"),
                                ring_infiltration_proportions = c(0.15)
                              ),
                              R2 = list(
                                name_of_cluster_cell = "Void",
                                infiltration_types = NULL,
                                infiltration_proportions = NULL,
                                shape = "Cylinder",
                                radius = 8,
                                start_loc = c(0, 0, 0),
                                end_loc   = c(20, 20 , 100),
                                name_of_ring_cell = "Endothelial",
                                ring_width = 5,
                                ring_infiltration_types = c("Others"),
                                ring_infiltration_proportions = c(0.15)
                              ),
                              R3 = list(
                                name_of_cluster_cell = "Tumour",
                                infiltration_types = c("Others"),
                                infiltration_proportions = c(0.05),
                                shape = "Ellipsoid",
                                x_radius = 10,
                                y_radius = 15,
                                z_radius = 20,
                                centre_loc = c(70, 70, 70),
                                x_y_rotation = 0,
                                x_z_rotation = 0,
                                y_z_rotation = 0,
                                name_of_ring_cell = "Immune",
                                ring_width = 5,
                                ring_infiltration_types = c("Others"),
                                ring_infiltration_proportions = c(0.15)
                              )
                            ),
                            plot_image = TRUE,
                            plot_categories = c("Others", "Tumour", "Immune", "Endothelial"),
                            plot_colours = c("lightgray", "orange", "skyblue", "#FF7F7F"))



###-------------------------------------------------------------------------###
### Heart with ring
###-------------------------------------------------------------------------###
bg_heart_ring <- simulate_rings3D(bg_sample = bg,
                                     n_ring = 2,
                                     ring_properties = list(
                                       R1 = list(
                                          name_of_cluster_cell = "Tumour",
                                          infiltration_types = c("Immune", "Others"),
                                          infiltration_proportions = c(0.1, 0.05),
                                          shape = "Ellipsoid",
                                          x_radius = 15,
                                          y_radius = 15,
                                          z_radius = 25,
                                          centre_loc = c(50, 50, 50),
                                          y_z_rotation = pi/4,
                                          x_z_rotation = 0,
                                          x_y_rotation = 0,
                                          name_of_ring_cell = "Immune",
                                          ring_width = 5,
                                          ring_infiltration_types = c("Others"),
                                          ring_infiltration_proportions = c(0.15)
                                       ), 
                                       R2 = list(
                                          name_of_cluster_cell = "Tumour",
                                          infiltration_types = c("Immune", "Others"),
                                          infiltration_proportions = c(0.1, 0.05),
                                          shape = "Ellipsoid",
                                          x_radius = 15,
                                          y_radius = 15,
                                          z_radius = 25,
                                          centre_loc = c(50, 66, 50),
                                          y_z_rotation = -pi/4,
                                          x_z_rotation = 0,
                                          x_y_rotation = 0,
                                          name_of_ring_cell = "Immune",
                                          ring_width = 5,
                                          ring_infiltration_types = c("Others"),
                                          ring_infiltration_proportions = c(0.15)
                                       )
                                     ),
                                     plot_image = TRUE,
                                     plot_categories = c("Others", "Immune", "Tumour"),
                                     plot_colours = c("lightgray", "skyblue", "red"))



###-------------------------------------------------------------------------###
### Hollow cylinder
###-------------------------------------------------------------------------###
bg_cylinder_ring <- simulate_rings3D(bg_sample = bg,
                                  n_ring = 1,
                                  ring_properties = list(
                                    R1 = list(
                                      name_of_cluster_cell = "Void",
                                      infiltration_types = NULL,
                                      infiltration_proportions = NULL,
                                      shape = "Cylinder",
                                      radius = 10,
                                      start_loc = c(0, 0, 0),
                                      end_loc = c(100, 80, 60),
                                      name_of_ring_cell = "Endothelial",
                                      ring_width = 5,
                                      ring_infiltration_types = c("Others"),
                                      ring_infiltration_proportions = c(0.15)
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_categories = c("Others", "Endothelial"),
                                  plot_colours = c("lightgray", "red"))




###-------------------------------------------------------------------------###
### Network with ring
###-------------------------------------------------------------------------###
bg_network <- simulate_rings3D(bg,
                               n_ring = 1,
                               ring_properties = list(
                                 R1 = list(
                                   shape = "Network",
                                   n_edges = 15,
                                   name_of_cluster_cell = "Immune1",
                                   infiltration_types = NULL,
                                   infiltration_proportions = NULL,
                                   width = 8,
                                   centre_loc = c(50, 50, 50), # Rough centre of network cluster
                                   radius = 50, # Rough radius spanned by the network cluster
                                   name_of_ring_cell = "Immune2",
                                   ring_width = 2,
                                   ring_infiltration_types = NULL,
                                   ring_infiltration_proportions = NULL
                                 )
                               ),
                               plot_image = TRUE,
                               plot_categories = c("Others", "Immune1", "Immune2"),
                               plot_colours = c("lightgray", "skyblue", "green"))



###-------------------------------------------------------------------------###
### Double rings
###-------------------------------------------------------------------------###
bg_dr <- simulate_double_rings3D(bg,
                                 n_dr = 1,
                                 dr_properties = list(
                                   D1 = list(
                                     name_of_cluster_cell = "Tumour",
                                     infiltration_types = c("Immune1", "Others"),
                                     infiltration_proportions = c(0.1, 0.05),
                                     shape = "Sphere",
                                     radius = 35,
                                     centre_loc = c(50, 50, 50),
                                     name_of_inner_ring_cell = "Immune1",
                                     inner_ring_width = 6,
                                     inner_ring_infiltration_types = c("Others"),
                                     inner_ring_infiltration_proportions = c(0.15),
                                     name_of_outer_ring_cell = "Immune2",
                                     outer_ring_width = 3,
                                     outer_ring_infiltration_types = c("Others"),
                                     outer_ring_infiltration_proportions = c(0.15)
                                   )
                                 ),
                                 plot_image = TRUE,
                                 plot_categories = c("Others", "Tumour", "Immune1", "Immune2"),
                                 plot_colours = c("lightgray", "orange", "blue", "green"))




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
                                     plot_colours = c("lightgray", "orange", "skyblue", "blue"))

