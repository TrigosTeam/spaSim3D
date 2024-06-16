### Add Cell.ID
chosen_bg <- bg_sphere
chosen_bg$Cell.ID <- (paste("Cell_", seq(nrow(chosen_bg)), sep="")) ## adding Cell.ID column
data <- chosen_bg


###-------------------------------------------------------------------------###
### Background
###-------------------------------------------------------------------------###
bg_r <- simulate_random_background_cells3D(n_cells = 10000,
                                          length = 100,
                                          width = 100,
                                          height = 100,
                                          minimum_distance_between_cells = 0.5,
                                          background_cell_type = "Others",
                                          plot_image = TRUE)

bg_n <- simulate_normal_background_cells3D(n_cells = 10000,
                                           length = 100,
                                           width = 100,
                                           height = 100,
                                           jitter_proportion = 0,
                                           background_cell_type = "Others",
                                           plot_image = TRUE)

###-------------------------------------------------------------------------###
### Mixing
###-------------------------------------------------------------------------###
bg_mix <- simulate_mixing3D(bg_n,
                            cell_types = c("Others", "Immune", "Tumour"),
                            cell_proportions = c(0.5, 0.25, 0.25),
                            plot_image = TRUE,
                            plot_cell_types = c("Others", "Immune", "Tumour"),
                            plot_colours = c("lightgray", "skyblue", "orange"))

###-------------------------------------------------------------------------###
### Clusters
###-------------------------------------------------------------------------###
bg_cluster <- simulate_clusters3D(bg_r,
                                  cluster_properties = list(
                                    C1 = list(
                                      shape = "Sphere",
                                      cluster_cell_types = c("Tumour", "Immune", "Others"),
                                      cluster_cell_proportions = c(0.55, 0.4, 0.05),
                                      radius = 25,
                                      centre_loc = c(40, 40, 40)
                                    ),
                                    C2 = list(
                                      shape = "Cylinder",
                                      cluster_cell_types = c("Endothelial", "Others"),
                                      cluster_cell_proportions = c(0.95, 0.05),
                                      radius = 10,
                                      start_loc = c(0, 0, 0),
                                      end_loc   = c(20, 20 , 100)
                                    ),
                                    C3 = list(
                                      shape = "Ellipsoid",
                                      cluster_cell_types = c("Tumour", "Immune", "Others"),
                                      cluster_cell_proportions = c(0.65, 0.3, 0.05),
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
                                  plot_cell_types = c("Others", "Immune", "Endothelial", "Tumour"),
                                  plot_colours = c("lightgray", "skyblue", "#FF7F7F", "orange"))

###-------------------------------------------------------------------------###
### Two separate spheres
###-------------------------------------------------------------------------###
bg_spheres <- simulate_clusters3D(bg_r,
                                  cluster_properties = list(
                                    C1 = list(
                                      shape = "Sphere",
                                      cluster_cell_types = c("Tumour"),
                                      cluster_cell_proportions = c(1),
                                      radius = 30,
                                      centre_loc = c(30, 30, 30)
                                    ),
                                    C2 = list(
                                      shape = "Sphere",
                                      cluster_cell_types = c("Immune"),
                                      cluster_cell_proportions = c(1),
                                      radius = 25,
                                      centre_loc = c(70, 70, 70)
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_cell_types = c("Others", "Tumour", "Immune"),
                                  plot_colours = c("lightgray", "orange", "skyblue"))

###-------------------------------------------------------------------------###
### One sphere with mixing
###-------------------------------------------------------------------------###
bg_sphere <-  simulate_clusters3D(bg_r,
                                  cluster_properties = list(
                                    C1 = list(
                                      shape = "Sphere",
                                      cluster_cell_types = c("Tumour", "Immune"),
                                      cluster_cell_proportions = c(0.5, 0.5),
                                      radius = 30,
                                      centre_loc = c(50, 50, 50)
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_cell_types = c("Others", "Tumour", "Immune"),
                                  plot_colours = c("lightgray", "orange", "skyblue"))


###-------------------------------------------------------------------------###
### Cylinder
###-------------------------------------------------------------------------###
bg_cylinder <- simulate_clusters3D(bg_r,
                                   cluster_properties = list(
                                     C1 = list(
                                       shape = "Cylinder",
                                       cluster_cell_types = c("Endothelial", "Immune", "Others"),
                                       cluster_cell_proportions = c(0.85, 0.1, 0.05),
                                       radius = 8,
                                       start_loc = c(0, 0, 0),
                                       end_loc   = c(30, 30 ,30)
                                     ),
                                     C2 = list(
                                       shape = "Cylinder",
                                       cluster_cell_types = c("Endothelial", "Immune", "Others"),
                                       cluster_cell_proportions = c(0.85, 0.1, 0.05),
                                       radius = 5,
                                       start_loc = c(30, 30, 30),
                                       end_loc   = c(30, 30 ,70)
                                     ),
                                     C3 = list(
                                       shape = "Cylinder",
                                       cluster_cell_types = c("Endothelial", "Immune", "Others"),
                                       cluster_cell_proportions = c(0.85, 0.1, 0.05),
                                       radius = 5,
                                       start_loc = c(30, 30, 30),
                                       end_loc   = c(60, 40 ,30)
                                     ),
                                     C4 = list(
                                       shape = "Cylinder",
                                       cluster_cell_types = c("Endothelial", "Immune", "Others"),
                                       cluster_cell_proportions = c(0.85, 0.1, 0.05),
                                       radius = 4,
                                       start_loc = c(60, 40, 30),
                                       end_loc   = c(100, 100 ,100)
                                     )
                                   ),
                                   plot_image = TRUE,
                                   plot_cell_types = c("Others", "Immune", "Endothelial"),
                                   plot_colours = c("lightgray", "skyblue", "tomato"))


###-------------------------------------------------------------------------###
### Void Cylinder
###-------------------------------------------------------------------------###
bg_void_cylinder <- simulate_clusters3D(bg_r,
                                        cluster_properties = list(
                                          C1 = list(
                                            shape = "Cylinder",
                                            cluster_cell_types = c("Void", "Immune"),
                                            cluster_cell_proportions = c(0.95, 0.05),
                                            radius = 20,
                                            start_loc = c(0, 0, 0),
                                            end_loc   = c(100, 100 , 100)
                                          )
                                        ),
                                        plot_image = TRUE,
                                        plot_cell_types = c("Others", "Immune"),
                                        plot_colours = c("lightgray", "skyblue"))


###-------------------------------------------------------------------------###
### Ellipsoid heart
###-------------------------------------------------------------------------###
bg_ellipsoid <- simulate_clusters3D(bg_r,
                                    cluster_properties = list(
                                      C1 = list(
                                        shape = "Ellipsoid",
                                        cluster_cell_types = c("Tumour", "Immune", "Others"),
                                        cluster_cell_proportions = c(0.85, 0.1, 0.05),
                                        x_radius = 20,
                                        y_radius = 20,
                                        z_radius = 30,
                                        centre_loc = c(50, 50, 50),
                                        y_z_rotation = pi/4,
                                        x_z_rotation = 0,
                                        x_y_rotation = 0
                                      ), 
                                      C2 = list(
                                        shape = "Ellipsoid",
                                        cluster_cell_types = c("Tumour", "Immune", "Others"),
                                        cluster_cell_proportions = c(0.85, 0.1, 0.05),
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
                                    plot_cell_types = c("Others", "Immune", "Tumour"),
                                    plot_colours = c("lightgray", "skyblue", "red"))




###-------------------------------------------------------------------------###
### Network
###-------------------------------------------------------------------------###
bg_network <- simulate_clusters3D(bg_r,
                                  cluster_properties = list(
                                    C1 = list(
                                      shape = "Network",
                                      cluster_cell_types = c("Immune1", "Immune2", "Immune3"),
                                      cluster_cell_proportions = c(0.85, 0.10, 0.05),
                                      n_edges = 15,
                                      width = 8,
                                      centre_loc = c(50, 50, 50), # Rough centre of network cluster
                                      radius = 50 # Rough radius spanned by the network cluster
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_cell_types = c("Others", "Immune1", "Immune2", "Immune3"),
                                  plot_colours = c("lightgray", "skyblue", "green", "tomato"))




###-------------------------------------------------------------------------###
### Ring
###-------------------------------------------------------------------------###
bg_ring <- simulate_rings3D(bg_r,
                            ring_properties = list(
                              R1 = list(
                                shape = "Sphere",
                                cluster_cell_types = c("Tumour", "Others"),
                                cluster_cell_proportions = c(0.95, 0.05),
                                radius = 20,
                                centre_loc = c(40, 40, 40),
                                ring_cell_types = c("Immune", "Others"),
                                ring_cell_proportions = c(0.85, 0.15),
                                ring_width = 5
                              ),
                              R2 = list(
                                shape = "Cylinder",
                                cluster_cell_types = c("Void"),
                                cluster_cell_proportions = c(1),
                                radius = 8,
                                start_loc = c(0, 0, 0),
                                end_loc   = c(20, 20 , 100),
                                ring_cell_types = c("Endothelial", "Others"),
                                ring_cell_proportions = c(0.85, 0.15),
                                ring_width = 5
                              ),
                              R3 = list(
                                shape = "Ellipsoid",
                                cluster_cell_types = c("Tumour", "Others"),
                                cluster_cell_proportions = c(0.95, 0.05),
                                x_radius = 10,
                                y_radius = 15,
                                z_radius = 20,
                                centre_loc = c(70, 70, 70),
                                x_y_rotation = 0,
                                x_z_rotation = 0,
                                y_z_rotation = 0,
                                ring_cell_types = c("Immune", "Others"),
                                ring_cell_proportions = c(0.85, 0.15),
                                ring_width = 5
                              )
                            ),
                            plot_image = TRUE,
                            plot_cell_types = c("Others", "Tumour", "Immune", "Endothelial"),
                            plot_colours = c("lightgray", "orange", "skyblue", "#FF7F7F"))



###-------------------------------------------------------------------------###
### Heart with ring
###-------------------------------------------------------------------------###
bg_heart_ring <- simulate_rings3D(bg_r,
                                  ring_properties = list(
                                    R1 = list(
                                      shape = "Ellipsoid",
                                      cluster_cell_types = c("Tumour", "Immune", "Others"),
                                      cluster_cell_proportions = c(0.85, 0.1, 0.05),
                                      x_radius = 15,
                                      y_radius = 15,
                                      z_radius = 25,
                                      centre_loc = c(50, 50, 50),
                                      y_z_rotation = pi/4,
                                      x_z_rotation = 0,
                                      x_y_rotation = 0,
                                      ring_cell_types = c("Immune", "Others"),
                                      ring_cell_proportions = c(0.85, 0.15),
                                      ring_width = 5
                                    ), 
                                    R2 = list(
                                      shape = "Ellipsoid",
                                      cluster_cell_types = c("Tumour", "Immune", "Others"),
                                      cluster_cell_proportions = c(0.85, 0.1, 0.05),
                                      x_radius = 15,
                                      y_radius = 15,
                                      z_radius = 25,
                                      centre_loc = c(50, 66, 50),
                                      y_z_rotation = -pi/4,
                                      x_z_rotation = 0,
                                      x_y_rotation = 0,
                                      ring_cell_types = c("Immune", "Others"),
                                      ring_cell_proportions = c(0.85, 0.15),
                                      ring_width = 5
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_cell_types = c("Others", "Immune", "Tumour"),
                                  plot_colours = c("lightgray", "skyblue", "red"))



###-------------------------------------------------------------------------###
### Hollow cylinder
###-------------------------------------------------------------------------###
bg_cylinder_ring <- simulate_rings3D(bg_r,
                                     ring_properties = list(
                                       R1 = list(
                                         shape = "Cylinder",
                                         cluster_cell_types = c("Void"),
                                         cluster_cell_proportions = c(1),
                                         radius = 10,
                                         start_loc = c(0, 0, 0),
                                         end_loc = c(100, 80, 60),
                                         ring_cell_types = c("Endothelial", "Others"),
                                         ring_cell_proportions = c(0.85, 0.15),
                                         ring_width = 5
                                       )
                                     ),
                                     plot_image = TRUE,
                                     plot_cell_types = c("Others", "Endothelial"),
                                     plot_colours = c("lightgray", "red"))




###-------------------------------------------------------------------------###
### Network with ring
###-------------------------------------------------------------------------###
bg_network <- simulate_rings3D(bg_r,
                               ring_properties = list(
                                 R1 = list(
                                   shape = "Network",
                                   cluster_cell_types = c("Immune1"),
                                   cluster_cell_proportions = c(1),
                                   n_edges = 15,
                                   width = 8,
                                   centre_loc = c(50, 50, 50), # Rough centre of network cluster
                                   radius = 50, # Rough radius spanned by the network cluster
                                   ring_cell_types = c("Immune2"),
                                   ring_cell_proportions = c(1),
                                   ring_width = 2
                                 )
                               ),
                               plot_image = TRUE,
                               plot_cell_types = c("Others", "Immune1", "Immune2"),
                               plot_colours = c("lightgray", "skyblue", "green"))



###-------------------------------------------------------------------------###
### Double rings
###-------------------------------------------------------------------------###
bg_dr <- simulate_double_rings3D(bg_r,
                                 dr_properties = list(
                                   D1 = list(
                                     shape = "Sphere",
                                     cluster_cell_types = c("Tumour", "Others"),
                                     cluster_cell_proportions = c(0.95, 0.05),
                                     radius = 20,
                                     centre_loc = c(40, 40, 40),
                                     inner_ring_cell_types = c("Immune1", "Others"),
                                     inner_ring_cell_proportions = c(0.85, 0.15),
                                     inner_ring_width = 5,
                                     outer_ring_cell_types = c("Immune2"),
                                     outer_ring_cell_proportions = c(1),
                                     outer_ring_width = 3
                                   ),
                                   D2 = list(
                                     shape = "Cylinder",
                                     cluster_cell_types = c("Void"),
                                     cluster_cell_proportions = c(1),
                                     radius = 8,
                                     start_loc = c(0, 0, 0),
                                     end_loc   = c(20, 20 , 100),
                                     inner_ring_cell_types = c("Endothelial", "Others"),
                                     inner_ring_cell_proportions = c(0.85, 0.15),
                                     inner_ring_width = 5,
                                     outer_ring_cell_types = c("Immune2"),
                                     outer_ring_cell_proportions = c(1),
                                     outer_ring_width = 3
                                   ),
                                   D3 = list(
                                     shape = "Ellipsoid",
                                     cluster_cell_types = c("Tumour", "Others"),
                                     cluster_cell_proportions = c(0.95, 0.05),
                                     x_radius = 10,
                                     y_radius = 15,
                                     z_radius = 20,
                                     centre_loc = c(70, 70, 70),
                                     x_y_rotation = 0,
                                     x_z_rotation = 0,
                                     y_z_rotation = 0,
                                     inner_ring_cell_types = c("Immune1", "Others"),
                                     inner_ring_cell_proportions = c(0.85, 0.15),
                                     inner_ring_width = 5,
                                     outer_ring_cell_types = c("Immune2"),
                                     outer_ring_cell_proportions = c(1),
                                     outer_ring_width = 3
                                   )
                                 ),
                                 plot_image = TRUE,
                                 plot_cell_types = c("Others", "Tumour", "Immune1", "Immune2", "Endothelial"),
                                 plot_colours = c("lightgray", "orange", "skyblue", "blue", "#FF7F7F")) 




###-------------------------------------------------------------------------###
### Heart with double rings
###-------------------------------------------------------------------------###
bg_heart_dr <- simulate_double_rings3D(bg_r,
                                       dr_properties = list(
                                       D1 = list(
                                         shape = "Ellipsoid",
                                         cluster_cell_types = c("Tumour", "Immune1", "Others"),
                                         cluster_cell_proportions = c(0.85, 0.1, 0.05),
                                         x_radius = 15,
                                         y_radius = 15,
                                         z_radius = 25,
                                         centre_loc = c(50, 50, 50),
                                         y_z_rotation = pi/4,
                                         x_z_rotation = 0,
                                         x_y_rotation = 0,
                                         inner_ring_cell_types = c("Immune1", "Others"),
                                         inner_ring_cell_proportions = c(0.85, 0.15),
                                         inner_ring_width = 3,
                                         outer_ring_cell_types = c("Immune2", "Others"),
                                         outer_ring_cell_proportions = c(0.85, 0.15),
                                         outer_ring_width = 2
                                       ), 
                                       D2 = list(
                                         shape = "Ellipsoid",
                                         cluster_cell_types = c("Tumour", "Immune1", "Others"),
                                         cluster_cell_proportions = c(0.85, 0.1, 0.05),
                                         x_radius = 15,
                                         y_radius = 15,
                                         z_radius = 25,
                                         centre_loc = c(50, 66, 50),
                                         y_z_rotation = -pi/4,
                                         x_z_rotation = 0,
                                         x_y_rotation = 0,
                                         inner_ring_cell_types = c("Immune1", "Others"),
                                         inner_ring_cell_proportions = c(0.85, 0.15),
                                         inner_ring_width = 3,
                                         outer_ring_cell_types = c("Immune2", "Others"),
                                         outer_ring_cell_proportions = c(0.85, 0.15),
                                         outer_ring_width = 2
                                       )
                                     ),
                                     plot_image = TRUE,
                                     plot_cell_types = c("Others", "Tumour", "Immune1", "Immune2"),
                                     plot_colours = c("lightgray", "orange", "skyblue", "blue"))



###-------------------------------------------------------------------------###
### Cylinder with double rings
###-------------------------------------------------------------------------###
bg_cylinder_dr <- simulate_double_rings3D(bg_n,
                                          dr_properties = list(
                                            D1 = list(
                                              shape = "Cylinder",
                                              cluster_cell_types = c("Tumour"),
                                              cluster_cell_proportions = c(1),
                                              radius = 10,
                                              start_loc = c(0, 0, 0),
                                              end_loc = c(100, 100, 100),
                                              inner_ring_cell_types = c("Endothelial"),
                                              inner_ring_cell_proportions = c(1),
                                              inner_ring_width = 8,
                                              outer_ring_cell_types = c("Immune"),
                                              outer_ring_cell_proportions = c(1),
                                              outer_ring_width = 5
                                            )
                                          ),
                                          plot_image = TRUE,
                                          plot_cell_types = c("Others", "Tumour", "Immune", "Endothelial"),
                                          plot_colours = c("lightgray", "orange", "skyblue", "tomato"))



###-------------------------------------------------------------------------###
### Network with double rings
###-------------------------------------------------------------------------###
bg_network_dr <- simulate_double_rings3D(bg_r,
                                         dr_properties = list(
                                           D1 = list(
                                             shape = "Network",
                                             cluster_cell_types = c("Tumour"),
                                             cluster_cell_proportions = c(1),
                                             n_edges = 15,
                                             width = 5,
                                             centre_loc = c(50, 50, 50),
                                             radius = 50,
                                             inner_ring_cell_types = c("Endothelial"),
                                             inner_ring_cell_proportions = c(1),
                                             inner_ring_width = 2,
                                             outer_ring_cell_types = c("Immune"),
                                             outer_ring_cell_proportions = c(1),
                                             outer_ring_width = 2
                                           )
                                         ),
                                         plot_image = TRUE,
                                         plot_cell_types = c("Others", "Tumour", "Immune", "Endothelial"),
                                         plot_colours = c("lightgray", "orange", "skyblue", "tomato"))

