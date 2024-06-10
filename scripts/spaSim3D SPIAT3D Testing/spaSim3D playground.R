### Background ----------------------------------------------------------------
bg <- simulate_background_cells3D(n_cells = 10000,
                                  length = 150,
                                  width  = 150,
                                  height = 150,
                                  method = "tumour",
                                  min_d = 6,
                                  oversampling_rate = 1.2,
                                  jitter_prop = 0,
                                  cell_type = "Others",
                                  plot_image = T)

### Mixing background ---------------------------------------------------------
bg_mix <- simulate_mixing3D(bg,
                            cell_types = c("Others", "Immune1"),
                            props = c(0.85, 0.15),
                            plot_image = TRUE,
                            plot_categories = c("Others", "Immune1"),
                            plot_colours = c("lightgray", "lightgreen"))


### Blood vessels -------------------------------------------------------------
bg_cylinder <- simulate_clusters3D(bg_sample = bg_mix,
                                   n_clusters = 3,
                                   cluster_properties = list(
                                     C1 = list(
                                       name_of_cluster_cell = "Blood vessel",
                                       infiltration_types = NULL,
                                       infiltration_proportions = NULL,
                                       shape = "Cylinder",
                                       radius = 15,
                                       start_loc = c(150, 0, 0),
                                       end_loc   = c(75, 100 ,50)
                                     ),
                                     C2 = list(
                                       name_of_cluster_cell = "Blood vessel",
                                       infiltration_types = NULL,
                                       infiltration_proportions = NULL,
                                       shape = "Cylinder",
                                       radius = 10,
                                       start_loc = c(75, 100, 50),
                                       end_loc   = c(100, 150 ,150)
                                     ),
                                     C3 = list(
                                       name_of_cluster_cell = "Blood vessel",
                                       infiltration_types = NULL,
                                       infiltration_proportions = NULL,
                                       shape = "Cylinder",
                                       radius = 10,
                                       start_loc = c(75, 100, 50),
                                       end_loc   = c(150, 150 ,0)
                                     )
                                   ),
                                   plot_image = TRUE,
                                   plot_categories = c("Others", "Immune1", "Blood vessel"),
                                   plot_colours = c("lightgray", "lightgreen", "tomato"))


### Tumour ring cluster ------------------------------------------------------
bg_tumour_ring <- simulate_rings3D(bg_sample = bg_cylinder,
                                   n_ring = 1,
                                   ring_properties = list(
                                     R1 = list(
                                       name_of_cluster_cell = "Tumour",
                                       infiltration_types = NULL,
                                       infiltration_proportions = NULL,
                                       shape = "Ellipsoid",
                                       x_radius = 30,
                                       y_radius = 30,
                                       z_radius = 40,
                                       centre_loc = c(75, 75, 100),
                                       y_z_rotation = -pi/4,
                                       x_z_rotation = 0,
                                       x_y_rotation = 0,
                                       name_of_ring_cell = "Immune2",
                                       ring_width = 10,
                                       ring_infiltration_types = NULL,
                                       ring_infiltration_proportions = NULL
                                     )
                                   ),
                                   plot_image = TRUE,
                                   plot_categories = c("Others", "Blood vessel", "Immune1", "Immune2", "Tumour"),
                                   plot_colours = c("lightgray", "red", "lightgreen", "skyblue", "orange"))




### Mini tumour clusters ------------------------------------------------------
bg_tumour_clusters <- simulate_clusters3D(bg_sample = bg_tumour_ring,
                                          n_clusters = 3,
                                          cluster_properties = list(
                                            C1 = list(
                                              name_of_cluster_cell = "Immune2",
                                              infiltration_types = NULL,
                                              infiltration_proportions = NULL,
                                              shape = "Sphere",
                                              radius = 20,
                                              centre_loc = c(75, 20, 50)
                                            ),
                                            C2 = list(
                                              name_of_cluster_cell = "Immune2",
                                              infiltration_types = NULL,
                                              infiltration_proportions = NULL,
                                              shape = "Sphere",
                                              radius = 20,
                                              centre_loc = c(75, 140, 90)
                                            ),
                                            C3 = list(
                                              name_of_cluster_cell = "Immune2",
                                              infiltration_types = NULL,
                                              infiltration_proportions = NULL,
                                              shape = "Sphere",
                                              radius = 20,
                                              centre_loc = c(140, 120, 130)
                                            )
                                          ),
                                          plot_image = TRUE,
                                          plot_categories = c("Others", "Blood vessel", "Immune1", "Immune2", "Tumour"),
                                          plot_colours = c("lightgray", "red", "lightgreen", "blue", "orange"))


### Plot ----------------------------------------------------------------------
plot_cell_categories3D(bg_tumour_clusters,
                       c("Others", "Blood vessel", "Immune1", "Immune2", "Tumour"),
                       c("lightgray", "red", "green", "blue", "orange"))




### Ellipsoid cluster with ring ---------------------------------------
bg <- simulate_background_cells3D(n_cells = 10000,
                                  length = 150,
                                  width  = 150,
                                  height = 150,
                                  method = "tumour",
                                  min_d = 6,
                                  oversampling_rate = 1.2,
                                  jitter_prop = 0,
                                  cell_type = "Others",
                                  plot_image = T)

bg_ellipsoid <- simulate_rings3D(bg_sample = bg,
                                 n_ring = 1,
                                 ring_properties = list(
                                   R1 = list(
                                     name_of_cluster_cell = "Tumour",
                                     infiltration_types = NULL,
                                     infiltration_proportions = NULL,
                                     shape = "Ellipsoid",
                                     x_radius = 60,
                                     y_radius = 60,
                                     z_radius = 75,
                                     centre_loc = c(75, 75, 75),
                                     y_z_rotation = -pi/4,
                                     x_z_rotation = 0,
                                     x_y_rotation = 0,
                                     name_of_ring_cell = "Immune",
                                     ring_width = 10,
                                     ring_infiltration_types = NULL,
                                     ring_infiltration_proportions = NULL
                                   )
                                 ),
                                 plot_image = TRUE,
                                 plot_categories = c("Others", "Tumour", "Immune"),
                                 plot_colours = c("lightgray", "orange", "skyblue"))


plot_cell_categories3D(bg_ellipsoid,
                       c("Others", "Tumour", "Immune"),
                       c("lightgray", "orange", "skyblue"))

### Mega network --------------------------
bg <- simulate_background_cells3D(n_cells = 10000,
                                  length = 150,
                                  width  = 150,
                                  height = 150,
                                  method = "tumour",
                                  min_d = 6,
                                  oversampling_rate = 1.2,
                                  jitter_prop = 0,
                                  cell_type = "Others",
                                  plot_image = T)

bg_network <- simulate_clusters3D(bg,
                                  n_clusters = 1,
                                  cluster_properties = list(
                                    C1 = list(
                                      shape = "Network",
                                      n_edges = 30,
                                      name_of_cluster_cell = "Tumour",
                                      infiltration_types = c("Immune"),
                                      infiltration_proportions = c(0.2),
                                      width = 17,
                                      centre_loc = c(70, 70, 70), # Rough centre of network cluster
                                      radius = 75 # Rough radius spanned by the network cluster
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_categories = c("Others", "Tumour", "Immune"),
                                  plot_colours = c("lightgray", "orange", "skyblue"))

plot_cell_categories3D(bg_network,
                       c("Others", "Tumour", "Immune"),
                       c("lightgray", "orange", "skyblue"))



### TLSN surrounded by small infiltrated tumour clusters -------------------
bg <- simulate_background_cells3D(n_cells = 10000,
                                  length = 150,
                                  width  = 150,
                                  height = 150,
                                  method = "tumour",
                                  min_d = 6,
                                  oversampling_rate = 1.2,
                                  jitter_prop = 0,
                                  cell_type = "Others",
                                  plot_image = T)

bg_mix <- simulate_mixing3D(bg,
                            cell_types = c("Others", "Immune1"),
                            props = c(0.85, 0.15),
                            plot_image = TRUE,
                            plot_categories = c("Others", "Immune1"),
                            plot_colours = c("lightgray", "green"))


bg_network <- simulate_clusters3D(bg_mix,
                                  n_clusters = 1,
                                  cluster_properties = list(
                                    C1 = list(
                                      shape = "Network",
                                      n_edges = 30,
                                      name_of_cluster_cell = "Immune2",
                                      infiltration_types = c("Immune3"),
                                      infiltration_proportions = c(0.5),
                                      width = 14,
                                      centre_loc = c(70, 70, 70), # Rough centre of network cluster
                                      radius = 75 # Rough radius spanned by the network cluster
                                    )
                                  ),
                                  plot_image = TRUE,
                                  plot_categories = c("Others", "Immune1", "Immune2", "Immune3"),
                                  plot_colours = c("lightgray", "green", "blue", "purple"))



bg_tumour_clusters <- simulate_clusters3D(bg_sample = bg_network,
                                          n_clusters = 3,
                                          cluster_properties = list(
                                            C1 = list(
                                              name_of_cluster_cell = "Tumour",
                                              infiltration_types = c("Immune3"),
                                              infiltration_proportions = c(0.4),
                                              shape = "Sphere",
                                              radius = 25,
                                              centre_loc = c(20, 130, 30)
                                            ),
                                            C2 = list(
                                              name_of_cluster_cell = "Tumour",
                                              infiltration_types = c("Immune3"),
                                              infiltration_proportions = c(0.4),
                                              shape = "Sphere",
                                              radius = 30,
                                              centre_loc = c(130, 90, 20)
                                            ),
                                            C3 = list(
                                              name_of_cluster_cell = "Tumour",
                                              infiltration_types = c("Immune3"),
                                              infiltration_proportions = c(0.4),
                                              shape = "Sphere",
                                              radius = 32,
                                              centre_loc = c(75, 75, 130)
                                            )
                                          ),
                                          plot_image = TRUE,
                                          plot_categories = c("Others", "Immune1", "Immune2", "Immune3", "Tumour"),
                                          plot_colours = c("lightgray", "green", "blue", "purple", "orange"))


plot_cell_categories3D(bg_tumour_clusters,
                       c("Others", "Immune1", "Immune2", "Immune3", "Tumour"),
                       c("lightgray", "green", "blue", "purple", "orange"))




### Sphere ??? -----
df <- data.frame(x=c(1,3,-3,-2), y=c(2,5,2,1),z=c(1,7,4,1))

library(rgl)
open3d()
plot3d(df,col=3,type="p", radius=0.5)
plot3d(df,col=rgb(1,0,0.3),alpha=0.5, add=T,type="s",radius=1)


sphere <- spheres3d(df$x, df$y, df$z, radius = 0.5, alpha = 0.5, color = "red")



icos <- icosahedron3d()
ids <- mfrow3d(1,2)
shade3d(icos, col = "red")
next3d()
shade3d(icos, col = "green")
