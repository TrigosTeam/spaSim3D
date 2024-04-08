###---------------------------------------------------------------------------
### 1. Background cells
###---------------------------------------------------------------------------
bg <- simulate_background_cells3D(n_cells = 10000,
                                  length = 100,
                                  width  = 100,
                                  height = 100,
                                  method = "tumour",
                                  min_d = 2,
                                  oversampling_rate = 1.2,
                                  jitter_prop = 0,
                                  cell_type = "Others",
                                  plot_image = F)

plot_cell_categories3D(bg,
                       cell_types_of_interest = "Others",
                       colour_vector = "#F0F0F0")



###---------------------------------------------------------------------------
### 2. Cell types: Random numbing sampling
###---------------------------------------------------------------------------
bg_mix <- simulate_mixing3D(bg,
                            idents = c("Tumour", "Immune", "Others"),
                            props = c(0.4, 0.4, 0.2),
                            plot_image = F)

plot_cell_categories3D(bg_mix,
                       cell_types_of_interest = c("Tumour", "Immune", "Others"),
                       colour_vector = c("orange", "skyblue", "#F0F0F0"))


###---------------------------------------------------------------------------
### 3. Cell aggregates: Geometric shapes
###---------------------------------------------------------------------------
bg_cluster <- simulate_clusters3D(bg,
                                  bg_type = "Others",
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
                                  plot_image = F)

### With 'Others' cell type
plot_cell_categories3D(bg_cluster,
                       cell_types_of_interest = c("Tumour", "Immune", "Endothelial", "Others"),
                       colour_vector = c("orange", "skyblue", "#FF7F7F", "#F0F0F0"))



###---------------------------------------------------------------------------
### 3. Immune rings: Concentric circles
###---------------------------------------------------------------------------
bg_ring <- simulate_rings3D(bg,
                            bg_type = "Others",
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
                            plot_image = F)

### With 'Others' cell type
plot_cell_categories3D(bg_ring,
                       cell_types_of_interest = c("Tumour", "Immune", "Endothelial", "Others"),
                       colour_vector = c("orange", "skyblue", "#FF7F7F", "#F0F0F0"))

