spe1 <- spe_clusters

# Get spe slice
spe1_z_coords <- spatialCoords(spe1)[ , "Cell.Z.Position"]
spe1_slice <- spe1[ , 40 < spe1_z_coords & spe1_z_coords < 60]

# Remove the z-coords
spe1_slice@int_colData@listData$spatialCoords <- spe1_slice@int_colData@listData$spatialCoords[ , c("Cell.X.Position", "Cell.Y.Position")]


### 1. Basic Metrics ----------------------------------------------------------

# Calculate Cell Proportions
cell_props1 <- calculate_cell_proportions2D(spe1,
                                            cell_types_of_interest = NULL,
                                            plot_image = TRUE)
print(cell_props1)

cell_props2 <- calculate_cell_proportions2D(spe1,
                                            cell_types_of_interest = c("Tumour", "Immune"),
                                            plot_image = TRUE)
print(cell_props2)


### 2. Colocalization metrics -------------------------------------------------

### Calculate Pairwise Distances between Cells
pairwise_distances <- calculate_pairwise_distances_between_cell_types2D(spe1,
                                                                        cell_types_of_interest = c("Tumour", "Immune"),
                                                                        plot_image = TRUE)


### Calculate Minimum Distances between cells
minimum_distances <- calculate_minimum_distances_between_cell_types2D(spe1,
                                                                      cell_types_of_interest = c("Tumour", "Immune"),
                                                                      plot_image = TRUE)





### Calculate Mixing Scores
mixing_scores <- calculate_mixing_scores2D(spe1,
                                           reference_cell_types = c("Tumour", "Immune"),
                                           target_cell_types = c("Tumour", "Immune"),
                                           radius = 20)
print(mixing_scores)

mixing_scores_gradient <- calculate_mixing_scores_gradient2D(spe1,
                                                             reference_cell_type = "Immune",
                                                             target_cell_type = "Tumour",
                                                             radii = 50)


### Calculate cells in the neighbourhood
neighbourhood_cells <- calculate_cells_in_neighbourhood2D(spe1,
                                                          reference_cell_type = "Tumour",
                                                          target_cell_types = c("Tumour", "Immune"),
                                                          radius = 30,
                                                          plot_image = F)

neighbourhood_cells_gradient <- calculate_cells_in_neighbourhood_gradient2D(spe1,
                                                                            reference_cell_type = "Tumour",
                                                                            target_cell_types = c("Tumour", "Immune"),
                                                                            radii = 30,
                                                                            plot_image = T)


## Calculate cell proportions in the neighbourhood
neighbourhood_cell_proportions <- calculate_cells_in_neighbourhood_proportions2D(spe1,
                                                                                 reference_cell_type = "Tumour",
                                                                                 target_cell_types = c("Tumour", "Immune"),
                                                                                 radius = 20)
print(neighbourhood_cell_proportions)

neighbourhood_cell_proportions_gradient <- calculate_cells_in_neighbourhood_proportions_gradient2D(spe1,
                                                                                                   reference_cell_type = "Tumour",
                                                                                                   target_cell_types = c("Tumour", "Immune"),
                                                                                                   radii = 50)


### Calculate cross-K function
cross_K <- calculate_cross_K2D(spe1,
                               reference_cell_type = "Tumour",
                               target_cell_type = "Immune",
                               radius = 20)
print(cross_K)


cross_K_gradient <- calculate_cross_K_gradient2D(spe1,
                                                 reference_cell_type = "Tumour",
                                                 target_cell_type = "Immune",
                                                 radii = 100)

plot_cross_K_gradient_ratio2D(cross_K_gradient_results = cross_K_gradient)

# Kcross_intersection <- calculate_Kcross_intersection3D(Kcross_results)
# 
# Kcross_AUC <- calculate_AUC_of_Kcross3D(Kcross_results)
# print(Kcross_AUC)


### Calculate entropy
entropy_background <- calculate_entropy_background2D(spe1,
                                                     cell_types_of_interest = c("Tumour", "Immune"))

print(entropy_background)

entropy_result <- calculate_entropy2D(spe1,
                                      radius = 20,
                                      reference_cell_type = "Tumour",
                                      target_cell_types = c("Tumour", "Immune"),
                                      plot_image = TRUE)





entropy_gradient <- calculate_entropy_gradient2D(spe1,
                                                 reference_cell_type = "Tumour",
                                                 target_cell_types = c("Tumour", "Immune"),
                                                 radii = 100,
                                                 plot_image = TRUE)



### 3. Spatial Heterogeneity metrics ------------------------------------------

### Determine entropy grid metrics
entropy_grid_metrics <- determine_entropy_grid_metrics2D(spe1,
                                                         n_splits = 10,
                                                         cell_types_of_interest = c("Tumour", "Immune"),
                                                         plot_image = TRUE)
plot_grid_metrics_discrete2D(entropy_grid_metrics, "entropy")



### Determine entropy prevalence
entropy_prevalence <- determine_prevalence2D(entropy_grid_metrics,
                                             metric_colname = "entropy",
                                             threshold = 0.5)
print(entropy_prevalence)

entropy_prevalence_gradient <- determine_prevalence_gradient2D(entropy_grid_metrics,
                                                               "entropy")

### Determine spatial autocorrelation
entropy_spatial_autocorrelation <- determine_spatial_autocorrelation2D(entropy_grid_metrics,
                                                                       metric_colname = "entropy",
                                                                       weight_method = "IDW")
print(entropy_spatial_autocorrelation)


### Determine cell proportion grid metrics
cell_proportion_grid_metrics <- determine_cell_proportion_grid_metrics2D(spe1,
                                                                         n_splits = 8,
                                                                         reference_cell_types = c("Tumour"),
                                                                         target_cell_types = c("Immune"),
                                                                         plot_image = TRUE)
plot_grid_metrics_discrete2D(cell_proportion_grid_metrics, "proportion")


### Determine cell proportion prevalence
cell_proportion_prevalence <- determine_prevalence2D(cell_proportion_grid_metrics,
                                                     metric_colname = "proportion",
                                                     threshold = 0.5)
print(cell_proportion_prevalence)

cell_proportion_prevalence_gradient <- determine_prevalence_gradient2D(cell_proportion_grid_metrics,
                                                                       metric_colname = "proportion")

## Determine spatial autocorrelation for cell proportions
cell_proportion_spatial_autocorrelation <- determine_spatial_autocorrelation2D(cell_proportion_grid_metrics,
                                                                               metric_colname = "proportion",
                                                                               weight_method = "binary")
print(cell_proportion_spatial_autocorrelation)



### 6. Plot data -------------------------------------------------------------
plot_cells2D(spe1,
             plot_cell_types = c("Tumour", "Immune", "Others"),
             plot_colours = c("orange", "skyblue", "lightgray"))
