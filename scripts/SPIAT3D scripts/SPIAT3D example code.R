spe1 <- spe_cluster


### 1. Basic Metrics ----------------------------------------------------------

# Calculate Cell Proportions
cell_props1 <- calculate_cell_proportions3D(spe1,
                                            cell_types_of_interest = NULL,
                                            plot_image = TRUE)
print(cell_props1)

cell_props2 <- calculate_cell_proportions3D(spe1,
                                            cell_types_of_interest = c("Tumour", "Immune"),
                                            plot_image = TRUE)
print(cell_props2)


### 2. Colocalization metrics -------------------------------------------------

### Calculate Pairwise Distances between Cells
pairwise_distances <- calculate_pairwise_distances_between_cell_types3D(spe1,
                                                                        cell_types_of_interest = c("Tumour", "Immune"),
                                                                        plot_image = TRUE)


### Calculate Minimum Distances between cells
minimum_distances <- calculate_minimum_distances_between_cell_types3D(spe1,
                                                                      cell_types_of_interest = c("Tumour", "Immune"),
                                                                      plot_image = TRUE)





### Calculate Mixing Scores
mixing_scores <- calculate_mixing_scores3D(spe1,
                                           reference_cell_types = c("Tumour", "Immune"),
                                           target_cell_types = c("Tumour", "Immune"),
                                           radius = 20)
print(mixing_scores)

mixing_scores_gradient <- calculate_mixing_scores_gradient3D(spe1,
                                                             reference_cell_type = "Tumour",
                                                             target_cell_type = "Immune1",
                                                             radii = 50)


### Calculate cells in the neighbourhood
neighbourhood_cells <- calculate_cells_in_neighbourhood3D(spe1,
                                                          reference_cell_type = "Tumour",
                                                          target_cell_types = c("Tumour", "Immune"),
                                                          radius = 30)



## Calculate cell proportions in the neighbourhood
neighbourhood_cell_proportions <- calculate_cells_in_neighbourhood_proportions3D(spe1,
                                                                                 reference_cell_type = "Tumour",
                                                                                 target_cell_types = c("Immune", "Immune1"),
                                                                                 radius = 20)
print(neighbourhood_cell_proportions)

neighbourhood_cell_proportions_gradient <- calculate_cells_in_neighbourhood_proportions_gradient3D(spe1,
                                                                                                   reference_cell_type = "Tumour",
                                                                                                   target_cell_types = c("Immune1", "Immune"),
                                                                                                   radii = 30)


### Calculate cross-K function
cross_K <- calculate_cross_K3D(spe1,
                               reference_cell_type = "Tumour",
                               target_cell_type = "Immune",
                               radius = 20)
print(cross_K)


cross_K_gradient <- calculate_cross_K_gradient3D(spe1,
                                                 reference_cell_type = "Tumour",
                                                 target_cell_type = "Immune",
                                                 radii = 50)

plot_cross_K_gradient_ratio3D(cross_K_gradient_results = cross_K_gradient)

# Kcross_intersection <- calculate_Kcross_intersection3D(Kcross_results)
# 
# Kcross_AUC <- calculate_AUC_of_Kcross3D(Kcross_results)
# print(Kcross_AUC)


### Calculate entropy
entropy_background <- calculate_entropy_background3D(spe1,
                                                     cell_types_of_interest = c("Tumour", "Immune"))

print(entropy_background)

entropy_result <- calculate_entropy3D(spe1,
                                      radius = 20,
                                      reference_cell_type = "Tumour",
                                      target_cell_types = c("Tumour", "Immune"),
                                      plot_image = TRUE)





entropy_gradient <- calculate_entropy_gradient3D(spe1,
                                                 reference_cell_type = "Tumour",
                                                 target_cell_types = c("Tumour", "Immune"),
                                                 radii = 65,
                                                 plot_image = TRUE)



### 3. Spatial Heterogeneity metrics ------------------------------------------

### Determine entropy grid metrics
entropy_grid_metrics <- determine_entropy_grid_metrics3D(spe1,
                                                         n_split = 8,
                                                         cell_types_of_interest = c("Tumour", "Immune"),
                                                         plot_image = TRUE)
plot_grid_metrics_discrete3D(entropy_grid_metrics, "entropy")


### Determine entropy prevalence
entropy_prevalence <- determine_prevalence3D(entropy_grid_metrics,
                                             metric_colname = "entropy",
                                             threshold = 0.5)
print(entropy_prevalence)

### Determine spatial autocorrelation
entropy_spatial_autocorrelation <- determine_spatial_autocorrelation3D(entropy_grid_metrics,
                                                                       metric_colname = "entropy",
                                                                       weight_method = "IDW")
print(entropy_spatial_autocorrelation)


### Determine cell proportion grid metrics
cell_proportion_grid_metrics <- determine_cell_proportion_grid_metrics3D(spe1,
                                                                         n_split = 8,
                                                                         reference_cell_types = c("Tumour"),
                                                                         target_cell_types = c("Immune", "Immune1"),
                                                                         plot_image = TRUE)
plot_grid_metrics_discrete3D(cell_proportion_grid_metrics, "proportion")


### Determine cell proportion prevalence
cell_proportion_prevalence <- determine_prevalence3D(cell_proportion_grid_metrics,
                                                     metric_colname = "proportion",
                                                     threshold = 0.5)
print(cell_proportion_prevalence)


## Determine spatial autocorrelation for cell proportions
cell_proportion_spatial_autocorrelation <- determine_spatial_autocorrelation3D(cell_proportion_grid_metrics,
                                                                               metric_colname = "proportion",
                                                                               weight_method = "binary")
print(cell_proportion_spatial_autocorrelation)



### 4. Margin of structure metrics --------------------------------------------
spe_alpha_hull <- determine_alpha_hull3D(spe1, c("Tumour", "Immune"), alpha = 3.85, minimum_cells_in_alpha_hull = 15)

plot_alpha_hull3D(spe_alpha_hull, c("Tumour", "Immune", "Immune1", "Others"), c("orange", "skyblue", "lightgreen", "lightgray"))

alpha_hull_props <- calculate_alpha_hull_cell_proportions3D(spe_alpha_hull)

alpha_hull_min_distances <- calculate_minimum_distances_to_alpha_hull3D(spe_alpha_hull, cell_types_of_interest = c("Tumour", "Immune", "Immune1"))

### 5. Presence of cluster metrics --------------------------------------------
ANNI_result <- average_nearest_neighbor_index3D(data,
                                                cell_types_of_interest = c("Tumour", "Immune"),
                                                n_simulations = 10)

library(ggplot2)
df <- data.frame(x = 0,
                 y = ANNI_result$ANNI_mean,
                 low = ANNI_result$ANNI_95confidence_interval[1],
                 up = ANNI_result$ANNI_95confidence_interval[2])

ggplot(df, aes(x, y)) + geom_point() + geom_errorbar(aes(ymin = low, ymax = up))



### 6. Plot data -------------------------------------------------------------
plot_cells3D(spe1,
             plot_cell_types = c("Tumour", "Immune", "Others"),
             plot_colours = c("orange", "skyblue", "lightgray"))
