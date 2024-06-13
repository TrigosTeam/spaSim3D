## Remember to add Cell.Id column if you haven't already
data <- chosen_bg


### 1. Basic Metrics ----------------------------------------------------------

# Calculate Cell Proportions
cell_props1 <- calculate_cell_proportions3D(data = data,
                                            reference_cell_types = NULL,
                                            cell_types_to_exclude = NULL)


cell_props2 <- calculate_cell_proportions3D(data = data,
                                            reference_cell_types = c("Tumour"),
                                            cell_types_to_exclude = c( "Others"))

plot_cell_percentages_bar3D(cell_props2)


### 2. Colocalization metrics -------------------------------------------------

### Calculate Pairwise Distances between Cells
pairwise_distances <- calculate_pairwise_distances_between_cell_types3D(data,
                                                                        c("Tumour", "Immune"))

plot_cell_distances_violin3D(pairwise_distances,
                             scales = "free_x")

pairwise_distances_summary <- summarise_distances_between_cell_types3D(pairwise_distances)

plot_cell_distances_summary_heatmap3D(pairwise_distances_summary,
                                      metric = "Mean")


### Calculate Minimum Distances between cells
minimum_distances <- calculate_minimum_distances_between_cell_types3D(data,
                                                                      c("Tumour", "Immune"))

plot_cell_distances_violin3D(minimum_distances,
                             scales = "free_x")

minimum_distances_summary <- summarise_distances_between_cell_types3D(minimum_distances)

plot_cell_distances_summary_heatmap3D(minimum_distances_summary,
                                      metric = "Mean")



### Calculate Mixing Scores
mixing_scores <- calculate_mixing_scores3D(data,
                                           reference_cell_types = c("Tumour", "Immune"),
                                           target_cell_types = c("Tumour", "Immune"),
                                           radius = 20)

mixing_scores_gradient <- calculate_mixing_scores_gradient3D(data,
                                                             reference_cell_type = "Immune",
                                                             target_cell_type = "Tumour",
                                                             radii = 30)


### Calculate Cells in the Neighborhood
neighborhood_cells <- calculate_cells_in_neighborhood3D(data,
                                                        reference_cell_types = c("Tumour", "Immune"),
                                                        target_cell_types = c("Tumour", "Immune"),
                                                        radius = 20)

plot_cells_in_neighborhood_violin3D(neighborhood_cells,
                                    scales = "free_x")

neighborhood_cells_summary <- summarise_cells_in_neighborhood3D(neighborhood_cells)


### Calculate cross-K function
Kcross_results <- calculate_Kcross3D(data,
                                     reference_cell_type = "Tumour",
                                     target_cell_type = "Immune",
                                     distance = 20)


Kcross_intersection <- calculate_Kcross_intersection3D(Kcross_results)

Kcross_AUC <- calculate_AUC_of_Kcross3D(Kcross_results)
print(Kcross_AUC)


### Calculate entropy
entropy_entire_image <- calculate_entropy3D(data,
                                            radius = NULL,
                                            reference_cell_type = NULL,
                                            target_cell_types = c("Tumour", "Immune", "Others"),
                                            log_base = NULL)
print(entropy_entire_image)

entropy_result <- calculate_entropy3D(data,
                                      radius = 20,
                                      reference_cell_type = "Tumour",
                                      target_cell_types = c("Tumour", "Immune", "Others"),
                                      log_base = NULL)

entropy_gradient <- calculate_entropy_gradient3D(data,
                                                 radii = 100,
                                                 reference_cell_type = "Tumour",
                                                 target_cell_types = c("Tumour", "Immune", "Others"))


entropy_gradient_aggregated <- calculate_entropy_gradient_aggregated3D(data,
                                                                       radii = 100,
                                                                       reference_cell_type = "Tumour",
                                                                       target_cell_types = c("Tumour", "Immune"))


### 3. Spatial Heterogeneity metrics ------------------------------------------

### Determine entropy grid metrics
entropy_grid_metrics <- determine_entropy_grid_metrics3D(data,
                                                         n_split = 8,
                                                         target_cell_types = c("Tumour", "Immune", "Others"),
                                                         size = 12,
                                                         plot_image = TRUE)

### Determine entropy prevalence
entropy_prevalence <- determine_prevalence3D(entropy_grid_metrics,
                                             metric_colname = "Entropy",
                                             threshold = 0.5)
print(entropy_prevalence)

### Determine spatial autocorrelation
entropy_spatial_autocorrelation <- determine_spatial_autocorrelation(entropy_grid_metrics,
                                                                     metric_colname = "Entropy",
                                                                     weight_method = "IDW")

print(entropy_spatial_autocorrelation)


### Determine cell proportion grid metrics
cell_proportion_grid_metrics <- determine_cell_proportion_grid_metrics3D(data,
                                                                         n_split = 8,
                                                                         reference_cell_types = c("Immune"),
                                                                         target_cell_types = c("Tumour"),
                                                                         size = 10,
                                                                         plot_image = TRUE)

### Determine cell proportion prevalence
cell_proportion_prevalence <- determine_prevalence3D(cell_proportion_grid_metrics,
                                                     metric_colname = "Proportion",
                                                     threshold = 0.5)
print(cell_proportion_prevalence)


## Determine spatial autocorrelation for cell proportions
cell_proportion_spatial_autocorrelation <- determine_spatial_autocorrelation(cell_proportion_grid_metrics,
                                                                             metric_colname = "Proportion",
                                                                             weight_method = "Binary")

print(cell_proportion_spatial_autocorrelation)



### 4. Margin of structure metrics --------------------------------------------




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
plot_cell_categories3D(data,
                       cell_types_of_interest = c("Tumour", "Immune", "Others"),
                       colour_vector = c("orange", "skyblue", "lightgray"))
