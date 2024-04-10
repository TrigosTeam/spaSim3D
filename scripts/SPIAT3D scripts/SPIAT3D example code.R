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

plot_cell_percentages_bar3D(cell_props1)


### 2. Colocalization metrics -------------------------------------------------

### Calculate Pairwise Distances between Cells
pairwise_distances <- calculate_pairwise_distances_between_cell_types3D(data,
                                                                        c("Tumour", "Immune"))

plot_cell_distances_violin3D(pairwise_distances,
                             scales = "free_x")

pairwise_distances_summary <- summarise_distances_between_cell_types3D(pairwise_distances)

plot_cell_distances_summary_heatmap3D(pairwise_distances_summary,
                                      metric = "Std.Dev")


### Calculate Minimum Distances between cells
minimum_distances <- calculate_minimum_distances_between_cell_types3D(data,
                                                                      c("Tumour", "Immune"))

plot_cell_distances_violin3D(minimum_distances,
                             scales = "free")

minimum_distances_summary <- summarise_distances_between_cell_types3D(minimum_distances)

plot_cell_distances_summary_heatmap3D(minimum_distances_summary,
                                      metric = "Mean")



### Calculate Mixing Scores
mixing_scores <- calculate_mixing_scores3D(data,
                                           reference_cell_types = c("Tumour", "Immune"),
                                           target_cell_types = c("Tumour", "Immune"),
                                           radius = 20)


### Calculate Cells in the Neighborhood
neighborhood_cells <- calculate_cells_in_neighborhood3D(data,
                                                        reference_cell_types = c("Tumour"),
                                                        target_cell_types = c("Tumour", "Immune", "Others"),
                                                        radius = 20)


neighborhood_cells_summary <- summarise_cells_in_neighborhood3D(neighborhood_cells)


### Calculate cross-K function
Kcross_results <- calculate_Kcross3D(data,
                                     reference_cell_type = "Tumour",
                                     target_cell_type = "Immune",
                                     distance = 35)


Kcross_intersection <- calculate_Kcross_intersection3D(Kcross_results)

Kcross_AUC <- calculate_AUC_of_Kcross3D(Kcross_results)

### Calculate entropy
entropy_entire_image <- calculate_entropy3D(data,
                                            radius = NULL,
                                            reference_cell_type = NULL,
                                            target_cell_types = c("Tumour", "Immune", "Others"),
                                            log_base = NULL)


entropy_result <- calculate_entropy3D(data,
                                      radius = 20,
                                      reference_cell_type = "Tumour",
                                      target_cell_types = c("Tumour", "Immune", "Others"),
                                      log_base = NULL)

entropy_gradient <- calculate_entropy_gradient3D(data,
                                                 radii = 20,
                                                 reference_cell_type = "Tumour",
                                                 target_cell_types = c("Tumour", "Immune", "Others"))


entropy_gradient_aggregated <- calculate_entropy_gradient_aggregated3D(data,
                                                                       radii = 40,
                                                                       reference_cell_type = "Tumour",
                                                                       target_cell_types = c("Tumour", "Immune", "Others"))


### 3. Spatial Heterogeneity metrics ------------------------------------------

### Determine entropy grid metrics
entropy_grid_metrics <- determine_entropy_grid_metrics3D(data,
                                                         n_split = 5,
                                                         target_cell_types = c("Tumour", "Immune", "Others"))

### Determine entropy prevalence
entropy_prevalence <- determine_entropy_prevalence3D(entropy_grid_metrics,
                                                     threshold = 0.5)

### Determine spatial autocorrelation
spatial_autocorrelation <- determine_spatial_autocorrelation(data,
                                                             entropy_grid_metrics,
                                                             5)




### 4. Margin of structure metrics --------------------------------------------




### 5. Presence of cluster metrics --------------------------------------------
