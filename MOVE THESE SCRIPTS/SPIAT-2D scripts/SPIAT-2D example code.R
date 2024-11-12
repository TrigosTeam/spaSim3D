spe1 <- mixed_spe_slices[[3]]


### 1. Basic Metrics ----------------------------------------------------------

# Calculate Cell Proportions
cell_props1 <- calculate_cell_proportions2D(spe1,
                                            cell_types_of_interest = NULL,
                                            plot_image = TRUE)
print(cell_props1)

cell_props2 <- calculate_cell_proportions2D(spe1,
                                            cell_types_of_interest = c("A", "B"),
                                            plot_image = TRUE)
print(cell_props2)


### 2. Colocalization metrics -------------------------------------------------

### Calculate Pairwise Distances between Cells
pairwise_distances <- calculate_pairwise_distances_between_cell_types2D(spe1,
                                                                        cell_types_of_interest = c("A", "B"),
                                                                        plot_image = TRUE)


### Calculate Minimum Distances between cells
minimum_distances <- calculate_minimum_distances_between_cell_types2D(spe1,
                                                                      cell_types_of_interest = c("A", "B"),
                                                                      plot_image = TRUE)





### Calculate Mixing Scores
mixing_scores <- calculate_mixing_scores2D(spe1,
                                           reference_cell_types = c("A", "B"),
                                           target_cell_types = c("A", "B"),
                                           radius = 20)
print(mixing_scores)

mixing_scores_gradient <- calculate_mixing_scores_gradient2D(spe1,
                                                             reference_cell_type = "A",
                                                             target_cell_type = "B",
                                                             radii = seq(10, 100, 10))


### Calculate cells in the neighbourhood
neighbourhood_cells <- calculate_cells_in_neighbourhood2D(spe1,
                                                          reference_cell_type = "A",
                                                          target_cell_types = c("A", "B"),
                                                          radius = 30,
                                                          plot_image = F)

neighbourhood_cells_gradient <- calculate_cells_in_neighbourhood_gradient2D(spe1,
                                                                            reference_cell_type = "A",
                                                                            target_cell_types = c("A", "B"),
                                                                            radii = seq(10, 100, 10),
                                                                            plot_image = T)


## Calculate cell proportions in the neighbourhood
neighbourhood_cell_proportions <- calculate_cells_in_neighbourhood_proportions2D(spe1,
                                                                                 reference_cell_type = "A",
                                                                                 target_cell_types = c("A", "B"),
                                                                                 radius = 20)
print(neighbourhood_cell_proportions)

neighbourhood_cell_proportions_gradient <- calculate_cells_in_neighbourhood_proportions_gradient2D(spe1,
                                                                                                   reference_cell_type = "A",
                                                                                                   target_cell_types = c("A", "B"),
                                                                                                   radii = seq(10, 100, 10))


### Calculate cross-K function
cross_K <- calculate_cross_K2D(spe1,
                               reference_cell_type = "A",
                               target_cell_type = "B",
                               radius = 20)
print(cross_K)


cross_K_gradient <- calculate_cross_K_gradient2D(spe1,
                                                 reference_cell_type = "A",
                                                 target_cell_type = "B",
                                                 radii = seq(10, 100, 10))


### Calculate entropy
entropy_background <- calculate_entropy_background2D(spe1,
                                                     cell_types_of_interest = c("A", "B"))

print(entropy_background)

entropy_result <- calculate_entropy2D(spe1,
                                      radius = 20,
                                      reference_cell_type = "A",
                                      target_cell_types = c("A", "B"))





entropy_gradient <- calculate_entropy_gradient2D(spe1,
                                                 reference_cell_type = "A",
                                                 target_cell_types = c("A", "B"),
                                                 radii = seq(10, 100, 10),
                                                 plot_image = FALSE)


### Using all_single_radius and all_gradient functions

all_single_radius_result <- calculate_all_single_radius_cc_metrics2D(spe1, "A", c("A", "B"), 20)

all_gradient_result <- calculate_all_gradient_cc_metrics2D(spe1, "A", c("A", "B"), seq(1, 50, 2))


### 3. Spatial Heterogeneity metrics ------------------------------------------

### calculate entropy grid metrics
entropy_grid_metrics <- calculate_entropy_grid_metrics2D(spe1,
                                                         n_splits = 8,
                                                         cell_types_of_interest = c("A", "B"),
                                                         plot_image = TRUE)


### calculate entropy prevalence
entropy_prevalence <- calculate_prevalence2D(entropy_grid_metrics,
                                             metric_colname = "entropy",
                                             threshold = 0.5)
print(entropy_prevalence)

entropy_prevalence_gradient <- calculate_prevalence_gradient2D(entropy_grid_metrics,
                                                               "entropy")

### calculate spatial autocorrelation
entropy_spatial_autocorrelation <- calculate_spatial_autocorrelation2D(entropy_grid_metrics,
                                                                       metric_colname = "entropy",
                                                                       weight_method = "IDW")
print(entropy_spatial_autocorrelation)


### calculate cell proportion grid metrics
cell_proportion_grid_metrics <- calculate_cell_proportion_grid_metrics2D(spe1,
                                                                         n_splits = 8,
                                                                         reference_cell_types = c("A"),
                                                                         target_cell_types = c("B"),
                                                                         plot_image = TRUE)


### calculate cell proportion prevalence
cell_proportion_prevalence <- calculate_prevalence2D(cell_proportion_grid_metrics,
                                                     metric_colname = "proportion",
                                                     threshold = 0.5)
print(cell_proportion_prevalence)

cell_proportion_prevalence_gradient <- calculate_prevalence_gradient2D(cell_proportion_grid_metrics,
                                                                       metric_colname = "proportion")

## calculate spatial autocorrelation for cell proportions
cell_proportion_spatial_autocorrelation <- calculate_spatial_autocorrelation2D(cell_proportion_grid_metrics,
                                                                               metric_colname = "proportion",
                                                                               weight_method = "rook")
print(cell_proportion_spatial_autocorrelation)



### Plot------
fig <- ggplot(data.frame(spatialCoords(spe1), "Cell.Type" = spe1$Cell.Type),  
              aes(Cell.X.Position, Cell.Y.Position, color = Cell.Type)) + 
  geom_point()
fig
