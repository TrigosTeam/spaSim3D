grid_m <- grid_metrics(formatted_spatial_data, calculate_entropy, 5,
                       cell_types_of_interest=c("Tumour","CD4 T"), feature_colname = "Cell.Type")

grid_percents <- calculate_percentage_of_grids(grid_m, threshold = 0.75, above = TRUE)
                     