library(SPIAT)

## Choose 3D spe object
spe_temp <- spe_cluster

## Get 2D slice from 3D spe object
z_coords <- spe_temp@int_colData$spatialCoords[ , 3]
spe_slice <- spe_temp[ , (z_coords > 35 & z_coords < 65)]
plot_cell_categories(spe_slice, "Tumour", "orange")

## Get border spe
border_spe <- identify_bordering_cells(spe_slice, "Tumour")

# margin_data <- calculate_distance_to_margin(border_spe)

## Transform border spe to a data frame
border_df <- data.frame(spatialCoords(border_spe), region = border_spe$Region)

## Plot
ggplot(border_df, aes(x = Cell.X.Position, y = Cell.Y.Position, color = region)) + 
  geom_point() + 
  theme_bw()
