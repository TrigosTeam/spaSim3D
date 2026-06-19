### Expectations:
# expect_equal()
# expect_true()
# expect_error()
# expect_warning()
# expect_snapshot()

### Creating various spe inputs: ----

# Normal spe with two different cell types

# Normal spe with four different cell types

# spe with two cells and one cell type

# spe with two cells and two different cell types
df_two_cells_two_cell_types <- data.frame("Cell.X.Position" = runif(2),
                                          "Cell.Y.Position" = runif(2),
                                          "Cell.Z.Position" = runif(2),
                                          "Cell.Type" = c('A', 'B'))

spe_two_cells_two_cell_types <- SpatialExperiment::SpatialExperiment(
  assay = matrix(data = NA, nrow = 0, ncol = nrow(df_two_cells_two_cell_types)),
  colData = df_two_cells_two_cell_types,
  spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"))


# spe with one cell
df_one_cell <- data.frame("Cell.X.Position" = runif(1),
                          "Cell.Y.Position" = runif(1),
                          "Cell.Z.Position" = runif(1),
                          "Cell.Type" = 'A')

spe_one_cell <- SpatialExperiment::SpatialExperiment(
  assay = matrix(data = NA, nrow = 0, ncol = nrow(df_one_cell)),
  colData = df_one_cell,
  spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"))



### Testing spe input: ----

# spe with two cells and two different cell types
test_that("Using spe with two cells and two different cell types for plot_cells3D works", {
  fig <- plot_cells3D(spe = spe_two_cells_two_cell_types,
                      plot_cell_types = c('A', 'B'),
                      plot_colours = c('red', 'blue'),
                      feature_colname = "Cell.Type")
  expect_error(fig, NA)
})

# spe with one cell
test_that("Using spe with one cell for plot_cells3D works", {
  fig <- plot_cells3D(spe = spe_one_cell,
                      plot_cell_types = 'A',
                      plot_colours = 'red',
                      feature_colname = "Cell.Type")
  expect_error(fig, NA)
})

# spe without correct spatialCoordsNames

### Testing plot_cell_types input: ----

### Testing plot_colours input: ----

### Testing feature_colname input: ----


