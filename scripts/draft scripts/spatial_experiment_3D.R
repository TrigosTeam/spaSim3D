library(SpatialExperiment)

# Data from spaSim3D_background_integrator: simulated_background_data

simulated_background_spe <- SpatialExperiment(
  assay = matrix(data = NA, nrow = nrow(simulated_background_data), ncol = nrow(simulated_background_data)),
  colData = simulated_background_data,
  spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"))

if (user_input_background == 1) {
  parameter_values <- append(list(background_type = "random"), parameter_values)
} else if (user_input_background == 2) {
  parameter_values <- append(list(background_type = "normal"), parameter_values)
}

simulated_background_spe@metadata <- list(background = list(parameter_values))



x <- list(b = 2, c = 3)
new_element_name <- "a"
new_element_value <- 1

x <- append(list(a = 1), x)
