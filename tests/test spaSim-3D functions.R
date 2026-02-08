# Testing all default clusters - should show no error
spe_metadata <- spe_metadata_background_template("random")
spe_metadata <- spe_metadata_cluster_template("regular", "sphere", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("regular", "ellipsoid", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("regular", "cylinder", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("regular", "network", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("ring", "sphere", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("ring", "ellipsoid", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("ring", "cylinder", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("ring", "network", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("double ring", "sphere", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("double ring", "ellipsoid", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("double ring", "cylinder", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("double ring", "network", spe_metadata)
spe_clusters <- simulate_spe_metadata3D(spe_metadata)



# Testing with dummy_values
spe_metadata <- spe_metadata_background_template("random")
spe_metadata <- spe_metadata_cluster_template("double ring", "sphere", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("double ring", "ellipsoid", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("double ring", "cylinder", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("double ring", "network", spe_metadata)

all_parameters <- unique(unlist(sapply(spe_metadata, names)))

dummy_values <- list(-3, # x_z_rotation, y_z_rotation, x_y_rotation
                     -4.5, # x_z_rotation, y_z_rotation, x_y_rotation
                     -0.3, # x_z_rotation, y_z_rotation, x_y_rotation
                     0, # x_z_rotation, y_z_rotation, x_y_rotation, minimum_distance_between_cells
                     0.3, # x_z_rotation, y_z_rotation, x_y_rotation, minimum_distance_between_cells, length, width, height, radius, x_radius, y_radius, z_radius, ring_width, inner_ring_width, outer_ring_width
                     3, # x_z_rotation, y_z_rotation, x_y_rotation, minimum_distance_between_cells, n_cells, n_edges, length, width, height, radius, x_radius, y_radius, z_radius, ring_width, inner_ring_width, outer_ring_width
                     4.5, # x_z_rotation, y_z_rotation, x_y_rotation, minimum_distance_between_cells, length, width, height, radius, x_radius, y_radius, z_radius, ring_width, inner_ring_width, outer_ring_width
                     c(1, 2, 3), # centre_loc, start_loc, end_loc
                     c(-2, 3, 0), # centre_loc, start_loc, end_loc
                     c(1.1, 1.2, 1.3), # centre_loc, start_loc, end_loc
                     c(0.4, 0.3, 0.3), # centre_loc, start_loc, end_loc, cluster_cell_proportions (only for sphere and ellipsoid)
                     "abba", # None
                     c("A", "B", "C"), # cluster_cell_types (only for sphere and ellipsoid)
                     c("A", "B", "C", "D"), # None
                     SpatialExperiment(), # None
                     list()) # None
names(dummy_values) <- as.character(seq(length(dummy_values)))

error_df <- data.frame(matrix(nrow = 2000, ncol = 5))
colnames(error_df) <- c("metadata", "parameter_name", "dummy_value", "error_message", "error_expected")
index <- 1


get_expected_error <- function(dummy_values_index, parameter_name, shape) {
  
  error_expected <- T
  if (j == 1 && parameter_name %in% c("x_z_rotation", "y_z_rotation", "x_y_rotation")) error_expected <- F
  if (j == 2 && parameter_name %in% c("x_z_rotation", "y_z_rotation", "x_y_rotation")) error_expected <- F
  if (j == 3 && parameter_name %in% c("x_z_rotation", "y_z_rotation", "x_y_rotation")) error_expected <- F
  if (j == 4 && parameter_name %in% c("x_z_rotation", "y_z_rotation", "x_y_rotation", "minimum_distance_between_cells")) error_expected <- F
  if (j == 5 && parameter_name %in% c("x_z_rotation", "y_z_rotation", "x_y_rotation", "minimum_distance_between_cells", "length", "width", "height", "radius", "x_radius", "y_radius", "z_radius", "ring_width", "inner_ring_width", "outer_ring_width")) error_expected <- F
  if (j == 6 && parameter_name %in% c("x_z_rotation", "y_z_rotation", "x_y_rotation", "minimum_distance_between_cells", 'n_cells', "n_edges", "length", "width", "height", "radius", "x_radius", "y_radius", "z_radius", "ring_width", "inner_ring_width", "outer_ring_width")) error_expected <- F
  if (j == 7 && parameter_name %in% c("x_z_rotation", "y_z_rotation", "x_y_rotation", "minimum_distance_between_cells", "length", "width", "height", "radius", "x_radius", "y_radius", "z_radius", "ring_width", "inner_ring_width", "outer_ring_width")) error_expected <- F
  if (j == 8 && parameter_name %in% c("centre_loc", "start_loc", "end_loc")) error_expected <- F
  if (j == 9 && parameter_name %in% c("centre_loc", "start_loc", "end_loc")) error_expected <- F
  if (j == 10 && parameter_name %in% c("centre_loc", "start_loc", "end_loc")) error_expected <- F
  
  if (j == 11 && parameter_name %in% c("centre_loc", "start_loc", "end_loc")) error_expected <- F
  if (j == 11 && parameter_name %in% c("cluster_cell_proportions") && shape %in% c("sphere", "ellipsoid")) error_expected <- F
  
  if (j == 12 && parameter_name %in% c()) error_expected <- F
  if (j == 13 && parameter_name %in% c("cluster_cell_types") && shape %in% c("sphere", "ellipsoid")) error_expected <- F
  if (j == 14 && parameter_name %in% c()) error_expected <- F
  if (j == 15 && parameter_name %in% c()) error_expected <- F
  
  return(error_expected)
}


for (i in seq(length(spe_metadata))) {
  for (parameter_name in names(spe_metadata[[i]])) {
    original_parameter_value <- spe_metadata[[i]][[parameter_name]]
    
    for (j in seq(length(dummy_values))) {
      
      error_expected <- get_expected_error(j, parameter_name, spe_metadata[[i]]$shape)
      
      spe_metadata[[i]][[parameter_name]] <- dummy_values[[j]]
      error_output <- tryCatch(
        expr = {
          simulate_spe_metadata3D(spe_metadata, plot_image = F)
        },
        error = function(e) {
          return(paste(e))
        }
      )
      if (is.character(error_output)) {
        error_df[index, ] <- c(names(spe_metadata)[i], parameter_name, names(dummy_values)[j], error_output, error_expected)
      } else {
        error_df[index, ] <- c(names(spe_metadata)[i], parameter_name, names(dummy_values)[j], "No error.", error_expected)        
      }
      index <- index + 1
    }
    
    spe_metadata[[i]][[parameter_name]] <- original_parameter_value
  }  
}
#
