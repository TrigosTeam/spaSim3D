simulate_sphere_cluster <- function(spe, cluster_properties) {
  
  # Check input parameters
  input_parameters <- cluster_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get sphere properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  radius <- cluster_properties$radius
  centre_loc <- cluster_properties$centre_loc
  
  # Change cell types in the sphere cluster
  spe_coords <- data.frame(spatialCoords(spe))
  
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                    (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                    (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= radius^2,
                                  sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                                  spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- cluster_properties
  
  return(spe)
}
