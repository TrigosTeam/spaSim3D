simulate_sphere_dr <- function(spe, dr_properties) {
  
  # Check input parameters
  input_parameters <- dr_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get sphere double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  radius <- dr_properties$radius
  centre_loc <- dr_properties$centre_loc
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  ## Change cell types in the sphere ringed cluster
  spe_coords <- data.frame(spatialCoords(spe))
  
  # Start with cells in outer ring  
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                    (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                    (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= (radius + inner_ring_width + outer_ring_width)^2,
                                  sample(outer_ring_cell_types, size = ncol(spe), replace = TRUE, prob = outer_ring_cell_proportions),
                                  spe[["Cell.Type"]])
  
  # Then do cells in inner ring  
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                    (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                    (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= (radius + inner_ring_width)^2,
                                  sample(inner_ring_cell_types, size = ncol(spe), replace = TRUE, prob = inner_ring_cell_proportions),
                                  spe[["Cell.Type"]])
  
  # Then do cells in the cluster 
  spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                    (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                    (spe_coords$Cell.Z.Position - centre_loc[3])^2 <= radius^2,
                                  sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                                  spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  return(spe)
}
