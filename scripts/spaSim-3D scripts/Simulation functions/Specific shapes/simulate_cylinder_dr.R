simulate_cylinder_dr <- function(spe, dr_properties) {
  
  # Check input parameters
  input_parameters <- dr_properties
  input_parameters[["spe"]] <- spe
  input_parameter_check_value <- check_input_parameters(input_parameters)
  if (!is.logical(input_parameter_check_value)) stop(input_parameter_error_message(input_parameter_check_value))
  
  # Get cylinder double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  radius <- dr_properties$radius
  start_loc <- dr_properties$start_loc
  end_loc <- dr_properties$end_loc
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  ## Check if start and end coordinates of the cylinder are the same
  if (identical(start_loc, end_loc)) warning("Start and end coordinates of the cylinder are the same.")
  
  ## Change cell types in the cylinder cluster
  spe_coords <- spatialCoords(spe)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  # Get vector between from each cell to start_loc
  v2 <- sweep(spe_coords, 2, end_loc, "-")

  # Start with cells in outer ring
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                     rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                    (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= (radius + inner_ring_width + outer_ring_width)^2), # Cell must be close enough to the cylinder line
                                  sample(outer_ring_cell_types, size = ncol(spe), replace = TRUE, prob = outer_ring_cell_proportions),
                                  spe[["Cell.Type"]])
    
  # Start with cells in inner ring
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                     rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                    (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= (radius + inner_ring_width)^2), # Cell must be close enough to the cylinder line
                                  sample(inner_ring_cell_types, size = ncol(spe), replace = TRUE, prob = inner_ring_cell_proportions),
                                  spe[["Cell.Type"]])
  
  # Then do cells in the cluster 
  spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                     rowSums(sweep(spe_coords, 2, v1, "*")) >= d1 & rowSums(sweep(spe_coords, 2, v1, "*")) <= d2) & # Cell must be between the planes
                                    (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) <= radius^2), # Cell must be close enough to the cylinder line
                                  sample(cluster_cell_types, size = ncol(spe), replace = TRUE, prob = cluster_cell_proportions),
                                  spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  spe@metadata[["simulation"]][[paste("cluster", length(spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  return(spe)
}
