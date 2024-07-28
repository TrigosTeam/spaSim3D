simulate_sphere_dr <- function(bg_spe, dr_properties) {
  
  # Get sphere double ring properties
  cluster_cell_types <- dr_properties$cluster_cell_types
  cluster_cell_proportions <- dr_properties$cluster_cell_proportions
  radius <- dr_properties$radius
  centre_loc <- dr_properties$centre_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
  
  inner_ring_cell_types <- dr_properties$inner_ring_cell_types
  inner_ring_cell_proportions <- dr_properties$inner_ring_cell_proportions
  inner_ring_width <- dr_properties$inner_ring_width
  
  ## Check number of inner ring cell types matches the number of inner ring cell proportions
  if (length(inner_ring_cell_types) != length(inner_ring_cell_proportions)) stop("Number of inner ring cell types doesn't match number of inner ring cell proportion.")
  
  ## Check inner ring cell proportions are not negative or greater than 1
  if (sum(inner_ring_cell_proportions < 0 | inner_ring_cell_proportions > 1) != 0) stop("Inner ring cell proportions cannot be negative or greater than 1")
  
  ## Check inner ring cell proportions add up to 1
  if (sum(inner_ring_cell_proportions) != 1) stop("Sum of inner ring cell proportions is NOT 1")
  
  outer_ring_cell_types <- dr_properties$outer_ring_cell_types
  outer_ring_cell_proportions <- dr_properties$outer_ring_cell_proportions
  outer_ring_width <- dr_properties$outer_ring_width
  
  ## Check number of outer ring cell types matches the number of outer ring cell proportions
  if (length(outer_ring_cell_types) != length(outer_ring_cell_proportions)) stop("Number of outer ring cell types doesn't match number of outer ring cell proportion.")
  
  ## Check outer ring cell proportions are not negative or greater than 1
  if (sum(outer_ring_cell_proportions < 0 | outer_ring_cell_proportions > 1) != 0) stop("Outer ring cell proportions cannot be negative or greater than 1")
  
  ## Check outer ring cell proportions add up to 1
  if (sum(outer_ring_cell_proportions) != 1) stop("Sum of outer ring cell proportions is NOT 1")
  
  ## Change cell types in the sphere ringed cluster
  spe_coords <- data.frame(spatialCoords(bg_spe))
  
  # Start with cells in outer ring  
  bg_spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                    (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                    (spe_coords$Cell.Z.Position - centre_loc[3])^2 < (radius + inner_ring_width + outer_ring_width)^2,
                                  sample(outer_ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = outer_ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Then do cells in inner ring  
  bg_spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                    (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                    (spe_coords$Cell.Z.Position - centre_loc[3])^2 < (radius + inner_ring_width)^2,
                                  sample(inner_ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = inner_ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Then do cells in the cluster 
  bg_spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                    (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                    (spe_coords$Cell.Z.Position - centre_loc[3])^2 < radius^2,
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(dr_properties$cluster_type)) dr_properties <- append(list(cluster_type = "double ring"), dr_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- dr_properties
  
  return(bg_spe)
}
