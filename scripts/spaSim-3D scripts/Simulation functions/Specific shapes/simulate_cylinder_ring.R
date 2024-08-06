simulate_cylinder_ring <- function(bg_spe, ring_properties) {
  
  # Get cylinder ring properties
  cluster_cell_types <- ring_properties$cluster_cell_types
  cluster_cell_proportions <- ring_properties$cluster_cell_proportions
  radius <- ring_properties$radius
  start_loc <- ring_properties$start_loc
  end_loc <- ring_properties$end_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (!all.equal(sum(cluster_cell_proportions), 1)) stop("Sum of cell proportions is NOT 1")
  
  ring_cell_types <- ring_properties$ring_cell_types
  ring_cell_proportions <- ring_properties$ring_cell_proportions
  ring_width <- ring_properties$ring_width
  
  ## Check number of ring cell types matches the number of cell proportions
  if (length(ring_cell_types) != length(ring_cell_proportions)) stop("Number of ring cell types doesn't match number of ring cell proportion.")
  
  ## Check ring cell proportions are not negative or greater than 1
  if (sum(ring_cell_proportions < 0 | ring_cell_proportions > 1) != 0) stop("Ring cell proportions cannot be negative or greater than 1")
  
  ## Check ring cell proportions add up to 1
  if (!all.equal(sum(ring_cell_proportions), 1)) stop("Sum of ring cell proportions is NOT 1")
  
  ## Check if start and end coordinates of the cylinder are the same
  if (identical(start_loc, end_loc)) warning("Start and end coordinates of the cylinder are the same.")
  
  ## Change cell types in the cylinder cluster
  spe_coords <- spatialCoords(bg_spe)
  
  # Get directional vector
  v1 <- end_loc - start_loc
  
  # Get 'd values of planes' at start_loc and end_loc
  d1 <- sum(v1 * start_loc)
  d2 <- sum(v1 * end_loc)
  
  # Get vector between from each cell to start_loc
  v2 <- sweep(spe_coords, 2, end_loc, "-")
  
  # Start with cells in ring
  bg_spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                     rowSums(sweep(spe_coords, 2, v1, "*")) > d1 & rowSums(sweep(spe_coords, 2, v1, "*")) < d2) & # Cell must be between the planes
                                    (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) < (radius + ring_width)^2), # Cell must be close enough to the cylinder line
                                  sample(ring_cell_types, size = ncol(bg_spe), replace = TRUE, prob = ring_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Then do cells in the cluster 
  bg_spe[["Cell.Type"]] <- ifelse((!(identical(start_loc, end_loc)) & # Start and end coordinates of the cylinder are the same
                                     rowSums(sweep(spe_coords, 2, v1, "*")) > d1 & rowSums(sweep(spe_coords, 2, v1, "*")) < d2) & # Cell must be between the planes
                                    (((v1[2]*v2[ , 3] - v1[3]*v2[ , 2])^2 + (v1[1]*v2[ , 3] - v1[3]*v2[ , 1])^2 + (v1[1]*v2[ , 2] - v1[2]*v2[ , 1])^2) / (v1[1]^2 + v1[2]^2 + v1[3]^2) < radius^2), # Cell must be close enough to the cylinder line
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  
  # Update current meta data
  if (is.null(ring_properties$cluster_type)) ring_properties <- append(list(cluster_type = "ring"), ring_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- ring_properties
  
  return(bg_spe)
}
