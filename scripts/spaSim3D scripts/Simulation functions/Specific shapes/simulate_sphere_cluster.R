simulate_sphere_cluster <- function(bg_spe, cluster_properties) {
  
  # Get sphere properties
  cluster_cell_types <- cluster_properties$cluster_cell_types
  cluster_cell_proportions <- cluster_properties$cluster_cell_proportions
  radius <- cluster_properties$radius
  centre_loc <- cluster_properties$centre_loc
  
  ## Check number of cell types matches the number of cell proportions
  if (length(cluster_cell_types) != length(cluster_cell_proportions)) stop("Number of cell types doesn't match number of cell proportion.")
  
  ## Check cell proportions are not negative or greater than 1
  if (sum(cluster_cell_proportions < 0 | cluster_cell_proportions > 1) != 0) stop("Cell proportions cannot be negative or greater than 1")
  
  ## Check cell proportions add up to 1
  if (sum(cluster_cell_proportions) != 1) stop("Sum of cell proportions is NOT 1")
  
  ## Change cell types in the sphere cluster
  spe_coords <- data.frame(spatialCoords(bg_spe))
  
  bg_spe[["Cell.Type"]] <- ifelse((spe_coords$Cell.X.Position - centre_loc[1])^2 +
                                    (spe_coords$Cell.Y.Position - centre_loc[2])^2 +
                                    (spe_coords$Cell.Z.Position - centre_loc[3])^2 < radius^2,
                                  sample(cluster_cell_types, size = ncol(bg_spe), replace = TRUE, prob = cluster_cell_proportions),
                                  bg_spe[["Cell.Type"]])
  
  # Update current meta data
  if (is.null(cluster_properties$cluster_type)) cluster_properties <- append(list(cluster_type = "regular"), cluster_properties)
  bg_spe@metadata[["simulation"]][[paste("cluster", length(bg_spe@metadata[["simulation"]]), sep="_")]] <- cluster_properties
  
  return(bg_spe)
}
