spe_metadata_cluster_template <- function(background_metadata, cluster_type, shape) {
  
  
  ### Get template for different shapes
  if (shape == "Sphere") {
    cluster_metadata <- list(shape = "Sphere",
                             cluster_cell_types = c("Tumour", "Immune", "Others"),
                             cluster_cell_proportions = c(0.8, 0.15, 0.05),
                             radius = 25,
                             centre_loc = c(40, 40, 40))
  }
  else if (shape == "Ellipsoid") {
    cluster_metadata <- list(shape = "Ellipsoid",
                             cluster_cell_types = c("Tumour", "Immune", "Others"),
                             cluster_cell_proportions = c(0.8, 0.15, 0.05),
                             x_radius = 15,
                             y_radius = 20,
                             z_radius = 25,
                             centre_loc = c(70, 70, 70),
                             x_y_rotation = 0,
                             x_z_rotation = 45,
                             y_z_rotation = 0)
  }
  else if (shape == "Cylinder") {
    cluster_metadata <- list(shape = "Cylinder",
                             cluster_cell_types = c("Endothelial", "Others"),
                             cluster_cell_proportions = c(0.95, 0.05),
                             radius = 10,
                             start_loc = c(0, 0, 0),
                             end_loc   = c(20, 20, 100)) 
  }
  else if (shape == "Network") {
    cluster_metadata <- list(shape = "Network",
                             cluster_cell_types = c("Tumour", "Others"),
                             cluster_cell_proportions = c(0.95, 0.05),
                             n_edges = 15,
                             width = 8,
                             centre_loc = c(50, 50, 50),
                             radius = 50)
  }
  else {
    stop("shape parameter must be 'Sphere', 'Ellipsoid', 'Cylinder' or 'Network'")
  }
  
  ### Add extra metadata for different cluster types
  if (cluster_type == "regular") {
    cluster_metadata <- append(list(cluster_type = "regular"), cluster_metadata)    
  }
  else if (cluster_type == "ring") {
    cluster_metadata <- append(list(cluster_type = "ring"), cluster_metadata)
    cluster_metadata$ring_cell_types <- c("Immune", "Others")
    cluster_metadata$ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$ring_width <- 5
  }
  else if (cluster_type == "double ring") {
    cluster_metadata <- append(list(cluster_type = "double ring"), cluster_metadata)
    cluster_metadata$inner_ring_cell_types <- c("Immune1", "Others")
    cluster_metadata$inner_ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$inner_ring_width <- 3
    cluster_metadata$outer_ring_cell_types <- c("Immune2", "Others")
    cluster_metadata$outer_ring_cell_proportions <- c(0.85, 0.15)
    cluster_metadata$outer_ring_width <- 2
  }
  else {
    stop("cluster_type parameter must be 'regular', 'ring' or 'double ring'")
  }
  
  background_metadata[[paste("cluster", length(background_metadata), sep="_")]] <- cluster_metadata
  
  return(background_metadata)
}
