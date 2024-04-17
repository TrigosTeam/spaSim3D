simulate_cylinder_tentacles <- function(bg_sample = bg,
                                        n_cylinders = 15,
                                        main_cell_type = "Immune1",
                                        infiltration_types = c("Immune2", "Immune3"),
                                        infiltration_proportions = c(0.10, 0.05),
                                        thickness = 8,
                                        cluster_centre = c(50, 50, 50), # Rough centre of cylinder tentacle cluster
                                        cluster_radius = 50, # Rough radius of cylinder tentacle cluster
                                        plot_image = T,
                                        plot_categories = c("Others", "Immune1", "Immune2", "Immune3"),
                                        plot_colours = c("lightgray", "skyblue", "green", "tomato")) {  
  
  
  ### Use graph theory language: edges and vertices
  n_edges <- n_cylinders
  
  n_vertices <- n_edges + 1 # number of vertices is always one more than the number of edges for the MST will we make
  
  ## Subset coordinate within the cluster_radius of the cluster_centre
  R <- cluster_radius^2
  
  D <- (bg_sample$Cell.X.Position - cluster_centre[1])^2 +
       (bg_sample$Cell.Y.Position - cluster_centre[2])^2 +
       (bg_sample$Cell.Z.Position - cluster_centre[3])^2
  
  chosen_cells <- bg_sample[D <= R, ]
  
  ## Subset further and pick n_vertices number of cells to represent the vertices
  chosen_rows <- sample(seq(nrow(chosen_cells)), n_vertices)
  chosen_cells <- chosen_cells[chosen_rows,
                               c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")]
  
  
  ## Get adjacency matrix from points (pairwise distance between points)
  # Assume all points have an edge between each other
  # Assume weight of each edge is equal to the distance between points
  adj_mat <- -1 * apcluster::negDistMat(chosen_cells)
  
  ## Use prim's algorithm to get edges (i.e. the cells connected by each edge)
  tree_edges <- prim(adj_mat)
  
  
  
  ### Determine thickness of cylinders so that cylinders further away are thinner
  tree_edges <- data.frame(tree_edges)
  colnames(tree_edges) <- c("Cell1", "Cell2")
  tree_edges$Depth <- NA # If tree_edge is not NA, we have already accounted for it
  
  # Get cells on the 'outskirts' of MST (i.e. leaf_vertices)
  tree_vertices <- c(tree_edges[ , 1], tree_edges[ , 2])
  
  leaf_vertices <- names(table(tree_vertices))[table(tree_vertices) == 1]
  leaf_vertices <- as.numeric(leaf_vertices)
  
  # Start with leaf_vertices
  curr_vertices <- leaf_vertices
  curr_depth <- 1
  
  while (NA %in% tree_edges$Depth) {
    
    # New vertices will be those adjacent to the current vertices
    new_vertices <- c()
    
    # Check each current vertex
    for (vertex in curr_vertices) {
      # Start with Cell1
      curr_edges <- which(tree_edges$Cell1 == vertex)
      
      tree_edges[curr_edges, "Depth"][is.na(tree_edges[curr_edges, "Depth"])] <- curr_depth
      
      new_vertices <- c(new_vertices, tree_edges[curr_edges, "Cell2"])
      
      # Then Cell2
      curr_edges <- which(tree_edges$Cell2 == vertex)
      
      tree_edges[curr_edges, "Depth"][is.na(tree_edges[curr_edges, "Depth"])] <- curr_depth
      
      new_vertices <- c(new_vertices, tree_edges[curr_edges, "Cell1"])
      
      # Only keep unique vertices
      new_vertices <- unique(new_vertices)
    }
    
    curr_depth <- curr_depth + 1
    curr_vertices <- new_vertices
  }
  
  
  
  
  ## Get cluster properties using edge data
  cluster_properties <- list()
  max_depth <- max(tree_edges[["Depth"]])
  
  for (i in seq(n_vertices - 1)) {
    start_loc <- as.numeric(chosen_cells[tree_edges[i, "Cell1"], ])
    end_loc <- as.numeric(chosen_cells[tree_edges[i, "Cell2"], ])
    curr_thickness <- (1 - 0.10 * (max_depth - tree_edges[i, "Depth"])) * thickness # 10% decrease with each depth
    
    # Very unlikely case when thickness is negative, just ignore these cylinders
    if (thickness < 0) {
      thickness <- 0
    }
    
    cluster_properties[[i]] <- list(name_of_cluster_cell = main_cell_type,
                                    infiltration_types = infiltration_types,
                                    infiltration_proportions = infiltration_proportions,
                                    shape = "Cylinder",
                                    radius = curr_thickness,
                                    start_loc = start_loc,
                                    end_loc = end_loc)
  }

  tentacles_bg <- simulate_clusters3D(bg,
                                      n_clusters = n_cylinders,
                                      cluster_properties = cluster_properties,
                                      plot_image = F)
  if (plot_image) {
    plot <- plot_cell_categories3D(tentacles_bg,
                                   plot_categories,
                                   plot_colours)  
    print(plot)
  }
  
  return (tentacles_bg)
  
}
