length <- width <- height <- 100
n_cells_inflated <- 12000
min_d <- 2

# Use poisson distribution to sample points
pois_df <- poisson_distribution3D(n_cells = n_cells_inflated, 
                                  length = length, 
                                  width = width, 
                                  height = height)


pois_df_distances <- dbscan::frNN(pois_df[ , c(1, 2, 3)], 
                                  eps = min_d,
                                  query = NULL, 
                                  sort = FALSE)

pois_df_distances_ids <- pois_df_distances$id

n_cells <- nrow(pois_df)
i <- 1

chosen_cells <- seq(length(pois_df_distances_ids))


while (i < n_cells) {

  cell <- chosen_cells[i]
  cells_to_remove <- pois_df_distances_ids[[cell]]
  n_cells_to_remove <- length(cells_to_remove)
  
  if (n_cells_to_remove != 0) {
    chosen_cells <- chosen_cells[-(cells_to_remove)]
    n_cells <- n_cells - length(pois_df_distances_ids[[cell]])
  }
  i <- i + 1 
}


pois_df <- pois_df[chosen_cells, ]
