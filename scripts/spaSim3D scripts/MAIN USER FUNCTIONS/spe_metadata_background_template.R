spe_metadata_background_template <- function(background_type) {
  
  if (background_type == "random") {
    background_metadata <- list(background = list(background_type = "random",
                                                  n_cells = 10000,
                                                  length = 100,
                                                  width = 100,
                                                  height = 100,
                                                  minimum_distance_between_cells = 2,
                                                  cell_types = c("Tumour", "Others"),
                                                  cell_proportions = c(0.05, 0.95)))
  }
  else if (background_type == "normal") {
    background_metadata <- list(background = list(background_type = "normal",
                                                  n_cells = 10000,
                                                  length = 100,
                                                  width = 100,
                                                  height = 100,
                                                  jitter_proportion = 0.25,
                                                  cell_types = c("Immune", "Others"),
                                                  cell_proportions = c(0.05, 0.95)))
  }
  else {
    stop("background_type parameter must be 'random' or 'normal'.")
  }
  
  return(background_metadata)
}
