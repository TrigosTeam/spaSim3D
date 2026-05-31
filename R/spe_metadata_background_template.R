#' @title Create spaSim3D background metadata template.
#'
#' @description This function creates a spaSim3D background metadata template
#'     with the predefined parameters required. Users can change the parameters
#'     after creating the template. Can be used as input for the following
#'     functions: spe_metadata_cluster_template, simulate_spe_metadata3D.
#' 
#' @param background_type Either "random" or "ordered". These are the two types
#'     of backgrounds available in spaSim3D, each with their own set of
#'     parameters.
#' @param original_spe_metadata A list containing the metadata for a spaSim3D
#'     SpatialExperiment object. This is generated from the output of the 
#'     spaSim3D function 'spe_metadata_cluster_template'. Defaults to NULL. If 
#'     NULL, the background metadata template will be returned. If original 
#'     metadata is inputted, the background element of the original metadata 
#'     will be replaced with the background metadata template. If the original 
#'     metadata didn't have a background element, a background element will be 
#'     added, containing the background metadata template.
#'
#' @return A list containing the background metadata template or the original
#'     metadata updated with the background metadata template.
#'
#' @examples
#' # Get background metadata
#' bg_metadata <- spe_metadata_background_template("random")
#' 
#' # Change background metadata
#' bg_metadata$background$n_cells <- 25000
#' 
#' Get spe from background metadata
#' spe_background <- simulate_spe_metadata3D(bg_metadata)
#' 
#' @export

spe_metadata_background_template <- function(background_type, 
                                             original_spe_metadata = NULL) {
  
  if (background_type == "random") {
    background_metadata <- list(background = list(background_type = "random",
                                                  n_cells = 20000,
                                                  length = 600,
                                                  width = 600,
                                                  height = 300,
                                                  minimum_distance_between_cells = 10,
                                                  cell_types = c("Tumour", "Others"),
                                                  cell_proportions = c(0.05, 0.95)))
  }
  else if (background_type == "ordered") {
    background_metadata <- list(background = list(background_type = "ordered",
                                                  n_cells = 20000,
                                                  length = 600,
                                                  width = 300,
                                                  height = 300,
                                                  jitter_proportion = 0.25,
                                                  cell_types = c("Immune", "Others"),
                                                  cell_proportions = c(0.05, 0.95)))
  }
  else {
    stop("background_type parameter must be 'random' or 'ordered'.")
  }
  
  
  # If original_spe_metadata input is not null, replace its background metadata with new background metadata
  if (!is.null(original_spe_metadata) && !is.null(original_spe_metadata[["background"]])) {
    original_spe_metadata[["background"]] <- background_metadata    
    return(original_spe_metadata)
  }
  else if (!is.null(original_spe_metadata) && is.null(original_spe_metadata[["background"]])) {
    original_spe_metadata <- c(background_metadata, original_spe_metadata)
    return(original_spe_metadata)
  }
  
  # Else, just return the background_metadata
  return(background_metadata)
}
