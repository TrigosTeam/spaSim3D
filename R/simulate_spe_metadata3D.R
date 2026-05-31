#' @title Simulate a SpatialExperiment object with spaSim3D metadata.
#'
#' @description This function generates a SpatialExperiment object with 
#'     parameters specified with inputted spaSim3D metadata. By doing so, it
#'     simulates a 3D tissue, beginning with the background cells of the tissue,
#'     then adding cell clusters to the tissue.
#' 
#' @param spe_metadata A list containing the metadata to generate the
#'     SpatialExperiment object. This is generated from the output of the 
#'     spaSim3D functions 'spe_metadata_background_template' and
#'    'spe_metadata_cluster_template'.
#' @param plot_image A logical indicating whether to plot 3D spatial data with 
#'     the added metadata. Defaults to TRUE.
#' @param plot_cell_types A string vector specifying the cell types to plot. If 
#'     NULL, all cell types in the `feature_colname` column will be considered. 
#'     Defaults to NULL.
#' @param plot_colours A string vector specifying the colours of the cell types
#'     when plotting. Must match the number of cell types specified in 
#'     `plot_cell_types`. If NULL, the viridis color pallete will be used. 
#'     Defaults to NULL.
#'
#' @return A SpatialExperiment object containing the cells of the simulated 3D
#'     tissue and including the inputted metadata.
#'
#' @examples
#' # Get background metadata
#' bg_metadata <- spe_metadata_background_template("random")
#' 
#' # Get cluster metadata (using background metadata as background)
#' cluster_metadata <- spe_metadata_cluster_template("regular", "sphere", bg_metadata)
#' cluster_metadata <- spe_metadata_cluster_template("ring", "ellipsoid", cluster_metadata)
#' 
#' # Get spe from cluster metadata
#' spe_clusters <- simulate_spe_metadata3D(cluster_metadata)
#' 
#' @export

simulate_spe_metadata3D <- function(spe_metadata, 
                                    plot_image = TRUE, 
                                    plot_cell_types = NULL,
                                    plot_colours = NULL) {
  
  # First element should contain background metadata
  bg_metadata <- spe_metadata[[1]]
  if (!(is.character(bg_metadata$background_type) && length(bg_metadata$background_type) == 1)) {
    stop("background_type parameter found in the metadata background list is not a character.")
  }
  
  if (bg_metadata$background_type == "random") {
    spe <- simulate_random_background_cells3D(bg_metadata$n_cells,
                                              bg_metadata$length,
                                              bg_metadata$width,
                                              bg_metadata$height,
                                              bg_metadata$minimum_distance_between_cells,
                                              plot_image = F)    
  }
  else if (bg_metadata$background_type == "ordered") {
    spe <- simulate_ordered_background_cells3D(bg_metadata$n_cells,
                                               bg_metadata$length,
                                               bg_metadata$width,
                                               bg_metadata$height,
                                               bg_metadata$jitter_proportion,
                                               plot_image = F) 
  }
  else {
    stop("background_type parameter found in the first list must be 'random' or 'ordered'.")
  }
  # Apply background mixing
  spe <- simulate_mixing3D(spe,
                           bg_metadata$cell_types,
                           bg_metadata$cell_proportions,
                           plot_image = F)
  
  ### If there is only background metadata, we are done
  if (length(spe_metadata) == 1) {
    
    # Plot
    if (plot_image) {
      fig <- plot_cells3D(spe,
                          plot_cell_types,
                          plot_colours)
      methods::show(fig)
    }
    
    return(spe)
  }
  
  ### All other elements should help to simulate clusters 
  for (i in 2:length(spe_metadata)) {
    cluster_metadata <- spe_metadata[[i]]
    
    if (!(is.character(cluster_metadata$cluster_type) && length(cluster_metadata$cluster_type) == 1)) {
      stop(paste("cluster_type parameter found in the metadata cluster list", i,"is not a character."))
    }
    
    if (cluster_metadata$cluster_type == "regular") {
      spe <- simulate_clusters3D(spe, list(cluster_metadata), plot_image = F)
    }
    else if (cluster_metadata$cluster_type == "ring") {
      spe <- simulate_rings3D(spe, list(cluster_metadata), plot_image = F)      
    }
    else if (cluster_metadata$cluster_type == "double ring") {
      spe <- simulate_double_rings3D(spe, list(cluster_metadata), plot_image = F)
    }
    else {
      stop("cluster_type parameter must be either 'regular', 'ring' or 'double ring'.")
    }
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_cells3D(spe,
                        plot_cell_types,
                        plot_colours)
    methods::show(fig)
  }
  
  return(spe)
}
