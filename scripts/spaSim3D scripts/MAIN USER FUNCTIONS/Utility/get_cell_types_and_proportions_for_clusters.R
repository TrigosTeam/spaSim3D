message_get_cell_types <- "Keep entering the name of cell types you would like (e.g. Tumour, Immune, etc.).\n    enter 'stop' to move on."



get_cell_types_and_proportions_for_clusters <- function(simulated_spe, simulate_function, properties, cell_type_option, cell_proportion_option, temp_cell_type) {
  
  ## Display the cell types currently found in simulated_spe to the user
  current_cell_types <- setdiff(unique(simulated_spe[["Cell.Type"]]), temp_cell_type)
  message("Your data currently has the following cell types:\n", paste(current_cell_types, collapse = ", "), "\n")
  
  ## Get cell types from user
  cell_types <- c()
  user_input <- ""
  message(message_get_cell_types)
  while (user_input != "stop") {
    
    user_input <- readline(prompt = "Enter a cell type, or enter 'stop': ")
    
    ## Ignore if user enters a blank string
    if (user_input == "") {
      
    }
    ## Add inputted cell type to cell_types vector
    else if (user_input != "stop") {
      cell_types <- c(cell_types, user_input)
      message(paste("Cell type added:", user_input))
    }
    ## User wants to stop but hasn't entered any cell types
    else if (user_input == "stop" && length(cell_types) == 0) {
      message("You have not entered any cell types. Try again\n")
      user_input <- ""
    }
    ## User wants to stop
    else {
      message(paste("Your cell types chosen are:", paste(cell_types, collapse = ", ")))
      
      ## Allow user to re-choose cell types
      message("Would like to re-choose these cell types?\n")
      user_input_y_or_n <- get_y_or_n_input()
      if (user_input_y_or_n == "y") {
        message("Your data currently has the following cell types:\n", paste(current_cell_types, collapse = ", "), "\n")
        cell_types <- c()
        message(message_get_cell_types)
        user_input <- ""
      }
    }
  }
  properties[[1]][[cell_type_option]] <- cell_types
  
  ## Get cell proportions from user
  cell_proportions <- c()
  max_proportion <- 1
  i <- 1
  message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
  while (i <= length(cell_types)) {
    
    ## For the last cell type, we can figure out what the cell proportion must be
    if (i == length(cell_types)) {
      cell_proportions <- c(cell_proportions, max_proportion)
      message("Cell proportion for ", cell_types[i], " must be ", max_proportion)
    }
    ## Add inputted cell proportion to cell_proportions vector
    else {
      cell_proportion <- get_numeric_between_input(paste("cell proportion of", cell_types[i], "cells"), 0, max_proportion)
      cell_proportions <- c(cell_proportions, cell_proportion)
      max_proportion <- 1 - sum(cell_proportions)
      message("Cell proportion for ", cell_types[i], " is ", cell_proportion)
    }
    i <- i + 1
    
    if (i > length(cell_types)) {
      properties[[1]][[cell_proportion_option]] <- cell_proportions
      
      ## Convert spe object to data frame
      df <- data.frame(spatialCoords(simulated_spe), "Cell.Type" = simulated_spe[["Cell.Type"]])
      
      ## Just change the cell type of the temp_cell_type, no need to actually re-simulate
      for (i in seq(nrow(df))) {
        if (df[i, "Cell.Type"] == temp_cell_type) {
          df[i, "Cell.Type"] <- get_random_cell_type(cell_types, cell_proportions)
        }
      }
      
      # Add Cell.ID column to data frame
      df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")
      
      # Update current meta data
      metadata <- simulated_spe@metadata
      metadata[[length(metadata)]][[cell_type_option]] <- cell_types
      metadata[[length(metadata)]][[cell_proportion_option]] <- cell_proportions
      
      # Convert data frame to spe object
      simulated_spe_new <- SpatialExperiment(
        assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
        colData = df,
        spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
        metadata = metadata)
      
      ## Generate simulation
      message("Generating simulation...")
      fig <- plot_cells3D(simulated_spe_new)
      print(fig)
      
      if (length(cell_types) == 1) break # If there is only one cell type, proportion is always 1
      
      ## Allow user to re-choose cell proportions  
      message("Would like to re-choose these cell proportions?\n")
      user_input_y_or_n <- get_y_or_n_input()
      if (user_input_y_or_n == "y") {
        cell_proportions <- c()
        max_proportion <- 1
        i <- 1
        message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
      }
    }
  }
  
  return(list(data = simulated_spe_new, properties = properties))
}