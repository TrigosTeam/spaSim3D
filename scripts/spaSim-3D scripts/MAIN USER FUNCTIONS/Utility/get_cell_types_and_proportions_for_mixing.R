message_get_cell_types <- "Keep entering the name of cell types you would like (e.g. Tumour, Immune, etc.).\n    enter 'stop' to move on."



get_cell_types_and_proportions_for_mixing <- function(simulated_spe) {
  
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
        cell_types <- c()
        message(message_get_cell_types)
        user_input <- ""
      }
    }
  }
  
  ## Get cell proportions from user
  cell_proportions <- c()
  max_proportion <- 1
  i <- 1
  message("For each cell type, choose their proportion in the simulation. They must add to 1.\n")
  while (i <= length(cell_types)) {
    
    ## For the last cell type, we can figure out what the cell proportion must be
    if (i == length(cell_types)) {
      cell_proportions <- c(cell_proportions, max_proportion)
      message("Cell proportion for ", cell_types[i], " must be ", round(max_proportion, 5))
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
      ## Generate simulation
      message("Generating simulation...")
      simulated_spe <- simulate_mixing3D(simulated_spe,
                                         cell_types,
                                         cell_proportions,
                                         plot_image = F)
      
      fig <- plot_cells3D(simulated_spe)
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
  
  return(simulated_spe)
}
