### Example code to run -------------------------------------------------------
simulated_background_data <- spaSim3D_background_integrator()

#

### Message strings------------------------------------------------------------
message_background <- paste("Hello spaSim-3D user, how do you want your background cells to look like?\n
          1. Random pattern\n
          2. Normal pattern\n\n",
                            "In a random pattern, cells are placed randomly...\n",
                            "In a normal pattern, cells follow a regularly spaced in a hexagonal grid\n",
                            "To choose, please enter 1 or 2.\n", sep = "")

message_background_random <- paste("We will need a few parameters before we can obtain the simulation\n",
                                   "    Window size - length, width and height\n",
                                   "    Number of cells\n",
                                   "    Minimum distance between cells\n",
                                   "If you want to change your inputs, you'll be able to at the end.\n", sep = "")

message_background_normal <- paste("We will need a few parameters before we can obtain the simulation\n",
                                   "    Window size - length, width and height\n",
                                   "    Number of cells\n",
                                   "    Amount of jitter (choose to give a bit of randomness)\n",
                                   "If you want to change your inputs, you'll be able to at the end.\n", sep = "")

message_mixing <- paste("Would you like to MIX the background cells with chosen cell types randomly?\n")

message_get_cell_types <- "Keep entering the name of cell types you would like.\n    enter 'stop' to move on."

### Functions -----------------------------------------------------------------
spaSim3D_background_integrator <- function() {
  
  # Ask if user wants a 'random' or 'normal' patterned background
  message(message_background)
  user_input_background <- get_integer_input_from_options(c(1, 2))
  
  ### Simulate random pattern
  if (user_input_background == 1) {
    
    # Get required parameters for a random background from user
    message(message_background_random)
    parameter_values <- list("length" = get_positive_numeric_input("length"),
                             "width" = get_positive_numeric_input("width"),
                             "height" = get_positive_numeric_input("height"),
                             "number of cells" = get_positive_numeric_input("number of cells"),
                             "minimum distance between cells" = get_positive_numeric_input("minimum distance between cells"))
    display_parameters(parameter_values)
    
    # Generate random background simulation using these parameters
    message("Generating simulation...")
    simulated_data <- simulate_random_background_cells3D(parameter_values[["number of cells"]],
                                                         parameter_values[["length"]],
                                                         parameter_values[["width"]],
                                                         parameter_values[["height"]],
                                                         parameter_values[["minimum distance between cells"]])
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["length"]] <- get_positive_numeric_input("length")
      if (user_input_parameter_choice == 2) parameter_values[["width"]] <- get_positive_numeric_input("width")
      if (user_input_parameter_choice == 3) parameter_values[["height"]] <- get_positive_numeric_input("height")
      if (user_input_parameter_choice == 4) parameter_values[["number of cells"]] <- get_positive_numeric_input("number of cells")
      if (user_input_parameter_choice == 5) parameter_values[["minimum distance between cells"]] <- get_positive_numeric_input("minimum distance between cells")
      
      # Generate random background simulation using updated parameters
      display_parameters(parameter_values)
      message("Generating simulation...")
      simulated_data <- simulate_random_background_cells3D(parameter_values[["number of cells"]],
                                                           parameter_values[["length"]],
                                                           parameter_values[["width"]],
                                                           parameter_values[["height"]],
                                                           parameter_values[["minimum distance between cells"]])
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Simulate normal pattern
  else if (user_input_background == 2) {

    # Get required parameters for a normal background from user
    message(message_background_normal)
    parameter_values <- list("length" = get_positive_numeric_input("length"),
                             "width" = get_positive_numeric_input("width"),
                             "height" = get_positive_numeric_input("height"),
                             "number of cells" = get_positive_numeric_input("number of cells"),
                             "amount of jitter" = get_numeric_between_input("amount of jitter", 0, 1))
    display_parameters(parameter_values)
    
    # Generate normal background simulation using these parameters
    message("Generating simulation...")
    simulated_data <- simulate_normal_background_cells3D(parameter_values[["number of cells"]],
                                                         parameter_values[["length"]],
                                                         parameter_values[["width"]],
                                                         parameter_values[["height"]],
                                                         parameter_values[["amount of jitter"]])
    
    # Allow user the option to change their input parameters
    message("Would you like to change your inputs?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {

      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values))) # 5 different parameters
      
      if (user_input_parameter_choice == 1) parameter_values[["length"]] <- get_positive_numeric_input("length")
      if (user_input_parameter_choice == 2) parameter_values[["width"]] <- get_positive_numeric_input("width")
      if (user_input_parameter_choice == 3) parameter_values[["height"]] <- get_positive_numeric_input("height")
      if (user_input_parameter_choice == 4) parameter_values[["number of cells"]] <- get_positive_numeric_input("number of cells")
      if (user_input_parameter_choice == 5) parameter_values[["amount of jitter"]] <- get_numeric_between_input("amount of jitter", 0, 1)
      
      # Generate normal background simulation using updated parameters
      display_parameters(parameter_values)
      message("Generating simulation...")
      simulated_data <- simulate_normal_background_cells3D(parameter_values[["number of cells"]],
                                                           parameter_values[["length"]],
                                                           parameter_values[["width"]],
                                                           parameter_values[["height"]],
                                                           parameter_values[["amount of jitter"]])
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  
  ### Simulate mixing
  message(message_mixing)
  simulated_data <- get_cell_types_and_proportions_for_mixing(simulated_data)
  
  message("All done!")
  
  return (simulated_data)
}

get_integer_input_from_options <- function(integer_options) {
    
  first_integers <- integer_options[1:(length(integer_options) - 1)]
  last_integer <- integer_options[length(integer_options)]
  integers_string <- paste(paste(first_integers, collapse = ", "), "or", last_integer)
  
  prompt <- paste("Enter either ", integers_string, ": ", sep = "")
  invalid_input_message <- paste("Invalid input. Please enter only", integers_string)
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = prompt)
    # Try converting to integer
    int_value <- tryCatch({as.integer(user_input)}, error = function(e) NA)
    
    # Check if conversion was successful and value is in integer_options
    if (!is.na(int_value) && int_value %in% integer_options) {
      valid_input <- TRUE
      message("Valid input received!")
    } else {
      message(invalid_input_message)
    }
  }

  return (int_value)
}

get_positive_numeric_input <- function(parameter) {
  
  prompt <- paste("Enter a positive numeric value for the ", parameter, ": ", sep = "")
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = prompt)
    # Try converting to numeric
    positive_numeric_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
    
    # Check if conversion was successful and the input is positive
    if (!is.na(positive_numeric_value) && positive_numeric_value > 0) {
      valid_input <- TRUE
      message("Valid input received!")
    } 
    # Non-positive input
    else if (!is.na(positive_numeric_value) && positive_numeric_value <= 0) {
      message("Non-positive input. Please enter a positive number")
    }
    # Non-numeric input
    else {
      message("Invalid input. Please enter a numeric value.")
    }
  }
  
  return (positive_numeric_value)
}

get_numeric_between_input <- function(parameter, lower, upper) {
  
  prompt <- paste("Enter a numeric value between ", lower, " and ", upper, " for the ", parameter, ": ", sep = "")
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = prompt)
    # Try converting to numeric
    numeric_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)
    
    # Check if conversion was successful and the input is positive
    if (!is.na(numeric_value) && numeric_value >= lower && numeric_value <= upper) {
      valid_input <- TRUE
      message("Valid input received!")
    } 
    # Negative input
    else if (!is.na(numeric_value) && (numeric_value < lower || numeric_value > upper)) {
      message("Out of bounds input. Please a number between ", lower, " and ", upper, ".", sep = "")
    }
    # Non-numeric input
    else {
      message("Invalid input. Please enter a numeric value.")
    }
  }
  
  return (numeric_value)
}

display_parameters <- function(parameter_values) {
  
  message("Your current inputs are:\n")
  
  display_message <- ""
  
  for (i in seq(length(parameter_values))) {
    display_message <- paste(display_message, "    ", i, ". ", names(parameter_values)[i], ": ", parameter_values[[i]], '\n', sep = "")
  }
  message(display_message)
}

get_y_or_n_input <- function() {
  
  valid_input <- FALSE
  while (!valid_input) {
    user_input <- readline(prompt = "Enter either y or n: ")
    
    if (user_input %in% c("y", "n")) {
      valid_input <- TRUE
      message("Valid input received!")
    }
    else {
      message("Invalid input. Please enter either y or n.")
    }
  }
  
  return (user_input)
}

get_cell_types_and_proportions_for_mixing <- function(simulated_data) {
  
  choose_cell_types_y_or_n <- get_y_or_n_input()
  if (choose_cell_types_y_or_n == "y") {
    
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
        ## Generate simulation
        message("Generating simulation...")
        simulated_data <- simulate_mixing3D(simulated_data,
                                            cell_types,
                                            cell_proportions,
                                            plot_image = F)
        
        fig <- plot_cell_categories3D(simulated_data)
        print(fig)
        
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
    
    return (simulated_data)
  }
}


