message_start <- paste("Hello spaSim-3D user, what would you like to simulate?\n
          1. Simulate background cells\n
          2. Simulate clusters\n",
                       "If you've just started, you'll first need background cells\n",
                       "Otherwise, you can add some clusters\n",
                       "To choose, please enter 1 or 2.\n", sep = "")

message_background <- paste("How do you want your background cells to look like?\n
          1. Random pattern\n
          2. Normal pattern\n",
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


spaSim3D_integrator <- function() {
  
  message(message_start)
  user_input_start <- get_integer_input_from_options(c(1, 2))
  
  ### Simulate background cells
  if (user_input_start == 1) {
    message(message_background)
    user_input_background <- get_integer_input_from_options(c(1, 2))
    
    ### Simulate random pattern
    if (user_input_background == 1) {
      message(message_background_random)
      
      parameter_values <- list("length" = get_positive_numeric_input("length"),
                               "width" = get_positive_numeric_input("width"),
                               "height" = get_positive_numeric_input("height"),
                               "number of cells" = get_positive_numeric_input("number of cells"),
                               "minimum distance between cells" = get_positive_numeric_input("minimum distance between cells"))
      
      display_parameters(parameter_values)
      message("Generating simulation...")
      simulated_data <- simulate_random_background_cells3D(parameter_values[["number of cells"]],
                                                           parameter_values[["length"]],
                                                           parameter_values[["width"]],
                                                           parameter_values[["height"]],
                                                           parameter_values[["minimum distance between cells"]])
      message("Would you like to change your inputs?")
      user_input_y_or_n <- get_y_or_n_input()
      
      while (user_input_y_or_n == "y") {
        user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values))) # 5 different parameters
        
        if (user_input_parameter_choice == 1) parameter_values[["length"]] <- get_positive_numeric_input("length")
        if (user_input_parameter_choice == 2) parameter_values[["width"]] <- get_positive_numeric_input("width")
        if (user_input_parameter_choice == 3) parameter_values[["height"]] <- get_positive_numeric_input("height")
        if (user_input_parameter_choice == 4) parameter_values[["number of cells"]] <- get_positive_numeric_input("number of cells")
        if (user_input_parameter_choice == 5) parameter_values[["minimum distance between cells"]] <- get_positive_numeric_input("minimum distance between cells")
        
        display_parameters(parameter_values)
        message("Generating simulation...")
        simulated_data <- simulate_random_background_cells3D(parameter_values[["number of cells"]],
                                                             parameter_values[["length"]],
                                                             parameter_values[["width"]],
                                                             parameter_values[["height"]],
                                                             parameter_values[["minimum distance between cells"]])
        message("Would you like to change your inputs?")
        user_input_y_or_n <- get_y_or_n_input()
      }
    }
    ### Simulate normal pattern
    else if (user_input_background == 2) {
      message(message_background_normal)
      
      parameter_values <- list("length" = get_positive_numeric_input("length"),
                               "width" = get_positive_numeric_input("width"),
                               "height" = get_positive_numeric_input("height"),
                               "number of cells" = get_positive_numeric_input("number of cells"),
                               "amount of jitter" = get_numeric_between_input("amount of jitter", 0, 1))
      
      display_parameters(parameter_values)
      message("Generating simulation...")
      simulated_data <- simulate_normal_background_cells3D(parameter_values[["number of cells"]],
                                                           parameter_values[["length"]],
                                                           parameter_values[["width"]],
                                                           parameter_values[["height"]],
                                                           parameter_values[["amount of jitter"]])
      message("Would you like to change your inputs?")
      user_input_y_or_n <- get_y_or_n_input()
      
      while (user_input_y_or_n == "y") {
        user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values))) # 5 different parameters
        
        if (user_input_parameter_choice == 1) parameter_values[["length"]] <- get_positive_numeric_input("length")
        if (user_input_parameter_choice == 2) parameter_values[["width"]] <- get_positive_numeric_input("width")
        if (user_input_parameter_choice == 3) parameter_values[["height"]] <- get_positive_numeric_input("height")
        if (user_input_parameter_choice == 4) parameter_values[["number of cells"]] <- get_positive_numeric_input("number of cells")
        if (user_input_parameter_choice == 5) parameter_values[["amount of jitter"]] <- get_numeric_between_input("amount of jitter", 0, 1)
        
        display_parameters(parameter_values)
        message("Generating simulation...")
        simulated_data <- simulate_normal_background_cells3D(parameter_values[["number of cells"]],
                                                             parameter_values[["length"]],
                                                             parameter_values[["width"]],
                                                             parameter_values[["height"]],
                                                             parameter_values[["amount of jitter"]])
        message("Would you like to change your inputs?")
        user_input_y_or_n <- get_y_or_n_input()
      }
    }
  }
  
  ### Simulate clusters
  else if (user_input_start == 2) {
    print("Still working on adding cluster functions")
  }
  
  return (simulated_data)
}









get_integer_input_from_options <- function(integer_options) {
    
  first_integers <- integer_options[1:(length(integer_options) - 1)]
  last_integer <- integer_options[length(integer_options)]
  integers_string <- paste(paste(first_integers, collapse = ", "), "or", last_integer)
  
  prompt <- paste("Enter either ", integers_string, ": \n", sep = "")
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



get_numeric_between_input <- function(parameter, lower = 0, upper = 1) {
  
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
  
  message("Your current inputs are: ")
  
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


