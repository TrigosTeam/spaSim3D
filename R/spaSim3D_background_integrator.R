spaSim3D_background_integrator <- function() {

  ### Message strings
  message_background <- paste("Hello spaSim-3D user, how do you want your background cells to look like?\n
          1. Random pattern\n
          2. Ordered pattern\n\n",
                              "In a random pattern, cells are placed randomly...\n",
                              "In a ordered pattern, cells follow a regularly spaced in a hexagonal grid\n",
                              "To choose, please enter 1 or 2.\n", sep = "")
  
  message_background_random <- paste("We will need a few parameters before we can obtain the simulation\n",
                                     "    Window size - length, width and height (e.g. 100 x 100 x 100)\n",
                                     "    Number of cells (e.g. 10000 cells)\n",
                                     "    Minimum distance between cells (e.g. minimum distance of 2)\n",
                                     "If you want to change your inputs, you'll be able to at the end.\n", sep = "")
  
  message_background_ordered <- paste("We will need a few parameters before we can obtain the simulation\n",
                                      "    Window size - length, width and height (e.g. 100 x 100 x 100)\n",
                                      "    Number of cells (e.g. 10000 cells)\n",
                                      "    Amount of jitter (choose to give a bit or a lot of randomness)\n",
                                      "If you want to change your inputs, you'll be able to at the end.\n", sep = "")
  
  message_mixing <- paste("Would you like to MIX the background cells with chosen cell types randomly?\n")
    
  # Ask if user wants a 'random' or 'ordered' patterned background
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
                             "minimum distance between cells" = get_non_negative_numeric_input("minimum distance between cells"))
    display_parameters(parameter_values)
    
    # Generate random background simulation using these parameters
    message("Generating simulation...")
    simulated_spe <- simulate_random_background_cells3D(parameter_values[["number of cells"]],
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
      if (user_input_parameter_choice == 5) parameter_values[["minimum distance between cells"]] <- get_non_negative_numeric_input("minimum distance between cells")
      
      # Generate random background simulation using updated parameters
      display_parameters(parameter_values)
      message("Generating simulation...")
      simulated_spe <- simulate_random_background_cells3D(parameter_values[["number of cells"]],
                                                          parameter_values[["length"]],
                                                          parameter_values[["width"]],
                                                          parameter_values[["height"]],
                                                          parameter_values[["minimum distance between cells"]])
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Simulate ordered pattern
  else if (user_input_background == 2) {

    # Get required parameters for a ordered background from user
    message(message_background_ordered)
    parameter_values <- list("length" = get_positive_numeric_input("length"),
                             "width" = get_positive_numeric_input("width"),
                             "height" = get_positive_numeric_input("height"),
                             "number of cells" = get_positive_numeric_input("number of cells"),
                             "amount of jitter" = get_numeric_between_input("amount of jitter", 0, 1))
    display_parameters(parameter_values)
    
    # Generate ordered background simulation using these parameters
    message("Generating simulation...")
    simulated_spe <- simulate_ordered_background_cells3D(parameter_values[["number of cells"]],
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
      
      # Generate ordered background simulation using updated parameters
      display_parameters(parameter_values)
      message("Generating simulation...")
      simulated_spe <- simulate_ordered_background_cells3D(parameter_values[["number of cells"]],
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
  choose_cell_types_y_or_n <- get_y_or_n_input()
  if (choose_cell_types_y_or_n == "y") {
    simulated_spe <- get_cell_types_and_proportions_for_mixing(simulated_spe) 
  }
  message("All done!")
  
  return(simulated_spe)
}
