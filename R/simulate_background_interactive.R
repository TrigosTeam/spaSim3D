#' @title Function to simulate background of 3D tissue interactively.
#'
#' @description This is an interactive function that allows users to simulate
#'     background cells of a 3D tissue. The function will ask you for your
#'     inputs, for which the user can enter them in the console.
#'
#' @return A 3D SpatialExperiment object with the background cells.
#'
#' @examples
#' # Simulate background
#' background_spe <- simulate_background_interactive()
#'
#' @export

simulate_background_interactive <- function() {

  ### Message strings
  message_background <- paste("Hello spaSim3D user, how do you want your background cells to look like?\n
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

  # Supportive functions
  get_integer_input_from_options <- function(integer_options) {

    first_integers <- integer_options[1:(length(integer_options) - 1)]
    last_integer <- integer_options[length(integer_options)]
    integers_string <- paste(paste(first_integers, collapse = ", "), "or", last_integer)

    prompt <- paste("Enter either ", integers_string, ": ", sep = "")

    valid_input <- FALSE
    while (!valid_input) {
      user_input <- readline(prompt = prompt)
      # Try converting to numeric
      integer_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)

      # Non-numeric input
      if (is.na(integer_value)) {
        message("Invalid input. Please enter a numeric integer value.")
      }
      # Numeric but not integer
      else if (integer_value%%1 != 0) {
        message("Non-integer input. Please enter an integer value.")
      }
      # Check if conversion was successful and value is in integer_options
      else if (integer_value %in% integer_options) {
        valid_input <- TRUE
        message("Valid input received!")
      }
      else {
        message(paste("Invalid input. Please enter only", integers_string))
      }
    }

    return(integer_value)
  }

  get_non_negative_numeric_input <- function(parameter) {

    prompt <- paste("Enter a non-negative numeric value for the ", parameter, ": ", sep = "")

    valid_input <- FALSE
    while (!valid_input) {
      user_input <- readline(prompt = prompt)
      # Try converting to numeric
      non_negative_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)

      # Non-numeric input
      if (is.na(non_negative_value)) {
        message("Invalid input. Please enter a numeric value.")
      }
      # Negative input
      else if (non_negative_value < 0) {
        message("Negative input. Please enter a non-negative number")
      }
      # Should be correct input
      else {
        valid_input <- TRUE
        message("Valid input received!")
      }
    }

    return(non_negative_value)
  }

  get_numeric_between_input <- function(parameter,
                                        lower,
                                        upper) {

    prompt <- paste("Enter a numeric value between ", lower, " and ", upper, " for the ", parameter, ": ", sep = "")

    valid_input <- FALSE
    while (!valid_input) {
      user_input <- readline(prompt = prompt)
      # Try converting to numeric
      numeric_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)

      # Non-numeric input
      if (is.na(numeric_value)) {
        message("Invalid input. Please enter a numeric value.")
      }
      # Out of bounds input
      else if (numeric_value < lower || numeric_value > upper) {
        message("Out of bounds input. Please a number between ", lower, " and ", upper, ".", sep = "")
      }
      # Should be correct
      else {
        valid_input <- TRUE
        message("Valid input received!")
      }
    }

    return(numeric_value)
  }

  get_positive_numeric_input <- function(parameter) {

    prompt <- paste("Enter a positive numeric value for the ", parameter, ": ", sep = "")

    valid_input <- FALSE
    while (!valid_input) {
      user_input <- readline(prompt = prompt)
      # Try converting to numeric
      positive_numeric_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)

      # Non-numeric input
      if (is.na(positive_numeric_value)) {
        message("Invalid input. Please enter a numeric value.")
      }
      # Non-positive input
      else if (positive_numeric_value <= 0) {
        message("Non-positive input. Please enter a positive number")
      }
      # Should be correct input
      else {
        valid_input <- TRUE
        message("Valid input received!")
      }
    }

    return(positive_numeric_value)
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

    return(user_input)
  }

  display_parameters <- function(parameter_values) {

    message("Your current inputs are:\n")

    display_message <- ""

    for (i in seq(length(parameter_values))) {
      display_message <- paste(display_message, "    ", i, ". ", names(parameter_values)[i], ": ", parameter_values[[i]], '\n', sep = "")
    }
    message(display_message)
  }

  get_cell_types_and_proportions_for_mixing <- function(simulated_spe) {

    message_get_cell_types <- "Keep entering the name of cell types you would like (e.g. Tumour, Immune, etc.).\n    enter 'stop' to move on."

    # Supportive functions
    get_numeric_between_input <- function(parameter,
                                          lower,
                                          upper) {

      prompt <- paste("Enter a numeric value between ", lower, " and ", upper, " for the ", parameter, ": ", sep = "")

      valid_input <- FALSE
      while (!valid_input) {
        user_input <- readline(prompt = prompt)
        # Try converting to numeric
        numeric_value <- tryCatch({as.numeric(user_input)}, error = function(e) NA)

        # Non-numeric input
        if (is.na(numeric_value)) {
          message("Invalid input. Please enter a numeric value.")
        }
        # Out of bounds input
        else if (numeric_value < lower || numeric_value > upper) {
          message("Out of bounds input. Please a number between ", lower, " and ", upper, ".", sep = "")
        }
        # Should be correct
        else {
          valid_input <- TRUE
          message("Valid input received!")
        }
      }

      return(numeric_value)
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

      return(user_input)
    }

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
