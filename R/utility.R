check_input_parameters <- function(input_parameters) {
  
  input_parameter_names <- names(input_parameters)
  
  check_value <- 0
  
  is_equal_with_tolerance <- function(x, y, tolerance = 1e-6) {
    abs(x - y) <= tolerance
  }
  
  for (input_parameter_name in input_parameter_names) {
    input_parameter <- input_parameters[[input_parameter_name]]
    
    # spe
    if (input_parameter_name == "spe" && class(input_parameter) != "SpatialExperiment") {
      check_value <- 1
      break
    }
    # Positive integer
    if (input_parameter_name %in% c("n_cells", "n_edges") && !(is.integer(input_parameter) && length(input_parameter) == 1 || (is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter > 0 && input_parameter%%1 == 0))) {
      check_value <- 2
      break
    }  
    # Positive numeric
    if (input_parameter_name %in% c("length", "width", "height", "radius", "x_radius", "y_radius", "z_radius", "ring_width", "inner_ring_width", "outer_ring_width") && !(is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter > 0)) {
      check_value <- 3
      break
    }
    # Non-negative numeric
    if (input_parameter_name %in% c("minimum_distance_between_cells") && !(is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter >= 0)) {
      check_value <- 4
      break
    }
    # Numeric between 0 and 1
    if (input_parameter_name %in% c("jitter_proportion") && !(is.numeric(input_parameter) && length(input_parameter) == 1 && input_parameter >= 0 && input_parameter <= 1)) {
      check_value <- 5
      break
    }
    # Character
    if (input_parameter_name %in% c("background_cell_type") && !(is.character(input_parameter)) && length(input_parameter) == 1) {
      check_value <- 6
      break
    }
    # Logical
    if (input_parameter_name %in% c("plot_image") && !(is.logical(input_parameter)) && length(input_parameter) == 1) {
      check_value <- 7
      break
    }
    # Character vector
    if (input_parameter_name %in% c("cell_types", "cluster_cell_types", "ring_cell_types", "inner_ring_cell_types", "outer_ring_cell_types") && 
        !(is.character(input_parameter))) {
      check_value <- 8
      break
    }
    # Numeric vector
    if (input_parameter_name %in% c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions") && 
        !(is.numeric(input_parameter))) {
      check_value <- 9
      break
    }
    # Numeric vector contains values between 0 and 1
    if (input_parameter_name %in% c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions") && 
        sum(input_parameter < 0 | input_parameter > 1) != 0) {
      check_value <- 10
      break
    }
    # Numeric vector contains values that sum to 1
    if (input_parameter_name %in% c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions") && 
        !is_equal_with_tolerance(sum(input_parameter), 1)) {
      check_value <- 11
      break
    }
    # Numeric vector of length 3
    if (input_parameter_name %in% c("centre_loc", "start_loc", "end_loc") && !(is.numeric(input_parameter) && length(input_parameter) == 3)) {
      check_value <- 12
      break
    }
    # Numeric
    if (input_parameter_name %in% c("y_z_rotation", "x_z_rotation", "x_y_rotation") && !(is.numeric(input_parameter) && length(input_parameter) == 1)) {
      check_value <- 13
      break
    }
  }
  
  # Two vectors match in length
  if (check_value == 0) {
    pairs <- data.frame(name1 = c("cell_types", "cluster_cell_types", "ring_cell_types", "inner_ring_cell_types", "outer_ring_cell_types"),
                        name2 = c("cell_proportions", "cluster_cell_proportions", "ring_cell_proportions", "inner_ring_cell_proportions", "outer_ring_cell_proportions"))  
    
    for (i in seq(nrow(pairs))) {
      name1 <- pairs[["name1"]][i]
      name2 <- pairs[["name2"]][i]
      if (name1 %in% input_parameter_names && name2 %in% input_parameter_names && length(input_parameters[[name1]]) != length(input_parameters[[name2]])) {
        check_value <- 14
        input_parameter_name <- c(name1, name2)
        break
      }
    }
  }
  
  # If check_value equals 0, all inputs are valid.
  if (check_value == 0) {
    return(TRUE)
  }
  # At least one input is not valid, return the first invalid input. 
  else {
    return(list(input_parameter_name = input_parameter_name, check_value = check_value)) 
  }
}

input_parameter_error_message <- function(input_parameter_check_value) {
  
  input_parameter_name <- input_parameter_check_value[[1]]
  check_value <- input_parameter_check_value[[2]]
  
  error_message <- switch(check_value,
                          "1" = paste(input_parameter_name, "is not a SpatialExperiment object."),
                          "2" = paste(input_parameter_name, "is not a positive integer."),
                          "3" = paste(input_parameter_name, "is not a positive numeric."),
                          "4" = paste(input_parameter_name, "is not a non-negative numeric."),
                          "5" = paste(input_parameter_name, "is not a numeric between 0 and 1."),
                          "6" = paste(input_parameter_name, "is not a character."),
                          "7" = paste(input_parameter_name, "is not a logical (TRUE or FALSE)."),
                          "8" = paste(input_parameter_name, "is not a character vector."),
                          "9" = paste(input_parameter_name, "is not a numeric vector."),
                          "10" = paste(input_parameter_name, "cannot be negative or greater than 1."),
                          "11" = paste(input_parameter_name, "does not sum to 1."),
                          "12" = paste(input_parameter_name, "is not a numeric vector of length 3."),
                          "13" = paste(input_parameter_name, "is not a numeric."),
                          "14" = paste(input_parameter_name[1], "and", input_parameter_name[2], "do not match in length."))
  
  return(error_message)
}