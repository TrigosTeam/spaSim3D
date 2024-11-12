poisson_distribution3D <- function(n_cells, length, width, height)  {
  
  # Choose lambda
  lambda <- 5
  
  # Set number of rows, columns and layers
  nRows <- nCols <- nLays <- round((n_cells/lambda)^(1/3))
  
  # Get number of cubes in grid
  nCubes <- nRows * nCols * nLays
  
  # Get pois vector
  pois <- rpois(nCubes, lambda)
  
  # Get points for each prism region
  x <- c()
  y <- c()
  z <- c()
  
  for (row in seq(nRows)) {
    
    for (col in seq(nCols)) {
      
      for (lay in seq(nLays)) {
        current_cube_index <- nRows^2 * (row - 1) + nCols * (col - 1) + lay
        
        x <- append(x, runif(pois[current_cube_index], row - 1, row))
        y <- append(y, runif(pois[current_cube_index], col - 1, col))
        z <- append(z, runif(pois[current_cube_index], lay - 1, lay))
      }
    }
  }
  x <- x * length / nRows
  y <- y * width / nCols
  z <- z * height / nLays
  
  df <- data.frame("Cell.X.Position" = x, 
                   "Cell.Y.Position" = y, 
                   "Cell.Z.Position" = z)
  
  return(df)
}

## Prim's algorithm function
# Input is the adjacency matrix of the graph (i.e. output from -1 * apcluster::negDistMat(df of coords))
prims_algorithm <- function(graph) {
  
  # Number of vertices is number of points
  num_vertices <- nrow(graph)
  
  # Start with no vertices selected except first
  selected <- rep(FALSE, num_vertices)
  selected[1] <- TRUE
  
  # Create tree_edge matrix. Currently zero, each row represents the two vertices the edge joins
  tree_edges <- matrix(0, 
                       nrow = num_vertices - 1,
                       ncol = 2)
  
  # Iterate until we select enough edges (one less than the number of vertices for a MST)
  num_edges <- 0
  while (num_edges < num_vertices - 1) {
    # Set initial temp values for weight and vertex
    min_weight <- Inf
    min_vertex <- -1
    
    # Iterate through each currently selected vertex
    for (i in seq(num_vertices)) {
      
      # Found a currently selected vertex
      if (selected[i] == TRUE) {
        
        # Iterate through each unselected vertex and find the nearest one
        for (j in seq(num_vertices)) {
          if (!selected[j] && graph[i, j] < min_weight) {
            min_weight <- graph[i, j]
            min_vertex <- j
            curr_vertex <- i
          }
        }
      }
    }
    
    # Current edge connects the min_vertex and curr_vertex
    tree_edges[num_edges + 1, ] <- c(min_vertex, curr_vertex)
    selected[min_vertex] <- TRUE
    num_edges <- num_edges + 1
  }
  return(tree_edges)
}

get_tree_depth <- function(tree_edges) {
  
  tree_edges <- data.frame(tree_edges)
  colnames(tree_edges) <- c("vertex1", "vertex2")
  
  # Set the initial depth of each tree_edge to be NA.
  tree_edges$depth <- NA
  
  # Get vertices on the 'outskirts' of MST (leaf_vertices which have a depth of 1)
  tree_vertices <- c(tree_edges[ , 1], tree_edges[ , 2])
  leaf_vertices <- as.numeric(names(table(tree_vertices))[table(tree_vertices) == 1])
  
  # Start with leaf_vertices
  curr_vertices <- leaf_vertices
  curr_depth <- 1
  
  while (NA %in% tree_edges$depth) {
    
    # New vertices will be those adjacent to the current vertices
    new_vertices <- c()
    
    # Check each current vertex
    for (vertex in curr_vertices) {
      # Start with vertex1
      curr_edges <- which(tree_edges$vertex1 == vertex)
      tree_edges[curr_edges, "depth"][is.na(tree_edges[curr_edges, "depth"])] <- curr_depth
      new_vertices <- c(new_vertices, tree_edges[curr_edges, "vertex2"])
      
      # Then vertex2
      curr_edges <- which(tree_edges$vertex2 == vertex)
      tree_edges[curr_edges, "depth"][is.na(tree_edges[curr_edges, "depth"])] <- curr_depth
      new_vertices <- c(new_vertices, tree_edges[curr_edges, "vertex1"])
      
      # Only keep unique vertices
      new_vertices <- unique(new_vertices)
    }
    
    curr_depth <- curr_depth + 1
    curr_vertices <- new_vertices
  }
  
  return(tree_edges)
}

check_input_parameters <- function(input_parameters) {
  
  input_parameter_names <- names(input_parameters)
  
  check_value <- 0
  
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

is_equal_with_tolerance <- function(x, y, tolerance = 1e-6) {
  abs(x - y) <= tolerance
}

get_integer_greater_than_or_equal_input <- function(parameter, lower) {
  
  prompt <- paste("Enter an integer in value greater than or equal to ", lower, " for the ", parameter, ": ", sep = "")
  
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
    # Integer but below the lower bound
    else if (integer_value < lower) {
      message("Out of bounds input. Please a number greater than or equal to ", lower)
    }
    else {
      valid_input <- TRUE
      message("Valid input received!")
    }
  }
  
  return(integer_value)
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
    } 
    else {
      message(invalid_input_message)
    }
  }
  
  return(int_value)
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

get_numeric_between_input <- function(parameter, lower, upper) {
  
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

get_cell_types_and_proportions_for_clusters <- function(simulated_spe, simulate_function, properties, cell_type_option, cell_proportion_option, temp_cell_type) {
  
  message_get_cell_types <- "Keep entering the name of cell types you would like (e.g. Tumour, Immune, etc.).\n    enter 'stop' to move on."
  
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
      properties[[1]][[cell_proportion_option]] <- cell_proportions
      
      ## Convert spe object to data frame
      df <- data.frame(spatialCoords(simulated_spe), "Cell.Type" = simulated_spe[["Cell.Type"]])
      
      ## Just change the cell type of the temp_cell_type, no need to actually re-simulate
      df[["Cell.Type"]] <- ifelse(df[["Cell.Type"]] == temp_cell_type, 
                                  sample(cell_types, size = length(df[["Cell.Type"]]), replace = TRUE, prob = cell_proportions), 
                                  df[["Cell.Type"]])
      
      # Add Cell.ID column to data frame
      df$Cell.ID <- paste("Cell", seq(nrow(df)), sep = "_")
      
      # Update current meta data
      metadata <- simulated_spe@metadata
      metadata[["simulation"]][[length(metadata[["simulation"]])]][[cell_type_option]] <- cell_types
      metadata[["simulation"]][[length(metadata[["simulation"]])]][[cell_proportion_option]] <- cell_proportions
      
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
