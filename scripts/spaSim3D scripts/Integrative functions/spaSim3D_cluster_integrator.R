### Example code to run -------------------------------------------------------
simulated_cluster_data <- spaSim3D_cluster_integrator(simulated_background_data)

#

### Message strings------------------------------------------------------------
message_no_simulated_data <- paste("Hello spaSim-3D user. Please input your simulated data into this function.\n",
                                   "If you don't have any, you can use the spaSim3D_background_integrator function")

message_shape_choice <- paste("Hello spaSim-3D user, what type of shape do you want your cluster to be?\n
          1. Sphere\n
          2. Ellipsoid\n
          3. Cylinder\n
          4. Network\n\n",
                              "To choose, please enter 1, 2, 3 or 4.\n", sep = "")

message_sphere_cluster <- paste("We will need a few parameters to generate a sphere cluster\n",
                                "    Radius\n",
                                "    Coordinates of sphere centre: x, y and z\n",
                                "If you want to change your inputs, you'll be able to at the end.\n", sep = "")

message_ellipsoid_cluster <- paste("We will need a few parameters to generate a ellipsoid cluster\n",
                                   "    Radii: x, y and z\n",
                                   "    Coordinates of ellipsoid centre: x, y and z\n",
                                   "    Angle of rotation in the x-axis, y-axis and z-axis\n",
                                   "If you want to change your inputs, you'll be able to at the end.\n", sep = "")

message_cylinder_cluster <- paste("We will need a few parameters to generate a cylinder cluster\n",
                                  "    Radius\n",
                                  "    Coordinates of the cylinder start point: x, y and z\n",
                                  "    Coordinates of the cylinder end point: x, y and z\n",
                                  "If you want to change your inputs, you'll be able to at the end.\n", sep = "")

message_network_cluster <- paste("We will need a few parameters to generate a network cluster\n",
                                 "    Number of branches\n",
                                 "    Width of each branch: x, y and z\n",
                                 "    Radius spanned by the whole network: x, y and z\n",
                                 "    Coordinates of network centre: x, y and z\n",
                                 "If you want to change your inputs, you'll be able to at the end.\n", sep = "")


message_get_cell_types <- "Keep entering the name of cell types you would like.\n    enter 'stop' to move on.\n"


message_cluster_choice <- paste("You can customise your cluster further if you'd like:\n
          1. Add a ring\n
          2. Add a double ring\n
          3. Continue\n\n",
                                "To choose, please enter 1, 2 or 3.\n", sep = "")


### Functions -----------------------------------------------------------------
spaSim3D_cluster_integrator <- function(simulated_data = NULL) {
  
  ## Start with checking if the user has inputted simulated_data (and not some garbage)
  if (is.null(simulated_data)) {
    message(message_no_simulated_data)
    return (NULL)
  }
  
  ## Plot the user's data so they can see what they already have
  fig <- plot_cell_categories3D(simulated_data)
  print(fig)
  
  ## Get user's choice for shape type (sphere, ellipsoid, cylinder or network)
  message(message_shape_choice)
  user_input_shape <- get_integer_input_from_options(1:4)
  
  ### Sphere
  if (user_input_shape == 1) {
    # Get required parameters for a sphere cluster from user
    message(message_sphere_cluster)
    parameter_values <- list("radius" = get_positive_numeric_input("radius"),
                             "centre x coordinate" = get_non_negative_numeric_input("centre x coordinate"),
                             "centre y coordinate" = get_non_negative_numeric_input("centre y coordinate"),
                             "centre z coordinate" = get_non_negative_numeric_input("centre z coordinate"))
    display_parameters(parameter_values)

    # Generate sphere cluster simulation using these parameters
    cluster_properties <- list(list(shape = "Sphere",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    radius = parameter_values[["radius"]],
                                    centre_loc = c(parameter_values[["centre x coordinate"]],
                                                   parameter_values[["centre y coordinate"]],
                                                   parameter_values[["centre z coordinate"]])))
    message("Generating simulation...")
    simulated_data_new <- simulate_clusters3D(simulated_data,
                                              1,
                                              cluster_properties,
                                              plot_image = TRUE,
                                              plot_categories = NULL,
                                              plot_colours = NULL)

    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {

      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))

      if (user_input_parameter_choice == 1) parameter_values[["radius"]] <- get_positive_numeric_input("radius")
      if (user_input_parameter_choice == 2) parameter_values[["centre x coordinate"]] <- get_non_negative_numeric_input("centre x coordinate")
      if (user_input_parameter_choice == 3) parameter_values[["centre y coordinate"]] <- get_non_negative_numeric_input("centre y coordinate")
      if (user_input_parameter_choice == 4) parameter_values[["centre z coordinate"]] <- get_non_negative_numeric_input("centre z coordinate")

      display_parameters(parameter_values)

      # Generate sphere cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "Sphere",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      radius = parameter_values[["radius"]],
                                      centre_loc = c(parameter_values[["centre x coordinate"]],
                                                     parameter_values[["centre y coordinate"]],
                                                     parameter_values[["centre z coordinate"]])))
      message("Generating simulation...")
      simulated_data_new <- simulate_clusters3D(simulated_data,
                                                1,
                                                cluster_properties,
                                                plot_image = TRUE,
                                                plot_categories = NULL,
                                                plot_colours = NULL)

      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Ellipsoid
  else if (user_input_shape == 2) {
    # Get required parameters for an ellipsoid cluster from user
    message(message_ellipsoid_cluster)
    parameter_values <- list("x radius" = get_positive_numeric_input("x radius"),
                             "y radius" = get_positive_numeric_input("y radius"),
                             "z radius" = get_positive_numeric_input("z radius"),
                             "centre x coordinate" = get_non_negative_numeric_input("centre x coordinate"),
                             "centre y coordinate" = get_non_negative_numeric_input("centre y coordinate"),
                             "centre z coordinate" = get_non_negative_numeric_input("centre z coordinate"),
                             "x-axis rotation angle" = get_non_negative_numeric_input("x-axis rotation angle"),
                             "y-axis rotation angle" = get_non_negative_numeric_input("y-axis rotation angle"),
                             "z-axis rotation angle" = get_non_negative_numeric_input("z-axis rotation angle"))
    display_parameters(parameter_values)
    
    # Generate ellipsoid cluster simulation using these parameters
    cluster_properties <- list(list(shape = "Ellipsoid",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    x_radius = parameter_values[["x radius"]],
                                    y_radius = parameter_values[["y radius"]],
                                    z_radius = parameter_values[["z radius"]],
                                    centre_loc = c(parameter_values[["centre x coordinate"]],
                                                   parameter_values[["centre y coordinate"]],
                                                   parameter_values[["centre z coordinate"]]),
                                    y_z_rotation = parameter_values[["x-axis rotation angle"]],
                                    x_z_rotation = parameter_values[["y-axis rotation angle"]],
                                    x_y_rotation = parameter_values[["z-axis rotation angle"]]))
    
    message("Generating simulation...")
    simulated_data_new <- simulate_clusters3D(simulated_data,
                                              1,
                                              cluster_properties,
                                              plot_image = TRUE,
                                              plot_categories = NULL,
                                              plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["x radius"]] <- get_positive_numeric_input("x radius")
      if (user_input_parameter_choice == 2) parameter_values[["y radius"]] <- get_positive_numeric_input("y radius")
      if (user_input_parameter_choice == 3) parameter_values[["z radius"]] <- get_positive_numeric_input("z radius")
      if (user_input_parameter_choice == 4) parameter_values[["centre x coordinate"]] <- get_non_negative_numeric_input("centre x coordinate")
      if (user_input_parameter_choice == 5) parameter_values[["centre y coordinate"]] <- get_non_negative_numeric_input("centre y coordinate")
      if (user_input_parameter_choice == 6) parameter_values[["centre z coordinate"]] <- get_non_negative_numeric_input("centre z coordinate")
      if (user_input_parameter_choice == 7) parameter_values[["x-axis rotation angle"]] <- get_non_negative_numeric_input("x-axis rotation angle")
      if (user_input_parameter_choice == 8) parameter_values[["y-axis rotation angle"]] <- get_non_negative_numeric_input("y-axis rotation angle")
      if (user_input_parameter_choice == 9) parameter_values[["z-axis rotation angle"]] <- get_non_negative_numeric_input("z-axis rotation angle")
      
      display_parameters(parameter_values)
      
      # Generate ellipsoid cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "Ellipsoid",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      x_radius = parameter_values[["x radius"]],
                                      y_radius = parameter_values[["y radius"]],
                                      z_radius = parameter_values[["z radius"]],
                                      centre_loc = c(parameter_values[["centre x coordinate"]],
                                                     parameter_values[["centre y coordinate"]],
                                                     parameter_values[["centre z coordinate"]]),
                                      y_z_rotation = parameter_values[["x-axis rotation angle"]],
                                      x_z_rotation = parameter_values[["y-axis rotation angle"]],
                                      x_y_rotation = parameter_values[["z-axis rotation angle"]]))
      
      message("Generating simulation...")
      simulated_data_new <- simulate_clusters3D(simulated_data,
                                                1,
                                                cluster_properties,
                                                plot_image = TRUE,
                                                plot_categories = NULL,
                                                plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Cylinder
  else if (user_input_shape == 3) {
    # Get required parameters for a cylinder cluster from user
    message(message_cylinder_cluster)
    parameter_values <- list("radius" = get_positive_numeric_input("radius"),
                             "start x coordinate" = get_non_negative_numeric_input("start x coordinate"),
                             "start y coordinate" = get_non_negative_numeric_input("start y coordinate"),
                             "start z coordinate" = get_non_negative_numeric_input("start z coordinate"),
                             "end x coordinate" = get_non_negative_numeric_input("end x coordinate"),
                             "end y coordinate" = get_non_negative_numeric_input("end y coordinate"),
                             "end z coordinate" = get_non_negative_numeric_input("end z coordinate"))
    display_parameters(parameter_values)
    
    # Generate cylinder cluster simulation using these parameters
    cluster_properties <- list(list(shape = "Cylinder",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    radius = parameter_values[["radius"]],
                                    start_loc = c(parameter_values[["start x coordinate"]],
                                                  parameter_values[["start y coordinate"]],
                                                  parameter_values[["start z coordinate"]]),
                                    end_loc = c(parameter_values[["end x coordinate"]],
                                                parameter_values[["end y coordinate"]],
                                                parameter_values[["end z coordinate"]])))
    message("Generating simulation...")
    simulated_data_new <- simulate_clusters3D(simulated_data,
                                              1,
                                              cluster_properties,
                                              plot_image = TRUE,
                                              plot_categories = NULL,
                                              plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))
      
      if (user_input_parameter_choice == 1) parameter_values[["radius"]] <- get_positive_numeric_input("radius")
      if (user_input_parameter_choice == 2) parameter_values[["start x coordinate"]] <- get_non_negative_numeric_input("start x coordinate")
      if (user_input_parameter_choice == 3) parameter_values[["start y coordinate"]] <- get_non_negative_numeric_input("start y coordinate")
      if (user_input_parameter_choice == 4) parameter_values[["start z coordinate"]] <- get_non_negative_numeric_input("start z coordinate")
      if (user_input_parameter_choice == 5) parameter_values[["end x coordinate"]] <- get_non_negative_numeric_input("end x coordinate")
      if (user_input_parameter_choice == 6) parameter_values[["end y coordinate"]] <- get_non_negative_numeric_input("end y coordinate")
      if (user_input_parameter_choice == 7) parameter_values[["end z coordinate"]] <- get_non_negative_numeric_input("end z coordinate")
      
      display_parameters(parameter_values)
      
      # Generate cylinder cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "Cylinder",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      radius = parameter_values[["radius"]],
                                      start_loc = c(parameter_values[["start x coordinate"]],
                                                    parameter_values[["start y coordinate"]],
                                                    parameter_values[["start z coordinate"]]),
                                      end_loc = c(parameter_values[["end x coordinate"]],
                                                  parameter_values[["end y coordinate"]],
                                                  parameter_values[["end z coordinate"]])))
      message("Generating simulation...")
      simulated_data_new <- simulate_clusters3D(simulated_data,
                                                1,
                                                cluster_properties,
                                                plot_image = TRUE,
                                                plot_categories = NULL,
                                                plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  ### Network
  else if (user_input_shape == 4) {
    # Get required parameters for a network cluster from user
    message(message_network_cluster)
    parameter_values <- list("number of branches" = get_integer_greater_than_or_equal_input("number of branches", 2),
                             "width of branch" = get_positive_numeric_input("width of branch"),
                             "radius spanned by network" = get_positive_numeric_input("radius spanned by network"),
                             "centre x coordinate" = get_non_negative_numeric_input("centre x coordinate"),
                             "centre y coordinate" = get_non_negative_numeric_input("centre y coordinate"),
                             "centre z coordinate" = get_non_negative_numeric_input("centre z coordinate"))
    display_parameters(parameter_values)
    
    # Generate network cluster simulation using these parameters
    cluster_properties <- list(list(shape = "Network",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    n_edges = parameter_values[["number of branches"]],
                                    width = parameter_values[["width of branch"]],
                                    radius = parameter_values[["radius spanned by network"]],
                                    centre_loc = c(parameter_values[["centre x coordinate"]],
                                                   parameter_values[["centre y coordinate"]],
                                                   parameter_values[["centre z coordinate"]])))
    message("Generating simulation...")
    simulated_data_new <- simulate_clusters3D(simulated_data,
                                              1,
                                              cluster_properties,
                                              plot_image = TRUE,
                                              plot_categories = NULL,
                                              plot_colours = NULL)
    
    # Allow user the option to change their input parameters
    message("Would you like to change your input parameters?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      user_input_parameter_choice <- get_integer_input_from_options(seq(length(parameter_values)))

      if (user_input_parameter_choice == 1) parameter_values[["number of branches"]] <- get_integer_greater_than_or_equal_input("number of branches", 2)
      if (user_input_parameter_choice == 2) parameter_values[["width of branch"]] <- get_positive_numeric_input("width of branch")
      if (user_input_parameter_choice == 3) parameter_values[["radius spanned by network"]] <- get_positive_numeric_input("radius spanned by network")
      if (user_input_parameter_choice == 4) parameter_values[["centre x coordinate"]] <- get_non_negative_numeric_input("centre x coordinate")
      if (user_input_parameter_choice == 5) parameter_values[["centre y coordinate"]] <- get_non_negative_numeric_input("centre y coordinate")
      if (user_input_parameter_choice == 6) parameter_values[["centre z coordinate"]] <- get_non_negative_numeric_input("centre z coordinate")
      
      display_parameters(parameter_values)
      
      # Generate sphere cluster simulation using updated parameters
      cluster_properties <- list(list(shape = "Network",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      n_edges = parameter_values[["number of branches"]],
                                      width = parameter_values[["width of network branch"]],
                                      radius = parameter_values[["radius spanned by network"]],
                                      centre_loc = c(parameter_values[["centre x coordinate"]],
                                                     parameter_values[["centre y coordinate"]],
                                                     parameter_values[["centre z coordinate"]])))
      message("Generating simulation...")
      simulated_data_new <- simulate_clusters3D(simulated_data,
                                                1,
                                                cluster_properties,
                                                plot_image = TRUE,
                                                plot_categories = NULL,
                                                plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }

  # Allow user to change the cell composition of the cluster
  message("Let's change the cell composition of this cluster")
  simulated_data_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_data_new,
                                                                                   simulate_clusters3D,
                                                                                   cluster_properties,
                                                                                   "cluster_cell_types",
                                                                                   "cluster_cell_proportions",
                                                                                   "Cluster")
  simulated_data_new <- simulated_data_new_and_properties[["data"]]
  cluster_properties <- simulated_data_new_and_properties[["properties"]]
  
  
  
  ## Get user's choice for cluster type (ringed, double ringed or continue)
  message(message_cluster_choice)
  user_input_cluster <- get_integer_input_from_options(1:3)
  
  ### Ring
  if (user_input_cluster == 1) {
    # Get width of ring from user
    message("For a single ring, we needs its width.\n")
    
    # Generate cluster with ring simulation using this width
    cluster_properties[[1]][["ring_width"]] <- get_positive_numeric_input("ring width")
    cluster_properties[[1]][["ring_cell_types"]] <- c("Ring")
    cluster_properties[[1]][["ring_cell_proportions"]] <- 1
    
    message("Generating simulation...")
    simulated_data_new <- simulate_rings3D(simulated_data,
                                           1,
                                           cluster_properties,
                                           plot_image = TRUE,
                                           plot_categories = NULL,
                                           plot_colours = NULL)
    
    # Allow user the option to change the ring width
    message("Would you like to change the ring width?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      cluster_properties[[1]][["ring_width"]] <- get_positive_numeric_input("ring width")
      
      message("Generating simulation...")
      simulated_data_new <- simulate_rings3D(simulated_data,
                                             1,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_categories = NULL,
                                             plot_colours = NULL)
      
      message("Would you like to the ring width?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
    
    # Allow user to change the cell composition of the ring
    message("Let's change the cell composition of the ring")
    
    simulated_data_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_data_new,
                                                                                     simulate_rings3D,
                                                                                     cluster_properties,
                                                                                     "ring_cell_types",
                                                                                     "ring_cell_proportions",
                                                                                     "Ring")
    
    simulated_data_new <- simulated_data_new_and_properties[["data"]]
  }
  ### Double ring
  else if (user_input_cluster == 2) {
    # Get width of inner and outer ring from user
    message("For a double ring, we needs the width of the inner and outer ring.\n")
    
    # Generate cluster with double ring simulation using both widths
    parameter_values <- list("inner ring width" = get_positive_numeric_input("inner ring width"),
                             "outer ring width" = get_positive_numeric_input("outer ring width"))
    display_parameters(parameter_values)
    
    cluster_properties[[1]][["inner_ring_width"]] <- parameter_values[["inner ring width"]]
    cluster_properties[[1]][["outer_ring_width"]] <- parameter_values[["outer ring width"]]
    cluster_properties[[1]][["inner_ring_cell_types"]] <- c("Inner ring")
    cluster_properties[[1]][["inner_ring_cell_proportions"]] <- 1
    cluster_properties[[1]][["outer_ring_cell_types"]] <- c("Outer ring")
    cluster_properties[[1]][["outer_ring_cell_proportions"]] <- 1
    
    message("Generating simulation...")
    simulated_data_new <- simulate_double_rings3D(simulated_data,
                                                  1,
                                                  cluster_properties,
                                                  plot_image = TRUE,
                                                  plot_categories = NULL,
                                                  plot_colours = NULL)
    
    # Allow user the option to change the widths of the inner or outer ring
    message("Would you like to change the widths of the inner or outer ring?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      if (user_input_parameter_choice == 1) parameter_values[["inner ring width"]] <- get_positive_numeric_input("inner ring width")
      if (user_input_parameter_choice == 2) parameter_values[["outer ring width"]] <- get_positive_numeric_input("outer ring width")
      
      cluster_properties[[1]][["inner_ring_width"]] <- parameter_values[["inner ring width"]]
      cluster_properties[[1]][["outer_ring_width"]] <- parameter_values[["outer ring width"]]
      
      display_parameters(parameter_values)
      
      message("Generating simulation...")
      simulated_data_new <- simulate_double_rings3D(simulated_data,
                                                    1,
                                                    cluster_properties,
                                                    plot_image = TRUE,
                                                    plot_categories = NULL,
                                                    plot_colours = NULL)
      
      message("Would you like to change the widths of the inner or outer ring?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
    
    # Allow user to change the cell composition of the inner ring
    message("Let's change the cell composition of the inner ring")
    simulated_data_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_data_new,
                                                                                     simulate_double_rings3D,
                                                                                     cluster_properties,
                                                                                     "inner_ring_cell_types",
                                                                                     "inner_ring_cell_proportions",
                                                                                     "Inner ring")
    
    simulated_data_new <- simulated_data_new_and_properties[["data"]]
    cluster_properties <- simulated_data_new_and_properties[["properties"]]
    
    # Allow user to change the cell composition of the outer ring
    message("Let's change the cell composition of the outer ring")
    simulated_data_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_data_new,
                                                                                     simulate_double_rings3D,
                                                                                     cluster_properties,
                                                                                     "outer_ring_cell_types",
                                                                                     "outer_ring_cell_proportions",
                                                                                     "Outer ring")
    
    simulated_data_new <- simulated_data_new_and_properties[["data"]]
  }
  ### Continue
  else if (user_input_cluster == 3) {
    
  }
 
  message("All done!")
  return (simulated_data_new) 
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
  
  return (int_value)
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
  
  return (positive_numeric_value)
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
  
  return (integer_value)
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
  
  return (non_negative_value)
}

get_cell_types_and_proportions_for_clusters <- function(simulated_data, simulate_function, properties, cell_type_option, cell_proportion_option, temp_cell_type) {
  
  ## Display the cell types currently found in simulated_data to the user
  current_cell_types <- setdiff(unique(simulated_data[["Cell.Type"]]), temp_cell_type)
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
      
      ## Generate simulation
      message("Generating simulation...")
      simulated_data_new <- simulated_data
      
      ## Just change the cell type of the temp_cell_type, no need to actually re-simulate
      for (i in seq(nrow(simulated_data_new))) {
        if (simulated_data_new[i, "Cell.Type"] == temp_cell_type) {
          simulated_data_new[i, "Cell.Type"] <- get_random_cell_type(cell_types, cell_proportions)
        }
      }
      fig <- plot_cell_categories3D(simulated_data_new)
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
  
  return (list(data = simulated_data_new, properties = properties))
}

get_random_cell_type <- function(cell_types, cell_proportions) {
  random <- runif(n = 1, min = 0, max = 1)
  
  i <- 1
  current_proportion <- 0
  
  while (i <= length(cell_types)){
    current_proportion <- current_proportion + cell_proportions[i]
    if (random <= current_proportion) {
      return (cell_types[i])
    }
    i <- i + 1
  }
}
