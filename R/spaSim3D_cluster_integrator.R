spaSim3D_cluster_integrator <- function(simulated_spe = NULL) {

  ### Message strings
  message_no_simulated_spe <- paste("Hello spaSim-3D user. Please input your simulated spe object into this function.\n",
                                    "If you don't have any, you can use the spaSim3D_background_integrator function")
  
  message_shape_choice <- paste("Hello spaSim-3D user, hopefully you can see a plot of your current spe object. What type of shape do you want your cluster to be?\n
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
  
  
  message_cluster_choice <- paste("You can customise your cluster further if you'd like:\n
          1. Add a ring\n
          2. Add a double ring\n
          3. Continue\n\n",
                                  "To choose, please enter 1, 2 or 3.\n", sep = "")
  
    
  ## Start with checking if the user has inputted spe object
  if (class(simulated_spe) != "SpatialExperiment") {
    stop(message_no_simulated_spe)
  }
  
  ## Plot the user's data so they can see what they already have
  fig <- plot_cells3D(simulated_spe)
  methods::show(fig)
  
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
    cluster_properties <- list(list(shape = "sphere",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    radius = parameter_values[["radius"]],
                                    centre_loc = c(parameter_values[["centre x coordinate"]],
                                                   parameter_values[["centre y coordinate"]],
                                                   parameter_values[["centre z coordinate"]])))
    message("Generating simulation...")
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
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
      cluster_properties <- list(list(shape = "sphere",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      radius = parameter_values[["radius"]],
                                      centre_loc = c(parameter_values[["centre x coordinate"]],
                                                     parameter_values[["centre y coordinate"]],
                                                     parameter_values[["centre z coordinate"]])))
      message("Generating simulation...")
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
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
    cluster_properties <- list(list(shape = "ellipsoid",
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
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
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
      cluster_properties <- list(list(shape = "ellipsoid",
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
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
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
    cluster_properties <- list(list(shape = "cylinder",
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
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
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
      cluster_properties <- list(list(shape = "cylinder",
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
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
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
    cluster_properties <- list(list(shape = "network",
                                    cluster_cell_types = "Cluster",
                                    cluster_cell_proportions = 1,
                                    n_edges = parameter_values[["number of branches"]],
                                    width = parameter_values[["width of branch"]],
                                    radius = parameter_values[["radius spanned by network"]],
                                    centre_loc = c(parameter_values[["centre x coordinate"]],
                                                   parameter_values[["centre y coordinate"]],
                                                   parameter_values[["centre z coordinate"]])))
    message("Generating simulation...")
    simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                             cluster_properties,
                                             plot_image = TRUE,
                                             plot_cell_types = NULL,
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
      cluster_properties <- list(list(shape = "network",
                                      cluster_cell_types = "Cluster",
                                      cluster_cell_proportions = 1,
                                      n_edges = parameter_values[["number of branches"]],
                                      width = parameter_values[["width of branch"]],
                                      radius = parameter_values[["radius spanned by network"]],
                                      centre_loc = c(parameter_values[["centre x coordinate"]],
                                                     parameter_values[["centre y coordinate"]],
                                                     parameter_values[["centre z coordinate"]])))
      
      message("Generating simulation...")
      simulated_spe_new <- simulate_clusters3D(simulated_spe,
                                               cluster_properties,
                                               plot_image = TRUE,
                                               plot_cell_types = NULL,
                                               plot_colours = NULL)
      
      message("Would you like to change your inputs?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
  }
  
  # Allow user to change the cell composition of the cluster
  message("Let's change the cell composition of this cluster")
  simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                  simulate_clusters3D,
                                                                                  cluster_properties,
                                                                                  "cluster_cell_types",
                                                                                  "cluster_cell_proportions",
                                                                                  "Cluster")
  simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
  cluster_properties <- simulated_spe_new_and_properties[["properties"]]
  
  
  
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
    simulated_spe_new <- simulate_rings3D(simulated_spe,
                                          cluster_properties,
                                          plot_image = TRUE,
                                          plot_cell_types = NULL,
                                          plot_colours = NULL)
    
    # Allow user the option to change the ring width
    message("Would you like to change the ring width?\n")
    change_input_parameters_y_or_n <- get_y_or_n_input()
    while (change_input_parameters_y_or_n == "y") {
      
      # Determine which parameter the user wants to change
      cluster_properties[[1]][["ring_width"]] <- get_positive_numeric_input("ring width")
      
      message("Generating simulation...")
      simulated_spe_new <- simulate_rings3D(simulated_spe,
                                            cluster_properties,
                                            plot_image = TRUE,
                                            plot_cell_types = NULL,
                                            plot_colours = NULL)
      
      message("Would you like to the ring width?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
    
    # Allow user to change the cell composition of the ring
    message("Let's change the cell composition of the ring")
    
    simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                    simulate_rings3D,
                                                                                    cluster_properties,
                                                                                    "ring_cell_types",
                                                                                    "ring_cell_proportions",
                                                                                    "Ring")
    
    simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
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
    simulated_spe_new <- simulate_double_rings3D(simulated_spe,
                                                 cluster_properties,
                                                 plot_image = TRUE,
                                                 plot_cell_types = NULL,
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
      simulated_spe_new <- simulate_double_rings3D(simulated_spe,
                                                   cluster_properties,
                                                   plot_image = TRUE,
                                                   plot_cell_types = NULL,
                                                   plot_colours = NULL)
      
      message("Would you like to change the widths of the inner or outer ring?\n")
      change_input_parameters_y_or_n <- get_y_or_n_input()
    }
    
    # Allow user to change the cell composition of the inner ring
    message("Let's change the cell composition of the inner ring")
    simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                    simulate_double_rings3D,
                                                                                    cluster_properties,
                                                                                    "inner_ring_cell_types",
                                                                                    "inner_ring_cell_proportions",
                                                                                    "Inner ring")
    
    simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
    cluster_properties <- simulated_spe_new_and_properties[["properties"]]
    
    # Allow user to change the cell composition of the outer ring
    message("Let's change the cell composition of the outer ring")
    simulated_spe_new_and_properties <- get_cell_types_and_proportions_for_clusters(simulated_spe_new,
                                                                                    simulate_double_rings3D,
                                                                                    cluster_properties,
                                                                                    "outer_ring_cell_types",
                                                                                    "outer_ring_cell_proportions",
                                                                                    "Outer ring")
    
    simulated_spe_new <- simulated_spe_new_and_properties[["data"]]
  }
  ### Continue
  else if (user_input_cluster == 3) {
    
  }
  
  message("All done!")
  return(simulated_spe_new) 
}
