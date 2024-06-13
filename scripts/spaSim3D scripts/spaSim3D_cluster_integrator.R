### Example code to run -------------------------------------------------------
simulated_cluster_data <- spaSim3D_cluster_integrator(simulated_background_data)

#

### Message strings------------------------------------------------------------
message_no_simulated_data <- paste("Hello spaSim-3D user. Please input a ___ into this function.\n",
                                   "If you don't have any, you can use the spaSim3D_background_integrator function")

message_cluster_choice <- paste("Hello spaSim-3D user, what type of cluster do you want to add?\n
          1. Regular cluster\n
          2. Ringed cluster\n
          3. Double ringed cluster\n\n",
                                "To choose, please enter 1, 2 or 3.\n", sep = "")

message_shape_choice <- paste("What type of shape do you want your cluster to be?\n
          1. Sphere\n
          2. Ellipsoid\n
          3. Cylinder\n
          4. Network\n\n",
                              "To choose, please enter 1, 2 or 3.\n", sep = "")


### Functions -----------------------------------------------------------------
spaSim3D_cluster_integrator <- function(simulated_data) {
  ## Start with checking if the user has inputted simulated_data (and not some garbage)
  
  if (FALSE) {
    message(message_no_simulated_data)
    return (NULL)
  }
  
  ## Get user's choice for cluster type (regular, ringed or double ringed)
  message(message_cluster_choice)
  user_input_cluster <- get_integer_input_from_options(1:3)
  
  ## Get user's choice for shape type (sphere, ellipsoid, cylinder or network)
  message(message_shape_choice)
  user_input_shape <- get_integer_input_from_options(1:4)
  
  
  ### Regular cluster
  if (user_input_cluster == 1) {
    ### Sphere
    if (user_input_shape == 1) {
      print("Sphere cluster is work in progress")    
    }
    ### Ellipsoid
    else if (user_input_shape == 2) {
      print("Ellipsoid cluster is work in progress")    
    }
    ### Cylinder
    else if (user_input_shape == 3) {
      print("Cylinder cluster is work in progress")    
    }
    ### Network
    else if (user_input_shape == 4) {
      print("Network cluster is work in progress")    
    }
  }
  
  ### Ringed cluster
  else if (user_input_cluster == 2) {
    print("Ringed cluster is work in progress")    
  }
  
  ### Double ringed cluster
  else if (user_input_cluster == 3) {
    print("Double ringed cluster is work in progress")
  }
  
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
