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
