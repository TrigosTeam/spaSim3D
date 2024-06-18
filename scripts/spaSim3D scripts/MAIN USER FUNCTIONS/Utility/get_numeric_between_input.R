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