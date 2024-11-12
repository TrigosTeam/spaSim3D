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
