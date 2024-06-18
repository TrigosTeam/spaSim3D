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