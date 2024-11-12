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
