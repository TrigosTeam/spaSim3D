calculate_cross_functions <- function(spe_object, 
                                      method = "Kcross", 
                                      cell_types_of_interest, 
                                      feature_colname, 
                                      plot_results = TRUE, 
                                      dist = NULL) {
  #CHECK
  formatted_data <-get_colData(spe_object)
  if (!all(cell_types_of_interest %in% formatted_data[[feature_colname]])) {
    stop("Cell type not found!")
  }
  
  # format spe to ppp object
  ppp_object <- format_spe_to_ppp(spe_object, feature_colname = feature_colname)
  
  ppp_object$marks <- as.factor(ppp_object$marks)
  
  
  
  # r
  if (is.null(dist)) {
    r <- NULL
  }
  
  else {
    r <- seq(0, dist, length.out = 100)
  }
  
  
  if (method == "Gcross"){
    p <- Gcross(ppp_object, cell_types_of_interest[1],
                cell_types_of_interest[2],correction = "border", r = r)
    if(plot_results){
      plot(p, main = paste("cross G function",attr(spe_object,"name")))
    }
  }
  
  
  else if (method == "Kcross"){
    p <- spatstat.explore::Kcross(ppp_object, 
                                  cell_types_of_interest[1], 
                                  cell_types_of_interest[2], 
                                  correction = "border", 
                                  r = r)
    
    
    if (plot_results) {
      if (is.null(dist)) {
        plot(p, main = paste("cross K function", attr(spe_object,"name")))
      }
      else {
        plot(p, main = paste("cross K function", attr(spe_object,"name")), xlim = c(0, dist))
      }
    }
  }
  
  
  else if (method == "Kcross.inhom"){
    p <- Kcross.inhom(ppp_object, cell_types_of_interest[1],
                      cell_types_of_interest[2],correction = "border", 
                      r = r)
    if(plot_results){
      if (is.null(dist)) plot(p, main = paste("cross K function",
                                              attr(spe_object,"name")))
      else plot(p, main = paste("cross K function",
                                attr(spe_object,"name")), 
                xlim = c(0,dist))
    }
  }
  else if (method == "Lcross"){
    p <- Lcross(ppp_object, cell_types_of_interest[1],
                cell_types_of_interest[2],correction = "border", r = r)
    if(plot_results){
      plot(p, main = paste("cross L function",attr(spe_object,"name")))
    }
  }
  else if (method == "Jcross"){
    p <- Jcross(ppp_object, cell_types_of_interest[1],
                cell_types_of_interest[2],correction = "border", r = r)
    if(plot_results){
      plot(p, main = paste("cross J function",attr(spe_object,"name")))
    }
  }
  
  
  return(p)
}