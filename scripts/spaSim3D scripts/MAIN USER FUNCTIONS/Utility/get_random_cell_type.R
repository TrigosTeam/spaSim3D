get_random_cell_type <- function(cell_types, cell_proportions) {
  
  random <- runif(n = 1, min = 0, max = 1)
  i <- 1
  current_proportion <- 0
  
  while (i <= length(cell_types)){
    current_proportion <- current_proportion + cell_proportions[i]
    if (random <= current_proportion) {
      return(cell_types[i])
    }
    i <- i + 1
  }
}