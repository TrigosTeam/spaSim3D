## Dimensions of study area
length <- 2000
width <- 2000

## Define grid dimensions
nrows <- 8
ncols <- 8
ngrid_squares <- nrows * ncols

grid_square_length <- length / nrows
grid_square_width <- width / ncols


### Algorithm to see proportion of cells in each grid square, show as histogram
data <- bgCluster

data <- determine_grid_square_nums(data = data,
                                   length = length, width = width,
                                   nrows = nrows, ncols = ncols)

cell_props <- find_grid_square_props(data = data,
                                     nrows = nrows, ncols = ncols)



library(ggplot2)

ggplot(data = data,
       mapping = aes(x = Cell.X.Position, y = Cell.Y.Position, color = Cell.Type)) + 
  geom_point() + 
  theme(
    panel.background = element_rect(fill = NA),
    panel.ontop = TRUE,
    panel.grid.major=element_line(color = "black"),
    panel.grid.minor=element_line(color = "black")
  ) + 
  scale_x_continuous(minor_breaks = seq(0, width, grid_square_width)) +
  scale_y_continuous(minor_breaks = seq(0, length, grid_square_length))


summary(cell_props[["Tumour"]][cell_props[["Tumour"]] > 0.0])
hist(cell_props[["Tumour"]][cell_props[["Tumour"]] > 0.0])



### Functions

determine_grid_square_nums <- function(data,
                                       length, width,
                                       nrows, ncols) {
  
  
  data$Cell.Num <- ncols * floor(data$Cell.Y.Position / (length / nrows)) +
                   floor(data$Cell.X.Position / (width / ncols)) + 1
  
  return(data)
}


find_grid_square_props <- function(data, 
                                   nrows, ncols) {

  cell_types <- unique(data$Cell.Type)
  cell_props <- vector(mode = 'list', length = length(cell_types))
  names(cell_props) <- cell_types
  
  for (i in seq(nrows * ncols)) {
    cells <- data$Cell.Type[data$Cell.Num == i]
    
    for (cell_type in cell_types) {
      cell_props[[cell_type]] <- append(cell_props[[cell_type]], sum(cells == cell_type) / length(cells)) 
      
    }
  }  
  
  return (cell_props)
}
