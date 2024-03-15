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
cell_dens <- find_grid_square_densities(data = data,
                                        nrows = nrows, ncols = ncols,
                                        length = grid_square_length,
                                        width = grid_square_width)



library(ggplot2)

## Plot all cells
plot_cells(data = data, length = length, width = width, nrows = nrows, ncols = ncols)



## Plot only showing a specific cell type
cell_type <- "Tumour"
data_specific_cell <- data[data$Cell.Type == cell_type, ]

plot_cells(data = data_specific_cell, length = length, width = width, nrows = nrows, ncols = ncols)




## Other stuff
summary(cell_props[["Tumour"]][cell_props[["Tumour"]] > 0.0])
hist(cell_props[["Tumour"]][cell_props[["Tumour"]] > 0.0])

summary(cell_dens[["Tumour"]][cell_dens[["Tumour"]] > 0.0])



## Plot squares
x <- c(15, 4, 3)
y <- c(8, 2, 10)
side_length <- c(2, 4, 6)
squares_data <- data.frame(x = x, y = y, side_length = side_length)

plot_squares(squares_data = squares_data,
             xmin = 0,
             xmax = 20,
             ymin = 0,
             ymax = 20)





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


find_grid_square_densities <- function(data, 
                                       nrows, ncols,
                                       length, width) {
  
  cell_types <- unique(data$Cell.Type)
  cell_dens <- vector(mode = 'list', length = length(cell_types))
  names(cell_dens) <- cell_types
  
  for (i in seq(nrows * ncols)) {
    cells <- data$Cell.Type[data$Cell.Num == i]
    
    for (cell_type in cell_types) {
      cell_dens[[cell_type]] <- append(cell_dens[[cell_type]], sum(cells == cell_type) / (length * width)) 
      
    }
  }  
  
  return (cell_dens)
}


plot_cells <- function(data, length, width, nrows, ncols) {
  
  ggplot(data = data,
         mapping = aes(x = Cell.X.Position, y = Cell.Y.Position, color = Cell.Type)) + 
    geom_point() + 
    theme(
      panel.background = element_rect(fill = NA),
      panel.ontop = TRUE,
      panel.grid.major=element_line(color = "black"),
      panel.grid.minor=element_line(color = "black")
    ) + 
    scale_x_continuous(minor_breaks = seq(0, width, width/ncols), limits = c(0, length)) +
    scale_y_continuous(minor_breaks = seq(0, length, length/nrows), limits = c(0, width)) 
  
}

plot_squares <- function(squares_data, xmin, xmax, ymin, ymax) {
  
  ggplot(squares_data, aes(x = x, y = y)) +
    geom_rect(aes(xmin = x, ymin = y, 
                  xmax = x + side_length, ymax = y + side_length),
              fill = "green",
              color = "black") +
    xlim(xmin, xmax) +
    ylim(ymin, ymax)
  
}

