## Dimensions of study area
length <- 2000
width <- 2000

## Define grid dimensions
nrows <- 8
ncols <- 8
ngrid_rect <- nrows * ncols

grid_rect_length <- length / nrows
grid_rect_width <- width / ncols


### Algorithm to see proportion of cells in each grid square, show as histogram
data <- bgCluster

data <- determine_grid_rect_nums(data = data,
                                 length = length, width = width,
                                 nrows = nrows, ncols = ncols)

cell_props <- find_grid_rect_props(data = data,
                                   nrows = nrows, ncols = ncols)
cell_dens <- find_grid_rect_densities(data = data,
                                      nrows = nrows, ncols = ncols,
                                      area = grid_rect_length * grid_rect_width)



library(ggplot2)

## Plot all cells
plot_cells(data = data, length = length, width = width, nrows = nrows, ncols = ncols)



## Plot only showing a specific cell type
cell_type <- "Tumour"
data_specific_cell <- data[data$Cell.Type == cell_type, ]

plot_cells(data = data_specific_cell, length = length, width = width, nrows = nrows, ncols = ncols)



## Other stuff
cell_type <- "Tumour"
summary(cell_props[[cell_type]][cell_props[[cell_type]] > 0.0])
hist(cell_props[[cell_type]][cell_props[[cell_type]] > 0.0])

summary(cell_dens[[cell_type]][cell_dens[[cell_type]] > 0.0])


## The cluser detection function

# rect_data <- list()
# rect_data[[length(rect_data) + 1]] <- c()

rect_data <- c()
cell_type <- "Tumour"

for (grid_rect in seq(ngrid_rect)) {
  
  if (cell_props[[cell_type]][grid_rect] > 0.0) {
    x_coord <- (grid_rect - 1)%%ncols * grid_rect_width
    y_coord <- floor((grid_rect - 1)/nrows) * grid_rect_length
    
    rect_data <- append(rect_data, check_grid_rect(data,
                                                   cell_type,
                                                   x_coord,
                                                   y_coord,
                                                   grid_rect_length,
                                                   grid_rect_width,
                                                   c()))
  }
}

## Plot rectangles
plot_rect(rect_data = rect_data,
          xmin = 0,
          xmax = width,
          ymin = 0,
          ymax = length)


# bottom left coord of rectangle: (x_coord, y_coord)
check_grid_rect <- function(data, cell_type, x_coord, y_coord, length, width, 
                            answer) {
  
  data_temp <- data$Cell.Type[data$Cell.X.Position >= x_coord &
                              data$Cell.X.Position < (x_coord + width) &
                              data$Cell.Y.Position >= y_coord &
                              data$Cell.Y.Position < (y_coord + length)]
  
  if (length(data_temp) == 0) {
    return (c())
  }
  
  cell_prop <- (sum(data_temp == cell_type)) / length(data_temp)
  
  if (length < 5 || width < 5) {
    return (c())
  }
  
  if (cell_prop > 0.85) {
    return (c(x_coord, y_coord, length, width))
  }
  
  else if (cell_prop > 0) {
    # Bottom Left
    answer <- append(answer, check_grid_rect(data,
                                             cell_type,
                                             x_coord,
                                             y_coord,
                                             length/2,
                                             width/2,
                                             answer))
    # Bottom Right
    answer <- append(answer, check_grid_rect(data,
                                             cell_type,
                                             x_coord + width/2,
                                             y_coord,
                                             length/2,
                                             width/2,
                                             answer))
    # Top Left
    answer <- append(answer, check_grid_rect(data,
                                             cell_type,
                                             x_coord,
                                             y_coord + length/2,
                                             length/2,
                                             width/2,
                                             answer))
    
    # Top Right
    answer <- append(answer, check_grid_rect(data,
                                             cell_type,
                                             x_coord + width/2,
                                             y_coord + length/2,
                                             length/2,
                                             width/2,
                                             answer))
    
    return (answer)
  }
  
  else {
    return (c())
  }
}




### Functions

determine_grid_rect_nums <- function(data,
                                     length, width,
                                     nrows, ncols) {
  
  
  data$Cell.Num <- ncols * floor(data$Cell.Y.Position / (length / nrows)) +
                   floor(data$Cell.X.Position / (width / ncols)) + 1
  
  return(data)
}


find_grid_rect_props <- function(data, 
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


find_grid_rect_densities <- function(data, 
                                     nrows, ncols,
                                     area) {
  
  cell_types <- unique(data$Cell.Type)
  cell_dens <- vector(mode = 'list', length = length(cell_types))
  names(cell_dens) <- cell_types
  
  for (i in seq(nrows * ncols)) {
    cells <- data$Cell.Type[data$Cell.Num == i]
    
    for (cell_type in cell_types) {
      cell_dens[[cell_type]] <- append(cell_dens[[cell_type]], sum(cells == cell_type) / area) 
      
    }
  }  
  
  return (cell_dens)
}


plot_cells <- function(data, length, width, nrows, ncols) {
  
  ggplot(data = data,
         aes(x = Cell.X.Position, y = Cell.Y.Position, color = Cell.Type)) + 
    geom_point() + 
    theme(
      panel.background = element_rect(fill = NA),
      panel.ontop = TRUE,
      panel.grid.major = element_line(color = "black"),
      panel.grid.minor = element_blank()
    ) + 
    scale_x_continuous(limits = c(0, width), breaks = round(seq(0, width, width/ncols))) +
    scale_y_continuous(limits = c(0, length), breaks = round(seq(0, length, length/nrows))) 
}


plot_rect <- function(rect_data, xmin, xmax, ymin, ymax) {
  
  df <- data.frame(x = rect_data[seq(1, length(rect_data), 4)],
                   y = rect_data[seq(2, length(rect_data), 4)],
                   length = rect_data[seq(3, length(rect_data), 4)],
                   width = rect_data[seq(4, length(rect_data), 4)])
  
  ggplot(df, aes(x = x, y = y)) +
    geom_rect(aes(xmin = x, ymin = y, 
                  xmax = x + width, ymax = y + length),
              fill = "green",
              color = "black") +
    xlim(xmin, xmax) +
    ylim(ymin, ymax)
  
}

