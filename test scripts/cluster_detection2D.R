## Dimensions of study area
length <- 2000
width <- 2000

## Define grid dimensions
nrows <- 8
ncols <- 8
ngrid_rect <- nrows * ncols

grid_rect_length <- length / nrows
grid_rect_width <- width / ncols

# Get data
data <- bgCluster # from spaSim script

# Determine the grid rectangle number each point belongs
# Bottom left is grid_rect_1, move right for grid_rect_2, 
# once at the end of the first row, go to the second row 
data <- determine_grid_rect_nums(data = data,
                                 length = length, width = width,
                                 nrows = nrows, ncols = ncols)

# Determine cell proportions for each cell type in each grid rectangle
cell_props <- find_grid_rect_props(data = data,
                                   nrows = nrows, ncols = ncols)

# Determine cell densities for each cell type in each grid rectangle
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
# hist(cell_props[[cell_type]][cell_props[[cell_type]] > 0.0])

summary(cell_dens[[cell_type]][cell_dens[[cell_type]] > 0.0])




# Determining the rectangles that represent the clusters
rect_data <- c()
cell_type <- "Tumour"

for (grid_rect in seq(ngrid_rect)) {
  
  # Only examine rectangles which contain the chosen cell_type
  if (cell_props[[cell_type]][grid_rect] > 0.0) {
    x_coord <- (grid_rect - 1)%%ncols * grid_rect_width
    y_coord <- floor((grid_rect - 1)/ncols) * grid_rect_length
    
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
          length = length,
          width = width,
          nrows = nrows,
          ncols = ncols)


# THE CLUSTER DETECTION FUNCTION
# bottom left coord of rectangle: (x_coord, y_coord)
check_grid_rect <- function(data, cell_type, x_coord, y_coord, length, width, 
                            answer) {
  
  # size of rectangle is getting too small
  if (length < 20 | width < 20) {
    return (c())
  }
  
  # look at cells only in the current rectangle
  data <- data[data$Cell.X.Position >= x_coord &
               data$Cell.X.Position < (x_coord + width) &
               data$Cell.Y.Position >= y_coord &
               data$Cell.Y.Position < (y_coord + length), ]
  
  # consider cell_types in the rectangle
  data_cell_type <- data$Cell.Type
  
  # no cells in the rectangle
  if (length(data_cell_type) == 0) {
    return (c())
  }

  # determine cell proportion for chosen cell_type in the rectangle
  cell_prop <- (sum(data_cell_type == cell_type)) / length(data_cell_type)
  
  # cell proportion in the rectangle must be more than the 75th percentile proportion
  if (cell_prop >= quantile(cell_props[[cell_type]][cell_props[[cell_type]] > 0.0], 0.75)) {
    return (c(x_coord, y_coord, length, width))
  }
  
  # some cell_type of interest still in the rectangle, check sub-rectangles
  else if (cell_prop > 0) {
    # Bottom Left
    answer <- append(answer, check_grid_rect(data,
                                             cell_type,
                                             x_coord,
                                             y_coord,
                                             length/2,
                                             width/2,
                                             c()))
    # Bottom Right
    answer <- append(answer, check_grid_rect(data,
                                             cell_type,
                                             x_coord + width/2,
                                             y_coord,
                                             length/2,
                                             width/2,
                                             c()))
    # Top Left
    answer <- append(answer, check_grid_rect(data,
                                             cell_type,
                                             x_coord,
                                             y_coord + length/2,
                                             length/2,
                                             width/2,
                                             c()))
    
    # Top Right
    answer <- append(answer, check_grid_rect(data,
                                             cell_type,
                                             x_coord + width/2,
                                             y_coord + length/2,
                                             length/2,
                                             width/2,
                                             c()))
    
    return (answer)
  }
  
  # cell proportion is zero
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
      
      answer <- ifelse(length(cells) > 0, 
                       sum(cells == cell_type) / length(cells), 0)
      
      cell_props[[cell_type]] <- append(cell_props[[cell_type]], answer) 
      
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
      answer <- ifelse(length(cells) > 0, 
                       sum(cells == cell_type) / area, 0)
      
      cell_dens[[cell_type]] <- append(cell_dens[[cell_type]], answer) 
      
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




plot_rect <- function(rect_data, length, width, nrows, ncols) {
  
  df <- data.frame(x = rect_data[seq(1, length(rect_data), 4)],
                   y = rect_data[seq(2, length(rect_data), 4)],
                   length = rect_data[seq(3, length(rect_data), 4)],
                   width = rect_data[seq(4, length(rect_data), 4)])
  
  ggplot(df, aes(x = x, y = y)) +
    geom_rect(aes(xmin = x, ymin = y, 
                  xmax = x + width, ymax = y + length),
              fill = "green",
              color = "black") +
    theme(
      panel.background = element_rect(fill = NA),
      panel.ontop = TRUE,
      panel.grid.major = element_line(color = "black"),
      panel.grid.minor = element_blank()
    ) + 
    scale_x_continuous(limits = c(0, width), breaks = round(seq(0, width, width/ncols))) +
    scale_y_continuous(limits = c(0, length), breaks = round(seq(0, length, length/nrows))) 
  
}



# x <- y <- c()
# density <- max(cell_dens[[cell_type]])
# 
# for (i in seq(nrow(df))) {
#   row <- df[i, ]
#   npoints <- round(density * row$length * row$width)
#   x <- append(x, runif(npoints, min = row$x, max = row$x + row$width))
#   y <- append(y, runif(npoints, min = row$y, max = row$y + row$length))
# }
# data_sim <- data.frame(Cell.X.Position = x, Cell.Y.Position = y, Cell.Type = cell_type)
# 
# plot_cells(data = data_sim, length = length, width = width, nrows = nrows, ncols = ncols)
