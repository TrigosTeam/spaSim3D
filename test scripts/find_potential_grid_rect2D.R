library(spaSim)

### Simulate 2D data ----------------------------------------------------------
bgHardcore <- simulate_background_cells(n_cells = 5000,
                                        width = 2000,
                                        height = 2000,
                                        method = "Hardcore",
                                        min_d = 10,
                                        #oversampling_rate,
                                        #jitter,
                                        #Cell.Type,
                                        plot_image = T)

cluster_prop <- list(
  C1 = list(
    name_of_cluster_cell = "Tumour",
    size = 500,
    shape = "Oval",
    centre_loc = data.frame("x" = 700, "y" = 800),
    infiltration_types = c("Immune1", "Others"),
    infiltration_proportions = c(0.10, 0.00)),
  C2 = list(
    name_of_cluster_cell = "Tumour",
    size = 400,
    shape = "Circle",
    centre_loc = data.frame("x" = 1500, "y" = 1500),
    infiltration_types = c("Immune1", "Others"),
    infiltration_proportions = c(0.5, 0.00)),
  C3 = list(
    name_of_cluster_cell = "Immune1",
    size = 250,
    shape = "Circle",
    centre_loc = data.frame("x" = 700, "y" = 700),
    infiltration_types = c("Tumour"),
    infiltration_proportions = c(0.4)
  )
)

bgCluster <- simulate_clusters(bg_sample = bgHardcore,
                               n_clusters = 3,
                               bg_type = "Others",
                               cluster_properties = cluster_prop,
                               plot_image = T,
                               plot_categories = c("Others", "Tumour", "Immune1"),
                               plot_colours = NULL)

data <- bgCluster

### Algorithm -----------------------------------------------------------------

## Dimensions of study area
length <- 2000
width <- 2000

## Define grid dimensions
nrows <- 8
ncols <- 8
ngrid_rect <- nrows * ncols

grid_rect_length <- length / nrows
grid_rect_width <- width / ncols

## Determine grid numbers for data
data <- determine_grid_rect_nums(data = data,
                                 length = length, width = width,
                                 nrows = nrows, ncols = ncols)


## Plot all cells
library(ggplot2)
plot_cells(data = data, length = length, width = width, nrows = nrows, ncols = ncols)


## Plot only showing a specific cell type
cell_type <- "Tumour"
data_specific_cell <- data[data$Cell.Type == cell_type, ]

plot_cells(data = data_specific_cell, length = length, width = width, nrows = nrows, ncols = ncols)





## Step 1. Create vector with number for each grid_rect
grid_rect_nums <- seq(ngrid_rect)


## Step 2. Determine cell proportions for each cell_type for each grid_rect
cell_props <- find_grid_rect_props(data = data,
                                   nrows = nrows, ncols = ncols)


## Step 3. Remove numbers from grid_rect_nums where cell_prop = 0
cell_type <- "Tumour"
grid_rect_nums <- grid_rect_nums[cell_props[[cell_type]] > 0]



## Step 4. Determine grid_rect with highest cell_prop
max_rect_num <- find_max_prop_grid_rect(cell_props = cell_props,
                                    grid_rect_nums = grid_rect_nums,
                                    cell_type = cell_type)

max_rect_prop <- cell_props[[cell_type]][max_rect_num]


## Step 6. Recursive algorithm. Determine if adjacent grid_rects have a
#          cell_prop > 25% of max_rect_prop. If they do, remove from
#          grid_rect_nums and append to cluster_rect_nums.

cluster_rect_nums <- grid_rect_nums

# REMOVE elements from grid_rect_nums which ARE included in the cluster
grid_rect_nums <- check_adjacent_grid_rects(max_rect_num, grid_rect_nums)

# Take the difference to get elements included in the cluster
cluster_rect_nums <- setdiff(cluster_rect_nums, grid_rect_nums)



check_adjacent_grid_rects <- function(curr_grid_rect_num, grid_rect_nums) {
  
  if (curr_grid_rect_num %in% grid_rect_nums == FALSE) {
    return (grid_rect_nums)
  }
  
  if (cell_props[[cell_type]][curr_grid_rect_num] > 0.25 * max_rect_prop) {
    
    grid_rect_nums <- grid_rect_nums[! grid_rect_nums %in% curr_grid_rect_num]
    
    ### CHECK LEFT, RIGHT, UP, DOWN
    
    ## Need to check if going left, right, up or down is possible
    
    # Left
    if (curr_grid_rect_num%%ncols != 1) {
      grid_rect_nums <- check_adjacent_grid_rects(curr_grid_rect_num - 1, 
                                                  grid_rect_nums)  
    }
    
    # Right
    if (curr_grid_rect_num%%ncols != 0) {
      grid_rect_nums <- check_adjacent_grid_rects(curr_grid_rect_num + 1, 
                                                  grid_rect_nums)  
    }
    
    # Up
    if (curr_grid_rect_num <= ncols * (nrows - 1)) {
      grid_rect_nums <- check_adjacent_grid_rects(curr_grid_rect_num + ncols, 
                                                         grid_rect_nums)
    }
    
    # Down
    if (curr_grid_rect_num > ncols) {
      grid_rect_nums <- check_adjacent_grid_rects(curr_grid_rect_num - ncols, 
                                                         grid_rect_nums)
    }
    
  }
  return (grid_rect_nums)
}

# x = (cluster_rect_nums - 1)%%ncols * grid_rect_width
# y = floor((cluster_rect_nums - 1)/ncols) * grid_rect_length
# rect_data = data.frame(x = x, y = y, length = grid_rect_length, width = grid_rect_width)
# plot_rect(rect_data, length, width, nrows, ncols)


## Functions
determine_grid_rect_nums <- function(data,
                                     length, width,
                                     nrows, ncols) {
  
  
  data$Cell.Num <- ncols * floor(data$Cell.Y.Position / (length / nrows)) +
    floor(data$Cell.X.Position / (width / ncols)) + 1
  
  return(data)
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


find_max_prop_grid_rect <- function(cell_props, grid_rect_nums, cell_type) {
  
  return (grid_rect_nums[which.max(cell_props[[cell_type]][grid_rect_nums])])
}
