length = 100
width = 100
height = 100
n_cells = 10000

library(plotly)

setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/Code/spaSim 3D/objects")
all_plots_data <- readRDS(file="all_plots_test_data.rda")
all_plots_meta_data <- readRDS(file="all_plots_meta_data.rda")

data <- all_plots_data[[33]]

s <- 150 # Size of window
delta <- -15 # Position of slice
thickness <- 5

# Plot in 3D
data3D <- data
data3D[0.5 * data$Cell.X.Position + 0.5 * data$Cell.Y.Position - data$Cell.Z.Position > -delta - thickness &
       0.5 * data$Cell.X.Position + 0.5 * data$Cell.Y.Position - data$Cell.Z.Position < -delta + thickness,
       "Cell.Type"] <- "Slice"
plot_cell_categories3D(data3D,
                       c("Others", "Tumour", "Immune", "Slice"),
                       c("lightgray", "orange", "skyblue", "tomato"))


## Plot in 2D
data2D <- data[0.5 * data$Cell.X.Position + 0.5 * data$Cell.Y.Position - data$Cell.Z.Position > -delta - thickness &
               0.5 * data$Cell.X.Position + 0.5 * data$Cell.Y.Position - data$Cell.Z.Position < -delta + thickness, ]

print(paste("Number of cells:", nrow(data2D)))

data2D$Cell.Type <- ordered(data2D$Cell.Type, levels = c("Others", "Tumour", "Immune"))
ggplot(data2D, aes(Cell.X.Position, Cell.Y.Position, color = Cell.Type)) +
  geom_point() +
  scale_colour_manual(values = c("lightgray", "orange", "skyblue"))
