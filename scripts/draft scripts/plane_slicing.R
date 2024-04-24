length = 100
width = 100
height = 100
n_cells = 10000

library(plotly)

setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/Code/spaSim 3D/objects")
all_plots_data <- readRDS(file="all_plots_test_data.rda")
all_plots_meta_data <- readRDS(file="all_plots_meta_data.rda")

data <- all_plots_data[[7]]

s <- 150 # Size of window
delta <- -15 # Position of slice
thickness <- 0.025

data[0.5 * data$Cell.X.Position + 0.5 * data$Cell.Y.Position - data$Cell.Z.Position > -delta - thickness * s &
     0.5 * data$Cell.X.Position + 0.5 * data$Cell.Y.Position - data$Cell.Z.Position < -delta + thickness * s,
     "Cell.Type"] <- "Slice"

# Get number of cells in the slice
sum(data$Cell.Type == "Slice")


plot_cell_categories3D(data,
                       c("Others", "Tumour", "Immune", "Slice"),
                       c("lightgray", "orange", "skyblue", "tomato"))


