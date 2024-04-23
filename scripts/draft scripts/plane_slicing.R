length = 100
width = 100
height = 100
n_cells = 10000

data <- data.frame(Cell.X.Position = runif(n_cells, 0, length),
                   Cell.Y.Position = runif(n_cells, 0, width),
                   Cell.Z.Position = runif(n_cells, 0, height),
                   Cell.Type = "Others")



s <- 100
delta <- 0

data[0.5 * data$Cell.X.Position + 0.5 * data$Cell.Y.Position - data$Cell.Z.Position > -delta - 0.05 * s &
     0.5 * data$Cell.X.Position + 0.5 * data$Cell.Y.Position - data$Cell.Z.Position < -delta + 0.05 * s,
     "Cell.Type"] <- "Tumour"



# s <- 100
# delta <- 20
# d <- 2 * delta * s

# if (delta >= 0) {
#   data[s * data$Cell.X.Position + s * data$Cell.Y.Position + (-2 * s + 4 * delta) * data$Cell.Z.Position > 0.8 * d &
#          s * data$Cell.X.Position + s * data$Cell.Y.Position + (-2 * s + 4 * delta) * data$Cell.Z.Position < 1.2 * d,
#        "Cell.Type"] <- "Tumour"  
# } else if (delta < 0) {
#   data[s * data$Cell.X.Position + s * data$Cell.Y.Position + (-2 * s + 4 * delta) * data$Cell.Z.Position < 0.8 * d &
#          s * data$Cell.X.Position + s * data$Cell.Y.Position + (-2 * s + 4 * delta) * data$Cell.Z.Position > 1.2 * d,
#        "Cell.Type"] <- "Tumour"
# }


library(plotly)
plot_cell_categories3D(data, 
                       cell_types_of_interest = c("Others", "Tumour"), 
                       colour_vector = c("lightgray", "orange"))

