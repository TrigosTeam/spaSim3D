library(plotly)
library(rgl)

### Colocalisation metrics

# APD
# Ref cells: (2, 1, 3), (3, 3, 1)
# Tar cells: (3, 1, 1), (1, 2, 4)

data <- data.frame(x = c(2, 3, 3, 1),
                   y = c(1, 3, 1, 2),
                   z = c(3, 1, 1, 4),
                   type = c("ref", "ref", "tar", "tar"),
                   line1 = c("y", "n", "y", "n"),
                   line2 = c("y", "n", "n", "y"))
data_ref <- data[data$type == "ref", ]
data_tar <- data[data$type == "tar", ]

fig <- plot_ly(type = 'scatter3d')

fig <- fig %>% layout(
  scene = list(aspectmode = 'cube',
               xaxis = list(title = 'x', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               yaxis = list(title = 'y', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               zaxis = list(title = 'z', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F))
  )

fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line1, mode = "lines", line = list(color = "gray", width = 5))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line2, mode = "lines", line = list(color = "gray", width = 5))
fig <- fig %>% add_markers(x = data_ref$x, y= data_ref$y, z = data_ref$z, marker = list(color = "red"), name = "ref")
fig <- fig %>% add_markers(x = data_tar$x, y= data_tar$y, z = data_tar$z, marker = list(color = "blue"), name = "tar")

fig



# AMD
# Ref cells: (2, 1, 3), (3, 3, 1), (1, 4, 2)
# Tar cells: (3, 1, 1), (1, 2, 4)
data <- data.frame(x = c(2, 3, 3, 1, 1),
                   y = c(1, 3, 1, 2, 4),
                   z = c(3, 1, 1, 4, 2),
                   type = c("ref", "ref", "tar", "tar", "ref"),
                   line1 = c("y", "n1", "n2", "y", "n3"),
                   line2 = c("n1", "y", "y", "n2", "n3"),
                   line3 = c("n1", "n2", "n3", "y", "y"))
data_ref <- data[data$type == "ref", ]
data_tar <- data[data$type == "tar", ]

fig <- plot_ly(type = 'scatter3d')

fig <- fig %>% layout(
  scene = list(aspectmode = 'cube',
               xaxis = list(title = 'x', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               yaxis = list(title = 'y', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               zaxis = list(title = 'z', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F))
)

fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line1, mode = "lines", line = list(color = "gray", width = 5))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line2, mode = "lines", line = list(color = "gray", width = 5))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line3, mode = "lines", line = list(color = "gray", width = 5))
fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z, marker = list(color = "red"), name = "ref")
fig <- fig %>% add_markers(x = data_tar$x, y = data_tar$y, z = data_tar$z, marker = list(color = "blue"), name = "tar")

fig


# MS & NMS
# Ref cells: (2, 1, 3), (3, 3, 2)
# Tar cells: (3, 3.5, 1), (0, 0, 4), (1, 1, 2.5)
data <- data.frame(x = c(2, 3, 3, 0, 1),
                   y = c(1, 3, 3.5, 0, 1),
                   z = c(3, 2, 1, 4, 2.5),
                   type = c("ref", "ref", "tar", "tar", "tar"),
                   line1 = c("y", "y", "n1", "n2", "n3"),
                   line2 = c("y", "n1", "n2", "n3", "y"),
                   line3 = c("y", "n1", "n2", "y", "n3"),
                   line4 = c("n1", "y", "n2", "n3", "y"),
                   line5 = c("n1", "y", "y", "n2", "n3"))
data_ref <- data[data$type == "ref", ]
data_tar <- data[data$type == "tar", ]

fig <- plot_ly(type = 'scatter3d')

fig <- fig %>% layout(
  scene = list(aspectmode = 'cube',
               xaxis = list(title = 'x', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               yaxis = list(title = 'y', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               zaxis = list(title = 'z', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F))
)

fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line1, mode = "lines", line = list(color = "tomato", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line2, mode = "lines", line = list(color = "skyblue", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line3, mode = "lines", line = list(color = "skyblue", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line4, mode = "lines", line = list(color = "skyblue", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line5, mode = "lines", line = list(color = "skyblue", width = 10))

fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z, marker = list(color = "red"), name = "ref")
fig <- fig %>% add_markers(x = data_tar$x, y = data_tar$y, z = data_tar$z, marker = list(color = "blue"), name = "tar")

fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z, 
                           marker = list(color = "lightgray", size = 145, opacity = 0.5))

fig


# CIN
# Ref cells: (2, 1, 3), (3, 3, 2)
# Tar cells blue: (3, 3.5, 1), (0, 0, 4), (1, 1, 2.5)
# Tar cells orange: (1.5, 2, 2), (3, 2, 2.5)
data <- data.frame(x = c(2, 3, 3, 0, 1, 1.5, 3.25),
                   y = c(1, 3, 3.5, 0, 1, 2, 3),
                   z = c(3, 2, 1, 4, 2.5, 2, 1.5),
                   type = c("ref", "ref", "tar1", "tar1", "tar1", "tar2", "tar2"),
                   line1 = c("y", "y", "n1", "n2", "n3", "n4", "n5"),
                   line2 = c("y", "n1", "n2", "n3", "y", "n4", "n5"),
                   line3 = c("y", "n1", "n2", "y", "n3", "n4", "n5"),
                   line4 = c("n1", "y", "n2", "n3", "y", "n4", "n5"),
                   line5 = c("n1", "y", "y", "n2", "n3", "n4", "n5"),
                   line6 = c("y", "n1", "n2", "n3", "n4", "y", "n5"),
                   line7 = c("n1", "y", "n2", "n3", "n4", "y", "n5"),
                   line8 = c("n1", "y", "n2", "n3", "n4", "n5", "y"))
data_ref <- data[data$type == "ref", ]
data_tar1 <- data[data$type == "tar1", ]
data_tar2 <- data[data$type == "tar2", ]

fig <- plot_ly(type = 'scatter3d')

fig <- fig %>% layout(
  scene = list(aspectmode = 'cube',
               xaxis = list(title = 'x', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               yaxis = list(title = 'y', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               zaxis = list(title = 'z', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F))
)

fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line2, mode = "lines", line = list(color = "skyblue", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line3, mode = "lines", line = list(color = "skyblue", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line4, mode = "lines", line = list(color = "skyblue", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line5, mode = "lines", line = list(color = "skyblue", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line5, mode = "lines", line = list(color = "skyblue", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line5, mode = "lines", line = list(color = "skyblue", width = 10))

fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line6, mode = "lines", line = list(color = "orchid", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line7, mode = "lines", line = list(color = "orchid", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line8, mode = "lines", line = list(color = "orchid", width = 10))

fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z, marker = list(color = "red"), name = "ref")
fig <- fig %>% add_markers(x = data_tar1$x, y = data_tar1$y, z = data_tar1$z, marker = list(color = "blue"), name = "tar1")
fig <- fig %>% add_markers(x = data_tar2$x, y = data_tar2$y, z = data_tar2$z, marker = list(color = "purple"), name = "tar2")

fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z, 
                           marker = list(color = "lightgray", size = 145, opacity = 0.5))

fig



# Entropy 1 (low entropy)
# Ref cells: (4, 4, 4)
# Tar cells: (3, 3.5, 1), (2, 2, 3), (1, 1, 2.5)
data <- data.frame(x = c(4, 3, 2, 1),
                   y = c(4, 3.5, 3.5, 1),
                   z = c(4, 1, 3, 2.5),
                   type = c("ref", "tar", "tar", "tar"),
                   line1 = c("y", "y", "n1", "n2"),
                   line2 = c("y", "n1", "y", "n2"),
                   line3 = c("y", "n1", "n2", "y"))
data_ref <- data[data$type == "ref", ]
data_tar <- data[data$type == "tar", ]

fig <- plot_ly(type = 'scatter3d')

fig <- fig %>% layout(
  scene = list(aspectmode = 'cube',
               xaxis = list(title = 'x', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               yaxis = list(title = 'y', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               zaxis = list(title = 'z', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F))
)

fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line1, mode = "lines", line = list(color = "skyblue", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line2, mode = "lines", line = list(color = "skyblue", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line3, mode = "lines", line = list(color = "skyblue", width = 10))


fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z, marker = list(color = "red"), name = "ref")
fig <- fig %>% add_markers(x = data_tar$x, y = data_tar$y, z = data_tar$z, marker = list(color = "blue"), name = "tar")

fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z, 
                           marker = list(color = "lightgray", size = 145, opacity = 0.5))

fig


# Entropy 2 (high entropy)
# Ref cells: (4, 4, 4)
# Tar1 cells: (3, 3.5, 1), (2, 2, 3)
# Tar 2 cells: (1, 1, 2.5), (1, 3, 3)
data <- data.frame(x = c(4, 3, 2, 1, 1),
                   y = c(4, 3.5, 3.5, 1, 3),
                   z = c(4, 1, 3, 2.5, 3),
                   type = c("ref", "tar1", "tar1", "tar2", "tar2"),
                   line1 = c("y", "y", "n1", "n2", "n3"),
                   line2 = c("y", "n1", "y", "n2", "n3"),
                   line3 = c("y", "n1", "n2", "y", "n3"),
                   line4 = c("y", "n1", "n2", "n3", "y"))
data_ref <- data[data$type == "ref", ]
data_tar1 <- data[data$type == "tar1", ]
data_tar2 <- data[data$type == "tar2", ]

fig <- plot_ly(type = 'scatter3d')

fig <- fig %>% layout(
  scene = list(aspectmode = 'cube',
               xaxis = list(title = 'x', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               yaxis = list(title = 'y', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               zaxis = list(title = 'z', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F))
)

fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line1, mode = "lines", line = list(color = "skyblue", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line2, mode = "lines", line = list(color = "skyblue", width = 10))

fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line3, mode = "lines", line = list(color = "orchid", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line4, mode = "lines", line = list(color = "orchid", width = 10))

fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z, marker = list(color = "red"), name = "ref")
fig <- fig %>% add_markers(x = data_tar1$x, y = data_tar1$y, z = data_tar1$z, marker = list(color = "blue"), name = "tar1")
fig <- fig %>% add_markers(x = data_tar2$x, y = data_tar2$y, z = data_tar2$z, marker = list(color = "purple"), name = "tar2")

fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z, 
                           marker = list(color = "lightgray", size = 145, opacity = 0.5))

fig


# Cross-K function 1 (high clustering)
# Ref cells: (4, 4, 4)
# Tar cells: (3, 3.5, 1), (2, 2.5, 3), (1, 3.5, 2.5) (3.5, 1, 3)
data <- data.frame(x = c(4, 3, 2, 1, 3.5),
                   y = c(4, 3.5, 3.5, 1, 1),
                   z = c(4, 1, 3, 2.5, 3),
                   type = c("ref", "tar", "tar", "tar", "tar"),
                   line1 = c("y", "y", "n1", "n2", "n3"),
                   line2 = c("y", "n1", "y", "n2", "n3"),
                   line3 = c("y", "n1", "n2", "y", "n3"),
                   line4 = c("y", "n1", "n2", "n3", "y"))
data_ref <- data[data$type == "ref", ]
data_tar <- data[data$type == "tar", ]

fig <- plot_ly(type = 'scatter3d')

fig <- fig %>% layout(
  scene = list(aspectmode = 'cube',
               xaxis = list(title = 'x', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               yaxis = list(title = 'y', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               zaxis = list(title = 'z', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F))
)

fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line1, mode = "lines", line = list(color = "skyblue", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line2, mode = "lines", line = list(color = "skyblue", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line3, mode = "lines", line = list(color = "skyblue", width = 10))
fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line4, mode = "lines", line = list(color = "skyblue", width = 10))


fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z, marker = list(color = "red"), name = "ref")
fig <- fig %>% add_markers(x = data_tar$x, y = data_tar$y, z = data_tar$z, marker = list(color = "blue"), name = "tar")

fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z, 
                           marker = list(color = "lightgray", size = 145, opacity = 0.5))

fig




# Cross-K function 1 (low clustering)
# Ref cells: (4, 4, 4)
# Tar cells: (3, 3.5, 1), (1, 5, 4), (1, 1, 5) (4.5, 1, 1)
data <- data.frame(x = c(4, 3, 1, 1, 4.5),
                   y = c(4, 3.5, 5, 1, 1),
                   z = c(4, 1, 4, 5, 1),
                   type = c("ref", "tar", "tar", "tar", "tar"),
                   line1 = c("y", "y", "n1", "n2", "n3"),
                   line2 = c("y", "n1", "y", "n2", "n3"),
                   line3 = c("y", "n1", "n2", "y", "n3"),
                   line4 = c("y", "n1", "n2", "n3", "y"))
data_ref <- data[data$type == "ref", ]
data_tar <- data[data$type == "tar", ]

fig <- plot_ly(type = 'scatter3d')

fig <- fig %>% layout(
  scene = list(aspectmode = 'cube',
               xaxis = list(title = 'x', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               yaxis = list(title = 'y', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               zaxis = list(title = 'z', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F))
)

fig <- fig %>% add_trace(x = data$x, y= data$y, z = data$z, color = data$line1, mode = "lines", line = list(color = "skyblue", width = 10))


fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z, marker = list(color = "red"), name = "ref")
fig <- fig %>% add_markers(x = data_tar$x, y = data_tar$y, z = data_tar$z, marker = list(color = "blue"), name = "tar")

fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z, 
                           marker = list(color = "lightgray", size = 145, opacity = 0.5))

fig




# gradient (increasing radii)
data_ref <- data.frame(x = c(1, 3, 4),
                       y = c(4, 2, 3),
                       z = c(3, 4 ,2))

data_tar <- data.frame(x = c(1, 2, 3, 2, 2, 3, 1, 4),
                       y = c(2, 4, 3, 1, 2, 3, 3, 3),
                       z = c(4, 1, 2, 3, 3, 2, 2, 1))


fig <- plot_ly(type = 'scatter3d')

fig <- fig %>% layout(
  scene = list(aspectmode = 'cube',
               xaxis = list(title = 'x', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               yaxis = list(title = 'y', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F),
               zaxis = list(title = 'z', range = c(0, 5), showgrid = T, showaxeslabels = T, showticklabels = F))
)

fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z, marker = list(color = "red"), name = "ref")
fig <- fig %>% add_markers(x = data_tar$x, y = data_tar$y, z = data_tar$z, marker = list(color = "blue"), name = "tar")

fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z,
                           marker = list(color = "lightgray", size = 50, opacity = 0.3))
fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z,
                           marker = list(color = "lightgray", size = 100, opacity = 0.3))
fig <- fig %>% add_markers(x = data_ref$x, y = data_ref$y, z = data_ref$z,
                           marker = list(color = "lightgray", size = 150, opacity = 0.3))

fig


