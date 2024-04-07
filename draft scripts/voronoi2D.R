library(deldir)
library(ggplot2)
library(ggvoronoi)


# Sample points
npoints <- 1000
length <- 200
width <- 200

x = runif(npoints, min = 0, max = length)
y = runif(npoints, min = 0, max = width)
points <- data.frame(x = x, 
                     y = y,
                     distance = sqrt((x - length/2)^2 + (y - width/2)^2))

circle <- data.frame(x = (length/2)*(1+cos(seq(0, 2*pi, length.out = 2500))),
                     y = (width/2)*(1+sin(seq(0, 2*pi, length.out = 2500))))

# Graph points to observe
ggplot(data = points, mapping = aes(x = x, y = y)) + geom_point()

# Plot the Voronoi diagram
ggplot(points) +
  geom_voronoi(aes(x = x, y = y, fill = distance), color = "black")

ggplot(points) +
  geom_voronoi(aes(x = x, y = y), fill = NA, color = "black")

ggplot(points) +
  geom_voronoi(aes(x = x, y = y, fill = distance), color = "black", outline = circle)
