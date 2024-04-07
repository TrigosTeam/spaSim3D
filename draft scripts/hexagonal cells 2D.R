# Define xy
length <- 50
width  <- 100

# Define number of points
n.points <- 1000

# Obtain distance between each point using MAGIC formula
s <- sqrt((2 * length * width)/(sqrt(3) * n.points))

# Get value for x.points (points in 1 row) and y.points (points in 1 column), rounded up
x.points <- ceiling(length/s)
y.points <- ceiling((2 * width)/(sqrt(3) * s))


# First, assume points are on a rectangular grid
x <- rep(1:x.points, y.points) * s
y <- rep(1:y.points, each = x.points) * ((sqrt(3)*s)/2)

plot(x, y, pch = 19)

# Next, shift every second row by (s/2)
x <- x + c(rep(0, x.points), rep(s/2, x.points))

plot(x, y, pch = 19)
