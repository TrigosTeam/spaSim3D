length <- 20
width <- 20
n_points <- 1000

x_radius <- 10
y_radius <- 15
centre_loc <- c(0, 0)
theta <- 3 * pi / 8

df <- data.frame(x = runif(n_points, -length, length), 
                 y = runif(n_points, -width, width), 
                 color = "lightgray")

for (i in seq(n_points)) {
  D <- ((cos(theta) * df$x[i] + sin(theta) * df$y[i] - centre_loc[1])/x_radius)^2 + 
       ((-sin(theta) * df$x[i] + cos(theta) * df$y[i] - centre_loc[2])/y_radius)^2
  if (D <= 1) {
    df$color[i] <- "blue"
  }
}

plot(df$x, df$y, pch = 19, col = df$color)
