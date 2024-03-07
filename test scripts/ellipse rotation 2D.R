length <- 20
width <- 20
n_points <- 1000

x_radius <- 15
y_radius <- 10
centre_loc <- c(5, 10)
theta <- 3 * pi / 16

df <- data.frame(x = runif(n_points, -length, length), 
                 y = runif(n_points, -width, width), 
                 color = "lightgray")

for (i in seq(n_points)) {
  x <- df$x[i] - centre_loc[1]
  y <- df$y[i] - centre_loc[2]
  
  D <- ((cos(theta) * x + sin(theta) * y)/x_radius)^2 + 
       ((-sin(theta) * x + cos(theta) * y)/y_radius)^2
  if (D <= 1) {
    df$color[i] <- "blue"
  }
}

df[n_points + 1, ] <- list(centre_loc[1],
                           centre_loc[2],
                           "red")

plot(df$x, df$y, pch = 19, col = df$color)

