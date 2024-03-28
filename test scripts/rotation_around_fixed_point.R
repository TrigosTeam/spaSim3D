## Define rectangle dimensions and centre
length = 40
width = 80
centre = c(200, 100) # x, y coords of centre

## Get random points inside this rectange
n_points <- 1000
x = runif(n_points, centre[1] - width/2, centre[1] + width/2)
y = runif(n_points, centre[2] - length/2, centre[2] + length/2)
df <- data.frame(x = x, y = y)


## Plot points
library(ggplot2)
ggplot(df, aes(x, y)) +
  geom_point(color = "red") +
  xlim(0, 300) +
  ylim(0, 300)


## Apply rotation to each point
theta = pi/4
x_new = cos(theta) * (x - centre[1]) - sin(theta) * (y - centre[2]) + centre[1]
y_new = sin(theta) * (x - centre[1]) + cos(theta) * (y - centre[2]) + centre[2]
df_new = data.frame(x = x_new, y = y_new)

## Plot rotated points
ggplot(df_new, aes(x, y)) +
  geom_point(color = "red") +
  xlim(0, 300) +
  ylim(0, 300)


## Apply a vertical and horizontal translation to each point
x_shift = 25
y_shift = 50
x_new1 = x_new + x_shift
y_new1 = y_new + y_shift
df_new1 = data.frame(x = x_new1, y = y_new1)

## Plot translated points
ggplot(df_new1, aes(x, y)) +
  geom_point(color = "red") +
  xlim(0, 300) +
  ylim(0, 300)



### --------------------------------------------------------------------------
### Slide alignment attempt
### --------------------------------------------------------------------------
## Define rectangle dimensions and centre
length <- 100
width <- 50
centre <- c(width/2, length/2)

## Get random points inside rectangle1
n_points <- 100
x = runif(n_points, centre[1] - width/2, centre[1] + width/2)
y = runif(n_points, centre[2] - length/2, centre[2] + length/2)
df1 <- data.frame(x = x, y = y)

## Add type
df1$type <- ifelse(df1$x > 25 & df1$y < 25, "typeA", "typeB")

## Plot points
library(ggplot2)
ggplot(df1, aes(x, y, color=type)) +
  geom_point() +
  xlim(-50, 150) +
  ylim(-50, 150)


## Get random points inside rectangle2
n_points <- 80
x = runif(n_points, centre[1] - width/2, centre[1] + width/2)
y = runif(n_points, centre[2] - length/2, centre[2] + length/2)
df2 <- data.frame(x = x, y = y)

## Add type
df2$type <- ifelse(df2$x > 30 & df2$y > 5 & df2$y < 30, "typeA", "typeB")

## Plot points
ggplot(df2, aes(x, y, color=type)) +
  geom_point() +
  xlim(-50, 150) +
  ylim(-50, 150)


###---------------------------------------------------------------------------
## Determine the translations and rotations required so that the min.distance
## between same cell types between the two rectangles is a minimum
###---------------------------------------------------------------------------

x_shifts <- seq(-5, 5, length.out = 10)
y_shifts <- seq(-5, 5, length.out = 10)
angle_shifts <- seq(-pi/6, pi/6, length.out = 10)

min_d = 100
for (x_shift in x_shifts) {
  for (y_shift in y_shifts) {
    for (angle_shift in angle_shifts) {
      x_new <- cos(angle_shift) * (df2$x - centre[1]) - 
               sin(angle_shift) * (df2$y - centre[2]) + centre[1] + x_shift
      
      y_new <- sin(angle_shift) * (df2$x - centre[1]) + 
               cos(angle_shift) * (df2$y - centre[2]) + centre[2] + y_shift
      
      df_new <- data.frame(x = x_new, y = y_new, type = df2$type)
      
      d <- average_minimum_distance(df1, df_new)
      
      if (d < min_d) {
        min_d <- d
        transformations <- c(x_shift, y_shift, angle_shift)
      }
    }
  }
}


## apply transformations to df2
x_shift <- transformations[1]
y_shift <- transformations[2]
angle_shift <- transformations[3]

x_new <- cos(angle_shift) * (df2$x - centre[1]) - 
  sin(angle_shift) * (df2$y - centre[2]) + centre[1] + x_shift

y_new <- sin(angle_shift) * (df2$x - centre[1]) + 
  cos(angle_shift) * (df2$y - centre[2]) + centre[2] + y_shift


df2_transformed <- data.frame(x = x_new, y = y_new, type = df2$type)

## Plot points
ggplot(df2_transformed, aes(x, y, color=type)) +
  geom_point() +
  xlim(-50, 150) +
  ylim(-50, 150)


## Plot df1 and df2_transformed together
df2_transformed[df2_transformed$type == "typeA", "type"] <- "typeA*" 
df2_transformed[df2_transformed$type == "typeB", "type"] <- "typeB*" 

df_combined <- data.frame(x = c(df1$x, df2_transformed$x),
                          y = c(df1$y, df2_transformed$y),
                          type = c(df1$type, df2_transformed$type))

ggplot(df_combined, aes(x, y, color=type)) +
  geom_point() +
  xlim(-50, 150) +
  ylim(-50, 150)




average_minimum_distance <- function(data1, data2) {
 
  ## Assume types are matching between data1 and data2 
  
  types <- unique(data1$type)
  sum <- 0
  
  for (type in types) {
    data1_temp <- data1[data1$type == type, ]
    data2_temp <- data2[data2$type == type, ]
    
    vec <- c()
    
    for (i in seq(nrow(data1_temp))) {
      min_d = 100
      point1 <- c(data1_temp[i, "x"], data1_temp[i, "y"])
      
      for (j in seq(nrow(data2_temp))) {
        point2 <- c(data2_temp[j, "x"], data2_temp[j, "y"])
        d <- get_distance(point1, point2)
        
        if (d < min_d) {
          min_d <- d
        }
        
      }
      vec <- append(vec, min_d)
    }
    sum <- sum + 
      (nrow(data1_temp) + nrow(data2_temp)) * mean(vec) / 
      (nrow(data1) + nrow(data2))
  }
  return (sum)
}


get_distance <- function(point1, point2) {
  
  return(sqrt((point1[1] - point2[1])^2 + (point1[2] - point2[2])^2))

}
