library(dbscan)

dbscan_clustering3D <- function(spe,
                                cell_types_of_interest,
                                radius,
                                minimum_cells_in_radius,
                                feature_colname = "Cell.Type") {
  
  spe_subset <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  spe_subset_coords <- spatialCoords(spe_subset)
  
  db <- dbscan::dbscan(spe_subset_coords, eps = radius, minPts = minimum_cells_in_radius)
  
  
  
  ## Convert spe object to data frame
  # df <- data.frame(spatialCoords(spe), 
  #                  "Cell.Type" = spe[[feature_colname]],
  #                  "Cell.ID" = spe[["Cell.ID"]])
  # 
  # df_cell_types_of_interest <- df[df$Cell.Type %in% cell_types_of_interest, ]
  # df_other_cell_types <- df[!(df$Cell.Type %in% cell_types_of_interest), ]
  # df_cell_types_of_interest$alpha_hull_number <- alpha_hull_numbers
  # df_other_cell_types$alpha_hull_number <- -1
  
}



cell_cords <- formatted_data[,c("Cell.X.Position", "Cell.Y.Position")]
#Use dbscan to generate clusters
db <- dbscan::dbscan(cell_cords, eps = radius, minPts = min_neighborhood_size)
#since dbscan outputs cluster 0 as noise, we add 1 to all cluster numbers to keep it consistent
formatted_data$Cluster <- factor(db$cluster + 1)




spe2 <- spe1[ , spe1$Cell.Type %in% c("Tumour", "Immune")]
coords <- spatialCoords(spe2)
kNNdistplot(coords, minPts = 25)
abline(h = 13, col = "red", lty = 2) # Noise starts around 13
db <- dbscan::dbscan(coords, eps = 13, minPts = 25, borderPoints = T)
db



### Example ------------------------------------------------------------------
## use dbscan on the iris data set
data(iris)
iris <- as.matrix(iris[, 1:4])
## Find suitable DBSCAN parameters:
## 1. We use minPts = dim + 1 = 5 for iris. A larger value can also be used.
## 2. We inspect the k-NN distance plot for k = minPts - 1 = 4
kNNdistplot(iris, minPts = 5)
## Noise seems to start around a 4-NN distance of .7
abline(h=.7, col = "red", lty = 2)
## Cluster with the chosen parameters
res <- dbscan(iris, eps = .7, minPts = 5)
res
pairs(iris, col = res$cluster + 1L)
## Use a precomputed frNN object
fr <- frNN(iris, eps = .7)
dbscan(fr, minPts = 5)
