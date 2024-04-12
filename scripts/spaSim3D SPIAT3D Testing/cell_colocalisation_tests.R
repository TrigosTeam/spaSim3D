## Get background cells
setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/Code/spaSim 3D/objects")
bg <- readRDS(file="bg.Rda")


## Nearing Immune Cluster 1 ------------------------------------------------
bg_cluster1 <- simulate_clusters3D(bg,
                                   n_clusters = 2,
                                   cluster_properties = list(
                                     C1 = list(
                                       name_of_cluster_cell = "Immune",
                                       infiltration_types = NULL,
                                       infiltration_proportions = NULL,
                                       shape = "Sphere",
                                       radius = 12,
                                       centre_loc = c(85, 85, 85)
                                     ),
                                     C2 = list(
                                       name_of_cluster_cell = "Tumour",
                                       infiltration_types = NULL,
                                       infiltration_proportions = NULL,
                                       shape = "Ellipsoid",
                                       x_radius = 15,
                                       y_radius = 20,
                                       z_radius = 25,
                                       centre_loc = c(50, 50, 50),
                                       x_y_rotation = 0,
                                       x_z_rotation = 0,
                                       y_z_rotation = 0
                                     )
                                   ),
                                   plot_image = F)
bg_cluster1$Cell.ID <- (paste("Cell_", seq(nrow(bg_cluster1)), sep=""))

## Nearing Immune Cluster 2 --------------------------------------------------
bg_cluster2 <- simulate_clusters3D(bg,
                                   n_clusters = 2,
                                   cluster_properties = list(
                                     C1 = list(
                                       name_of_cluster_cell = "Immune",
                                       infiltration_types = NULL,
                                       infiltration_proportions = NULL,
                                       shape = "Sphere",
                                       radius = 12,
                                       centre_loc = c(75, 75, 75)
                                     ),
                                     C2 = list(
                                       name_of_cluster_cell = "Tumour",
                                       infiltration_types = NULL,
                                       infiltration_proportions = NULL,
                                       shape = "Ellipsoid",
                                       x_radius = 15,
                                       y_radius = 20,
                                       z_radius = 25,
                                       centre_loc = c(50, 50, 50),
                                       x_y_rotation = 0,
                                       x_z_rotation = 0,
                                       y_z_rotation = 0
                                     )
                                   ),
                                   plot_image = F)
bg_cluster2$Cell.ID <- (paste("Cell_", seq(nrow(bg_cluster2)), sep=""))

## Nearing Immune Cluster 3 --------------------------------------------------
bg_cluster3 <- simulate_clusters3D(bg,
                                   n_clusters = 2,
                                   cluster_properties = list(
                                     C1 = list(
                                       name_of_cluster_cell = "Immune",
                                       infiltration_types = NULL,
                                       infiltration_proportions = NULL,
                                       shape = "Sphere",
                                       radius = 12,
                                       centre_loc = c(65, 65, 65)
                                     ),
                                     C2 = list(
                                       name_of_cluster_cell = "Tumour",
                                       infiltration_types = NULL,
                                       infiltration_proportions = NULL,
                                       shape = "Ellipsoid",
                                       x_radius = 15,
                                       y_radius = 20,
                                       z_radius = 25,
                                       centre_loc = c(50, 50, 50),
                                       x_y_rotation = 0,
                                       x_z_rotation = 0,
                                       y_z_rotation = 0
                                     )
                                   ),
                                   plot_image = F)
bg_cluster3$Cell.ID <- (paste("Cell_", seq(nrow(bg_cluster3)), sep=""))

## Immune Ring ---------------------------------------------------------------
bg_cluster4 <- simulate_rings3D(bg,
                                n_ring = 1,
                                ring_properties = list(
                                  R1 = list(
                                    name_of_cluster_cell = "Tumour",
                                    infiltration_types = NULL,
                                    infiltration_proportions = NULL,
                                    shape = "Ellipsoid",
                                    x_radius = 15,
                                    y_radius = 20,
                                    z_radius = 25,
                                    centre_loc = c(50, 50, 50),
                                    x_y_rotation = 0,
                                    x_z_rotation = 0,
                                    y_z_rotation = 0,
                                    name_of_ring_cell = "Immune",
                                    ring_width = 5,
                                    ring_infiltration_types = NULL,
                                    ring_infiltration_proportions = NULL
                                  )
                                ),
                                plot_image = F)
bg_cluster4$Cell.ID <- (paste("Cell_", seq(nrow(bg_cluster4)), sep=""))

## Increasing Immune Infiltration 1 ------------------------------------------
bg_cluster5 <- simulate_clusters3D(bg,
                                   n_clusters = 1,
                                   cluster_properties = list(
                                     C1 = list(
                                       name_of_cluster_cell = "Tumour",
                                       infiltration_types = "Immune",
                                       infiltration_proportions = 0.2,
                                       shape = "Ellipsoid",
                                       x_radius = 15,
                                       y_radius = 20,
                                       z_radius = 25,
                                       centre_loc = c(50, 50, 50),
                                       x_y_rotation = 0,
                                       x_z_rotation = 0,
                                       y_z_rotation = 0
                                     )
                                   ),
                                   plot_image = F)
bg_cluster5$Cell.ID <- (paste("Cell_", seq(nrow(bg_cluster5)), sep=""))

## Increasing Immune Infiltration 2 ------------------------------------------
bg_cluster6 <- simulate_clusters3D(bg,
                                   n_clusters = 1,
                                   cluster_properties = list(
                                     C1 = list(
                                       name_of_cluster_cell = "Tumour",
                                       infiltration_types = "Immune",
                                       infiltration_proportions = 0.4,
                                       shape = "Ellipsoid",
                                       x_radius = 15,
                                       y_radius = 20,
                                       z_radius = 25,
                                       centre_loc = c(50, 50, 50),
                                       x_y_rotation = 0,
                                       x_z_rotation = 0,
                                       y_z_rotation = 0
                                     )
                                   ),
                                   plot_image = F)
bg_cluster6$Cell.ID <- (paste("Cell_", seq(nrow(bg_cluster6)), sep=""))

## Increasing Immune Infiltration 3 ------------------------------------------
bg_cluster7 <- simulate_clusters3D(bg,
                                   n_clusters = 1,
                                   cluster_properties = list(
                                     C1 = list(
                                       name_of_cluster_cell = "Tumour",
                                       infiltration_types = "Immune",
                                       infiltration_proportions = 0.6,
                                       shape = "Ellipsoid",
                                       x_radius = 15,
                                       y_radius = 20,
                                       z_radius = 25,
                                       centre_loc = c(50, 50, 50),
                                       x_y_rotation = 0,
                                       x_z_rotation = 0,
                                       y_z_rotation = 0
                                     )
                                   ),
                                   plot_image = F)
bg_cluster7$Cell.ID <- (paste("Cell_", seq(nrow(bg_cluster7)), sep=""))



## All Data together ---------------------------------------------------------
# bg_clusters <- list(bg_cluster1,
#                     bg_cluster2,
#                     bg_cluster3,
#                     bg_cluster4,
#                     bg_cluster5,
#                     bg_cluster6,
#                     bg_cluster7)


setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/Code/spaSim 3D/objects")
# saveRDS(bg_clusters, file = "bg_clusters.rda")
bg_clusters <- readRDS(file="bg_clusters.rda")


### Using SPIAT3D for basic metrics ------------------------------------------
metrics_list <- vector(mode = "list", length = 7)
names(metrics_list) <- c("APD", "AMD", "MS", "NMS", "CKI", "CKAUC", "CIN")

## Assume Tumour is reference, Immune is target
reference_cell_type <- "Tumour"
target_cell_type <- "Immune"

## Test APD and AMD as they don't required radius input
for (cluster_data in bg_clusters) {
  
  ## Get pairwise distance data
  APD_data <- calculate_pairwise_distances_between_cell_types3D(cluster_data,
                                                                cell_types_of_interest = c(reference_cell_type, target_cell_type))
  metrics_list$APD <- append(metrics_list$APD, mean(APD_data[APD_data$Pair %in% c("Tumour/Immune", "Immune/Tumour"), "Distance"]))
  
  
  ## Get minimum distance data
  AMD_data <- calculate_minimum_distances_between_cell_types3D(cluster_data,
                                                               cell_types_of_interest = c(reference_cell_type, target_cell_type))
  metrics_list$AMD <- append(metrics_list$AMD, mean(AMD_data[AMD_data$Pair == "Tumour/Immune", "Distance"]))
  
}


## Test the remaining metrics
radius <- 30

for (cluster_data in bg_clusters) {
  
  ## Get mixing score data
  MS_data <- calculate_mixing_scores3D(cluster_data,
                                       reference_cell_types = reference_cell_type,
                                       target_cell_types = target_cell_type,
                                       radius = radius)
  
  metrics_list$MS <- append(metrics_list$MS, MS_data[["Mixing_score"]])
  metrics_list$NMS <- append(metrics_list$NMS, MS_data[["Normalised_mixing_score"]])
  
  ## Get cross K data
  cross_k_data <- calculate_Kcross3D(cluster_data,
                                     reference_cell_type = reference_cell_type,
                                     target_cell_type = target_cell_type,
                                     distance = radius,
                                     plot_results = FALSE)
  
  CKI <- calculate_Kcross_intersection3D(cross_k_data)
  if (length(CKI) == 2) {
    CKI <- CKI$Distance[1] ## Choose the first distance when crossing occurs
  }
  
  CKAUC <- calculate_AUC_of_Kcross3D(cross_k_data)
  
  metrics_list$CKI <- append(metrics_list$CKI, CKI)
  metrics_list$CKAUC <- append(metrics_list$CKAUC, CKAUC)
  
  ## Get cells in neighborhood data
  CIN_data <- calculate_cells_in_neighborhood3D(cluster_data,
                                                reference_cell_types = reference_cell_type,
                                                target_cell_types = target_cell_type,
                                                radius = radius)
  
  metrics_list$CIN <- append(metrics_list$CIN, mean(CIN_data[[reference_cell_type]][[target_cell_type]]))
}

metrics_list$MS <- as.numeric(metrics_list$MS)
metrics_list$NMS <- as.numeric(metrics_list$NMS)


### Plotting 
library(ggplot2)

result <- reshape2::melt(metrics_list)
colnames(result) <- c('Value', 'Metric')
result$Cluster_number <- rep(1:length(metrics_list[[1]]), length(metrics_list))
result$Value <- as.numeric(result$Value)
result$Metric <- factor(result$Metric, levels = names(metrics_list))

ggplot(result, aes(Cluster_number, Value)) +
  geom_line(color = "red") +
  geom_point() +
  facet_wrap(~ Metric, nrow = 7, scales = "free_y") +
  scale_x_continuous("Cluster Number", breaks = 1:7)


### Using SPIAT3D for gradient metrics ---------------------------------------

gradient_metrics_list <- vector(mode = "list", length = 3)
names(gradient_metrics_list) <- c("MSG", "CKG", "EGA")

radii <- 50
gradient_metrics_list[["MSG"]] <- data.frame(matrix(nrow = 0, ncol = radii))


for (cluster_data in bg_clusters) {
 
  ## Get mixing score gradient
  MSG_data <- calculate_mixing_scores_gradient3D(cluster_data,
                                                 reference_cell_type = "Immune",
                                                 target_cell_type = "Tumour",
                                                 radii = radii,
                                                 plot_image = F)
  MSG_data <- MSG_data$Normalised_mixing_score
  gradient_metrics_list[["MSG"]] <- rbind(gradient_metrics_list[["MSG"]], MSG_data)
  
}

colnames(gradient_metrics_list[["MSG"]]) <- seq(radii)



## Graph Mixing score gradient
df <- t(gradient_metrics_list[["MSG"]])
df <- reshape2::melt(df)
colnames(df) <- c("Radius", "Cluster_number", "Value")
df$Value <- as.numeric(df$Value)

ggplot(df, aes(Radius, Value)) +
  geom_line(color = "blue") +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  facet_wrap(~ Cluster_number, nrow = 7, scales = "free_y") +
  scale_x_continuous("Radius", breaks = 1:radii, labels = NULL)




gradient_metrics_list[["EGA"]] <- data.frame(matrix(nrow = 0, ncol = radii))


for (cluster_data in bg_clusters) {
  
  ## Get Entropy gradient aggregated
  EGA_data <- calculate_entropy_gradient_aggregated3D(cluster_data,
                                                      radii = radii,
                                                      reference_cell_type = "Tumour",
                                                      target_cell_types = c("Tumour", "Immune"), #Need both
                                                      plot_image = FALSE)
  EGA_data <- EGA_data$Entropy
  gradient_metrics_list[["EGA"]] <- rbind(gradient_metrics_list[["EGA"]], EGA_data)  
}

colnames(gradient_metrics_list[["EGA"]]) <- seq(radii)


## Graph Entropy gradient aggregated
df <- t(gradient_metrics_list[["EGA"]])
df <- reshape2::melt(df)
colnames(df) <- c("Radius", "Cluster_number", "Value")
df$Value <- as.numeric(df$Value)

ggplot(df, aes(Radius, Value)) +
  geom_line(color = "blue") +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  facet_wrap(~ Cluster_number, nrow = 7, scales = "free_y") +
  scale_x_continuous("Radius", breaks = 1:radii, labels = NULL)






gradient_metrics_list[["CKG"]] <- list()
radii <- 30
i <- 1
for (cluster_data in bg_clusters) {
  
  # Get Cross-K gradient
  CKG_data <- calculate_Kcross3D(cluster_data,
                                 reference_cell_type = reference_cell_type,
                                 target_cell_type = target_cell_type,
                                 distance = radii,
                                 plot_results = F)

  gradient_metrics_list[["CKG"]][[i]] <- CKG_data
  i <- i + 1
  
}



