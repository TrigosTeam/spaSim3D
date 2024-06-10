library(SPIAT)
library(plotly)
library(ggplot2)

setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/Code/spaSim 3D/objects")
all_plots_data <- readRDS(file="all_plots_test_data.rda")
all_plots_meta_data <- readRDS(file="all_plots_meta_data.rda")


data_index <- 279
## Plot 3D data normally
plot_cell_categories3D(all_plots_data[[data_index]],
                       c("Others", "Tumour", "Immune"),
                       c("lightgray", "orange", "skyblue"))


## Bar plot for comparison between 3D and a single 2D slice
metrics_data1 <- metrics_data[metrics_data$name %in% c("3D", "slice4"), ]
metrics_data1$name[2] <- "2D"
colnames(metrics_data1)[6] <- "ACIN"
metrics_list_names <- c("APD", "AMD", "MS", "NMS", "ACIN", "CKAUC")

plot_result <- reshape2::melt(metrics_data1, id.vars = "name", mesaure.vars = metrics_list_names)
colnames(plot_result) <- c("name", "metric", "value")

plot_result$metric <- factor(plot_result$metric, levels = metrics_list_names)
plot <- ggplot(plot_result, aes(name, value, fill = name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~metric, scales = "free_y", ncol = 3)
plot

# 
# ggplot(metrics_data1, aes(name, APD, fill = name, ymin = APD - 2, ymax = APD + 2)) +
#   geom_bar(stat = "identity") +
#   scale_y_continuous(limits = c(0, 50)) +
#   labs(title = "APD", y= "value", x = "") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   geom_errorbar(width = 0.2)




metrics_data <- get_all_data(all_plots_data = all_plots_data,
                             all_plots_meta_data = all_plots_meta_data,
                             data_index = data_index,
                             reference_cell_type = "Tumour",
                             plot_image = TRUE)

plot_slice_3D(all_plots_data = all_plots_data,
              all_plots_meta_data = all_plots_meta_data,
              data_index = data_index,
              slice_num = 4)

plot_slices_2D(all_plots_data = all_plots_data,
               all_plots_meta_data = all_plots_meta_data,
               data_index = data_index)



get_all_data <- function(all_plots_data, 
                         all_plots_meta_data, 
                         data_index,
                         reference_cell_type,
                         plot_image = TRUE) {
  
  plots_data <- all_plots_data[[data_index]]
  plots_meta_data <- all_plots_meta_data[data_index, ]
  
  ## Get slice data for 7 slices
  slice_data <- vector(mode = "list", length = 7)
  
  delta <- -30 # Position of slice
  thickness <- 5 # Thickness of slice
  
  for (i in 1:7) {
    ## Separate clusters
    if (substr(plots_meta_data$cluster_type, 1, 1) == "S") {
      slice_data[[i]] <-  plots_data[0.5 * plots_data$Cell.X.Position + 0.5 * plots_data$Cell.Y.Position - plots_data$Cell.Z.Position > -delta - thickness &
                                     0.5 * plots_data$Cell.X.Position + 0.5 * plots_data$Cell.Y.Position - plots_data$Cell.Z.Position < -delta + thickness, ]  
    }
    ## Ring or Mixed clusters
    else {
      slice_data[[i]] <-  plots_data[plots_data$Cell.Z.Position > 75 + delta - thickness &
                                     plots_data$Cell.Z.Position < 75 + delta + thickness, ]   
    }
    delta <- delta + 2 * thickness
  }
  
  
  
  ### Get metrics data for 3D data and 2D slices
  metrics_data <- data.frame(matrix(nrow = 8, ncol = 6))
  metrics_list_names <- c("APD", "AMD", "MS", "NMS", "CKAUC", "CIN")
  colnames(metrics_data) <-   metrics_list_names
  
  if (reference_cell_type == "Tumour") {
    target_cell_type <- "Immune"
  } 
  else {
    target_cell_type <- "Tumour"
  }
  
  radius <- 30
  
  for (i in 1:8) {
    
    ## Get metrics data for 3D data
    if (i == 1) {
      result <- get_metrics_data_3D(plots_data,
                                    reference_cell_type,
                                    target_cell_type,
                                    radius)
      
      metrics_data[i, ] <- result
    }
    ## Get metrics data for 2D slices
    else {
      result <- get_metrics_data_2D(slice_data[[i - 1]],
                                    reference_cell_type,
                                    target_cell_type,
                                    radius)
      
      metrics_data[i, ] <- result
    }
  }
  
  # Plot
  if (plot_image == TRUE) {
    metrics_data$name <- paste("slice", 0:7, sep = "")
    metrics_data[1, "name"] <- "3D"
    
    plot_result <- reshape2::melt(metrics_data, id.vars = "name", mesaure.vars = c("APD", "AMD", "MS", "NMS", "CKAUC", "CIN"))
    colnames(plot_result) <- c("name", "metric", "value")
    
    plot <- ggplot(plot_result, aes(name, value, group = 1)) +
            geom_line(color = "red") +
            geom_point() +
            facet_wrap(~metric, nrow = 7, scales = "free_y")
    print(plot)
  }
  
  print(plots_meta_data)
  
  return (metrics_data)
}


get_metrics_data_3D <- function(plots_data,
                                reference_cell_type,
                                target_cell_type,
                                radius) {
  
  ## Order: c("APD", "AMD", "MS", "NMS", "CKAUC", "CIN")
  answer <- c()
  
  ## Get average pairwise distance data
  APD_data <- calculate_pairwise_distances_between_cell_types3D(plots_data, cell_types_of_interest = c("Tumour", "Immune"))
  answer <- append(answer, mean(APD_data[APD_data$Pair %in% c("Tumour/Immune", "Immune/Tumour"), "Distance"]))
  
  ## Get average minimum distance data
  AMD_data <- calculate_minimum_distances_between_cell_types3D(plots_data, cell_types_of_interest = c("Tumour", "Immune"))
  answer <- append(answer, mean(AMD_data[AMD_data$Pair == "Tumour/Immune", "Distance"]))
  
  ## Get mixing score and normalised mixing score data
  MS_data <- calculate_mixing_scores3D(plots_data,
                                       reference_cell_types = reference_cell_type,
                                       target_cell_types = target_cell_type,
                                       radius = radius)
  
  answer <- append(answer, MS_data[["Mixing_score"]])
  answer <- append(answer, MS_data[["Normalised_mixing_score"]])
  
  ## Get cross K AUC data
  cross_k_data <- calculate_Kcross3D(plots_data,
                                     reference_cell_type = reference_cell_type,
                                     target_cell_type = target_cell_type,
                                     distance = radius,
                                     plot_results = FALSE)
  
  CKAUC_data <- calculate_AUC_of_Kcross3D(cross_k_data)
  
  answer <- append(answer, CKAUC_data)
  
  ## Get cells in neighborhood data
  CIN_data <- calculate_cells_in_neighborhood3D(plots_data,
                                                reference_cell_types = reference_cell_type,
                                                target_cell_types = target_cell_type,
                                                radius = radius)
  
  answer <- append(answer, mean(CIN_data[[reference_cell_type]][[target_cell_type]]))
  

  return (answer)
}


get_metrics_data_2D <- function(plots_data,
                                reference_cell_type,
                                target_cell_type,
                                radius) {
  
  rownames(plots_data) <- plots_data$Cell.ID
  plots_data <- plots_data[ , c("Cell.X.Position", "Cell.Y.Position", "Cell.Type")]
  plots_data <- format_colData_to_spe(plots_data)
  
  
  ## Order: c("APD", "AMD", "MS", "NMS", "CKAUC", "CIN")
  answer <- c()
  
  ## Get average pairwise distance data
  if (sum(plots_data$Cell.Type == reference_cell_type) > 2 && sum(plots_data$Cell.Type == target_cell_type) > 2) {
    APD_data <- calculate_pairwise_distances_between_celltypes(plots_data, 
                                                               cell_types_of_interest = c("Tumour", "Immune"),
                                                               feature_colname = "Cell.Type")
    answer <- append(answer, mean(APD_data[APD_data$Pair %in% c("Tumour/Immune", "Immune/Tumour"), "Distance"]))
    
    ## Get average minimum distance data
    AMD_data <- calculate_minimum_distances_between_celltypes(plots_data, 
                                                              cell_types_of_interest = c("Tumour", "Immune"),
                                                              feature_colname = "Cell.Type")
    answer <- append(answer, mean(AMD_data[AMD_data$Pair == "Tumour/Immune", "Distance"]))
    
    ## Get mixing score and normalised mixing score data
    MS_data <- mixing_score_summary(plots_data,
                                    reference_celltype = reference_cell_type,
                                    target_celltype = target_cell_type,
                                    radius = radius,
                                    feature_colname = "Cell.Type")
    
    answer <- append(answer, MS_data[["Mixing_score"]])
    answer <- append(answer, MS_data[["Normalised_mixing_score"]])
    
    ## Get cross K AUC data
    cross_k_data <- calculate_cross_functions(plots_data,
                                              cell_types_of_interest = c(reference_cell_type, target_cell_type),
                                              feature_colname = "Cell.Type",
                                              dist = 30,
                                              plot_results = FALSE)
    
    CKAUC_data <- AUC_of_cross_function(cross_k_data)
    answer <- append(answer, CKAUC_data)
    
    ## Get cells in neighborhood data
    CIN_data <- number_of_cells_within_radius(plots_data,
                                              reference_celltype = reference_cell_type,
                                              target_celltype = target_cell_type,
                                              radius = radius,
                                              feature_colname = "Cell.Type")
    answer <- append(answer, mean(CIN_data[[reference_cell_type]][[target_cell_type]]))
  }
  else {
    answer <- c(NA, NA, NA, NA, NA, NA)
  }
  
  return (answer)
}




## Build function that can take slice of 3D data, and plot 3D data for chosen slice number (1 - 7)
plot_slice_3D <- function(all_plots_data = all_plots_data,
                          all_plots_meta_data = all_plots_meta_data,
                          data_index = data_index,
                          slice_nums) {
  
  plots_data <- all_plots_data[[data_index]]
  plots_meta_data <- all_plots_meta_data[data_index, ]
  
  thickness <- 5 # Thickness of slice
  
  for (slice_num in slice_nums) {
    delta <- -30 + (slice_num - 1) * (2 * thickness) # Position of slice
    
    ## Separate clusters
    if (substr(plots_meta_data$cluster_type, 1, 1) == "S") {
      plots_data[0.5 * plots_data$Cell.X.Position + 0.5 * plots_data$Cell.Y.Position - plots_data$Cell.Z.Position > -delta - thickness &
                   0.5 * plots_data$Cell.X.Position + 0.5 * plots_data$Cell.Y.Position - plots_data$Cell.Z.Position < -delta + thickness, "Cell.Type"] <- "Slice"  
    }
    ## Ring or Mixed clusters
    else {
      plots_data[plots_data$Cell.Z.Position > 75 + delta - thickness &
                   plots_data$Cell.Z.Position < 75 + delta + thickness, "Cell.Type"] <- "Slice"   
    }
  }
  plot_cell_categories3D(plots_data, 
                         c("Tumour", "Immune", "Slice", "Others"),
                         c("orange", "skyblue", "tomato", "lightgray"))
}





## Build another function that can plot all 7 2D slice data at once.

plot_slices_2D <- function(all_plots_data = all_plots_data,
                           all_plots_meta_data = all_plots_meta_data,
                           data_index = data_index) {
  
  plots_data <- all_plots_data[[data_index]]
  plots_meta_data <- all_plots_meta_data[data_index, ]
  
  ## Get slice data for 7 slices
  slice_data <- vector(mode = "list", length = 7)
  
  delta <- -30 # Position of slice
  thickness <- 5 # Thickness of slice
  
  for (i in 1:7) {
    ## Separate clusters
    if (substr(plots_meta_data$cluster_type, 1, 1) == "S") {
      slice_data[[i]] <-  plots_data[0.5 * plots_data$Cell.X.Position + 0.5 * plots_data$Cell.Y.Position - plots_data$Cell.Z.Position > -delta - thickness &
                                     0.5 * plots_data$Cell.X.Position + 0.5 * plots_data$Cell.Y.Position - plots_data$Cell.Z.Position < -delta + thickness, 
                                     c("Cell.X.Position", "Cell.Y.Position", "Cell.Type")]  
    }
    ## Ring or Mixed clusters
    else {
      slice_data[[i]] <-  plots_data[plots_data$Cell.Z.Position > 75 + delta - thickness &
                                     plots_data$Cell.Z.Position < 75 + delta + thickness,
                                     c("Cell.X.Position", "Cell.Y.Position", "Cell.Type")]   
    }
    
    delta <- delta + 2 * thickness
  }
  
  
  result <- reshape2::melt(slice_data,
                           measure.vars = "Cell.Type")
  
  result[["variable"]] <- NULL
  colnames(result) <- c("Cell.X.Position", "Cell.Y.Position", "Cell.Type", "Slice.Num")
  result$Cell.Type <- ordered(result$Cell.Type, c("Others", "Tumour", "Immune"))
  
  p <- ggplot(result, aes(Cell.X.Position, Cell.Y.Position, color = Cell.Type)) +
        geom_point() +
        scale_colour_manual(values = c("lightgray","orange", "skyblue")) +
        facet_wrap(~Slice.Num, nrow = 2, ncol = 4)
  
  return (p)
}
  




### Extra stuff -------------------------------------------------------------

## Plot in 2D
data2D <- slice_data[[7]]
data2D$Cell.Type <- ordered(data2D$Cell.Type, levels = c("Others", "Tumour", "Immune"))
ggplot(data2D, aes(Cell.X.Position, Cell.Y.Position, color = Cell.Type)) +
  geom_point() +
  scale_colour_manual(values = c("lightgray", "orange", "skyblue"))



## Get Tumour and Immune cell count data (should differ)
for (i in 1:7) {
  nTumour <- sum(slice_data[[i]]$Cell.Type == "Tumour")
  nImmune <- sum(slice_data[[i]]$Cell.Type == "Immune")
  print(paste("nTumour:", nTumour, "nImmune:", nImmune))
}
