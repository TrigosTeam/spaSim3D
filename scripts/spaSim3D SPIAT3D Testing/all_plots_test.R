setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/Code/spaSim 3D/objects")
all_plots_data <- readRDS(file="all_plots_test_data.rda")
all_plots_meta_data <- readRDS(file="all_plots_meta_data.rda")


# For shape, size, reference_cell_type, cluster_type,
# One of these parameters must be set to the string: 'test'
# The other parameters must be fixed
# shape: "S", "E" or "N"
# size: "s", "m", or "l"
# reference_cell_type: "Tumour", or "Immune"
# cluster_type: "S1", "S2", "S3", "R1", "R2", "R3", "M1", "M2", "M3"

library(ggplot2)

## 'shape' is the variable to test
result <- get_metrics_list(all_plots_data,
                           all_plots_meta_data,
                           shape = "test",
                           size = "l",
                           reference_cell_type = "Tumour",
                           cluster_type = 'M2')


## 'size' is the variable to test
result <- get_metrics_list(all_plots_data,
                           all_plots_meta_data,
                           shape = "N",
                           size = "test",
                           reference_cell_type = "Tumour",
                           cluster_type = 'S3')


## 'cluster_type' is the variable to test
result <- get_metrics_list(all_plots_data,
                           all_plots_meta_data,
                           shape = "S",
                           size = "s",
                           reference_cell_type = "Tumour",
                           cluster_type = 'test')


## Get plot_data and metrics_list
plot_data <- result[[1]]
metrics_list <- result[[2]]

library(plotly)
plot_cell_categories3D(plot_data[[3]],
                       c("Others", "Tumour", "Immune"),
                       c("lightgray", "orange", "skyblue"))



get_metrics_list <- function(plots_data, 
                             plots_meta_data,
                             shape,
                             size,
                             reference_cell_type,
                             cluster_type,
                             plot_image = TRUE) {
  
  n_metrics <- 7
  metrics_list <- vector(mode = "list", length = n_metrics)
  names(metrics_list) <- c("APD", "AMD", "MS", "NMS", "CKI", "CKAUC", "CIN")
  
  # Get target cell type
  if (reference_cell_type == "Tumour") {
    target_cell_type <- "Immune"
  }
  else {
   target_cell_type <- "Tumour" 
  }
  
  
  ### Select data for input parameters
  if (shape == 'test') {
    chosen_data <- (plots_meta_data$size == size | plots_meta_data$size == paste(size, size, sep = "")) &
                   (plots_meta_data$cluster_type == cluster_type)
    
    plots_data <- plots_data[chosen_data]
    
    variable <- "shape"
    variable_names <- plots_meta_data[chosen_data, variable]
  }
  
  else if (size == 'test') {
    chosen_data <- (plots_meta_data$shape == shape | plots_meta_data$shape == paste(shape, shape, sep = "")) &
                   (plots_meta_data$cluster_type == cluster_type)
    
    plots_data <- plots_data[chosen_data]
    
    variable <- "size"
    variable_names <- plots_meta_data[chosen_data, variable]
  }
  
  else if (cluster_type == 'test') {
    chosen_data <- (plots_meta_data$shape == shape | plots_meta_data$shape == paste(shape, shape, sep = "")) &
                   (plots_meta_data$size == size | plots_meta_data$size == paste(size, size, sep = "")) # Include SS and ss
    
    plots_data <- plots_data[chosen_data]
    
    variable <- "cluster_type"
    variable_names <- plots_meta_data[chosen_data, variable]
  }
  
  ### Get metrics data
  radius <- 30
  
  for (plot_data in plots_data) {
    ## Get average pairwise distance data
    APD_data <- calculate_pairwise_distances_between_cell_types3D(plot_data, cell_types_of_interest = c("Tumour", "Immune"))
    metrics_list$APD <- append(metrics_list$APD, mean(APD_data[APD_data$Pair %in% c("Tumour/Immune", "Immune/Tumour"), "Distance"]))
    
    ## Get average minimum distance data
    AMD_data <- calculate_minimum_distances_between_cell_types3D(plot_data, cell_types_of_interest = c("Tumour", "Immune"))
    metrics_list$AMD <- append(metrics_list$AMD, mean(AMD_data[AMD_data$Pair == "Tumour/Immune", "Distance"]))
    
    ## Get mixing score and normalised mixing score data
    MS_data <- calculate_mixing_scores3D(plot_data,
                                         reference_cell_types = reference_cell_type,
                                         target_cell_types = target_cell_type,
                                         radius = radius)
    
    metrics_list$MS <- append(metrics_list$MS, MS_data[["Mixing_score"]])
    metrics_list$NMS <- append(metrics_list$NMS, MS_data[["Normalised_mixing_score"]])
    
    ## Get cross K intersection and cross K AUC data
    cross_k_data <- calculate_Kcross3D(plot_data,
                                       reference_cell_type = reference_cell_type,
                                       target_cell_type = target_cell_type,
                                       distance = radius,
                                       plot_results = FALSE)
    
    CKI_data <- calculate_Kcross_intersection3D(cross_k_data)
    if (length(CKI_data) == 2) {
      CKI_data <- CKI_data$Distance[1] ## Choose the first distance when crossing occurs
    }
    
    CKAUC_data <- calculate_AUC_of_Kcross3D(cross_k_data)
    
    metrics_list$CKI <- append(metrics_list$CKI, CKI_data)
    metrics_list$CKAUC <- append(metrics_list$CKAUC, CKAUC_data)
    
    ## Get cells in neighborhood data
    CIN_data <- calculate_cells_in_neighborhood3D(plot_data,
                                                  reference_cell_types = reference_cell_type,
                                                  target_cell_types = target_cell_type,
                                                  radius = radius)
    
    metrics_list$CIN <- append(metrics_list$CIN, mean(CIN_data[[reference_cell_type]][[target_cell_type]]))
  }
  
  metrics_list$MS <- as.numeric(metrics_list$MS)
  metrics_list$NMS <- as.numeric(metrics_list$NMS)
  
  # Plot results
  if (plot_image == TRUE) {
    
    result <- reshape2::melt(metrics_list)
    colnames(result) <- c('value', 'metric')
    result[ , "variable"] <- rep(variable_names, n_metrics)
    
    result$value <- as.numeric(result$value)
    result$variable <- ordered(result$variable, levels = variable_names)
    result$metric <- factor(result$metric, levels = names(metrics_list))
    
    plot <- ggplot(result, aes(x = factor(variable), y = value, group = 1)) +
              geom_line(color = "red") +
              facet_wrap(~metric, nrow = n_metrics, scales = "free_y")
              #scale_x_continuous(variable, breaks = 1:length(variable_names))
    
    print(plot)
  }
  
  ## return both chosen plots_data and metrics_list
  answer <- list()
  answer[[1]] <- plots_data
  answer[[2]] <- metrics_list

  return (answer)
}





### Create a new function which makes two of the parameters variables, and the other one constant
## e.g. cluster_type <- 'test1', shape <- 'test2' and size <- 's'
## cluster_type is the variable on the x-axis, shape is the variable that controls the colour of the line
