setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/Code/spaSim 3D/objects")
all_plots_data <- readRDS(file="all_plots_test_data.rda")
all_plots_meta_data <- readRDS(file="all_plots_meta_data.rda")


## Make SS -> S, EE -> E, NN -> N
## Make ss -> s, mm -> m, ll -> l
all_plots_meta_data[all_plots_meta_data$shape == "SS", "shape"] <- "S"
all_plots_meta_data[all_plots_meta_data$shape == "EE", "shape"] <- "E"
all_plots_meta_data[all_plots_meta_data$shape == "NN", "shape"] <- "N"
all_plots_meta_data[all_plots_meta_data$size == "ss", "size"] <- "s"
all_plots_meta_data[all_plots_meta_data$size == "mm", "size"] <- "m"
all_plots_meta_data[all_plots_meta_data$size == "ll", "size"] <- "l"

all_plots_meta_data <-
  all_plots_meta_data[all_plots_meta_data$shape == "S" |
                        all_plots_meta_data$shape == "E" |
                        all_plots_meta_data$shape == "N", ]

all_plots_meta_data <-
  all_plots_meta_data[all_plots_meta_data$size == "s" |
                        all_plots_meta_data$size == "m" |
                        all_plots_meta_data$size == "l", ]

all_plots_data <- all_plots_data[as.numeric(rownames(all_plots_meta_data))]



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
                           shape = "S",
                           size = "test",
                           reference_cell_type = "Tumour",
                           cluster_type = 'S3')


## 'cluster_type' is the variable to test
result <- get_metrics_list(all_plots_data,
                           all_plots_meta_data,
                           shape = "N",
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
      CKI_data <- CKI_data$Distance[nrow(CKI_data)] ## Choose the largest distance when crossing occurs
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
















library(ggplot2)
## 'shape' and 'size' is the variable to test
result <- get_metrics_list1(all_plots_data,
                            all_plots_meta_data,
                            shape = "test",
                            size = "test",
                            reference_cell_type = "Tumour",
                            cluster_type = 'M2')


## 'size' and 'cluster_type' is the variable to test
result <- get_metrics_list1(all_plots_data,
                            all_plots_meta_data,
                            shape = "S",
                            size = "test",
                            reference_cell_type = "Tumour",
                            cluster_type = 'test')


## 'shape' and 'cluster_type' is the variable to test
result <- get_metrics_list1(all_plots_data,
                            all_plots_meta_data,
                            shape = "test",
                            size = "l",
                            reference_cell_type = "Tumour",
                            cluster_type = 'test')


## Get plot_data and metrics_list
plot_data <- result[[1]]
metrics_list <- result[[2]]


### Create a new function which makes two of the parameters variables, and the other one constant
## e.g. cluster_type <- 'test1', shape <- 'test2' and size <- 's'
## cluster_type is the variable on the x-axis, shape is the variable that controls the colour of the line

get_metrics_list1 <- function(plots_data, 
                              plots_meta_data,
                              shape,
                              size,
                              reference_cell_type,
                              cluster_type,
                              plot_image = TRUE) {
  
  n_metrics <- 7
  metrics_list <- vector(mode = "list", length = n_metrics)
  metrics_list_names <- c("APD", "AMD", "MS", "NMS", "CKI", "CKAUC", "CIN")
  names(metrics_list) <- metrics_list_names
  
  # Get target cell type
  if (reference_cell_type == "Tumour") {
    target_cell_type <- "Immune"
  }
  else {
    target_cell_type <- "Tumour" 
  }
  
  
  ### Select data for input parameters
  if (shape == 'test' && size == 'test') {
    chosen_data <- (plots_meta_data$cluster_type == cluster_type)
    plots_data <- plots_data[chosen_data]
    plots_meta_data <- plots_meta_data[chosen_data, ]
    
    var1 <- "shape"
    var1_names <- unique(plots_meta_data[["shape"]])
    
    var2 <- "size"
    var2_names <- unique(plots_meta_data[["size"]])
  }
  
  else if (cluster_type == 'test' && size == 'test') {
    chosen_data <- (plots_meta_data$shape == shape | plots_meta_data$shape == paste(shape, shape, sep = ""))
    plots_data <- plots_data[chosen_data]
    plots_meta_data <- plots_meta_data[chosen_data, ]
    
    var1 <- "cluster_type"
    var1_names <- unique(plots_meta_data[["cluster_type"]])
    
    var2 <- "size"
    var2_names <- unique(plots_meta_data[["size"]])
  }
  
  else if (cluster_type == 'test' && shape == 'test') {
    chosen_data <- (plots_meta_data$size == size | plots_meta_data$size == paste(size, size, sep = ""))
    plots_data <- plots_data[chosen_data]
    plots_meta_data <- plots_meta_data[chosen_data, ]
    
    var1 <- "cluster_type"
    var1_names <- unique(plots_meta_data[["cluster_type"]])
    
    var2 <- "shape"
    var2_names <- unique(plots_meta_data[["shape"]])
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
  
  ## Combine meta data and metrics list data
  plots_meta_data <- cbind(plots_meta_data, metrics_list)
  
  # Plot results
  if (plot_image == TRUE) {
    
    ## Melt
    result <- reshape2::melt(plots_meta_data, measure.vars = metrics_list_names)
    
    colnames(result) <- c('cluster_type', 'shape', 'size', 'metric', 'value')
    
    result$value <- as.numeric(result$value)
    result[[var1]] <- ordered(result[[var1]], levels = var1_names)
    result[[var2]] <- ordered(result[[var2]], levels = var2_names)
    result$metric <- factor(result$metric, levels = names(metrics_list))
    
    plot <- ggplot(result, aes(x = factor(.data[[var1]]), y = value, group = .data[[var2]])) +
      facet_wrap(~metric, nrow = n_metrics, scales = "free_y") +
      geom_line(aes(color = .data[[var2]])) + 
      scale_color_manual(values = RColorBrewer::brewer.pal(length(var2_names), "Set1"))
    
    print(plot)
  }
  
  
  ## return both chosen plots_data and metrics_list
  answer <- list()
  answer[[1]] <- plots_data
  answer[[2]] <- metrics_list
  
  return (answer)
}





















## 'shape' and 'cluster_type' is the variable to test
result <- get_metrics_list2(all_plots_data,
                            all_plots_meta_data,
                            reference_cell_type = "Immune")


## Get plot_data and metrics_list
plot_data <- result[[1]]
metrics_list <- result[[2]]


### Create a new function which compares all parameters: shape size & cluster_type

get_metrics_list2 <- function(plots_data, 
                              plots_meta_data,
                              reference_cell_type,
                              plot_image = TRUE) {
  
  n_metrics <- 7
  metrics_list <- vector(mode = "list", length = n_metrics)
  metrics_list_names <- c("APD", "AMD", "MS", "NMS", "CKI", "CKAUC", "CIN")
  names(metrics_list) <- metrics_list_names
  
  # Get target cell type
  if (reference_cell_type == "Tumour") {
    target_cell_type <- "Immune"
  }
  else {
    target_cell_type <- "Tumour" 
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
  
  ## Combine meta data and metrics list data
  plots_meta_data <- cbind(plots_meta_data, metrics_list)
  
  # Plot results
  if (plot_image == TRUE) {
    
    ## Melt
    result <- reshape2::melt(plots_meta_data, measure.vars = metrics_list_names)
    
    colnames(result) <- c('cluster_type', 'shape', 'size', 'metric', 'value')
    
    result$cluster_type <- ordered(result$cluster_type, levels = unique(result$cluster_type))
    result$shape <- ordered(result$shape, levels = unique(result$shape))
    result$size <- ordered(result$size, levels = unique(result$size))
    result$metric <- ordered(result$metric, levels = names(metrics_list))
    
    plot <- ggplot(result, aes(x = cluster_type, y = value, group = size)) +
      facet_grid(rows = vars(metric), cols = vars(shape), scales = "free_y") +
      geom_line(aes(color = size)) + 
      scale_color_manual(values = RColorBrewer::brewer.pal(length(unique(result$size)), "Set1"))
    
    print(plot)
  }
  
  
  ## return both chosen plots_data and metrics_list
  answer <- list()
  answer[[1]] <- plots_data
  answer[[2]] <- metrics_list
  
  return (answer)
}
