library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)

### Utility function to get metric cell types -----
get_metric_cell_types <- function(metric) {
  # Get metric_cell_types
  if (metric %in% c("AMD", "ACIN", "CKR", "ACIN_AUC", "CKR_AUC")) {
    metric_cell_types <- data.frame(ref = c("A"), tar = c("B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("MS", "NMS", "MS_AUC", "NMS_AUC")) {
    metric_cell_types <- data.frame(ref = c("A"), tar = c("B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("ACINP", "ACINP_AUC")) {
    metric_cell_types <- data.frame(ref = c("A"), tar = c("B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("AE", "AE_AUC")) {
    metric_cell_types <- data.frame(ref = c("A"), tar = c("A,B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("prop_SAC", "prop_prev", "prop_AUC")) {
    metric_cell_types <- data.frame(ref = c("A"), tar = c("B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("entropy_SAC", "entropy_prev", "entropy_AUC")) {
    metric_cell_types <- data.frame(cell_types = c("A,B"))
  }
  else {
    stop("metric not found")
  }
  return(metric_cell_types)
}


### Utility function to subset metric_df -----
subset_metric_df <- function(metric,
                             metric_df,
                             metric_cell_types,
                             index) {
  if (metric %in% c("AMD", "ACIN", "CKR", "MS", "NMS", "ACIN_AUC", "CKR_AUC", "MS_AUC", "NMS_AUC", "prop_SAC", "prop_prev", "prop_AUC")) {
    metric_df <- metric_df[metric_df$reference == metric_cell_types[index, "ref"] & metric_df$target == metric_cell_types[index, "tar"], ] 
  }
  else if (metric %in% c("ACINP", "AE", "ACINP_AUC", "AE_AUC")) {
    metric_df <- metric_df[metric_df$reference == metric_cell_types[index, "ref"], ] 
  }
  else if (metric %in% c("entropy_SAC", "entropy_prev", "entropy_AUC")) {
    metric_df <- metric_df[metric_df$cell_types == metric_cell_types[index, "cell_types"], ]
  }
  else {
    stop("metric not found")
  }
  return(metric_df)
}

### Utility function to duplicate df by *rows* -----
duplicate_df <- function(df, n_times) {
  df <- df %>%
    mutate(row_num = row_number())
  df <- do.call(bind_rows, replicate(n_times, df, simplify = FALSE)) %>%
    arrange(row_num)
  df$row_num <- NULL
  
  return(df)
}
### Function to compare 3D vs 2D and error vs 2D (not annotating for arrangement or shape) random slice ---------------
plot_3D_and_error_vs_2D_metric_random_slice_no_annotating <- function(metric, 
                                                                      metric_df3D,
                                                                      metric_df2D,
                                                                      plots_metadata) {
  
  # P-value multiplying factor
  p_value_factor <- 11 # number of metrics
  
  # Get metric_cell_types
  metric_cell_types <- get_metric_cell_types(metric)
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- round(as.numeric(x_sci[ , 1]), 1)
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot_3D_vs_2D <- function(data, x_aes, y_aes, title = "") {
    
    wilcox_test  <- wilcox.test(data[[x_aes]], data[[y_aes]], paired = TRUE)
    p_value <- wilcox_test$p.value * p_value_factor
    
    if (p_value > 1) p_value <- 1
    if (p_value == 0) p_value <- 2.2e-308
    if (0 < p_value && p_value < 1e-3)  {
      p_value <- formatCustomSci(p_value)
    }
    else {
      p_value  <- round(p_value, 3)
    }
    title <- paste("p =", p_value)
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes)) +
      labs(x = x_aes, y = y_aes) +
      theme_bw() +
      ggtitle(title) +
      theme(plot.title = element_text(size = 9)) +
      xlim(min(c(data[[x_aes]], data[[y_aes]]), na.rm = T), max(c(data[[x_aes]], data[[y_aes]]), na.rm = T)) +
      ylim(min(c(data[[x_aes]], data[[y_aes]]), na.rm = T), max(c(data[[x_aes]], data[[y_aes]]), na.rm = T)) +
      geom_point(size = 0.5) +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "longdash")
    
    return(plot)
  }
  create_plot_error_vs_2D <- function(data, x_aes, y_aes, title = "") {
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes)) +
      labs(title = title, x = x_aes, y = paste(metric, "error (%)")) +
      theme_bw() +
      # ylim(min(c(0, data[[y_aes]]), na.rm = T), max(c(0, data[[y_aes]]), na.rm = T)) +
      ylim(-200, 1000) +
      geom_point(size = 0.5) +
      geom_abline(intercept = 0, slope = 0, color = "red", linetype = "longdash")
    
    return(plot)
  }
  
  
  # Put plots into an organised list
  plots_list_3D_vs_2D <- list()
  plots_list_error_vs_2D <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    metric_df3D_subset <- subset_metric_df(metric, metric_df3D, metric_cell_types, i)
    metric_df2D_subset <- subset_metric_df(metric, metric_df2D, metric_cell_types, i)
    
    plot_df <- data.frame(row.names = rownames(metric_df3D_subset))
    plot_df[[paste(metric, "3D", sep = "_")]] <- metric_df3D_subset[[metric]]
    
    # Choose a random slice from metric_df2D_subset
    n_slices <- length(unique(metric_df2D[["slice"]]))
    
    metric_df2D_subset$key <- paste(metric_df2D_subset[["spe"]], metric_df2D_subset[["slice"]], sep = "_")
    plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D_subset[metric_df2D_subset$key %in% paste(unique(metric_df2D_subset[["spe"]]), sample(seq(n_slices), nrow(metric_df3D_subset), replace = TRUE), sep = "_"), 
                                                                    metric]
    
    # Add error
    plot_df[[paste(metric, "error", sep = "_")]] <- 100 * (plot_df[[paste(metric, "2D", sep = "_")]] - plot_df[[paste(metric, "3D", sep = "_")]]) / (plot_df[[paste(metric, "3D", sep = "_")]])
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list_3D_vs_2D[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- paste(metric, "3D", sep = "_")
      y_aes <- paste(metric, "2D", sep = "_")
      title <- plot_def$title
      plot <- create_plot_3D_vs_2D(data = plot_df, x_aes = x_aes, y_aes = y_aes, title = title)
      return(plot)
    })
    
    plots_list_error_vs_2D[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- paste(metric, "3D", sep = "_")
      y_aes <- paste(metric, "error", sep = "_")
      title <- plot_def$title
      plot <- create_plot_error_vs_2D(data = plot_df, x_aes = x_aes, y_aes = y_aes, title = title)
      return(plot)
    })
  }
  
  # Combine the plots together using metric_cell_types
  combined_plots_list_3D_vs_2D <- list()
  combined_plots_list_error_vs_2D <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    ## Part 1
    # Remove legend from base plots
    for (j in seq(length(plots_list_3D_vs_2D[[metric_cell_types[i, ncol(metric_cell_types)]]]))) {
      plots_list_3D_vs_2D[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] <- 
        plots_list_3D_vs_2D[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] + theme(legend.position = "none")
    }
    
    # Getting current set of cell types from metric_cell_types
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list_3D_vs_2D[[cells]], nrow = 1, ncol = length(plots_list_3D_vs_2D[[cells]]))
    
    combined_plots_list_3D_vs_2D[[cells]] <- plots
    
    ## Part 2
    # Remove legend from base plots
    for (j in seq(length(plots_list_error_vs_2D[[metric_cell_types[i, ncol(metric_cell_types)]]]))) {
      plots_list_error_vs_2D[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] <- 
        plots_list_error_vs_2D[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] + theme(legend.position = "none")
    }
    
    # Getting current set of cell types from metric_cell_types
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list_error_vs_2D[[cells]], nrow = 1, ncol = length(plots_list_error_vs_2D[[cells]]))
    
    combined_plots_list_error_vs_2D[[cells]] <- plots
  }
  
  # Combine the combined plots into one big plot
  combined_plots_3D_vs_2D <- plot_grid(plotlist = combined_plots_list_3D_vs_2D,
                                       nrow = length(combined_plots_list_3D_vs_2D), 
                                       ncol = 1)
  
  combined_plots_error_vs_2D <- plot_grid(plotlist = combined_plots_list_error_vs_2D,
                                          nrow = length(combined_plots_list_error_vs_2D), 
                                          ncol = 1)
  
  return(list(plots_3D_vs_2D = combined_plots_3D_vs_2D,
              plots_error_vs_2D = combined_plots_error_vs_2D))
}

library(cowplot)
library(ggplot2)
