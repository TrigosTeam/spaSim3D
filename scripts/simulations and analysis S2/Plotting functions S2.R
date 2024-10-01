library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)

### Function for non-gradient output ------------------------------------------
plot_non_gradient_metric <- function(spes_table, 
                                     metric, 
                                     metric_df, 
                                     arrangement, 
                                     plots_metadata) {
  
  ### Modify plots_metadata
  # Change plots_metadata arrangement to inputted arrangement
  plots_metadata$arrangement$x_aes <- arrangement
  
  # Change plots_metadata y_aes to inputted metric
  for (i in seq_along(plots_metadata)) {
    # Modify the y_aes element
    plots_metadata[[i]]$y_aes <- metric
  }
  
  
  # Get metric_cell_types
  if (metric == "AMD") {
    metric_cell_types <- data.frame(ref = c("A", "A", "B", "B"), tar = c("A", "B", "A", "B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("prop_SAC", "prop_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("entropy_SAC", "entropy_AUC")) {
    metric_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  }
  
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, title = "") {
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw()
    
    # Use scientific notation for ellipsoid volume
    if (x_aes == "E_volume") {
      plot <- plot + 
        geom_point() +
        scale_x_continuous(labels = formatCustomSci)
    }
    else if (typeof(data[[x_aes]]) == "double") {
      breaks <- pretty(c(min(data[[x_aes]]), max(data[[x_aes]])), n = 2)
      plot <- plot + 
        geom_point() + 
        scale_x_continuous(breaks = breaks)
    }
    # Factored character is an integer
    else if (typeof(data[[x_aes]]) %in% c("integer", "character")) {
      plot <- plot + 
        geom_violin()
    }
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    if (metric %in% c("AMD", "prop_SAC", "prop_AUC")) {
      plot_df <- metric_df[metric_df$reference == metric_cell_types[i, "ref"] & metric_df$target == metric_cell_types[i, "tar"], ] 
    }
    else if (metric %in% c("entropy_SAC", "entropy_AUC")) {
      plot_df <- metric_df[metric_df$cell_types == metric_cell_types[i, "cell_types"], ]
    }
    
    # Combine spes_table and metric_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- plot_def$x_aes
      y_aes <- plot_def$y_aes
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, title = title)
      return(plot)
    })
  }
  
  # Combine the plots together using metric_cell_types
  combined_plots_list <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    # Get final column
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list[[cells]], nrow = 1, ncol = length(plots_list[[cells]]))
    
    if (metric %in% c("AMD", "prop_SAC", "prop_AUC")) {
      title <- ggdraw() +
        draw_label(paste("Reference:", metric_cell_types[i, "ref"], "Target:", metric_cell_types[i, "tar"]),
                   fontface = 'bold')
    }
    else if (metric %in% c("entropy_SAC", "entropy_AUC")) {
      title <- ggdraw() + 
        draw_label(paste("Cell types of interest:", cells), 
                   fontface='bold')
    }
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    combined_plots_list[[cells]] <- fig
  }
  
  # Combine the combined plots into one big plot
  non_gradient_metric_plot <- plot_grid(plotlist = combined_plots_list,
                                        nrow = length(combined_plots_list), 
                                        ncol = 1)
  
  return(non_gradient_metric_plot)
}

### Function for gradient output ----------------------------------------------
plot_gradient_metric <- function(spes_table, 
                                 metric, 
                                 metric_df, 
                                 arrangement, 
                                 gradient_type,
                                 plots_metadata) {
  ### Modify plots_metadata
  # Change plots_metadata arrangement to inputted arrangement
  plots_metadata$arrangement$color_aes <- arrangement
  
  # Change plots_metadata x_aes to gradient_type and y_aes to inputted metric
  for (i in seq_along(plots_metadata)) {
    plots_metadata[[i]]$x_aes <- gradient_type
    plots_metadata[[i]]$y_aes <- metric
  }
  
  # Get gradient/threshold values
  if (gradient_type == "radius") {
    gradient <- seq(20, 100, 10)
    gradient_colnames <- paste("r", gradient, sep = "")    
  }
  else if (gradient_type == "threshold") {
    gradient <- seq(0.01, 1, 0.01)
    gradient_colnames <- paste("t", gradient, sep = "")
  }
  
  # Get metric_cell_types
  if (metric %in% c("MS", "NMS")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("B", "A"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  if (metric %in% c("ACINP")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("A", "A"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  if (metric %in% c("AE")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("A,B", "A,B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("ACIN", "CKR")) {
    metric_cell_types <- data.frame(ref = c("A", "A", "B", "B"), tar = c("A", "B", "A", "B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("prop_prevalence")) {
    metric_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("entropy_prevalence")) {
    metric_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  }
  
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  create_plot <- function(data, x_aes, y_aes, color_aes, group_aes, title = "") {
    
    breaks <- pretty(c(min(data[[color_aes]]), max(data[[color_aes]])), n = 3)
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes, group = group_aes, color = color_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw() +
      geom_line()
    
    # Use scientific notation for ellipsoid volume
    if (color_aes == "E_volume") {
      plot <- plot + scale_color_continuous(labels = formatCustomSci)
    }
    else {
      plot <- plot + scale_color_continuous(breaks = breaks)
    }
    
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    if (metric %in% c("MS", "NMS", "ACIN", "CKR", "prop_prevalence")) {
      plot_df <- metric_df[metric_df$reference == metric_cell_types[i, "ref"] & metric_df$target == metric_cell_types[i, "tar"], ] 
    }
    else if (metric %in% c("ACINP", "AE")) {
      plot_df <- metric_df[metric_df$reference == metric_cell_types[i, "ref"], ]
    }
    else if (metric %in% c("entropy_prevalence")) {
      plot_df <- metric_df[metric_df$cell_types == metric_cell_types[i, "cell_types"], ]
    }
    else {
      stop("Invalid metric")
    }
    # Combine spes_table and metric_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , gradient_colnames)
    
    # Change last 2 column names
    colnames(plot_df)[c(ncol(plot_df) - 1, ncol(plot_df))] <- c(gradient_type, metric)
    
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df[[gradient_type]] <- unfactor(plot_df[[gradient_type]])
    plot_df[[gradient_type]] <- as.numeric(substr(plot_df[[gradient_type]], 2, nchar(plot_df[[gradient_type]])))
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    if (!is.null(plot_df$slice)) {
      plot_df$slice <- as.character(plot_df$slice)
      plot_df$key <- paste(plot_df$spe, plot_df$slice, sep = "_")
      group_aes = "key"
    }
    else {
      group_aes = "spe"
    }
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- plot_def$x_aes
      y_aes <- plot_def$y_aes
      color_aes <- plot_def$color_aes
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, group_aes = group_aes, color_aes = color_aes, title = title)
      return(plot)
    })
  }
  
  # Extract legends from first set of plots
  legends_list <- lapply(plots_list[[1]], function(plot) {
    plot_legend <- get_legend(plot + theme(legend.direction = "horizontal"))
    return(plot_legend)
  })
  legends <- plot_grid(plotlist = legends_list, nrow = 1)
  
  # Combine the plots together using metric_cell_types
  combined_plots_list <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    # Remove legend from base plots
    for (j in seq(length(plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]]))) {
      plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] <- 
        plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] + theme(legend.position = "none")
    }
    
    # Getting current set of cell types from metric_cell_types
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list[[cells]], nrow = 1, ncol = length(plots_list[[cells]]))
    
    if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_prevalence")) {
      title <- ggdraw() +
        draw_label(paste("Reference:", metric_cell_types[i, "ref"], "Target:", metric_cell_types[i, "tar"]),
                   fontface = 'bold')
    }
    else if (metric %in% c("entropy_prevalence")) {
      title <- ggdraw() + 
        draw_label(paste("Cell types of interest:", cells), 
                   fontface='bold')
    }
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    combined_plots_list[[cells]] <- fig
  }
  
  # Combine the combined plots into one big plot
  gradient_metric_plot <- plot_grid(plotlist = combined_plots_list,
                                    nrow = length(combined_plots_list), ncol = 1)
  
  # Add legends
  gradient_metric_with_legends_plot <- plot_grid(gradient_metric_plot, legends,
                                                 nrow = 2, ncol = 1,
                                                 rel_heights = c(1, 0.1))
  
  return(gradient_metric_with_legends_plot)
}


### Function to compare 3D vs 2D one slice ------------------------------------------
plot_3D_vs_2D_metric_one_slice <- function(spes_table, 
                                           metric, 
                                           metric_df3D,
                                           metric_df2D,
                                           arrangement, 
                                           plots_metadata) {
  
  ### Modify plots_metadata
  # Change plots_metadata arrangement to inputted arrangement
  plots_metadata$arrangement$color_aes <- arrangement
  
  # Get metric_cell_types
  if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "A", "B", "B"), tar = c("A", "B", "A", "B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("MS_AUC", "NMS_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("B", "A"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("ACINP_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("A", "A"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("AE_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("A,B", "A,B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("prop_SAC", "prop_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("entropy_SAC", "entropy_AUC")) {
    metric_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  }
  
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, color_aes, title = "") {
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes, color = color_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw() +
      xlim(min(c(data[[x_aes]], data[[y_aes]])), max(c(data[[x_aes]], data[[y_aes]]))) +
      ylim(min(c(data[[x_aes]], data[[y_aes]])), max(c(data[[x_aes]], data[[y_aes]]))) +
      geom_abline(intercept = 0, slope = 1, color = "black", linetype = "longdash")
    
    
    # Use scientific notation for ellipsoid volume
    if (color_aes == "E_volume") {
      plot <- plot + 
        geom_point() +
        scale_color_continuous(labels = formatCustomSci)
    }
    else if (typeof(data[[color_aes]]) == "double") {
      breaks <- pretty(c(min(data[[color_aes]]), max(data[[color_aes]])), n = 3)
      plot <- plot + 
        geom_point() + 
        scale_color_continuous(breaks = breaks)
    }
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC", "MS_AUC", "NMS_AUC", "prop_SAC", "prop_AUC")) {
      plot_df <- metric_df3D[metric_df3D$reference == metric_cell_types[i, "ref"] & metric_df3D$target == metric_cell_types[i, "tar"], ] 
      colnames(plot_df)[colnames(plot_df) == metric] <- paste(metric, "3D", sep = "_")
      plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D[metric_df2D$reference == metric_cell_types[i, "ref"] & metric_df2D$target == metric_cell_types[i, "tar"], metric] 
    }
    else if (metric %in% c("ACINP_AUC", "AE_AUC")) {
      plot_df <- metric_df3D[metric_df3D$reference == metric_cell_types[i, "ref"], ] 
      colnames(plot_df)[colnames(plot_df) == metric] <- paste(metric, "3D", sep = "_")
      plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D[metric_df2D$reference == metric_cell_types[i, "ref"], metric] 
      
    }
    else if (metric %in% c("entropy_SAC", "entropy_AUC")) {
      plot_df <- metric_df3D[metric_df3D$cell_types == metric_cell_types[i, "cell_types"], ]
      colnames(plot_df)[colnames(plot_df) == metric] <- paste(metric, "3D", sep = "_")
      plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D[metric_df2D$cell_types == metric_cell_types[i, "cell_types"], metric] 
    }
    
    # Combine spes_table and metric_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- paste(metric, "3D", sep = "_")
      y_aes <- paste(metric, "2D", sep = "_")
      color_aes <- plot_def$color_aes
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, color_aes = color_aes, title = title)
      return(plot)
    })
  }
  
  # Extract legends from first set of plots
  legends_list <- lapply(plots_list[[1]], function(plot) {
    plot_legend <- get_legend(plot + theme(legend.direction = "horizontal"))
    return(plot_legend)
  })
  legends <- plot_grid(plotlist = legends_list, nrow = 1)
  
  
  # Combine the plots together using metric_cell_types
  combined_plots_list <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    # Remove legend from base plots
    for (j in seq(length(plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]]))) {
      plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] <- 
        plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] + theme(legend.position = "none")
    }
    
    # Getting current set of cell types from metric_cell_types
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list[[cells]], nrow = 1, ncol = length(plots_list[[cells]]))
    
    if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "prop_SAC", "prop_AUC")) {
      title <- ggdraw() +
        draw_label(paste("Reference:", metric_cell_types[i, "ref"], "Target:", metric_cell_types[i, "tar"]),
                   fontface = 'bold')
    }
    else if (metric %in% c("entropy_SAC", "entropy_AUC")) {
      title <- ggdraw() + 
        draw_label(paste("Cell types of interest:", cells), 
                   fontface='bold')
    }
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    combined_plots_list[[cells]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plots <- plot_grid(plotlist = combined_plots_list,
                              nrow = length(combined_plots_list), 
                              ncol = 1)
  
  combined_plots_with_legends <- plot_grid(combined_plots, 
                                           legends,
                                           nrow = 2, ncol = 1,
                                           rel_heights = c(1, 0.1))
  return(combined_plots_with_legends)
}


### Function to compare 3D vs 2D all slices ------------------------------------------
plot_3D_vs_2D_metric_all_slices <- function(spes_table, 
                                            metric, 
                                            metric_df3D,
                                            metric_df2D,
                                            arrangement, 
                                            plots_metadata) {
  
  ### Modify plots_metadata
  # Change plots_metadata arrangement to inputted arrangement
  plots_metadata$arrangement$label <- arrangement  
  
  # Get metric_cell_types
  if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "A", "B", "B"), tar = c("A", "B", "A", "B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("MS_AUC", "NMS_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("B", "A"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("ACINP_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("A", "A"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("AE_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("A,B", "A,B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("prop_SAC", "prop_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("entropy_SAC", "entropy_AUC")) {
    metric_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  }
  
  # Get number of slices
  n_slices <- length(unique(metric_df2D[["slice"]]))
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, color_aes, label, title = "") {
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes, color = color_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw() +
      xlim(min(c(data[[x_aes]], data[[y_aes]]), na.rm = T), max(c(data[[x_aes]], data[[y_aes]]), na.rm = T)) +
      ylim(min(c(data[[x_aes]], data[[y_aes]]), na.rm = T), max(c(data[[x_aes]], data[[y_aes]]), na.rm = T)) +
      geom_abline(intercept = 0, slope = 1, color = "black", linetype = "longdash") +
      geom_point()
    
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC", "MS_AUC", "NMS_AUC", "prop_SAC", "prop_AUC")) {
      plot_df <- metric_df3D[metric_df3D$reference == metric_cell_types[i, "ref"] & metric_df3D$target == metric_cell_types[i, "tar"], ] 
      colnames(plot_df)[colnames(plot_df) == metric] <- paste(metric, "3D", sep = "_")
      plot_df <- cbind(spes_table, plot_df)
      
      plot_df <- plot_df %>%
        mutate(row_num = row_number())
      plot_df <- do.call(bind_rows, replicate(n_slices, plot_df, simplify = FALSE)) %>%
        arrange(row_num)
      plot_df$row_num <- NULL
      
      plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D[metric_df2D$reference == metric_cell_types[i, "ref"] & metric_df2D$target == metric_cell_types[i, "tar"], metric]
      plot_df[["slice"]] <- metric_df2D[metric_df2D$reference == metric_cell_types[i, "ref"] & metric_df2D$target == metric_cell_types[i, "tar"], "slice"]
    }
    else if (metric %in% c("ACINP_AUC", "AE_AUC")) {
      plot_df <- metric_df3D[metric_df3D$reference == metric_cell_types[i, "ref"], ] 
      colnames(plot_df)[colnames(plot_df) == metric] <- paste(metric, "3D", sep = "_")
      plot_df <- cbind(spes_table, plot_df)
      
      plot_df <- plot_df %>%
        mutate(row_num = row_number())
      plot_df <- do.call(bind_rows, replicate(n_slices, plot_df, simplify = FALSE)) %>%
        arrange(row_num)
      plot_df$row_num <- NULL
      
      plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D[metric_df2D$reference == metric_cell_types[i, "ref"], metric] 
      plot_df[["slice"]] <- metric_df2D[metric_df2D$reference == metric_cell_types[i, "ref"], "slice"] 
      
    }
    else if (metric %in% c("entropy_SAC", "entropy_AUC")) {
      plot_df <- metric_df3D[metric_df3D$cell_types == metric_cell_types[i, "cell_types"], ]
      colnames(plot_df)[colnames(plot_df) == metric] <- paste(metric, "3D", sep = "_")
      plot_df <- cbind(spes_table, plot_df)
      
      plot_df <- plot_df %>%
        mutate(row_num = row_number())
      plot_df <- do.call(bind_rows, replicate(n_slices, plot_df, simplify = FALSE)) %>%
        arrange(row_num)
      plot_df$row_num <- NULL
      
      plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D[metric_df2D$cell_types == metric_cell_types[i, "cell_types"], metric]
      plot_df[["slice"]] <- metric_df2D[metric_df2D$cell_types == metric_cell_types[i, "cell_types"], "slice"]
    }
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- paste(metric, "3D", sep = "_")
      y_aes <- paste(metric, "2D", sep = "_")
      color_aes <- plot_def$color_aes
      label <- plot_def$label
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, color_aes = color_aes, label = label, title = title)
      return(plot)
    })
  }
  
  # Extract legends from first set of plots
  legends_list <- lapply(plots_list[[1]], function(plot) {
    plot_legend <- get_legend(plot + theme(legend.direction = "horizontal"))
    return(plot_legend)
  })
  legends <- plot_grid(plotlist = legends_list, nrow = 1)
  
  
  # Combine the plots together using metric_cell_types
  combined_plots_list <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    # Remove legend from base plots
    for (j in seq(length(plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]]))) {
      plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] <- 
        plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] + theme(legend.position = "none")
    }
    
    # Getting current set of cell types from metric_cell_types
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list[[cells]], nrow = 1, ncol = length(plots_list[[cells]]))
    
    if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "prop_SAC", "prop_AUC")) {
      title <- ggdraw() +
        draw_label(paste("Reference:", metric_cell_types[i, "ref"], "Target:", metric_cell_types[i, "tar"]),
                   fontface = 'bold')
    }
    else if (metric %in% c("entropy_SAC", "entropy_AUC")) {
      title <- ggdraw() + 
        draw_label(paste("Cell types of interest:", cells), 
                   fontface='bold')
    }
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    combined_plots_list[[cells]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plots <- plot_grid(plotlist = combined_plots_list,
                              nrow = length(combined_plots_list), 
                              ncol = 1)
  
  # Get labels, if specified
  labels_vec <- unlist(lapply(plots_metadata, function(x) {
    return(x$label)
  }))
  labels <- list()
  for (label in labels_vec) {
    label <- ggdraw() +
      draw_label(label, fontface = 'bold')
    labels[[length(labels) + 1]] <- label
  }
  labels <- plot_grid(plotlist = labels, nrow = 1)
  combined_plots_with_legends_and_labels <- plot_grid(combined_plots, 
                                                      legends,
                                                      labels,
                                                      nrow = 3, ncol = 1,
                                                      rel_heights = c(1, 0.1, 0.1))
  return(combined_plots_with_legends_and_labels)
}

### Function to compare 3D vs 2D (not indexing for arrangement or shape) all slices ---------------
plot_3D_vs_2D_metric_all_slices1 <- function(metric, 
                                             metric_df3D,
                                             metric_df2D,
                                             plots_metadata) {
  
  # Get metric_cell_types
  if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "A", "B", "B"), tar = c("A", "B", "A", "B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("MS_AUC", "NMS_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("B", "A"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("ACINP_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("A", "A"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("AE_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "B"), tar = c("A,B", "A,B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("prop_SAC", "prop_AUC")) {
    metric_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("entropy_SAC", "entropy_AUC")) {
    metric_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  }
  
  # Get number of slices
  n_slices <- length(unique(metric_df2D[["slice"]]))
  
  # Define plotting function
  formatCustomSci <- function(x) {
    x_sci <- str_split_fixed(formatC(x, format = "e"), "e", 2)
    alpha <- as.numeric(x_sci[ , 1])
    power <- as.integer(x_sci[ , 2])
    paste(alpha, power, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, color_aes, title = "") {
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes, color = color_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw() +
      xlim(min(c(data[[x_aes]], data[[y_aes]]), na.rm = T), max(c(data[[x_aes]], data[[y_aes]]), na.rm = T)) +
      ylim(min(c(data[[x_aes]], data[[y_aes]]), na.rm = T), max(c(data[[x_aes]], data[[y_aes]]), na.rm = T)) +
      geom_abline(intercept = 0, slope = 1, color = "black", linetype = "longdash") +
      geom_point()
    
    return(plot)
  }
  
  # Put plots into an organised list
  plots_list <- list()
  
  for (i in seq(nrow(metric_cell_types))) {
    
    # Subset metric_df for chosen pair/cell types
    if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC", "MS_AUC", "NMS_AUC", "prop_SAC", "prop_AUC")) {
      plot_df <- metric_df3D[metric_df3D$reference == metric_cell_types[i, "ref"] & metric_df3D$target == metric_cell_types[i, "tar"], ] 
      colnames(plot_df)[colnames(plot_df) == metric] <- paste(metric, "3D", sep = "_")

      plot_df <- plot_df %>%
        mutate(row_num = row_number())
      plot_df <- do.call(bind_rows, replicate(n_slices, plot_df, simplify = FALSE)) %>%
        arrange(row_num)
      plot_df$row_num <- NULL
      
      plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D[metric_df2D$reference == metric_cell_types[i, "ref"] & metric_df2D$target == metric_cell_types[i, "tar"], metric]
      plot_df[["slice"]] <- metric_df2D[metric_df2D$reference == metric_cell_types[i, "ref"] & metric_df2D$target == metric_cell_types[i, "tar"], "slice"]
    }
    else if (metric %in% c("ACINP_AUC", "AE_AUC")) {
      plot_df <- metric_df3D[metric_df3D$reference == metric_cell_types[i, "ref"], ] 
      colnames(plot_df)[colnames(plot_df) == metric] <- paste(metric, "3D", sep = "_")
      
      plot_df <- plot_df %>%
        mutate(row_num = row_number())
      plot_df <- do.call(bind_rows, replicate(n_slices, plot_df, simplify = FALSE)) %>%
        arrange(row_num)
      plot_df$row_num <- NULL
      
      plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D[metric_df2D$reference == metric_cell_types[i, "ref"], metric] 
      plot_df[["slice"]] <- metric_df2D[metric_df2D$reference == metric_cell_types[i, "ref"], "slice"] 
      
    }
    else if (metric %in% c("entropy_SAC", "entropy_AUC")) {
      plot_df <- metric_df3D[metric_df3D$cell_types == metric_cell_types[i, "cell_types"], ]
      colnames(plot_df)[colnames(plot_df) == metric] <- paste(metric, "3D", sep = "_")
      
      plot_df <- plot_df %>%
        mutate(row_num = row_number())
      plot_df <- do.call(bind_rows, replicate(n_slices, plot_df, simplify = FALSE)) %>%
        arrange(row_num)
      plot_df$row_num <- NULL
      
      plot_df[[paste(metric, "2D", sep = "_")]] <- metric_df2D[metric_df2D$cell_types == metric_cell_types[i, "cell_types"], metric]
      plot_df[["slice"]] <- metric_df2D[metric_df2D$cell_types == metric_cell_types[i, "cell_types"], "slice"]
    }
    
    # Factor
    if (!is.null(plot_df$shape)) plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    if (!is.null(plot_df$slice)) plot_df$slice <- as.character(plot_df$slice)
    
    # Generate plots based on plots_metadata, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plots_metadata, function(plot_def) {
      x_aes <- paste(metric, "3D", sep = "_")
      y_aes <- paste(metric, "2D", sep = "_")
      color_aes <- plot_def$color_aes
      title <- plot_def$title
      plot <- create_plot(data = plot_df, x_aes = x_aes, y_aes = y_aes, color_aes = color_aes, title = title)
      return(plot)
    })
  }
  
  # Extract legends from first set of plots
  legends_list <- lapply(plots_list[[1]], function(plot) {
    plot_legend <- get_legend(plot + theme(legend.direction = "horizontal"))
    return(plot_legend)
  })
  legends <- plot_grid(plotlist = legends_list, nrow = 1)
  
  
  # Combine the plots together using metric_cell_types
  combined_plots_list <- list()
  for (i in seq(nrow(metric_cell_types))) {
    
    # Remove legend from base plots
    for (j in seq(length(plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]]))) {
      plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] <- 
        plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]][[j]] + theme(legend.position = "none")
    }
    
    # Getting current set of cell types from metric_cell_types
    cells <- metric_cell_types[i, ncol(metric_cell_types)]
    
    plots <- plot_grid(plotlist = plots_list[[cells]], nrow = 1, ncol = length(plots_list[[cells]]))
    
    if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "prop_SAC", "prop_AUC")) {
      title <- ggdraw() +
        draw_label(paste("Reference:", metric_cell_types[i, "ref"], "Target:", metric_cell_types[i, "tar"]),
                   fontface = 'bold')
    }
    else if (metric %in% c("entropy_SAC", "entropy_AUC")) {
      title <- ggdraw() + 
        draw_label(paste("Cell types of interest:", cells), 
                   fontface='bold')
    }
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    combined_plots_list[[cells]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plots <- plot_grid(plotlist = combined_plots_list,
                              nrow = length(combined_plots_list), 
                              ncol = 1)
  
  combined_plots_with_legends <- plot_grid(combined_plots, 
                                           legends,
                                           nrow = 2, ncol = 1,
                                           rel_heights = c(1, 0.1))
  return(combined_plots_with_legends)
}
