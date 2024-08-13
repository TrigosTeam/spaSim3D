library(cowplot)
library(ggplot2)
mixed_arrangements <- c("M1", "M2", "M3")
ringed_arrangements <- c("R1", "R2", "R3")
all_arrangements <- c("M2", "R2", "S2")

### 1.1.1. Function to get plot for APD ----------------------------------------
### 1.1.2. Function to get plot for AMD -------------------------------------
plot_AMD_metric <- function(spes_table, AMD_df, bg_AMD_df, arrangements) {
  
  # AMD pairs are A/A, A/B, B/A, B/B
  AMD_pairs <- data.frame(cell1 = c("A", "A", "B", "B"),
                          cell2 = c("A", "B", "A", "B"))
  AMD_pairs$pair <- paste(AMD_pairs$cell1, AMD_pairs$cell2, sep = "/")
  
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(nrow(AMD_pairs))) {
    
    # Subset mixed_AMD_df for chosen pair
    plot_df <- AMD_df[AMD_df$reference == AMD_pairs[i, "cell1"] & AMD_df$target == AMD_pairs[i, "cell2"], ]
    
    # Combine mixed_spes_table and mixed_AMD_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Create a 'key' column which groups simulations if they have the same bg_type, shape and size (but not arrangement)
    plot_df$key <- paste(plot_df$bg_type, plot_df$shape, plot_df$size, sep = "_")
    
    # Get plot_bg_df from bg_AMD_df
    plot_bg_df <- data.frame(arrangement = arrangements, 
                             AMD = bg_AMD_df[bg_AMD_df$reference == AMD_pairs[i, "cell1"] & bg_AMD_df$target == AMD_pairs[i, "cell2"], "AMD"],
                             key = "bg")
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shape <- factor(plot_df$shape, c("Sphere", "Ellipsoid", "Network"))
    plot_df$size <- factor(plot_df$size, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    fig_bg_type <- ggplot(plot_df, aes(arrangement, AMD, group = key, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_line(data = plot_bg_df, aes(arrangement, AMD), col = "black", linetype = 2)
    
    fig_shape <- ggplot(plot_df, aes(arrangement, AMD, group = key, col = shape)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_line(data = plot_bg_df, aes(arrangement, AMD), col = "black", linetype = 2)
    
    fig_size <- ggplot(plot_df, aes(arrangement, AMD, group = key, col = size)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_line(data = plot_bg_df, aes(arrangement, AMD), col = "black", linetype = 2)
    
    all_plots_list[[AMD_pairs[i, "pair"]]] <- list(bg_type = fig_bg_type, shape = fig_shape, size = fig_size)
  }
  
  
  # Combine the plots together by pairs
  library(cowplot)
  plots_pair_list <- list()
  
  for (i in seq(nrow(AMD_pairs))) {
    pair <- AMD_pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$bg_type, all_plots_list[[pair]]$shape, all_plots_list[[pair]]$size, nrow = 1, ncol = 3)
    
    title <- ggdraw() + 
      draw_label(paste("Reference:", AMD_pairs[i, "cell1"], "Target:", AMD_pairs[i, "cell2"]), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_pair_list[[pair]] <- fig
  }
  
  # Combine the combined plots into one big plot
  AMD_plot <- plot_grid(plots_pair_list$`A/A`, 
                        plots_pair_list$`A/B`, 
                        plots_pair_list$`B/A`, 
                        plots_pair_list$`B/B`, 
                        nrow = 2, ncol = 2, scale = 0.9)
  
  methods::show(AMD_plot)
  
  return(AMD_plot)
}


### 1.2.1. Function to get plot for MS, NMS, ACINP, AE gradient metrics ------------

plot_gradient_metrics_type1 <- function(spes_table, gradient_metric_df, bg_gradient_metric_df, metric, arrangements) {
  
  # Constants
  cell_types <- c("A", "B") # Use A as reference, and B as target, and vice versa
  radii <- 50
  radii_colnames <- paste("r", seq(radii), sep = "")
  
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(length(cell_types))) {
    # Subset gradient_metric_df for current reference cell
    plot_df <- gradient_metric_df[gradient_metric_df$reference == cell_types[i], ]
    
    # Combine spes_table and mixed_AMD_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , radii_colnames)
  
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df$variable <- unfactor(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shape <- factor(plot_df$shape, c("Sphere", "Ellipsoid", "Network"))
    plot_df$size <- factor(plot_df$size, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    # Do the same for the bg_gradient_metric_df
    plot_bg_df <- bg_gradient_metric_df[bg_gradient_metric_df$reference == cell_types[i], ]
    plot_bg_df <- reshape2::melt(plot_bg_df, , radii_colnames)
    plot_bg_df$variable <- unfactor(plot_bg_df$variable)
    plot_bg_df$variable <- as.numeric(substr(plot_bg_df$variable, 2, nchar(plot_bg_df$variable)))
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = arrangement)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2)

    fig_bg_type <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2)
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_size <- ggplot(plot_df, aes(variable, value, group = spe, col = size)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2)
    
    all_plots_list[[cell_types[i]]] <- list(arrangement = fig_arrangement, bg_type = fig_bg_type, shape = fig_shape, size = fig_size)
    
  }
  
  
  # Combine the plots together by reference cell type
  plots_ref_list <- list()
  
  for (i in seq(length(cell_types))) {
    reference_cell_type <- cell_types[i]
    
    plots <- plot_grid(all_plots_list[[reference_cell_type]]$arrangement,
                       all_plots_list[[reference_cell_type]]$bg_type, 
                       all_plots_list[[reference_cell_type]]$shape, 
                       all_plots_list[[reference_cell_type]]$size, nrow = 1, ncol = 4)
    
    title <- ggdraw() + 
      draw_label(paste("Reference:", reference_cell_type), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_ref_list[[reference_cell_type]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plot <- plot_grid(plots_ref_list$A,
                             plots_ref_list$B, 
                             nrow = 2, ncol = 1, scale = 0.9)
  
  methods::show(combined_plot)
  
  return(combined_plot)
}






### 1.2.2. Function to get plot for ACIN, CKR gradient metrics ------------------

plot_gradient_metrics_type2 <- function(spes_table, gradient_metric_df, bg_gradient_metric_df, metric, arrangements, min_radius, max_radius) {
  
  # Constants
  pairs <- data.frame(cell1 = c("A", "A", "B", "B"),
                          cell2 = c("A", "B", "A", "B"))
  pairs$pair <- paste(pairs$cell1, pairs$cell2, sep = "/")
  
  radii <- 50
  radii_colnames <- paste("r", seq(radii), sep = "")
  
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(nrow(pairs))) {
    # Subset gradient_metric_df for current reference cell
    plot_df <- gradient_metric_df[gradient_metric_df$reference == pairs[i, "cell1"] &
                                    gradient_metric_df$target == pairs[i, "cell2"], ]
    
    # Combine spes_table and mixed_AMD_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , radii_colnames)
    
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df$variable <- unfactor(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    
    plot_df <- plot_df[plot_df$variable >= min_radius & plot_df$variable <= max_radius, ]
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shape <- factor(plot_df$shape, c("Sphere", "Ellipsoid", "Network"))
    plot_df$size <- factor(plot_df$size, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    # Do the same for the bg_gradient_metric_df
    plot_bg_df <- bg_gradient_metric_df[bg_gradient_metric_df$reference == pairs[i, "cell1"] &
                                          bg_gradient_metric_df$target == pairs[i, "cell2"], ]
    plot_bg_df <- reshape2::melt(plot_bg_df, , radii_colnames)
    plot_bg_df$variable <- unfactor(plot_bg_df$variable)
    plot_bg_df$variable <- as.numeric(substr(plot_bg_df$variable, 2, nchar(plot_bg_df$variable)))
    
    plot_bg_df <- plot_bg_df[plot_bg_df$variable >= min_radius & plot_bg_df$variable <= max_radius, ]
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = arrangement)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_bg_type <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_size <- ggplot(plot_df, aes(variable, value, group = spe, col = size)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    all_plots_list[[pairs[i, "pair"]]] <- list(arrangement = fig_arrangement, bg_type = fig_bg_type, shape = fig_shape, size = fig_size)
    
  }
  
  
  # Combine the plots together by reference cell type
  plots_pair_list <- list()
  
  for (i in seq(nrow(pairs))) {
    pair <- pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$arrangement,
                       all_plots_list[[pair]]$bg_type, 
                       all_plots_list[[pair]]$shape, 
                       all_plots_list[[pair]]$size, nrow = 1, ncol = 4)
    
    title <- ggdraw() + 
      draw_label(paste("Reference:", pairs[i, "cell1"], "Target:", pairs[i, "cell2"]), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_pair_list[[pair]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plot <- plot_grid(plots_pair_list$`A/A`, 
                             plots_pair_list$`A/B`, 
                             plots_pair_list$`B/A`, 
                             plots_pair_list$`B/B`, 
                             nrow = 2, ncol = 2, scale = 0.9)
  
  methods::show(combined_plot)
  
  return(combined_plot)
}


### 1.3.1. Function to get plot for proportion SAC ----------------------------------------
plot_proportion_SAC <- function(spes_table, SAC_df, bg_SAC_df, arrangements) {
  
  # Get possible reference and target cell combinations
  prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
  prop_cell_types$pair <- paste(prop_cell_types$ref, prop_cell_types$tar, sep = "/")
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    
    # Subset for current reference (and target) cell
    plot_df <- SAC_df[SAC_df$reference == prop_cell_types$ref[i], ]
    
    # Combine spes_table and SAC_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Create a 'key' column which groups simulations if they have the same bg_type, shape and size (but not arrangement)
    plot_df$key <- paste(plot_df$bg_type, plot_df$shape, plot_df$size, sep = "_")

    # Get plot_bg_df from bg_SAC_df
    plot_bg_df <- data.frame(arrangement = arrangements, 
                             proportion = bg_SAC_df[bg_SAC_df$reference == prop_cell_types$ref[i], "proportion"],
                             key = "bg")
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shape <- factor(plot_df$shape, c("Sphere", "Ellipsoid", "Network"))
    plot_df$size <- factor(plot_df$size, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    fig_bg_type <- ggplot(plot_df, aes(arrangement, proportion, group = key, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    fig_shape <- ggplot(plot_df, aes(arrangement, proportion, group = key, col = shape)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    fig_size <- ggplot(plot_df, aes(arrangement, proportion, group = key, col = size)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(bg_type = fig_bg_type, shape = fig_shape, size = fig_size)
  }
  
  
  # Combine the plots together by reference target pairs
  library(cowplot)
  plots_pair_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$bg_type, all_plots_list[[pair]]$shape, all_plots_list[[pair]]$size, nrow = 1, ncol = 3)
    
    title <- ggdraw() + 
      draw_label(paste("Reference/Target:", pair), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_pair_list[[pair]] <- fig
  }
  
  # Combine the combined plots into one big plot
  SAC_plot <- plot_grid(plots_pair_list[[prop_cell_types$pair[1]]], 
                        plots_pair_list[[prop_cell_types$pair[2]]],
                        nrow = 2, ncol = 1, scale = 0.9)
  
  methods::show(SAC_plot)
  
  return(SAC_plot)
}


### 1.3.2. Function to get plot for entropy SAC ----------------------------------------
plot_entropy_SAC <- function(spes_table, SAC_df, bg_SAC_df, arrangements) {
  
  # Get possible cell type of interest combinations
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    
    # Subset for current cell type of interest combintation
    plot_df <- SAC_df[SAC_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Combine spes_table and SAC_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Create a 'key' column which groups simulations if they have the same bg_type, shape and size (but not arrangement)
    plot_df$key <- paste(plot_df$bg_type, plot_df$shape, plot_df$size, sep = "_")
    
    # Get plot_bg_df from bg_SAC_df
    plot_bg_df <- data.frame(arrangement = arrangements, 
                             proportion = bg_SAC_df[bg_SAC_df$cell_types == entropy_cell_types$cell_types[i], "entropy"],
                             key = "bg")
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shape <- factor(plot_df$shape, c("Sphere", "Ellipsoid", "Network"))
    plot_df$size <- factor(plot_df$size, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    fig_bg_type <- ggplot(plot_df, aes(arrangement, entropy, group = key, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    fig_shape <- ggplot(plot_df, aes(arrangement, entropy, group = key, col = shape)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    fig_size <- ggplot(plot_df, aes(arrangement, entropy, group = key, col = size)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(bg_type = fig_bg_type, shape = fig_shape, size = fig_size)
  }
  
  # Combine the plots together by cell types of interest
  library(cowplot)
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$bg_type, all_plots_list[[cell_types]]$shape, all_plots_list[[cell_types]]$size, nrow = 1, ncol = 3)
    
    title <- ggdraw() + 
      draw_label(paste("Cell types of interest:", cell_types), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_cell_types_list[[cell_types]] <- fig
  }
  
  # Combine the combined plots into one big plot
  SAC_plot <- plot_grid(plots_cell_types_list[[entropy_cell_types$cell_types[1]]], 
                        plots_cell_types_list[[entropy_cell_types$cell_types[2]]],
                        nrow = 2, ncol = 1, scale = 0.9)
  
  methods::show(SAC_plot)
  
  return(SAC_plot)
}
### 1.3.3. Function to get plot for proportion prevalence ----------------------------------
plot_proportion_prevalence <- function(spes_table, prevalence_df, bg_df, arrangements) {
  
  # Constants
  prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
  prop_cell_types$pair <- paste(prop_cell_types$ref, prop_cell_types$tar, sep = "/")
  
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    # Subset prevalence_df for current metric
    plot_df <- prevalence_df[prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Combine spes_table and prevalence_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , threshold_colnames)
    
    # Extract threshold value from threshold strings (t0.01 -> 0.01...)
    plot_df$variable <- as.character(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shape <- factor(plot_df$shape, c("Sphere", "Ellipsoid", "Network"))
    plot_df$size <- factor(plot_df$size, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    # Do the same for the bg_df
    plot_bg_df <- bg_df[bg_df$reference == prop_cell_types$ref[i], ]
    plot_bg_df <- reshape2::melt(plot_bg_df, , threshold_colnames)
    plot_bg_df$variable <- unfactor(plot_bg_df$variable)
    plot_bg_df$variable <- as.numeric(substr(plot_bg_df$variable, 2, nchar(plot_bg_df$variable)))
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = arrangement)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_bg_type <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_size <- ggplot(plot_df, aes(variable, value, group = spe, col = size)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(arrangement = fig_arrangement, bg_type = fig_bg_type, shape = fig_shape, size = fig_size)
    
  }
  
  
  # Combine the plots together by reference (and target) cell pairs
  plots_pair_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$arrangement,
                       all_plots_list[[pair]]$bg_type, 
                       all_plots_list[[pair]]$shape, 
                       all_plots_list[[pair]]$size, nrow = 1, ncol = 4)
    
    title <- ggdraw() + 
      draw_label(paste("Reference/Target:", pair), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_pair_list[[pair]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plot <- plot_grid(plots_pair_list[[prop_cell_types$pair[1]]],
                             plots_pair_list[[prop_cell_types$pair[2]]], 
                             nrow = 2, ncol = 1, scale = 0.9)
  
  methods::show(combined_plot)
  
  return(combined_plot)
}



### 1.3.4. Function to get plot for entropy prevalence ----------------------------------
plot_entropy_prevalence <- function(spes_table, prevalence_df, bg_df, arrangements) {
  
  # Constants
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    # Subset prevalence_df for current cell type of interest
    plot_df <- prevalence_df[prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Combine spes_table and prevalence_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , threshold_colnames)
    
    # Extract threshold value from threshold strings (t0.01 -> 0.01...)
    plot_df$variable <- as.character(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shape <- factor(plot_df$shape, c("Sphere", "Ellipsoid", "Network"))
    plot_df$size <- factor(plot_df$size, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    # Do the same for the bg_df
    plot_bg_df <- bg_df[bg_df$cell_types == entropy_cell_types$cell_types[i], ]
    plot_bg_df <- reshape2::melt(plot_bg_df, , threshold_colnames)
    plot_bg_df$variable <- unfactor(plot_bg_df$variable)
    plot_bg_df$variable <- as.numeric(substr(plot_bg_df$variable, 2, nchar(plot_bg_df$variable)))
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = arrangement)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_bg_type <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_size <- ggplot(plot_df, aes(variable, value, group = spe, col = size)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(arrangement = fig_arrangement, bg_type = fig_bg_type, shape = fig_shape, size = fig_size)
    
  }
  
  # Combine the plots together by cell types of interest
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$arrangement,
                       all_plots_list[[cell_types]]$bg_type, 
                       all_plots_list[[cell_types]]$shape, 
                       all_plots_list[[cell_types]]$size, nrow = 1, ncol = 4)
    
    title <- ggdraw() + 
      draw_label(paste("Cell types of interest:", cell_types), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_cell_types_list[[cell_types]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plot <- plot_grid(plots_cell_types_list[[entropy_cell_types$cell_types[1]]],
                             plots_cell_types_list[[entropy_cell_types$cell_types[2]]], 
                             nrow = 2, ncol = 1, scale = 0.9)
  
  methods::show(combined_plot)
  
  return(combined_plot)
}



### 1.3.5. Function to get plot for proportion prevalence AUC -----------------
plot_proportion_prevalence_AUC <- function(spes_table, prevalence_df, bg_df, arrangements) {
  
  # Constants
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
  prop_cell_types$pair <- paste(prop_cell_types$ref, prop_cell_types$tar, sep = "/")

  # Get AUC for each prevalence gradient
  prevalence_df$AUC <- apply(prevalence_df[ , threshold_colnames], 1, sum) * 0.01
  prevalence_df <- prevalence_df[ , c("spe", "reference", "target", "AUC")]
  
  # Do the same for bg
  bg_df$AUC <- apply(bg_df[ , threshold_colnames], 1, sum) * 0.01
  bg_df <- bg_df[ , c("spe", "reference", "target", "AUC")]
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    
    # Subset for current reference (and target) cell
    plot_df <- prevalence_df[prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Combine spes_table and updated prevalence_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Create a 'key' column which groups simulations if they have the same bg_type, shape and size (but not arrangement)
    plot_df$key <- paste(plot_df$bg_type, plot_df$shape, plot_df$size, sep = "_")
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shape <- factor(plot_df$shape, c("Sphere", "Ellipsoid", "Network"))
    plot_df$size <- factor(plot_df$size, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    # Get plot_bg_df from bg_df
    plot_bg_df <- data.frame(arrangement = arrangements,
                             AUC = bg_df[bg_df$reference == prop_cell_types$ref[i], "AUC"],
                             key = "bg")
    
    fig_bg_type <- ggplot(plot_df, aes(arrangement, AUC, group = key, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("AUC") +
      geom_line(data = plot_bg_df, aes(arrangement, AUC), col = "black", linetype = 2)
    
    fig_shape <- ggplot(plot_df, aes(arrangement, AUC, group = key, col = shape)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("AUC") +
      geom_line(data = plot_bg_df, aes(arrangement, AUC), col = "black", linetype = 2)
    
    fig_size <- ggplot(plot_df, aes(arrangement, AUC, group = key, col = size)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("AUC") +
      geom_line(data = plot_bg_df, aes(arrangement, AUC), col = "black", linetype = 2)
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(bg_type = fig_bg_type, shape = fig_shape, size = fig_size)
  }
  
  
  # Combine the plots together by reference target pairs
  library(cowplot)
  plots_pair_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$bg_type, all_plots_list[[pair]]$shape, all_plots_list[[pair]]$size, nrow = 1, ncol = 3)
    
    title <- ggdraw() + 
      draw_label(paste("Reference/Target:", pair), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_pair_list[[pair]] <- fig
  }
  
  # Combine the combined plots into one big plot
  AUC_plot <- plot_grid(plots_pair_list[[prop_cell_types$pair[1]]], 
                        plots_pair_list[[prop_cell_types$pair[2]]],
                        nrow = 2, ncol = 1, scale = 0.9)
  
  methods::show(AUC_plot)
  
  return(AUC_plot)
}

### 1.3.6. Function to get plot for entropy prevalence AUC ----------------------------------------
plot_entropy_prevalence_AUC <- function(spes_table, prevalence_df, bg_df, arrangements) {
  
  # Constants
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))

  # Get AUC for each prevalence gradient
  prevalence_df$AUC <- apply(prevalence_df[ , threshold_colnames], 1, sum) * 0.01
  prevalence_df <- prevalence_df[ , c("spe", "cell_types", "AUC")]
  
  # Do the same for bg
  bg_df$AUC <- apply(bg_df[ , threshold_colnames], 1, sum) * 0.01
  bg_df <- bg_df[ , c("spe", "cell_types", "AUC")]
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    
    # Subset for current cell type of interest combintation
    plot_df <- prevalence_df[prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Combine spes_table and SAC_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Create a 'key' column which groups simulations if they have the same bg_type, shape and size (but not arrangement)
    plot_df$key <- paste(plot_df$bg_type, plot_df$shape, plot_df$size, sep = "_")
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shape <- factor(plot_df$shape, c("Sphere", "Ellipsoid", "Network"))
    plot_df$size <- factor(plot_df$size, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    # Get plot_bg_df from bg_df
    plot_bg_df <- data.frame(arrangement = arrangements,
                             AUC = bg_df[bg_df$cell_types == entropy_cell_types$cell_types[i], "AUC"],
                             key = "bg")
    
    fig_bg_type <- ggplot(plot_df, aes(arrangement, AUC, group = key, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("AUC") +
      geom_line(data = plot_bg_df, aes(arrangement, AUC), col = "black", linetype = 2)
    
    fig_shape <- ggplot(plot_df, aes(arrangement, AUC, group = key, col = shape)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("AUC") +
      geom_line(data = plot_bg_df, aes(arrangement, AUC), col = "black", linetype = 2)
    
    fig_size <- ggplot(plot_df, aes(arrangement, AUC, group = key, col = size)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("AUC") +
      geom_line(data = plot_bg_df, aes(arrangement, AUC), col = "black", linetype = 2)
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(bg_type = fig_bg_type, shape = fig_shape, size = fig_size)
  }
  
  # Combine the plots together by cell types of interest
  library(cowplot)
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$bg_type, all_plots_list[[cell_types]]$shape, all_plots_list[[cell_types]]$size, nrow = 1, ncol = 3)
    
    title <- ggdraw() + 
      draw_label(paste("Cell types of interest:", cell_types), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_cell_types_list[[cell_types]] <- fig
  }
  
  # Combine the combined plots into one big plot
  AUC_plot <- plot_grid(plots_cell_types_list[[entropy_cell_types$cell_types[1]]], 
                        plots_cell_types_list[[entropy_cell_types$cell_types[2]]],
                        nrow = 2, ncol = 1, scale = 0.9)
  
  methods::show(AUC_plot)
  
  return(AUC_plot)
}


### 2.1. mixed spes APD ------------------------------------------------------
### 2.2. mixed spes AMD ------------------------------------------------------

# Read mixed_spes_table
setwd("~/Objects/supervised/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Read mixed_AMD_df
setwd("~/Objects/supervised/mixed_spes/analysis_3D")
mixed_AMD_df <- read.table("mixed_AMD_df.csv")

# Read bg_AMD_df
setwd("~/Objects/supervised/background_spe")
bg_AMD_df <- read.table("bg_AMD_df.csv")

mixed_AMD_plot <- plot_AMD_metric(mixed_spes_table, mixed_AMD_df, bg_AMD_df, mixed_arrangements)

setwd("~/Objects/mixed_spes/analysis_3D/plots")
# saveRDS(mixed_AMD_plot, "mixed_AMD_plot.rds")




### 2.3. mixed spes MS, NMS, ACINP, AE -----------------------------------

# Read mixed_spes_table
setwd("~/Objects/supervised/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Read mixed MS, NMS, ACINP, AE dfs
setwd("~/Objects/supervised/mixed_spes/analysis_3D")
mixed_MS_df <- read.table("mixed_MS_df.csv")
mixed_NMS_df <- read.table("mixed_NMS_df.csv")
mixed_ACINP_df <- read.table("mixed_ACINP_df.csv")
mixed_AE_df <- read.table("mixed_AE_df.csv")

# Read bg_dfs
setwd("~/Objects/supervised/background_spe")
bg_MS_df <- read.table("bg_MS_df.csv")
bg_NMS_df <- read.table("bg_NMS_df.csv")
bg_ACINP_df <- read.table("bg_ACINP_df.csv")
bg_AE_df <- read.table("bg_AE_df.csv")

mixed_MS_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_MS_df, bg_MS_df, "MS", mixed_arrangements)
mixed_NMS_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_NMS_df, bg_NMS_df, "NMS", mixed_arrangements)
mixed_ACINP_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_ACINP_df, bg_ACINP_df, "ACINP", mixed_arrangements)
mixed_AE_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_AE_df, bg_AE_df, "AE", mixed_arrangements)

setwd("~/Objects/mixed_spes/analysis_3D/plots")
# saveRDS()

### 2.4. mixed spes ACIN, CKR ------------------------------------------------

# Read mixed_spes_table
setwd("~/Objects/supervised/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Read mixed ACIN, CKR
setwd("~/Objects/supervised/mixed_spes/analysis_3D")
mixed_ACIN_df <- read.table("mixed_ACIN_df.csv")
mixed_CKR_df <- read.table("mixed_CKR_df.csv")

# Read bg_dfs
setwd("~/Objects/supervised/background_spe")
bg_ACIN_df <- read.table("bg_ACIN_df.csv")
bg_CKR_df <- read.table("bg_CKR_df.csv")

# Get plots
mixed_ACIN_plot <- plot_gradient_metrics_type2(mixed_spes_table, mixed_ACIN_df, bg_ACIN_df, "ACIN", mixed_arrangements, 0, 50)

mixed_CKR_plot <- plot_gradient_metrics_type2(mixed_spes_table, mixed_CKR_df, bg_CKR_df, "CKR", mixed_arrangements, 15, 50)


### 2.5. mixed spes SAC ------------------------------------------------------

# Read mixed_spes_table
setwd("~/Objects/supervised/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Read mixed_SAC_df
setwd("~/Objects/supervised/mixed_spes/analysis_3D")
mixed_prop_SAC_df <- read.table("mixed_prop_SAC_df.csv")
mixed_entropy_SAC_df <- read.table("mixed_entropy_SAC_df.csv")

# Read bg_SAC_df
setwd("~/Objects/supervised/background_spe")
bg_prop_SAC_df <- read.table("bg_prop_SAC_df.csv")
bg_entropy_SAC_df <- read.table("bg_entropy_SAC_df.csv")

mixed_prop_SAC_plot <- plot_proportion_SAC(mixed_spes_table, mixed_prop_SAC_df, bg_prop_SAC_df, mixed_arrangements)
mixed_entropy_SAC_plot <- plot_entropy_SAC(mixed_spes_table, mixed_entropy_SAC_df, bg_entropy_SAC_df, mixed_arrangements)

# setwd("~/Objects/mixed_spes/analysis_3D/plots")
# saveRDS(mixed_SAC_plot, "mixed_SAC_plot.rds")

### 2.6. mixed spes prevalence ------------------------------------------------

# Read mixed_spes_table
setwd("~/Objects/supervised/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Read mixed prevalence dfs
setwd("~/Objects/supervised/mixed_spes/analysis_3D")
mixed_prop_prevalence_df <- read.table("mixed_prop_prevalence_df.csv")
mixed_entropy_prevalence_df <- read.table("mixed_entropy_prevalence_df.csv")

# Read bg prevalence dfs
setwd("~/Objects/supervised/background_spe")
bg_prop_prevalence_df <- read.table("bg_prop_prevalence_df.csv")
bg_entropy_prevalence_df <- read.table("bg_entropy_prevalence_df.csv")

mixed_prop_prevalence_plot <- plot_proportion_prevalence(mixed_spes_table, mixed_prop_prevalence_df, bg_prop_prevalence_df, mixed_arrangements)
mixed_entropy_prevalence_plot <- plot_entropy_prevalence(mixed_spes_table, mixed_entropy_prevalence_df, bg_entropy_prevalence_df, mixed_arrangements)

mixed_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(mixed_spes_table, mixed_prop_prevalence_df, bg_prop_prevalence_df, mixed_arrangements)
mixed_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(mixed_spes_table, mixed_entropy_prevalence_df, bg_entropy_prevalence_df, mixed_arrangements)

setwd("~/Objects/mixed_spes/analysis_3D/plots")
# saveRDS(mixed_prevalence_plot, "mixed_prevalence_plot.rds")



### 3.1. ringed spes APD ------------------------------------------------------
### 3.2. ringed spes AMD ------------------------------------------------------

# Read ringed_spes_table
setwd("~/Objects/supervised/spes_table")
ringed_spes_table <- read.table("ringed_spes_table.csv")

# Read ringed_AMD_df
setwd("~/Objects/supervised/ringed_spes/analysis_3D")
ringed_AMD_df <- read.table("ringed_AMD_df.csv")

# Read bg_AMD_df
setwd("~/Objects/supervised/background_spe")
bg_AMD_df <- read.table("bg_AMD_df.csv")

ringed_AMD_plot <- plot_AMD_metric(ringed_spes_table, ringed_AMD_df, bg_AMD_df, ringed_arrangements)

# setwd("~/Objects/supervised/ringed_spes/analysis_3D/plots")
# saveRDS(ringed_AMD_plot, "ringed_AMD_plot.rds")




### 3.3. ringed spes MS, NMS, ACINP, AE -----------------------------------

# Read ringed_spes_table
setwd("~/Objects/supervised/spes_table")
ringed_spes_table <- read.table("ringed_spes_table.csv")

# Read ringed MS, NMS, ACINP, AE dfs
setwd("~/Objects/supervised/ringed_spes/analysis_3D")
ringed_MS_df <- read.table("ringed_MS_df.csv")
ringed_NMS_df <- read.table("ringed_NMS_df.csv")
ringed_ACINP_df <- read.table("ringed_ACINP_df.csv")
ringed_AE_df <- read.table("ringed_AE_df.csv")

# Read bg_dfs
setwd("~/Objects/supervised/background_spe")
bg_MS_df <- read.table("bg_MS_df.csv")
bg_NMS_df <- read.table("bg_NMS_df.csv")
bg_ACINP_df <- read.table("bg_ACINP_df.csv")
bg_AE_df <- read.table("bg_AE_df.csv")

ringed_MS_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_MS_df, bg_MS_df, "MS", ringed_arrangements)
ringed_NMS_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_NMS_df, bg_NMS_df, "NMS", ringed_arrangements)
ringed_ACINP_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_ACINP_df, bg_ACINP_df, "ACINP", ringed_arrangements)
ringed_AE_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_AE_df, bg_AE_df, "AE", ringed_arrangements)

setwd("~/Objects/ringed_spes/analysis_3D/plots")
# saveRDS()

### 3.4. ringed spes ACIN, CKR ------------------------------------------------

# Read ringed_spes_table
setwd("~/Objects/supervised/spes_table")
ringed_spes_table <- read.table("ringed_spes_table.csv")

# Read ringed ACIN, CKR
setwd("~/Objects/supervised/ringed_spes/analysis_3D")
ringed_ACIN_df <- read.table("ringed_ACIN_df.csv")
ringed_CKR_df <- read.table("ringed_CKR_df.csv")

# Read bg_dfs
setwd("~/Objects/supervised/background_spe")
bg_ACIN_df <- read.table("bg_ACIN_df.csv")
bg_CKR_df <- read.table("bg_CKR_df.csv")

# Get plots
ringed_ACIN_plot <- plot_gradient_metrics_type2(ringed_spes_table, ringed_ACIN_df, bg_ACIN_df, "ACIN", ringed_arrangements, 0, 50)

ringed_CKR_plot <- plot_gradient_metrics_type2(ringed_spes_table, ringed_CKR_df, bg_CKR_df, "CKR", ringed_arrangements, 15, 50)


### 3.5. ringed spes SAC ------------------------------------------------------

# Read ringed_spes_table
setwd("~/Objects/supervised/spes_table")
ringed_spes_table <- read.table("ringed_spes_table.csv")

# Read ringed_SAC_df
setwd("~/Objects/supervised/ringed_spes/analysis_3D")
ringed_prop_SAC_df <- read.table("ringed_prop_SAC_df.csv")
ringed_entropy_SAC_df <- read.table("ringed_entropy_SAC_df.csv")

# Read bg_SAC_df
setwd("~/Objects/supervised/background_spe")
bg_prop_SAC_df <- read.table("bg_prop_SAC_df.csv")
bg_entropy_SAC_df <- read.table("bg_entropy_SAC_df.csv")

ringed_prop_SAC_plot <- plot_proportion_SAC(ringed_spes_table, ringed_prop_SAC_df, bg_prop_SAC_df, ringed_arrangements)
ringed_entropy_SAC_plot <- plot_entropy_SAC(ringed_spes_table, ringed_entropy_SAC_df, bg_entropy_SAC_df, ringed_arrangements)

# setwd("~/Objects/ringed_spes/analysis_3D/plots")
# saveRDS(ringed_SAC_plot, "ringed_SAC_plot.rds")

### 3.6. ringed spes prevalence ------------------------------------------------

# Read ringed_spes_table
setwd("~/Objects/supervised/spes_table")
ringed_spes_table <- read.table("ringed_spes_table.csv")

# Read ringed prevalence dfs
setwd("~/Objects/supervised/ringed_spes/analysis_3D")
ringed_prop_prevalence_df <- read.table("ringed_prop_prevalence_df.csv")
ringed_entropy_prevalence_df <- read.table("ringed_entropy_prevalence_df.csv")

# Read bg prevalence dfs
setwd("~/Objects/supervised/background_spe")
bg_prop_prevalence_df <- read.table("bg_prop_prevalence_df.csv")
bg_entropy_prevalence_df <- read.table("bg_entropy_prevalence_df.csv")

ringed_prop_prevalence_plot <- plot_proportion_prevalence(ringed_spes_table, ringed_prop_prevalence_df, bg_prop_prevalence_df, ringed_arrangements)
ringed_entropy_prevalence_plot <- plot_entropy_prevalence(ringed_spes_table, ringed_entropy_prevalence_df, bg_entropy_prevalence_df, ringed_arrangements)

ringed_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(ringed_spes_table, ringed_prop_prevalence_df, bg_prop_prevalence_df, ringed_arrangements)
ringed_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(ringed_spes_table, ringed_entropy_prevalence_df, bg_entropy_prevalence_df, ringed_arrangements)

setwd("~/Objects/ringed_spes/analysis_3D/plots")
# saveRDS(ringed_prevalence_plot, "ringed_prevalence_plot.rds")



### 4.1. all spes APD ------------------------------------------------------
### 4.2. all spes AMD ------------------------------------------------------

# Read spes_table
setwd("~/Objects/supervised/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")
ringed_spes_table <- read.table("ringed_spes_table.csv")
separated_spes_table <- read.table("separated_spes_table.csv")

mixed_spes_table$spe <- paste("mixed_spe_", seq(nrow(mixed_spes_table)), sep = "")
ringed_spes_table$spe <- paste("ringed_spe_", seq(nrow(ringed_spes_table)), sep = "")
separated_spes_table$spe <- paste("separated_spe_", seq(nrow(separated_spes_table)), sep = "")

# Subset separated_spes_table for when the shapes are the same, and the sizes are the same
# Subset all spes_table for "M2", "R2" and "S2" arrangement
mixed_spes_table <- mixed_spes_table[mixed_spes_table$arrangement == "M2", ]
ringed_spes_table <- ringed_spes_table[ringed_spes_table$arrangement == "R2", ]
separated_spes_table <- separated_spes_table[separated_spes_table$shapeA == separated_spes_table$shapeB & 
                                               separated_spes_table$sizeA == separated_spes_table$sizeB &
                                               separated_spes_table$arrangement == "S2", c("bg_type", "shapeA", "sizeA", "arrangement", "spe")]
colnames(separated_spes_table) <- c("bg_type", "shape", "size", "arrangement", "spe")

# Combine all spes tables together
all_spes_table <- rbind(mixed_spes_table, ringed_spes_table, separated_spes_table)


# Read AMD dfs
setwd("~/Objects/supervised/mixed_spes/analysis_3D")
mixed_AMD_df <- read.table("mixed_AMD_df.csv")
setwd("~/Objects/supervised/ringed_spes/analysis_3D")
ringed_AMD_df <- read.table("ringed_AMD_df.csv")
setwd("~/Objects/supervised/separated_spes/analysis_3D")
separated_AMD_df <- read.table("separated_AMD_df.csv")

# Subset all AMD df to align with spes_table
mixed_AMD_df <- mixed_AMD_df[mixed_AMD_df$spe %in% mixed_spes_table$spe, ]
ringed_AMD_df <- ringed_AMD_df[ringed_AMD_df$spe %in% ringed_spes_table$spe, ]
separated_AMD_df <- separated_AMD_df[separated_AMD_df$spe %in% separated_spes_table$spe, ]

# Remove spe column in spes_table
all_spes_table$spe <- NULL

# Combine all AMD dfs together
all_AMD_df <- rbind(mixed_AMD_df, ringed_AMD_df, separated_AMD_df)

# Read bg_AMD_df
setwd("~/Objects/supervised/background_spe")
bg_AMD_df <- read.table("bg_AMD_df.csv")

all_AMD_plot <- plot_AMD_metric(all_spes_table, all_AMD_df, bg_AMD_df, all_arrangements)

# setwd("~/Objects/supervised/ringed_spes/analysis_3D/plots")
# saveRDS(ringed_AMD_plot, "ringed_AMD_plot.rds")




### 4.3. all spes MS, NMS, ACINP, AE -----------------------------------

# Read spes_table
setwd("~/Objects/supervised/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")
ringed_spes_table <- read.table("ringed_spes_table.csv")
separated_spes_table <- read.table("separated_spes_table.csv")

mixed_spes_table$spe <- paste("mixed_spe_", seq(nrow(mixed_spes_table)), sep = "")
ringed_spes_table$spe <- paste("ringed_spe_", seq(nrow(ringed_spes_table)), sep = "")
separated_spes_table$spe <- paste("separated_spe_", seq(nrow(separated_spes_table)), sep = "")

# Subset separated_spes_table for when the shapes are the same, and the sizes are the same
# Subset all spes_table for "M2", "R2" and "S2" arrangement
mixed_spes_table <- mixed_spes_table[mixed_spes_table$arrangement == "M2", ]
ringed_spes_table <- ringed_spes_table[ringed_spes_table$arrangement == "R2", ]
separated_spes_table <- separated_spes_table[separated_spes_table$shapeA == separated_spes_table$shapeB & 
                                               separated_spes_table$sizeA == separated_spes_table$sizeB &
                                               separated_spes_table$arrangement == "S2", c("bg_type", "shapeA", "sizeA", "arrangement", "spe")]
colnames(separated_spes_table) <- c("bg_type", "shape", "size", "arrangement", "spe")

# Combine all spes tables together
all_spes_table <- rbind(mixed_spes_table, ringed_spes_table, separated_spes_table)



# Read ringed MS, NMS, ACINP, AE dfs
setwd("~/Objects/supervised/mixed_spes/analysis_3D")
mixed_MS_df <- read.table("mixed_MS_df.csv")
mixed_NMS_df <- read.table("mixed_NMS_df.csv")
mixed_ACINP_df <- read.table("mixed_ACINP_df.csv")
mixed_AE_df <- read.table("mixed_AE_df.csv")

setwd("~/Objects/supervised/ringed_spes/analysis_3D")
ringed_MS_df <- read.table("ringed_MS_df.csv")
ringed_NMS_df <- read.table("ringed_NMS_df.csv")
ringed_ACINP_df <- read.table("ringed_ACINP_df.csv")
ringed_AE_df <- read.table("ringed_AE_df.csv")

setwd("~/Objects/supervised/separated_spes/analysis_3D")
separated_MS_df <- read.table("separated_MS_df.csv")
separated_NMS_df <- read.table("separated_NMS_df.csv")
separated_ACINP_df <- read.table("separated_ACINP_df.csv")
separated_AE_df <- read.table("separated_AE_df.csv")



# Subset gradient_df to allign with spes_table
mixed_MS_df <- mixed_MS_df[mixed_MS_df$spe %in% mixed_spes_table$spe, ]
mixed_NMS_df <- mixed_NMS_df[mixed_NMS_df$spe %in% mixed_spes_table$spe, ]
mixed_ACINP_df <- mixed_ACINP_df[mixed_ACINP_df$spe %in% mixed_spes_table$spe, ]
mixed_AE_df <- mixed_AE_df[mixed_AE_df$spe %in% mixed_spes_table$spe, ]

ringed_MS_df <- ringed_MS_df[ringed_MS_df$spe %in% ringed_spes_table$spe, ]
ringed_NMS_df <- ringed_NMS_df[ringed_NMS_df$spe %in% ringed_spes_table$spe, ]
ringed_ACINP_df <- ringed_ACINP_df[ringed_ACINP_df$spe %in% ringed_spes_table$spe, ]
ringed_AE_df <- ringed_AE_df[ringed_AE_df$spe %in% ringed_spes_table$spe, ]

separated_MS_df <- separated_MS_df[separated_MS_df$spe %in% separated_spes_table$spe, ]
separated_NMS_df <- separated_NMS_df[separated_NMS_df$spe %in% separated_spes_table$spe, ]
separated_ACINP_df <- separated_ACINP_df[separated_ACINP_df$spe %in% separated_spes_table$spe, ]
separated_AE_df <- separated_AE_df[separated_AE_df$spe %in% separated_spes_table$spe, ]

# Remove spe column in spes_table
all_spes_table$spe <- NULL

# Combine all MS, NMS, ACINP, AE dfs together
all_MS_df <- rbind(mixed_MS_df, ringed_MS_df, separated_MS_df)
all_NMS_df <- rbind(mixed_NMS_df, ringed_NMS_df, separated_NMS_df)
all_ACINP_df <- rbind(mixed_ACINP_df, ringed_ACINP_df, separated_ACINP_df)
all_AE_df <- rbind(mixed_AE_df, ringed_AE_df, separated_AE_df)


# Read bg_dfs
setwd("~/Objects/supervised/background_spe")
bg_MS_df <- read.table("bg_MS_df.csv")
bg_NMS_df <- read.table("bg_NMS_df.csv")
bg_ACINP_df <- read.table("bg_ACINP_df.csv")
bg_AE_df <- read.table("bg_AE_df.csv")

all_MS_plot <- plot_gradient_metrics_type1(all_spes_table, all_MS_df, bg_MS_df, "MS", all_arrangements)
all_NMS_plot <- plot_gradient_metrics_type1(all_spes_table, all_NMS_df, bg_NMS_df, "NMS", all_arrangements)
all_ACINP_plot <- plot_gradient_metrics_type1(all_spes_table, all_ACINP_df, bg_ACINP_df, "ACINP", all_arrangements)
all_AE_plot <- plot_gradient_metrics_type1(all_spes_table, all_AE_df, bg_AE_df, "AE", all_arrangements)

setwd("~/Objects/supervised/ringed_spes/analysis_3D/plots")
# saveRDS()

### 4.4. all spes ACIN, CKR ------------------------------------------------

# Read spes_table
setwd("~/Objects/supervised/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")
ringed_spes_table <- read.table("ringed_spes_table.csv")
separated_spes_table <- read.table("separated_spes_table.csv")

mixed_spes_table$spe <- paste("mixed_spe_", seq(nrow(mixed_spes_table)), sep = "")
ringed_spes_table$spe <- paste("ringed_spe_", seq(nrow(ringed_spes_table)), sep = "")
separated_spes_table$spe <- paste("separated_spe_", seq(nrow(separated_spes_table)), sep = "")

# Subset separated_spes_table for when the shapes are the same, and the sizes are the same
# Subset all spes_table for "M2", "R2" and "S2" arrangement
mixed_spes_table <- mixed_spes_table[mixed_spes_table$arrangement == "M2", ]
ringed_spes_table <- ringed_spes_table[ringed_spes_table$arrangement == "R2", ]
separated_spes_table <- separated_spes_table[separated_spes_table$shapeA == separated_spes_table$shapeB & 
                                               separated_spes_table$sizeA == separated_spes_table$sizeB &
                                               separated_spes_table$arrangement == "S2", c("bg_type", "shapeA", "sizeA", "arrangement", "spe")]
colnames(separated_spes_table) <- c("bg_type", "shape", "size", "arrangement", "spe")

# Combine all spes tables together
all_spes_table <- rbind(mixed_spes_table, ringed_spes_table, separated_spes_table)


# Read ringed CKR and ACIN
setwd("~/Objects/supervised/mixed_spes/analysis_3D")
mixed_ACIN_df <- read.table("mixed_ACIN_df.csv")
mixed_CKR_df <- read.table("mixed_CKR_df.csv")

setwd("~/Objects/supervised/ringed_spes/analysis_3D")
ringed_ACIN_df <- read.table("ringed_ACIN_df.csv")
ringed_CKR_df <- read.table("ringed_CKR_df.csv")

setwd("~/Objects/supervised/separated_spes/analysis_3D")
separated_ACIN_df <- read.table("separated_ACIN_df.csv")
separated_CKR_df <- read.table("separated_CKR_df.csv")



# Subset gradient_df to allign with spes_table
mixed_ACIN_df <- mixed_ACIN_df[mixed_ACIN_df$spe %in% mixed_spes_table$spe, ]
mixed_CKR_df <- mixed_CKR_df[mixed_CKR_df$spe %in% mixed_spes_table$spe, ]

ringed_ACIN_df <- ringed_ACIN_df[ringed_ACIN_df$spe %in% ringed_spes_table$spe, ]
ringed_CKR_df <- ringed_CKR_df[ringed_CKR_df$spe %in% ringed_spes_table$spe, ]

separated_ACIN_df <- separated_ACIN_df[separated_ACIN_df$spe %in% separated_spes_table$spe, ]
separated_CKR_df <- separated_CKR_df[separated_CKR_df$spe %in% separated_spes_table$spe, ]

# Remove spe column in spes_table
all_spes_table$spe <- NULL

# Combine all MS, NMS, ACINP, AE dfs together
all_ACIN_df <- rbind(mixed_ACIN_df, ringed_ACIN_df, separated_ACIN_df)
all_CKR_df <- rbind(mixed_CKR_df, ringed_CKR_df, separated_CKR_df)

# Read bg_dfs
setwd("~/Objects/supervised/background_spe")
bg_ACIN_df <- read.table("bg_ACIN_df.csv")
bg_CKR_df <- read.table("bg_CKR_df.csv")

# Get plots
all_ACIN_plot <- plot_gradient_metrics_type2(all_spes_table, all_ACIN_df, bg_ACIN_df, "ACIN", all_arrangements, 0, 50)

all_CKR_plot <- plot_gradient_metrics_type2(all_spes_table, all_CKR_df, bg_CKR_df, "CKR", all_arrangements, 15, 50)


### 4.5. all spes SAC ------------------------------------------------------

# Read spes_table
setwd("~/Objects/supervised/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")
ringed_spes_table <- read.table("ringed_spes_table.csv")
separated_spes_table <- read.table("separated_spes_table.csv")

mixed_spes_table$spe <- paste("mixed_spe_", seq(nrow(mixed_spes_table)), sep = "")
ringed_spes_table$spe <- paste("ringed_spe_", seq(nrow(ringed_spes_table)), sep = "")
separated_spes_table$spe <- paste("separated_spe_", seq(nrow(separated_spes_table)), sep = "")

# Subset separated_spes_table for when the shapes are the same, and the sizes are the same
# Subset all spes_table for "M2", "R2" and "S2" arrangement
mixed_spes_table <- mixed_spes_table[mixed_spes_table$arrangement == "M2", ]
ringed_spes_table <- ringed_spes_table[ringed_spes_table$arrangement == "R2", ]
separated_spes_table <- separated_spes_table[separated_spes_table$shapeA == separated_spes_table$shapeB & 
                                               separated_spes_table$sizeA == separated_spes_table$sizeB &
                                               separated_spes_table$arrangement == "S2", c("bg_type", "shapeA", "sizeA", "arrangement", "spe")]
colnames(separated_spes_table) <- c("bg_type", "shape", "size", "arrangement", "spe")

# Combine all spes tables together
all_spes_table <- rbind(mixed_spes_table, ringed_spes_table, separated_spes_table)



# Read SAC dfs
setwd("~/Objects/supervised/mixed_spes/analysis_3D")
mixed_prop_SAC_df <- read.table("mixed_prop_SAC_df.csv")
mixed_entropy_SAC_df <- read.table("mixed_entropy_SAC_df.csv")
setwd("~/Objects/supervised/ringed_spes/analysis_3D")
ringed_prop_SAC_df <- read.table("ringed_prop_SAC_df.csv")
ringed_entropy_SAC_df <- read.table("ringed_entropy_SAC_df.csv")
setwd("~/Objects/supervised/separated_spes/analysis_3D")
separated_prop_SAC_df <- read.table("separated_prop_SAC_df.csv")
separated_entropy_SAC_df <- read.table("separated_entropy_SAC_df.csv")

# Subset all prop_SAC df to align with spes_table
mixed_prop_SAC_df <- mixed_prop_SAC_df[mixed_prop_SAC_df$spe %in% mixed_spes_table$spe, ]
ringed_prop_SAC_df <- ringed_prop_SAC_df[ringed_prop_SAC_df$spe %in% ringed_spes_table$spe, ]
separated_prop_SAC_df <- separated_prop_SAC_df[separated_prop_SAC_df$spe %in% separated_spes_table$spe, ]
mixed_entropy_SAC_df <- mixed_entropy_SAC_df[mixed_entropy_SAC_df$spe %in% mixed_spes_table$spe, ]
ringed_entropy_SAC_df <- ringed_entropy_SAC_df[ringed_entropy_SAC_df$spe %in% ringed_spes_table$spe, ]
separated_entropy_SAC_df <- separated_entropy_SAC_df[separated_entropy_SAC_df$spe %in% separated_spes_table$spe, ]

# Remove spe column in spes_table
all_spes_table$spe <- NULL

# Combine all SAC dfs together
all_prop_SAC_df <- rbind(mixed_prop_SAC_df, ringed_prop_SAC_df, separated_prop_SAC_df)
all_entropy_SAC_df <- rbind(mixed_entropy_SAC_df, ringed_entropy_SAC_df, separated_entropy_SAC_df)

# Read bg_SAC_df
setwd("~/Objects/supervised/background_spe")
bg_prop_SAC_df <- read.table("bg_prop_SAC_df.csv")
bg_entropy_SAC_df <- read.table("bg_entropy_SAC_df.csv")

all_prop_SAC_plot <- plot_proportion_SAC(all_spes_table, all_prop_SAC_df, bg_prop_SAC_df, all_arrangements)
all_entropy_SAC_plot <- plot_entropy_SAC(all_spes_table, all_entropy_SAC_df, bg_entropy_SAC_df, all_arrangements)

# setwd("~/Objects/ringed_spes/analysis_3D/plots")
# saveRDS(ringed_SAC_plot, "ringed_SAC_plot.rds")

### 4.6. all spes prevalence ------------------------------------------------

# Read spes_table
setwd("~/Objects/supervised/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")
ringed_spes_table <- read.table("ringed_spes_table.csv")
separated_spes_table <- read.table("separated_spes_table.csv")

mixed_spes_table$spe <- paste("mixed_spe_", seq(nrow(mixed_spes_table)), sep = "")
ringed_spes_table$spe <- paste("ringed_spe_", seq(nrow(ringed_spes_table)), sep = "")
separated_spes_table$spe <- paste("separated_spe_", seq(nrow(separated_spes_table)), sep = "")

# Subset separated_spes_table for when the shapes are the same, and the sizes are the same
# Subset all spes_table for "M2", "R2" and "S2" arrangement
mixed_spes_table <- mixed_spes_table[mixed_spes_table$arrangement == "M2", ]
ringed_spes_table <- ringed_spes_table[ringed_spes_table$arrangement == "R2", ]
separated_spes_table <- separated_spes_table[separated_spes_table$shapeA == separated_spes_table$shapeB & 
                                               separated_spes_table$sizeA == separated_spes_table$sizeB &
                                               separated_spes_table$arrangement == "S2", c("bg_type", "shapeA", "sizeA", "arrangement", "spe")]
colnames(separated_spes_table) <- c("bg_type", "shape", "size", "arrangement", "spe")

# Combine all spes tables together
all_spes_table <- rbind(mixed_spes_table, ringed_spes_table, separated_spes_table)



# Read prevalence dfs
setwd("~/Objects/supervised/mixed_spes/analysis_3D")
mixed_prop_prevalence_df <- read.table("mixed_prop_prevalence_df.csv")
mixed_entropy_prevalence_df <- read.table("mixed_entropy_prevalence_df.csv")
setwd("~/Objects/supervised/ringed_spes/analysis_3D")
ringed_prop_prevalence_df <- read.table("ringed_prop_prevalence_df.csv")
ringed_entropy_prevalence_df <- read.table("ringed_entropy_prevalence_df.csv")
setwd("~/Objects/supervised/separated_spes/analysis_3D")
separated_prop_prevalence_df <- read.table("separated_prop_prevalence_df.csv")
separated_entropy_prevalence_df <- read.table("separated_entropy_prevalence_df.csv")

# Subset all prop_prevalence df to align with spes_table
mixed_prop_prevalence_df <- mixed_prop_prevalence_df[mixed_prop_prevalence_df$spe %in% mixed_spes_table$spe, ]
ringed_prop_prevalence_df <- ringed_prop_prevalence_df[ringed_prop_prevalence_df$spe %in% ringed_spes_table$spe, ]
separated_prop_prevalence_df <- separated_prop_prevalence_df[separated_prop_prevalence_df$spe %in% separated_spes_table$spe, ]
mixed_entropy_prevalence_df <- mixed_entropy_prevalence_df[mixed_entropy_prevalence_df$spe %in% mixed_spes_table$spe, ]
ringed_entropy_prevalence_df <- ringed_entropy_prevalence_df[ringed_entropy_prevalence_df$spe %in% ringed_spes_table$spe, ]
separated_entropy_prevalence_df <- separated_entropy_prevalence_df[separated_entropy_prevalence_df$spe %in% separated_spes_table$spe, ]

# Remove spe column in spes_table
all_spes_table$spe <- NULL

# Combine all prevalence dfs together
all_prop_prevalence_df <- rbind(mixed_prop_prevalence_df, ringed_prop_prevalence_df, separated_prop_prevalence_df)
all_entropy_prevalence_df <- rbind(mixed_entropy_prevalence_df, ringed_entropy_prevalence_df, separated_entropy_prevalence_df)



# Read bg prevalence dfs
setwd("~/Objects/supervised/background_spe")
bg_prop_prevalence_df <- read.table("bg_prop_prevalence_df.csv")
bg_entropy_prevalence_df <- read.table("bg_entropy_prevalence_df.csv")

all_prop_prevalence_plot <- plot_proportion_prevalence(all_spes_table, all_prop_prevalence_df, bg_prop_prevalence_df, all_arrangements)
all_entropy_prevalence_plot <- plot_entropy_prevalence(all_spes_table, all_entropy_prevalence_df, bg_entropy_prevalence_df, all_arrangements)

all_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(all_spes_table, all_prop_prevalence_df, bg_prop_prevalence_df, all_arrangements)
all_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(all_spes_table, all_entropy_prevalence_df, bg_entropy_prevalence_df, all_arrangements)

setwd("~/Objects/ringed_spes/analysis_3D/plots")
# saveRDS(ringed_prevalence_plot, "ringed_prevalence_plot.rds")


