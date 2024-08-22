library(cowplot)
library(ggplot2)
separated_arrangements <- c("S1", "S2", "S3")

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
    plot_df$key <- paste(plot_df$bg_type, plot_df$shapeA, plot_df$sizeA, plot_df$shapeB, plot_df$sizeB, sep = "_")
    
    # Get plot_bg_df from bg_AMD_df
    plot_bg_df <- data.frame(arrangement = arrangements, 
                             AMD = bg_AMD_df[bg_AMD_df$reference == AMD_pairs[i, "cell1"] & bg_AMD_df$target == AMD_pairs[i, "cell2"], "AMD"],
                             key = "bg")
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shapeA <- factor(plot_df$shapeA, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeA <- factor(plot_df$sizeA, c("Small", "Medium", "Large"))
    plot_df$shapeB <- factor(plot_df$shapeB, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeB <- factor(plot_df$sizeB, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    fig_bg_type <- ggplot(plot_df, aes(arrangement, AMD, group = key, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_line(data = plot_bg_df, aes(arrangement, AMD), col = "black", linetype = 2)
    
    fig_shapeA <- ggplot(plot_df, aes(arrangement, AMD, group = key, col = shapeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_line(data = plot_bg_df, aes(arrangement, AMD), col = "black", linetype = 2)
    
    fig_sizeA <- ggplot(plot_df, aes(arrangement, AMD, group = key, col = sizeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_line(data = plot_bg_df, aes(arrangement, AMD), col = "black", linetype = 2)
    
    fig_shapeB <- ggplot(plot_df, aes(arrangement, AMD, group = key, col = shapeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_line(data = plot_bg_df, aes(arrangement, AMD), col = "black", linetype = 2)
    
    fig_sizeB <- ggplot(plot_df, aes(arrangement, AMD, group = key, col = sizeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      geom_line(data = plot_bg_df, aes(arrangement, AMD), col = "black", linetype = 2)
    
    all_plots_list[[AMD_pairs[i, "pair"]]] <- list(bg_type = fig_bg_type, 
                                                   shapeA = fig_shapeA, sizeA = fig_sizeA,
                                                   shapeB = fig_shapeB, sizeB = fig_sizeB)
  }
  
  
  # Combine the plots together by pairs
  library(cowplot)
  plots_pair_list <- list()
  
  for (i in seq(nrow(AMD_pairs))) {
    pair <- AMD_pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$bg_type, 
                       all_plots_list[[pair]]$shapeA, all_plots_list[[pair]]$sizeA, 
                       all_plots_list[[pair]]$shapeB, all_plots_list[[pair]]$sizeB, 
                       nrow = 1, ncol = 5)
    
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
    plot_df$shapeA <- factor(plot_df$shapeA, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeA <- factor(plot_df$sizeA, c("Small", "Medium", "Large"))
    plot_df$shapeB <- factor(plot_df$shapeB, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeB <- factor(plot_df$sizeB, c("Small", "Medium", "Large"))
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
    
    fig_shapeA <- ggplot(plot_df, aes(variable, value, group = spe, col = shapeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_sizeA <- ggplot(plot_df, aes(variable, value, group = spe, col = sizeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2)
    
    fig_shapeB <- ggplot(plot_df, aes(variable, value, group = spe, col = shapeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_sizeB <- ggplot(plot_df, aes(variable, value, group = spe, col = sizeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2)
    
    all_plots_list[[cell_types[i]]] <- list(arrangement = fig_arrangement, bg_type = fig_bg_type, 
                                            shapeA = fig_shapeA, sizeA = fig_sizeA,
                                            shapeB = fig_shapeB, sizeB = fig_sizeB)
    
  }
  
  
  # Combine the plots together by reference cell type
  plots_ref_list <- list()
  
  for (i in seq(length(cell_types))) {
    reference_cell_type <- cell_types[i]
    
    plots <- plot_grid(all_plots_list[[reference_cell_type]]$arrangement,
                       all_plots_list[[reference_cell_type]]$bg_type, 
                       all_plots_list[[reference_cell_type]]$shapeA, 
                       all_plots_list[[reference_cell_type]]$sizeA, 
                       all_plots_list[[reference_cell_type]]$shapeB, 
                       all_plots_list[[reference_cell_type]]$sizeB,
                       nrow = 1, ncol = 6)
    
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
    plot_df$shapeA <- factor(plot_df$shapeA, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeA <- factor(plot_df$sizeA, c("Small", "Medium", "Large"))
    plot_df$shapeB <- factor(plot_df$shapeB, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeB <- factor(plot_df$sizeB, c("Small", "Medium", "Large"))
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
    
    fig_shapeA <- ggplot(plot_df, aes(variable, value, group = spe, col = shapeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_sizeA <- ggplot(plot_df, aes(variable, value, group = spe, col = sizeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_shapeB <- ggplot(plot_df, aes(variable, value, group = spe, col = shapeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_sizeB <- ggplot(plot_df, aes(variable, value, group = spe, col = sizeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric) +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    all_plots_list[[pairs[i, "pair"]]] <- list(arrangement = fig_arrangement, bg_type = fig_bg_type, 
                                               shapeA = fig_shapeA, sizeA = fig_sizeA,
                                               shapeB = fig_shapeB, sizeB = fig_sizeB)
    
  }
  
  
  # Combine the plots together by reference cell type
  plots_pair_list <- list()
  
  for (i in seq(nrow(pairs))) {
    pair <- pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$arrangement,
                       all_plots_list[[pair]]$bg_type, 
                       all_plots_list[[pair]]$shapeA, 
                       all_plots_list[[pair]]$sizeA,
                       all_plots_list[[pair]]$shapeB, 
                       all_plots_list[[pair]]$sizeB,
                       nrow = 1, ncol = 6)
    
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
    plot_df$key <- paste(plot_df$bg_type, plot_df$shapeA, plot_df$sizeA, plot_df$shapeB, plot_df$sizeB, sep = "_")

    # Get plot_bg_df from bg_SAC_df
    plot_bg_df <- data.frame(arrangement = arrangements, 
                             proportion = bg_SAC_df[bg_SAC_df$reference == prop_cell_types$ref[i], "proportion"],
                             key = "bg")
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shapeA <- factor(plot_df$shapeA, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeA <- factor(plot_df$sizeA, c("Small", "Medium", "Large"))
    plot_df$shapeB <- factor(plot_df$shapeB, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeB <- factor(plot_df$sizeB, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    fig_bg_type <- ggplot(plot_df, aes(arrangement, proportion, group = key, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    fig_shapeA <- ggplot(plot_df, aes(arrangement, proportion, group = key, col = shapeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    fig_sizeA <- ggplot(plot_df, aes(arrangement, proportion, group = key, col = sizeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    fig_shapeB <- ggplot(plot_df, aes(arrangement, proportion, group = key, col = shapeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    fig_sizeB <- ggplot(plot_df, aes(arrangement, proportion, group = key, col = sizeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(bg_type = fig_bg_type, 
                                                      shapeA = fig_shapeA, sizeA = fig_sizeA, 
                                                      shapeB = fig_shapeB, sizeB = fig_sizeB)
  }
  
  
  # Combine the plots together by reference target pairs
  library(cowplot)
  plots_pair_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$bg_type, 
                       all_plots_list[[pair]]$shapeA, all_plots_list[[pair]]$sizeA,
                       all_plots_list[[pair]]$shapeB, all_plots_list[[pair]]$sizeB,
                       nrow = 1, ncol = 5)
    
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
    plot_df$key <- paste(plot_df$bg_type, plot_df$shape, plot_df$size, plot_df$shapeB, plot_df$sizeB, sep = "_")
    
    # Get plot_bg_df from bg_SAC_df
    plot_bg_df <- data.frame(arrangement = arrangements, 
                             proportion = bg_SAC_df[bg_SAC_df$cell_types == entropy_cell_types$cell_types[i], "entropy"],
                             key = "bg")
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shapeA <- factor(plot_df$shapeA, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeA <- factor(plot_df$sizeA, c("Small", "Medium", "Large"))
    plot_df$shapeB <- factor(plot_df$shapeB, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeB <- factor(plot_df$sizeB, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    size = 1
    alpha = 1
    
    fig_bg_type <- ggplot(plot_df, aes(arrangement, entropy, group = key, col = bg_type)) +
      geom_point(size = size, alpha = alpha, aes(shape = bg_type)) +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    fig_shapeA <- ggplot(plot_df, aes(arrangement, entropy, group = key, col = shapeA)) +
      geom_point(size = size, alpha = alpha, aes(shape = shapeA)) +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    fig_sizeA <- ggplot(plot_df, aes(arrangement, entropy, group = key, col = sizeA)) +
      geom_point(size = size, alpha = alpha, aes(shape = sizeA)) +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    fig_shapeB <- ggplot(plot_df, aes(arrangement, entropy, group = key, col = shapeB)) +
      geom_point(size = size, alpha = alpha, aes(shape = shapeB)) +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    fig_sizeB <- ggplot(plot_df, aes(arrangement, entropy, group = key, col = sizeB)) +
      geom_point(size = size, alpha = alpha, aes(shape = sizeB)) +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC") +
      geom_line(data = plot_bg_df, aes(arrangement, proportion), col = "black", linetype = 2)
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(bg_type = fig_bg_type, 
                                                               shapeA = fig_shapeA, sizeA = fig_sizeA,
                                                               shapeB = fig_shapeB, sizeB = fig_sizeB)
  }
  
  # Combine the plots together by cell types of interest
  library(cowplot)
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$bg_type, 
                       all_plots_list[[cell_types]]$shapeA, all_plots_list[[cell_types]]$sizeA,
                       all_plots_list[[cell_types]]$shapeB, all_plots_list[[cell_types]]$sizeB,
                       nrow = 1, ncol = 5)
    
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
    plot_df$shapeA <- factor(plot_df$shapeA, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeA <- factor(plot_df$sizeA, c("Small", "Medium", "Large"))
    plot_df$shapeB <- factor(plot_df$shapeB, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeB <- factor(plot_df$sizeB, c("Small", "Medium", "Large"))
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
    
    fig_shapeA <- ggplot(plot_df, aes(variable, value, group = spe, col = shapeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_sizeA <- ggplot(plot_df, aes(variable, value, group = spe, col = sizeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_shapeB <- ggplot(plot_df, aes(variable, value, group = spe, col = shapeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_sizeB <- ggplot(plot_df, aes(variable, value, group = spe, col = sizeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(arrangement = fig_arrangement, bg_type = fig_bg_type, 
                                                      shapeA = fig_shapeA, sizeA = fig_sizeA,
                                                      shapeB = fig_shapeB, sizeB = fig_sizeB)
    
  }
  
  
  # Combine the plots together by reference (and target) cell pairs
  plots_pair_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$arrangement,
                       all_plots_list[[pair]]$bg_type, 
                       all_plots_list[[pair]]$shapeA, 
                       all_plots_list[[pair]]$sizeA, 
                       all_plots_list[[pair]]$shapeB, 
                       all_plots_list[[pair]]$sizeB, 
                       nrow = 1, ncol = 6)
    
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
    plot_df$shapeA <- factor(plot_df$shapeA, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeA <- factor(plot_df$sizeA, c("Small", "Medium", "Large"))
    plot_df$shapeB <- factor(plot_df$shapeB, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeB <- factor(plot_df$sizeB, c("Small", "Medium", "Large"))
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
    
    fig_shapeA <- ggplot(plot_df, aes(variable, value, group = spe, col = shapeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_sizeA <- ggplot(plot_df, aes(variable, value, group = spe, col = sizeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_shapeB <- ggplot(plot_df, aes(variable, value, group = spe, col = shapeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    fig_sizeB <- ggplot(plot_df, aes(variable, value, group = spe, col = sizeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence") +
      geom_line(data = plot_bg_df, aes(variable, value), col = "black", linetype = 2) 
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(arrangement = fig_arrangement, bg_type = fig_bg_type, 
                                                               shapeA = fig_shapeA, sizeA = fig_sizeA,
                                                               shapeB = fig_shapeB, sizeB = fig_sizeB)
    
  }
  
  # Combine the plots together by cell types of interest
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$arrangement,
                       all_plots_list[[cell_types]]$bg_type, 
                       all_plots_list[[cell_types]]$shapeA, 
                       all_plots_list[[cell_types]]$sizeA,
                       all_plots_list[[cell_types]]$shapeB, 
                       all_plots_list[[cell_types]]$sizeB,
                       nrow = 1, ncol = 6)
    
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
    plot_df$key <- paste(plot_df$bg_type, plot_df$shapeA, plot_df$sizeA, plot_df$shapeB, plot_df$sizeB, sep = "_")
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shapeA <- factor(plot_df$shapeA, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeA <- factor(plot_df$sizeA, c("Small", "Medium", "Large"))
    plot_df$shapeB <- factor(plot_df$shapeB, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeB <- factor(plot_df$sizeB, c("Small", "Medium", "Large"))
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
    
    fig_shapeA <- ggplot(plot_df, aes(arrangement, AUC, group = key, col = shapeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("AUC") +
      geom_line(data = plot_bg_df, aes(arrangement, AUC), col = "black", linetype = 2)
    
    fig_sizeA <- ggplot(plot_df, aes(arrangement, AUC, group = key, col = sizeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("AUC") +
      geom_line(data = plot_bg_df, aes(arrangement, AUC), col = "black", linetype = 2)
    
    fig_shapeB <- ggplot(plot_df, aes(arrangement, AUC, group = key, col = shapeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("AUC") +
      geom_line(data = plot_bg_df, aes(arrangement, AUC), col = "black", linetype = 2)
    
    fig_sizeB <- ggplot(plot_df, aes(arrangement, AUC, group = key, col = sizeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("AUC") +
      geom_line(data = plot_bg_df, aes(arrangement, AUC), col = "black", linetype = 2)
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(bg_type = fig_bg_type, 
                                                      shapeA = fig_shapeA, sizeA = fig_sizeA,
                                                      shapeB = fig_shapeB, sizeB = fig_sizeB)
  }
  
  
  # Combine the plots together by reference target pairs
  library(cowplot)
  plots_pair_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$bg_type, 
                       all_plots_list[[pair]]$shapeA, all_plots_list[[pair]]$sizeA,
                       all_plots_list[[pair]]$shapeB, all_plots_list[[pair]]$sizeB,
                       nrow = 1, ncol = 5)
    
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
    plot_df$key <- paste(plot_df$bg_type, plot_df$shapeA, plot_df$sizeA, plot_df$shapeB, plot_df$sizeB, sep = "_")
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shapeA <- factor(plot_df$shapeA, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeA <- factor(plot_df$sizeA, c("Small", "Medium", "Large"))
    plot_df$shapeB <- factor(plot_df$shapeB, c("Sphere", "Ellipsoid", "Network"))
    plot_df$sizeB <- factor(plot_df$sizeB, c("Small", "Medium", "Large"))
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
    
    fig_shapeA <- ggplot(plot_df, aes(arrangement, AUC, group = key, col = shapeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("AUC") +
      geom_line(data = plot_bg_df, aes(arrangement, AUC), col = "black", linetype = 2)
    
    fig_sizeA <- ggplot(plot_df, aes(arrangement, AUC, group = key, col = sizeA)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("AUC") +
      geom_line(data = plot_bg_df, aes(arrangement, AUC), col = "black", linetype = 2)
    
    fig_shapeB <- ggplot(plot_df, aes(arrangement, AUC, group = key, col = shapeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("AUC") +
      geom_line(data = plot_bg_df, aes(arrangement, AUC), col = "black", linetype = 2)
    
    fig_sizeB <- ggplot(plot_df, aes(arrangement, AUC, group = key, col = sizeB)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("AUC") +
      geom_line(data = plot_bg_df, aes(arrangement, AUC), col = "black", linetype = 2)
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(bg_type = fig_bg_type, 
                                                               shapeA = fig_shapeA, sizeA = fig_sizeA,
                                                               shapeB = fig_shapeB, sizeB = fig_sizeB)
  }
  
  # Combine the plots together by cell types of interest
  library(cowplot)
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$bg_type, 
                       all_plots_list[[cell_types]]$shapeA, all_plots_list[[cell_types]]$sizeA,
                       all_plots_list[[cell_types]]$shapeB, all_plots_list[[cell_types]]$sizeB,
                       nrow = 1, ncol = 5)
    
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


### 2.1. separated spes APD ------------------------------------------------------
### 2.2. separated spes AMD ------------------------------------------------------

# Read separated_spes_table
setwd("~/Objects/supervised/spes_table")
separated_spes_table <- read.table("separated_spes_table.csv")

# Read separated_AMD_df
setwd("~/Objects/supervised/separated_spes/analysis_3D")
separated_AMD_df <- read.table("separated_AMD_df.csv")

# Read bg_AMD_df
setwd("~/Objects/supervised/background_spe")
bg_AMD_df <- read.table("bg_AMD_df.csv")

separated_AMD_plot <- plot_AMD_metric(separated_spes_table, separated_AMD_df, bg_AMD_df, separated_arrangements)

setwd("~/Objects/separated_spes/analysis_3D/plots")
# saveRDS(separated_AMD_plot, "separated_AMD_plot.rds")




### 2.3. separated spes MS, NMS, ACINP, AE -----------------------------------

# Read separated_spes_table
setwd("~/Objects/supervised/spes_table")
separated_spes_table <- read.table("separated_spes_table.csv")

# Read separated MS, NMS, ACINP, AE dfs
setwd("~/Objects/supervised/separated_spes/analysis_3D")
separated_MS_df <- read.table("separated_MS_df.csv")
separated_NMS_df <- read.table("separated_NMS_df.csv")
separated_ACINP_df <- read.table("separated_ACINP_df.csv")
separated_AE_df <- read.table("separated_AE_df.csv")

# Read bg_dfs
setwd("~/Objects/supervised/background_spe")
bg_MS_df <- read.table("bg_MS_df.csv")
bg_NMS_df <- read.table("bg_NMS_df.csv")
bg_ACINP_df <- read.table("bg_ACINP_df.csv")
bg_AE_df <- read.table("bg_AE_df.csv")

separated_MS_plot <- plot_gradient_metrics_type1(separated_spes_table, separated_MS_df, bg_MS_df, "MS", separated_arrangements)
separated_NMS_plot <- plot_gradient_metrics_type1(separated_spes_table, separated_NMS_df, bg_NMS_df, "NMS", separated_arrangements)
separated_ACINP_plot <- plot_gradient_metrics_type1(separated_spes_table, separated_ACINP_df, bg_ACINP_df, "ACINP", separated_arrangements)
separated_AE_plot <- plot_gradient_metrics_type1(separated_spes_table, separated_AE_df, bg_AE_df, "AE", separated_arrangements)

setwd("~/Objects/separated_spes/analysis_3D/plots")
# saveRDS()

### 2.4. separated spes ACIN, CKR ------------------------------------------------

# Read separated_spes_table
setwd("~/Objects/supervised/spes_table")
separated_spes_table <- read.table("separated_spes_table.csv")

# Read separated ACIN, CKR
setwd("~/Objects/supervised/separated_spes/analysis_3D")
separated_ACIN_df <- read.table("separated_ACIN_df.csv")
separated_CKR_df <- read.table("separated_CKR_df.csv")

# Read bg_dfs
setwd("~/Objects/supervised/background_spe")
bg_ACIN_df <- read.table("bg_ACIN_df.csv")
bg_CKR_df <- read.table("bg_CKR_df.csv")

# Get plots
separated_ACIN_plot <- plot_gradient_metrics_type2(separated_spes_table, separated_ACIN_df, bg_ACIN_df, "ACIN", separated_arrangements, 0, 50)

separated_CKR_plot <- plot_gradient_metrics_type2(separated_spes_table, separated_CKR_df, bg_CKR_df, "CKR", separated_arrangements, 15, 50)


### 2.5. separated spes SAC ------------------------------------------------------

# Read separated_spes_table
setwd("~/Objects/supervised/spes_table")
separated_spes_table <- read.table("separated_spes_table.csv")

# Read separated_SAC_df
setwd("~/Objects/supervised/separated_spes/analysis_3D")
separated_prop_SAC_df <- read.table("separated_prop_SAC_df.csv")
separated_entropy_SAC_df <- read.table("separated_entropy_SAC_df.csv")

# Read bg_SAC_df
setwd("~/Objects/supervised/background_spe")
bg_prop_SAC_df <- read.table("bg_prop_SAC_df.csv")
bg_entropy_SAC_df <- read.table("bg_entropy_SAC_df.csv")

separated_prop_SAC_plot <- plot_proportion_SAC(separated_spes_table, separated_prop_SAC_df, bg_prop_SAC_df, separated_arrangements)
separated_entropy_SAC_plot <- plot_entropy_SAC(separated_spes_table, separated_entropy_SAC_df, bg_entropy_SAC_df, separated_arrangements)

# setwd("~/Objects/separated_spes/analysis_3D/plots")
# saveRDS(separated_SAC_plot, "separated_SAC_plot.rds")

### 2.6. separated spes prevalence ------------------------------------------------

# Read separated_spes_table
setwd("~/Objects/supervised/spes_table")
separated_spes_table <- read.table("separated_spes_table.csv")

# Read separated prevalence dfs
setwd("~/Objects/supervised/separated_spes/analysis_3D")
separated_prop_prevalence_df <- read.table("separated_prop_prevalence_df.csv")
separated_entropy_prevalence_df <- read.table("separated_entropy_prevalence_df.csv")

# Read bg prevalence dfs
setwd("~/Objects/supervised/background_spe")
bg_prop_prevalence_df <- read.table("bg_prop_prevalence_df.csv")
bg_entropy_prevalence_df <- read.table("bg_entropy_prevalence_df.csv")

separated_prop_prevalence_plot <- plot_proportion_prevalence(separated_spes_table, separated_prop_prevalence_df, bg_prop_prevalence_df, separated_arrangements)
separated_entropy_prevalence_plot <- plot_entropy_prevalence(separated_spes_table, separated_entropy_prevalence_df, bg_entropy_prevalence_df, separated_arrangements)

separated_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(separated_spes_table, separated_prop_prevalence_df, bg_prop_prevalence_df, separated_arrangements)
separated_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(separated_spes_table, separated_entropy_prevalence_df, bg_entropy_prevalence_df, separated_arrangements)

setwd("~/Objects/separated_spes/analysis_3D/plots")
# saveRDS(separated_prevalence_plot, "separated_prevalence_plot.rds")



