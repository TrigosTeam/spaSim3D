library(cowplot)
library(ggplot2)

### Read tables --------------------------------------------------------------
# Read mixed_spes_table
setwd("~/Objects/unsupervised/spes_table")
mixed_spes_table <- read.table("mixed_spes_table_unsupervised.csv")

# Read ringed_spes_table
setwd("~/Objects/unsupervised/spes_table")
ringed_spes_table <- read.table("ringed_spes_table_unsupervised.csv")

# Read separated_spes_table
setwd("~/Objects/unsupervised/spes_table")
separated_spes_table <- read.table("separated_spes_table_unsupervised.csv")
separated_spes_table$distance <- separated_spes_table$centre_x_coord_B - separated_spes_table$centre_x_coord_A

# Start with cluster_A
separated_A_spes_table <- separated_spes_table
colnames(separated_A_spes_table)[3:7] <- c("shape", "radius_x_E", "radius_y_E", "radius_z_E", "width_N")

# Then do cluster_B
separated_B_spes_table <- separated_spes_table
colnames(separated_B_spes_table)[9:13] <- c("shape", "radius_x_E", "radius_y_E", "radius_z_E", "width_N")

### 1.1.2. Function to get plot for AMD -------------------------------------
plot_AMD_metric <- function(spes_table, AMD_df, arrangement_colname) {
  
  # AMD pairs are A/A, A/B, B/A, B/B
  AMD_pairs <- data.frame(cell1 = c("A", "A", "B", "B"),
                          cell2 = c("A", "B", "A", "B"))
  AMD_pairs$pair <- paste(AMD_pairs$cell1, AMD_pairs$cell2, sep = "/")
  
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(nrow(AMD_pairs))) {
    
    # Subset AMD_df for chosen pair
    plot_df <- AMD_df[AMD_df$reference == AMD_pairs[i, "cell1"] & AMD_df$target == AMD_pairs[i, "cell2"], ]
    
    # Combine spes_table and AMD_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Remove spheres
    plot_df <- plot_df[plot_df$shape != "Sphere", ]
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    fig_bg_prop_A <- ggplot(plot_df, aes(!!sym(arrangement_colname), AMD, col = bg_prop_A)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_bg_prop_B <- ggplot(plot_df, aes(!!sym(arrangement_colname), AMD, col = bg_prop_B)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), AMD, col = shape)) +
      geom_point() +
      theme_bw()

    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AMD, col = variation_E)) +
      geom_point() +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AMD, col = volume_E)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), AMD, col = width_N)) +
      geom_point() +
      theme_bw()
    
    all_plots_list[[AMD_pairs[i, "pair"]]] <- list(bg_prop_A = fig_bg_prop_A + theme(legend.position="none"), 
                                                   bg_prop_B = fig_bg_prop_B + theme(legend.position="none"),
                                                   shape = fig_shape + theme(legend.position="none"),
                                                   variation_E = fig_variation_E + theme(legend.position="none"),
                                                   volume_E = fig_volume_E + theme(legend.position="none"),
                                                   width_N = fig_width_N + theme(legend.position="none"))
  }

  # Get legends
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_bg_prop_a, 
                       legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by pairs
  plots_pair_list <- list()
  
  for (i in seq(nrow(AMD_pairs))) {
    pair <- AMD_pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$bg_prop_A,
                       all_plots_list[[pair]]$bg_prop_B,
                       all_plots_list[[pair]]$shape, 
                       all_plots_list[[pair]]$variation_E,
                       all_plots_list[[pair]]$volume_E,
                       all_plots_list[[pair]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[pair]]))
    
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
                        legends,
                        nrow = 5, ncol = 1, 
                        rel_heights = c(1, 1, 1, 1, 0.5))
  
  methods::show(AMD_plot)
  
  return(AMD_plot)
}


### 1.2.1. Function to get plot for MS, NMS, ACINP, AE gradient metrics ------------
plot_gradient_metrics_type1 <- function(spes_table, gradient_metric_df, metric, arrangement_colname) {
  
  # Constants
  cell_types <- c("A", "B") # Use A as reference, and B as target, and vice versa
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(length(cell_types))) {
    # Subset gradient_metric_df for current reference cell
    plot_df <- gradient_metric_df[gradient_metric_df$reference == cell_types[i], ]
    
    # Combine spes_table and gradient_metric_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , radii_colnames)
    
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df$variable <- unfactor(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    
    # Remove spheres
    plot_df <- plot_df[plot_df$shape != "Sphere", ]
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))

    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    fig_bg_prop_A <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_prop_A)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_bg_prop_B <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_prop_B)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = variation_E)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = volume_E)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = spe, col = width_N)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    all_plots_list[[cell_types[i]]] <- list(arrangement = fig_arrangement + theme(legend.position = "none"), 
                                            bg_prop_A = fig_bg_prop_A + theme(legend.position = "none"), 
                                            bg_prop_B = fig_bg_prop_B + theme(legend.position = "none"),
                                            shape = fig_shape + theme(legend.position = "none"),
                                            variation_E = fig_variation_E + theme(legend.position = "none"),
                                            volume_E = fig_volume_E + theme(legend.position = "none"),
                                            width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  # Get legends
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_arrangement,
                       legend_bg_prop_a, 
                       legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference cell type
  plots_ref_list <- list()
  
  for (i in seq(length(cell_types))) {
    reference_cell_type <- cell_types[i]
    
    plots <- plot_grid(all_plots_list[[reference_cell_type]]$arrangement,
                       all_plots_list[[reference_cell_type]]$bg_prop_A,
                       all_plots_list[[reference_cell_type]]$bg_prop_B,
                       all_plots_list[[reference_cell_type]]$shape, 
                       all_plots_list[[reference_cell_type]]$variation_E,
                       all_plots_list[[reference_cell_type]]$volume_E,
                       all_plots_list[[reference_cell_type]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[reference_cell_type]]))
    
    title <- ggdraw() + 
      draw_label(paste("Reference:", reference_cell_type), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_ref_list[[reference_cell_type]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plot <- plot_grid(plots_ref_list$A,
                             plots_ref_list$B, 
                             legends,
                             nrow = 3, ncol = 1,
                             rel_heights = c(1, 1, 0.5))
  
  methods::show(combined_plot)
  
  return(combined_plot)
}


### 1.2.2. Function to get plot for ACIN, CKR gradient metrics ------------------
plot_gradient_metrics_type2 <- function(spes_table, gradient_metric_df, metric, arrangement_colname, min_radius, max_radius) {
  
  # Constants
  pairs <- data.frame(cell1 = c("A", "A", "B", "B"),
                      cell2 = c("A", "B", "A", "B"))
  pairs$pair <- paste(pairs$cell1, pairs$cell2, sep = "/")
  
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  
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
    
    # Remove spheres
    plot_df <- plot_df[plot_df$shape != "Sphere", ]
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    fig_bg_prop_A <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_prop_A)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_bg_prop_B <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_prop_B)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      labs(x = "radius", y = metric) +
      geom_line() +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = variation_E)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = volume_E)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = spe, col = width_N)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    all_plots_list[[pairs[i, "pair"]]] <- list(arrangement = fig_arrangement + theme(legend.position = "none"), 
                                            bg_prop_A = fig_bg_prop_A + theme(legend.position = "none"), 
                                            bg_prop_B = fig_bg_prop_B + theme(legend.position = "none"),
                                            shape = fig_shape + theme(legend.position = "none"),
                                            variation_E = fig_variation_E + theme(legend.position = "none"),
                                            volume_E = fig_volume_E + theme(legend.position = "none"),
                                            width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  
  # Get legends
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_arrangement,
                       legend_bg_prop_a, 
                       legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference cell type
  plots_pair_list <- list()
  
  for (i in seq(nrow(pairs))) {
    pair <- pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$arrangement,
                       all_plots_list[[pair]]$bg_prop_A,
                       all_plots_list[[pair]]$bg_prop_B,
                       all_plots_list[[pair]]$shape, 
                       all_plots_list[[pair]]$variation_E,
                       all_plots_list[[pair]]$volume_E,
                       all_plots_list[[pair]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[pair]]))
    
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
                             legends,
                             nrow = 5, ncol = 1,
                             rel_heights = c(1, 1, 1, 1, 0.5))
  
  methods::show(combined_plot)
  
  return(combined_plot)
}





### 1.3.1. Function to get plot for proportion SAC ----------------------------------------
plot_proportion_SAC <- function(spes_table, SAC_df, arrangement_colname) {
  
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
    
    # Remove spheres
    plot_df <- plot_df[plot_df$shape != "Sphere", ]
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    
    fig_bg_prop_A <- ggplot(plot_df, aes(!!sym(arrangement_colname), proportion, col = bg_prop_A)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_bg_prop_B <- ggplot(plot_df, aes(!!sym(arrangement_colname), proportion, col = bg_prop_B)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), proportion, col = shape)) +
      geom_point() +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), proportion, col = variation_E)) +
      geom_point() +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), proportion, col = volume_E)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), proportion, col = width_N)) +
      geom_point() +
      theme_bw()
    
    all_plots_list[[prop_cell_types[i, "pair"]]] <- list(bg_prop_A = fig_bg_prop_A + theme(legend.position="none"), 
                                                         bg_prop_B = fig_bg_prop_B + theme(legend.position="none"),
                                                         shape = fig_shape + theme(legend.position="none"),
                                                         variation_E = fig_variation_E + theme(legend.position="none"),
                                                         volume_E = fig_volume_E + theme(legend.position="none"),
                                                         width_N = fig_width_N + theme(legend.position="none"))
  }
  
  
  # Get legends
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_bg_prop_a, 
                       legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by pairs
  
  plots_pair_list <- list()
  
  for (i in seq(nrow(prop_cell_types))) {
    pair <- prop_cell_types[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$bg_prop_A,
                       all_plots_list[[pair]]$bg_prop_B,
                       all_plots_list[[pair]]$shape, 
                       all_plots_list[[pair]]$variation_E,
                       all_plots_list[[pair]]$volume_E,
                       all_plots_list[[pair]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[pair]]))
    
    title <- ggdraw() + 
      draw_label(paste("Reference/Target:", pair), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_pair_list[[pair]] <- fig
  }
  
  # Combine the combined plots into one big plot
  SAC_plot <- plot_grid(plots_pair_list[[prop_cell_types$pair[1]]], 
                        plots_pair_list[[prop_cell_types$pair[2]]],
                        legends,
                        nrow = 3, ncol = 1, 
                        rel_heights = c(1, 1, 0.5))
  
  methods::show(SAC_plot)
  
  return(SAC_plot)
}



### 1.3.2. Function to get plot for entropy SAC ----------------------------------------
plot_entropy_SAC <- function(spes_table, SAC_df, arrangement_colname) {
  
  # Get possible cell type of interest combinations
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    
    # Subset for current cell type of interest combination
    plot_df <- SAC_df[SAC_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Combine spes_table and SAC_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Remove spheres
    plot_df <- plot_df[plot_df$shape != "Sphere", ]
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    fig_bg_prop_A <- ggplot(plot_df, aes(!!sym(arrangement_colname), entropy, col = bg_prop_A)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_bg_prop_B <- ggplot(plot_df, aes(!!sym(arrangement_colname), entropy, col = bg_prop_B)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), entropy, col = shape)) +
      geom_point() +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), entropy, col = variation_E)) +
      geom_point() +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), entropy, col = volume_E)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), entropy, col = width_N)) +
      geom_point() +
      theme_bw()
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(bg_prop_A = fig_bg_prop_A + theme(legend.position="none"), 
                                                               bg_prop_B = fig_bg_prop_B + theme(legend.position="none"),
                                                               shape = fig_shape + theme(legend.position="none"),
                                                               variation_E = fig_variation_E + theme(legend.position="none"),
                                                               volume_E = fig_volume_E + theme(legend.position="none"),
                                                               width_N = fig_width_N + theme(legend.position="none"))
  }
  
  # Get legends
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_bg_prop_a, 
                       legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by cell types of interest
  
  plots_cell_types_list <- list()
  
  for (i in seq(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$bg_prop_A,
                       all_plots_list[[cell_types]]$bg_prop_B,
                       all_plots_list[[cell_types]]$shape, 
                       all_plots_list[[cell_types]]$variation_E,
                       all_plots_list[[cell_types]]$volume_E,
                       all_plots_list[[cell_types]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[cell_types]]))
    
    title <- ggdraw() + 
      draw_label(paste("Cell types of interest:", cell_types), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_cell_types_list[[cell_types]] <- fig
  }
  
  
  # Combine the combined plots into one big plot
  SAC_plot <- plot_grid(plots_cell_types_list[[entropy_cell_types$cell_types[1]]], 
                        plots_cell_types_list[[entropy_cell_types$cell_types[2]]],
                        legends,
                        nrow = 3, ncol = 1,
                        rel_heights = c(1, 1, 0.5))
  
  methods::show(SAC_plot)

  return(SAC_plot)
}
### 1.3.3. Function to get plot for proportion prevalence ----------------------------------
plot_proportion_prevalence <- function(spes_table, prevalence_df, arrangement_colname) {
  
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
    
    # Remove spheres
    plot_df <- plot_df[plot_df$shape != "Sphere", ]
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    fig_bg_prop_A <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_prop_A)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_bg_prop_B <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_prop_B)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      labs(x = "threshold", y = "prevalence") +
      geom_line() +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = variation_E)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = volume_E)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = spe, col = width_N)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(arrangement = fig_arrangement + theme(legend.position = "none"), 
                                                      bg_prop_A = fig_bg_prop_A + theme(legend.position = "none"), 
                                                      bg_prop_B = fig_bg_prop_B + theme(legend.position = "none"),
                                                      shape = fig_shape + theme(legend.position = "none"),
                                                      variation_E = fig_variation_E + theme(legend.position = "none"),
                                                      volume_E = fig_volume_E + theme(legend.position = "none"),
                                                      width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  # Get legends
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_arrangement,
                       legend_bg_prop_a, 
                       legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  
  # Combine the plots together by reference (and target) cell pairs
  plots_pair_list <- list()

  for (i in seq(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$arrangement,
                       all_plots_list[[pair]]$bg_prop_A,
                       all_plots_list[[pair]]$bg_prop_B,
                       all_plots_list[[pair]]$shape, 
                       all_plots_list[[pair]]$variation_E,
                       all_plots_list[[pair]]$volume_E,
                       all_plots_list[[pair]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[pair]]))
    
    title <- ggdraw() + 
      draw_label(paste("Reference/Target:", pair), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_pair_list[[pair]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plot <- plot_grid(plots_pair_list[[prop_cell_types$pair[1]]],
                             plots_pair_list[[prop_cell_types$pair[2]]], 
                             legends,
                             nrow = 3, ncol = 1,
                             rel_heights = c(1, 1, 0.5))
  
  methods::show(combined_plot)
  
  return(combined_plot)
}



### 1.3.4. Function to get plot for entropy prevalence ----------------------------------
plot_entropy_prevalence <- function(spes_table, prevalence_df, arrangement_colname) {
  
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
    
    # Remove spheres
    plot_df <- plot_df[plot_df$shape != "Sphere", ]
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    fig_bg_prop_A <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_prop_A)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_bg_prop_B <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_prop_B)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      labs(x = "threshold", y = "prevalence") +
      geom_line() +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = variation_E)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = volume_E)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = spe, col = width_N)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(arrangement = fig_arrangement + theme(legend.position = "none"), 
                                                               bg_prop_A = fig_bg_prop_A + theme(legend.position = "none"), 
                                                               bg_prop_B = fig_bg_prop_B + theme(legend.position = "none"),
                                                               shape = fig_shape + theme(legend.position = "none"),
                                                               variation_E = fig_variation_E + theme(legend.position = "none"),
                                                               volume_E = fig_volume_E + theme(legend.position = "none"),
                                                               width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  # Get legends
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_arrangement,
                       legend_bg_prop_a, 
                       legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  
  # Combine the plots together by reference (and target) cell pairs
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$arrangement,
                       all_plots_list[[cell_types]]$bg_prop_A,
                       all_plots_list[[cell_types]]$bg_prop_B,
                       all_plots_list[[cell_types]]$shape, 
                       all_plots_list[[cell_types]]$variation_E,
                       all_plots_list[[cell_types]]$volume_E,
                       all_plots_list[[cell_types]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[cell_types]]))
    
    title <- ggdraw() + 
      draw_label(paste("Cell types of interest:", cell_types), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_cell_types_list[[cell_types]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plot <- plot_grid(plots_cell_types_list[[entropy_cell_types$cell_types[1]]],
                             plots_cell_types_list[[entropy_cell_types$cell_types[2]]], 
                             legends,
                             nrow = 3, ncol = 1,
                             rel_heights = c(1, 1, 0.5))
  
  methods::show(combined_plot)
  
  return(combined_plot)
}
### 1.3.5. Function to get plot for proportion prevalence AUC -----------------
plot_proportion_prevalence_AUC <- function(spes_table, prevalence_df, arrangement_colname) {
  
  # Constants
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
  prop_cell_types$pair <- paste(prop_cell_types$ref, prop_cell_types$tar, sep = "/")
  
  # Get AUC for each prevalence gradient
  prevalence_df$AUC <- apply(prevalence_df[ , threshold_colnames], 1, sum) * 0.01
  prevalence_df <- prevalence_df[ , c("spe", "reference", "target", "AUC")]
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    
    # Subset for current reference (and target) cell
    plot_df <- prevalence_df[prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Combine spes_table and updated prevalence_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Remove spheres
    plot_df <- plot_df[plot_df$shape != "Sphere", ]
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    fig_bg_prop_A <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = bg_prop_A)) +
      geom_point() +
      ylab("AUC") +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_bg_prop_B <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = bg_prop_B)) +
      geom_point() +
      ylab("AUC") +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = shape)) +
      geom_point() +
      ylab("AUC") +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AUC, col = variation_E)) +
      geom_point() +
      ylab("AUC") +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AUC, col = volume_E)) +
      geom_point() +
      ylab("AUC") +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), AUC, col = width_N)) +
      geom_point() +
      ylab("AUC") +
      theme_bw()
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(bg_prop_A = fig_bg_prop_A + theme(legend.position="none"), 
                                                      bg_prop_B = fig_bg_prop_B + theme(legend.position="none"),
                                                      shape = fig_shape + theme(legend.position="none"),
                                                      variation_E = fig_variation_E + theme(legend.position="none"),
                                                      volume_E = fig_volume_E + theme(legend.position="none"),
                                                      width_N = fig_width_N + theme(legend.position="none"))
  }
  
  # Get legends
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_bg_prop_a, 
                       legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference target pairs
  plots_pair_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$bg_prop_A,
                       all_plots_list[[pair]]$bg_prop_B,
                       all_plots_list[[pair]]$shape, 
                       all_plots_list[[pair]]$variation_E,
                       all_plots_list[[pair]]$volume_E,
                       all_plots_list[[pair]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[pair]]))
    
    title <- ggdraw() + 
      draw_label(paste("Reference/Target:", pair), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_pair_list[[pair]] <- fig
  }
  
  # Combine the combined plots into one big plot
  AUC_plot <- plot_grid(plots_pair_list[[prop_cell_types$pair[1]]], 
                        plots_pair_list[[prop_cell_types$pair[2]]],    
                        legends,
                        nrow = 3, ncol = 1, 
                        rel_heights = c(1, 1, 0.5))
  
  methods::show(AUC_plot)
  
  return(AUC_plot)
}

### 1.3.6. Function to get plot for entropy prevalence AUC -----------------
plot_entropy_prevalence_AUC <- function(spes_table, prevalence_df, arrangement_colname) {
  
  # Constants
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  
  # Get AUC for each prevalence gradient
  prevalence_df$AUC <- apply(prevalence_df[ , threshold_colnames], 1, sum) * 0.01
  prevalence_df <- prevalence_df[ , c("spe", "cell_types", "AUC")]
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    
    # Subset for current cell type of interest combintation
    plot_df <- prevalence_df[prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Combine spes_table and updated prevalence_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Remove spheres
    plot_df <- plot_df[plot_df$shape != "Sphere", ]
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    fig_bg_prop_A <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = bg_prop_A)) +
      geom_point() +
      ylab("AUC") +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_bg_prop_B <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = bg_prop_B)) +
      geom_point() +
      ylab("AUC") +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = shape)) +
      geom_point() +
      ylab("AUC") +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AUC, col = variation_E)) +
      geom_point() +
      ylab("AUC") +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AUC, col = volume_E)) +
      geom_point() +
      ylab("AUC") +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), AUC, col = width_N)) +
      geom_point() +
      ylab("AUC") +
      theme_bw()
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(bg_prop_A = fig_bg_prop_A + theme(legend.position="none"), 
                                                               bg_prop_B = fig_bg_prop_B + theme(legend.position="none"),
                                                               shape = fig_shape + theme(legend.position="none"),
                                                               variation_E = fig_variation_E + theme(legend.position="none"),
                                                               volume_E = fig_volume_E + theme(legend.position="none"),
                                                               width_N = fig_width_N + theme(legend.position="none"))
  }
  
  # Get legends
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_bg_prop_a, 
                       legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference target pairs
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$bg_prop_A,
                       all_plots_list[[cell_types]]$bg_prop_B,
                       all_plots_list[[cell_types]]$shape, 
                       all_plots_list[[cell_types]]$variation_E,
                       all_plots_list[[cell_types]]$volume_E,
                       all_plots_list[[cell_types]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[cell_types]]))
    
    title <- ggdraw() + 
      draw_label(paste("Cell types of interest:", cell_types), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_cell_types_list[[cell_types]] <- fig
  }
  
  # Combine the combined plots into one big plot
  AUC_plot <- plot_grid(plots_cell_types_list[[entropy_cell_types$cell_types[1]]], 
                        plots_cell_types_list[[entropy_cell_types$cell_types[2]]],    
                        legends,
                        nrow = 3, ncol = 1, 
                        rel_heights = c(1, 1, 0.5))
  
  methods::show(AUC_plot)
  
  return(AUC_plot)
}
### 2.2. mixed spes AMD ------------------------------------------------------

# Read mixed_AMD_df
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_AMD_df <- read.table("mixed_AMD_df.csv")

mixed_AMD_plot <- plot_AMD_metric(mixed_spes_table, mixed_AMD_df, "cluster_prop_B")

setwd("~/Objects/unsupervised/plots/original")
saveRDS(mixed_AMD_plot, "mixed_AMD_plot.RDS")

### 2.3. mixed spes MS, NMS, ACINP, AE -----------------------------------

# Read mixed MS, NMS, ACINP, AE dfs
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_MS_df <- read.table("mixed_MS_df.csv")
mixed_NMS_df <- read.table("mixed_NMS_df.csv")
mixed_ACINP_df <- read.table("mixed_ACINP_df.csv")
mixed_AE_df <- read.table("mixed_AE_df.csv")

mixed_MS_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_MS_df, "MS", "cluster_prop_B")
mixed_NMS_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_NMS_df, "NMS", "cluster_prop_B")
mixed_ACINP_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_ACINP_df, "ACINP", "cluster_prop_B")
mixed_AE_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_AE_df, "AE", "cluster_prop_B")

setwd("~/Objects/unsupervised/plots/original")
saveRDS(mixed_MS_plot, "mixed_MS_plot.RDS")
saveRDS(mixed_NMS_plot, "mixed_NMS_plot.RDS")
saveRDS(mixed_ACINP_plot, "mixed_ACINP_plot.RDS")
saveRDS(mixed_AE_plot, "mixed_AE_plot.RDS")

### 2.4. mixed spes ACIN, CKR ------------------------------------------------

# Read mixed ACIN, CKR
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_ACIN_df <- read.table("mixed_ACIN_df.csv")
mixed_CKR_df <- read.table("mixed_CKR_df.csv")

# Get plots
mixed_ACIN_plot <- plot_gradient_metrics_type2(mixed_spes_table, mixed_ACIN_df, "ACIN", "cluster_prop_B", 20, 100)
mixed_CKR_plot <- plot_gradient_metrics_type2(mixed_spes_table, mixed_CKR_df, "CKR", "cluster_prop_B", 20, 100)

setwd("~/Objects/unsupervised/plots/original")
saveRDS(mixed_ACIN_plot, "mixed_ACIN_plot.RDS")
saveRDS(mixed_CKR_plot, "mixed_CKR_plot.RDS")

### 2.5. mixed spes SAC ------------------------------------------------------

# Read mixed_SAC_df
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_prop_SAC_df <- read.table("mixed_prop_SAC_df.csv")
mixed_entropy_SAC_df <- read.table("mixed_entropy_SAC_df.csv")

mixed_prop_SAC_plot <- plot_proportion_SAC(mixed_spes_table, mixed_prop_SAC_df, "cluster_prop_B")
mixed_entropy_SAC_plot <- plot_entropy_SAC(mixed_spes_table, mixed_entropy_SAC_df, "cluster_prop_B")

setwd("~/Objects/unsupervised/plots/original")
saveRDS(mixed_prop_SAC_plot, "mixed_prop_SAC_plot.RDS")
saveRDS(mixed_entropy_SAC_plot, "mixed_entropy_SAC_plot.RDS")

### 2.6. mixed spes prevalence ------------------------------------------------

# Read mixed prevalence dfs
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_prop_prevalence_df <- read.table("mixed_prop_prevalence_df.csv")
mixed_entropy_prevalence_df <- read.table("mixed_entropy_prevalence_df.csv")

mixed_prop_prevalence_plot <- plot_proportion_prevalence(mixed_spes_table, mixed_prop_prevalence_df, "cluster_prop_B")
mixed_entropy_prevalence_plot <- plot_entropy_prevalence(mixed_spes_table, mixed_entropy_prevalence_df, "cluster_prop_B")
mixed_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(mixed_spes_table, mixed_prop_prevalence_df, "cluster_prop_B")
mixed_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(mixed_spes_table, mixed_entropy_prevalence_df, "cluster_prop_B")

setwd("~/Objects/unsupervised/plots/original")
saveRDS(mixed_prop_prevalence_plot, "mixed_prop_prevalence_plot.RDS")
saveRDS(mixed_entropy_prevalence_plot, "mixed_entropy_prevalence_plot.RDS")
saveRDS(mixed_prop_prevalence_AUC_plot, "mixed_prop_prevalence_AUC_plot.RDS")
saveRDS(mixed_entropy_prevalence_AUC_plot, "mixed_entropy_prevalence_AUC_plot.RDS")



### 3.2. ringed spes AMD ------------------------------------------------------

# Read ringed_AMD_df
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_AMD_df <- read.table("ringed_AMD_df.csv")

ringed_AMD_plot <- plot_AMD_metric(ringed_spes_table, ringed_AMD_df, "width_ring_factor")

setwd("~/Objects/unsupervised/plots/original")
saveRDS(ringed_AMD_plot, "ringed_AMD_plot.RDS")

### 3.3. ringed spes MS, NMS, ACINP, AE -----------------------------------

# Read ringed MS, NMS, ACINP, AE dfs
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_MS_df <- read.table("ringed_MS_df.csv")
ringed_NMS_df <- read.table("ringed_NMS_df.csv")
ringed_ACINP_df <- read.table("ringed_ACINP_df.csv")
ringed_AE_df <- read.table("ringed_AE_df.csv")

ringed_MS_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_MS_df, "MS", "width_ring_factor")
ringed_NMS_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_NMS_df, "NMS", "width_ring_factor")
ringed_ACINP_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_ACINP_df, "ACINP", "width_ring_factor")
ringed_AE_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_AE_df, "AE", "width_ring_factor")

setwd("~/Objects/unsupervised/plots/original")
saveRDS(ringed_MS_plot, "ringed_MS_plot.RDS")
saveRDS(ringed_NMS_plot, "ringed_NMS_plot.RDS")
saveRDS(ringed_ACINP_plot, "ringed_ACINP_plot.RDS")
saveRDS(ringed_AE_plot, "ringed_AE_plot.RDS")

### 3.4. ringed spes ACIN, CKR ------------------------------------------------

# Read ringed ACIN, CKR
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_ACIN_df <- read.table("ringed_ACIN_df.csv")
ringed_CKR_df <- read.table("ringed_CKR_df.csv")

# Get plots
ringed_ACIN_plot <- plot_gradient_metrics_type2(ringed_spes_table, ringed_ACIN_df, "ACIN", "width_ring_factor", 20, 100)
ringed_CKR_plot <- plot_gradient_metrics_type2(ringed_spes_table, ringed_CKR_df, "CKR", "width_ring_factor", 20, 100)

setwd("~/Objects/unsupervised/plots/original")
saveRDS(ringed_ACIN_plot, "ringed_ACIN_plot.RDS")
saveRDS(ringed_CKR_plot, "ringed_CKR_plot.RDS")

### 3.5. ringed spes SAC ------------------------------------------------------

# Read ringed_SAC_df
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_prop_SAC_df <- read.table("ringed_prop_SAC_df.csv")
ringed_entropy_SAC_df <- read.table("ringed_entropy_SAC_df.csv")

ringed_prop_SAC_plot <- plot_proportion_SAC(ringed_spes_table, ringed_prop_SAC_df, "width_ring_factor")
ringed_entropy_SAC_plot <- plot_entropy_SAC(ringed_spes_table, ringed_entropy_SAC_df, "width_ring_factor")

setwd("~/Objects/unsupervised/plots/original")
saveRDS(ringed_prop_SAC_plot, "ringed_prop_SAC_plot.RDS")
saveRDS(ringed_entropy_SAC_plot, "ringed_entropy_SAC_plot.RDS")

### 3.6. ringed spes prevalence ------------------------------------------------

# Read ringed prevalence dfs
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_prop_prevalence_df <- read.table("ringed_prop_prevalence_df.csv")
ringed_entropy_prevalence_df <- read.table("ringed_entropy_prevalence_df.csv")

ringed_prop_prevalence_plot <- plot_proportion_prevalence(ringed_spes_table, ringed_prop_prevalence_df, "width_ring_factor")
ringed_entropy_prevalence_plot <- plot_entropy_prevalence(ringed_spes_table, ringed_entropy_prevalence_df, "width_ring_factor")
ringed_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(ringed_spes_table, ringed_prop_prevalence_df, "width_ring_factor")
ringed_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(ringed_spes_table, ringed_entropy_prevalence_df, "width_ring_factor")

setwd("~/Objects/unsupervised/plots/original")
saveRDS(ringed_prop_prevalence_plot, "ringed_prop_prevalence_plot.RDS")
saveRDS(ringed_entropy_prevalence_plot, "ringed_entropy_prevalence_plot.RDS")
saveRDS(ringed_prop_prevalence_AUC_plot, "ringed_prop_prevalence_AUC_plot.RDS")
saveRDS(ringed_entropy_prevalence_AUC_plot, "ringed_entropy_prevalence_AUC_plot.RDS")



### 4.2. separated spes AMD ------------------------------------------------------

# Read separated_AMD_df
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_AMD_df <- read.table("separated_AMD_df.csv")


separated_A_AMD_plot <- plot_AMD_metric(separated_A_spes_table, separated_AMD_df, "distance")
separated_B_AMD_plot <- plot_AMD_metric(separated_B_spes_table, separated_AMD_df, "distance")

setwd("~/Objects/unsupervised/plots/original")
saveRDS(separated_A_AMD_plot, "separated_A_AMD_plot.RDS")
saveRDS(separated_B_AMD_plot, "separated_B_AMD_plot.RDS")

### 4.3. separated spes MS, NMS, ACINP, AE -----------------------------------

# Read separated MS, NMS, ACINP, AE dfs
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_MS_df <- read.table("separated_MS_df.csv")
separated_NMS_df <- read.table("separated_NMS_df.csv")
separated_ACINP_df <- read.table("separated_ACINP_df.csv")
separated_AE_df <- read.table("separated_AE_df.csv")

# Start with cluster_A
separated_A_MS_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_MS_df, "MS", "distance")
separated_A_NMS_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_NMS_df, "NMS", "distance")
separated_A_ACINP_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_ACINP_df, "ACINP", "distance")
separated_A_AE_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_AE_df, "AE", "distance")

# Then do cluster_B
separated_B_MS_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_MS_df, "MS", "distance")
separated_B_NMS_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_NMS_df, "NMS", "distance")
separated_B_ACINP_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_ACINP_df, "ACINP", "distance")
separated_B_AE_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_AE_df, "AE", "distance")

setwd("~/Objects/unsupervised/plots/original")
saveRDS(separated_A_MS_plot, "separated_A_MS_plot.RDS")
saveRDS(separated_A_NMS_plot, "separated_A_NMS_plot.RDS")
saveRDS(separated_A_ACINP_plot, "separated_A_ACINP_plot.RDS")
saveRDS(separated_A_AE_plot, "separated_A_AE_plot.RDS")

saveRDS(separated_B_MS_plot, "separated_B_MS_plot.RDS")
saveRDS(separated_B_NMS_plot, "separated_B_NMS_plot.RDS")
saveRDS(separated_B_ACINP_plot, "separated_B_ACINP_plot.RDS")
saveRDS(separated_B_AE_plot, "separated_B_AE_plot.RDS")

### 4.4. separated spes ACIN, CKR ------------------------------------------------

# Read separated ACIN, CKR
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_ACIN_df <- read.table("separated_ACIN_df.csv")
separated_CKR_df <- read.table("separated_CKR_df.csv")

# Start with cluster_A
separated_A_ACIN_plot <- plot_gradient_metrics_type2(separated_A_spes_table, separated_ACIN_df, "ACIN", "distance", 20, 100)
separated_A_CKR_plot <- plot_gradient_metrics_type2(separated_A_spes_table, separated_CKR_df, "CKR", "distance", 20, 100)

# Then do cluster_B
separated_B_ACIN_plot <- plot_gradient_metrics_type2(separated_B_spes_table, separated_ACIN_df, "ACIN", "distance", 20, 100)
separated_B_CKR_plot <- plot_gradient_metrics_type2(separated_B_spes_table, separated_CKR_df, "CKR", "distance", 20, 100)

setwd("~/Objects/unsupervised/plots/original")
saveRDS(separated_A_ACIN_plot, "separated_A_ACIN_plot.RDS")
saveRDS(separated_A_CKR_plot, "separated_A_CKR_plot.RDS")

saveRDS(separated_B_ACIN_plot, "separated_B_ACIN_plot.RDS")
saveRDS(separated_B_CKR_plot, "separated_B_CKR_plot.RDS")

### 4.5. separated spes SAC ------------------------------------------------------

# Read separated_SAC_df
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_prop_SAC_df <- read.table("separated_prop_SAC_df.csv")
separated_entropy_SAC_df <- read.table("separated_entropy_SAC_df.csv")

# Start with cluster_A
separated_A_prop_SAC_plot <- plot_proportion_SAC(separated_A_spes_table, separated_prop_SAC_df, "distance")
separated_A_entropy_SAC_plot <- plot_entropy_SAC(separated_A_spes_table, separated_entropy_SAC_df, "distance")

# Then do cluster_B
separated_B_prop_SAC_plot <- plot_proportion_SAC(separated_B_spes_table, separated_prop_SAC_df, "distance")
separated_B_entropy_SAC_plot <- plot_entropy_SAC(separated_B_spes_table, separated_entropy_SAC_df, "distance")

setwd("~/Objects/unsupervised/plots/original")
saveRDS(separated_A_prop_SAC_plot, "separated_A_prop_SAC_plot.RDS")
saveRDS(separated_A_entropy_SAC_plot, "separated_A_entropy_SAC_plot.RDS")

saveRDS(separated_B_prop_SAC_plot, "separated_B_prop_SAC_plot.RDS")
saveRDS(separated_B_entropy_SAC_plot, "separated_B_entropy_SAC_plot.RDS")

### 4.6. separated spes prevalence ------------------------------------------------

# Read separated prevalence dfs
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_prop_prevalence_df <- read.table("separated_prop_prevalence_df.csv")
separated_entropy_prevalence_df <- read.table("separated_entropy_prevalence_df.csv")

# Start with cluster_A
separated_A_prop_prevalence_plot <- plot_proportion_prevalence(separated_A_spes_table, separated_prop_prevalence_df, "distance")
separated_A_entropy_prevalence_plot <- plot_entropy_prevalence(separated_A_spes_table, separated_entropy_prevalence_df, "distance")
separated_A_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(separated_A_spes_table, separated_prop_prevalence_df, "distance")
separated_A_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(separated_A_spes_table, separated_entropy_prevalence_df, "distance")

# Then do cluster_B
separated_B_prop_prevalence_plot <- plot_proportion_prevalence(separated_B_spes_table, separated_prop_prevalence_df, "distance")
separated_B_entropy_prevalence_plot <- plot_entropy_prevalence(separated_B_spes_table, separated_entropy_prevalence_df, "distance")
separated_B_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(separated_B_spes_table, separated_prop_prevalence_df, "distance")
separated_B_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(separated_B_spes_table, separated_entropy_prevalence_df, "distance")

setwd("~/Objects/unsupervised/plots/original")
saveRDS(separated_A_prop_prevalence_plot, "separated_A_prop_prevalence_plot.RDS")
saveRDS(separated_A_entropy_prevalence_plot, "separated_A_entropy_prevalence_plot.RDS")
saveRDS(separated_A_prop_prevalence_AUC_plot, "separated_A_prop_prevalence_AUC_plot.RDS")
saveRDS(separated_A_entropy_prevalence_AUC_plot, "separated_A_entropy_prevalence_AUC_plot.RDS")

saveRDS(separated_B_prop_prevalence_plot, "separated_B_prop_prevalence_plot.RDS")
saveRDS(separated_B_entropy_prevalence_plot, "separated_B_entropy_prevalence_plot.RDS")
saveRDS(separated_B_prop_prevalence_AUC_plot, "separated_B_prop_prevalence_AUC_plot.RDS")
saveRDS(separated_B_entropy_prevalence_AUC_plot, "separated_B_entropy_prevalence_AUC_plot.RDS")


### 5.0. pdf with all plots -----------------------------
setwd("~/Objects/unsupervised/plots/original")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_prevalence_AUC", "entropy_SAC", "entropy_prevalence", "entropy_prevalence_AUC")
arrangements <- c("mixed", "ringed", "separated_A", "separated_B")

pdf("plots.pdf", width = 15, height = 10)

for (metric in metrics) {
  for (arrangement in arrangements) {
    plot_file_name <- paste(arrangement, "_", metric, "_plot.RDS", sep = "")
    print(readRDS(plot_file_name))
  }
}
dev.off()

### Without background noise cells functions ---------------------------------
plot_AMD_metric <- function(spes_table, AMD_df, arrangement_colname) {
  
  # AMD pairs are A/A, A/B, B/A, B/B
  AMD_pairs <- data.frame(cell1 = c("A", "A", "B", "B"),
                          cell2 = c("A", "B", "A", "B"))
  AMD_pairs$pair <- paste(AMD_pairs$cell1, AMD_pairs$cell2, sep = "/")
  
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(nrow(AMD_pairs))) {
    
    # Subset AMD_df for chosen pair
    plot_df <- AMD_df[AMD_df$reference == AMD_pairs[i, "cell1"] & AMD_df$target == AMD_pairs[i, "cell2"], ]
    
    # Combine spes_table and AMD_df
    plot_df <- cbind(spes_table, plot_df)

    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), AMD, col = shape)) +
      geom_point() +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AMD, col = variation_E)) +
      geom_point() +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AMD, col = volume_E)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), AMD, col = width_N)) +
      geom_point() +
      theme_bw()
    
    all_plots_list[[AMD_pairs[i, "pair"]]] <- list(shape = fig_shape + theme(legend.position="none"),
                                                   variation_E = fig_variation_E + theme(legend.position="none"),
                                                   volume_E = fig_volume_E + theme(legend.position="none"),
                                                   width_N = fig_width_N + theme(legend.position="none"))
  }
  
  # Get legends
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by pairs
  plots_pair_list <- list()
  
  for (i in seq(nrow(AMD_pairs))) {
    pair <- AMD_pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$shape, 
                       all_plots_list[[pair]]$variation_E,
                       all_plots_list[[pair]]$volume_E,
                       all_plots_list[[pair]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[pair]]))
    
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
                        legends,
                        nrow = 5, ncol = 1, 
                        rel_heights = c(1, 1, 1, 1, 0.5))
  
  methods::show(AMD_plot)
  
  return(AMD_plot)
}
plot_gradient_metrics_type1 <- function(spes_table, gradient_metric_df, metric, arrangement_colname) {
  
  # Constants
  cell_types <- c("A", "B") # Use A as reference, and B as target, and vice versa
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(length(cell_types))) {
    # Subset gradient_metric_df for current reference cell
    plot_df <- gradient_metric_df[gradient_metric_df$reference == cell_types[i], ]
    
    # Combine spes_table and gradient_metric_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , radii_colnames)
    
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df$variable <- unfactor(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    
    # Remove spheres
    plot_df <- plot_df[plot_df$shape != "Sphere", ]
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = variation_E)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = volume_E)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = spe, col = width_N)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    all_plots_list[[cell_types[i]]] <- list(arrangement = fig_arrangement + theme(legend.position = "none"),
                                            shape = fig_shape + theme(legend.position = "none"),
                                            variation_E = fig_variation_E + theme(legend.position = "none"),
                                            volume_E = fig_volume_E + theme(legend.position = "none"),
                                            width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  # Get legends
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_arrangement,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference cell type
  plots_ref_list <- list()
  
  for (i in seq(length(cell_types))) {
    reference_cell_type <- cell_types[i]
    
    plots <- plot_grid(all_plots_list[[reference_cell_type]]$arrangement,
                       all_plots_list[[reference_cell_type]]$shape, 
                       all_plots_list[[reference_cell_type]]$variation_E,
                       all_plots_list[[reference_cell_type]]$volume_E,
                       all_plots_list[[reference_cell_type]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[reference_cell_type]]))
    
    title <- ggdraw() + 
      draw_label(paste("Reference:", reference_cell_type), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_ref_list[[reference_cell_type]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plot <- plot_grid(plots_ref_list$A,
                             plots_ref_list$B, 
                             legends,
                             nrow = 3, ncol = 1,
                             rel_heights = c(1, 1, 0.5))
  
  methods::show(combined_plot)
  
  return(combined_plot)
}
plot_gradient_metrics_type2 <- function(spes_table, gradient_metric_df, metric, arrangement_colname, min_radius, max_radius) {
  
  # Constants
  pairs <- data.frame(cell1 = c("A", "A", "B", "B"),
                      cell2 = c("A", "B", "A", "B"))
  pairs$pair <- paste(pairs$cell1, pairs$cell2, sep = "/")
  
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  
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
    
    # Remove spheres
    plot_df <- plot_df[plot_df$shape != "Sphere", ]
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      labs(x = "radius", y = metric) +
      geom_line() +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = variation_E)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = volume_E)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = spe, col = width_N)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    all_plots_list[[pairs[i, "pair"]]] <- list(arrangement = fig_arrangement + theme(legend.position = "none"),
                                               shape = fig_shape + theme(legend.position = "none"),
                                               variation_E = fig_variation_E + theme(legend.position = "none"),
                                               volume_E = fig_volume_E + theme(legend.position = "none"),
                                               width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  
  # Get legends
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_arrangement,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference cell type
  plots_pair_list <- list()
  
  for (i in seq(nrow(pairs))) {
    pair <- pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$arrangement,
                       all_plots_list[[pair]]$shape, 
                       all_plots_list[[pair]]$variation_E,
                       all_plots_list[[pair]]$volume_E,
                       all_plots_list[[pair]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[pair]]))
    
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
                             legends,
                             nrow = 5, ncol = 1,
                             rel_heights = c(1, 1, 1, 1, 0.5))
  
  methods::show(combined_plot)
  
  return(combined_plot)
}
plot_proportion_SAC <- function(spes_table, SAC_df, arrangement_colname) {
  
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
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), proportion, col = shape)) +
      geom_point() +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), proportion, col = variation_E)) +
      geom_point() +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), proportion, col = volume_E)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), proportion, col = width_N)) +
      geom_point() +
      theme_bw()
    
    all_plots_list[[prop_cell_types[i, "pair"]]] <- list(shape = fig_shape + theme(legend.position="none"),
                                                         variation_E = fig_variation_E + theme(legend.position="none"),
                                                         volume_E = fig_volume_E + theme(legend.position="none"),
                                                         width_N = fig_width_N + theme(legend.position="none"))
  }
  
  
  # Get legends
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by pairs
  
  plots_pair_list <- list()
  
  for (i in seq(nrow(prop_cell_types))) {
    pair <- prop_cell_types[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$shape, 
                       all_plots_list[[pair]]$variation_E,
                       all_plots_list[[pair]]$volume_E,
                       all_plots_list[[pair]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[pair]]))
    
    title <- ggdraw() + 
      draw_label(paste("Reference/Target:", pair), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_pair_list[[pair]] <- fig
  }
  
  # Combine the combined plots into one big plot
  SAC_plot <- plot_grid(plots_pair_list[[prop_cell_types$pair[1]]], 
                        plots_pair_list[[prop_cell_types$pair[2]]],
                        legends,
                        nrow = 3, ncol = 1, 
                        rel_heights = c(1, 1, 0.5))
  
  methods::show(SAC_plot)
  
  return(SAC_plot)
}
plot_entropy_SAC <- function(spes_table, SAC_df, arrangement_colname) {
  
  # Get possible cell type of interest combinations
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    
    # Subset for current cell type of interest combination
    plot_df <- SAC_df[SAC_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Combine spes_table and SAC_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), entropy, col = shape)) +
      geom_point() +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), entropy, col = variation_E)) +
      geom_point() +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), entropy, col = volume_E)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), entropy, col = width_N)) +
      geom_point() +
      theme_bw()
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(shape = fig_shape + theme(legend.position="none"),
                                                               variation_E = fig_variation_E + theme(legend.position="none"),
                                                               volume_E = fig_volume_E + theme(legend.position="none"),
                                                               width_N = fig_width_N + theme(legend.position="none"))
  }
  
  # Get legends
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by cell types of interest
  
  plots_cell_types_list <- list()
  
  for (i in seq(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$shape, 
                       all_plots_list[[cell_types]]$variation_E,
                       all_plots_list[[cell_types]]$volume_E,
                       all_plots_list[[cell_types]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[cell_types]]))
    
    title <- ggdraw() + 
      draw_label(paste("Cell types of interest:", cell_types), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_cell_types_list[[cell_types]] <- fig
  }
  
  
  # Combine the combined plots into one big plot
  SAC_plot <- plot_grid(plots_cell_types_list[[entropy_cell_types$cell_types[1]]], 
                        plots_cell_types_list[[entropy_cell_types$cell_types[2]]],
                        legends,
                        nrow = 3, ncol = 1,
                        rel_heights = c(1, 1, 0.5))
  
  methods::show(SAC_plot)
  
  return(SAC_plot)
}
plot_proportion_prevalence <- function(spes_table, prevalence_df, arrangement_colname) {
  
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
    
    # Remove spheres
    plot_df <- plot_df[plot_df$shape != "Sphere", ]
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      labs(x = "threshold", y = "prevalence") +
      geom_line() +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = variation_E)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = volume_E)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = spe, col = width_N)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(arrangement = fig_arrangement + theme(legend.position = "none"),
                                                      shape = fig_shape + theme(legend.position = "none"),
                                                      variation_E = fig_variation_E + theme(legend.position = "none"),
                                                      volume_E = fig_volume_E + theme(legend.position = "none"),
                                                      width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  # Get legends
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_arrangement,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  
  # Combine the plots together by reference (and target) cell pairs
  plots_pair_list <- list()
  
  for (i in seq(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$arrangement,
                       all_plots_list[[pair]]$shape, 
                       all_plots_list[[pair]]$variation_E,
                       all_plots_list[[pair]]$volume_E,
                       all_plots_list[[pair]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[pair]]))
    
    title <- ggdraw() + 
      draw_label(paste("Reference/Target:", pair), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_pair_list[[pair]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plot <- plot_grid(plots_pair_list[[prop_cell_types$pair[1]]],
                             plots_pair_list[[prop_cell_types$pair[2]]], 
                             legends,
                             nrow = 3, ncol = 1,
                             rel_heights = c(1, 1, 0.5))
  
  methods::show(combined_plot)
  
  return(combined_plot)
}
plot_entropy_prevalence <- function(spes_table, prevalence_df, arrangement_colname) {
  
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
    
    # Remove spheres
    plot_df <- plot_df[plot_df$shape != "Sphere", ]
    
    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      labs(x = "threshold", y = "prevalence") +
      geom_line() +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = variation_E)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = volume_E)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = spe, col = width_N)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(arrangement = fig_arrangement + theme(legend.position = "none"), 
                                                               shape = fig_shape + theme(legend.position = "none"),
                                                               variation_E = fig_variation_E + theme(legend.position = "none"),
                                                               volume_E = fig_volume_E + theme(legend.position = "none"),
                                                               width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  # Get legends
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_arrangement,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  
  # Combine the plots together by reference (and target) cell pairs
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$arrangement,
                       all_plots_list[[cell_types]]$shape, 
                       all_plots_list[[cell_types]]$variation_E,
                       all_plots_list[[cell_types]]$volume_E,
                       all_plots_list[[cell_types]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[cell_types]]))
    
    title <- ggdraw() + 
      draw_label(paste("Cell types of interest:", cell_types), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_cell_types_list[[cell_types]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plot <- plot_grid(plots_cell_types_list[[entropy_cell_types$cell_types[1]]],
                             plots_cell_types_list[[entropy_cell_types$cell_types[2]]], 
                             legends,
                             nrow = 3, ncol = 1,
                             rel_heights = c(1, 1, 0.5))
  
  methods::show(combined_plot)
  
  return(combined_plot)
}
plot_proportion_prevalence_AUC <- function(spes_table, prevalence_df, arrangement_colname) {
  
  # Constants
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
  prop_cell_types$pair <- paste(prop_cell_types$ref, prop_cell_types$tar, sep = "/")
  
  # Get AUC for each prevalence gradient
  prevalence_df$AUC <- apply(prevalence_df[ , threshold_colnames], 1, sum) * 0.01
  prevalence_df <- prevalence_df[ , c("spe", "reference", "target", "AUC")]
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    
    # Subset for current reference (and target) cell
    plot_df <- prevalence_df[prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Combine spes_table and updated prevalence_df
    plot_df <- cbind(spes_table, plot_df)

    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = shape)) +
      geom_point() +
      ylab("AUC") +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AUC, col = variation_E)) +
      geom_point() +
      ylab("AUC") +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AUC, col = volume_E)) +
      geom_point() +
      ylab("AUC") +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), AUC, col = width_N)) +
      geom_point() +
      ylab("AUC") +
      theme_bw()
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(shape = fig_shape + theme(legend.position="none"),
                                                      variation_E = fig_variation_E + theme(legend.position="none"),
                                                      volume_E = fig_volume_E + theme(legend.position="none"),
                                                      width_N = fig_width_N + theme(legend.position="none"))
  }
  
  # Get legends
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference target pairs
  plots_pair_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$shape, 
                       all_plots_list[[pair]]$variation_E,
                       all_plots_list[[pair]]$volume_E,
                       all_plots_list[[pair]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[pair]]))
    
    title <- ggdraw() + 
      draw_label(paste("Reference/Target:", pair), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_pair_list[[pair]] <- fig
  }
  
  # Combine the combined plots into one big plot
  AUC_plot <- plot_grid(plots_pair_list[[prop_cell_types$pair[1]]], 
                        plots_pair_list[[prop_cell_types$pair[2]]],    
                        legends,
                        nrow = 3, ncol = 1, 
                        rel_heights = c(1, 1, 0.5))
  
  methods::show(AUC_plot)
  
  return(AUC_plot)
}
plot_entropy_prevalence_AUC <- function(spes_table, prevalence_df, arrangement_colname) {
  
  # Constants
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  
  # Get AUC for each prevalence gradient
  prevalence_df$AUC <- apply(prevalence_df[ , threshold_colnames], 1, sum) * 0.01
  prevalence_df <- prevalence_df[ , c("spe", "cell_types", "AUC")]
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    
    # Subset for current cell type of interest combintation
    plot_df <- prevalence_df[prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Combine spes_table and updated prevalence_df
    plot_df <- cbind(spes_table, plot_df)

    # Factor
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = shape)) +
      geom_point() +
      ylab("AUC") +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    
    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AUC, col = variation_E)) +
      geom_point() +
      ylab("AUC") +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AUC, col = volume_E)) +
      geom_point() +
      ylab("AUC") +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)
    
    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), AUC, col = width_N)) +
      geom_point() +
      ylab("AUC") +
      theme_bw()
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(shape = fig_shape + theme(legend.position="none"),
                                                               variation_E = fig_variation_E + theme(legend.position="none"),
                                                               volume_E = fig_volume_E + theme(legend.position="none"),
                                                               width_N = fig_width_N + theme(legend.position="none"))
  }
  
  # Get legends
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference target pairs
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$shape, 
                       all_plots_list[[cell_types]]$variation_E,
                       all_plots_list[[cell_types]]$volume_E,
                       all_plots_list[[cell_types]]$width_N, 
                       nrow = 1, ncol = length(all_plots_list[[cell_types]]))
    
    title <- ggdraw() + 
      draw_label(paste("Cell types of interest:", cell_types), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_cell_types_list[[cell_types]] <- fig
  }
  
  # Combine the combined plots into one big plot
  AUC_plot <- plot_grid(plots_cell_types_list[[entropy_cell_types$cell_types[1]]], 
                        plots_cell_types_list[[entropy_cell_types$cell_types[2]]],    
                        legends,
                        nrow = 3, ncol = 1, 
                        rel_heights = c(1, 1, 0.5))
  
  methods::show(AUC_plot)
  
  return(AUC_plot)
}

### Without background noise cells set up ---------------------------------------------
mixed_spes_table <- mixed_spes_table[mixed_spes_table$bg_prop_A == 0 & mixed_spes_table$bg_prop_B == 0, ]
ringed_spes_table <- ringed_spes_table[ringed_spes_table$bg_prop_A == 0 & ringed_spes_table$bg_prop_B == 0, ]
separated_spes_table <- separated_spes_table[separated_spes_table$bg_prop_A == 0 & separated_spes_table$bg_prop_B == 0, ]
separated_A_spes_table <- separated_A_spes_table[separated_A_spes_table$bg_prop_A == 0 & separated_A_spes_table$bg_prop_B == 0, ]
separated_B_spes_table <- separated_B_spes_table[separated_B_spes_table$bg_prop_A == 0 & separated_B_spes_table$bg_prop_B == 0, ]

### Without background noise cells analysis ---------------------------------------------

# Read mixed_AMD_df
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_AMD_df <- read.table("mixed_AMD_df.csv")
mixed_AMD_df <- mixed_AMD_df[mixed_AMD_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]

mixed_AMD_plot <- plot_AMD_metric(mixed_spes_table, mixed_AMD_df, "cluster_prop_B")

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(mixed_AMD_plot, "mixed_AMD_plot_no_background_noise.RDS")



# Read mixed MS, NMS, ACINP, AE dfs
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_MS_df <- read.table("mixed_MS_df.csv")
mixed_NMS_df <- read.table("mixed_NMS_df.csv")
mixed_ACINP_df <- read.table("mixed_ACINP_df.csv")
mixed_AE_df <- read.table("mixed_AE_df.csv")

mixed_MS_df <- mixed_MS_df[mixed_MS_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_NMS_df <- mixed_NMS_df[mixed_NMS_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_ACINP_df <- mixed_ACINP_df[mixed_ACINP_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_AE_df <- mixed_AE_df[mixed_AE_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]

mixed_MS_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_MS_df, "MS", "cluster_prop_B")
mixed_NMS_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_NMS_df, "NMS", "cluster_prop_B")
mixed_ACINP_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_ACINP_df, "ACINP", "cluster_prop_B")
mixed_AE_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_AE_df, "AE", "cluster_prop_B")

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(mixed_MS_plot, "mixed_MS_plot_no_background_noise.RDS")
saveRDS(mixed_NMS_plot, "mixed_NMS_plot_no_background_noise.RDS")
saveRDS(mixed_ACINP_plot, "mixed_ACINP_plot_no_background_noise.RDS")
saveRDS(mixed_AE_plot, "mixed_AE_plot_no_background_noise.RDS")



# Read mixed ACIN, CKR
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_ACIN_df <- read.table("mixed_ACIN_df.csv")
mixed_CKR_df <- read.table("mixed_CKR_df.csv")

mixed_ACIN_df <- mixed_ACIN_df[mixed_ACIN_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_CKR_df <- mixed_CKR_df[mixed_CKR_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]

# Get plots
mixed_ACIN_plot <- plot_gradient_metrics_type2(mixed_spes_table, mixed_ACIN_df, "ACIN", "cluster_prop_B", 20, 100)
mixed_CKR_plot <- plot_gradient_metrics_type2(mixed_spes_table, mixed_CKR_df, "CKR", "cluster_prop_B", 20, 100)

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(mixed_ACIN_plot, "mixed_ACIN_plot_no_background_noise.RDS")
saveRDS(mixed_CKR_plot, "mixed_CKR_plot_no_background_noise.RDS")



# Read mixed_SAC_df
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_prop_SAC_df <- read.table("mixed_prop_SAC_df.csv")
mixed_entropy_SAC_df <- read.table("mixed_entropy_SAC_df.csv")

mixed_prop_SAC_df <- mixed_prop_SAC_df[mixed_prop_SAC_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_entropy_SAC_df <- mixed_entropy_SAC_df[mixed_entropy_SAC_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]

mixed_prop_SAC_plot <- plot_proportion_SAC(mixed_spes_table, mixed_prop_SAC_df, "cluster_prop_B")
mixed_entropy_SAC_plot <- plot_entropy_SAC(mixed_spes_table, mixed_entropy_SAC_df, "cluster_prop_B")

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(mixed_prop_SAC_plot, "mixed_prop_SAC_plot_no_background_noise.RDS")
saveRDS(mixed_entropy_SAC_plot, "mixed_entropy_SAC_plot_no_background_noise.RDS")



# Read mixed prevalence dfs
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_prop_prevalence_df <- read.table("mixed_prop_prevalence_df.csv")
mixed_entropy_prevalence_df <- read.table("mixed_entropy_prevalence_df.csv")

mixed_prop_prevalence_df <- mixed_prop_prevalence_df[mixed_prop_prevalence_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_entropy_prevalence_df <- mixed_entropy_prevalence_df[mixed_entropy_prevalence_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]

mixed_prop_prevalence_plot <- plot_proportion_prevalence(mixed_spes_table, mixed_prop_prevalence_df, "cluster_prop_B")
mixed_entropy_prevalence_plot <- plot_entropy_prevalence(mixed_spes_table, mixed_entropy_prevalence_df, "cluster_prop_B")
mixed_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(mixed_spes_table, mixed_prop_prevalence_df, "cluster_prop_B")
mixed_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(mixed_spes_table, mixed_entropy_prevalence_df, "cluster_prop_B")

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(mixed_prop_prevalence_plot, "mixed_prop_prevalence_plot_no_background_noise.RDS")
saveRDS(mixed_entropy_prevalence_plot, "mixed_entropy_prevalence_plot_no_background_noise.RDS")
saveRDS(mixed_prop_prevalence_AUC_plot, "mixed_prop_prevalence_AUC_plot_no_background_noise.RDS")
saveRDS(mixed_entropy_prevalence_AUC_plot, "mixed_entropy_prevalence_AUC_plot_no_background_noise.RDS")





# Read ringed_AMD_df
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_AMD_df <- read.table("ringed_AMD_df.csv")
ringed_AMD_df <- ringed_AMD_df[ringed_AMD_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]

ringed_AMD_plot <- plot_AMD_metric(ringed_spes_table, ringed_AMD_df, "width_ring_factor")

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(ringed_AMD_plot, "ringed_AMD_plot_no_background_noise.RDS")



# Read ringed MS, NMS, ACINP, AE dfs
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_MS_df <- read.table("ringed_MS_df.csv")
ringed_NMS_df <- read.table("ringed_NMS_df.csv")
ringed_ACINP_df <- read.table("ringed_ACINP_df.csv")
ringed_AE_df <- read.table("ringed_AE_df.csv")

ringed_MS_df <- ringed_MS_df[ringed_MS_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_NMS_df <- ringed_NMS_df[ringed_NMS_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_ACINP_df <- ringed_ACINP_df[ringed_ACINP_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_AE_df <- ringed_AE_df[ringed_AE_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]

ringed_MS_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_MS_df, "MS", "width_ring_factor")
ringed_NMS_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_NMS_df, "NMS", "width_ring_factor")
ringed_ACINP_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_ACINP_df, "ACINP", "width_ring_factor")
ringed_AE_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_AE_df, "AE", "width_ring_factor")

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(ringed_MS_plot, "ringed_MS_plot_no_background_noise.RDS")
saveRDS(ringed_NMS_plot, "ringed_NMS_plot_no_background_noise.RDS")
saveRDS(ringed_ACINP_plot, "ringed_ACINP_plot_no_background_noise.RDS")
saveRDS(ringed_AE_plot, "ringed_AE_plot_no_background_noise.RDS")


# Read ringed ACIN, CKR
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_ACIN_df <- read.table("ringed_ACIN_df.csv")
ringed_CKR_df <- read.table("ringed_CKR_df.csv")

ringed_ACIN_df <- ringed_ACIN_df[ringed_ACIN_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_CKR_df <- ringed_CKR_df[ringed_CKR_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]

# Get plots
ringed_ACIN_plot <- plot_gradient_metrics_type2(ringed_spes_table, ringed_ACIN_df, "ACIN", "width_ring_factor", 20, 100)
ringed_CKR_plot <- plot_gradient_metrics_type2(ringed_spes_table, ringed_CKR_df, "CKR", "width_ring_factor", 20, 100)

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(ringed_ACIN_plot, "ringed_ACIN_plot_no_background_noise.RDS")
saveRDS(ringed_CKR_plot, "ringed_CKR_plot_no_background_noise.RDS")



# Read ringed_SAC_df
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_prop_SAC_df <- read.table("ringed_prop_SAC_df.csv")
ringed_entropy_SAC_df <- read.table("ringed_entropy_SAC_df.csv")

ringed_prop_SAC_df <- ringed_prop_SAC_df[ringed_prop_SAC_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_entropy_SAC_df <- ringed_entropy_SAC_df[ringed_entropy_SAC_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]

ringed_prop_SAC_plot <- plot_proportion_SAC(ringed_spes_table, ringed_prop_SAC_df, "width_ring_factor")
ringed_entropy_SAC_plot <- plot_entropy_SAC(ringed_spes_table, ringed_entropy_SAC_df, "width_ring_factor")

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(ringed_prop_SAC_plot, "ringed_prop_SAC_plot_no_background_noise.RDS")
saveRDS(ringed_entropy_SAC_plot, "ringed_entropy_SAC_plot_no_background_noise.RDS")


# Read ringed prevalence dfs
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_prop_prevalence_df <- read.table("ringed_prop_prevalence_df.csv")
ringed_entropy_prevalence_df <- read.table("ringed_entropy_prevalence_df.csv")

ringed_prop_prevalence_df <- ringed_prop_prevalence_df[ringed_prop_prevalence_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_entropy_prevalence_df <- ringed_entropy_prevalence_df[ringed_entropy_prevalence_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]

ringed_prop_prevalence_plot <- plot_proportion_prevalence(ringed_spes_table, ringed_prop_prevalence_df, "width_ring_factor")
ringed_entropy_prevalence_plot <- plot_entropy_prevalence(ringed_spes_table, ringed_entropy_prevalence_df, "width_ring_factor")
ringed_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(ringed_spes_table, ringed_prop_prevalence_df, "width_ring_factor")
ringed_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(ringed_spes_table, ringed_entropy_prevalence_df, "width_ring_factor")

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(ringed_prop_prevalence_plot, "ringed_prop_prevalence_plot_no_background_noise.RDS")
saveRDS(ringed_entropy_prevalence_plot, "ringed_entropy_prevalence_plot_no_background_noise.RDS")
saveRDS(ringed_prop_prevalence_AUC_plot, "ringed_prop_prevalence_AUC_plot_no_background_noise.RDS")
saveRDS(ringed_entropy_prevalence_AUC_plot, "ringed_entropy_prevalence_AUC_plot_no_background_noise.RDS")




# Read separated_AMD_df
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_AMD_df <- read.table("separated_AMD_df.csv")
separated_AMD_df <- separated_AMD_df[separated_AMD_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]


separated_A_AMD_plot <- plot_AMD_metric(separated_A_spes_table, separated_AMD_df, "distance")
separated_B_AMD_plot <- plot_AMD_metric(separated_B_spes_table, separated_AMD_df, "distance")

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(separated_A_AMD_plot, "separated_A_AMD_plot_no_background_noise.RDS")
saveRDS(separated_B_AMD_plot, "separated_B_AMD_plot_no_background_noise.RDS")


# Read separated MS, NMS, ACINP, AE dfs
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_MS_df <- read.table("separated_MS_df.csv")
separated_NMS_df <- read.table("separated_NMS_df.csv")
separated_ACINP_df <- read.table("separated_ACINP_df.csv")
separated_AE_df <- read.table("separated_AE_df.csv")

separated_MS_df <- separated_MS_df[separated_MS_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_NMS_df <- separated_NMS_df[separated_NMS_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_ACINP_df <- separated_ACINP_df[separated_ACINP_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_AE_df <- separated_AE_df[separated_AE_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]

# Start with cluster_A
separated_A_MS_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_MS_df, "MS", "distance")
separated_A_NMS_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_NMS_df, "NMS", "distance")
separated_A_ACINP_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_ACINP_df, "ACINP", "distance")
separated_A_AE_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_AE_df, "AE", "distance")

# Then do cluster_B
separated_B_MS_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_MS_df, "MS", "distance")
separated_B_NMS_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_NMS_df, "NMS", "distance")
separated_B_ACINP_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_ACINP_df, "ACINP", "distance")
separated_B_AE_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_AE_df, "AE", "distance")

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(separated_A_MS_plot, "separated_A_MS_plot_no_background_noise.RDS")
saveRDS(separated_A_NMS_plot, "separated_A_NMS_plot_no_background_noise.RDS")
saveRDS(separated_A_ACINP_plot, "separated_A_ACINP_plot_no_background_noise.RDS")
saveRDS(separated_A_AE_plot, "separated_A_AE_plot_no_background_noise.RDS")

saveRDS(separated_B_MS_plot, "separated_B_MS_plot_no_background_noise.RDS")
saveRDS(separated_B_NMS_plot, "separated_B_NMS_plot_no_background_noise.RDS")
saveRDS(separated_B_ACINP_plot, "separated_B_ACINP_plot_no_background_noise.RDS")
saveRDS(separated_B_AE_plot, "separated_B_AE_plot_no_background_noise.RDS")


# Read separated ACIN, CKR
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_ACIN_df <- read.table("separated_ACIN_df.csv")
separated_CKR_df <- read.table("separated_CKR_df.csv")

separated_ACIN_df <- separated_ACIN_df[separated_ACIN_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_CKR_df <- separated_CKR_df[separated_CKR_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]

# Start with cluster_A
separated_A_ACIN_plot <- plot_gradient_metrics_type2(separated_A_spes_table, separated_ACIN_df, "ACIN", "distance", 20, 100)
separated_A_CKR_plot <- plot_gradient_metrics_type2(separated_A_spes_table, separated_CKR_df, "CKR", "distance", 20, 100)

# Then do cluster_B
separated_B_ACIN_plot <- plot_gradient_metrics_type2(separated_B_spes_table, separated_ACIN_df, "ACIN", "distance", 20, 100)
separated_B_CKR_plot <- plot_gradient_metrics_type2(separated_B_spes_table, separated_CKR_df, "CKR", "distance", 20, 100)

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(separated_A_ACIN_plot, "separated_A_ACIN_plot_no_background_noise.RDS")
saveRDS(separated_A_CKR_plot, "separated_A_CKR_plot_no_background_noise.RDS")

saveRDS(separated_B_ACIN_plot, "separated_B_ACIN_plot_no_background_noise.RDS")
saveRDS(separated_B_CKR_plot, "separated_B_CKR_plot_no_background_noise.RDS")


# Read separated_SAC_df
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_prop_SAC_df <- read.table("separated_prop_SAC_df.csv")
separated_entropy_SAC_df <- read.table("separated_entropy_SAC_df.csv")

separated_prop_SAC_df <- separated_prop_SAC_df[separated_prop_SAC_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_entropy_SAC_df <- separated_entropy_SAC_df[separated_entropy_SAC_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]

# Start with cluster_A
separated_A_prop_SAC_plot <- plot_proportion_SAC(separated_A_spes_table, separated_prop_SAC_df, "distance")
separated_A_entropy_SAC_plot <- plot_entropy_SAC(separated_A_spes_table, separated_entropy_SAC_df, "distance")

# Then do cluster_B
separated_B_prop_SAC_plot <- plot_proportion_SAC(separated_B_spes_table, separated_prop_SAC_df, "distance")
separated_B_entropy_SAC_plot <- plot_entropy_SAC(separated_B_spes_table, separated_entropy_SAC_df, "distance")

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(separated_A_prop_SAC_plot, "separated_A_prop_SAC_plot_no_background_noise.RDS")
saveRDS(separated_A_entropy_SAC_plot, "separated_A_entropy_SAC_plot_no_background_noise.RDS")

saveRDS(separated_B_prop_SAC_plot, "separated_B_prop_SAC_plot_no_background_noise.RDS")
saveRDS(separated_B_entropy_SAC_plot, "separated_B_entropy_SAC_plot_no_background_noise.RDS")


# Read separated prevalence dfs
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_prop_prevalence_df <- read.table("separated_prop_prevalence_df.csv")
separated_entropy_prevalence_df <- read.table("separated_entropy_prevalence_df.csv")

separated_prop_prevalence_df <- separated_prop_prevalence_df[separated_prop_prevalence_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_entropy_prevalence_df <- separated_entropy_prevalence_df[separated_entropy_prevalence_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]

# Start with cluster_A
separated_A_prop_prevalence_plot <- plot_proportion_prevalence(separated_A_spes_table, separated_prop_prevalence_df, "distance")
separated_A_entropy_prevalence_plot <- plot_entropy_prevalence(separated_A_spes_table, separated_entropy_prevalence_df, "distance")
separated_A_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(separated_A_spes_table, separated_prop_prevalence_df, "distance")
separated_A_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(separated_A_spes_table, separated_entropy_prevalence_df, "distance")

# Then do cluster_B
separated_B_prop_prevalence_plot <- plot_proportion_prevalence(separated_B_spes_table, separated_prop_prevalence_df, "distance")
separated_B_entropy_prevalence_plot <- plot_entropy_prevalence(separated_B_spes_table, separated_entropy_prevalence_df, "distance")
separated_B_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(separated_B_spes_table, separated_prop_prevalence_df, "distance")
separated_B_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(separated_B_spes_table, separated_entropy_prevalence_df, "distance")

setwd("~/Objects/unsupervised/plots/original_no_background_noise")
saveRDS(separated_A_prop_prevalence_plot, "separated_A_prop_prevalence_plot_no_background_noise.RDS")
saveRDS(separated_A_entropy_prevalence_plot, "separated_A_entropy_prevalence_plot_no_background_noise.RDS")
saveRDS(separated_A_prop_prevalence_AUC_plot, "separated_A_prop_prevalence_AUC_plot_no_background_noise.RDS")
saveRDS(separated_A_entropy_prevalence_AUC_plot, "separated_A_entropy_prevalence_AUC_plot_no_background_noise.RDS")

saveRDS(separated_B_prop_prevalence_plot, "separated_B_prop_prevalence_plot_no_background_noise.RDS")
saveRDS(separated_B_entropy_prevalence_plot, "separated_B_entropy_prevalence_plot_no_background_noise.RDS")
saveRDS(separated_B_prop_prevalence_AUC_plot, "separated_B_prop_prevalence_AUC_plot_no_background_noise.RDS")
saveRDS(separated_B_entropy_prevalence_AUC_plot, "separated_B_entropy_prevalence_AUC_plot_no_background_noise.RDS")




### Without background noise pdf with all plots -------------------------
setwd("~/Objects/unsupervised/plots/original_no_background_noise")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_prevalence_AUC", "entropy_SAC", "entropy_prevalence", "entropy_prevalence_AUC")
arrangements <- c("mixed", "ringed", "separated_A", "separated_B")

pdf("plots_no_background_noise.pdf", width = 15, height = 10)

for (metric in metrics) {
  for (arrangement in arrangements) {
    plot_file_name <- paste(arrangement, "_", metric, "_plot_no_background_noise.RDS", sep = "")
    print(readRDS(plot_file_name))
  }
}
dev.off()
