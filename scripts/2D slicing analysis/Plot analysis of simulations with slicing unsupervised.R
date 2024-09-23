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
plot_AMD_metric <- function(spes_table, AMD_df, slices_AMD_df, arrangement_colname) {
  
  # AMD pairs are A/A, A/B, B/A, B/B
  AMD_pairs <- data.frame(cell1 = c("A", "A", "B", "B"),
                          cell2 = c("A", "B", "A", "B"))
  AMD_pairs$pair <- paste(AMD_pairs$cell1, AMD_pairs$cell2, sep = "/")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(nrow(AMD_pairs))) {
    
    # Subset AMD_df for chosen pair
    AMD_df_subset <- AMD_df[AMD_df$reference == AMD_pairs[i, "cell1"] & AMD_df$target == AMD_pairs[i, "cell2"], ]
    
    # Subset slices_AMD_df for chosen pair
    slices_AMD_df_subset <- slices_AMD_df[slices_AMD_df$reference == AMD_pairs[i, "cell1"] & slices_AMD_df$target == AMD_pairs[i, "cell2"], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_AMD_df_subset, AMD_df_subset, "spe", suffix = c("_2D", "_3D"))
  
    slices_AMD_df_subset$AMD <- (joint_df$AMD_2D - joint_df$AMD_3D) / joint_df$AMD_3D
    
    # Combine spes_table and AMD_df
    plot_df <- cbind(spes_table, slices_AMD_df_subset)

    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)

    fig_slice <- ggplot(plot_df, aes(!!sym(arrangement_colname), AMD, col = slice)) +
      geom_point() +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))

    fig_bg_prop_A <- ggplot(plot_df, aes(!!sym(arrangement_colname), AMD, col = bg_prop_A)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))

    fig_bg_prop_B <- ggplot(plot_df, aes(!!sym(arrangement_colname), AMD, col = bg_prop_B)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))

    # fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), AMD, col = shape)) +
    #   geom_point() +
    #   theme_bw()

    # radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    # plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    # plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    # 
    # fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AMD, col = variation_E)) +
    #   geom_point() +
    #   theme_bw()
    
    # fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AMD, col = volume_E)) +
    #   geom_point() +
    #   theme_bw() +
    #   scale_color_continuous(n.breaks = 4)
    # 
    # fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), AMD, col = width_N)) +
    #   geom_point() +
    #   theme_bw()
    
    all_plots_list[[AMD_pairs[i, "pair"]]] <- list(slice = fig_slice + theme(legend.position="none"),
                                                   bg_prop_A = fig_bg_prop_A + theme(legend.position="none"),
                                                   bg_prop_B = fig_bg_prop_B + theme(legend.position="none"))
                                                   # shape = fig_shape + theme(legend.position="none")
                                                   # variation_E = fig_variation_E + theme(legend.position="none"),
                                                   # volume_E = fig_volume_E + theme(legend.position="none"),
                                                   # width_N = fig_width_N + theme(legend.position="none")
  }
  
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  # legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  # legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  # legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  # legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       legend_bg_prop_a,
                       legend_bg_prop_B,
                       # legend_shape,
                       # legend_variation_E,
                       # legend_volume_E,
                       # legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by pairs
  plots_pair_list <- list()
  
  for (i in seq(nrow(AMD_pairs))) {
    pair <- AMD_pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$slice,
                       all_plots_list[[pair]]$bg_prop_A,
                       all_plots_list[[pair]]$bg_prop_B,
                       # all_plots_list[[pair]]$shape,
                       # all_plots_list[[pair]]$variation_E,
                       # all_plots_list[[pair]]$volume_E,
                       # all_plots_list[[pair]]$width_N,
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

plot_AMD_density_metric <- function(spes_table, AMD_df, slices_AMD_df, arrangement_colname) {
  
  # AMD pairs are A/A, A/B, B/A, B/B
  AMD_pairs <- data.frame(cell1 = c("A", "A", "B", "B"),
                          cell2 = c("A", "B", "A", "B"))
  AMD_pairs$pair <- paste(AMD_pairs$cell1, AMD_pairs$cell2, sep = "/")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(nrow(AMD_pairs))) {
    
    # Subset AMD_df for chosen pair
    AMD_df_subset <- AMD_df[AMD_df$reference == AMD_pairs[i, "cell1"] & AMD_df$target == AMD_pairs[i, "cell2"], ]
    
    # Subset slices_AMD_df for chosen pair
    slices_AMD_df_subset <- slices_AMD_df[slices_AMD_df$reference == AMD_pairs[i, "cell1"] & slices_AMD_df$target == AMD_pairs[i, "cell2"], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_AMD_df_subset, AMD_df_subset, "spe", suffix = c("_2D", "_3D"))
    
    slices_AMD_df_subset$AMD <- NULL
    slices_AMD_df_subset$AMD_truth <- joint_df$AMD_3D
    slices_AMD_df_subset$AMD_error <- (joint_df$AMD_2D - joint_df$AMD_3D) / joint_df$AMD_3D
    
    # Combine spes_table and AMD_df
    plot_df <- cbind(spes_table, slices_AMD_df_subset)
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    
    fig_density <- ggplot(plot_df, aes(AMD_truth, AMD_error)) +
      geom_point() +
      theme_bw()
    
    all_plots_list[[AMD_pairs[i, "pair"]]] <- list(density = fig_density)

  }
  
  
  # Combine the plots together by pairs
  plots_pair_list <- list()
  
  for (i in seq(nrow(AMD_pairs))) {
    pair <- AMD_pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$density,
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
                        nrow = 5, ncol = 1, 
                        rel_heights = c(1, 1, 1, 1, 0.5))
  
  methods::show(AMD_plot)
  
  return(AMD_plot)
}


### 1.2.1. Function to get plot for MS, NMS, ACINP, AE gradient metrics ------------
plot_gradient_metrics_type1 <- function(spes_table, gradient_metric_df, slices_gradient_metric_df, metric, arrangement_colname) {
  
  # Constants
  cell_types <- c("A", "B") # Use A as reference, and B as target, and vice versa
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(length(cell_types))) {
    # Subset gradient_metric_df for current reference cell
    gradient_metric_df_subset <- gradient_metric_df[gradient_metric_df$reference == cell_types[i], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_gradient_metric_df_subset <- slices_gradient_metric_df[slices_gradient_metric_df$reference == cell_types[i], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_gradient_metric_df_subset, gradient_metric_df_subset, by = "spe", suffix = c("_2D", "_3D"))

    for (radii_colname in radii_colnames) {
      slices_gradient_metric_df_subset[ , radii_colname] <- 
        (joint_df[ , paste(radii_colname, "_2D", sep = "")] - joint_df[ , paste(radii_colname, "_3D", sep = "")]) / joint_df[ , paste(radii_colname, "_3D", sep = "")]
    }

    # Combine spes_table and slices_gradient_metric_df_subset
    plot_df <- cbind(spes_table, slices_gradient_metric_df_subset)

    # Melt
    plot_df <- reshape2::melt(plot_df, , radii_colnames)

    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    plot_df$key <- paste(plot_df$spe, plot_df$slice, sep = "_")
    
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df$variable <- unfactor(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    
    fig_slice <- ggplot(plot_df, aes(variable, value, group = key, col = slice)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = key, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    fig_bg_prop_A <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_A)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))

    fig_bg_prop_B <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_B)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))

    # fig_shape <- ggplot(plot_df, aes(variable, value, group = key, col = shape)) +
    #   geom_line() +
    #   labs(x = "radius", y = metric) +
    #   theme_bw()

    # radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    # plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    # plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    # 
    # fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = variation_E)) +
    #   geom_line() +
    #   labs(x = "radius", y = metric) +
    #   theme_bw()
    
    # fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = volume_E)) +
    #   geom_line() +
    #   labs(x = "radius", y = metric) +
    #   theme_bw() +
    #   scale_color_continuous(n.breaks = 4)
    # 
    # fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = spe, col = width_N)) +
    #   geom_line() +
    #   labs(x = "radius", y = metric) +
    #   theme_bw()
    
    all_plots_list[[cell_types[i]]] <- list(slice = fig_slice + theme(legend.position = "none"),
                                            arrangement = fig_arrangement + theme(legend.position = "none"), 
                                            bg_prop_A = fig_bg_prop_A + theme(legend.position = "none"),
                                            bg_prop_B = fig_bg_prop_B + theme(legend.position = "none"))
                                            # shape = fig_shape + theme(legend.position = "none")
                                            # variation_E = fig_variation_E + theme(legend.position = "none"),
                                            # volume_E = fig_volume_E + theme(legend.position = "none"),
                                            # width_N = fig_width_N + theme(legend.position = "none")
                                            
    
  }
  
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  # legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  # # legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  # legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  # legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       legend_arrangement,
                       legend_bg_prop_a,
                       legend_bg_prop_B,
                       # legend_shape,
                       # legend_variation_E,
                       # legend_volume_E,
                       # legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference cell type
  plots_ref_list <- list()
  
  for (i in seq(length(cell_types))) {
    
    reference_cell_type <- cell_types[i]
    
    plots <- plot_grid(all_plots_list[[reference_cell_type]]$slice,
                       all_plots_list[[reference_cell_type]]$arrangement,
                       all_plots_list[[reference_cell_type]]$bg_prop_A,
                       all_plots_list[[reference_cell_type]]$bg_prop_B,
                       # all_plots_list[[reference_cell_type]]$shape,
                       # all_plots_list[[reference_cell_type]]$variation_E,
                       # all_plots_list[[reference_cell_type]]$volume_E,
                       # all_plots_list[[reference_cell_type]]$width_N, 
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

plot_gradient_metrics_type1_boxplot <- function(spes_table, gradient_metric_df, slices_gradient_metric_df, metric, arrangement_colname) {
  
  # Constants
  cell_types <- c("A", "B") # Use A as reference, and B as target, and vice versa
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(length(cell_types))) {
    # Subset gradient_metric_df for current reference cell
    gradient_metric_df_subset <- gradient_metric_df[gradient_metric_df$reference == cell_types[i], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_gradient_metric_df_subset <- slices_gradient_metric_df[slices_gradient_metric_df$reference == cell_types[i], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_gradient_metric_df_subset, gradient_metric_df_subset, by = "spe", suffix = c("_2D", "_3D"))
    
    for (radii_colname in radii_colnames) {
      slices_gradient_metric_df_subset[ , radii_colname] <- 
        (joint_df[ , paste(radii_colname, "_2D", sep = "")] - joint_df[ , paste(radii_colname, "_3D", sep = "")]) / joint_df[ , paste(radii_colname, "_3D", sep = "")]
    }
    
    # Combine spes_table and slices_gradient_metric_df_subset
    plot_df <- cbind(spes_table, slices_gradient_metric_df_subset)
    
    plot_df$slice <- as.character(plot_df$slice)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , radii_colnames)
    
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df$variable <- unfactor(plot_df$variable)
    plot_df$variable <- substr(plot_df$variable, 2, nchar(plot_df$variable))
    plot_df$variable <- factor(plot_df$variable, as.character(radii))
    
    fig_boxplot <- ggplot(plot_df, aes(variable, value, col = slice)) +
      geom_boxplot() +
      theme_bw() +
      facet_wrap(~slice, ncol = 5) +
      scale_color_manual(values = viridis::viridis(5)) +
      labs(x = "radius", y = metric)
    
    all_plots_list[[cell_types[i]]] <- list(fig_boxplot = fig_boxplot)
  }
  
  # Combine the plots together by reference cell type
  plots_ref_list <- list()
  
  for (i in seq(length(cell_types))) {
    reference_cell_type <- cell_types[i]
    
    plots <- plot_grid(all_plots_list[[reference_cell_type]]$fig_boxplot,
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
                             nrow = 2, ncol = 1)
  
  methods::show(combined_plot)
  
  return(combined_plot)
}


### 1.2.2. Function to get plot for ACIN, CKR gradient metrics ------------------
plot_gradient_metrics_type2 <- function(spes_table, gradient_metric_df, slices_gradient_metric_df, metric, arrangement_colname, min_radius, max_radius) {
  
  # Constants
  pairs <- data.frame(cell1 = c("A", "A", "B", "B"),
                      cell2 = c("A", "B", "A", "B"))
  pairs$pair <- paste(pairs$cell1, pairs$cell2, sep = "/")
  
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(nrow(pairs))) {

    # Subset gradient_metric_df for current reference cell
    gradient_metric_df_subset <- gradient_metric_df[gradient_metric_df$reference == pairs[i, "cell1"] &
                                                      gradient_metric_df$target == pairs[i, "cell2"], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_gradient_metric_df_subset <- slices_gradient_metric_df[slices_gradient_metric_df$reference == pairs[i, "cell1"] &
                                                                    slices_gradient_metric_df$target == pairs[i, "cell2"], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_gradient_metric_df_subset, gradient_metric_df_subset, by = "spe", suffix = c("_2D", "_3D"))
    
    for (radii_colname in radii_colnames) {
      slices_gradient_metric_df_subset[ , radii_colname] <- 
        (joint_df[ , paste(radii_colname, "_2D", sep = "")] - joint_df[ , paste(radii_colname, "_3D", sep = "")]) / joint_df[ , paste(radii_colname, "_3D", sep = "")]
    }
    
    
    # Combine spes_table and mixed_AMD_df
    plot_df <- cbind(spes_table, slices_gradient_metric_df_subset)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , radii_colnames)
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    plot_df$key <- paste(plot_df$spe, plot_df$slice, sep = "_")
    
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df$variable <- unfactor(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    
    plot_df <- plot_df[plot_df$variable >= min_radius & plot_df$variable <= max_radius, ]
  
    fig_slice <- ggplot(plot_df, aes(variable, value, group = key, col = slice)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
      
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = key, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
  
    fig_bg_prop_A <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_A)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))

    fig_bg_prop_B <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_B)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))

    # fig_shape <- ggplot(plot_df, aes(variable, value, group = key, col = shape)) +
    #   geom_line() +
    #   theme_bw()
    # 
    # radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    # plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    # plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    # 
    # fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = variation_E)) +
    #   geom_line() +
    #   labs(x = "radius", y = metric) +
    #   theme_bw()
    # 
    # fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = volume_E)) +
    #   geom_line() +
    #   labs(x = "radius", y = metric) +
    #   theme_bw() +
    #   scale_color_continuous(n.breaks = 4)
    # 
    # fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = spe, col = width_N)) +
    #   geom_line() +
    #   labs(x = "radius", y = metric) +
    #   theme_bw()
    
    all_plots_list[[pairs[i, "pair"]]] <- list(slice = fig_slice + theme(legend.position = "none"),
                                               arrangement = fig_arrangement + theme(legend.position = "none"),
                                               bg_prop_A = fig_bg_prop_A + theme(legend.position = "none"),
                                               bg_prop_B = fig_bg_prop_B + theme(legend.position = "none"))
                                               # shape = fig_shape + theme(legend.position = "none"))
                                               # variation_E = fig_variation_E + theme(legend.position = "none"),
                                               # volume_E = fig_volume_E + theme(legend.position = "none"),
                                               # width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  # legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  # legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  # legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  # legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       legend_arrangement,
                       legend_bg_prop_a,
                       legend_bg_prop_B,
                       # legend_shape,
                       # legend_variation_E,
                       # legend_volume_E,
                       # legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference cell type
  plots_pair_list <- list()
  
  for (i in seq(nrow(pairs))) {
    pair <- pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$slice,
                       all_plots_list[[pair]]$arrangement,
                       all_plots_list[[pair]]$bg_prop_A,
                       all_plots_list[[pair]]$bg_prop_B,
                       # all_plots_list[[pair]]$shape,
                       # all_plots_list[[pair]]$variation_E,
                       # all_plots_list[[pair]]$volume_E,
                       # all_plots_list[[pair]]$width_N, 
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

plot_gradient_metrics_type2_boxplot <- function(spes_table, gradient_metric_df, slices_gradient_metric_df, metric, arrangement_colname, min_radius, max_radius) {
  
  # Constants
  pairs <- data.frame(cell1 = c("A", "A", "B", "B"),
                      cell2 = c("A", "B", "A", "B"))
  pairs$pair <- paste(pairs$cell1, pairs$cell2, sep = "/")
  
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(nrow(pairs))) {
    
    # Subset gradient_metric_df for current reference cell
    gradient_metric_df_subset <- gradient_metric_df[gradient_metric_df$reference == pairs[i, "cell1"] &
                                                      gradient_metric_df$target == pairs[i, "cell2"], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_gradient_metric_df_subset <- slices_gradient_metric_df[slices_gradient_metric_df$reference == pairs[i, "cell1"] &
                                                                    slices_gradient_metric_df$target == pairs[i, "cell2"], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_gradient_metric_df_subset, gradient_metric_df_subset, by = "spe", suffix = c("_2D", "_3D"))
    
    for (radii_colname in radii_colnames) {
      slices_gradient_metric_df_subset[ , radii_colname] <- 
        (joint_df[ , paste(radii_colname, "_2D", sep = "")] - joint_df[ , paste(radii_colname, "_3D", sep = "")]) / joint_df[ , paste(radii_colname, "_3D", sep = "")]
    }
    
    
    # Combine spes_table and mixed_AMD_df
    plot_df <- cbind(spes_table, slices_gradient_metric_df_subset)
    
    plot_df$slice <- as.character(plot_df$slice)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , radii_colnames)
    
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df$variable <- unfactor(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    plot_df <- plot_df[plot_df$variable >= min_radius & plot_df$variable <= max_radius, ]
    plot_df$variable <- factor(as.character(plot_df$variable), radii)
    
    fig_boxplot <- ggplot(plot_df, aes(variable, value, col = slice)) +
      geom_boxplot() +
      theme_bw() +
      facet_wrap(~slice, ncol = 5) +
      scale_color_manual(values = viridis::viridis(5)) +
      labs(x = "radius", y = metric)
    
    all_plots_list[[pairs[i, "pair"]]] <- list(fig_boxplot = fig_boxplot)
  }
  # Combine the plots together by reference cell type
  plots_pair_list <- list()
  
  for (i in seq(nrow(pairs))) {
    pair <- pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$fig_boxplot,
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
                             nrow = 4, ncol = 1)
  
  methods::show(combined_plot)
  
  return(combined_plot)
}
### 1.3.1. Function to get plot for proportion SAC ----------------------------------------
plot_proportion_SAC <- function(spes_table, SAC_df, slices_SAC_df, arrangement_colname) {
  
  # Get possible reference and target cell combinations
  prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
  prop_cell_types$pair <- paste(prop_cell_types$ref, prop_cell_types$tar, sep = "/")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    
    # Subset SAC_df for chosen pair
    SAC_df_subset <- SAC_df[SAC_df$reference == prop_cell_types$ref[i], ]
    
    # Subset slices_SAC_df for chosen pair
    slices_SAC_df_subset <- slices_SAC_df[slices_SAC_df$reference == prop_cell_types$ref[i], ]
    
    # Get difference between SAC values in 3D and 2D slices.
    joint_df <- full_join(slices_SAC_df_subset, SAC_df_subset, "spe", suffix = c("_2D", "_3D"))
    
    slices_SAC_df_subset$SAC <- (joint_df$proportion_2D - joint_df$proportion_3D) / joint_df$proportion_3D
    
    # Combine spes_table and SAC_df
    plot_df <- cbind(spes_table, slices_SAC_df_subset)
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    
    fig_slice <- ggplot(plot_df, aes(!!sym(arrangement_colname), proportion, col = slice)) +
      geom_point() +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
    fig_bg_prop_A <- ggplot(plot_df, aes(!!sym(arrangement_colname), proportion, col = bg_prop_A)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_bg_prop_B <- ggplot(plot_df, aes(!!sym(arrangement_colname), proportion, col = bg_prop_B)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    # fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), proportion, col = shape)) +
    #   geom_point() +
    #   theme_bw()
    # 
    # radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    # plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    # plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    # 
    # fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), proportion, col = variation_E)) +
    #   geom_point() +
    #   theme_bw()
    
    # fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), proportion, col = volume_E)) +
    #   geom_point() +
    #   theme_bw() +
    #   scale_color_continuous(n.breaks = 4)
    # 
    # fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), proportion, col = width_N)) +
    #   geom_point() +
    #   theme_bw()
    
    all_plots_list[[prop_cell_types[i, "pair"]]] <- list(slice = fig_slice + theme(legend.position="none"),
                                                         bg_prop_A = fig_bg_prop_A + theme(legend.position="none"), 
                                                         bg_prop_B = fig_bg_prop_B + theme(legend.position="none"))
                                                         # shape = fig_shape + theme(legend.position="none"))
                                                         # variation_E = fig_variation_E + theme(legend.position="none"),
                                                         # volume_E = fig_volume_E + theme(legend.position="none"),
                                                         # width_N = fig_width_N + theme(legend.position="none"))
  }
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  # legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  # legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  # legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  # legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       legend_bg_prop_a, 
                       legend_bg_prop_B,
                       # legend_shape,
                       # legend_variation_E,
                       # legend_volume_E,
                       # legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by pairs
  
  plots_pair_list <- list()
  
  for (i in seq(nrow(prop_cell_types))) {
    pair <- prop_cell_types[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$slice,
                       all_plots_list[[pair]]$bg_prop_A,
                       all_plots_list[[pair]]$bg_prop_B,
                       # all_plots_list[[pair]]$shape, 
                       # all_plots_list[[pair]]$variation_E,
                       # all_plots_list[[pair]]$volume_E,
                       # all_plots_list[[pair]]$width_N, 
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
plot_entropy_SAC <- function(spes_table, SAC_df, slices_SAC_df, arrangement_colname) {
  
  # Get possible cell type of interest combinations
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    
    # Subset SAC_df for chosen pair
    SAC_df_subset <- SAC_df[SAC_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Subset slices_SAC_df for chosen pair
    slices_SAC_df_subset <- slices_SAC_df[slices_SAC_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Get difference between SAC values in 3D and 2D slices.
    joint_df <- full_join(slices_SAC_df_subset, SAC_df_subset, "spe", suffix = c("_2D", "_3D"))
    
    slices_SAC_df_subset$SAC <- (joint_df$entropy_2D - joint_df$entropy_3D) / joint_df$entropy_3D
    
    # Combine spes_table and SAC_df
    plot_df <- cbind(spes_table, slices_SAC_df_subset)
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    
    fig_slice <- ggplot(plot_df, aes(!!sym(arrangement_colname), entropy, col = slice)) +
      geom_point() +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
    fig_bg_prop_A <- ggplot(plot_df, aes(!!sym(arrangement_colname), entropy, col = bg_prop_A)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_bg_prop_B <- ggplot(plot_df, aes(!!sym(arrangement_colname), entropy, col = bg_prop_B)) +
      geom_point() +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    # fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), entropy, col = shape)) +
    #   geom_point() +
    #   theme_bw()
    # 
    # radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    # plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    # plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    # 
    # fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), entropy, col = variation_E)) +
    #   geom_point() +
    #   theme_bw()
    # 
    # fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), entropy, col = volume_E)) +
    #   geom_point() +
    #   theme_bw() +
    #   scale_color_continuous(n.breaks = 4)
    # 
    # fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), entropy, col = width_N)) +
    #   geom_point() +
    #   theme_bw()
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(slice = fig_slice + theme(legend.position="none"),
                                                               bg_prop_A = fig_bg_prop_A + theme(legend.position="none"), 
                                                               bg_prop_B = fig_bg_prop_B + theme(legend.position="none"))
                                                               # shape = fig_shape + theme(legend.position="none"))
                                                               # variation_E = fig_variation_E + theme(legend.position="none"),
                                                               # volume_E = fig_volume_E + theme(legend.position="none"),
                                                               # width_N = fig_width_N + theme(legend.position="none"))
  }
  
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  # legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  # legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  # legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  # legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       legend_bg_prop_a, 
                       legend_bg_prop_B,
                       # legend_shape,
                       # legend_variation_E,
                       # legend_volume_E,
                       # legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by cell types of interest
  
  plots_cell_types_list <- list()
  
  for (i in seq(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$slice,
                       all_plots_list[[cell_types]]$bg_prop_A,
                       all_plots_list[[cell_types]]$bg_prop_B,
                       # all_plots_list[[cell_types]]$shape, 
                       # all_plots_list[[cell_types]]$variation_E,
                       # all_plots_list[[cell_types]]$volume_E,
                       # all_plots_list[[cell_types]]$width_N, 
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
plot_proportion_prevalence <- function(spes_table, prevalence_df, slices_prevalence_df, arrangement_colname) {
  
  # Constants
  prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
  prop_cell_types$pair <- paste(prop_cell_types$ref, prop_cell_types$tar, sep = "/")
  
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    
    # Subset gradient_metric_df for current reference cell
    prevalence_df_subset <- prevalence_df[prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_prevalence_df_subset <- slices_prevalence_df[slices_prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_prevalence_df_subset, prevalence_df_subset, by = "spe", suffix = c("_2D", "_3D"))

    for (threshold_colname in threshold_colnames) {
      slices_prevalence_df_subset[ , threshold_colname] <- 
        (joint_df[ , paste(threshold_colname, "_2D", sep = "")] - joint_df[ , paste(threshold_colname, "_3D", sep = "")]) / joint_df[ , paste(threshold_colname, "_3D", sep = "")]
    }
    
    # Combine spes_table and prevalence_df
    plot_df <- cbind(spes_table, slices_prevalence_df_subset)

    # Melt
    plot_df <- reshape2::melt(plot_df, , threshold_colnames)
    
    # Extract threshold value from threshold strings (t0.01 -> 0.01...)
    plot_df$variable <- as.character(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    plot_df$key <- paste(plot_df$spe, plot_df$slice, sep = "_")
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = key, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    fig_slice <- ggplot(plot_df, aes(variable, value, group = key, col = slice)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
    fig_bg_prop_A <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_A)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))

    fig_bg_prop_B <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_B)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))

    # fig_shape <- ggplot(plot_df, aes(variable, value, group = key, col = shape)) +
    #   labs(x = "threshold", y = "prevalence") +
    #   geom_line() +
    #   theme_bw()
    # 
    # radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    # plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    # plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    # 
    # fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = variation_E)) +
    #   geom_line() +
    #   labs(x = "threshold", y = "prevalence") +
    #   theme_bw()
    # 
    # fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = volume_E)) +
    #   geom_line() +
    #   labs(x = "threshold", y = "prevalence") +
    #   theme_bw() +
    #   scale_color_continuous(n.breaks = 4)
    # 
    # fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = spe, col = width_N)) +
    #   geom_line() +
    #   labs(x = "threshold", y = "prevalence") +
    #   theme_bw()
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(slice = fig_slice + theme(legend.position = "none"),
                                                      arrangement = fig_arrangement + theme(legend.position = "none"), 
                                                      bg_prop_A = fig_bg_prop_A + theme(legend.position = "none"),
                                                      bg_prop_B = fig_bg_prop_B + theme(legend.position = "none"))
                                                      # shape = fig_shape + theme(legend.position = "none"))
                                                      # variation_E = fig_variation_E + theme(legend.position = "none"),
                                                      # volume_E = fig_volume_E + theme(legend.position = "none"),
                                                      # width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  # legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  # legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  # legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  # legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       legend_arrangement,
                       legend_bg_prop_a,
                       legend_bg_prop_B,
                       # legend_shape,
                       # legend_variation_E,
                       # legend_volume_E,
                       # legend_width_N,
                       nrow = 1)
  
  
  # Combine the plots together by reference (and target) cell pairs
  plots_pair_list <- list()
  
  for (i in seq(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$slice,
                       all_plots_list[[pair]]$arrangement,
                       all_plots_list[[pair]]$bg_prop_A,
                       all_plots_list[[pair]]$bg_prop_B,
                       # all_plots_list[[pair]]$shape,
                       # all_plots_list[[pair]]$variation_E,
                       # all_plots_list[[pair]]$volume_E,
                       # all_plots_list[[pair]]$width_N, 
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

plot_proportion_prevalence_boxplot <- function(spes_table, prevalence_df, slices_prevalence_df, arrangement_colname) {
  
  # Constants
  prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
  prop_cell_types$pair <- paste(prop_cell_types$ref, prop_cell_types$tar, sep = "/")
  
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")

  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    
    # Subset gradient_metric_df for current reference cell
    prevalence_df_subset <- prevalence_df[prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_prevalence_df_subset <- slices_prevalence_df[slices_prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_prevalence_df_subset, prevalence_df_subset, by = "spe", suffix = c("_2D", "_3D"))
    
    for (threshold_colname in threshold_colnames) {
      slices_prevalence_df_subset[ , threshold_colname] <- 
        (joint_df[ , paste(threshold_colname, "_2D", sep = "")] - joint_df[ , paste(threshold_colname, "_3D", sep = "")]) / joint_df[ , paste(threshold_colname, "_3D", sep = "")]
    }
    
    # Combine spes_table and prevalence_df
    plot_df <- cbind(spes_table, slices_prevalence_df_subset)
    
    plot_df$slice <- as.character(plot_df$slice)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , threshold_colnames)
    
    # Extract threshold value from threshold strings (t0.01 -> 0.01...)
    plot_df$variable <- as.character(plot_df$variable)
    plot_df$variable <- substr(plot_df$variable, 2, nchar(plot_df$variable))
    plot_df$variable <- factor(plot_df$variable, as.character(thresholds))
    
    fig_boxplot <- ggplot(plot_df, aes(variable, value, col = slice)) +
      geom_boxplot() +
      theme_bw() +
      facet_wrap(~slice, ncol = 5) +
      scale_color_manual(values = viridis::viridis(5)) +
      labs(x = "threshold", y = "prevalence") +
      scale_x_discrete(breaks = c(0.01, 1), labels = c("0.01", "1"))
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(boxplot = fig_boxplot)
  }

    # Combine the plots together by reference (and target) cell pairs
  plots_pair_list <- list()
  
  for (i in seq(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$boxplot,
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
                             nrow = 2, ncol = 1,
                             rel_heights = c(1, 1, 0.5))
  
  methods::show(combined_plot)
  
  return(combined_plot)
}

### 1.3.4. Function to get plot for entropy prevalence ----------------------------------
plot_entropy_prevalence <- function(spes_table, prevalence_df, slices_prevalence_df, arrangement_colname) {
  
  # Constants
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    
    # Subset gradient_metric_df for current reference cell
    prevalence_df_subset <- prevalence_df[prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_prevalence_df_subset <- slices_prevalence_df[slices_prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_prevalence_df_subset, prevalence_df_subset, by = "spe", suffix = c("_2D", "_3D"))
    
    for (threshold_colname in threshold_colnames) {
      slices_prevalence_df_subset[ , threshold_colname] <- 
        (joint_df[ , paste(threshold_colname, "_2D", sep = "")] - joint_df[ , paste(threshold_colname, "_3D", sep = "")]) / joint_df[ , paste(threshold_colname, "_3D", sep = "")]
    }
    
    # Combine spes_table and prevalence_df
    plot_df <- cbind(spes_table, slices_prevalence_df_subset)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , threshold_colnames)
    
    # Extract threshold value from threshold strings (t0.01 -> 0.01...)
    plot_df$variable <- as.character(plot_df$variable)
    plot_df$variable <- substr(plot_df$variable, 2, nchar(plot_df$variable))
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    plot_df$key <- paste(plot_df$spe, plot_df$slice, sep = "_")
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = key, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    fig_slice <- ggplot(plot_df, aes(variable, value, group = key, col = slice)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))

    fig_bg_prop_A <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_A)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))

    fig_bg_prop_B <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_B)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_continuous(breaks = c(0.0, 0.05, 0.1))

    # fig_shape <- ggplot(plot_df, aes(variable, value, group = key, col = shape)) +
    #   labs(x = "threshold", y = "prevalence") +
    #   geom_line() +
    #   theme_bw()
    #
    # radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    # plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    # plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    # 
    # fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = variation_E)) +
    #   geom_line() +
    #   labs(x = "threshold", y = "prevalence") +
    #   theme_bw()
    # 
    # fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = spe, col = volume_E)) +
    #   geom_line() +
    #   labs(x = "threshold", y = "prevalence") +
    #   theme_bw() +
    #   scale_color_continuous(n.breaks = 4)
    # 
    # fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = spe, col = width_N)) +
    #   geom_line() +
    #   labs(x = "threshold", y = "prevalence") +
    #   theme_bw()
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(slice = fig_slice + theme(legend.position = "none"),
                                                               arrangement = fig_arrangement + theme(legend.position = "none"), 
                                                               bg_prop_A = fig_bg_prop_A + theme(legend.position = "none"),
                                                               bg_prop_B = fig_bg_prop_B + theme(legend.position = "none"))
                                                               # shape = fig_shape + theme(legend.position = "none"))
                                                               # variation_E = fig_variation_E + theme(legend.position = "none"),
                                                               # volume_E = fig_volume_E + theme(legend.position = "none"),
                                                               # width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  # legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  # legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  # legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  # legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       legend_arrangement,
                       legend_bg_prop_a,
                       legend_bg_prop_B,
                       # legend_shape,
                       # legend_variation_E,
                       # legend_volume_E,
                       # legend_width_N,
                       nrow = 1)
  
  
  # Combine the plots together by reference (and target) cell pairs
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$slice,
                       all_plots_list[[cell_types]]$arrangement,
                       all_plots_list[[cell_types]]$bg_prop_A,
                       all_plots_list[[cell_types]]$bg_prop_B,
                       # all_plots_list[[cell_types]]$shape,
                       # all_plots_list[[cell_types]]$variation_E,
                       # all_plots_list[[cell_types]]$volume_E,
                       # all_plots_list[[cell_types]]$width_N, 
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

plot_entropy_prevalence_boxplot <- function(spes_table, prevalence_df, slices_prevalence_df, arrangement_colname) {
  
  # Constants
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    
    # Subset gradient_metric_df for current reference cell
    prevalence_df_subset <- prevalence_df[prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_prevalence_df_subset <- slices_prevalence_df[slices_prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_prevalence_df_subset, prevalence_df_subset, by = "spe", suffix = c("_2D", "_3D"))
    
    for (threshold_colname in threshold_colnames) {
      slices_prevalence_df_subset[ , threshold_colname] <- 
        (joint_df[ , paste(threshold_colname, "_2D", sep = "")] - joint_df[ , paste(threshold_colname, "_3D", sep = "")]) / joint_df[ , paste(threshold_colname, "_3D", sep = "")]
    }
    
    # Combine spes_table and prevalence_df
    plot_df <- cbind(spes_table, slices_prevalence_df_subset)
    
    plot_df$slice <- as.character(plot_df$slice)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , threshold_colnames)
    
    # Extract threshold value from threshold strings (t0.01 -> 0.01...)
    plot_df$variable <- as.character(plot_df$variable)
    plot_df$variable <- substr(plot_df$variable, 2, nchar(plot_df$variable))
    plot_df$variable <- factor(plot_df$variable, as.character(thresholds))
    
    fig_boxplot <- ggplot(plot_df, aes(variable, value, col = slice)) +
      geom_boxplot() +
      theme_bw() +
      facet_wrap(~slice, ncol = 5) +
      scale_color_manual(values = viridis::viridis(5)) +
      labs(x = "threshold", y = "prevalence") +
      scale_x_discrete(breaks = c(0.01, 1), labels = c("0.01", "1"))
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(boxplot = fig_boxplot)
  }

    # Combine the plots together by reference (and target) cell pairs
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$boxplot,
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
                             nrow = 2, ncol = 1,
                             rel_heights = c(1, 1, 0.5))
  
  methods::show(combined_plot)
  
  return(combined_plot)
}


### 1.3.5. Function to get plot for proportion prevalence AUC -----------------
plot_proportion_prevalence_AUC <- function(spes_table, prevalence_df, slices_prevalence_df, arrangement_colname) {
  
  # Constants
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
  prop_cell_types$pair <- paste(prop_cell_types$ref, prop_cell_types$tar, sep = "/")
  
  # Get AUC for each prevalence gradient
  prevalence_df$AUC <- apply(prevalence_df[ , threshold_colnames], 1, sum) * 0.01
  prevalence_df <- prevalence_df[ , c("spe", "reference", "target", "AUC")]
  
  # Get AUC for each prevalence gradient (slices)
  slices_prevalence_df$AUC <- apply(slices_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
  slices_prevalence_df <- slices_prevalence_df[ , c("spe", "slice", "reference", "target", "AUC")]
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    
    # Subset prevalence_df for chosen pair
    prevalence_df_subset <- prevalence_df[prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Subset slices_prevalence_df for chosen pair
    slices_prevalence_df_subset <- slices_prevalence_df[slices_prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Get difference between prevalence values in 3D and 2D slices.
    joint_df <- full_join(slices_prevalence_df_subset, prevalence_df_subset, "spe", suffix = c("_2D", "_3D"))
    
    slices_prevalence_df_subset$AUC <- (joint_df$AUC_2D - joint_df$AUC_3D) / joint_df$AUC_3D
    
    # Combine spes_table and SAC_df
    plot_df <- cbind(spes_table, slices_prevalence_df_subset)

    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    
    fig_slice <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = slice)) +
      geom_point() +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
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
    
    # fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = shape)) +
    #   geom_point() +
    #   ylab("AUC") +
    #   theme_bw()
    #
    # radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    # plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    # plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    # 
    # fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AUC, col = variation_E)) +
    #   geom_point() +
    #   ylab("AUC") +
    #   theme_bw()
    # 
    # fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AUC, col = volume_E)) +
    #   geom_point() +
    #   ylab("AUC") +
    #   theme_bw() +
    #   scale_color_continuous(n.breaks = 4)
    # 
    # fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), AUC, col = width_N)) +
    #   geom_point() +
    #   ylab("AUC") +
    #   theme_bw()
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(slice = fig_slice + theme(legend.position="none"), 
                                                      bg_prop_A = fig_bg_prop_A + theme(legend.position="none"), 
                                                      bg_prop_B = fig_bg_prop_B + theme(legend.position="none"))
                                                      # shape = fig_shape + theme(legend.position="none"))
                                                      # variation_E = fig_variation_E + theme(legend.position="none"),
                                                      # volume_E = fig_volume_E + theme(legend.position="none"),
                                                      # width_N = fig_width_N + theme(legend.position="none"))
  }
  
  # Get legends
  legend_slice <-  get_legend(fig_slice + theme(legend.direction = "horizontal"))
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  # legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  # legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  # legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  # legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       legend_bg_prop_a, 
                       legend_bg_prop_B,
                       # legend_shape,
                       # legend_variation_E,
                       # legend_volume_E,
                       # legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference target pairs
  plots_pair_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$slice,
                       all_plots_list[[pair]]$bg_prop_A,
                       all_plots_list[[pair]]$bg_prop_B,
                       # all_plots_list[[pair]]$shape, 
                       # all_plots_list[[pair]]$variation_E,
                       # all_plots_list[[pair]]$volume_E,
                       # all_plots_list[[pair]]$width_N, 
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
plot_entropy_prevalence_AUC <- function(spes_table, prevalence_df, slices_prevalence_df, arrangement_colname) {
  
  # Constants
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  
  # Get AUC for each prevalence gradient
  prevalence_df$AUC <- apply(prevalence_df[ , threshold_colnames], 1, sum) * 0.01
  prevalence_df <- prevalence_df[ , c("spe", "cell_types", "AUC")]
  
  # Get AUC for each prevalence gradient (slices)
  slices_prevalence_df$AUC <- apply(slices_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
  slices_prevalence_df <- slices_prevalence_df[ , c("spe", "slice", "cell_types", "AUC")]
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    
    # Subset prevalence_df for chosen pair
    prevalence_df_subset <- prevalence_df[prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Subset slices_prevalence_df for chosen pair
    slices_prevalence_df_subset <- slices_prevalence_df[slices_prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Get difference between prevalence values in 3D and 2D slices.
    joint_df <- full_join(slices_prevalence_df_subset, prevalence_df_subset, "spe", suffix = c("_2D", "_3D"))
    
    slices_prevalence_df_subset$AUC <- (joint_df$AUC_2D - joint_df$AUC_3D) / joint_df$AUC_3D
    
    # Combine spes_table and SAC_df
    plot_df <- cbind(spes_table, slices_prevalence_df_subset)
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    
    fig_slice <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = slice)) +
      geom_point() +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
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
    
    # fig_shape <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = shape)) +
    #   geom_point() +
    #   ylab("AUC") +
    #   theme_bw()
    # 
    # radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    # plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    # plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
    # 
    # fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AUC, col = variation_E)) +
    #   geom_point() +
    #   ylab("AUC") +
    #   theme_bw()
    # 
    # fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(!!sym(arrangement_colname), AUC, col = volume_E)) +
    #   geom_point() +
    #   ylab("AUC") +
    #   theme_bw() +
    #   scale_color_continuous(n.breaks = 4)
    # 
    # fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(!!sym(arrangement_colname), AUC, col = width_N)) +
    #   geom_point() +
    #   ylab("AUC") +
    #   theme_bw()
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(slice = fig_slice + theme(legend.position="none"),
                                                               bg_prop_A = fig_bg_prop_A + theme(legend.position="none"), 
                                                               bg_prop_B = fig_bg_prop_B + theme(legend.position="none"))
                                                               # shape = fig_shape + theme(legend.position="none"))
                                                               # variation_E = fig_variation_E + theme(legend.position="none"),
                                                               # volume_E = fig_volume_E + theme(legend.position="none"),
                                                               # width_N = fig_width_N + theme(legend.position="none"))
  }
  
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  # legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  # legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  # legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  # legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice, 
                       legend_bg_prop_a, 
                       legend_bg_prop_B,
                       # legend_shape,
                       # legend_variation_E,
                       # legend_volume_E,
                       # legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference target pairs
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$slice,
                       all_plots_list[[cell_types]]$bg_prop_A,
                       all_plots_list[[cell_types]]$bg_prop_B,
                       # all_plots_list[[cell_types]]$shape, 
                       # all_plots_list[[cell_types]]$variation_E,
                       # all_plots_list[[cell_types]]$volume_E,
                       # all_plots_list[[cell_types]]$width_N, 
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

# Read mixed_slices_AMD_df
setwd("~/Objects/unsupervised/mixed_spes/analysis_2D")
mixed_slices_AMD_df <- read.table("mixed_slices_AMD_df.csv")

mixed_AMD_plot <- plot_AMD_metric(mixed_spes_table, mixed_AMD_df, mixed_slices_AMD_df, "cluster_prop_B")

mixed_AMD_density_plot <- plot_AMD_density_metric(mixed_spes_table, mixed_AMD_df, mixed_slices_AMD_df, "cluster_prop_B")

setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(mixed_AMD_plot, "mixed_AMD_plot_slicing.RDS")

### 2.3. mixed spes MS, NMS, ACINP, AE -----------------------------------

# Read mixed MS, NMS, ACINP, AE dfs
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_MS_df <- read.table("mixed_MS_df.csv")
mixed_NMS_df <- read.table("mixed_NMS_df.csv")
mixed_ACINP_df <- read.table("mixed_ACINP_df.csv")
mixed_AE_df <- read.table("mixed_AE_df.csv")

# Read mixed MS, NMS, ACINP, AE dfs (slices)
setwd("~/Objects/unsupervised/mixed_spes/analysis_2D")
mixed_slices_MS_df <- read.table("mixed_slices_MS_df.csv")
mixed_slices_NMS_df <- read.table("mixed_slices_NMS_df.csv")
mixed_slices_ACINP_df <- read.table("mixed_slices_ACINP_df.csv")
mixed_slices_AE_df <- read.table("mixed_slices_AE_df.csv")

# Lines
mixed_MS_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_MS_df, mixed_slices_MS_df, "MS", "cluster_prop_B")
mixed_NMS_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_NMS_df, mixed_slices_NMS_df, "NMS", "cluster_prop_B")
mixed_ACINP_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_ACINP_df, mixed_slices_ACINP_df, "ACINP", "cluster_prop_B")
mixed_AE_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_AE_df, mixed_slices_AE_df, "AE", "cluster_prop_B")

# Boxplots
mixed_MS_box_plot <- plot_gradient_metrics_type1_boxplot(mixed_spes_table, mixed_MS_df, mixed_slices_MS_df, "MS", "cluster_prop_B")
mixed_NMS_box_plot <- plot_gradient_metrics_type1_boxplot(mixed_spes_table, mixed_NMS_df, mixed_slices_NMS_df, "NMS", "cluster_prop_B")
mixed_ACINP_box_plot <- plot_gradient_metrics_type1_boxplot(mixed_spes_table, mixed_ACINP_df, mixed_slices_ACINP_df, "ACINP", "cluster_prop_B")
mixed_AE_box_plot <- plot_gradient_metrics_type1_boxplot(mixed_spes_table, mixed_AE_df, mixed_slices_AE_df, "AE", "cluster_prop_B")


setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(mixed_MS_plot, "mixed_MS_plot_slicing.RDS")
saveRDS(mixed_NMS_plot, "mixed_NMS_plot_slicing.RDS")
saveRDS(mixed_ACINP_plot, "mixed_ACINP_plot_slicing.RDS")
saveRDS(mixed_AE_plot, "mixed_AE_plot_slicing.RDS")

saveRDS(mixed_MS_box_plot, "mixed_MS_box_plot_slicing.RDS")
saveRDS(mixed_NMS_box_plot, "mixed_NMS_box_plot_slicing.RDS")
saveRDS(mixed_ACINP_box_plot, "mixed_ACINP_box_plot_slicing.RDS")
saveRDS(mixed_AE_box_plot, "mixed_AE_box_plot_slicing.RDS")

### 2.4. mixed spes ACIN, CKR ------------------------------------------------

# Read mixed ACIN, CKR
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_ACIN_df <- read.table("mixed_ACIN_df.csv")
mixed_CKR_df <- read.table("mixed_CKR_df.csv")

# Read mixed ACIN, CKR (slices)
setwd("~/Objects/unsupervised/mixed_spes/analysis_2D")
mixed_slices_ACIN_df <- read.table("mixed_slices_ACIN_df.csv")
mixed_slices_CKR_df <- read.table("mixed_slices_CKR_df.csv")

# Lines
mixed_ACIN_plot <- plot_gradient_metrics_type2(mixed_spes_table, mixed_ACIN_df, mixed_slices_ACIN_df, "ACIN", "cluster_prop_B", 20, 100)
mixed_CKR_plot <- plot_gradient_metrics_type2(mixed_spes_table, mixed_CKR_df, mixed_slices_CKR_df, "CKR", "cluster_prop_B", 20, 100)

# Box plots
mixed_ACIN_box_plot <- plot_gradient_metrics_type2_boxplot(mixed_spes_table, mixed_ACIN_df, mixed_slices_ACIN_df, "ACIN", "cluster_prop_B", 20, 100)
mixed_CKR_box_plot <- plot_gradient_metrics_type2_boxplot(mixed_spes_table, mixed_CKR_df, mixed_slices_CKR_df, "CKR", "cluster_prop_B", 20, 100)

setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(mixed_ACIN_plot, "mixed_ACIN_plot_slicing.RDS")
saveRDS(mixed_CKR_plot, "mixed_CKR_plot_slicing.RDS")

saveRDS(mixed_ACIN_box_plot, "mixed_ACIN_box_plot_slicing.RDS")
saveRDS(mixed_CKR_box_plot, "mixed_CKR_box_plot_slicing.RDS")

### 2.5. mixed spes SAC ------------------------------------------------------

# Read mixed_SAC_df
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_prop_SAC_df <- read.table("mixed_prop_SAC_df.csv")
mixed_entropy_SAC_df <- read.table("mixed_entropy_SAC_df.csv")

# Read mixed_SAC_df (slices)
setwd("~/Objects/unsupervised/mixed_spes/analysis_2D")
mixed_slices_prop_SAC_df <- read.table("mixed_slices_prop_SAC_df.csv")
mixed_slices_entropy_SAC_df <- read.table("mixed_slices_entropy_SAC_df.csv")

mixed_prop_SAC_plot <- plot_proportion_SAC(mixed_spes_table, mixed_prop_SAC_df, mixed_slices_prop_SAC_df, "cluster_prop_B")
mixed_entropy_SAC_plot <- plot_entropy_SAC(mixed_spes_table, mixed_entropy_SAC_df, mixed_slices_entropy_SAC_df, "cluster_prop_B")

setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(mixed_prop_SAC_plot, "mixed_prop_SAC_plot_slicing.RDS")
saveRDS(mixed_entropy_SAC_plot, "mixed_entropy_SAC_plot_slicing.RDS")



### 2.6. mixed spes prevalence ------------------------------------------------

# Read mixed prevalence dfs
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_prop_prevalence_df <- read.table("mixed_prop_prevalence_df.csv")
mixed_entropy_prevalence_df <- read.table("mixed_entropy_prevalence_df.csv")

# Read mixed prevalence dfs (slices)
setwd("~/Objects/unsupervised/mixed_spes/analysis_2D")
mixed_slices_prop_prevalence_df <- read.table("mixed_slices_prop_prevalence_df.csv")
mixed_slices_entropy_prevalence_df <- read.table("mixed_slices_entropy_prevalence_df.csv")

# Lines
mixed_prop_prevalence_plot <- plot_proportion_prevalence(mixed_spes_table, mixed_prop_prevalence_df, mixed_slices_prop_prevalence_df, "cluster_prop_B")
mixed_entropy_prevalence_plot <- plot_entropy_prevalence(mixed_spes_table, mixed_entropy_prevalence_df, mixed_slices_entropy_prevalence_df, "cluster_prop_B")

# Boxplots
mixed_prop_prevalence_box_plot <- plot_proportion_prevalence_boxplot(mixed_spes_table, mixed_prop_prevalence_df, mixed_slices_prop_prevalence_df, "cluster_prop_B")
mixed_entropy_prevalence_box_plot <- plot_entropy_prevalence_boxplot(mixed_spes_table, mixed_entropy_prevalence_df, mixed_slices_entropy_prevalence_df, "cluster_prop_B")

# AUC
mixed_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(mixed_spes_table, mixed_prop_prevalence_df, mixed_slices_prop_prevalence_df, "cluster_prop_B")
mixed_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(mixed_spes_table, mixed_entropy_prevalence_df, mixed_slices_entropy_prevalence_df, "cluster_prop_B")

setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(mixed_prop_prevalence_plot, "mixed_prop_prevalence_plot_slicing.RDS")
saveRDS(mixed_entropy_prevalence_plot, "mixed_entropy_prevalence_plot_slicing.RDS")

saveRDS(mixed_prop_prevalence_box_plot, "mixed_prop_prevalence_box_plot_slicing.RDS")
saveRDS(mixed_entropy_prevalence_box_plot, "mixed_entropy_prevalence_box_plot_slicing.RDS")

saveRDS(mixed_prop_prevalence_AUC_plot, "mixed_prop_prevalence_AUC_plot_slicing.RDS")
saveRDS(mixed_entropy_prevalence_AUC_plot, "mixed_entropy_prevalence_AUC_plot_slicing.RDS")


### 3.2. ringed spes AMD ------------------------------------------------------

# Read ringed_AMD_df
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_AMD_df <- read.table("ringed_AMD_df.csv")

# Read ringed_slices_AMD_df
setwd("~/Objects/unsupervised/ringed_spes/analysis_2D")
ringed_slices_AMD_df <- read.table("ringed_slices_AMD_df.csv")

ringed_AMD_plot <- plot_AMD_metric(ringed_spes_table, ringed_AMD_df, ringed_slices_AMD_df, "width_ring_factor")

setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(ringed_AMD_plot, "ringed_AMD_plot_slicing.RDS")

### 3.3. ringed spes MS, NMS, ACINP, AE -----------------------------------

# Read ringed MS, NMS, ACINP, AE dfs
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_MS_df <- read.table("ringed_MS_df.csv")
ringed_NMS_df <- read.table("ringed_NMS_df.csv")
ringed_ACINP_df <- read.table("ringed_ACINP_df.csv")
ringed_AE_df <- read.table("ringed_AE_df.csv")

# Read ringed MS, NMS, ACINP, AE dfs (slices)
setwd("~/Objects/unsupervised/ringed_spes/analysis_2D")
ringed_slices_MS_df <- read.table("ringed_slices_MS_df.csv")
ringed_slices_NMS_df <- read.table("ringed_slices_NMS_df.csv")
ringed_slices_ACINP_df <- read.table("ringed_slices_ACINP_df.csv")
ringed_slices_AE_df <- read.table("ringed_slices_AE_df.csv")

# Lines
ringed_MS_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_MS_df, ringed_slices_MS_df, "MS", "width_ring_factor")
ringed_NMS_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_NMS_df, ringed_slices_NMS_df, "NMS", "width_ring_factor")
ringed_ACINP_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_ACINP_df, ringed_slices_ACINP_df, "ACINP", "width_ring_factor")
ringed_AE_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_AE_df, ringed_slices_AE_df, "AE", "width_ring_factor")

# Boxplots
ringed_MS_box_plot <- plot_gradient_metrics_type1_boxplot(ringed_spes_table, ringed_MS_df, ringed_slices_MS_df, "MS", "width_ring_factor")
ringed_NMS_box_plot <- plot_gradient_metrics_type1_boxplot(ringed_spes_table, ringed_NMS_df, ringed_slices_NMS_df, "NMS", "width_ring_factor")
ringed_ACINP_box_plot <- plot_gradient_metrics_type1_boxplot(ringed_spes_table, ringed_ACINP_df, ringed_slices_ACINP_df, "ACINP", "width_ring_factor")
ringed_AE_box_plot <- plot_gradient_metrics_type1_boxplot(ringed_spes_table, ringed_AE_df, ringed_slices_AE_df, "AE", "width_ring_factor")

setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(ringed_MS_plot, "ringed_MS_plot_slicing.RDS")
saveRDS(ringed_NMS_plot, "ringed_NMS_plot_slicing.RDS")
saveRDS(ringed_ACINP_plot, "ringed_ACINP_plot_slicing.RDS")
saveRDS(ringed_AE_plot, "ringed_AE_plot_slicing.RDS")

saveRDS(ringed_MS_box_plot, "ringed_MS_box_plot_slicing.RDS")
saveRDS(ringed_NMS_box_plot, "ringed_NMS_box_plot_slicing.RDS")
saveRDS(ringed_ACINP_box_plot, "ringed_ACINP_box_plot_slicing.RDS")
saveRDS(ringed_AE_box_plot, "ringed_AE_box_plot_slicing.RDS")

### 3.4. ringed spes ACIN, CKR ------------------------------------------------

# Read ringed ACIN, CKR
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_ACIN_df <- read.table("ringed_ACIN_df.csv")
ringed_CKR_df <- read.table("ringed_CKR_df.csv")

# Read ringed ACIN, CKR (slices)
setwd("~/Objects/unsupervised/ringed_spes/analysis_2D")
ringed_slices_ACIN_df <- read.table("ringed_slices_ACIN_df.csv")
ringed_slices_CKR_df <- read.table("ringed_slices_CKR_df.csv")

# Lines
ringed_ACIN_plot <- plot_gradient_metrics_type2(ringed_spes_table, ringed_ACIN_df, ringed_slices_ACIN_df, "ACIN", "width_ring_factor", 20, 100)
ringed_CKR_plot <- plot_gradient_metrics_type2(ringed_spes_table, ringed_CKR_df, ringed_slices_CKR_df, "CKR", "width_ring_factor", 20, 100)

# Box plots
ringed_ACIN_box_plot <- plot_gradient_metrics_type2_boxplot(ringed_spes_table, ringed_ACIN_df, ringed_slices_ACIN_df, "ACIN", "width_ring_factor", 20, 100)
ringed_CKR_box_plot <- plot_gradient_metrics_type2_boxplot(ringed_spes_table, ringed_CKR_df, ringed_slices_CKR_df, "CKR", "width_ring_factor", 20, 100)

setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(ringed_ACIN_plot, "ringed_ACIN_plot_slicing.RDS")
saveRDS(ringed_CKR_plot, "ringed_CKR_plot_slicing.RDS")

saveRDS(ringed_ACIN_box_plot, "ringed_ACIN_box_plot_slicing.RDS")
saveRDS(ringed_CKR_box_plot, "ringed_CKR_box_plot_slicing.RDS")

### 3.5. ringed spes SAC ------------------------------------------------------

# Read ringed_SAC_df
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_prop_SAC_df <- read.table("ringed_prop_SAC_df.csv")
ringed_entropy_SAC_df <- read.table("ringed_entropy_SAC_df.csv")

# Read ringed_SAC_df (slices)
setwd("~/Objects/unsupervised/ringed_spes/analysis_2D")
ringed_slices_prop_SAC_df <- read.table("ringed_slices_prop_SAC_df.csv")
ringed_slices_entropy_SAC_df <- read.table("ringed_slices_entropy_SAC_df.csv")

ringed_prop_SAC_plot <- plot_proportion_SAC(ringed_spes_table, ringed_prop_SAC_df, ringed_slices_prop_SAC_df, "width_ring_factor")
ringed_entropy_SAC_plot <- plot_entropy_SAC(ringed_spes_table, ringed_entropy_SAC_df, ringed_slices_entropy_SAC_df, "width_ring_factor")

setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(ringed_prop_SAC_plot, "ringed_prop_SAC_plot_slicing.RDS")
saveRDS(ringed_entropy_SAC_plot, "ringed_entropy_SAC_plot_slicing.RDS")



### 3.6. ringed spes prevalence ------------------------------------------------

# Read ringed prevalence dfs
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_prop_prevalence_df <- read.table("ringed_prop_prevalence_df.csv")
ringed_entropy_prevalence_df <- read.table("ringed_entropy_prevalence_df.csv")

# Read ringed prevalence dfs (slices)
setwd("~/Objects/unsupervised/ringed_spes/analysis_2D")
ringed_slices_prop_prevalence_df <- read.table("ringed_slices_prop_prevalence_df.csv")
ringed_slices_entropy_prevalence_df <- read.table("ringed_slices_entropy_prevalence_df.csv")

# Lines
ringed_prop_prevalence_plot <- plot_proportion_prevalence(ringed_spes_table, ringed_prop_prevalence_df, ringed_slices_prop_prevalence_df, "width_ring_factor")
ringed_entropy_prevalence_plot <- plot_entropy_prevalence(ringed_spes_table, ringed_entropy_prevalence_df, ringed_slices_entropy_prevalence_df, "width_ring_factor")

# Boxplots
ringed_prop_prevalence_box_plot <- plot_proportion_prevalence_boxplot(ringed_spes_table, ringed_prop_prevalence_df, ringed_slices_prop_prevalence_df, "width_ring_factor")
ringed_entropy_prevalence_box_plot <- plot_entropy_prevalence_boxplot(ringed_spes_table, ringed_entropy_prevalence_df, ringed_slices_entropy_prevalence_df, "width_ring_factor")

# AUC
ringed_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(ringed_spes_table, ringed_prop_prevalence_df, ringed_slices_prop_prevalence_df, "width_ring_factor")
ringed_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(ringed_spes_table, ringed_entropy_prevalence_df, ringed_slices_entropy_prevalence_df, "width_ring_factor")

setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(ringed_prop_prevalence_plot, "ringed_prop_prevalence_plot_slicing.RDS")
saveRDS(ringed_entropy_prevalence_plot, "ringed_entropy_prevalence_plot_slicing.RDS")
saveRDS(ringed_prop_prevalence_box_plot, "ringed_prop_prevalence_box_plot_slicing.RDS")
saveRDS(ringed_entropy_prevalence_box_plot, "ringed_entropy_prevalence_box_plot_slicing.RDS")
saveRDS(ringed_prop_prevalence_AUC_plot, "ringed_prop_prevalence_AUC_plot_slicing.RDS")
saveRDS(ringed_entropy_prevalence_AUC_plot, "ringed_entropy_prevalence_AUC_plot_slicing.RDS")



### 4.2. separated spes AMD ------------------------------------------------------

# Read separated_AMD_df
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_AMD_df <- read.table("separated_AMD_df.csv")

# Read separated_slices_AMD_df
setwd("~/Objects/unsupervised/separated_spes/analysis_2D")
separated_slices_AMD_df <- read.table("separated_slices_AMD_df.csv")

# Plots
separated_A_AMD_plot <- plot_AMD_metric(separated_A_spes_table, separated_AMD_df, separated_slices_AMD_df, "distance")

# separated_B_AMD_plot <- plot_AMD_metric(separated_B_spes_table, separated_AMD_df, separated_slices_AMD_df, "distance")

setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(separated_A_AMD_plot, "separated_A_AMD_plot_slicing.RDS")

# saveRDS(separated_B_AMD_plot, "separated_B_AMD_plot_slicing.RDS")

### 4.3. separated spes MS, NMS, ACINP, AE -----------------------------------

# Read separated MS, NMS, ACINP, AE dfs
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_MS_df <- read.table("separated_MS_df.csv")
separated_NMS_df <- read.table("separated_NMS_df.csv")
separated_ACINP_df <- read.table("separated_ACINP_df.csv")
separated_AE_df <- read.table("separated_AE_df.csv")

# Read separated MS, NMS, ACINP, AE dfs (slices)
setwd("~/Objects/unsupervised/separated_spes/analysis_2D")
separated_slices_MS_df <- read.table("separated_slices_MS_df.csv")
separated_slices_NMS_df <- read.table("separated_slices_NMS_df.csv")
separated_slices_ACINP_df <- read.table("separated_slices_ACINP_df.csv")
separated_slices_AE_df <- read.table("separated_slices_AE_df.csv")

# Plots
separated_A_MS_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_MS_df, separated_slices_MS_df, "MS", "distance")
separated_A_NMS_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_NMS_df, separated_slices_NMS_df, "NMS", "distance")
separated_A_ACINP_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_ACINP_df, separated_slices_ACINP_df, "ACINP", "distance")
separated_A_AE_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_AE_df, separated_slices_AE_df, "AE", "distance")

separated_A_MS_box_plot <- plot_gradient_metrics_type1_boxplot(separated_A_spes_table, separated_MS_df, separated_slices_MS_df, "MS", "distance")
separated_A_NMS_box_plot <- plot_gradient_metrics_type1_boxplot(separated_A_spes_table, separated_NMS_df, separated_slices_NMS_df, "NMS", "distance")
separated_A_ACINP_box_plot <- plot_gradient_metrics_type1_boxplot(separated_A_spes_table, separated_ACINP_df, separated_slices_ACINP_df, "ACINP", "distance")
separated_A_AE_box_plot <- plot_gradient_metrics_type1_boxplot(separated_A_spes_table, separated_AE_df, separated_slices_AE_df, "AE", "distance")

# separated_B_MS_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_MS_df, separated_slices_MS_df, "MS", "distance")
# separated_B_NMS_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_NMS_df, separated_slices_NMS_df, "NMS", "distance")
# separated_B_ACINP_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_ACINP_df, separated_slices_ACINP_df, "ACINP", "distance")
# separated_B_AE_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_AE_df, separated_slices_AE_df, "AE", "distance")
# 
# separated_B_MS_box_plot <- plot_gradient_metrics_type1_boxplot(separated_B_spes_table, separated_MS_df, separated_slices_MS_df, "MS", "distance")
# separated_B_NMS_box_plot <- plot_gradient_metrics_type1_boxplot(separated_B_spes_table, separated_NMS_df, separated_slices_NMS_df, "NMS", "distance")
# separated_B_ACINP_box_plot <- plot_gradient_metrics_type1_boxplot(separated_B_spes_table, separated_ACINP_df, separated_slices_ACINP_df, "ACINP", "distance")
# separated_B_AE_box_plot <- plot_gradient_metrics_type1_boxplot(separated_B_spes_table, separated_AE_df, separated_slices_AE_df, "AE", "distance")

setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(separated_A_MS_plot, "separated_A_MS_plot_slicing.RDS")
saveRDS(separated_A_NMS_plot, "separated_A_NMS_plot_slicing.RDS")
saveRDS(separated_A_ACINP_plot, "separated_A_ACINP_plot_slicing.RDS")
saveRDS(separated_A_AE_plot, "separated_A_AE_plot_slicing.RDS")

saveRDS(separated_A_MS_box_plot, "separated_A_MS_box_plot_slicing.RDS")
saveRDS(separated_A_NMS_box_plot, "separated_A_NMS_box_plot_slicing.RDS")
saveRDS(separated_A_ACINP_box_plot, "separated_A_ACINP_box_plot_slicing.RDS")
saveRDS(separated_A_AE_box_plot, "separated_A_AE_box_plot_slicing.RDS")

# saveRDS(separated_B_MS_plot, "separated_B_MS_plot_slicing.RDS")
# saveRDS(separated_B_NMS_plot, "separated_B_NMS_plot_slicing.RDS")
# saveRDS(separated_B_ACINP_plot, "separated_B_ACINP_plot_slicing.RDS")
# saveRDS(separated_B_AE_plot, "separated_B_AE_plot_slicing.RDS")
# 
# saveRDS(separated_B_MS_box_plot, "separated_B_MS_box_plot_slicing.RDS")
# saveRDS(separated_B_NMS_box_plot, "separated_B_NMS_box_plot_slicing.RDS")
# saveRDS(separated_B_ACINP_box_plot, "separated_B_ACINP_box_plot_slicing.RDS")
# saveRDS(separated_B_AE_box_plot, "separated_B_AE_box_plot_slicing.RDS")

### 4.4. separated spes ACIN, CKR ------------------------------------------------

# Read separated ACIN, CKR
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_ACIN_df <- read.table("separated_ACIN_df.csv")
separated_CKR_df <- read.table("separated_CKR_df.csv")

# Read separated ACIN, CKR (slices)
setwd("~/Objects/unsupervised/separated_spes/analysis_2D")
separated_slices_ACIN_df <- read.table("separated_slices_ACIN_df.csv")
separated_slices_CKR_df <- read.table("separated_slices_CKR_df.csv")

# Plots
separated_A_ACIN_plot <- plot_gradient_metrics_type2(separated_A_spes_table, separated_ACIN_df, separated_slices_ACIN_df, "ACIN", "distance", 20, 100)
separated_A_CKR_plot <- plot_gradient_metrics_type2(separated_A_spes_table, separated_CKR_df, separated_slices_CKR_df, "CKR", "distance", 20, 100)

separated_A_ACIN_box_plot <- plot_gradient_metrics_type2_boxplot(separated_A_spes_table, separated_ACIN_df, separated_slices_ACIN_df, "ACIN", "distance", 20, 100)
separated_A_CKR_box_plot <- plot_gradient_metrics_type2_boxplot(separated_A_spes_table, separated_CKR_df, separated_slices_CKR_df, "CKR", "distance", 20, 100)

# separated_B_ACIN_plot <- plot_gradient_metrics_type2(separated_B_spes_table, separated_ACIN_df, separated_slices_ACIN_df, "ACIN", "distance", 20, 100)
# separated_B_CKR_plot <- plot_gradient_metrics_type2(separated_B_spes_table, separated_CKR_df, separated_slices_CKR_df, "CKR", "distance", 20, 100)
# 
# separated_B_ACIN_box_plot <- plot_gradient_metrics_type2_boxplot(separated_B_spes_table, separated_ACIN_df, separated_slices_ACIN_df, "ACIN", "distance", 20, 100)
# separated_B_CKR_box_plot <- plot_gradient_metrics_type2_boxplot(separated_B_spes_table, separated_CKR_df, separated_slices_CKR_df, "CKR", "distance", 20, 100)

setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(separated_A_ACIN_plot, "separated_A_ACIN_plot_slicing.RDS")
saveRDS(separated_A_CKR_plot, "separated_A_CKR_plot_slicing.RDS")

saveRDS(separated_A_ACIN_box_plot, "separated_A_ACIN_box_plot_slicing.RDS")
saveRDS(separated_A_CKR_box_plot, "separated_A_CKR_box_plot_slicing.RDS")

# saveRDS(separated_B_ACIN_plot, "separated_B_ACIN_plot_slicing.RDS")
# saveRDS(separated_B_CKR_plot, "separated_B_CKR_plot_slicing.RDS")
# 
# saveRDS(separated_B_ACIN_box_plot, "separated_B_ACIN_plot_slicing.RDS")
# saveRDS(separated_B_CKR_box_plot, "separated_B_CKR_plot_slicing.RDS")

### 4.5. separated spes SAC ------------------------------------------------------

# Read separated_SAC_df
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_prop_SAC_df <- read.table("separated_prop_SAC_df.csv")
separated_entropy_SAC_df <- read.table("separated_entropy_SAC_df.csv")

# Read separated_SAC_df (slices)
setwd("~/Objects/unsupervised/separated_spes/analysis_2D")
separated_slices_prop_SAC_df <- read.table("separated_slices_prop_SAC_df.csv")
separated_slices_entropy_SAC_df <- read.table("separated_slices_entropy_SAC_df.csv")

# Plots
separated_A_prop_SAC_plot <- plot_proportion_SAC(separated_A_spes_table, separated_prop_SAC_df, separated_slices_prop_SAC_df, "distance")
separated_A_entropy_SAC_plot <- plot_entropy_SAC(separated_A_spes_table, separated_entropy_SAC_df, separated_slices_entropy_SAC_df, "distance")

# separated_B_prop_SAC_plot <- plot_proportion_SAC(separated_B_spes_table, separated_prop_SAC_df, separated_slices_prop_SAC_df, "distance")
# separated_B_entropy_SAC_plot <- plot_entropy_SAC(separated_B_spes_table, separated_entropy_SAC_df, separated_slices_entropy_SAC_df, "distance")

setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(separated_A_prop_SAC_plot, "separated_A_prop_SAC_plot_slicing.RDS")
saveRDS(separated_A_entropy_SAC_plot, "separated_A_entropy_SAC_plot_slicing.RDS")

# saveRDS(separated_B_prop_SAC_plot, "separated_B_prop_SAC_plot_slicing.RDS")
# saveRDS(separated_B_entropy_SAC_plot, "separated_B_entropy_SAC_plot_slicing.RDS")

### 4.6. separated spes prevalence ------------------------------------------------

# Read separated prevalence dfs
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_prop_prevalence_df <- read.table("separated_prop_prevalence_df.csv")
separated_entropy_prevalence_df <- read.table("separated_entropy_prevalence_df.csv")

# Read separated prevalence dfs (slices)
setwd("~/Objects/unsupervised/separated_spes/analysis_2D")
separated_slices_prop_prevalence_df <- read.table("separated_slices_prop_prevalence_df.csv")
separated_slices_entropy_prevalence_df <- read.table("separated_slices_entropy_prevalence_df.csv")

# Plots
separated_A_prop_prevalence_plot <- plot_proportion_prevalence(separated_A_spes_table, separated_prop_prevalence_df, separated_slices_prop_prevalence_df, "distance")
separated_A_entropy_prevalence_plot <- plot_entropy_prevalence(separated_A_spes_table, separated_entropy_prevalence_df, separated_slices_entropy_prevalence_df, "distance")

separated_A_prop_prevalence_box_plot <- plot_proportion_prevalence_boxplot(separated_A_spes_table, separated_prop_prevalence_df, separated_slices_prop_prevalence_df, "distance")
separated_A_entropy_prevalence_box_plot <- plot_entropy_prevalence_boxplot(separated_A_spes_table, separated_entropy_prevalence_df, separated_slices_entropy_prevalence_df, "distance")

separated_A_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(separated_A_spes_table, separated_prop_prevalence_df, separated_slices_prop_prevalence_df, "distance")
separated_A_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(separated_A_spes_table, separated_entropy_prevalence_df, separated_slices_entropy_prevalence_df, "distance")

# separated_B_prop_prevalence_plot <- plot_proportion_prevalence(separated_B_spes_table, separated_prop_prevalence_df, separated_slices_prop_prevalence_df, "distance")
# separated_B_entropy_prevalence_plot <- plot_entropy_prevalence(separated_B_spes_table, separated_entropy_prevalence_df, separated_slices_entropy_prevalence_df, "distance")
# 
# separated_B_prop_prevalence_box_plot <- plot_proportion_prevalence_boxplot(separated_B_spes_table, separated_prop_prevalence_df, separated_slices_prop_prevalence_df, "distance")
# separated_B_entropy_prevalence_box_plot <- plot_entropy_prevalence_boxplot(separated_B_spes_table, separated_entropy_prevalence_df, separated_slices_entropy_prevalence_df, "distance")
# 
# separated_B_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(separated_B_spes_table, separated_prop_prevalence_df, separated_slices_prop_prevalence_df, "distance")
# separated_B_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(separated_B_spes_table, separated_entropy_prevalence_df, separated_slices_entropy_prevalence_df, "distance")

setwd("~/Objects/unsupervised/plots/slicing")
saveRDS(separated_A_prop_prevalence_plot, "separated_A_prop_prevalence_plot_slicing.RDS")
saveRDS(separated_A_entropy_prevalence_plot, "separated_A_entropy_prevalence_plot_slicing.RDS")
saveRDS(separated_A_prop_prevalence_box_plot, "separated_A_prop_prevalence_box_plot_slicing.RDS")
saveRDS(separated_A_entropy_prevalence_box_plot, "separated_A_entropy_prevalence_box_plot_slicing.RDS")
saveRDS(separated_A_prop_prevalence_AUC_plot, "separated_A_prop_prevalence_AUC_plot_slicing.RDS")
saveRDS(separated_A_entropy_prevalence_AUC_plot, "separated_A_entropy_prevalence_AUC_plot_slicing.RDS")

# saveRDS(separated_B_prop_prevalence_plot, "separated_B_prop_prevalence_plot_slicing.RDS")
# saveRDS(separated_B_entropy_prevalence_plot, "separated_B_entropy_prevalence_plot_slicing.RDS")
# saveRDS(separated_B_prop_prevalence_box_plot, "separated_B_prop_prevalence_box_plot_slicing.RDS")
# saveRDS(separated_B_entropy_prevalence_box_plot, "separated_B_entropy_prevalence_box_plot_slicing.RDS")
# saveRDS(separated_B_prop_prevalence_AUC_plot, "separated_B_prop_prevalence_AUC_plot_slicing.RDS")
# saveRDS(separated_B_entropy_prevalence_AUC_plot, "separated_B_entropy_prevalence_AUC_plot_slicing.RDS")

### 5.0. pdf with all plots -----------------------------
setwd("~/Objects/unsupervised/plots/slicing")
metrics <- c("AMD", 
             "MS", "MS_box", "NMS", "NMS_box", "ACINP", "ACINP_box", "AE", "AE_box", 
             "ACIN", "ACIN_box", "CKR", "CKR_box", 
             "prop_SAC", "prop_prevalence", "prop_prevalence_box", "prop_prevalence_AUC", 
             "entropy_SAC", "entropy_prevalence", "entropy_prevalence_box","entropy_prevalence_AUC")
arrangements <- c("mixed", "ringed", "separated_A") #ignoring separated_B

pdf("plots.pdf", width = 15, height = 10)

for (metric in metrics) {
  for (arrangement in arrangements) {
    plot_file_name <- paste(arrangement, "_", metric, "_plot_slicing.RDS", sep = "")
    print(readRDS(plot_file_name))
  }
}
dev.off()

### Without background noise cells functions ---------------------------------
plot_AMD_metric <- function(spes_table, AMD_df, slices_AMD_df, arrangement_colname) {
  
  # AMD pairs are A/A, A/B, B/A, B/B
  AMD_pairs <- data.frame(cell1 = c("A", "A", "B", "B"),
                          cell2 = c("A", "B", "A", "B"))
  AMD_pairs$pair <- paste(AMD_pairs$cell1, AMD_pairs$cell2, sep = "/")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(nrow(AMD_pairs))) {
    
    # Subset AMD_df for chosen pair
    AMD_df_subset <- AMD_df[AMD_df$reference == AMD_pairs[i, "cell1"] & AMD_df$target == AMD_pairs[i, "cell2"], ]
    
    # Subset slices_AMD_df for chosen pair
    slices_AMD_df_subset <- slices_AMD_df[slices_AMD_df$reference == AMD_pairs[i, "cell1"] & slices_AMD_df$target == AMD_pairs[i, "cell2"], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_AMD_df_subset, AMD_df_subset, "spe", suffix = c("_2D", "_3D"))
    
    slices_AMD_df_subset$AMD <- (joint_df$AMD_2D - joint_df$AMD_3D) / joint_df$AMD_3D
    
    # Combine spes_table and AMD_df
    plot_df <- cbind(spes_table, slices_AMD_df_subset)
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    
    fig_slice <- ggplot(plot_df, aes(!!sym(arrangement_colname), AMD, col = slice)) +
      geom_point() +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
    # fig_bg_prop_A <- ggplot(plot_df, aes(!!sym(arrangement_colname), AMD, col = bg_prop_A)) +
    #   geom_point() +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    # 
    # fig_bg_prop_B <- ggplot(plot_df, aes(!!sym(arrangement_colname), AMD, col = bg_prop_B)) +
    #   geom_point() +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
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
    
    all_plots_list[[AMD_pairs[i, "pair"]]] <- list(slice = fig_slice + theme(legend.position="none"),
                                                   # bg_prop_A = fig_bg_prop_A + theme(legend.position="none"),
                                                   # bg_prop_B = fig_bg_prop_B + theme(legend.position="none"),
                                                   shape = fig_shape + theme(legend.position="none"),
                                                   variation_E = fig_variation_E + theme(legend.position="none"),
                                                   volume_E = fig_volume_E + theme(legend.position="none"),
                                                   width_N = fig_width_N + theme(legend.position="none"))
  }
  
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  # legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  # legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       # legend_bg_prop_a,
                       # legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by pairs
  plots_pair_list <- list()
  
  for (i in seq(nrow(AMD_pairs))) {
    pair <- AMD_pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$slice,
                       # all_plots_list[[pair]]$bg_prop_A,
                       # all_plots_list[[pair]]$bg_prop_B,
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

plot_gradient_metrics_type1 <- function(spes_table, gradient_metric_df, slices_gradient_metric_df, metric, arrangement_colname) {
  
  # Constants
  cell_types <- c("A", "B") # Use A as reference, and B as target, and vice versa
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(length(cell_types))) {
    # Subset gradient_metric_df for current reference cell
    gradient_metric_df_subset <- gradient_metric_df[gradient_metric_df$reference == cell_types[i], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_gradient_metric_df_subset <- slices_gradient_metric_df[slices_gradient_metric_df$reference == cell_types[i], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_gradient_metric_df_subset, gradient_metric_df_subset, by = "spe", suffix = c("_2D", "_3D"))
    
    for (radii_colname in radii_colnames) {
      slices_gradient_metric_df_subset[ , radii_colname] <- 
        (joint_df[ , paste(radii_colname, "_2D", sep = "")] - joint_df[ , paste(radii_colname, "_3D", sep = "")]) / joint_df[ , paste(radii_colname, "_3D", sep = "")]
    }
    
    # Combine spes_table and slices_gradient_metric_df_subset
    plot_df <- cbind(spes_table, slices_gradient_metric_df_subset)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , radii_colnames)
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    plot_df$key <- paste(plot_df$spe, plot_df$slice, sep = "_")
    
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df$variable <- unfactor(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    
    fig_slice <- ggplot(plot_df, aes(variable, value, group = key, col = slice)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = key, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    # fig_bg_prop_A <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_A)) +
    #   geom_line() +
    #   labs(x = "radius", y = metric) +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    # 
    # fig_bg_prop_B <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_B)) +
    #   geom_line() +
    #   labs(x = "radius", y = metric) +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = key, col = shape)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100

    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = key, col = variation_E)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = key, col = volume_E)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)

    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = key, col = width_N)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    all_plots_list[[cell_types[i]]] <- list(slice = fig_slice + theme(legend.position = "none"),
                                            arrangement = fig_arrangement + theme(legend.position = "none"), 
                                            # bg_prop_A = fig_bg_prop_A + theme(legend.position = "none"),
                                            # bg_prop_B = fig_bg_prop_B + theme(legend.position = "none"),
                                            shape = fig_shape + theme(legend.position = "none"),
                                            variation_E = fig_variation_E + theme(legend.position = "none"),
                                            volume_E = fig_volume_E + theme(legend.position = "none"),
                                            width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  # legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  # legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       legend_arrangement,
                       # legend_bg_prop_a,
                       # legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference cell type
  plots_ref_list <- list()
  
  for (i in seq(length(cell_types))) {
    
    reference_cell_type <- cell_types[i]
    
    plots <- plot_grid(all_plots_list[[reference_cell_type]]$slice,
                       all_plots_list[[reference_cell_type]]$arrangement,
                       # all_plots_list[[reference_cell_type]]$bg_prop_A,
                       # all_plots_list[[reference_cell_type]]$bg_prop_B,
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

plot_gradient_metrics_type1_boxplot <- function(spes_table, gradient_metric_df, slices_gradient_metric_df, metric, arrangement_colname) {
  
  # Constants
  cell_types <- c("A", "B") # Use A as reference, and B as target, and vice versa
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(length(cell_types))) {
    # Subset gradient_metric_df for current reference cell
    gradient_metric_df_subset <- gradient_metric_df[gradient_metric_df$reference == cell_types[i], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_gradient_metric_df_subset <- slices_gradient_metric_df[slices_gradient_metric_df$reference == cell_types[i], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_gradient_metric_df_subset, gradient_metric_df_subset, by = "spe", suffix = c("_2D", "_3D"))
    
    for (radii_colname in radii_colnames) {
      slices_gradient_metric_df_subset[ , radii_colname] <- 
        (joint_df[ , paste(radii_colname, "_2D", sep = "")] - joint_df[ , paste(radii_colname, "_3D", sep = "")]) / joint_df[ , paste(radii_colname, "_3D", sep = "")]
    }
    
    # Combine spes_table and slices_gradient_metric_df_subset
    plot_df <- cbind(spes_table, slices_gradient_metric_df_subset)
    
    plot_df$slice <- as.character(plot_df$slice)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , radii_colnames)
    
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df$variable <- unfactor(plot_df$variable)
    plot_df$variable <- substr(plot_df$variable, 2, nchar(plot_df$variable))
    plot_df$variable <- factor(plot_df$variable, as.character(radii))
    
    fig_boxplot <- ggplot(plot_df, aes(variable, value, col = slice)) +
      geom_boxplot() +
      theme_bw() +
      facet_wrap(~slice, ncol = 5) +
      scale_color_manual(values = viridis::viridis(5)) +
      labs(x = "radius", y = metric)
    
    all_plots_list[[cell_types[i]]] <- list(fig_boxplot = fig_boxplot)
  }
  
  # Combine the plots together by reference cell type
  plots_ref_list <- list()
  
  for (i in seq(length(cell_types))) {
    reference_cell_type <- cell_types[i]
    
    plots <- plot_grid(all_plots_list[[reference_cell_type]]$fig_boxplot,
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
                             nrow = 2, ncol = 1)
  
  methods::show(combined_plot)
  
  return(combined_plot)
}

plot_gradient_metrics_type2 <- function(spes_table, gradient_metric_df, slices_gradient_metric_df, metric, arrangement_colname, min_radius, max_radius) {
  
  # Constants
  pairs <- data.frame(cell1 = c("A", "A", "B", "B"),
                      cell2 = c("A", "B", "A", "B"))
  pairs$pair <- paste(pairs$cell1, pairs$cell2, sep = "/")
  
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(nrow(pairs))) {
    
    # Subset gradient_metric_df for current reference cell
    gradient_metric_df_subset <- gradient_metric_df[gradient_metric_df$reference == pairs[i, "cell1"] &
                                                      gradient_metric_df$target == pairs[i, "cell2"], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_gradient_metric_df_subset <- slices_gradient_metric_df[slices_gradient_metric_df$reference == pairs[i, "cell1"] &
                                                                    slices_gradient_metric_df$target == pairs[i, "cell2"], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_gradient_metric_df_subset, gradient_metric_df_subset, by = "spe", suffix = c("_2D", "_3D"))
    
    for (radii_colname in radii_colnames) {
      slices_gradient_metric_df_subset[ , radii_colname] <- 
        (joint_df[ , paste(radii_colname, "_2D", sep = "")] - joint_df[ , paste(radii_colname, "_3D", sep = "")]) / joint_df[ , paste(radii_colname, "_3D", sep = "")]
    }
    
    
    # Combine spes_table and mixed_AMD_df
    plot_df <- cbind(spes_table, slices_gradient_metric_df_subset)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , radii_colnames)
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    plot_df$key <- paste(plot_df$spe, plot_df$slice, sep = "_")
    
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df$variable <- unfactor(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    
    plot_df <- plot_df[plot_df$variable >= min_radius & plot_df$variable <= max_radius, ]
    
    fig_slice <- ggplot(plot_df, aes(variable, value, group = key, col = slice)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = key, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    # fig_bg_prop_A <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_A)) +
    #   geom_line() +
    #   labs(x = "radius", y = metric) +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    # 
    # fig_bg_prop_B <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_B)) +
    #   geom_line() +
    #   labs(x = "radius", y = metric) +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = key, col = shape)) +
      geom_line() +
      theme_bw()

    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100

    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = key, col = variation_E)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()

    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = key, col = volume_E)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)

    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = key, col = width_N)) +
      geom_line() +
      labs(x = "radius", y = metric) +
      theme_bw()
    
    all_plots_list[[pairs[i, "pair"]]] <- list(slice = fig_slice + theme(legend.position = "none"),
                                               arrangement = fig_arrangement + theme(legend.position = "none"),
                                               # bg_prop_A = fig_bg_prop_A + theme(legend.position = "none"),
                                               # bg_prop_B = fig_bg_prop_B + theme(legend.position = "none"),
                                               shape = fig_shape + theme(legend.position = "none"),
                                               variation_E = fig_variation_E + theme(legend.position = "none"),
                                               volume_E = fig_volume_E + theme(legend.position = "none"),
                                               width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  # legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  # legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       legend_arrangement,
                       # legend_bg_prop_a,
                       # legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference cell type
  plots_pair_list <- list()
  
  for (i in seq(nrow(pairs))) {
    pair <- pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$slice,
                       all_plots_list[[pair]]$arrangement,
                       # all_plots_list[[pair]]$bg_prop_A,
                       # all_plots_list[[pair]]$bg_prop_B,
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

plot_gradient_metrics_type2_boxplot <- function(spes_table, gradient_metric_df, slices_gradient_metric_df, metric, arrangement_colname, min_radius, max_radius) {
  
  # Constants
  pairs <- data.frame(cell1 = c("A", "A", "B", "B"),
                      cell2 = c("A", "B", "A", "B"))
  pairs$pair <- paste(pairs$cell1, pairs$cell2, sep = "/")
  
  radii <- seq(20, 100, 10)
  radii_colnames <- paste("r", radii, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(nrow(pairs))) {
    
    # Subset gradient_metric_df for current reference cell
    gradient_metric_df_subset <- gradient_metric_df[gradient_metric_df$reference == pairs[i, "cell1"] &
                                                      gradient_metric_df$target == pairs[i, "cell2"], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_gradient_metric_df_subset <- slices_gradient_metric_df[slices_gradient_metric_df$reference == pairs[i, "cell1"] &
                                                                    slices_gradient_metric_df$target == pairs[i, "cell2"], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_gradient_metric_df_subset, gradient_metric_df_subset, by = "spe", suffix = c("_2D", "_3D"))
    
    for (radii_colname in radii_colnames) {
      slices_gradient_metric_df_subset[ , radii_colname] <- 
        (joint_df[ , paste(radii_colname, "_2D", sep = "")] - joint_df[ , paste(radii_colname, "_3D", sep = "")]) / joint_df[ , paste(radii_colname, "_3D", sep = "")]
    }
    
    
    # Combine spes_table and mixed_AMD_df
    plot_df <- cbind(spes_table, slices_gradient_metric_df_subset)
    
    plot_df$slice <- as.character(plot_df$slice)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , radii_colnames)
    
    # Extract radius value from radius strings (r1 -> 1, r2 -> 2...)
    plot_df$variable <- unfactor(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    plot_df <- plot_df[plot_df$variable >= min_radius & plot_df$variable <= max_radius, ]
    plot_df$variable <- factor(as.character(plot_df$variable), radii)
    
    fig_boxplot <- ggplot(plot_df, aes(variable, value, col = slice)) +
      geom_boxplot() +
      theme_bw() +
      facet_wrap(~slice, ncol = 5) +
      scale_color_manual(values = viridis::viridis(5)) +
      labs(x = "radius", y = metric)
    
    all_plots_list[[pairs[i, "pair"]]] <- list(fig_boxplot = fig_boxplot)
  }
  # Combine the plots together by reference cell type
  plots_pair_list <- list()
  
  for (i in seq(nrow(pairs))) {
    pair <- pairs[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$fig_boxplot,
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
                             nrow = 4, ncol = 1)
  
  methods::show(combined_plot)
  
  return(combined_plot)
}

plot_proportion_SAC <- function(spes_table, SAC_df, slices_SAC_df, arrangement_colname) {
  
  # Get possible reference and target cell combinations
  prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
  prop_cell_types$pair <- paste(prop_cell_types$ref, prop_cell_types$tar, sep = "/")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    
    # Subset SAC_df for chosen pair
    SAC_df_subset <- SAC_df[SAC_df$reference == prop_cell_types$ref[i], ]
    
    # Subset slices_SAC_df for chosen pair
    slices_SAC_df_subset <- slices_SAC_df[slices_SAC_df$reference == prop_cell_types$ref[i], ]
    
    # Get difference between SAC values in 3D and 2D slices.
    joint_df <- full_join(slices_SAC_df_subset, SAC_df_subset, "spe", suffix = c("_2D", "_3D"))
    
    slices_SAC_df_subset$SAC <- (joint_df$proportion_2D - joint_df$proportion_3D) / joint_df$proportion_3D
    
    # Combine spes_table and SAC_df
    plot_df <- cbind(spes_table, slices_SAC_df_subset)
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    
    fig_slice <- ggplot(plot_df, aes(!!sym(arrangement_colname), proportion, col = slice)) +
      geom_point() +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
    # fig_bg_prop_A <- ggplot(plot_df, aes(!!sym(arrangement_colname), proportion, col = bg_prop_A)) +
    #   geom_point() +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    # 
    # fig_bg_prop_B <- ggplot(plot_df, aes(!!sym(arrangement_colname), proportion, col = bg_prop_B)) +
    #   geom_point() +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
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
    
    all_plots_list[[prop_cell_types[i, "pair"]]] <- list(slice = fig_slice + theme(legend.position="none"),
                                                         # bg_prop_A = fig_bg_prop_A + theme(legend.position="none"),
                                                         # bg_prop_B = fig_bg_prop_B + theme(legend.position="none"),
                                                         shape = fig_shape + theme(legend.position="none"),
                                                         variation_E = fig_variation_E + theme(legend.position="none"),
                                                         volume_E = fig_volume_E + theme(legend.position="none"),
                                                         width_N = fig_width_N + theme(legend.position="none"))
  }
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  # legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  # legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       # legend_bg_prop_a, 
                       # legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by pairs
  
  plots_pair_list <- list()
  
  for (i in seq(nrow(prop_cell_types))) {
    pair <- prop_cell_types[i, "pair"]
    
    plots <- plot_grid(all_plots_list[[pair]]$slice,
                       # all_plots_list[[pair]]$bg_prop_A,
                       # all_plots_list[[pair]]$bg_prop_B,
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

plot_entropy_SAC <- function(spes_table, SAC_df, slices_SAC_df, arrangement_colname) {
  
  # Get possible cell type of interest combinations
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    
    # Subset SAC_df for chosen pair
    SAC_df_subset <- SAC_df[SAC_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Subset slices_SAC_df for chosen pair
    slices_SAC_df_subset <- slices_SAC_df[slices_SAC_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Get difference between SAC values in 3D and 2D slices.
    joint_df <- full_join(slices_SAC_df_subset, SAC_df_subset, "spe", suffix = c("_2D", "_3D"))
    
    slices_SAC_df_subset$SAC <- (joint_df$entropy_2D - joint_df$entropy_3D) / joint_df$entropy_3D
    
    # Combine spes_table and SAC_df
    plot_df <- cbind(spes_table, slices_SAC_df_subset)
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    
    fig_slice <- ggplot(plot_df, aes(!!sym(arrangement_colname), entropy, col = slice)) +
      geom_point() +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
    # fig_bg_prop_A <- ggplot(plot_df, aes(!!sym(arrangement_colname), entropy, col = bg_prop_A)) +
    #   geom_point() +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    # 
    # fig_bg_prop_B <- ggplot(plot_df, aes(!!sym(arrangement_colname), entropy, col = bg_prop_B)) +
    #   geom_point() +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
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
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(slice = fig_slice + theme(legend.position="none"),
                                                               # bg_prop_A = fig_bg_prop_A + theme(legend.position="none"), 
                                                               # bg_prop_B = fig_bg_prop_B + theme(legend.position="none"),
                                                               shape = fig_shape + theme(legend.position="none"),
                                                               variation_E = fig_variation_E + theme(legend.position="none"),
                                                               volume_E = fig_volume_E + theme(legend.position="none"),
                                                               width_N = fig_width_N + theme(legend.position="none"))
  }
  
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  # legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  # legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       # legend_bg_prop_a, 
                       # legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by cell types of interest
  
  plots_cell_types_list <- list()
  
  for (i in seq(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$slice,
                       # all_plots_list[[cell_types]]$bg_prop_A,
                       # all_plots_list[[cell_types]]$bg_prop_B,
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

plot_proportion_prevalence <- function(spes_table, prevalence_df, slices_prevalence_df, arrangement_colname) {
  
  # Constants
  prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
  prop_cell_types$pair <- paste(prop_cell_types$ref, prop_cell_types$tar, sep = "/")
  
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    
    # Subset gradient_metric_df for current reference cell
    prevalence_df_subset <- prevalence_df[prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_prevalence_df_subset <- slices_prevalence_df[slices_prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_prevalence_df_subset, prevalence_df_subset, by = "spe", suffix = c("_2D", "_3D"))
    
    for (threshold_colname in threshold_colnames) {
      slices_prevalence_df_subset[ , threshold_colname] <- 
        (joint_df[ , paste(threshold_colname, "_2D", sep = "")] - joint_df[ , paste(threshold_colname, "_3D", sep = "")]) / joint_df[ , paste(threshold_colname, "_3D", sep = "")]
    }
    
    # Combine spes_table and prevalence_df
    plot_df <- cbind(spes_table, slices_prevalence_df_subset)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , threshold_colnames)
    
    # Extract threshold value from threshold strings (t0.01 -> 0.01...)
    plot_df$variable <- as.character(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    plot_df$key <- paste(plot_df$spe, plot_df$slice, sep = "_")
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = key, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    fig_slice <- ggplot(plot_df, aes(variable, value, group = key, col = slice)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
    # fig_bg_prop_A <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_A)) +
    #   geom_line() +
    #   labs(x = "threshold", y = "prevalence") +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    # 
    # fig_bg_prop_B <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_B)) +
    #   geom_line() +
    #   labs(x = "threshold", y = "prevalence") +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = key, col = shape)) +
      labs(x = "threshold", y = "prevalence") +
      geom_line() +
      theme_bw()

    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100

    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = key, col = variation_E)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()

    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = key, col = volume_E)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)

    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = key, col = width_N)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(slice = fig_slice + theme(legend.position = "none"),
                                                      arrangement = fig_arrangement + theme(legend.position = "none"), 
                                                      # bg_prop_A = fig_bg_prop_A + theme(legend.position = "none"),
                                                      # bg_prop_B = fig_bg_prop_B + theme(legend.position = "none"),
                                                      shape = fig_shape + theme(legend.position = "none"),
                                                      variation_E = fig_variation_E + theme(legend.position = "none"),
                                                      volume_E = fig_volume_E + theme(legend.position = "none"),
                                                      width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  # legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  # legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       legend_arrangement,
                       # legend_bg_prop_a,
                       # legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  
  # Combine the plots together by reference (and target) cell pairs
  plots_pair_list <- list()
  
  for (i in seq(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$slice,
                       all_plots_list[[pair]]$arrangement,
                       # all_plots_list[[pair]]$bg_prop_A,
                       # all_plots_list[[pair]]$bg_prop_B,
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

plot_proportion_prevalence_boxplot <- function(spes_table, prevalence_df, slices_prevalence_df, arrangement_colname) {
  
  # Constants
  prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
  prop_cell_types$pair <- paste(prop_cell_types$ref, prop_cell_types$tar, sep = "/")
  
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    
    # Subset gradient_metric_df for current reference cell
    prevalence_df_subset <- prevalence_df[prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_prevalence_df_subset <- slices_prevalence_df[slices_prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_prevalence_df_subset, prevalence_df_subset, by = "spe", suffix = c("_2D", "_3D"))
    
    for (threshold_colname in threshold_colnames) {
      slices_prevalence_df_subset[ , threshold_colname] <- 
        (joint_df[ , paste(threshold_colname, "_2D", sep = "")] - joint_df[ , paste(threshold_colname, "_3D", sep = "")]) / joint_df[ , paste(threshold_colname, "_3D", sep = "")]
    }
    
    # Combine spes_table and prevalence_df
    plot_df <- cbind(spes_table, slices_prevalence_df_subset)
    
    plot_df$slice <- as.character(plot_df$slice)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , threshold_colnames)
    
    # Extract threshold value from threshold strings (t0.01 -> 0.01...)
    plot_df$variable <- as.character(plot_df$variable)
    plot_df$variable <- substr(plot_df$variable, 2, nchar(plot_df$variable))
    plot_df$variable <- factor(plot_df$variable, as.character(thresholds))
    
    fig_boxplot <- ggplot(plot_df, aes(variable, value, col = slice)) +
      geom_boxplot() +
      theme_bw() +
      facet_wrap(~slice, ncol = 5) +
      scale_color_manual(values = viridis::viridis(5)) +
      labs(x = "threshold", y = "prevalence") +
      scale_x_discrete(breaks = c(0.01, 1), labels = c("0.01", "1"))
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(boxplot = fig_boxplot)
  }
  
  # Combine the plots together by reference (and target) cell pairs
  plots_pair_list <- list()
  
  for (i in seq(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$boxplot,
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
                             nrow = 2, ncol = 1,
                             rel_heights = c(1, 1, 0.5))
  
  methods::show(combined_plot)
  
  return(combined_plot)
}

plot_entropy_prevalence <- function(spes_table, prevalence_df, slices_prevalence_df, arrangement_colname) {
  
  # Constants
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    
    # Subset gradient_metric_df for current reference cell
    prevalence_df_subset <- prevalence_df[prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_prevalence_df_subset <- slices_prevalence_df[slices_prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_prevalence_df_subset, prevalence_df_subset, by = "spe", suffix = c("_2D", "_3D"))
    
    for (threshold_colname in threshold_colnames) {
      slices_prevalence_df_subset[ , threshold_colname] <- 
        (joint_df[ , paste(threshold_colname, "_2D", sep = "")] - joint_df[ , paste(threshold_colname, "_3D", sep = "")]) / joint_df[ , paste(threshold_colname, "_3D", sep = "")]
    }
    
    # Combine spes_table and prevalence_df
    plot_df <- cbind(spes_table, slices_prevalence_df_subset)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , threshold_colnames)
    
    # Extract threshold value from threshold strings (t0.01 -> 0.01...)
    plot_df$variable <- as.character(plot_df$variable)
    plot_df$variable <- substr(plot_df$variable, 2, nchar(plot_df$variable))
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    plot_df$key <- paste(plot_df$spe, plot_df$slice, sep = "_")
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = key, col = !!sym(arrangement_colname))) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    fig_slice <- ggplot(plot_df, aes(variable, value, group = key, col = slice)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
    # fig_bg_prop_A <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_A)) +
    #   geom_line() +
    #   labs(x = "threshold", y = "prevalence") +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    # 
    # fig_bg_prop_B <- ggplot(plot_df, aes(variable, value, group = key, col = bg_prop_B)) +
    #   geom_line() +
    #   labs(x = "threshold", y = "prevalence") +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))

    fig_shape <- ggplot(plot_df, aes(variable, value, group = key, col = shape)) +
      labs(x = "threshold", y = "prevalence") +
      geom_line() +
      theme_bw()

    radii_E_df <- plot_df[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
    plot_df$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * plot_df$radius_z_E
    plot_df$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100

    fig_variation_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = key, col = variation_E)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()

    fig_volume_E <- ggplot(plot_df %>% filter(shape == "Ellipsoid"), aes(variable, value, group = key, col = volume_E)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw() +
      scale_color_continuous(n.breaks = 4)

    fig_width_N <- ggplot(plot_df %>% filter(!is.na(width_N)), aes(variable, value, group = key, col = width_N)) +
      geom_line() +
      labs(x = "threshold", y = "prevalence") +
      theme_bw()
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(slice = fig_slice + theme(legend.position = "none"),
                                                               arrangement = fig_arrangement + theme(legend.position = "none"), 
                                                               # bg_prop_A = fig_bg_prop_A + theme(legend.position = "none"),
                                                               # bg_prop_B = fig_bg_prop_B + theme(legend.position = "none"),
                                                               shape = fig_shape + theme(legend.position = "none"),
                                                               variation_E = fig_variation_E + theme(legend.position = "none"),
                                                               volume_E = fig_volume_E + theme(legend.position = "none"),
                                                               width_N = fig_width_N + theme(legend.position = "none"))
    
  }
  
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  legend_arrangement <- get_legend(fig_arrangement + theme(legend.direction = "horizontal"))
  # legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  # legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       legend_arrangement,
                       # legend_bg_prop_a,
                       # legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  
  # Combine the plots together by reference (and target) cell pairs
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$slice,
                       all_plots_list[[cell_types]]$arrangement,
                       # all_plots_list[[cell_types]]$bg_prop_A,
                       # all_plots_list[[cell_types]]$bg_prop_B,
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

plot_entropy_prevalence_boxplot <- function(spes_table, prevalence_df, slices_prevalence_df, arrangement_colname) {
  
  # Constants
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    
    # Subset gradient_metric_df for current reference cell
    prevalence_df_subset <- prevalence_df[prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Subset slices_gradient_metric_df for current reference cell
    slices_prevalence_df_subset <- slices_prevalence_df[slices_prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Get difference between AMD values in 3D and 2D slices.
    joint_df <- full_join(slices_prevalence_df_subset, prevalence_df_subset, by = "spe", suffix = c("_2D", "_3D"))
    
    for (threshold_colname in threshold_colnames) {
      slices_prevalence_df_subset[ , threshold_colname] <- 
        (joint_df[ , paste(threshold_colname, "_2D", sep = "")] - joint_df[ , paste(threshold_colname, "_3D", sep = "")]) / joint_df[ , paste(threshold_colname, "_3D", sep = "")]
    }
    
    # Combine spes_table and prevalence_df
    plot_df <- cbind(spes_table, slices_prevalence_df_subset)
    
    plot_df$slice <- as.character(plot_df$slice)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , threshold_colnames)
    
    # Extract threshold value from threshold strings (t0.01 -> 0.01...)
    plot_df$variable <- as.character(plot_df$variable)
    plot_df$variable <- substr(plot_df$variable, 2, nchar(plot_df$variable))
    plot_df$variable <- factor(plot_df$variable, as.character(thresholds))
    
    fig_boxplot <- ggplot(plot_df, aes(variable, value, col = slice)) +
      geom_boxplot() +
      theme_bw() +
      facet_wrap(~slice, ncol = 5) +
      scale_color_manual(values = viridis::viridis(5)) +
      labs(x = "threshold", y = "prevalence") +
      scale_x_discrete(breaks = c(0.01, 1), labels = c("0.01", "1"))
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(boxplot = fig_boxplot)
  }
  
  # Combine the plots together by reference (and target) cell pairs
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$boxplot,
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
                             nrow = 2, ncol = 1,
                             rel_heights = c(1, 1, 0.5))
  
  methods::show(combined_plot)
  
  return(combined_plot)
}

plot_proportion_prevalence_AUC <- function(spes_table, prevalence_df, slices_prevalence_df, arrangement_colname) {
  
  # Constants
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
  prop_cell_types$pair <- paste(prop_cell_types$ref, prop_cell_types$tar, sep = "/")
  
  # Get AUC for each prevalence gradient
  prevalence_df$AUC <- apply(prevalence_df[ , threshold_colnames], 1, sum) * 0.01
  prevalence_df <- prevalence_df[ , c("spe", "reference", "target", "AUC")]
  
  # Get AUC for each prevalence gradient (slices)
  slices_prevalence_df$AUC <- apply(slices_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
  slices_prevalence_df <- slices_prevalence_df[ , c("spe", "slice", "reference", "target", "AUC")]
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    
    # Subset prevalence_df for chosen pair
    prevalence_df_subset <- prevalence_df[prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Subset slices_prevalence_df for chosen pair
    slices_prevalence_df_subset <- slices_prevalence_df[slices_prevalence_df$reference == prop_cell_types$ref[i], ]
    
    # Get difference between prevalence values in 3D and 2D slices.
    joint_df <- full_join(slices_prevalence_df_subset, prevalence_df_subset, "spe", suffix = c("_2D", "_3D"))
    
    slices_prevalence_df_subset$AUC <- (joint_df$AUC_2D - joint_df$AUC_3D) / joint_df$AUC_3D
    
    # Combine spes_table and SAC_df
    plot_df <- cbind(spes_table, slices_prevalence_df_subset)
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    
    fig_slice <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = slice)) +
      geom_point() +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
    # fig_bg_prop_A <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = bg_prop_A)) +
    #   geom_point() +
    #   ylab("AUC") +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    # 
    # fig_bg_prop_B <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = bg_prop_B)) +
    #   geom_point() +
    #   ylab("AUC") +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
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
    
    all_plots_list[[prop_cell_types$pair[i]]] <- list(slice = fig_slice + theme(legend.position="none"), 
                                                      # bg_prop_A = fig_bg_prop_A + theme(legend.position="none"), 
                                                      # bg_prop_B = fig_bg_prop_B + theme(legend.position="none"),
                                                      shape = fig_shape + theme(legend.position="none"),
                                                      variation_E = fig_variation_E + theme(legend.position="none"),
                                                      volume_E = fig_volume_E + theme(legend.position="none"),
                                                      width_N = fig_width_N + theme(legend.position="none"))
  }
  
  # Get legends
  legend_slice <-  get_legend(fig_slice + theme(legend.direction = "horizontal"))
  # legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  # legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice,
                       # legend_bg_prop_a, 
                       # legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference target pairs
  plots_pair_list <- list()
  
  for (i in seq_len(nrow(prop_cell_types))) {
    pair <- prop_cell_types$pair[i]
    
    plots <- plot_grid(all_plots_list[[pair]]$slice,
                       # all_plots_list[[pair]]$bg_prop_A,
                       # all_plots_list[[pair]]$bg_prop_B,
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

plot_entropy_prevalence_AUC <- function(spes_table, prevalence_df, slices_prevalence_df, arrangement_colname) {
  
  # Constants
  thresholds <- seq(0.01, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))
  
  # Get AUC for each prevalence gradient
  prevalence_df$AUC <- apply(prevalence_df[ , threshold_colnames], 1, sum) * 0.01
  prevalence_df <- prevalence_df[ , c("spe", "cell_types", "AUC")]
  
  # Get AUC for each prevalence gradient (slices)
  slices_prevalence_df$AUC <- apply(slices_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
  slices_prevalence_df <- slices_prevalence_df[ , c("spe", "slice", "cell_types", "AUC")]
  
  # Duplicate spes_table 5 times (5 slices)
  spes_table <- spes_table %>%
    mutate(row_num = row_number())
  spes_table <- do.call(bind_rows, replicate(5, spes_table, simplify = FALSE)) %>%
    arrange(row_num)
  spes_table$row_num <- NULL
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    
    # Subset prevalence_df for chosen pair
    prevalence_df_subset <- prevalence_df[prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Subset slices_prevalence_df for chosen pair
    slices_prevalence_df_subset <- slices_prevalence_df[slices_prevalence_df$cell_types == entropy_cell_types$cell_types[i], ]
    
    # Get difference between prevalence values in 3D and 2D slices.
    joint_df <- full_join(slices_prevalence_df_subset, prevalence_df_subset, "spe", suffix = c("_2D", "_3D"))
    
    slices_prevalence_df_subset$AUC <- (joint_df$AUC_2D - joint_df$AUC_3D) / joint_df$AUC_3D
    
    # Combine spes_table and SAC_df
    plot_df <- cbind(spes_table, slices_prevalence_df_subset)
    
    # Slight changes
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    plot_df$slice <- as.character(plot_df$slice)
    
    fig_slice <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = slice)) +
      geom_point() +
      theme_bw() +
      scale_color_manual(values = viridis::viridis(5))
    
    # fig_bg_prop_A <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = bg_prop_A)) +
    #   geom_point() +
    #   ylab("AUC") +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    # 
    # fig_bg_prop_B <- ggplot(plot_df, aes(!!sym(arrangement_colname), AUC, col = bg_prop_B)) +
    #   geom_point() +
    #   ylab("AUC") +
    #   theme_bw() +
    #   scale_color_continuous(breaks = c(0.0, 0.05, 0.1))
    
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
    
    all_plots_list[[entropy_cell_types$cell_types[i]]] <- list(slice = fig_slice + theme(legend.position="none"),
                                                               # bg_prop_A = fig_bg_prop_A + theme(legend.position="none"), 
                                                               # bg_prop_B = fig_bg_prop_B + theme(legend.position="none"),
                                                               shape = fig_shape + theme(legend.position="none"),
                                                               variation_E = fig_variation_E + theme(legend.position="none"),
                                                               volume_E = fig_volume_E + theme(legend.position="none"),
                                                               width_N = fig_width_N + theme(legend.position="none"))
  }
  
  # Get legends
  legend_slice <- get_legend(fig_slice + theme(legend.direction = "horizontal"))
  # legend_bg_prop_a <- get_legend(fig_bg_prop_A + theme(legend.direction = "horizontal"))
  # legend_bg_prop_B <- get_legend(fig_bg_prop_B + theme(legend.direction = "horizontal"))
  legend_shape <- get_legend(fig_shape + theme(legend.direction = "horizontal"))
  legend_variation_E <- get_legend(fig_variation_E + theme(legend.direction = "horizontal"))
  legend_volume_E <- get_legend(fig_volume_E + theme(legend.direction = "horizontal"))
  legend_width_N <- get_legend(fig_width_N + theme(legend.direction = "horizontal"))
  
  legends <- plot_grid(legend_slice, 
                       # legend_bg_prop_a, 
                       # legend_bg_prop_B,
                       legend_shape,
                       legend_variation_E,
                       legend_volume_E,
                       legend_width_N,
                       nrow = 1)
  
  # Combine the plots together by reference target pairs
  plots_cell_types_list <- list()
  
  for (i in seq_len(nrow(entropy_cell_types))) {
    cell_types <- entropy_cell_types$cell_types[i]
    
    plots <- plot_grid(all_plots_list[[cell_types]]$slice,
                       # all_plots_list[[cell_types]]$bg_prop_A,
                       # all_plots_list[[cell_types]]$bg_prop_B,
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

# Read mixed_slices_AMD_df
setwd("~/Objects/unsupervised/mixed_spes/analysis_2D")
mixed_slices_AMD_df <- read.table("mixed_slices_AMD_df.csv")
mixed_slices_AMD_df <- mixed_slices_AMD_df[mixed_slices_AMD_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]

mixed_AMD_plot <- plot_AMD_metric(mixed_spes_table, mixed_AMD_df, mixed_slices_AMD_df, "cluster_prop_B")

setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(mixed_AMD_plot, "mixed_AMD_plot_slicing_no_background_noise.RDS")




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

# Read mixed MS, NMS, ACINP, AE dfs (slices)
setwd("~/Objects/unsupervised/mixed_spes/analysis_2D")
mixed_slices_MS_df <- read.table("mixed_slices_MS_df.csv")
mixed_slices_NMS_df <- read.table("mixed_slices_NMS_df.csv")
mixed_slices_ACINP_df <- read.table("mixed_slices_ACINP_df.csv")
mixed_slices_AE_df <- read.table("mixed_slices_AE_df.csv")

mixed_slices_MS_df <- mixed_slices_MS_df[mixed_slices_MS_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_slices_NMS_df <- mixed_slices_NMS_df[mixed_slices_NMS_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_slices_ACINP_df <- mixed_slices_ACINP_df[mixed_slices_ACINP_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_slices_AE_df <- mixed_slices_AE_df[mixed_slices_AE_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]

# Lines
mixed_MS_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_MS_df, mixed_slices_MS_df, "MS", "cluster_prop_B")
mixed_NMS_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_NMS_df, mixed_slices_NMS_df, "NMS", "cluster_prop_B")
mixed_ACINP_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_ACINP_df, mixed_slices_ACINP_df, "ACINP", "cluster_prop_B")
mixed_AE_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_AE_df, mixed_slices_AE_df, "AE", "cluster_prop_B")

# Boxplots
mixed_MS_box_plot <- plot_gradient_metrics_type1_boxplot(mixed_spes_table, mixed_MS_df, mixed_slices_MS_df, "MS", "cluster_prop_B")
mixed_NMS_box_plot <- plot_gradient_metrics_type1_boxplot(mixed_spes_table, mixed_NMS_df, mixed_slices_NMS_df, "NMS", "cluster_prop_B")
mixed_ACINP_box_plot <- plot_gradient_metrics_type1_boxplot(mixed_spes_table, mixed_ACINP_df, mixed_slices_ACINP_df, "ACINP", "cluster_prop_B")
mixed_AE_box_plot <- plot_gradient_metrics_type1_boxplot(mixed_spes_table, mixed_AE_df, mixed_slices_AE_df, "AE", "cluster_prop_B")


setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(mixed_MS_plot, "mixed_MS_plot_slicing_no_background_noise.RDS")
saveRDS(mixed_NMS_plot, "mixed_NMS_plot_slicing_no_background_noise.RDS")
saveRDS(mixed_ACINP_plot, "mixed_ACINP_plot_slicing_no_background_noise.RDS")
saveRDS(mixed_AE_plot, "mixed_AE_plot_slicing_no_background_noise.RDS")

saveRDS(mixed_MS_box_plot, "mixed_MS_box_plot_slicing_no_background_noise.RDS")
saveRDS(mixed_NMS_box_plot, "mixed_NMS_box_plot_slicing_no_background_noise.RDS")
saveRDS(mixed_ACINP_box_plot, "mixed_ACINP_box_plot_slicing_no_background_noise.RDS")
saveRDS(mixed_AE_box_plot, "mixed_AE_box_plot_slicing_no_background_noise.RDS")


# Read mixed ACIN, CKR
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_ACIN_df <- read.table("mixed_ACIN_df.csv")
mixed_CKR_df <- read.table("mixed_CKR_df.csv")

mixed_ACIN_df <- mixed_ACIN_df[mixed_ACIN_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_CKR_df <- mixed_CKR_df[mixed_CKR_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]

# Read mixed ACIN, CKR (slices)
setwd("~/Objects/unsupervised/mixed_spes/analysis_2D")
mixed_slices_ACIN_df <- read.table("mixed_slices_ACIN_df.csv")
mixed_slices_CKR_df <- read.table("mixed_slices_CKR_df.csv")

mixed_slices_ACIN_df <- mixed_slices_ACIN_df[mixed_slices_ACIN_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_slices_CKR_df <- mixed_slices_CKR_df[mixed_slices_CKR_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]

# Lines
mixed_ACIN_plot <- plot_gradient_metrics_type2(mixed_spes_table, mixed_ACIN_df, mixed_slices_ACIN_df, "ACIN", "cluster_prop_B", 20, 100)
mixed_CKR_plot <- plot_gradient_metrics_type2(mixed_spes_table, mixed_CKR_df, mixed_slices_CKR_df, "CKR", "cluster_prop_B", 20, 100)

# Box plots
mixed_ACIN_box_plot <- plot_gradient_metrics_type2_boxplot(mixed_spes_table, mixed_ACIN_df, mixed_slices_ACIN_df, "ACIN", "cluster_prop_B", 20, 100)
mixed_CKR_box_plot <- plot_gradient_metrics_type2_boxplot(mixed_spes_table, mixed_CKR_df, mixed_slices_CKR_df, "CKR", "cluster_prop_B", 20, 100)

setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(mixed_ACIN_plot, "mixed_ACIN_plot_slicing_no_background_noise.RDS")
saveRDS(mixed_CKR_plot, "mixed_CKR_plot_slicing_no_background_noise.RDS")

saveRDS(mixed_ACIN_box_plot, "mixed_ACIN_box_plot_slicing_no_background_noise.RDS")
saveRDS(mixed_CKR_box_plot, "mixed_CKR_box_plot_slicing_no_background_noise.RDS")


# Read mixed_SAC_df
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_prop_SAC_df <- read.table("mixed_prop_SAC_df.csv")
mixed_entropy_SAC_df <- read.table("mixed_entropy_SAC_df.csv")

mixed_prop_SAC_df <- mixed_prop_SAC_df[mixed_prop_SAC_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_entropy_SAC_df <- mixed_entropy_SAC_df[mixed_entropy_SAC_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]

# Read mixed_SAC_df (slices)
setwd("~/Objects/unsupervised/mixed_spes/analysis_2D")
mixed_slices_prop_SAC_df <- read.table("mixed_slices_prop_SAC_df.csv")
mixed_slices_entropy_SAC_df <- read.table("mixed_slices_entropy_SAC_df.csv")

mixed_slices_prop_SAC_df <- mixed_slices_prop_SAC_df[mixed_slices_prop_SAC_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_slices_entropy_SAC_df <- mixed_slices_entropy_SAC_df[mixed_slices_entropy_SAC_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]

mixed_prop_SAC_plot <- plot_proportion_SAC(mixed_spes_table, mixed_prop_SAC_df, mixed_slices_prop_SAC_df, "cluster_prop_B")
mixed_entropy_SAC_plot <- plot_entropy_SAC(mixed_spes_table, mixed_entropy_SAC_df, mixed_slices_entropy_SAC_df, "cluster_prop_B")

setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(mixed_prop_SAC_plot, "mixed_prop_SAC_plot_slicing_no_background_noise.RDS")
saveRDS(mixed_entropy_SAC_plot, "mixed_entropy_SAC_plot_slicing_no_background_noise.RDS")


# Read mixed prevalence dfs
setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
mixed_prop_prevalence_df <- read.table("mixed_prop_prevalence_df.csv")
mixed_entropy_prevalence_df <- read.table("mixed_entropy_prevalence_df.csv")

mixed_prop_prevalence_df <- mixed_prop_prevalence_df[mixed_prop_prevalence_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_entropy_prevalence_df <- mixed_entropy_prevalence_df[mixed_entropy_prevalence_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]

# Read mixed prevalence dfs (slices)
setwd("~/Objects/unsupervised/mixed_spes/analysis_2D")
mixed_slices_prop_prevalence_df <- read.table("mixed_slices_prop_prevalence_df.csv")
mixed_slices_entropy_prevalence_df <- read.table("mixed_slices_entropy_prevalence_df.csv")

mixed_slices_prop_prevalence_df <- mixed_slices_prop_prevalence_df[mixed_slices_prop_prevalence_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
mixed_slices_entropy_prevalence_df <- mixed_slices_entropy_prevalence_df[mixed_slices_entropy_prevalence_df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]

# Lines
mixed_prop_prevalence_plot <- plot_proportion_prevalence(mixed_spes_table, mixed_prop_prevalence_df, mixed_slices_prop_prevalence_df, "cluster_prop_B")
mixed_entropy_prevalence_plot <- plot_entropy_prevalence(mixed_spes_table, mixed_entropy_prevalence_df, mixed_slices_entropy_prevalence_df, "cluster_prop_B")

# Boxplots
mixed_prop_prevalence_box_plot <- plot_proportion_prevalence_boxplot(mixed_spes_table, mixed_prop_prevalence_df, mixed_slices_prop_prevalence_df, "cluster_prop_B")
mixed_entropy_prevalence_box_plot <- plot_entropy_prevalence_boxplot(mixed_spes_table, mixed_entropy_prevalence_df, mixed_slices_entropy_prevalence_df, "cluster_prop_B")

# AUC
mixed_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(mixed_spes_table, mixed_prop_prevalence_df, mixed_slices_prop_prevalence_df, "cluster_prop_B")
mixed_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(mixed_spes_table, mixed_entropy_prevalence_df, mixed_slices_entropy_prevalence_df, "cluster_prop_B")

setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(mixed_prop_prevalence_plot, "mixed_prop_prevalence_plot_slicing_no_background_noise.RDS")
saveRDS(mixed_entropy_prevalence_plot, "mixed_entropy_prevalence_plot_slicing_no_background_noise.RDS")

saveRDS(mixed_prop_prevalence_box_plot, "mixed_prop_prevalence_box_plot_slicing_no_background_noise.RDS")
saveRDS(mixed_entropy_prevalence_box_plot, "mixed_entropy_prevalence_box_plot_slicing_no_background_noise.RDS")

saveRDS(mixed_prop_prevalence_AUC_plot, "mixed_prop_prevalence_AUC_plot_slicing_no_background_noise.RDS")
saveRDS(mixed_entropy_prevalence_AUC_plot, "mixed_entropy_prevalence_AUC_plot_slicing_no_background_noise.RDS")



# Read ringed_AMD_df
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_AMD_df <- read.table("ringed_AMD_df.csv")
ringed_AMD_df <- ringed_AMD_df[ringed_AMD_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]

# Read ringed_slices_AMD_df
setwd("~/Objects/unsupervised/ringed_spes/analysis_2D")
ringed_slices_AMD_df <- read.table("ringed_slices_AMD_df.csv")
ringed_slices_AMD_df <- ringed_slices_AMD_df[ringed_slices_AMD_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]

ringed_AMD_plot <- plot_AMD_metric(ringed_spes_table, ringed_AMD_df, ringed_slices_AMD_df, "width_ring_factor")

setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(ringed_AMD_plot, "ringed_AMD_plot_slicing_no_background_noise.RDS")



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

# Read ringed MS, NMS, ACINP, AE dfs (slices)
setwd("~/Objects/unsupervised/ringed_spes/analysis_2D")
ringed_slices_MS_df <- read.table("ringed_slices_MS_df.csv")
ringed_slices_NMS_df <- read.table("ringed_slices_NMS_df.csv")
ringed_slices_ACINP_df <- read.table("ringed_slices_ACINP_df.csv")
ringed_slices_AE_df <- read.table("ringed_slices_AE_df.csv")

ringed_slices_MS_df <- ringed_slices_MS_df[ringed_slices_MS_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_slices_NMS_df <- ringed_slices_NMS_df[ringed_slices_NMS_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_slices_ACINP_df <- ringed_slices_ACINP_df[ringed_slices_ACINP_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_slices_AE_df <- ringed_slices_AE_df[ringed_slices_AE_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]

# Lines
ringed_MS_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_MS_df, ringed_slices_MS_df, "MS", "width_ring_factor")
ringed_NMS_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_NMS_df, ringed_slices_NMS_df, "NMS", "width_ring_factor")
ringed_ACINP_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_ACINP_df, ringed_slices_ACINP_df, "ACINP", "width_ring_factor")
ringed_AE_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_AE_df, ringed_slices_AE_df, "AE", "width_ring_factor")

# Boxplots
ringed_MS_box_plot <- plot_gradient_metrics_type1_boxplot(ringed_spes_table, ringed_MS_df, ringed_slices_MS_df, "MS", "width_ring_factor")
ringed_NMS_box_plot <- plot_gradient_metrics_type1_boxplot(ringed_spes_table, ringed_NMS_df, ringed_slices_NMS_df, "NMS", "width_ring_factor")
ringed_ACINP_box_plot <- plot_gradient_metrics_type1_boxplot(ringed_spes_table, ringed_ACINP_df, ringed_slices_ACINP_df, "ACINP", "width_ring_factor")
ringed_AE_box_plot <- plot_gradient_metrics_type1_boxplot(ringed_spes_table, ringed_AE_df, ringed_slices_AE_df, "AE", "width_ring_factor")

setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(ringed_MS_plot, "ringed_MS_plot_slicing_no_background_noise.RDS")
saveRDS(ringed_NMS_plot, "ringed_NMS_plot_slicing_no_background_noise.RDS")
saveRDS(ringed_ACINP_plot, "ringed_ACINP_plot_slicing_no_background_noise.RDS")
saveRDS(ringed_AE_plot, "ringed_AE_plot_slicing_no_background_noise.RDS")

saveRDS(ringed_MS_box_plot, "ringed_MS_box_plot_slicing_no_background_noise.RDS")
saveRDS(ringed_NMS_box_plot, "ringed_NMS_box_plot_slicing_no_background_noise.RDS")
saveRDS(ringed_ACINP_box_plot, "ringed_ACINP_box_plot_slicing_no_background_noise.RDS")
saveRDS(ringed_AE_box_plot, "ringed_AE_box_plot_slicing_no_background_noise.RDS")


# Read ringed ACIN, CKR
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_ACIN_df <- read.table("ringed_ACIN_df.csv")
ringed_CKR_df <- read.table("ringed_CKR_df.csv")

ringed_ACIN_df <- ringed_ACIN_df[ringed_ACIN_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_CKR_df <- ringed_CKR_df[ringed_CKR_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]

# Read ringed ACIN, CKR (slices)
setwd("~/Objects/unsupervised/ringed_spes/analysis_2D")
ringed_slices_ACIN_df <- read.table("ringed_slices_ACIN_df.csv")
ringed_slices_CKR_df <- read.table("ringed_slices_CKR_df.csv")

ringed_slices_ACIN_df <- ringed_slices_ACIN_df[ringed_slices_ACIN_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_slices_CKR_df <- ringed_slices_CKR_df[ringed_slices_CKR_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]

# Lines
ringed_ACIN_plot <- plot_gradient_metrics_type2(ringed_spes_table, ringed_ACIN_df, ringed_slices_ACIN_df, "ACIN", "width_ring_factor", 20, 100)
ringed_CKR_plot <- plot_gradient_metrics_type2(ringed_spes_table, ringed_CKR_df, ringed_slices_CKR_df, "CKR", "width_ring_factor", 20, 100)

# Box plots
ringed_ACIN_box_plot <- plot_gradient_metrics_type2_boxplot(ringed_spes_table, ringed_ACIN_df, ringed_slices_ACIN_df, "ACIN", "width_ring_factor", 20, 100)
ringed_CKR_box_plot <- plot_gradient_metrics_type2_boxplot(ringed_spes_table, ringed_CKR_df, ringed_slices_CKR_df, "CKR", "width_ring_factor", 20, 100)

setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(ringed_ACIN_plot, "ringed_ACIN_plot_slicing_no_background_noise.RDS")
saveRDS(ringed_CKR_plot, "ringed_CKR_plot_slicing_no_background_noise.RDS")

saveRDS(ringed_ACIN_box_plot, "ringed_ACIN_box_plot_slicing_no_background_noise.RDS")
saveRDS(ringed_CKR_box_plot, "ringed_CKR_box_plot_slicing_no_background_noise.RDS")


# Read ringed_SAC_df
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_prop_SAC_df <- read.table("ringed_prop_SAC_df.csv")
ringed_entropy_SAC_df <- read.table("ringed_entropy_SAC_df.csv")

ringed_prop_SAC_df <- ringed_prop_SAC_df[ringed_prop_SAC_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_entropy_SAC_df <- ringed_entropy_SAC_df[ringed_entropy_SAC_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]

# Read ringed_SAC_df (slices)
setwd("~/Objects/unsupervised/ringed_spes/analysis_2D")
ringed_slices_prop_SAC_df <- read.table("ringed_slices_prop_SAC_df.csv")
ringed_slices_entropy_SAC_df <- read.table("ringed_slices_entropy_SAC_df.csv")

ringed_slices_prop_SAC_df <- ringed_slices_prop_SAC_df[ringed_slices_prop_SAC_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_slices_entropy_SAC_df <- ringed_slices_entropy_SAC_df[ringed_slices_entropy_SAC_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]

ringed_prop_SAC_plot <- plot_proportion_SAC(ringed_spes_table, ringed_prop_SAC_df, ringed_slices_prop_SAC_df, "width_ring_factor")
ringed_entropy_SAC_plot <- plot_entropy_SAC(ringed_spes_table, ringed_entropy_SAC_df, ringed_slices_entropy_SAC_df, "width_ring_factor")

setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(ringed_prop_SAC_plot, "ringed_prop_SAC_plot_slicing_no_background_noise.RDS")
saveRDS(ringed_entropy_SAC_plot, "ringed_entropy_SAC_plot_slicing_no_background_noise.RDS")



# Read ringed prevalence dfs
setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
ringed_prop_prevalence_df <- read.table("ringed_prop_prevalence_df.csv")
ringed_entropy_prevalence_df <- read.table("ringed_entropy_prevalence_df.csv")

ringed_prop_prevalence_df <- ringed_prop_prevalence_df[ringed_prop_prevalence_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_entropy_prevalence_df <- ringed_entropy_prevalence_df[ringed_entropy_prevalence_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]

# Read ringed prevalence dfs (slices)
setwd("~/Objects/unsupervised/ringed_spes/analysis_2D")
ringed_slices_prop_prevalence_df <- read.table("ringed_slices_prop_prevalence_df.csv")
ringed_slices_entropy_prevalence_df <- read.table("ringed_slices_entropy_prevalence_df.csv")

ringed_slices_prop_prevalence_df <- ringed_slices_prop_prevalence_df[ringed_slices_prop_prevalence_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
ringed_slices_entropy_prevalence_df <- ringed_slices_entropy_prevalence_df[ringed_slices_entropy_prevalence_df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]

# Lines
ringed_prop_prevalence_plot <- plot_proportion_prevalence(ringed_spes_table, ringed_prop_prevalence_df, ringed_slices_prop_prevalence_df, "width_ring_factor")
ringed_entropy_prevalence_plot <- plot_entropy_prevalence(ringed_spes_table, ringed_entropy_prevalence_df, ringed_slices_entropy_prevalence_df, "width_ring_factor")

# Boxplots
ringed_prop_prevalence_box_plot <- plot_proportion_prevalence_boxplot(ringed_spes_table, ringed_prop_prevalence_df, ringed_slices_prop_prevalence_df, "width_ring_factor")
ringed_entropy_prevalence_box_plot <- plot_entropy_prevalence_boxplot(ringed_spes_table, ringed_entropy_prevalence_df, ringed_slices_entropy_prevalence_df, "width_ring_factor")

# AUC
ringed_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(ringed_spes_table, ringed_prop_prevalence_df, ringed_slices_prop_prevalence_df, "width_ring_factor")
ringed_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(ringed_spes_table, ringed_entropy_prevalence_df, ringed_slices_entropy_prevalence_df, "width_ring_factor")

setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(ringed_prop_prevalence_plot, "ringed_prop_prevalence_plot_slicing_no_background_noise.RDS")
saveRDS(ringed_entropy_prevalence_plot, "ringed_entropy_prevalence_plot_slicing_no_background_noise.RDS")
saveRDS(ringed_prop_prevalence_box_plot, "ringed_prop_prevalence_box_plot_slicing_no_background_noise.RDS")
saveRDS(ringed_entropy_prevalence_box_plot, "ringed_entropy_prevalence_box_plot_slicing_no_background_noise.RDS")
saveRDS(ringed_prop_prevalence_AUC_plot, "ringed_prop_prevalence_AUC_plot_slicing_no_background_noise.RDS")
saveRDS(ringed_entropy_prevalence_AUC_plot, "ringed_entropy_prevalence_AUC_plot_slicing_no_background_noise.RDS")



# Read separated_AMD_df
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_AMD_df <- read.table("separated_AMD_df.csv")
separated_AMD_df <- separated_AMD_df[separated_AMD_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]

# Read separated_slices_AMD_df
setwd("~/Objects/unsupervised/separated_spes/analysis_2D")
separated_slices_AMD_df <- read.table("separated_slices_AMD_df.csv")
separated_slices_AMD_df <- separated_slices_AMD_df[separated_slices_AMD_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]

# Plots
separated_A_AMD_plot <- plot_AMD_metric(separated_A_spes_table, separated_AMD_df, separated_slices_AMD_df, "distance")

separated_B_AMD_plot <- plot_AMD_metric(separated_B_spes_table, separated_AMD_df, separated_slices_AMD_df, "distance")

setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(separated_A_AMD_plot, "separated_A_AMD_plot_slicing_no_background_noise.RDS")

saveRDS(separated_B_AMD_plot, "separated_B_AMD_plot_slicing_no_background_noise.RDS")


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

# Read separated MS, NMS, ACINP, AE dfs (slices)
setwd("~/Objects/unsupervised/separated_spes/analysis_2D")
separated_slices_MS_df <- read.table("separated_slices_MS_df.csv")
separated_slices_NMS_df <- read.table("separated_slices_NMS_df.csv")
separated_slices_ACINP_df <- read.table("separated_slices_ACINP_df.csv")
separated_slices_AE_df <- read.table("separated_slices_AE_df.csv")

separated_slices_MS_df <- separated_slices_MS_df[separated_slices_MS_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_slices_NMS_df <- separated_slices_NMS_df[separated_slices_NMS_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_slices_ACINP_df <- separated_slices_ACINP_df[separated_slices_ACINP_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_slices_AE_df <- separated_slices_AE_df[separated_slices_AE_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]

# Plots
separated_A_MS_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_MS_df, separated_slices_MS_df, "MS", "distance")
separated_A_NMS_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_NMS_df, separated_slices_NMS_df, "NMS", "distance")
separated_A_ACINP_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_ACINP_df, separated_slices_ACINP_df, "ACINP", "distance")
separated_A_AE_plot <- plot_gradient_metrics_type1(separated_A_spes_table, separated_AE_df, separated_slices_AE_df, "AE", "distance")

separated_A_MS_box_plot <- plot_gradient_metrics_type1_boxplot(separated_A_spes_table, separated_MS_df, separated_slices_MS_df, "MS", "distance")
separated_A_NMS_box_plot <- plot_gradient_metrics_type1_boxplot(separated_A_spes_table, separated_NMS_df, separated_slices_NMS_df, "NMS", "distance")
separated_A_ACINP_box_plot <- plot_gradient_metrics_type1_boxplot(separated_A_spes_table, separated_ACINP_df, separated_slices_ACINP_df, "ACINP", "distance")
separated_A_AE_box_plot <- plot_gradient_metrics_type1_boxplot(separated_A_spes_table, separated_AE_df, separated_slices_AE_df, "AE", "distance")

separated_B_MS_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_MS_df, separated_slices_MS_df, "MS", "distance")
separated_B_NMS_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_NMS_df, separated_slices_NMS_df, "NMS", "distance")
separated_B_ACINP_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_ACINP_df, separated_slices_ACINP_df, "ACINP", "distance")
separated_B_AE_plot <- plot_gradient_metrics_type1(separated_B_spes_table, separated_AE_df, separated_slices_AE_df, "AE", "distance")

separated_B_MS_box_plot <- plot_gradient_metrics_type1_boxplot(separated_B_spes_table, separated_MS_df, separated_slices_MS_df, "MS", "distance")
separated_B_NMS_box_plot <- plot_gradient_metrics_type1_boxplot(separated_B_spes_table, separated_NMS_df, separated_slices_NMS_df, "NMS", "distance")
separated_B_ACINP_box_plot <- plot_gradient_metrics_type1_boxplot(separated_B_spes_table, separated_ACINP_df, separated_slices_ACINP_df, "ACINP", "distance")
separated_B_AE_box_plot <- plot_gradient_metrics_type1_boxplot(separated_B_spes_table, separated_AE_df, separated_slices_AE_df, "AE", "distance")

setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(separated_A_MS_plot, "separated_A_MS_plot_slicing_no_background_noise.RDS")
saveRDS(separated_A_NMS_plot, "separated_A_NMS_plot_slicing_no_background_noise.RDS")
saveRDS(separated_A_ACINP_plot, "separated_A_ACINP_plot_slicing_no_background_noise.RDS")
saveRDS(separated_A_AE_plot, "separated_A_AE_plot_slicing_no_background_noise.RDS")

saveRDS(separated_A_MS_box_plot, "separated_A_MS_box_plot_slicing_no_background_noise.RDS")
saveRDS(separated_A_NMS_box_plot, "separated_A_NMS_box_plot_slicing_no_background_noise.RDS")
saveRDS(separated_A_ACINP_box_plot, "separated_A_ACINP_box_plot_slicing_no_background_noise.RDS")
saveRDS(separated_A_AE_box_plot, "separated_A_AE_box_plot_slicing_no_background_noise.RDS")

saveRDS(separated_B_MS_plot, "separated_B_MS_plot_slicing_no_background_noise.RDS")
saveRDS(separated_B_NMS_plot, "separated_B_NMS_plot_slicing_no_background_noise.RDS")
saveRDS(separated_B_ACINP_plot, "separated_B_ACINP_plot_slicing_no_background_noise.RDS")
saveRDS(separated_B_AE_plot, "separated_B_AE_plot_slicing_no_background_noise.RDS")

saveRDS(separated_B_MS_box_plot, "separated_B_MS_box_plot_slicing_no_background_noise.RDS")
saveRDS(separated_B_NMS_box_plot, "separated_B_NMS_box_plot_slicing_no_background_noise.RDS")
saveRDS(separated_B_ACINP_box_plot, "separated_B_ACINP_box_plot_slicing_no_background_noise.RDS")
saveRDS(separated_B_AE_box_plot, "separated_B_AE_box_plot_slicing_no_background_noise.RDS")


# Read separated ACIN, CKR
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_ACIN_df <- read.table("separated_ACIN_df.csv")
separated_CKR_df <- read.table("separated_CKR_df.csv")

separated_ACIN_df <- separated_ACIN_df[separated_ACIN_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_CKR_df <- separated_CKR_df[separated_CKR_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]

# Read separated ACIN, CKR (slices)
setwd("~/Objects/unsupervised/separated_spes/analysis_2D")
separated_slices_ACIN_df <- read.table("separated_slices_ACIN_df.csv")
separated_slices_CKR_df <- read.table("separated_slices_CKR_df.csv")

separated_slices_ACIN_df <- separated_slices_ACIN_df[separated_slices_ACIN_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_slices_CKR_df <- separated_slices_CKR_df[separated_slices_CKR_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]

# Plots
separated_A_ACIN_plot <- plot_gradient_metrics_type2(separated_A_spes_table, separated_ACIN_df, separated_slices_ACIN_df, "ACIN", "distance", 20, 100)
separated_A_CKR_plot <- plot_gradient_metrics_type2(separated_A_spes_table, separated_CKR_df, separated_slices_CKR_df, "CKR", "distance", 20, 100)

separated_A_ACIN_box_plot <- plot_gradient_metrics_type2_boxplot(separated_A_spes_table, separated_ACIN_df, separated_slices_ACIN_df, "ACIN", "distance", 20, 100)
separated_A_CKR_box_plot <- plot_gradient_metrics_type2_boxplot(separated_A_spes_table, separated_CKR_df, separated_slices_CKR_df, "CKR", "distance", 20, 100)

separated_B_ACIN_plot <- plot_gradient_metrics_type2(separated_B_spes_table, separated_ACIN_df, separated_slices_ACIN_df, "ACIN", "distance", 20, 100)
separated_B_CKR_plot <- plot_gradient_metrics_type2(separated_B_spes_table, separated_CKR_df, separated_slices_CKR_df, "CKR", "distance", 20, 100)

separated_B_ACIN_box_plot <- plot_gradient_metrics_type2_boxplot(separated_B_spes_table, separated_ACIN_df, separated_slices_ACIN_df, "ACIN", "distance", 20, 100)
separated_B_CKR_box_plot <- plot_gradient_metrics_type2_boxplot(separated_B_spes_table, separated_CKR_df, separated_slices_CKR_df, "CKR", "distance", 20, 100)

setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(separated_A_ACIN_plot, "separated_A_ACIN_plot_slicing_no_background_noise.RDS")
saveRDS(separated_A_CKR_plot, "separated_A_CKR_plot_slicing_no_background_noise.RDS")

saveRDS(separated_A_ACIN_box_plot, "separated_A_ACIN_box_plot_slicing_no_background_noise.RDS")
saveRDS(separated_A_CKR_box_plot, "separated_A_CKR_box_plot_slicing_no_background_noise.RDS")

saveRDS(separated_B_ACIN_plot, "separated_B_ACIN_plot_slicing_no_background_noise.RDS")
saveRDS(separated_B_CKR_plot, "separated_B_CKR_plot_slicing_no_background_noise.RDS")

saveRDS(separated_B_ACIN_box_plot, "separated_B_ACIN_plot_slicing_no_background_noise.RDS")
saveRDS(separated_B_CKR_box_plot, "separated_B_CKR_plot_slicing_no_background_noise.RDS")


# Read separated_SAC_df
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_prop_SAC_df <- read.table("separated_prop_SAC_df.csv")
separated_entropy_SAC_df <- read.table("separated_entropy_SAC_df.csv")

separated_prop_SAC_df <- separated_prop_SAC_df[separated_prop_SAC_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_entropy_SAC_df <- separated_entropy_SAC_df[separated_entropy_SAC_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]

# Read separated_SAC_df (slices)
setwd("~/Objects/unsupervised/separated_spes/analysis_2D")
separated_slices_prop_SAC_df <- read.table("separated_slices_prop_SAC_df.csv")
separated_slices_entropy_SAC_df <- read.table("separated_slices_entropy_SAC_df.csv")

separated_slices_prop_SAC_df <- separated_slices_prop_SAC_df[separated_slices_prop_SAC_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_slices_entropy_SAC_df <- separated_slices_entropy_SAC_df[separated_slices_entropy_SAC_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]

# Plots
separated_A_prop_SAC_plot <- plot_proportion_SAC(separated_A_spes_table, separated_prop_SAC_df, separated_slices_prop_SAC_df, "distance")
separated_A_entropy_SAC_plot <- plot_entropy_SAC(separated_A_spes_table, separated_entropy_SAC_df, separated_slices_entropy_SAC_df, "distance")

separated_B_prop_SAC_plot <- plot_proportion_SAC(separated_B_spes_table, separated_prop_SAC_df, separated_slices_prop_SAC_df, "distance")
separated_B_entropy_SAC_plot <- plot_entropy_SAC(separated_B_spes_table, separated_entropy_SAC_df, separated_slices_entropy_SAC_df, "distance")

setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(separated_A_prop_SAC_plot, "separated_A_prop_SAC_plot_slicing_no_background_noise.RDS")
saveRDS(separated_A_entropy_SAC_plot, "separated_A_entropy_SAC_plot_slicing_no_background_noise.RDS")

saveRDS(separated_B_prop_SAC_plot, "separated_B_prop_SAC_plot_slicing_no_background_noise.RDS")
saveRDS(separated_B_entropy_SAC_plot, "separated_B_entropy_SAC_plot_slicing_no_background_noise.RDS")


# Read separated prevalence dfs
setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
separated_prop_prevalence_df <- read.table("separated_prop_prevalence_df.csv")
separated_entropy_prevalence_df <- read.table("separated_entropy_prevalence_df.csv")

separated_prop_prevalence_df <- separated_prop_prevalence_df[separated_prop_prevalence_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_entropy_prevalence_df <- separated_entropy_prevalence_df[separated_entropy_prevalence_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]

# Read separated prevalence dfs (slices)
setwd("~/Objects/unsupervised/separated_spes/analysis_2D")
separated_slices_prop_prevalence_df <- read.table("separated_slices_prop_prevalence_df.csv")
separated_slices_entropy_prevalence_df <- read.table("separated_slices_entropy_prevalence_df.csv")

separated_slices_prop_prevalence_df <- separated_slices_prop_prevalence_df[separated_slices_prop_prevalence_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
separated_slices_entropy_prevalence_df <- separated_slices_entropy_prevalence_df[separated_slices_entropy_prevalence_df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]

# Plots
separated_A_prop_prevalence_plot <- plot_proportion_prevalence(separated_A_spes_table, separated_prop_prevalence_df, separated_slices_prop_prevalence_df, "distance")
separated_A_entropy_prevalence_plot <- plot_entropy_prevalence(separated_A_spes_table, separated_entropy_prevalence_df, separated_slices_entropy_prevalence_df, "distance")

separated_A_prop_prevalence_box_plot <- plot_proportion_prevalence_boxplot(separated_A_spes_table, separated_prop_prevalence_df, separated_slices_prop_prevalence_df, "distance")
separated_A_entropy_prevalence_box_plot <- plot_entropy_prevalence_boxplot(separated_A_spes_table, separated_entropy_prevalence_df, separated_slices_entropy_prevalence_df, "distance")

separated_A_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(separated_A_spes_table, separated_prop_prevalence_df, separated_slices_prop_prevalence_df, "distance")
separated_A_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(separated_A_spes_table, separated_entropy_prevalence_df, separated_slices_entropy_prevalence_df, "distance")

separated_B_prop_prevalence_plot <- plot_proportion_prevalence(separated_B_spes_table, separated_prop_prevalence_df, separated_slices_prop_prevalence_df, "distance")
separated_B_entropy_prevalence_plot <- plot_entropy_prevalence(separated_B_spes_table, separated_entropy_prevalence_df, separated_slices_entropy_prevalence_df, "distance")

separated_B_prop_prevalence_box_plot <- plot_proportion_prevalence_boxplot(separated_B_spes_table, separated_prop_prevalence_df, separated_slices_prop_prevalence_df, "distance")
separated_B_entropy_prevalence_box_plot <- plot_entropy_prevalence_boxplot(separated_B_spes_table, separated_entropy_prevalence_df, separated_slices_entropy_prevalence_df, "distance")

separated_B_prop_prevalence_AUC_plot <- plot_proportion_prevalence_AUC(separated_B_spes_table, separated_prop_prevalence_df, separated_slices_prop_prevalence_df, "distance")
separated_B_entropy_prevalence_AUC_plot <- plot_entropy_prevalence_AUC(separated_B_spes_table, separated_entropy_prevalence_df, separated_slices_entropy_prevalence_df, "distance")

setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
saveRDS(separated_A_prop_prevalence_plot, "separated_A_prop_prevalence_plot_slicing_no_background_noise.RDS")
saveRDS(separated_A_entropy_prevalence_plot, "separated_A_entropy_prevalence_plot_slicing_no_background_noise.RDS")
saveRDS(separated_A_prop_prevalence_box_plot, "separated_A_prop_prevalence_box_plot_slicing_no_background_noise.RDS")
saveRDS(separated_A_entropy_prevalence_box_plot, "separated_A_entropy_prevalence_box_plot_slicing_no_background_noise.RDS")
saveRDS(separated_A_prop_prevalence_AUC_plot, "separated_A_prop_prevalence_AUC_plot_slicing_no_background_noise.RDS")
saveRDS(separated_A_entropy_prevalence_AUC_plot, "separated_A_entropy_prevalence_AUC_plot_slicing_no_background_noise.RDS")

saveRDS(separated_B_prop_prevalence_plot, "separated_B_prop_prevalence_plot_slicing_no_background_noise.RDS")
saveRDS(separated_B_entropy_prevalence_plot, "separated_B_entropy_prevalence_plot_slicing_no_background_noise.RDS")
saveRDS(separated_B_prop_prevalence_box_plot, "separated_B_prop_prevalence_box_plot_slicing_no_background_noise.RDS")
saveRDS(separated_B_entropy_prevalence_box_plot, "separated_B_entropy_prevalence_box_plot_slicing_no_background_noise.RDS")
saveRDS(separated_B_prop_prevalence_AUC_plot, "separated_B_prop_prevalence_AUC_plot_slicing_no_background_noise.RDS")
saveRDS(separated_B_entropy_prevalence_AUC_plot, "separated_B_entropy_prevalence_AUC_plot_slicing_no_background_noise.RDS")



### Without background noise pdf with all plots -------------------------
setwd("~/Objects/unsupervised/plots/slicing_no_background_noise")
metrics <- c("AMD", 
             "MS", "MS_box", "NMS", "NMS_box", "ACINP", "ACINP_box", "AE", "AE_box", 
             "ACIN", "ACIN_box", "CKR", "CKR_box", 
             "prop_SAC", "prop_prevalence", "prop_prevalence_box", "prop_prevalence_AUC", 
             "entropy_SAC", "entropy_prevalence", "entropy_prevalence_box","entropy_prevalence_AUC")
arrangements <- c("mixed", "ringed", "separated_A", "separated_B")

pdf("plots.pdf", width = 15, height = 10)

for (metric in metrics) {
  for (arrangement in arrangements) {
    plot_file_name <- paste(arrangement, "_", metric, "_plot_slicing_no_background_noise.RDS", sep = "")
    print(readRDS(plot_file_name))
  }
}
dev.off()
