library(cowplot)

### 1.1. Function to get plot for APD ----------------------------------------
### 1.2. Function to get plot for AMD -------------------------------------
plot_AMD_metric <- function(spes_table, AMD_df, arrangements) {
  
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
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shape <- factor(plot_df$shape, c("Sphere", "Ellipsoid", "Network"))
    plot_df$size <- factor(plot_df$size, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    fig_bg_type <- ggplot(plot_df, aes(arrangement, AMD, group = key, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom")
    
    fig_shape <- ggplot(plot_df, aes(arrangement, AMD, group = key, col = shape)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom")
    
    fig_size <- ggplot(plot_df, aes(arrangement, AMD, group = key, col = size)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom")
    
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


### 1.3. Function to get plot for MS, NMS, ACINP, AE gradient metrics ------------

plot_gradient_metrics_type1 <- function(spes_table, gradient_metric_df, metric, arrangements) {
  
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
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = arrangement)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric)
    
    fig_bg_type <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric)
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric)
    
    fig_size <- ggplot(plot_df, aes(variable, value, group = spe, col = size)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric)
    
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






### 1.4. Function to get plot for ACIN, CKR gradient metrics ------------------

plot_gradient_metrics_type2 <- function(spes_table, gradient_metric_df, metric, arrangements, min_radius, max_radius) {
  
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
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = arrangement)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric)
    
    fig_bg_type <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric)
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric)
    
    fig_size <- ggplot(plot_df, aes(variable, value, group = spe, col = size)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "radius", y = metric)
    
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


### 1.5. Function to get plot for SAC ----------------------------------------
plot_SAC_metric <- function(spes_table, SAC_df, arrangements) {
  
  # SAC applies to proportion and entropy
  metrics <- c("proportion", "entropy")
  
  # Combine spes_table and SAC_df
  plot_df <- cbind(spes_table, SAC_df)
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(length(metrics))) {
    
    # Create a 'key' column which groups simulations if they have the same bg_type, shape and size (but not arrangement)
    plot_df$key <- paste(plot_df$bg_type, plot_df$shape, plot_df$size, sep = "_")
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shape <- factor(plot_df$shape, c("Sphere", "Ellipsoid", "Network"))
    plot_df$size <- factor(plot_df$size, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    fig_bg_type <- ggplot(plot_df, aes(arrangement, !!sym(metrics[i]), group = key, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC")
    
    fig_shape <- ggplot(plot_df, aes(arrangement, !!sym(metrics[i]), group = key, col = shape)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC")
    
    fig_size <- ggplot(plot_df, aes(arrangement, !!sym(metrics[i]), group = key, col = size)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      ylab("SAC")
    
    all_plots_list[[metrics[i]]] <- list(bg_type = fig_bg_type, shape = fig_shape, size = fig_size)
  }
  
  
  # Combine the plots together by pairs
  library(cowplot)
  plots_metric_list <- list()
  
  for (i in seq(length(metrics))) {
    metric <- metrics[i]
    
    plots <- plot_grid(all_plots_list[[metric]]$bg_type, all_plots_list[[metric]]$shape, all_plots_list[[metric]]$size, nrow = 1, ncol = 3)
    
    title <- ggdraw() + 
      draw_label(paste("Metric:", metric), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_metric_list[[metric]] <- fig
  }
  
  # Combine the combined plots into one big plot
  SAC_plot <- plot_grid(plots_metric_list$proportion, 
                        plots_metric_list$entropy,
                        nrow = 2, ncol = 1, scale = 0.9)
  
  methods::show(SAC_plot)
  
  return(SAC_plot)
}


### 1.6. Function to get plot for prevalence ----------------------------------
plot_prevalence <- function(spes_table, prevalence_df, arrangements) {
  
  # Constants
  metrics <- c("proportion", "entropy")
  thresholds <- seq(0, 1, 0.01)
  threshold_colnames <- paste("t", thresholds, sep = "")
  
  
  # Put all plots into an organised list
  all_plots_list <- list()
  
  for (i in seq(length(metrics))) {
    # Subset prevalence_df for current metric
    plot_df <- prevalence_df[prevalence_df$metric == metrics[i], ]
    
    # Combine spes_table and mixed_AMD_df
    plot_df <- cbind(spes_table, plot_df)
    
    # Melt
    plot_df <- reshape2::melt(plot_df, , threshold_colnames)
    
    # Extract threshold value from threshold strings (t0 -> 0, t0.01 -> 0.01...)
    plot_df$variable <- unfactor(plot_df$variable)
    plot_df$variable <- as.numeric(substr(plot_df$variable, 2, nchar(plot_df$variable)))
    
    # Factor
    plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
    plot_df$shape <- factor(plot_df$shape, c("Sphere", "Ellipsoid", "Network"))
    plot_df$size <- factor(plot_df$size, c("Small", "Medium", "Large"))
    plot_df$arrangement <- factor(plot_df$arrangement, arrangements)
    
    fig_arrangement <- ggplot(plot_df, aes(variable, value, group = spe, col = arrangement)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence")
    
    fig_bg_type <- ggplot(plot_df, aes(variable, value, group = spe, col = bg_type)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence")
    
    fig_shape <- ggplot(plot_df, aes(variable, value, group = spe, col = shape)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence")
    
    fig_size <- ggplot(plot_df, aes(variable, value, group = spe, col = size)) +
      geom_line() +
      theme_bw() +
      theme(legend.position="bottom") +
      labs(x = "threshold", y = "prevalence")
    
    all_plots_list[[metrics[i]]] <- list(arrangement = fig_arrangement, bg_type = fig_bg_type, shape = fig_shape, size = fig_size)
    
  }
  
  
  # Combine the plots together by metric
  plots_metric_list <- list()
  
  for (i in seq(length(metrics))) {
    metric <- metrics[i]
    
    plots <- plot_grid(all_plots_list[[metric]]$arrangement,
                       all_plots_list[[metric]]$bg_type, 
                       all_plots_list[[metric]]$shape, 
                       all_plots_list[[metric]]$size, nrow = 1, ncol = 4)
    
    title <- ggdraw() + 
      draw_label(paste("Metric:", metric), 
                 fontface='bold')
    
    fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
    
    plots_metric_list[[metric]] <- fig
  }
  
  # Combine the combined plots into one big plot
  combined_plot <- plot_grid(plots_metric_list$proportion,
                             plots_metric_list$entropy, 
                             nrow = 2, ncol = 1, scale = 0.9)
  
  methods::show(combined_plot)
  
  return(combined_plot)
}



### 2.1. Mixed spes APD ------------------------------------------------------
### 2.2. Mixed spes AMD ------------------------------------------------------

# Read mixed_spes_table
setwd("~/Objects/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Read mixed_AMD_df
setwd("~/Objects/mixed_spes/analysis_3D")
mixed_AMD_df <- read.table("mixed_AMD_df.csv")

mixed_AMD_plot <- plot_AMD_metric(mixed_spes_table, mixed_AMD_df, c("M1", "M2", "M3"))

setwd("~/Objects/mixed_spes/analysis_3D/plots")
saveRDS(mixed_AMD_plot, "mixed_AMD_plot.rds")




### 2.3. Mixed spes MS, NMS, ACINP, AE -----------------------------------

# Read mixed_spes_table
setwd("~/Objects/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Read mixed MS, NMS, ACINP, AE dfs
setwd("~/Objects/mixed_spes/analysis_3D")
mixed_MS_df <- read.table("mixed_MS_df.csv")
mixed_NMS_df <- read.table("mixed_NMS_df.csv")
mixed_ACINP_df <- read.table("mixed_ACINP_df.csv")
mixed_AE_df <- read.table("mixed_AE_df.csv")

mixed_MS_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_MS_df, "MS", c("M1", "M2", "M3"))
mixed_NMS_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_NMS_df, "NMS", c("M1", "M2", "M3"))
mixed_ACINP_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_ACINP_df, "ACINP", c("M1", "M2", "M3"))
mixed_AE_plot <- plot_gradient_metrics_type1(mixed_spes_table, mixed_AE_df, "AE", c("M1", "M2", "M3"))

setwd("~/Objects/mixed_spes/analysis_3D/plots")
# saveRDS()

### 2.4. Mixed spes ACIN, CKR ------------------------------------------------

# Read mixed_spes_table
setwd("~/Objects/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Read mixed ACIN, CKR
setwd("~/Objects/mixed_spes/analysis_3D")
mixed_ACIN_df <- read.table("mixed_ACIN_df.csv")
mixed_CKR_df <- read.table("mixed_CKR_df.csv")

# Get plots
mixed_ACIN_plot <- plot_gradient_metrics_type2(mixed_spes_table, mixed_ACIN_df, "ACIN", c("M1", "M2", "M3"), 0, 50)

mixed_CKR_plot <- plot_gradient_metrics_type2(mixed_spes_table, mixed_CKR_df, "CKR", c("M1", "M2", "M3"), 15, 50)


### 2.5. Mixed spes SAC ------------------------------------------------------

# Read mixed_spes_table
setwd("~/Objects/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Read mixed_SAC_df
setwd("~/Objects/mixed_spes/analysis_3D")
mixed_SAC_df <- read.table("mixed_SAC_df.csv")

mixed_SAC_plot <- plot_SAC_metric(mixed_spes_table, mixed_SAC_df, c("M1", "M2", "M3"))

setwd("~/Objects/mixed_spes/analysis_3D/plots")
saveRDS(mixed_SAC_plot, "mixed_SAC_plot.rds")

### 2.6. Mixed spes prevalence ------------------------------------------------
# Read mixed_spes_table
setwd("~/Objects/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Read mixed_AMD_df
setwd("~/Objects/mixed_spes/analysis_3D")
mixed_prevalence_df <- read.table("mixed_prevalence_df.csv")

mixed_prevalence_plot <- plot_prevalence(mixed_spes_table, mixed_prevalence_df, c("M1", "M2", "M3"))

setwd("~/Objects/mixed_spes/analysis_3D/plots")


### 3.1. Ringed spes AMD -------------------------------------------------

# Read ringed_spes_table
setwd("~/Objects/spes_table")
ringed_spes_table <- read.table("ringed_spes_table.csv")

# Read ringed_AMD_df
setwd("~/Objects/ringed_spes/analysis_3D")
ringed_AMD_df <- read.table("ringed_AMD_df.csv")

ringed_AMD_plot <- plot_AMD_metric(ringed_spes_table, ringed_AMD_df, c("R1", "R2", "R3"))

setwd("~/Objects/ringed_spes/analysis_3D/plots")
saveRDS(ringed_AMD_plot, "ringed_AMD_plot.rds")


### 3.2. Ringed spes MS, NMS, ACINP, AE ---------------------------------------
# Read ringed_spes_table
setwd("~/Objects/spes_table")
ringed_spes_table <- read.table("ringed_spes_table.csv")

# Read ringed MS, NMS, ACINP, AE dfs
setwd("~/Objects/ringed_spes/analysis_3D")
ringed_MS_df <- read.table("ringed_MS_df.csv")
ringed_NMS_df <- read.table("ringed_NMS_df.csv")
ringed_ACINP_df <- read.table("ringed_ACINP_df.csv")
ringed_AE_df <- read.table("ringed_AE_df.csv")

ringed_MS_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_MS_df, "MS", c("R1", "R2", "R3"))
ringed_NMS_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_NMS_df, "NMS", c("R1", "R2", "R3"))
ringed_ACINP_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_ACINP_df, "ACINP", c("R1", "R2", "R3"))
ringed_AE_plot <- plot_gradient_metrics_type1(ringed_spes_table, ringed_AE_df, "AE", c("R1", "R2", "R3"))

setwd("~/Objects/ringed_spes/analysis_3D/plots")
# saveRDS()

