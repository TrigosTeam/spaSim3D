library(cowplot)

### 1.1. Function to get plot for APD, AMD -------------------------------------

### 1.2. Function to get plot for MS, NMS, ACINP, AE gradient metrics ------------

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






### 2.1. Mixed spes APD, AMD ------------------------------------------------------

# Read mixed_spes_table
setwd("~/Objects/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Read mixed_AMD_df
setwd("~/Objects/mixed_spes/analysis_3D")
mixed_AMD_df <- read.table("mixed_AMD_df.csv")


# AMD pairs are A/A, A/B, B/A, B/B
AMD_pairs <- data.frame(cell1 = c("A", "A", "B", "B"),
                        cell2 = c("A", "B", "A", "B"))
AMD_pairs$pair <- paste(AMD_pairs$cell1, AMD_pairs$cell2, sep = "/")


# Put all plots into an organised list
mixed_AMD_all_plots_list <- list()

for (i in seq(nrow(AMD_pairs))) {
  
  # Subset mixed_AMD_df for chosen pair
  plot_df <- mixed_AMD_df[mixed_AMD_df$reference == AMD_pairs[i, "cell1"] & mixed_AMD_df$target == AMD_pairs[i, "cell2"], ]
  
  # Combine mixed_spes_table and mixed_AMD_df
  plot_df <- cbind(mixed_spes_table, plot_df)
 
   # Create a 'key' column which groups simulations if they have the same bg_type, shape and size (but not arrangement)
  plot_df$key <- paste(plot_df$bg_type, plot_df$shape, plot_df$size, sep = "_")
  
  # Factor
  plot_df$bg_type <- factor(plot_df$bg_type, c("O", "A", "B", "AB"))
  plot_df$shape <- factor(plot_df$shape, c("Sphere", "Ellipsoid", "Network"))
  plot_df$size <- factor(plot_df$size, c("Small", "Medium", "Large"))
  plot_df$arrangement <- factor(plot_df$arrangement, c("M1", "M2", "M3"))
  
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
  
  mixed_AMD_all_plots_list[[AMD_pairs[i, "pair"]]] <- list(bg_type = fig_bg_type, shape = fig_shape, size = fig_size)
}


# Combine the plots together by pairs
library(cowplot)
mixed_AMD_plots_pair_list <- list()

for (i in seq(nrow(AMD_pairs))) {
  pair <- AMD_pairs[i, "pair"]
  
  plots <- plot_grid(mixed_AMD_all_plots_list[[pair]]$bg_type, mixed_AMD_all_plots_list[[pair]]$shape, mixed_AMD_all_plots_list[[pair]]$size, nrow = 1, ncol = 3)
  
  title <- ggdraw() + 
    draw_label(paste("Reference:", AMD_pairs[i, "cell1"], "Target:", AMD_pairs[i, "cell2"]), 
               fontface='bold')
  
  fig <- plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1))
  
  mixed_AMD_plots_pair_list[[pair]] <- fig
}

# Combine the combined plots into one big plot
mixed_AMD_plot <- plot_grid(mixed_AMD_plots_pair_list$`A/A`, 
                            mixed_AMD_plots_pair_list$`A/B`, 
                            mixed_AMD_plots_pair_list$`B/A`, 
                            mixed_AMD_plots_pair_list$`B/B`, 
                            nrow = 2, ncol = 2, scale = 0.9)

setwd("~/Objects/mixed_spes/analysis_3D/plots")
saveRDS(mixed_AMD_plot, "mixed_AMD_plot.rds")

methods::show(mixed_AMD_plot)




### 2.2.  Mixed spes MS, NMS, ACINP, AE -----------------------------------

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

### 3.1. Ringed spes APD, AMD -------------------------------------------------

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

