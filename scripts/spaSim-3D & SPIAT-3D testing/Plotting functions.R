### Function for non-gradient output ------------------------------------------
plot_non_gradient_metric <- function(spes_table, 
                                     metric, 
                                     metric_cell_types, 
                                     metric_df, 
                                     arrangement, 
                                     plot_list) {
  
  ### Modify plot_list
  # Change plot_list arrangement to inputted arrangement
  plot_list$arrangement$x_aes <- arrangement
  
  # Change plot_list y_aes to inputted metric
  for (i in seq_along(plot_list)) {
    # Modify the y_aes element
    plot_list[[i]]$y_aes <- metric
  }
  

  # Define plotting function
  create_plot <- function(data, x_aes, y_aes, title = "") {
    data <- data[!is.na(data[[x_aes]]), ]
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw()
    
    if (typeof(data[[x_aes]]) == "double") {
      plot <- plot + geom_point()
    }
    # Factored character is an integer
    else if (typeof(data[[x_aes]]) == "integer") {
      plot <- plot + geom_violin()
    }
    return(plot)
  }
  
  # Add Ellipsoid variation and volume to spes_table
  radii_E_df <- spes_table[ , c("radius_x_E", "radius_y_E", "radius_z_E")]
  spes_table$volume_E <- radii_E_df$radius_x_E * radii_E_df$radius_y_E * radii_E_df$radius_z_E
  spes_table$variation_E <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
  
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
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
    
    # Generate plots based on plot_list, use final column of metric_cell_types
    plots_list[[metric_cell_types[i, ncol(metric_cell_types)]]] <- lapply(plot_list, function(plot_def) {
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
                                        nrow = length(combined_plots_list), ncol = 1)
  
  methods::show(non_gradient_metric_plot)
  
  return(non_gradient_metric_plot)
}

plot_gradient_metric <- function(spes_table, 
                                 metric, 
                                 metric_cell_types, 
                                 metric_df, 
                                 arrangement_colname, 
                                 plot_list) {
  
}


### Test for non-gradient metrics  ----------------------------------------------------------------------
# Read mixed_spes_table
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/spe_tables")
mixed_spes_table <- read.table("mixed_spes_table.csv")

### Set up plot lists
non_gradient_plot_list <- list(
  arrangement = list(x_aes = "temp", y_aes = "temp"),
  bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "temp"),
  bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "temp"),
  shape = list(x_aes = "shape", y_aes = "temp"),
  variation_E = list(x_aes = "variation_E", y_aes = "temp"),
  volume_E = list(x_aes = "volume_E", y_aes = "temp"),
  width_N = list(x_aes = "width_N", y_aes = "temp")
)





# Read mixed_AMD_df
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/mixed")
mixed_AMD_df <- read.table("mixed_AMD_df.csv")

# AMD pairs are A/A, A/B, B/A, B/B
AMD_cell_types <- data.frame(ref = c("A", "A", "B", "B"), tar = c("A", "B", "A", "B"))
AMD_cell_types$pair <- paste(AMD_cell_types$ref, AMD_cell_types$tar, sep = "/")

mixed_AMD_plot <- plot_non_gradient_metric(mixed_spes_table, 
                                           "AMD", 
                                           AMD_cell_types,
                                           mixed_AMD_df, 
                                           "cluster_prop_B", 
                                           non_gradient_plot_list)


setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/mixed")
mixed_prop_SAC_df <- read.table("mixed_prop_SAC_df.csv")
mixed_entropy_SAC_df <- read.table("mixed_entropy_SAC_df.csv")

prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))
prop_cell_types$pair <- paste(prop_cell_types$ref, prop_cell_types$tar, sep = "/")

entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))

mixed_prop_SAC_plot <- plot_non_gradient_metric(mixed_spes_table, 
                                                "prop_SAC", 
                                                prop_cell_types,
                                                mixed_prop_SAC_df, 
                                                "cluster_prop_B", 
                                                non_gradient_plot_list)

mixed_entropy_SAC_plot <- plot_non_gradient_metric(mixed_spes_table, 
                                                   "entropy_SAC", 
                                                   entropy_cell_types,
                                                   mixed_entropy_SAC_df, 
                                                   "cluster_prop_B", 
                                                   non_gradient_plot_list)

# Read mixed prevalence dfs
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/mixed")
mixed_prop_prevalence_df <- read.table("mixed_prop_prevalence_df.csv")
mixed_entropy_prevalence_df <- read.table("mixed_entropy_prevalence_df.csv")

thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")

mixed_prop_prevalence_df$prop_AUC <- apply(mixed_prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
mixed_prop_prevalence_df <- mixed_prop_prevalence_df[ , c("spe", "reference", "target", "prop_AUC")]

mixed_entropy_prevalence_df$entropy_AUC <- apply(mixed_entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
mixed_entropy_prevalence_df <- mixed_entropy_prevalence_df[ , c("spe", "cell_types", "entropy_AUC")]

mixed_prop_prevalence_AUC_plot <- plot_non_gradient_metric(mixed_spes_table, 
                                                           "prop_AUC", 
                                                           prop_cell_types,
                                                           mixed_prop_prevalence_df, 
                                                           "cluster_prop_B", 
                                                           non_gradient_plot_list)

mixed_entropy_prevaelnce_AUC_plot <- plot_non_gradient_metric(mixed_spes_table, 
                                                              "entropy_AUC", 
                                                              entropy_cell_types,
                                                              mixed_entropy_prevalence_df, 
                                                              "cluster_prop_B", 
                                                              non_gradient_plot_list)

