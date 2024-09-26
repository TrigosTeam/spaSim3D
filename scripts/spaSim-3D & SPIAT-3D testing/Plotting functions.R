library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)

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
    paste(alpha * 10, power - 1L, sep = "e")
  }
  
  create_plot <- function(data, x_aes, y_aes, title = "") {
    data <- data[!is.na(data[[x_aes]]), ]
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw()
    
    if (typeof(data[[x_aes]]) == "double") {
      plot <- plot + geom_point()
    }
    # Factored character is an integer
    else if (typeof(data[[x_aes]]) %in% c("integer", "character")) {
      plot <- plot + geom_violin()
    }
    # Use scientific notation for ellipsoid volume
    if (x_aes == "E_volume") {
      plot <- plot + scale_x_continuous(labels = formatCustomSci)
    }
    return(plot)
  }
  
  # Add Ellipsoid variation and volume to spes_table
  radii_E_df <- spes_table[ , c("E_radius_x", "E_radius_y", "E_radius_z")]
  spes_table$E_volume <- radii_E_df$E_radius_x * radii_E_df$E_radius_y * radii_E_df$E_radius_z
  spes_table$E_radii_CV <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
  
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
    if (!is.null(plot_df$slice)) {
      plot_df$slice <- as.character(plot_df$slice)
    }
    
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
                                        nrow = length(combined_plots_list), ncol = 1)
  
  methods::show(non_gradient_metric_plot)
  
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
    paste(alpha * 10, power - 1L, sep = "e")
  }
  create_plot <- function(data, x_aes, y_aes, color_aes, group_aes, title = "") {
    
    data <- data[!is.na(data[[color_aes]]), ]
    
    plot <- ggplot(data, aes_string(x = x_aes, y = y_aes, group = group_aes, color = color_aes)) +
      labs(title = title, x = x_aes, y = y_aes) +
      theme_bw() +
      geom_line()
    
    # Use scientific notation for ellipsoid volume
    if (color_aes == "E_volume") {
      plot <- plot + scale_x_continuous(labels = formatCustomSci)
    }
    
    return(plot)
  }
  
  # Add Ellipsoid variation and volume to spes_table
  radii_E_df <- spes_table[ , c("E_radius_x", "E_radius_y", "E_radius_z")]
  spes_table$E_volume <- radii_E_df$E_radius_x * radii_E_df$E_radius_y * radii_E_df$E_radius_z
  spes_table$E_radii_CV <- (apply(radii_E_df, 1, sd) / rowMeans(radii_E_df)) * 100
  
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
    plot_df$shape <- factor(plot_df$shape, c("Ellipsoid", "Network"))
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
    
    # Get final column
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
  
  methods::show(gradient_metric_with_legends_plot)
  
  return(gradient_metric_with_legends_plot)
}


### Test for non-gradient metrics  ----------------------------------------------------------------------

### Set up plot lists
non_gradient_plots_metadata <- list(
  arrangement = list(x_aes = "temp_arrangement", y_aes = "metric"),
  bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "metric"),
  bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "metric"),
  shape = list(x_aes = "shape", y_aes = "metric"),
  E_radii_CV = list(x_aes = "E_radii_CV", y_aes = "metric"),
  E_volume = list(x_aes = "E_volume", y_aes = "metric"),
  N_width = list(x_aes = "N_width", y_aes = "metric")
)

# Read mixed_spes_table
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/spe_tables")
mixed_spes_table <- read.table("mixed_spes_table.csv")





# Read mixed_AMD_df
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/mixed")
mixed_AMD_df <- read.table("mixed_AMD_df.csv")

mixed_AMD_plot <- plot_non_gradient_metric(mixed_spes_table, 
                                           "AMD", 
                                           mixed_AMD_df, 
                                           "cluster_prop_B", 
                                           non_gradient_plots_metadata)


setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/mixed")
mixed_prop_SAC_df <- read.table("mixed_prop_SAC_df.csv")
mixed_entropy_SAC_df <- read.table("mixed_entropy_SAC_df.csv")

mixed_prop_SAC_plot <- plot_non_gradient_metric(mixed_spes_table, 
                                                "prop_SAC", 
                                                mixed_prop_SAC_df, 
                                                "cluster_prop_B", 
                                                non_gradient_plots_metadata)

mixed_entropy_SAC_plot <- plot_non_gradient_metric(mixed_spes_table, 
                                                   "entropy_SAC", 
                                                   mixed_entropy_SAC_df, 
                                                   "cluster_prop_B", 
                                                   non_gradient_plots_metadata)

# Read mixed prevalence dfs
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/mixed")
mixed_prop_prevalence_df <- read.table("mixed_prop_prevalence_df.csv")
mixed_entropy_prevalence_df <- read.table("mixed_entropy_prevalence_df.csv")

thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")

mixed_prop_prevalence_df$prop_AUC <- apply(mixed_prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
mixed_prop_AUC_df <- mixed_prop_prevalence_df[ , c("spe", "reference", "target", "prop_AUC")]

mixed_entropy_prevalence_df$entropy_AUC <- apply(mixed_entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
mixed_entropy_AUC_df <- mixed_entropy_prevalence_df[ , c("spe", "cell_types", "entropy_AUC")]

mixed_prop_prevalence_AUC_plot <- plot_non_gradient_metric(mixed_spes_table, 
                                                           "prop_AUC", 
                                                           mixed_prop_AUC_df, 
                                                           "cluster_prop_B", 
                                                           non_gradient_plots_metadata)

mixed_entropy_prevaelnce_AUC_plot <- plot_non_gradient_metric(mixed_spes_table, 
                                                              "entropy_AUC", 
                                                              mixed_entropy_AUC_df, 
                                                              "cluster_prop_B", 
                                                              non_gradient_plots_metadata)



### Test for gradient metrics  ----------------------------------------------------------------------

### Set up plot lists
gradient_plots_metadata <- list(
  arrangement = list(x_aes = "gradient", y_aes = "metric", color_aes = "temp_arrangement"),
  bg_prop_A = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_A"),
  bg_prop_B = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_B"),
  shape = list(x_aes = "gradient", y_aes = "metric", color_aes = "shape"),
  E_radii_CV = list(x_aes = "gradient", y_aes = "metric", color_aes = "E_radii_CV"),
  E_volume = list(x_aes = "gradient", y_aes = "metric", color_aes = "E_volume"),
  N_width = list(x_aes = "gradient", y_aes = "metric", color_aes = "N_width")
)

# Read mixed_spes_table
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/spe_tables")
mixed_spes_table <- read.table("mixed_spes_table.csv")

setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/mixed")
mixed_MS_df <- read.table("mixed_MS_df.csv")

mixed_MS_plot <- plot_gradient_metric(mixed_spes_table, 
                                      "MS",
                                      mixed_MS_df, 
                                      "cluster_prop_B", 
                                      "radius",
                                      gradient_plots_metadata)


# Read mixed ACIN, CKR
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/mixed")
mixed_ACIN_df <- read.table("mixed_ACIN_df.csv")
mixed_CKR_df <- read.table("mixed_CKR_df.csv")

mixed_ACIN_plot <- plot_gradient_metric(mixed_spes_table, 
                                        "ACIN",
                                        mixed_ACIN_df, 
                                        "cluster_prop_B", 
                                        "radius",
                                        gradient_plots_metadata)

# Read mixed prevalence dfs
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/mixed")
mixed_prop_prevalence_df <- read.table("mixed_prop_prevalence_df.csv")
mixed_entropy_prevalence_df <- read.table("mixed_entropy_prevalence_df.csv")


mixed_prop_prevalence_plot <- plot_gradient_metric(mixed_spes_table, 
                                                   "prop_prevalence",
                                                   mixed_prop_prevalence_df, 
                                                   "cluster_prop_B",
                                                   "threshold",
                                                   gradient_plots_metadata)

mixed_entropy_prevalence_plot <- plot_gradient_metric(mixed_spes_table, 
                                                      "entropy_prevalence",
                                                      mixed_entropy_prevalence_df, 
                                                      "cluster_prop_B",
                                                      "threshold",
                                                      gradient_plots_metadata)



