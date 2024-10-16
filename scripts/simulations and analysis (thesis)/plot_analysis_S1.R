library(cowplot)
library(ggplot2)
library(S4Vectors)
library(stringr)
library(dplyr)

### Utility functions -------------------------
get_gradient <- function(metric) {
  if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR")) return("radius")
  return("threshold")
}

### Read metric_df_lists --------------------------------------------------------------
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
metric_df_lists3D <- readRDS("metric_df_lists3D.RDS")
metric_df_lists2D <- readRDS("metric_df_lists2D.RDS")

### Turn gradient radii metrics into AUC and add to metric_df list ------------------
get_AUC_for_radii_gradient_metrics <- function(y) {
  x <- radii
  h <- diff(x)[1]
  n <- length(x)
  
  AUC <- (h / 2) * (y[1] + 2 * sum(y[2:(n - 1)]) + y[n])
  
  return(AUC)
}

arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

radii <- seq(20, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

gradient_radii_metrics <- c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR")


for (arrangement in arrangements) {
  
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (metric in gradient_radii_metrics) {
      metric_AUC_name <- paste(metric, "AUC", sep = "_")
      
      if (metric %in% c("MS", "NMS", "ACIN", "CKR")) {
        subset_colnames <- c("spe", "reference", "target", metric_AUC_name)
      }
      else {
        subset_colnames <- c("spe", "reference", metric_AUC_name)
      }
      
      # 3D
      df <- metric_df_lists3D[[spes_metadata_index]][[metric]]
      df[[metric_AUC_name]] <- apply(df[ , radii_colnames], 1, get_AUC_for_radii_gradient_metrics)
      
      df <- df[ , subset_colnames]
      metric_df_lists3D[[spes_metadata_index]][[metric_AUC_name]] <- df
      
      # 2D
      subset_colnames <- c(subset_colnames, "slice")
      df <- metric_df_lists2D[[spes_metadata_index]][[metric]]
      df[[metric_AUC_name]] <- apply(df[ , radii_colnames], 1, get_AUC_for_radii_gradient_metrics)
      
      df <- df[ , subset_colnames]
      metric_df_lists2D[[spes_metadata_index]][[metric_AUC_name]] <- df
      
    }
  }
}


### Turn threshold radii metrics into AUC and add to metric_df list--------------
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")

for (arrangement in arrangements) {
  
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    # prop_AUC 3D
    prop_prev_df <- metric_df_lists3D[[spes_metadata_index]][["prop_prev"]]
    prop_prev_df$prop_AUC <- apply(prop_prev_df[ , threshold_colnames], 1, sum) * 0.01
    prop_AUC_df <- prop_prev_df[ , c("spe", "reference", "target", "prop_AUC")]
    metric_df_lists3D[[spes_metadata_index]][["prop_AUC"]] <- prop_AUC_df
    
    # entropy_AUC 3D
    entropy_prev_df <- metric_df_lists3D[[spes_metadata_index]][["entropy_prev"]]
    entropy_prev_df$entropy_AUC <- apply(entropy_prev_df[ , threshold_colnames], 1, sum) * 0.01
    entropy_AUC_df <- entropy_prev_df[ , c("spe", "cell_types", "entropy_AUC")]
    metric_df_lists3D[[spes_metadata_index]][["entropy_AUC"]] <- entropy_AUC_df
    
    # prop_AUC 2D
    prop_prev_df <- metric_df_lists2D[[spes_metadata_index]][["prop_prev"]]
    prop_prev_df$prop_AUC <- apply(prop_prev_df[ , threshold_colnames], 1, sum) * 0.01
    prop_AUC_df <- prop_prev_df[ , c("spe", "slice", "reference", "target", "prop_AUC")]
    metric_df_lists2D[[spes_metadata_index]][["prop_AUC"]] <- prop_AUC_df
    
    # entropy_AUC 2D
    entropy_prev_df <- metric_df_lists2D[[spes_metadata_index]][["entropy_prev"]]
    entropy_prev_df$entropy_AUC <- apply(entropy_prev_df[ , threshold_colnames], 1, sum) * 0.01
    entropy_AUC_df <- entropy_prev_df[ , c("spe", "slice", "cell_types", "entropy_AUC")]
    metric_df_lists2D[[spes_metadata_index]][["entropy_AUC"]] <- entropy_AUC_df
  }
}


### Get plots for 3D metric analysis (non-gradient) ------------------------------------------
# Read spes_table
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B"
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 

# Set up plots metadata
non_gradient_plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "temp_arrangement", y_aes = "metric"),
    bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "metric"),
    bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "metric"),
    E_volume = list(x_aes = "E_volume", y_aes = "metric")
  ),
  network = list(
    arrangement = list(x_aes = "temp_arrangement", y_aes = "metric"),
    bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "metric"),
    bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "metric"),
    N_width = list(x_aes = "N_width", y_aes = "metric")
  )
)

# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))


metric_plots3D_non_gradient <- list(mixed_ellipsoid = list(),
                                    mixed_network = list(),
                                    ringed_ellipsoid = list(),
                                    ringed_network = list(),
                                    separated_ellipsoid = list(),
                                    separated_network = list())
# list for plot labels
plot_labels <- list()
i <- 1
for (shape in shapes) {
  for (arrangement in arrangements) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    plot_labels[[spes_metadata_index]] <- LETTERS[i:(i + 3)]
    i <- i + 4
  }
}

for (shape in shapes) {
  for (arrangement in arrangements) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      metric_plots3D_non_gradient[[spes_metadata_index]][[metric]] <- plot_non_gradient_metric(spes_table_subset, 
                                                                                               metric, 
                                                                                               metric_df_lists3D[[spes_metadata_index]][[metric]], 
                                                                                               arrangement_parameters[[arrangement]], 
                                                                                               non_gradient_plots_metadata[[shape]],
                                                                                               plot_labels[[spes_metadata_index]])
      
    }
  }
}

# Put plots into a pdf
setwd("~/R/thesis_plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")
pdf("plots3D_non_gradient.pdf", width = 16, height = 8)

for (metric in metrics) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots3D_non_gradient[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(10, 10, 10, 10))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(arrangements),
                    ncol = length(shapes),
                    byrow = FALSE)
  print(plot)
}
dev.off()


### Get plots for 3D metric analysis (gradient) ------------------------------------------
# Read spes_table
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B"
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 

gradient_plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "gradient", y_aes = "metric", color_aes = "temp_arrangement"),
    bg_prop_A = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_A"),
    bg_prop_B = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_B"),
    E_volume = list(x_aes = "gradient", y_aes = "metric", color_aes = "E_volume")
  ),
  network = list(
    arrangement = list(x_aes = "gradient", y_aes = "metric", color_aes = "temp_arrangement"),
    bg_prop_A = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_A"),
    bg_prop_B = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_B"),
    N_width = list(x_aes = "gradient", y_aes = "metric", color_aes = "N_width")
  )
)

# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("ACIN", "ACINP", "AE", "MS", "NMS", "CKR", "prop_prev", "entropy_prev")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))


metric_plots3D_gradient <- list(mixed_ellipsoid = list(),
                                mixed_network = list(),
                                ringed_ellipsoid = list(),
                                ringed_network = list(),
                                separated_ellipsoid = list(),
                                separated_network = list())

# list for plot labels
plot_labels <- list()
i <- 1
for (shape in shapes) {
  for (arrangement in arrangements) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    plot_labels[[spes_metadata_index]] <- LETTERS[i:(i + 3)]
    i <- i + 4
  }
}

for (shape in shapes) {
  for (arrangement in arrangements) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      
      metric_plots3D_gradient[[spes_metadata_index]][[metric]] <- plot_gradient_metric(spes_table_subset, 
                                                                                       metric,
                                                                                       metric_df_lists3D[[spes_metadata_index]][[metric]], 
                                                                                       arrangement_parameters[[arrangement]], 
                                                                                       get_gradient(metric),
                                                                                       gradient_plots_metadata[[shape]],
                                                                                       plot_labels[[spes_metadata_index]])
      
    }
  }
}




# Put plots into a pdf
setwd("~/R/thesis_plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

metrics <- c("ACIN", "ACINP", "AE", "MS", "NMS", "CKR", "prop_prev", "entropy_prev")

pdf("plots3D_gradient.pdf", width = 17.5, height = 8.5)

for (metric in metrics) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots3D_gradient[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(arrangements),
                    ncol = length(shapes),
                    byrow = FALSE)
  print(plot)
}
dev.off()

### Get plots for 2D metric analysis (non-gradient, all slices with ground truth) -----
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B" 
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 




# Set up plots metadata
plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "temp_arrangement", y_aes = "2D", color_aes = "slice"),
    bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "2D", color_aes = "slice"),
    bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "2D", color_aes = "slice"),
    E_volume = list(x_aes = "E_volume", y_aes = "2D", color_aes = "slice")
  ),
  network = list(
    arrangement = list(x_aes = "temp_arrangement", y_aes = "2D", color_aes = "slice"),
    bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "2D", color_aes = "slice"),
    bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "2D", color_aes = "slice"),
    N_width = list(x_aes = "N_width", y_aes = "2D", color_aes = "slice")
  )
)



# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))

metric_plots2D_non_gradient_all_slices_ground_truth <- list(mixed_ellipsoid = list(),
                                                            mixed_network = list(),
                                                            ringed_ellipsoid = list(),
                                                            ringed_network = list(),
                                                            separated_ellipsoid = list(),
                                                            separated_network = list())

# list for plot labels
plot_labels <- list()
i <- 1
for (shape in shapes) {
  for (arrangement in arrangements) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    plot_labels[[spes_metadata_index]] <- LETTERS[i:(i + 3)]
    i <- i + 4
  }
}

for (shape in shapes) {
  for (arrangement in arrangements) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      metric_plots2D_non_gradient_all_slices_ground_truth[[spes_metadata_index]][[metric]] <- plot_non_gradient_metric_all_slices_ground_truth(spes_table_subset, 
                                                                                                                                               metric, 
                                                                                                                                               metric_df_lists3D[[spes_metadata_index]][[metric]], 
                                                                                                                                               metric_df_lists2D[[spes_metadata_index]][[metric]], 
                                                                                                                                               arrangement_parameters[[arrangement]], 
                                                                                                                                               plots_metadata[[shape]],
                                                                                                                                               plot_labels[[spes_metadata_index]])
    }
  }
}


# Put plots into a pdf
setwd("~/R/thesis_plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")
pdf("plots2D_all_slices_ground_truth.pdf", width = 16, height = 8)

for (metric in metrics) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots2D_non_gradient_all_slices_ground_truth[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(10, 10, 10, 10))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(arrangements),
                    ncol = length(shapes),
                    byrow = FALSE)
  print(plot)
}
dev.off()




### Get plots for error for non-gradient metrics ----------------------------------
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B" 
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 


# Set up plots metadata
plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "temp_arrangement", y_aes = "error", color_aes = "slice"),
    bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "error", color_aes = "slice"),
    bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "error", color_aes = "slice"),
    E_volume = list(x_aes = "E_volume", y_aes = "error", color_aes = "slice")
  ),
  network = list(
    arrangement = list(x_aes = "temp_arrangement", y_aes = "error", color_aes = "slice"),
    bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "error", color_aes = "slice"),
    bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "error", color_aes = "slice"),
    N_width = list(x_aes = "N_width", y_aes = "error", color_aes = "slice")
  )
)



# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))

metric_plots_error_non_gradient <- list(mixed_ellipsoid = list(),
                                        mixed_network = list(),
                                        ringed_ellipsoid = list(),
                                        ringed_network = list(),
                                        separated_ellipsoid = list(),
                                        separated_network = list())

# list for plot labels
plot_labels <- list()
i <- 1
for (shape in shapes) {
  for (arrangement in arrangements) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    plot_labels[[spes_metadata_index]] <- LETTERS[i:(i + 3)]
    i <- i + 4
  }
}


for (shape in shapes) {
  for (arrangement in arrangements) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      metric_plots_error_non_gradient[[spes_metadata_index]][[metric]] <- plot_error_non_gradient_metric(spes_table_subset, 
                                                                                                         metric, 
                                                                                                         metric_df_lists3D[[spes_metadata_index]][[metric]],
                                                                                                         metric_df_lists2D[[spes_metadata_index]][[metric]], 
                                                                                                         arrangement_parameters[[arrangement]], 
                                                                                                         plots_metadata[[shape]],
                                                                                                         plot_labels[[spes_metadata_index]])
    }
  }
}


# Put plots into a pdf
setwd("~/R/thesis_plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")
pdf("plots_error_non_gradient_all_slices.pdf", width = 15.5, height = 8)

for (metric in metrics) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots_error_non_gradient[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(10, 10, 10, 10))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(arrangements),
                    ncol = length(shapes),
                    byrow = FALSE)
  print(plot)
}
dev.off()




### Get plots for violin plots for slice (all slices) -----
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
spes_table <- read.table("spes_table.csv")
spes_table$cluster_prop_B <- 1 - spes_table$cluster_prop_A
spes_table[spes_table$variable_parameter == "cluster_prop_A", "variable_parameter"] <- "cluster_prop_B" 
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 




# Set up plots metadata
plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "slice", y_aes = "metric", label = "temp_arrangement"),
    bg_prop_A = list(x_aes = "slice", y_aes = "metric", label = "bg_prop_A"),
    bg_prop_B = list(x_aes = "slice", y_aes = "metric", label = "bg_prop_B"),
    E_volume = list(x_aes = "slice", y_aes = "metric", label = "E_volume")
  ),
  network = list(
    arrangement = list(x_aes = "slice", y_aes = "metric", label = "temp_arrangement"),
    bg_prop_A = list(x_aes = "slice", y_aes = "metric", label = "bg_prop_A"),
    bg_prop_B = list(x_aes = "slice", y_aes = "metric", label = "bg_prop_B"),
    N_width = list(x_aes = "slice", y_aes = "metric", label = "N_width")
  )
)



# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))

metric_plots_violin_all_slices <- list(mixed_ellipsoid = list(),
                                       mixed_network = list(),
                                       ringed_ellipsoid = list(),
                                       ringed_network = list(),
                                       separated_ellipsoid = list(),
                                       separated_network = list())


# list for plot labels
plot_labels <- list()
i <- 1
for (shape in shapes) {
  for (arrangement in arrangements) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    plot_labels[[spes_metadata_index]] <- LETTERS[i:(i + 3)]
    i <- i + 4
  }
}

for (shape in shapes) {
  for (arrangement in arrangements) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      metric_plots_violin_all_slices[[spes_metadata_index]][[metric]] <- plot_violin_all_slices(spes_table_subset, 
                                                                                                metric, 
                                                                                                metric_df_lists2D[[spes_metadata_index]][[metric]], 
                                                                                                arrangement_parameters[[arrangement]], 
                                                                                                plots_metadata[[shape]],
                                                                                                plot_labels[[spes_metadata_index]])
    }
  }
}




# Put plots into a pdf
setwd("~/R/thesis_plots/S1")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")
pdf("metric_plots_violin_all_slices.pdf", width = 16, height = 8)

for (metric in metrics) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots_violin_all_slices[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(10, 10, 10, 10))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(arrangements),
                    ncol = length(shapes),
                    byrow = FALSE)
  print(plot)
}
dev.off()

