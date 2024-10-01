library(cowplot)
library(ggplot2)

### Utility functions -------------------------
get_gradient <- function(metric) {
  if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR")) return("radius")
  return("threshold")
}

### Read metric_df_lists --------------------------------------------------------------
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
metric_df_lists3D <- readRDS("metric_df_lists3D.RDS")
metric_df_lists2D <- readRDS("metric_df_lists2D.RDS")

### Put prop_AUC and entropy_AUC into metric_df list--------------
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")

thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")

for (arrangement in arrangements) {
  
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    # prop_AUC 3D
    prop_prevalence_df <- metric_df_lists3D[[spes_metadata_index]][["prop_prevalence"]]
    prop_prevalence_df$prop_AUC <- apply(prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    prop_AUC_df <- prop_prevalence_df[ , c("spe", "reference", "target", "prop_AUC")]
    metric_df_lists3D[[spes_metadata_index]][["prop_AUC"]] <- prop_AUC_df
    
    # entropy_AUC 3D
    entropy_prevalence_df <- metric_df_lists3D[[spes_metadata_index]][["entropy_prevalence"]]
    entropy_prevalence_df$entropy_AUC <- apply(entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    entropy_AUC_df <- entropy_prevalence_df[ , c("spe", "cell_types", "entropy_AUC")]
    metric_df_lists3D[[spes_metadata_index]][["entropy_AUC"]] <- entropy_AUC_df
    
    # prop_AUC 2D
    prop_prevalence_df <- metric_df_lists2D[[spes_metadata_index]][["prop_prevalence"]]
    prop_prevalence_df$prop_AUC <- apply(prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    prop_AUC_df <- prop_prevalence_df[ , c("spe", "slice", "reference", "target", "prop_AUC")]
    metric_df_lists2D[[spes_metadata_index]][["prop_AUC"]] <- prop_AUC_df
    
    # entropy_AUC 2D
    entropy_prevalence_df <- metric_df_lists2D[[spes_metadata_index]][["entropy_prevalence"]]
    entropy_prevalence_df$entropy_AUC <- apply(entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    entropy_AUC_df <- entropy_prevalence_df[ , c("spe", "slice", "cell_types", "entropy_AUC")]
    metric_df_lists2D[[spes_metadata_index]][["entropy_AUC"]] <- entropy_AUC_df
  }
}


### Get plots for 3D metric analysis ------------------------------------------
# Read spes_table
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
spes_table <- read.table("spes_table.csv")
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
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))


metric_plots3D <- list(mixed_ellipsoid = list(),
                       mixed_network = list(),
                       ringed_ellipsoid = list(),
                       ringed_network = list(),
                       separated_ellipsoid = list(),
                       separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      if (metric %in% c("AMD", "prop_SAC", "entropy_SAC", "prop_AUC", "entropy_AUC")) {
        metric_plots3D[[spes_metadata_index]][[metric]] <- plot_non_gradient_metric(spes_table_subset, 
                                                                                    metric, 
                                                                                    metric_df_lists3D[[spes_metadata_index]][[metric]], 
                                                                                    arrangement_parameters[[arrangement]], 
                                                                                    non_gradient_plots_metadata[[shape]])
      }
      else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_prevalence", "entropy_prevalence")) {
        metric_plots3D[[spes_metadata_index]][[metric]] <- plot_gradient_metric(spes_table_subset, 
                                                                                metric,
                                                                                metric_df_lists3D[[spes_metadata_index]][[metric]], 
                                                                                arrangement_parameters[[arrangement]], 
                                                                                get_gradient(metric),
                                                                                gradient_plots_metadata[[shape]])
      }
    }
  }
}

# Put plots into a pdf
setwd("~/R/plots/V2")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
# metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")

metrics_set1 <- c("AMD", "ACIN", "CKR")
metrics_set2 <- c("MS", "NMS", "ACINP", "AE", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")

pdf("plots3D.pdf", width = 25, height = 12)

for (metric in metrics_set1) {
  for (shape in shapes) {
    curr_metric_plots <- list()
    for (arrangement in arrangements) {
        spes_metadata_index <- paste(arrangement, shape, sep = "_")
        curr_metric_plots[[arrangement]] <- metric_plots3D[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
      plot <- plot_grid(plotlist = curr_metric_plots,
                        nrow = 1, 
                        ncol = length(arrangements))
      print(plot)
  }
}

for (metric in metrics_set2) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots3D[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(shapes), 
                    ncol = length(arrangements))
  print(plot)
}


dev.off()

### Get plots for 2D metric analysis (middle slice only) ------------------------------------------
# Read spes_table
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
spes_table <- read.table("spes_table.csv")
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 

# Subset metric_df_lists2D to only include the middle slice
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    curr_list <- metric_df_lists2D[[spes_metadata_index]]
    
    for (i in seq(length(curr_list))) {
      curr_df <- curr_list[[i]]
      curr_df <- curr_df[curr_df$slice == 3, ]
      curr_list[[i]] <- curr_df
    }
    metric_df_lists2D[[spes_metadata_index]] <- curr_list
  }
}

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
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))


metric_plots2D <- list(mixed_ellipsoid = list(),
                       mixed_network = list(),
                       ringed_ellipsoid = list(),
                       ringed_network = list(),
                       separated_ellipsoid = list(),
                       separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      if (metric %in% c("AMD", "prop_SAC", "entropy_SAC", "prop_AUC", "entropy_AUC")) {
        metric_plots2D[[spes_metadata_index]][[metric]] <- plot_non_gradient_metric(spes_table_subset, 
                                                                                    metric, 
                                                                                    metric_df_lists2D[[spes_metadata_index]][[metric]], 
                                                                                    arrangement_parameters[[arrangement]], 
                                                                                    non_gradient_plots_metadata[[shape]])
      }
      else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_prevalence", "entropy_prevalence")) {
        metric_plots2D[[spes_metadata_index]][[metric]] <- plot_gradient_metric(spes_table_subset, 
                                                                                metric,
                                                                                metric_df_lists2D[[spes_metadata_index]][[metric]], 
                                                                                arrangement_parameters[[arrangement]], 
                                                                                get_gradient(metric),
                                                                                gradient_plots_metadata[[shape]])
      }
    }
  }
}

# Put plots into a pdf
setwd("~/R/plots/V2")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
# metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")

metrics_set1 <- c("AMD", "ACIN", "CKR")
metrics_set2 <- c("MS", "NMS", "ACINP", "AE", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")

pdf("plots2D_middle_slice.pdf", width = 25, height = 12)

for (metric in metrics_set1) {
  for (shape in shapes) {
    curr_metric_plots <- list()
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[arrangement]] <- metric_plots2D[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
    plot <- plot_grid(plotlist = curr_metric_plots,
                      nrow = 1, 
                      ncol = length(arrangements))
    print(plot)
  }
}

for (metric in metrics_set2) {
  curr_metric_plots <- list()
  for (shape in shapes) {
    for (arrangement in arrangements) {
      spes_metadata_index <- paste(arrangement, shape, sep = "_")
      curr_metric_plots[[spes_metadata_index]] <- metric_plots2D[[spes_metadata_index]][[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
    }
  }
  plot <- plot_grid(plotlist = curr_metric_plots,
                    nrow = length(shapes), 
                    ncol = length(arrangements))
  print(plot)
}

dev.off()

### Get plots with 2D (middle slice only) on the x-axis and 3D on the y-axis ----------------
setwd("~/R/spaSim-3D/scripts/simulations and analysis S1/S1 data")
spes_table <- read.table("spes_table.csv")
spes_table$distance <- 450 - spes_table$cluster_x_coord
spes_table[spes_table$variable_parameter == "cluster_x_coord", "variable_parameter"] <- "distance" 
spes_table$E_volume <- (4/3) * pi * spes_table$E_radius_x * spes_table$E_radius_y * spes_table$E_radius_z
spes_table[spes_table$variable_parameter == "E_radius_z", "variable_parameter"] <- "E_volume" 



# Subset metric_df_lists2D to only include the middle slice
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    curr_list <- metric_df_lists2D[[spes_metadata_index]]
    
    for (i in seq(length(curr_list))) {
      curr_df <- curr_list[[i]]
      curr_df <- curr_df[curr_df$slice == 3, ]
      curr_list[[i]] <- curr_df
    }
    metric_df_lists2D[[spes_metadata_index]] <- curr_list
  }
}




# Set up plots metadata
plots_metadata <- list(
  ellipsoid = list(
    arrangement = list(x_aes = "2D", y_aes = "3D", color_aes = "temp_arrangement"),
    bg_prop_A = list(x_aes = "2D", y_aes = "3D", color_aes = "bg_prop_A"),
    bg_prop_B = list(x_aes = "2D", y_aes = "3D", color_aes = "bg_prop_B"),
    E_volume = list(x_aes = "2D", y_aes = "3D", color_aes = "E_volume")
  ),
  network = list(
    arrangement = list(x_aes = "2D", y_aes = "3D", color_aes = "temp_arrangement"),
    bg_prop_A = list(x_aes = "2D", y_aes = "3D", color_aes = "bg_prop_A"),
    bg_prop_B = list(x_aes = "2D", y_aes = "3D", color_aes = "bg_prop_B"),
    N_width = list(x_aes = "2D", y_aes = "3D", color_aes = "N_width")
  )
)



# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")

background_parameters <- c("bg_prop_A", "bg_prop_B")

arrangement_parameters <- list(mixed = "cluster_prop_B",
                               ringed = "ring_width_factor",
                               separated = "distance")

shape_parameters <- list(ellipsoid = c("E_volume"),
                         network = c("N_width"))

metric_plots2D_vs_3D <- list(mixed_ellipsoid = list(),
                             mixed_network = list(),
                             ringed_ellipsoid = list(),
                             ringed_network = list(),
                             separated_ellipsoid = list(),
                             separated_network = list())

for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    spes_table_subset <- spes_table[spes_table$variable_parameter %in% c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]]), 
                                    c(background_parameters, shape_parameters[[shape]], arrangement_parameters[[arrangement]], "variable_parameter")]
    
    for (metric in metrics) {
      if (metric %in% c("AMD", "prop_SAC", "entropy_SAC", "prop_AUC", "entropy_AUC")) {
        metric_plots2D[[spes_metadata_index]][[metric]] <- plot_3D_vs_2D_metric(spes_table_subset, 
                                                                                metric, 
                                                                                metric_df_lists3D[[spes_metadata_index]][[metric]],
                                                                                metric_df_lists2D[[spes_metadata_index]][[metric]], 
                                                                                arrangement_parameters[[arrangement]], 
                                                                                plots_metadata[[shape]])
      }
      else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_prevalence", "entropy_prevalence")) {
        # metric_plots2D[[spes_metadata_index]][[metric]] <- plot_gradient_metric(spes_table_subset, 
        #                                                                         metric,
        #                                                                         metric_df_lists2D[[spes_metadata_index]][[metric]], 
        #                                                                         arrangement_parameters[[arrangement]], 
        #                                                                         get_gradient(metric),
        #                                                                         gradient_plots_metadata[[shape]])
      }
    }
  }
}









