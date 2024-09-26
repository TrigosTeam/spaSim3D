library(cowplot)
library(ggplot2)

### Read tables --------------------------------------------------------------
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/spe_tables")
mixed_spes_table <- read.table("mixed_spes_table.csv")
ringed_spes_table <- read.table("ringed_spes_table.csv")
separated_spes_table <- read.table("separated_spes_table.csv")
separated_spes_table$distance <- 450 - separated_spes_table$cluster_x_coord

arrangement_tables <- list(mixed = mixed_spes_table,
                           ringed = ringed_spes_table,
                           separated = separated_spes_table)

### Set up plots metadata --------------------------------------------------------
non_gradient_plots_metadata <- list(
  arrangement = list(x_aes = "temp_arrangement", y_aes = "metric"),
  bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "metric"),
  bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "metric"),
  shape = list(x_aes = "shape", y_aes = "metric"),
  variation_E = list(x_aes = "variation_E", y_aes = "metric"),
  volume_E = list(x_aes = "volume_E", y_aes = "metric"),
  width_N = list(x_aes = "width_N", y_aes = "metric")
)

gradient_plots_metadata <- list(
  arrangement = list(x_aes = "gradient", y_aes = "metric", color_aes = "temp_arrangement"),
  bg_prop_A = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_A"),
  bg_prop_B = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_B"),
  shape = list(x_aes = "gradient", y_aes = "metric", color_aes = "shape"),
  variation_E = list(x_aes = "gradient", y_aes = "metric", color_aes = "variation_E"),
  volume_E = list(x_aes = "gradient", y_aes = "metric", color_aes = "volume_E"),
  width_N = list(x_aes = "gradient", y_aes = "metric", color_aes = "width_N")
)




### Utility functions -------------------------
get_gradient <- function(metric) {
  if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR")) return("radius")
  return("threshold")
}
### Put plots into list--------------
thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")
arrangements <- c("mixed", "ringed", "separated")
arrangement_parameters <- list(mixed = "cluster_prop_B", ringed = "ring_width_factor", separated = "distance")

metric_dfs <- list(mixed = list(),
                   ringed = list(),
                   separated = list())

for (arrangement in arrangements) {
  setwd(paste("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/", arrangement,  sep = ""))
  for (metric in metrics) {
    if (metric == "prop_AUC") {
      prop_prevalence_df <- metric_dfs[[arrangement]][["prop_prevalence"]]
      prop_prevalence_df$prop_AUC <- apply(prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
      prop_AUC_df <- prop_prevalence_df[ , c("spe", "reference", "target", "prop_AUC")]
      metric_dfs[[arrangement]][[metric]] <- prop_AUC_df
    }
    else if (metric == "entropy_AUC") {
      entropy_prevalence_df <- metric_dfs[[arrangement]][["entropy_prevalence"]]
      entropy_prevalence_df$entropy_AUC <- apply(entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
      entropy_AUC_df <- entropy_prevalence_df[ , c("spe", "cell_types", "entropy_AUC")]
      metric_dfs[[arrangement]][[metric]] <- entropy_AUC_df
    }
    else {
      metric_dfs[[arrangement]][[metric]] <- read.table(paste(arrangement, metric, "df.csv", sep = "_")) 
    }
  }
}


metric_plots <- list(mixed = list(),
                     ringed = list(),
                     separated = list())

for (arrangement in arrangements) {
  for (metric in metrics) {
    if (metric %in% c("AMD", "prop_SAC", "entropy_SAC", "prop_AUC", "entropy_AUC")) {
      metric_plots[[arrangement]][[metric]] <- plot_non_gradient_metric(arrangement_tables[[arrangement]], 
                                                                        metric, 
                                                                        metric_dfs[[arrangement]][[metric]], 
                                                                        arrangement_parameters[[arrangement]], 
                                                                        non_gradient_plots_metadata)
    }
    else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_prevalence", "entropy_prevalence")) {
      metric_plots[[arrangement]][[metric]] <- plot_gradient_metric(arrangement_tables[[arrangement]], 
                                                                    metric,
                                                                    metric_dfs[[arrangement]][[metric]], 
                                                                    arrangement_parameters[[arrangement]], 
                                                                    get_gradient(metric),
                                                                    gradient_plots_metadata)
    }
  }
}

setwd("~/R/plots/3D")
saveRDS(metric_plots, "metric_plots_3D.RDS")



### Get pdf plot ------------------------------------------------------

setwd("~/R/plots/3D")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")
arrangements <- c("mixed", "ringed", "separated")
metric_plots <- readRDS("metric_plots_3D.RDS")

pdf("plots.pdf", width = 20, height = 12)

for (metric in metrics) {
  for (arrangement in arrangements) {
    print(metric_plots[[arrangement]][[metric]])
  }
}
dev.off()

### Without background noise -------------------------------------------------
# Read tables
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/spe_tables")
mixed_spes_table <- read.table("mixed_spes_table.csv")
ringed_spes_table <- read.table("ringed_spes_table.csv")
separated_spes_table <- read.table("separated_spes_table.csv")
separated_spes_table$distance <- 450 - separated_spes_table$cluster_x_coord
mixed_spes_table <- mixed_spes_table[mixed_spes_table$bg_prop_A == 0 & mixed_spes_table$bg_prop_B == 0, ]
ringed_spes_table <- ringed_spes_table[ringed_spes_table$bg_prop_A == 0 & ringed_spes_table$bg_prop_B == 0, ]
separated_spes_table <- separated_spes_table[separated_spes_table$bg_prop_A == 0 & separated_spes_table$bg_prop_B == 0, ]

arrangement_tables <- list(mixed = mixed_spes_table,
                           ringed = ringed_spes_table,
                           separated = separated_spes_table)

# Set up plots metadata
non_gradient_plots_metadata <- list(
  arrangement = list(x_aes = "temp_arrangement", y_aes = "metric"),
  # bg_prop_A = list(x_aes = "bg_prop_A", y_aes = "metric"),
  # bg_prop_B = list(x_aes = "bg_prop_B", y_aes = "metric"),
  shape = list(x_aes = "shape", y_aes = "metric"),
  variation_E = list(x_aes = "variation_E", y_aes = "metric"),
  volume_E = list(x_aes = "volume_E", y_aes = "metric"),
  width_N = list(x_aes = "width_N", y_aes = "metric")
)

gradient_plots_metadata <- list(
  arrangement = list(x_aes = "gradient", y_aes = "metric", color_aes = "temp_arrangement"),
  # bg_prop_A = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_A"),
  # bg_prop_B = list(x_aes = "gradient", y_aes = "metric", color_aes = "bg_prop_B"),
  shape = list(x_aes = "gradient", y_aes = "metric", color_aes = "shape"),
  variation_E = list(x_aes = "gradient", y_aes = "metric", color_aes = "variation_E"),
  volume_E = list(x_aes = "gradient", y_aes = "metric", color_aes = "volume_E"),
  width_N = list(x_aes = "gradient", y_aes = "metric", color_aes = "width_N")
)


# Utility functions
get_gradient <- function(metric) {
  if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR")) return("radius")
  return("threshold")
}



# Put plots into list
thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")
arrangements <- c("mixed", "ringed", "separated")
arrangement_parameters <- list(mixed = "cluster_prop_B", ringed = "ring_width_factor", separated = "distance")

metric_dfs <- list(mixed = list(),
                   ringed = list(),
                   separated = list())

for (arrangement in arrangements) {
  setwd(paste("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/", arrangement,  sep = ""))
  for (metric in metrics) {
    if (metric == "prop_AUC") {
      prop_prevalence_df <- metric_dfs[[arrangement]][["prop_prevalence"]]
      prop_prevalence_df$prop_AUC <- apply(prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
      prop_AUC_df <- prop_prevalence_df[ , c("spe", "reference", "target", "prop_AUC")]
      metric_dfs[[arrangement]][[metric]] <- prop_AUC_df
    }
    else if (metric == "entropy_AUC") {
      entropy_prevalence_df <- metric_dfs[[arrangement]][["entropy_prevalence"]]
      entropy_prevalence_df$entropy_AUC <- apply(entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
      entropy_AUC_df <- entropy_prevalence_df[ , c("spe", "cell_types", "entropy_AUC")]
      metric_dfs[[arrangement]][[metric]] <- entropy_AUC_df
    }
    else {
      df <- read.table(paste(arrangement, metric, "df.csv", sep = "_")) 
      df <- df[df$spe %in% paste(arrangement, "_spe_", rownames(arrangement_tables[[arrangement]]), sep = ""), ]
      metric_dfs[[arrangement]][[metric]] <- df
    }
  }
}


metric_plots <- list(mixed = list(),
                     ringed = list(),
                     separated = list())

for (arrangement in arrangements) {
  for (metric in metrics) {
    if (metric %in% c("AMD", "prop_SAC", "entropy_SAC", "prop_AUC", "entropy_AUC")) {
      metric_plots[[arrangement]][[metric]] <- plot_non_gradient_metric(arrangement_tables[[arrangement]], 
                                                                        metric, 
                                                                        metric_dfs[[arrangement]][[metric]], 
                                                                        arrangement_parameters[[arrangement]], 
                                                                        non_gradient_plots_metadata)
    }
    else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_prevalence", "entropy_prevalence")) {
      metric_plots[[arrangement]][[metric]] <- plot_gradient_metric(arrangement_tables[[arrangement]], 
                                                                    metric,
                                                                    metric_dfs[[arrangement]][[metric]], 
                                                                    arrangement_parameters[[arrangement]], 
                                                                    get_gradient(metric),
                                                                    gradient_plots_metadata)
    }
  }
}
setwd("~/R/plots/3D_no_bg_noise")
saveRDS(metric_plots, "metric_plots_3D_no_bg_noise.RDS")



# Get pdf plot
setwd("~/R/plots/3D_no_bg_noise")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")
arrangements <- c("mixed", "ringed", "separated")

metric_plots <- readRDS("metric_plots_3D_no_bg_noise.RDS")

pdf("plots.pdf", width = 15, height = 12)

for (metric in metrics) {
  for (arrangement in arrangements) {
    print(metric_plots[[arrangement]][[metric]])
  }
}
dev.off()


