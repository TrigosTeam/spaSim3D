library(cowplot)
library(ggplot2)

### Read tables --------------------------------------------------------------
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/spe_tables")
mixed_spes_table <- read.table("mixed_spes_table.csv")
ringed_spes_table <- read.table("ringed_spes_table.csv")
separated_spes_table <- read.table("separated_spes_table.csv")
separated_spes_table$distance <- separated_spes_table$centre_x_coord_B - separated_spes_table$centre_x_coord_A

### Set up plot lists --------------------------------------------------------

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




### Mixed spes --------------
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/mixed")

thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")
mixed_dfs <- list()

for (metric in metrics) {
  if (metric == "prop_AUC") {
    mixed_prop_prevalence_df <- mixed_dfs[["prop_prevalence"]]
    mixed_prop_prevalence_df$prop_AUC <- apply(mixed_prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    mixed_prop_AUC_df <- mixed_prop_prevalence_df[ , c("spe", "reference", "target", "prop_AUC")]
  }
  else if (metric == "entropy_AUC") {
    mixed_entropy_prevalence_df <- mixed_dfs[["entropy_prevalence"]]
    mixed_entropy_prevalence_df$entropy_AUC <- apply(mixed_entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    mixed_entropy_AUC_df <- mixed_entropy_prevalence_df[ , c("spe", "cell_types", "entropy_AUC")]
  }
  else {
    mixed_dfs[[metric]] <- read.table(paste("mixed_", metric, "_df.csv", sep = "")) 
  }
}

mixed_plots <- list()

for (metric in metrics) {
  if (metric %in% c("AMD", "prop_SAC", "entropy_SAC", "prop_AUC", "entropy_AUC")) {
    mixed_plots[[metric]] <- plot_non_gradient_metric(mixed_spes_table, 
                                                      metric, 
                                                      mixed_dfs[[metric]], 
                                                      "cluster_prop_B", 
                                                      non_gradient_plots_metadata)
  }
  else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "entropy_SAC", "prop_prevalence", "entropy_prevalence")) {
    mixed_plots[[metric]] <- plot_gradient_metric(mixed_spes_table, 
                                                  metric,
                                                  mixed_dfs[[metric]], 
                                                  "cluster_prop_B", 
                                                  "radius",
                                                  gradient_plots_metadata)
  }
}

setwd("~/R/plots/3D")
saveRDS(mixed_plots, "mixed_plots_3D.RDS")

### Ringed spes --------------
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/ringed")

thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")
ringed_dfs <- list()

for (metric in metrics) {
  if (metric == "prop_AUC") {
    ringed_prop_prevalence_df <- ringed_dfs[["prop_prevalence"]]
    ringed_prop_prevalence_df$prop_AUC <- apply(ringed_prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    ringed_prop_AUC_df <- ringed_prop_prevalence_df[ , c("spe", "reference", "target", "prop_AUC")]
  }
  else if (metric == "entropy_AUC") {
    ringed_entropy_prevalence_df <- ringed_dfs[["entropy_prevalence"]]
    ringed_entropy_prevalence_df$entropy_AUC <- apply(ringed_entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    ringed_entropy_AUC_df <- ringed_entropy_prevalence_df[ , c("spe", "cell_types", "entropy_AUC")]
  }
  else {
    ringed_dfs[[metric]] <- read.table(paste("ringed_", metric, "_df.csv", sep = "")) 
  }
}

ringed_plots <- list()

for (metric in metrics) {
  if (metric %in% c("AMD", "prop_SAC", "entropy_SAC", "prop_AUC", "entropy_AUC")) {
    ringed_plots[[metric]] <- plot_non_gradient_metric(ringed_spes_table, 
                                                      metric, 
                                                      ringed_dfs[[metric]], 
                                                      "ring_width_factor", 
                                                      non_gradient_plots_metadata)
  }
  else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "entropy_SAC", "prop_prevalence", "entropy_prevalence")) {
    ringed_plots[[metric]] <- plot_gradient_metric(ringed_spes_table, 
                                                  metric,
                                                  ringed_dfs[[metric]], 
                                                  "ring_width_factor", 
                                                  "radius",
                                                  gradient_plots_metadata)
  }
}

setwd("~/R/plots/3D")
saveRDS(ringed_plots, "ringed_plots_3D.RDS")

### Separated spes --------------
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/separated")

thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")
separated_dfs <- list()

for (metric in metrics) {
  if (metric == "prop_AUC") {
    separated_prop_prevalence_df <- separated_dfs[["prop_prevalence"]]
    separated_prop_prevalence_df$prop_AUC <- apply(separated_prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    separated_prop_AUC_df <- separated_prop_prevalence_df[ , c("spe", "reference", "target", "prop_AUC")]
  }
  else if (metric == "entropy_AUC") {
    separated_entropy_prevalence_df <- separated_dfs[["entropy_prevalence"]]
    separated_entropy_prevalence_df$entropy_AUC <- apply(separated_entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    separated_entropy_AUC_df <- separated_entropy_prevalence_df[ , c("spe", "cell_types", "entropy_AUC")]
  }
  else {
    separated_dfs[[metric]] <- read.table(paste("separated_", metric, "_df.csv", sep = "")) 
  }
}

separated_plots <- list()

for (metric in metrics) {
  if (metric %in% c("AMD", "prop_SAC", "entropy_SAC", "prop_AUC", "entropy_AUC")) {
    separated_plots[[metric]] <- plot_non_gradient_metric(separated_spes_table, 
                                                      metric, 
                                                      separated_dfs[[metric]], 
                                                      "distance", 
                                                      non_gradient_plots_metadata)
  }
  else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "entropy_SAC", "prop_prevalence", "entropy_prevalence")) {
    separated_plots[[metric]] <- plot_gradient_metric(separated_spes_table, 
                                                  metric,
                                                  separated_dfs[[metric]], 
                                                  "distance", 
                                                  "radius",
                                                  gradient_plots_metadata)
  }
}

setwd("~/R/plots/3D")
saveRDS(separated_plots, "separated_plots_3D.RDS")

### Get pdf plot ------------------------------------------------------

setwd("~/R/plots/3D")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")
arrangements <- c("mixed", "ringed", "separated")

mixed_plots <- readRDS("mixed_plots_3D.RDS")
ringed_plots <- readRDS("ringed_plots_3D.RDS")
separated_plots <- readRDS("separated_plots_3D.RDS")

pdf("plots.pdf", width = 15, height = 10)

for (metric in metrics) {
  for (arrangement in arrangements) {
    if (arrangement == "mixed") {
      print(mixed_plots[[metric]])
    }
    else if (arrangement == "ringed") {
      print(ringed_plots[[metric]])
    }
    else if (arrangement == "separated") {
      print(separated_plots[[metric]])
    }
  }
}
dev.off()

### Without background noise -------------------------------------------------
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/spe_tables")
mixed_spes_table <- read.table("mixed_spes_table.csv")
ringed_spes_table <- read.table("ringed_spes_table.csv")
separated_spes_table <- read.table("separated_spes_table.csv")
separated_spes_table$distance <- 450 - separated_spes_table$centre_x_coord
mixed_spes_table <- mixed_spes_table[mixed_spes_table$bg_prop_A == 0 & mixed_spes_table$bg_prop_B == 0, ]
ringed_spes_table <- ringed_spes_table[ringed_spes_table$bg_prop_A == 0 & ringed_spes_table$bg_prop_B == 0, ]
separated_spes_table <- separated_spes_table[separated_spes_table$bg_prop_A == 0 & separated_spes_table$bg_prop_B == 0, ]

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


### Mixed spes
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/mixed")

thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")
mixed_dfs <- list()

for (metric in metrics) {
  if (metric == "prop_AUC") {
    mixed_prop_prevalence_df <- mixed_dfs[["prop_prevalence"]]
    mixed_prop_prevalence_df$prop_AUC <- apply(mixed_prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    mixed_prop_AUC_df <- mixed_prop_prevalence_df[ , c("spe", "reference", "target", "prop_AUC")]
  }
  else if (metric == "entropy_AUC") {
    mixed_entropy_prevalence_df <- mixed_dfs[["entropy_prevalence"]]
    mixed_entropy_prevalence_df$entropy_AUC <- apply(mixed_entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    mixed_entropy_AUC_df <- mixed_entropy_prevalence_df[ , c("spe", "cell_types", "entropy_AUC")]
  }
  else {
    df <- read.table(paste("mixed_", metric, "_df.csv", sep = "")) 
    df <- df[df$spe %in% paste("mixed_spe_", rownames(mixed_spes_table), sep = ""), ]
    mixed_dfs[[metric]] <- df
  }
}

mixed_plots <- list()

for (metric in metrics) {
  if (metric %in% c("AMD", "prop_SAC", "entropy_SAC", "prop_AUC", "entropy_AUC")) {
    mixed_plots[[metric]] <- plot_non_gradient_metric(mixed_spes_table, 
                                                      metric, 
                                                      mixed_dfs[[metric]], 
                                                      "cluster_prop_B", 
                                                      non_gradient_plots_metadata)
  }
  else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "entropy_SAC", "prop_prevalence", "entropy_prevalence")) {
    mixed_plots[[metric]] <- plot_gradient_metric(mixed_spes_table, 
                                                  metric,
                                                  mixed_dfs[[metric]], 
                                                  "cluster_prop_B", 
                                                  "radius",
                                                  gradient_plots_metadata)
  }
}

setwd("~/R/plots/3D_no_bg_noise")
saveRDS(mixed_plots, "mixed_plots_3D_no_bg_noise.RDS")

### Ringed spes
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/ringed")

thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")
ringed_dfs <- list()

for (metric in metrics) {
  if (metric == "prop_AUC") {
    ringed_prop_prevalence_df <- ringed_dfs[["prop_prevalence"]]
    ringed_prop_prevalence_df$prop_AUC <- apply(ringed_prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    ringed_prop_AUC_df <- ringed_prop_prevalence_df[ , c("spe", "reference", "target", "prop_AUC")]
  }
  else if (metric == "entropy_AUC") {
    ringed_entropy_prevalence_df <- ringed_dfs[["entropy_prevalence"]]
    ringed_entropy_prevalence_df$entropy_AUC <- apply(ringed_entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    ringed_entropy_AUC_df <- ringed_entropy_prevalence_df[ , c("spe", "cell_types", "entropy_AUC")]
  }
  else {
    df <- read.table(paste("ringed_", metric, "_df.csv", sep = "")) 
    df <- df[df$spe %in% paste("ringed_spe_", rownames(ringed_spes_table), sep = ""), ]
    ringed_dfs[[metric]] <- df
  }
}

ringed_plots <- list()

for (metric in metrics) {
  if (metric %in% c("AMD", "prop_SAC", "entropy_SAC", "prop_AUC", "entropy_AUC")) {
    ringed_plots[[metric]] <- plot_non_gradient_metric(ringed_spes_table, 
                                                       metric, 
                                                       ringed_dfs[[metric]], 
                                                       "ring_width_factor", 
                                                       non_gradient_plots_metadata)
  }
  else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "entropy_SAC", "prop_prevalence", "entropy_prevalence")) {
    ringed_plots[[metric]] <- plot_gradient_metric(ringed_spes_table, 
                                                   metric,
                                                   ringed_dfs[[metric]], 
                                                   "ring_width_factor", 
                                                   "radius",
                                                   gradient_plots_metadata)
  }
}

setwd("~/R/plots/3D_no_bg_noise")
saveRDS(ringed_plots, "ringed_plots_3D_no_bg_noise.RDS")

### Separated spes
setwd("~/R/spaSim-3D/scripts/spaSim-3D & SPIAT-3D testing/analysis3D_tables/separated")

thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")
separated_dfs <- list()

for (metric in metrics) {
  if (metric == "prop_AUC") {
    separated_prop_prevalence_df <- separated_dfs[["prop_prevalence"]]
    separated_prop_prevalence_df$prop_AUC <- apply(separated_prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    separated_prop_AUC_df <- separated_prop_prevalence_df[ , c("spe", "reference", "target", "prop_AUC")]
  }
  else if (metric == "entropy_AUC") {
    separated_entropy_prevalence_df <- separated_dfs[["entropy_prevalence"]]
    separated_entropy_prevalence_df$entropy_AUC <- apply(separated_entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
    separated_entropy_AUC_df <- separated_entropy_prevalence_df[ , c("spe", "cell_types", "entropy_AUC")]
  }
  else {
    df <- read.table(paste("separated_", metric, "_df.csv", sep = "")) 
    df <- df[df$spe %in% paste("separated_spe_", rownames(separated_spes_table), sep = ""), ]
    separated_dfs[[metric]] <- df 
  }
}

separated_plots <- list()

for (metric in metrics) {
  if (metric %in% c("AMD", "prop_SAC", "entropy_SAC", "prop_AUC", "entropy_AUC")) {
    separated_plots[[metric]] <- plot_non_gradient_metric(separated_spes_table, 
                                                          metric, 
                                                          separated_dfs[[metric]], 
                                                          "distance", 
                                                          non_gradient_plots_metadata)
  }
  else if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "entropy_SAC", "prop_prevalence", "entropy_prevalence")) {
    separated_plots[[metric]] <- plot_gradient_metric(separated_spes_table, 
                                                      metric,
                                                      separated_dfs[[metric]], 
                                                      "distance", 
                                                      "radius",
                                                      gradient_plots_metadata)
  }
}

setwd("~/R/plots/3D_no_bg_noise")
saveRDS(separated_plots, "separated_plots_3D_no_bg_noise.RDS")


setwd("~/R/plots/3D_no_bg_noise")
metrics <- c("AMD", "MS", "NMS", "ACINP", "AE", "ACIN", "CKR", "prop_SAC", "prop_prevalence", "prop_AUC", "entropy_SAC", "entropy_prevalence", "entropy_AUC")
arrangements <- c("mixed", "ringed", "separated")

mixed_plots <- readRDS("mixed_plots_3D_no_bg_noise.RDS")
ringed_plots <- readRDS("ringed_plots_3D_no_bg_noise.RDS")
separated_plots <- readRDS("separated_plots_3D_no_bg_noise.RDS")

pdf("plots.pdf", width = 15, height = 10)

for (metric in metrics) {
  for (arrangement in arrangements) {
    if (arrangement == "mixed") {
      print(mixed_plots[[metric]])
    }
    else if (arrangement == "ringed") {
      print(ringed_plots[[metric]])
    }
    else if (arrangement == "separated") {
      print(separated_plots[[metric]])
    }
  }
}
dev.off()




