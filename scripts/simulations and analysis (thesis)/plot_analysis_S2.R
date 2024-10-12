library(cowplot)
library(ggplot2)

### Utility functions -------------------------
get_gradient <- function(metric) {
  if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR")) return("radius")
  return("threshold")
}

### Read metric_df_lists --------------------------------------------------------------
setwd("~/R/spaSim-3D/scripts/simulations and analysis S2/S2 data")
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


### Get plots with 2D (all slices) on the x-axis and 3D on the y-axis (not annotating for arrangement or shape and choosing random slice) ----------------

# Set up plots metadata
plots_metadata <- list(
  temp <- list(x_aes = "3D", y_aes = "2D")
)

# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")


# Merge lists in metric_lists
metric_df_lists3D_merged <- list()
metric_df_lists2D_merged <- list()

i <- 1
for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (metric in metrics) {
      if (i == 1)  {
        metric_df_lists3D_merged[[metric]] <- data.frame()
        metric_df_lists2D_merged[[metric]] <- data.frame()
      }
      if (i > 1) {
        temp <- nrow(metric_df_lists3D[[spes_metadata_index]][[metric]])
        n_slices <- length(unique(metric_df_lists2D[[spes_metadata_index]][[metric]][["slice"]]))
        if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC")) {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp/4) * (i - 1) + (temp/4)), each = 4), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp)/4 * (i - 1) + (temp/4)), each = 4 * n_slices), sep = "_")
        }
        else {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2 * n_slices), sep = "_")
        }
      }
      metric_df_lists3D_merged[[metric]] <- rbind(metric_df_lists3D_merged[[metric]], metric_df_lists3D[[spes_metadata_index]][[metric]])
      metric_df_lists2D_merged[[metric]] <- rbind(metric_df_lists2D_merged[[metric]], metric_df_lists2D[[spes_metadata_index]][[metric]])
    }
    
    i <- i + 1
  }
}

metric_plots_3D_vs_2D_random_slice <- list()

for (metric in metrics) {
  metric_plots_3D_vs_2D_random_slice[[metric]] <- plot_3D_vs_2D_metric_random_slice_no_annotating(metric, 
                                                                                                  metric_df_lists3D_merged[[metric]],
                                                                                                  metric_df_lists2D_merged[[metric]], 
                                                                                                  plots_metadata)
}



# Put plots into a pdf
setwd("~/R/thesis_plots/S2")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

plot_list <- list()

pdf("plots_3D_vs_2D_random_slice.pdf", width = 10, height = 8)

for (metric in metrics) {
  plot_list[[metric]] <- metric_plots_3D_vs_2D_random_slice[[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
}

plots <- plot_grid(plotlist = plot_list,
                   nrow = 3,
                   ncol = 4,
                   labels = LETTERS[1:13])

print(plots)
dev.off()



### Get plots with 2D (all slices) on the x-axis and 3D and ERROR on the y-axis (not annotating for arrangement or shape and choosing random slice) ----------------

# Set up plots metadata
plots_metadata <- list(
  temp <- list(x_aes = "3D", y_aes = "error")
)

# Generate plots and plots into a list
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")


# Merge lists in metric_lists
metric_df_lists3D_merged <- list()
metric_df_lists2D_merged <- list()

i <- 1
for (arrangement in arrangements) {
  for (shape in shapes) {
    spes_metadata_index <- paste(arrangement, shape, sep = "_")
    
    for (metric in metrics) {
      if (i == 1)  {
        metric_df_lists3D_merged[[metric]] <- data.frame()
        metric_df_lists2D_merged[[metric]] <- data.frame()
      }
      if (i > 1) {
        temp <- nrow(metric_df_lists3D[[spes_metadata_index]][[metric]])
        n_slices <- length(unique(metric_df_lists2D[[spes_metadata_index]][[metric]][["slice"]]))
        if (metric %in% c("AMD", "ACIN_AUC", "CKR_AUC")) {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp/4) * (i - 1) + (temp/4)), each = 4), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/4) * (i - 1) + 1, (temp)/4 * (i - 1) + (temp/4)), each = 4 * n_slices), sep = "_")
        }
        else {
          metric_df_lists3D[[spes_metadata_index]][[metric]][["spe"]] <- 
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2), sep = "_")
          metric_df_lists2D[[spes_metadata_index]][[metric]][["spe"]] <-
            paste("spe", rep(seq((temp/2) * (i - 1) + 1, (temp/2) * (i - 1) + (temp/2)), each = 2 * n_slices), sep = "_")
        }
      }
      metric_df_lists3D_merged[[metric]] <- rbind(metric_df_lists3D_merged[[metric]], metric_df_lists3D[[spes_metadata_index]][[metric]])
      metric_df_lists2D_merged[[metric]] <- rbind(metric_df_lists2D_merged[[metric]], metric_df_lists2D[[spes_metadata_index]][[metric]])
    }
    
    i <- i + 1
  }
}

metric_plots_3D_vs_2D_random_slice <- list()
metric_plots_error_vs_2D_random_slice <- list()

for (metric in metrics) {
  temp_plots <- plot_3D_and_error_vs_2D_metric_random_slice_no_annotating(metric, 
                                                                          metric_df_lists3D_merged[[metric]],
                                                                          metric_df_lists2D_merged[[metric]], 
                                                                          plots_metadata)
  
  metric_plots_3D_vs_2D_random_slice[[metric]] <- temp_plots$plots_3D_vs_2D
  metric_plots_error_vs_2D_random_slice[[metric]] <- temp_plots$plots_error_vs_2D
}



# Put plots into a pdf
setwd("~/R/thesis_plots/S2")
arrangements <- c("mixed", "ringed", "separated")
shapes <- c("ellipsoid", "network")
metrics <- c("AMD", "ACIN_AUC", "ACINP_AUC", "AE_AUC", "MS_AUC", "NMS_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

## 3D vs 2D
plot_list <- list()

pdf("plots_3D_vs_2D_random_slice.pdf", width = 10, height = 8)

for (metric in metrics) {
  plot_list[[metric]] <- metric_plots_3D_vs_2D_random_slice[[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
}

plots <- plot_grid(plotlist = plot_list,
                   nrow = 3,
                   ncol = 4,
                   labels = LETTERS[1:13])
print(plots)
dev.off()


## Error vs 2D
plot_list <- list()

pdf("plots_error_vs_2D_random_slice.pdf", width = 10, height = 8)

for (metric in metrics) {
  plot_list[[metric]] <- metric_plots_error_vs_2D_random_slice[[metric]] + theme(plot.margin = margin(15, 15, 15, 15))  
}

plots <- plot_grid(plotlist = plot_list,
                   nrow = 3,
                   ncol = 4,
                   labels = LETTERS[1:13])
print(plots)
dev.off()





