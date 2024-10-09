### Read data -----

# In raw data, each cell is defined by a number corresponding to its cell type
# Use the cell type dictionary to determine the cell type of a cell from its number
setwd("~/R/Lin et al - human colorectal cancer/other_data")
cell_type_dictionary <- read.csv("cell_type_dictionary.csv")
specific_cell_types <- unique(cell_type_dictionary$Type_Name) # E.g. Tumour/Epi, Ki67+ Tumour/Epi...
generic_cell_types <- unique(cell_type_dictionary$Category) # E.g. Tumour, Stroma, Immune

setwd("~/R/Lin et al - human colorectal cancer/CRC1_data_updated/raw_slice_data")

slice_nums <- c("002", "007", "014", "020", "025", "029",
                "034", "039", "044", "049", "050", "051",
                "052", "054", "059", "064", "069", "074",
                "078", "084", "086", "091", "097")
file_names <- paste("CRC01-", slice_nums, ".csv", sep = "")

# Get z coord of each slice in micrometers using the slice number
get_slice_z_coord <- function(slice_num) {
  
  slice_num <- as.integer(slice_num)
  
  if (slice_num <= 85) {
    return(slice_num * 5)
  }
  
  else {
    return(85 * 5 + (slice_num - 85) * 4)
  }
}

# Get 2D data
# Store each slice as a data.frame, all in a big list
list_dfs2D <- list()

# Loop through each slice
for (i in seq(length(slice_nums))) {
  
  file_name <- file_names[i]
  
  curr_slice_df <- read.csv(file_name, sep = ",")
  curr_slice_df <- curr_slice_df[ , c("NewType", "Xtt", "Ytt")]
  curr_slice_df$Ztt <- get_slice_z_coord(slice_nums[i])
  
  curr_slice_df$Cell.Type.Specific <- specific_cell_types[curr_slice_df$NewType]
  curr_slice_df$Cell.Type.Generic <- generic_cell_types[ifelse(curr_slice_df$NewType <= 3, 1, 
                                                               ifelse(curr_slice_df$NewType <= 6, 2, 3))]
  
  curr_slice_df$NewType <- NULL
  
  list_dfs2D[[i]] <- curr_slice_df
}

# Get 3D data
# Store as one data frame
df3D <- do.call(rbind, list_dfs2D)

# setwd("~/R/Lin et al - human colorectal cancer/CRC1_data_final")
# saveRDS(list_dfs2D, "list_dfs2D.RDS")
# saveRDS(df3D, "df3D.RDS")


### Set up dfs to contain analysis results -----
setwd("~/R/Lin et al - human colorectal cancer/CRC1_data_final")
df3D <- readRDS("df3D.RDS")
colnames(df3D)[1:3] <- c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")
n_slices <- length(unique(df3D$Cell.Z.Position))

# Define AMD data frames as well as constants
cell_types <- c("Tumour", "Immune")
AMD_pairs <- c("Tumour/Tumour", "Tumour/Immune", "Immune/Tumour", "Immune/Immune")

AMD_df_colnames <- c("slice", "reference", "target", "AMD")
AMD_df <- data.frame(matrix(nrow = length(AMD_pairs) * (n_slices + 1), ncol = length(AMD_df_colnames)))
colnames(AMD_df) <- AMD_df_colnames

# Define MS, NMS, ACIN, ACINP, CKR, AE data frames as well as constants
radii <- seq(20, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

MS_df_colnames <- c("slice", "reference", "target", radii_colnames)
MS_df <- data.frame(matrix(nrow = length(cell_types) * (n_slices + 1), ncol = length(MS_df_colnames)))
colnames(MS_df) <- MS_df_colnames

NMS_df_colnames <- c("slice", "reference", "target", radii_colnames)
NMS_df <- data.frame(matrix(nrow = length(cell_types) * (n_slices + 1), ncol = length(MS_df_colnames)))
colnames(NMS_df) <- NMS_df_colnames

# Target is always Tumour and Immune together
# Only choose prop(Tumour) as prop(Immune) = 1 - prop(Tumour) always
ACINP_df_colnames <- c("slice", "reference", radii_colnames)
ACINP_df <- data.frame(matrix(nrow = length(cell_types) * (n_slices + 1), ncol = length(ACINP_df_colnames)))
colnames(ACINP_df) <- ACINP_df_colnames

# Target is always Tumour and Immune together
AE_df_colnames <- c("slice", "reference", radii_colnames)
AE_df <- data.frame(matrix(nrow = length(cell_types) * (n_slices + 1), ncol = length(AE_df_colnames)))
colnames(AE_df) <- AE_df_colnames

## ACIN and CKR are twice as large
# (ref Tumour and tar Tumour or Immune) OR (ref Immune and tar Immune or Tumour)
ACIN_df_colnames <- c("slice", "reference", "target", radii_colnames)
ACIN_df <- data.frame(matrix(nrow = length(cell_types)^2 * (n_slices + 1), ncol = length(ACIN_df_colnames)))
colnames(ACIN_df) <- ACIN_df_colnames

# (ref Tumour and tar Tumour or Immune) OR (ref Immune and tar Immune or Tumour)
CKR_df_colnames <- c("slice", "reference", "target", radii_colnames)
CKR_df <- data.frame(matrix(nrow = length(cell_types)^2 * (n_slices + 1), ncol = length(CKR_df_colnames)))
colnames(CKR_df) <- CKR_df_colnames

# Define SAC and prevalence data frames as well as constants
n_splits <- 10
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

prop_cell_types <- data.frame(ref = c("Tumour"), tar = c("Immune"))

prop_SAC_df_colnames <- c("slice", "reference", "target", "prop_SAC")
prop_SAC_df <- data.frame(matrix(nrow = nrow(prop_cell_types) * (n_slices + 1), ncol = length(prop_SAC_df_colnames)))
colnames(prop_SAC_df) <- prop_SAC_df_colnames

prop_prevalence_df_colnames <- c("slice", "reference", "target", thresholds_colnames)
prop_prevalence_df <- data.frame(matrix(nrow = nrow(prop_cell_types) * (n_slices + 1), ncol = length(prop_prevalence_df_colnames)))
colnames(prop_prevalence_df) <- prop_prevalence_df_colnames


entropy_cell_types <- data.frame(cell_types = c("Tumour,Immune"))

entropy_SAC_df_colnames <- c("slice", "cell_types", "entropy_SAC")
entropy_SAC_df <- data.frame(matrix(nrow = nrow(entropy_cell_types) * (n_slices + 1), ncol = length(entropy_SAC_df_colnames)))
colnames(entropy_SAC_df) <- entropy_SAC_df_colnames

entropy_prevalence_df_colnames <- c("slice", "cell_types", thresholds_colnames)
entropy_prevalence_df <- data.frame(matrix(nrow = nrow(entropy_cell_types) * (n_slices + 1), ncol = length(entropy_prevalence_df_colnames)))
colnames(entropy_prevalence_df) <- entropy_prevalence_df_colnames


# Add all to list:
metric_df_list <- list(AMD = AMD_df,
                       MS = MS_df,
                       NMS = NMS_df,
                       ACINP = ACINP_df,
                       AE = AE_df,
                       ACIN = ACIN_df,
                       CKR = CKR_df,
                       prop_SAC = prop_SAC_df,
                       prop_prevalence = prop_prevalence_df,
                       entropy_SAC = entropy_SAC_df,
                       entropy_prevalence = entropy_prevalence_df)


### Analyse 2D and 3D data -------
slice_z_coords <- unique(df3D$Cell.Z.Position)
feature_colname <- "Cell.Type.Generic"
df3D[df3D$Cell.Type.Generic == "Tumor", "Cell.Type.Generic"] <- "Tumour"


for (i in seq(n_slices + 1)) {
  
  # i represents the current slice index
  if (i != n_slices + 1) {
    df <- df3D[df3D$Cell.Z.Position == slice_z_coords[i], ]
    dimension <- "2D"
    slice <- i
  }
  # if i == n_slices + 1, analyse in 3D instead
  else {
    df <- df3D
    dimension <- "3D"
    slice <- 0
  }
  
  minimum_distance_data <- calculate_minimum_distances_between_cell_types(df,
                                                                          cell_types,
                                                                          feature_colname,
                                                                          dimension = dimension)
  
  minimum_distance_data_summary <- summarise_distances_between_cell_types(minimum_distance_data)
  ## Fill in 4 rows at a time for AMD df (as we have Tumour/Tumour, Tumour/Immune, Immune/Tumour, Immune/IMmune)
  index <- 4 * (i - 1) + 1 # index is 1, 5, 9, 13
  
  metric_df_list[["AMD"]][index:(index + 3), "slice"] <- slice
  metric_df_list[["AMD"]][index:(index + 3), "reference"] <- minimum_distance_data_summary$reference
  metric_df_list[["AMD"]][index:(index + 3), "target"] <- minimum_distance_data_summary$target
  metric_df_list[["AMD"]][index:(index + 3), "AMD"] <- minimum_distance_data_summary$mean
  
  
  
  index1 <- 2 * (i - 1) + 1 # index1 is 1, 3, 5, ...
  index2 <- 4 * (i - 1) + 1 # index2 is 1, 5, 9, 13...
  for (reference_cell_type in cell_types) {
    gradient_data <- calculate_all_gradient_cc_metrics(df,
                                                       reference_cell_type,
                                                       cell_types,
                                                       radii,
                                                       feature_colname,
                                                       dimension = dimension)
    
    target_cell_type <- setdiff(cell_types, reference_cell_type)
    
    if (!is.null(gradient_data)) {
      metric_df_list[["MS"]][index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
      metric_df_list[["NMS"]][index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
      metric_df_list[["ACINP"]][index1, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][["Tumour"]]
      metric_df_list[["AE"]][index1, radii_colnames] <- gradient_data[["entropy"]]$entropy
    }
    else {
      metric_df_list[["MS"]][index1, radii_colnames] <- NA
      metric_df_list[["NMS"]][index1, radii_colnames] <- NA
      metric_df_list[["ACINP"]][index1, radii_colnames] <- NA
      metric_df_list[["AE"]][index1, radii_colnames] <- NA
    }
    
    metric_df_list[["MS"]][index1, "slice"] <- slice
    metric_df_list[["MS"]][index1, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
    
    metric_df_list[["NMS"]][index1, "slice"] <- slice
    metric_df_list[["NMS"]][index1, c("reference", "target")] <- c(reference_cell_type, target_cell_type)
    
    metric_df_list[["ACINP"]][index1, "slice"] <- slice
    metric_df_list[["ACINP"]][index1, c("reference")] <- c(reference_cell_type)
    
    metric_df_list[["AE"]][index1, "slice"] <- slice
    metric_df_list[["AE"]][index1, c("reference")] <- c(reference_cell_type)

    index1 <- index1 + 1
    
    for (target_cell_type in cell_types) {
      ## Calculate ACIN and CKR as target cell type can also be the reference cell type
      
      # ACIN
      metric_df_list[["ACIN"]][index2, "slice"] <- slice
      metric_df_list[["ACIN"]][index2, c("reference", "target")] <- c(reference_cell_type, target_cell_type)

      # CKR
      metric_df_list[["CKR"]][index2, "slice"] <- slice
      metric_df_list[["CKR"]][index2, c("reference", "target")] <- c(reference_cell_type, target_cell_type)

      
      if (!is.null(gradient_data)) {
        metric_df_list[["ACIN"]][index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
        metric_df_list[["CKR"]][index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]]$cross_K_ratio
      }
      else {
        metric_df_list[["ACIN"]][index2, radii_colnames] <- NA
        metric_df_list[["CKR"]][index2, radii_colnames] <- NA
      }
      
      index2 <- index2 + 1
    }
  }
  
  
  # Get proportion grid metrics
  for (j in seq_len(nrow(prop_cell_types))) {
    proportion_grid_metrics <- calculate_cell_proportion_grid_metrics(df, 
                                                                      n_splits,
                                                                      strsplit(prop_cell_types$ref[j], ",")[[1]], 
                                                                      strsplit(prop_cell_types$tar[j], ",")[[1]],
                                                                      feature_colname,
                                                                      dimension = dimension)
    
    
    if (!is.null(proportion_grid_metrics)) {
      proportion_SAC <- calculate_spatial_autocorrelation(proportion_grid_metrics, 
                                                          "proportion",
                                                          weight_method = 0.1,
                                                          dimension = dimension)
      
      proportion_prevalence_df <- calculate_prevalence_gradient(proportion_grid_metrics,
                                                                "proportion")
    }
    else {
      proportion_SAC <- NA
      proportion_prevalence_df <- data.frame(threshold = seq(0.01, 1, 0.01), prevalence = NA)
    }

    index <- nrow(prop_cell_types) * (i - 1) + j
    metric_df_list[["prop_SAC"]][index, "slice"] <- slice
    metric_df_list[["prop_SAC"]][index, c("reference", "target")] <- c(prop_cell_types$ref[j], prop_cell_types$tar[j])
    metric_df_list[["prop_SAC"]][index, "prop_SAC"] <- proportion_SAC
    
    metric_df_list[["prop_prevalence"]][index, "slice"] <- slice
    metric_df_list[["prop_prevalence"]][index, c("reference", "target")] <- c(prop_cell_types$ref[j], prop_cell_types$tar[j])
    metric_df_list[["prop_prevalence"]][index, thresholds_colnames] <- proportion_prevalence_df$prevalence
  }
  
  # Get entropy grid metrics
  for (j in seq_len(nrow(entropy_cell_types))) {
    entropy_grid_metrics <- calculate_entropy_grid_metrics(df, 
                                                           n_splits,
                                                           strsplit(entropy_cell_types$cell_types[j], ",")[[1]],
                                                           feature_colname,
                                                           dimension = dimension)
    
    if (!is.null(entropy_grid_metrics)) {
      entropy_SAC <- calculate_spatial_autocorrelation(entropy_grid_metrics, 
                                                       "entropy",
                                                       weight_method = 0.1,
                                                       dimension = dimension)
      
      entropy_prevalence_df <- calculate_prevalence_gradient(entropy_grid_metrics,
                                                             "entropy")
      
    }
    else {
      entropy_SAC <- NA
      entropy_prevalence_df <- data.frame(threshold = seq(0.01, 1, 0.01), prevalence = NA)
    }
    
    index <- nrow(entropy_cell_types) * (i - 1) + j
    metric_df_list[["entropy_SAC"]][index, "slice"] <- slice
    metric_df_list[["entropy_SAC"]][index, c("cell_types")] <- c(entropy_cell_types$cell_types[j])
    metric_df_list[["entropy_SAC"]][index, "entropy_SAC"] <- entropy_SAC
    
    metric_df_list[["entropy_prevalence"]][index, "slice"] <- slice
    metric_df_list[["entropy_prevalence"]][index, c("cell_types")] <- c(entropy_cell_types$cell_types[j])
    metric_df_list[["entropy_prevalence"]][index, thresholds_colnames] <- entropy_prevalence_df$prevalence
  }  
}


# setwd("~/R/Lin et al - human colorectal cancer/CRC1_data_final")
# saveRDS(metric_df_list, "metric_df_list.RDS")

### Plot analysis of 2D and 3D data -----
setwd("~/R/Lin et al - human colorectal cancer/CRC1_data_final")
metric_df_list <- readRDS("metric_df_list.RDS")

get_gradient <- function(metric) {
  if (metric %in% c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR")) {
    return("radius")
  }
  else if (metric %in% c("prop_prevalence", "entropy_prevalence")) {
    return("threshold")  
  }
  else {
    stop("Invalid metric. Must be gradient-based")
  }
}


## Turn gradient radii metrics into AUC and add to metric_df list
get_AUC_for_radii_gradient_metrics <- function(y) {
  x <- radii
  h <- diff(x)[1]
  n <- length(x)
  
  AUC <- (h / 2) * (y[1] + 2 * sum(y[2:(n - 1)]) + y[n])
  
  return(AUC)
}


radii <- seq(20, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

gradient_radii_metrics <- c("MS", "NMS", "ACINP", "AE", "ACIN", "CKR")


for (metric in gradient_radii_metrics) {
  metric_AUC_name <- paste(metric, "AUC", sep = "_")
  
  if (metric %in% c("MS", "NMS", "ACIN", "CKR")) {
    subset_colnames <- c("slice", "reference", "target", metric_AUC_name)
  }
  else {
    subset_colnames <- c("slice", "reference", metric_AUC_name)
  }
  
  df <- metric_df_list[[metric]]
  df[[metric_AUC_name]] <- apply(df[ , radii_colnames], 1, get_AUC_for_radii_gradient_metrics)
  metric_df_list[[metric_AUC_name]] <- df
  
}

## Turn threshold radii metrics into AUC and add to metric_df list
thresholds <- seq(0.01, 1, 0.01)
threshold_colnames <- paste("t", thresholds, sep = "")

# prop_AUC 3D
prop_prevalence_df <- metric_df_list[["prop_prevalence"]]
prop_prevalence_df$prop_AUC <- apply(prop_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
prop_AUC_df <- prop_prevalence_df[ , c("slice", "reference", "target", "prop_AUC")]
metric_df_list[["prop_AUC"]] <- prop_AUC_df

# entropy_AUC 3D
entropy_prevalence_df <- metric_df_list[["entropy_prevalence"]]
entropy_prevalence_df$entropy_AUC <- apply(entropy_prevalence_df[ , threshold_colnames], 1, sum) * 0.01
entropy_AUC_df <- entropy_prevalence_df[ , c("slice", "cell_types", "entropy_AUC")]
metric_df_list[["entropy_AUC"]] <- entropy_AUC_df


## Functions to plot
# Utility function to get metric cell types
get_metric_cell_types <- function(metric) {
  # Get metric_cell_types
  if (metric %in% c("AMD", "ACIN", "CKR", "ACIN_AUC", "CKR_AUC")) {
    metric_cell_types <- data.frame(ref = c("Tumour"), tar = c("Immune"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("MS", "NMS", "MS_AUC", "NMS_AUC")) {
    metric_cell_types <- data.frame(ref = c("Tumour"), tar = c("Immune"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("ACINP", "ACINP_AUC")) {
    metric_cell_types <- data.frame(ref = c("Tumour"), tar = c("Tumour"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("AE", "AE_AUC")) {
    metric_cell_types <- data.frame(ref = c("Tumour"), tar = c("Tumour,Immune"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("prop_SAC", "prop_prevalence", "prop_AUC")) {
    metric_cell_types <- data.frame(ref = c("Tumour"), tar = c("Immune"))
    metric_cell_types$pair <- paste(metric_cell_types$ref, metric_cell_types$tar, sep = "/")
  }
  else if (metric %in% c("entropy_SAC", "entropy_prevalence", "entropy_AUC")) {
    metric_cell_types <- data.frame(cell_types = c("Tumour,Immune"))
  }
  else {
    stop("metric not found")
  }
  return(metric_cell_types)
}

# Utility function to subset metric_df
subset_metric_df <- function(metric_df,
                             metric_cell_types,
                             metric,
                             index) {
  
  if (metric %in% c("AMD", "ACIN", "CKR", "MS", "NMS", "ACIN_AUC", "CKR_AUC", "MS_AUC", "NMS_AUC", "prop_SAC", "prop_prevalence", "prop_AUC")) {
    metric_df_subset <- metric_df[metric_df$reference == metric_cell_types[index, "ref"] & metric_df$target == metric_cell_types[index, "tar"], ] 
  }
  else if (metric %in% c("ACINP", "AE", "ACINP_AUC", "AE_AUC")) {
    metric_df_subset <- metric_df[metric_df$reference == metric_cell_types[index, "ref"], ] 
  }
  else if (metric %in% c("entropy_SAC", "entropy_prevalence", "entropy_AUC")) {
    metric_df_subset <- metric_df[metric_df$cell_types == metric_cell_types[index, "cell_types"], ]
  }
  else {
    stop("metric not found")
  }
  
  return(metric_df_subset)
}



plotting_function <- function(metric_df_list,
                              metrics) {
  
  # Get all dfs for each metric and put them into a combined df
  combined_df <- data.frame()
  
  for (metric in metrics) {

    # Get metric_df for current metric
    metric_df <- metric_df_list[[metric]]

    # Get metric cell types for current metric (should only be one set/ one row)
    metric_cell_types <- get_metric_cell_types(metric)
    
    # Subset metric_df
    metric_df_subset <- subset_metric_df(metric_df,
                                         metric_cell_types,
                                         metric,
                                         1) # Always first row
    
    # Change and further subset columns of metric_df_subset
    colnames(metric_df_subset)[colnames(metric_df_subset) == metric] <- "value"
    metric_df_subset$metric <- metric
    metric_df_subset <- metric_df_subset[ , c("slice", "value", "metric")]
    
    # Add to combined df
    combined_df <- rbind(combined_df, metric_df_subset)
  }
  combined_df$dummy <- "dummy"
  
  # Create the dot plot, highlighting the maximum slice points with a star shape and using facets
  fig <- ggplot(combined_df[combined_df$slice != 0, ], aes(x = dummy, y = value)) +
    geom_jitter(width = 0.2, height = 0, size = 2) +
    geom_point(data = combined_df[combined_df$slice == 0, ], color = "red", shape = 8, size = 6) +
    labs(x = "Metric", y = "Value") +
    facet_wrap(~ metric, strip.position = "bottom", scales = "free_y", ncol = 4) +
    theme_bw()  +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  
  return(fig)
  
}

## Get plot
metrics <- c("AMD", "MS_AUC", "NMS_AUC", "ACINP_AUC", "AE_AUC", "ACIN_AUC", "CKR_AUC", "prop_SAC", "prop_AUC", "entropy_SAC", "entropy_AUC")

plot <- plotting_function(metric_df_list,
                          metrics)
methods::show(plot)

setwd("~/R/Lin et al - human colorectal cancer/CRC1_data_final")
pdf("lin_et_al_3D_vs_2D.pdf", width = 10, height = 8)

print(plot)

dev.off()



