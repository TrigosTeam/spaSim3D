### 1.1. Mixed spes - cc distance metrics ----------------------------------------

# Get mixed_spes_table
setwd("~/Objects/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Get number of mixed spes
n_mixed_spes <- nrow(mixed_spes_table)

# Define APD and AMD data frames as well as constants
cell_types <- c("A", "B")
APD_pairs <- c("A/A", "A/B", "B/B")
AMD_pairs <- c("A/A", "A/B", "B/A", "B/B")

mixed_APD_df <- data.frame(matrix(nrow = n_mixed_spes * length(APD_pairs), ncol = 3))
colnames(mixed_APD_df) <- c("spe", "pair", "APD")

mixed_AMD_df <- data.frame(matrix(nrow = n_mixed_spes * length(AMD_pairs), ncol = 4))
colnames(mixed_AMD_df) <- c("spe", "reference", "target", "AMD")


# Loop through each mixed spes and get APD and AMD
setwd("~/Objects/mixed_spes")
for (i in seq(n_mixed_spes)) {

  # Read in current mixed spe
  mixed_spe_name <- paste("mixed_spe_", i, sep = "")
  mixed_spe_file_name <- paste(mixed_spe_name, ".rds", sep = "")
  mixed_spe <- readRDS(mixed_spe_file_name)

  # pairwise_distance_data <- calculate_pairwise_distances_between_cell_types3D(mixed_spe,
  #                                                                             cell_types,
  #                                                                             show_summary = F,
  #                                                                             plot_image = F)
  # pairwise_distance_data_summary <- summarise_distances_between_cell_types3D(pairwise_distance_data)
  # 
  # ## Fill in 3 rows at a time for APD df (as we have A/A, A/B, B/B)
  # index <- 3 * (i - 1) + 1 # index is 1, 4, 7, 10 ...
  # mixed_APD_df[index:(index + 2), "spe"] <- mixed_spe_name
  # mixed_APD_df[index:(index + 2), "pair"] <- pairwise_distance_data_summary$pair
  # mixed_APD_df[index:(index + 2), "APD"] <- pairwise_distance_data_summary$mean
  
  minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(mixed_spe,
                                                                            cell_types,
                                                                            show_summary = F,
                                                                            plot_image = F)
  
  minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)
  ## Fill in 4 rows at a time for AMD df (as we have A/A, A/B, B/A, B/B)
  index <- 4 * (i - 1) + 1 # index is 1, 5, 9, 13
  mixed_AMD_df[index:(index + 3), "spe"] <- mixed_spe_name
  mixed_AMD_df[index:(index + 3), "reference"] <- minimum_distance_data_summary$reference
  mixed_AMD_df[index:(index + 3), "target"] <- minimum_distance_data_summary$target
  mixed_AMD_df[index:(index + 3), "AMD"] <- minimum_distance_data_summary$mean
  
}

setwd("~/Objects/mixed_spes/analysis_3D")
write.table(mixed_APD_df, file = "mixed_APD_df.csv")
write.table(mixed_AMD_df, file = "mixed_AMD_df.csv")

### 1.2. Mixed spes - cc gradient based metrics ----------------------------------

# Get mixed_spes_table
setwd("~/Objects/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Get number of mixed spes
n_mixed_spes <- nrow(mixed_spes_table)

# Define MS, NMS, ACIN, ACINP, CKR, AE data frames as well as constants
cell_types <- c("A", "B") # Use A as reference, and B as target, and vice versa
radii <- 50
radii_colnames <- paste("r", seq(radii), sep = "")

mixed_MS_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types), ncol = 3 + radii))
colnames(mixed_MS_df) <- c("spe", "reference", "target", radii_colnames)

mixed_NMS_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types), ncol = 3 + radii))
colnames(mixed_NMS_df) <- c("spe", "reference", "target", radii_colnames)

# Target is always A and B together
# Only choose prop(A) as prop(B) = 1 - prop(A) always
mixed_ACINP_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types), ncol = 2 + radii))
colnames(mixed_ACINP_df) <- c("spe", "reference", radii_colnames)

# Target is always A and B together
mixed_AE_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types), ncol = 2 + radii))
colnames(mixed_AE_df) <- c("spe", "reference", radii_colnames)


## ACIN and CKR are twice as large
# (ref A and tar A or B) OR (ref B and tar B or A)
mixed_ACIN_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types)^2, ncol = 3 + radii))
colnames(mixed_ACIN_df) <- c("spe", "reference", "target", radii_colnames)

# (ref A and tar A or B) OR (ref B and tar B or A)
mixed_CKR_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types)^2, ncol = 3 + radii))
colnames(mixed_CKR_df) <- c("spe", "reference", "target", radii_colnames)



# Loop through each mixed spes and get gradient-based metrics
setwd("~/Objects/mixed_spes")
for (i in seq(n_mixed_spes)) {
  
  # Read in current mixed spe
  mixed_spe_name <- paste("mixed_spe_", i, sep = "")
  mixed_spe_file_name <- paste(mixed_spe_name, ".rds", sep = "")
  mixed_spe <- readRDS(mixed_spe_file_name)
  
  
  index1 <- 2 * (i - 1) + 1 # index1 is 1, 3, 5, ...
  index2 <- 4 * (i - 1) + 1 # index2 is 1, 5, 9, 13...
  for (reference_cell_type in cell_types) {
    ## Calculate MS, NMS, ACINP and AE first
    target_cell_type <- setdiff(cell_types, reference_cell_type)
    
    # MS and NMS
    MS_NMS_data <- calculate_mixing_scores_gradient3D(mixed_spe,
                                                      reference_cell_type,
                                                      target_cell_type,
                                                      radii,
                                                      plot_image = F)
    mixed_MS_df[index1, c("spe", "reference", "target")] <- c(mixed_spe_name, reference_cell_type, target_cell_type)
    mixed_MS_df[index1, radii_colnames] <- MS_NMS_data$mixing_score
    
    mixed_NMS_df[index1, c("spe", "reference", "target")] <- c(mixed_spe_name, reference_cell_type, target_cell_type)
    mixed_NMS_df[index1, radii_colnames] <- MS_NMS_data$normalised_mixing_score
    
    
    # ACINP
    ACINP_data <- calculate_cells_in_neighbourhood_proportions_gradient3D(mixed_spe,
                                                                          reference_cell_type,
                                                                          cell_types, # Use both cell types, but keep prop for A
                                                                          radii,
                                                                          plot_image = F)
    mixed_ACINP_df[index1, c("spe", "reference")] <- c(mixed_spe_name, reference_cell_type)
    mixed_ACINP_df[index1, radii_colnames] <- ACINP_data[["A"]]
    
    
    # AE
    AE_data <- calculate_entropy_gradient3D(mixed_spe,
                                            reference_cell_type,
                                            cell_types, # Use both A and B
                                            radii,
                                            plot_image = F)
    mixed_AE_df[index1, c("spe", "reference")] <- c(mixed_spe_name, reference_cell_type)
    mixed_AE_df[index1, radii_colnames] <- AE_data$entropy
    
    index1 <- index1 + 1
    
    
    for (target_cell_type in cell_types) {
      ## Calculate ACIN and CKR as target cell type can also be the reference cell type
      
      # ACIN
      ACIN_data <- calculate_cells_in_neighbourhood_gradient3D(mixed_spe,
                                                               reference_cell_type,
                                                               target_cell_type,
                                                               radii,
                                                               plot_image = F)
      mixed_ACIN_df[index2, c("spe", "reference", "target")] <- c(mixed_spe_name, reference_cell_type, target_cell_type)
      mixed_ACIN_df[index2, radii_colnames] <- ACIN_data[[target_cell_type]]
      
      
      # CKR
      CK_data <- calculate_cross_K_gradient3D(mixed_spe,
                                              reference_cell_type,
                                              target_cell_type,
                                              radii,
                                              plot_image = F)
      mixed_CKR_df[index2, c("spe", "reference", "target")] <- c(mixed_spe_name, reference_cell_type, target_cell_type)
      mixed_CKR_df[index2, radii_colnames] <- CK_data$observed_cross_K / CK_data$expected_cross_K
      
      index2 <- index2 + 1
    }
  }
}

setwd("~/Objects/mixed_spes/analysis_3D")
write.table(mixed_MS_df, file = "mixed_MS_df.csv")
write.table(mixed_NMS_df, file = "mixed_NMS_df.csv")
write.table(mixed_ACINP_df, file = "mixed_ACINP_df.csv")
write.table(mixed_AE_df, file = "mixed_AE_df.csv")
write.table(mixed_ACIN_df, file = "mixed_ACIN_df.csv")
write.table(mixed_CKR_df, file = "mixed_CKR_df.csv")


### 1.3. Mixed spes - heterogeneity metrics --------------------------------------
# Get mixed_spes_table
setwd("~/Objects/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Get number of mixed spes
n_mixed_spes <- nrow(mixed_spes_table)

# Define SAC and prevalence data frames as well as constants
n_splits <- 10
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))

mixed_prop_SAC_df_colnames <- c("spe", "reference", "target", "proportion")
mixed_prop_SAC_df <- data.frame(matrix(nrow = n_mixed_spes * nrow(prop_cell_types), ncol = length(mixed_prop_SAC_df_colnames)))
colnames(mixed_prop_SAC_df) <- mixed_prop_SAC_df_colnames

mixed_prop_prevalence_df_colnames <- c("spe", "reference", "target", thresholds_colnames)
mixed_prop_prevalence_df <- data.frame(matrix(nrow = n_mixed_spes * nrow(prop_cell_types), ncol = length(mixed_prop_prevalence_df_colnames)))
colnames(mixed_prop_prevalence_df) <- mixed_prop_prevalence_df_colnames


entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))

mixed_entropy_SAC_df_colnames <- c("spe", "cell_types", "entropy")
mixed_entropy_SAC_df <- data.frame(matrix(nrow = n_mixed_spes * nrow(entropy_cell_types), ncol = length(mixed_entropy_SAC_df_colnames)))
colnames(mixed_entropy_SAC_df) <- mixed_entropy_SAC_df_colnames

mixed_entropy_prevalence_df_colnames <- c("spe", "cell_types", thresholds_colnames)
mixed_entropy_prevalence_df <- data.frame(matrix(nrow = n_mixed_spes * nrow(entropy_cell_types), ncol = length(mixed_entropy_prevalence_df_colnames)))
colnames(mixed_entropy_prevalence_df) <- mixed_entropy_prevalence_df_colnames


# Loop through each mixed spes and get SAC and prevalence
setwd("~/Objects/mixed_spes")
for (i in seq(n_mixed_spes)) {
  
  # Read in current mixed spe
  mixed_spe_name <- paste("mixed_spe_", i, sep = "")
  mixed_spe_file_name <- paste(mixed_spe_name, ".rds", sep = "")
  mixed_spe <- readRDS(mixed_spe_file_name)
  
  # Get proportion grid metrics
  for (j in seq_len(nrow(prop_cell_types))) {
    proportion_grid_metrics <- determine_cell_proportion_grid_metrics3D(mixed_spe, 
                                                                        n_splits,
                                                                        strsplit(prop_cell_types$ref[j], ",")[[1]], 
                                                                        strsplit(prop_cell_types$tar[j], ",")[[1]],
                                                                        plot_image = F)
    
    proportion_SAC <- determine_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                          "proportion",
                                                          "binary")
    
    proportion_prevalence_df <- determine_prevalence_gradient3D(proportion_grid_metrics,
                                                                "proportion",
                                                                show_AUC = T,
                                                                plot_image = F)
    
    index <- nrow(prop_cell_types) * (i - 1) + j
    mixed_prop_SAC_df[index, c("spe", "reference", "target")] <- c(mixed_spe_name, prop_cell_types$ref[j], prop_cell_types$tar[j])
    mixed_prop_SAC_df[index, "proportion"] <- proportion_SAC
    
    mixed_prop_prevalence_df[index, c("spe", "reference", "target")] <- c(mixed_spe_name, prop_cell_types$ref[j], prop_cell_types$tar[j])
    mixed_prop_prevalence_df[index, thresholds_colnames] <- proportion_prevalence_df$prevalence
  }
  
  # Get entropy grid metrics
  for (j in seq_len(nrow(entropy_cell_types))) {
    entropy_grid_metrics <- determine_entropy_grid_metrics3D(mixed_spe, 
                                                             n_splits,
                                                             strsplit(entropy_cell_types$cell_types[j], ",")[[1]], 
                                                             plot_image = F)
    
    entropy_SAC <- determine_spatial_autocorrelation3D(entropy_grid_metrics, 
                                                       "entropy",
                                                       "binary")
    
    entropy_prevalence_df <- determine_prevalence_gradient3D(entropy_grid_metrics,
                                                             "entropy",
                                                             show_AUC = T,
                                                             plot_image = F)
    
    index <- nrow(entropy_cell_types) * (i - 1) + j
    mixed_entropy_SAC_df[index, c("spe", "cell_types")] <- c(mixed_spe_name, entropy_cell_types$cell_types[j])
    mixed_entropy_SAC_df[index, "entropy"] <- entropy_SAC
    
    mixed_entropy_prevalence_df[index, c("spe", "cell_types")] <- c(mixed_spe_name, entropy_cell_types$cell_types[j])
    mixed_entropy_prevalence_df[index, thresholds_colnames] <- entropy_prevalence_df$prevalence
  }
}

setwd("~/Objects/mixed_spes/analysis_3D")
# write.table(mixed_prop_SAC_df, file = "mixed_prop_SAC_df.csv")
# write.table(mixed_prop_prevalence_df, file = "mixed_prop_prevalence_df.csv")
# write.table(mixed_entropy_SAC_df, file = "mixed_entropy_SAC_df.csv")
# write.table(mixed_entropy_prevalence_df, file = "mixed_entropy_prevalence_df.csv")


### 2.1. Ringed spes - cc distance metrics ----------------------------------------

# Get ringed_spes_table
setwd("~/Objects/spes_table")
ringed_spes_table <- read.table("ringed_spes_table.csv")

# Get number of ringed spes
n_ringed_spes <- nrow(ringed_spes_table)

# Define APD and AMD data frames as well as constants
cell_types <- c("A", "B")
APD_pairs <- c("A/A", "A/B", "B/B")
AMD_pairs <- c("A/A", "A/B", "B/A", "B/B")

ringed_APD_df <- data.frame(matrix(nrow = n_ringed_spes * length(APD_pairs), ncol = 3))
colnames(ringed_APD_df) <- c("spe", "pair", "APD")

ringed_AMD_df <- data.frame(matrix(nrow = n_ringed_spes * length(AMD_pairs), ncol = 4))
colnames(ringed_AMD_df) <- c("spe", "reference", "target", "AMD")


# Loop through each ringed spes and get APD and AMD
setwd("~/Objects/ringed_spes")
for (i in seq(n_ringed_spes)) {
  
  # Read in current ringed spe
  ringed_spe_name <- paste("ringed_spe_", i, sep = "")
  ringed_spe_file_name <- paste(ringed_spe_name, ".rds", sep = "")
  ringed_spe <- readRDS(ringed_spe_file_name)
  
  # pairwise_distance_data <- calculate_pairwise_distances_between_cell_types3D(ringed_spe,
  #                                                                             cell_types,
  #                                                                             show_summary = F,
  #                                                                             plot_image = F)
  # pairwise_distance_data_summary <- summarise_distances_between_cell_types3D(pairwise_distance_data)
  # 
  # ## Fill in 3 rows at a time for APD df (as we have A/A, A/B, B/B)
  # index <- 3 * (i - 1) + 1 # index is 1, 4, 7, 10 ...
  # ringed_APD_df[index:(index + 2), "spe"] <- ringed_spe_name
  # ringed_APD_df[index:(index + 2), "pair"] <- pairwise_distance_data_summary$pair
  # ringed_APD_df[index:(index + 2), "APD"] <- pairwise_distance_data_summary$mean
  
  minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(ringed_spe,
                                                                            cell_types,
                                                                            show_summary = F,
                                                                            plot_image = F)
  
  minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)
  ## Fill in 4 rows at a time for AMD df (as we have A/A, A/B, B/A, B/B)
  index <- 4 * (i - 1) + 1 # index is 1, 5, 9, 13
  ringed_AMD_df[index:(index + 3), "spe"] <- ringed_spe_name
  ringed_AMD_df[index:(index + 3), "reference"] <- minimum_distance_data_summary$reference
  ringed_AMD_df[index:(index + 3), "target"] <- minimum_distance_data_summary$target
  ringed_AMD_df[index:(index + 3), "AMD"] <- minimum_distance_data_summary$mean
}

setwd("~/Objects/ringed_spes/analysis_3D")
write.table(ringed_APD_df, file = "ringed_APD_df.csv")
write.table(ringed_AMD_df, file = "ringed_AMD_df.csv")

### 2.2. Ringed spes - cc gradient based metrics ----------------------------------

# Get ringed_spes_table
setwd("~/Objects/spes_table")
ringed_spes_table <- read.table("ringed_spes_table.csv")

# Get number of ringed spes
n_ringed_spes <- nrow(ringed_spes_table)

# Define MS, NMS, ACIN, ACINP, CKR, AE data frames as well as constants
cell_types <- c("A", "B") # Use A as reference, and B as target, and vice versa
radii <- 50
radii_colnames <- paste("r", seq(radii), sep = "")

ringed_MS_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types), ncol = 3 + radii))
colnames(ringed_MS_df) <- c("spe", "reference", "target", radii_colnames)

ringed_NMS_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types), ncol = 3 + radii))
colnames(ringed_NMS_df) <- c("spe", "reference", "target", radii_colnames)

# Target is always A and B together
# Only choose prop(A) as prop(B) = 1 - prop(A) always
ringed_ACINP_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types), ncol = 2 + radii))
colnames(ringed_ACINP_df) <- c("spe", "reference", radii_colnames)

# Target is always A and B together
ringed_AE_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types), ncol = 2 + radii))
colnames(ringed_AE_df) <- c("spe", "reference", radii_colnames)


## ACIN and CKR are twice as large
# (ref A and tar A or B) OR (ref B and tar B or A)
ringed_ACIN_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types)^2, ncol = 3 + radii))
colnames(ringed_ACIN_df) <- c("spe", "reference", "target", radii_colnames)

# (ref A and tar A or B) OR (ref B and tar B or A)
ringed_CKR_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types)^2, ncol = 3 + radii))
colnames(ringed_CKR_df) <- c("spe", "reference", "target", radii_colnames)



# Loop through each ringed spes and get gradient-based metrics
setwd("~/Objects/ringed_spes")
for (i in seq(n_ringed_spes)) {
  
  # Read in current ringed spe
  ringed_spe_name <- paste("ringed_spe_", i, sep = "")
  ringed_spe_file_name <- paste(ringed_spe_name, ".rds", sep = "")
  ringed_spe <- readRDS(ringed_spe_file_name)
  
  
  index1 <- 2 * (i - 1) + 1 # index1 is 1, 3, 5, ...
  index2 <- 4 * (i - 1) + 1 # index2 is 1, 5, 9, 13...
  for (reference_cell_type in cell_types) {
    ## Calculate MS, NMS, ACINP and AE first
    target_cell_type <- setdiff(cell_types, reference_cell_type)
    
    # MS and NMS
    MS_NMS_data <- calculate_mixing_scores_gradient3D(ringed_spe,
                                                      reference_cell_type,
                                                      target_cell_type,
                                                      radii,
                                                      plot_image = F)
    ringed_MS_df[index1, c("spe", "reference", "target")] <- c(ringed_spe_name, reference_cell_type, target_cell_type)
    ringed_MS_df[index1, radii_colnames] <- MS_NMS_data$mixing_score
    
    ringed_NMS_df[index1, c("spe", "reference", "target")] <- c(ringed_spe_name, reference_cell_type, target_cell_type)
    ringed_NMS_df[index1, radii_colnames] <- MS_NMS_data$normalised_mixing_score
    
    
    # ACINP
    ACINP_data <- calculate_cells_in_neighbourhood_proportions_gradient3D(ringed_spe,
                                                                          reference_cell_type,
                                                                          cell_types, # Use both cell types, but keep prop for A
                                                                          radii,
                                                                          plot_image = F)
    ringed_ACINP_df[index1, c("spe", "reference")] <- c(ringed_spe_name, reference_cell_type)
    ringed_ACINP_df[index1, radii_colnames] <- ACINP_data[["A"]]
    
    
    # AE
    AE_data <- calculate_entropy_gradient3D(ringed_spe,
                                            reference_cell_type,
                                            cell_types, # Use both A and B
                                            radii,
                                            plot_image = F)
    ringed_AE_df[index1, c("spe", "reference")] <- c(ringed_spe_name, reference_cell_type)
    ringed_AE_df[index1, radii_colnames] <- AE_data$entropy
    
    index1 <- index1 + 1
    
    
    for (target_cell_type in cell_types) {
      ## Calculate ACIN and CKR as target cell type can also be the reference cell type
      
      # ACIN
      ACIN_data <- calculate_cells_in_neighbourhood_gradient3D(ringed_spe,
                                                               reference_cell_type,
                                                               target_cell_type,
                                                               radii,
                                                               plot_image = F)
      ringed_ACIN_df[index2, c("spe", "reference", "target")] <- c(ringed_spe_name, reference_cell_type, target_cell_type)
      ringed_ACIN_df[index2, radii_colnames] <- ACIN_data[[target_cell_type]]
      
      
      # CKR
      CK_data <- calculate_cross_K_gradient3D(ringed_spe,
                                              reference_cell_type,
                                              target_cell_type,
                                              radii,
                                              plot_image = F)
      ringed_CKR_df[index2, c("spe", "reference", "target")] <- c(ringed_spe_name, reference_cell_type, target_cell_type)
      ringed_CKR_df[index2, radii_colnames] <- CK_data$observed_cross_K / CK_data$expected_cross_K
      
      index2 <- index2 + 1
    }
  }
}

setwd("~/Objects/ringed_spes/analysis_3D")
write.table(ringed_MS_df, file = "ringed_MS_df.csv")
write.table(ringed_NMS_df, file = "ringed_NMS_df.csv")
write.table(ringed_ACINP_df, file = "ringed_ACINP_df.csv")
write.table(ringed_AE_df, file = "ringed_AE_df.csv")
write.table(ringed_ACIN_df, file = "ringed_ACIN_df.csv")
write.table(ringed_CKR_df, file = "ringed_CKR_df.csv")


### 2.3. Ringed spes - heterogeneity metrics --------------------------------------
# Get ringed_spes_table
setwd("~/Objects/spes_table")
ringed_spes_table <- read.table("ringed_spes_table.csv")

# Get number of ringed spes
n_ringed_spes <- nrow(ringed_spes_table)

# Define SAC and prevalence data frames as well as constants
n_splits <- 10
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))

ringed_prop_SAC_df_colnames <- c("spe", "reference", "target", "proportion")
ringed_prop_SAC_df <- data.frame(matrix(nrow = n_ringed_spes * nrow(prop_cell_types), ncol = length(ringed_prop_SAC_df_colnames)))
colnames(ringed_prop_SAC_df) <- ringed_prop_SAC_df_colnames

ringed_prop_prevalence_df_colnames <- c("spe", "reference", "target", thresholds_colnames)
ringed_prop_prevalence_df <- data.frame(matrix(nrow = n_ringed_spes * nrow(prop_cell_types), ncol = length(ringed_prop_prevalence_df_colnames)))
colnames(ringed_prop_prevalence_df) <- ringed_prop_prevalence_df_colnames


entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))

ringed_entropy_SAC_df_colnames <- c("spe", "cell_types", "entropy")
ringed_entropy_SAC_df <- data.frame(matrix(nrow = n_ringed_spes * nrow(entropy_cell_types), ncol = length(ringed_entropy_SAC_df_colnames)))
colnames(ringed_entropy_SAC_df) <- ringed_entropy_SAC_df_colnames

ringed_entropy_prevalence_df_colnames <- c("spe", "cell_types", thresholds_colnames)
ringed_entropy_prevalence_df <- data.frame(matrix(nrow = n_ringed_spes * nrow(entropy_cell_types), ncol = length(ringed_entropy_prevalence_df_colnames)))
colnames(ringed_entropy_prevalence_df) <- ringed_entropy_prevalence_df_colnames


# Loop through each ringed spes and get SAC and prevalence
setwd("~/Objects/ringed_spes")
for (i in seq(n_ringed_spes)) {
  
  # Read in current ringed spe
  ringed_spe_name <- paste("ringed_spe_", i, sep = "")
  ringed_spe_file_name <- paste(ringed_spe_name, ".rds", sep = "")
  ringed_spe <- readRDS(ringed_spe_file_name)
  
  # Get proportion grid metrics
  for (j in seq_len(nrow(prop_cell_types))) {
    proportion_grid_metrics <- determine_cell_proportion_grid_metrics3D(ringed_spe, 
                                                                        n_splits,
                                                                        strsplit(prop_cell_types$ref[j], ",")[[1]], 
                                                                        strsplit(prop_cell_types$tar[j], ",")[[1]],
                                                                        plot_image = F)
    
    proportion_SAC <- determine_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                          "proportion",
                                                          "binary")
    
    proportion_prevalence_df <- determine_prevalence_gradient3D(proportion_grid_metrics,
                                                                "proportion",
                                                                show_AUC = T,
                                                                plot_image = F)
    
    index <- nrow(prop_cell_types) * (i - 1) + j
    ringed_prop_SAC_df[index, c("spe", "reference", "target")] <- c(ringed_spe_name, prop_cell_types$ref[j], prop_cell_types$tar[j])
    ringed_prop_SAC_df[index, "proportion"] <- proportion_SAC
    
    ringed_prop_prevalence_df[index, c("spe", "reference", "target")] <- c(ringed_spe_name, prop_cell_types$ref[j], prop_cell_types$tar[j])
    ringed_prop_prevalence_df[index, thresholds_colnames] <- proportion_prevalence_df$prevalence
  }
  
  # Get entropy grid metrics
  for (j in seq_len(nrow(entropy_cell_types))) {
    entropy_grid_metrics <- determine_entropy_grid_metrics3D(ringed_spe, 
                                                             n_splits,
                                                             strsplit(entropy_cell_types$cell_types[j], ",")[[1]], 
                                                             plot_image = F)
    
    entropy_SAC <- determine_spatial_autocorrelation3D(entropy_grid_metrics, 
                                                       "entropy",
                                                       "binary")
    
    entropy_prevalence_df <- determine_prevalence_gradient3D(entropy_grid_metrics,
                                                             "entropy",
                                                             show_AUC = T,
                                                             plot_image = F)
    
    index <- nrow(entropy_cell_types) * (i - 1) + j
    ringed_entropy_SAC_df[index, c("spe", "cell_types")] <- c(ringed_spe_name, entropy_cell_types$cell_types[j])
    ringed_entropy_SAC_df[index, "entropy"] <- entropy_SAC
    
    ringed_entropy_prevalence_df[index, c("spe", "cell_types")] <- c(ringed_spe_name, entropy_cell_types$cell_types[j])
    ringed_entropy_prevalence_df[index, thresholds_colnames] <- entropy_prevalence_df$prevalence
  }
}

setwd("~/Objects/ringed_spes/analysis_3D")
write.table(ringed_prop_SAC_df, file = "ringed_prop_SAC_df.csv")
write.table(ringed_prop_prevalence_df, file = "ringed_prop_prevalence_df.csv")
write.table(ringed_entropy_SAC_df, file = "ringed_entropy_SAC_df.csv")
write.table(ringed_entropy_prevalence_df, file = "ringed_entropy_prevalence_df.csv")


### 3.1. Separated spes - cc distance metrics ----------------------------------------

# Get separated_spes_table
setwd("~/Objects/spes_table")
separated_spes_table <- read.table("separated_spes_table.csv")

# Get number of separated spes
n_separated_spes <- nrow(separated_spes_table)

# Define APD and AMD data frames as well as constants
cell_types <- c("A", "B")
APD_pairs <- c("A/A", "A/B", "B/B")
AMD_pairs <- c("A/A", "A/B", "B/A", "B/B")

separated_APD_df <- data.frame(matrix(nrow = n_separated_spes * length(APD_pairs), ncol = 3))
colnames(separated_APD_df) <- c("spe", "pair", "APD")

separated_AMD_df <- data.frame(matrix(nrow = n_separated_spes * length(AMD_pairs), ncol = 4))
colnames(separated_AMD_df) <- c("spe", "reference", "target", "AMD")


# Loop through each separated spes and get APD and AMD
setwd("~/Objects/separated_spes")
for (i in seq(n_separated_spes)) {
  
  # if (separated_spes_table[i, "shapeA"] != separated_spes_table[i, "shapeB"] ||
  #     separated_spes_table[i, "sizeA"] != separated_spes_table[i, "sizeB"] ||
  #     separated_spes_table[i, "arrangement"] != "S2") {
  #   i <- i + 1
  #   next
  # }
  
  
  # Read in current separated spe
  separated_spe_name <- paste("separated_spe_", i, sep = "")
  separated_spe_file_name <- paste(separated_spe_name, ".rds", sep = "")
  separated_spe <- readRDS(separated_spe_file_name)
  
  # pairwise_distance_data <- calculate_pairwise_distances_between_cell_types3D(separated_spe,
  #                                                                             cell_types,
  #                                                                             show_summary = F,
  #                                                                             plot_image = F)
  # pairwise_distance_data_summary <- summarise_distances_between_cell_types3D(pairwise_distance_data)

  # ## Fill in 3 rows at a time for APD df (as we have A/A, A/B, B/B)
  # index <- 3 * (i - 1) + 1 # index is 1, 4, 7, 10 ...
  # separated_APD_df[index:(index + 2), "spe"] <- separated_spe_name
  # separated_APD_df[index:(index + 2), "pair"] <- pairwise_distance_data_summary$pair
  # separated_APD_df[index:(index + 2), "APD"] <- pairwise_distance_data_summary$mean

  minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(separated_spe,
                                                                            cell_types,
                                                                            show_summary = F,
                                                                            plot_image = F)
  
  minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)
  ## Fill in 4 rows at a time for AMD df (as we have A/A, A/B, B/A, B/B)
  index <- 4 * (i - 1) + 1 # index is 1, 5, 9, 13
  separated_AMD_df[index:(index + 3), "spe"] <- separated_spe_name
  separated_AMD_df[index:(index + 3), "reference"] <- minimum_distance_data_summary$reference
  separated_AMD_df[index:(index + 3), "target"] <- minimum_distance_data_summary$target
  separated_AMD_df[index:(index + 3), "AMD"] <- minimum_distance_data_summary$mean
}

setwd("~/Objects/separated_spes/analysis_3D")
write.table(separated_APD_df, file = "separated_APD_df.csv")
write.table(separated_AMD_df, file = "separated_AMD_df.csv")

### 3.2. Separated spes - cc gradient based metrics ----------------------------------

# Get separated_spes_table
setwd("~/Objects/spes_table")
separated_spes_table <- read.table("separated_spes_table.csv")

# Get number of separated spes
n_separated_spes <- nrow(separated_spes_table)

# Define MS, NMS, ACIN, ACINP, CKR, AE data frames as well as constants
cell_types <- c("A", "B") # Use A as reference, and B as target, and vice versa
radii <- 50
radii_colnames <- paste("r", seq(radii), sep = "")

separated_MS_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types), ncol = 3 + radii))
colnames(separated_MS_df) <- c("spe", "reference", "target", radii_colnames)

separated_NMS_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types), ncol = 3 + radii))
colnames(separated_NMS_df) <- c("spe", "reference", "target", radii_colnames)

# Target is always A and B together
# Only choose prop(A) as prop(B) = 1 - prop(A) always
separated_ACINP_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types), ncol = 2 + radii))
colnames(separated_ACINP_df) <- c("spe", "reference", radii_colnames)

# Target is always A and B together
separated_AE_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types), ncol = 2 + radii))
colnames(separated_AE_df) <- c("spe", "reference", radii_colnames)


## ACIN and CKR are twice as large
# (ref A and tar A or B) OR (ref B and tar B or A)
separated_ACIN_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types)^2, ncol = 3 + radii))
colnames(separated_ACIN_df) <- c("spe", "reference", "target", radii_colnames)

# (ref A and tar A or B) OR (ref B and tar B or A)
separated_CKR_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types)^2, ncol = 3 + radii))
colnames(separated_CKR_df) <- c("spe", "reference", "target", radii_colnames)



# Loop through each separated spes and get gradient-based metrics
setwd("~/Objects/separated_spes")
for (i in seq(n_separated_spes)) {
  
  # Read in current separated spe
  separated_spe_name <- paste("separated_spe_", i, sep = "")
  separated_spe_file_name <- paste(separated_spe_name, ".rds", sep = "")
  separated_spe <- readRDS(separated_spe_file_name)
  
  
  index1 <- 2 * (i - 1) + 1 # index1 is 1, 3, 5, ...
  index2 <- 4 * (i - 1) + 1 # index2 is 1, 5, 9, 13...
  for (reference_cell_type in cell_types) {
    ## Calculate MS, NMS, ACINP and AE first
    target_cell_type <- setdiff(cell_types, reference_cell_type)
    
    # MS and NMS
    MS_NMS_data <- calculate_mixing_scores_gradient3D(separated_spe,
                                                      reference_cell_type,
                                                      target_cell_type,
                                                      radii,
                                                      plot_image = F)
    separated_MS_df[index1, c("spe", "reference", "target")] <- c(separated_spe_name, reference_cell_type, target_cell_type)
    separated_MS_df[index1, radii_colnames] <- MS_NMS_data$mixing_score
    
    separated_NMS_df[index1, c("spe", "reference", "target")] <- c(separated_spe_name, reference_cell_type, target_cell_type)
    separated_NMS_df[index1, radii_colnames] <- MS_NMS_data$normalised_mixing_score
    
    
    # ACINP
    ACINP_data <- calculate_cells_in_neighbourhood_proportions_gradient3D(separated_spe,
                                                                          reference_cell_type,
                                                                          cell_types, # Use both cell types, but keep prop for A
                                                                          radii,
                                                                          plot_image = F)
    separated_ACINP_df[index1, c("spe", "reference")] <- c(separated_spe_name, reference_cell_type)
    separated_ACINP_df[index1, radii_colnames] <- ACINP_data[["A"]]
    
    
    # AE
    AE_data <- calculate_entropy_gradient3D(separated_spe,
                                            reference_cell_type,
                                            cell_types, # Use both A and B
                                            radii,
                                            plot_image = F)
    separated_AE_df[index1, c("spe", "reference")] <- c(separated_spe_name, reference_cell_type)
    separated_AE_df[index1, radii_colnames] <- AE_data$entropy
    
    index1 <- index1 + 1
    
    
    for (target_cell_type in cell_types) {
      ## Calculate ACIN and CKR as target cell type can also be the reference cell type
      
      # ACIN
      ACIN_data <- calculate_cells_in_neighbourhood_gradient3D(separated_spe,
                                                               reference_cell_type,
                                                               target_cell_type,
                                                               radii,
                                                               plot_image = F)
      separated_ACIN_df[index2, c("spe", "reference", "target")] <- c(separated_spe_name, reference_cell_type, target_cell_type)
      separated_ACIN_df[index2, radii_colnames] <- ACIN_data[[target_cell_type]]
      
      
      # CKR
      CK_data <- calculate_cross_K_gradient3D(separated_spe,
                                              reference_cell_type,
                                              target_cell_type,
                                              radii,
                                              plot_image = F)
      separated_CKR_df[index2, c("spe", "reference", "target")] <- c(separated_spe_name, reference_cell_type, target_cell_type)
      separated_CKR_df[index2, radii_colnames] <- CK_data$observed_cross_K / CK_data$expected_cross_K
      
      index2 <- index2 + 1
    }
  }
}

setwd("~/Objects/separated_spes/analysis_3D")
write.table(separated_MS_df, file = "separated_MS_df.csv")
write.table(separated_NMS_df, file = "separated_NMS_df.csv")
write.table(separated_ACINP_df, file = "separated_ACINP_df.csv")
write.table(separated_AE_df, file = "separated_AE_df.csv")
write.table(separated_ACIN_df, file = "separated_ACIN_df.csv")
write.table(separated_CKR_df, file = "separated_CKR_df.csv")


### 3.3. Separated spes - heterogeneity metrics --------------------------------------
# Get separated_spes_table
setwd("~/Objects/spes_table")
separated_spes_table <- read.table("separated_spes_table.csv")

# Get number of separated spes
n_separated_spes <- nrow(separated_spes_table)

# Define SAC and prevalence data frames as well as constants
n_splits <- 10
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))

separated_prop_SAC_df_colnames <- c("spe", "reference", "target", "proportion")
separated_prop_SAC_df <- data.frame(matrix(nrow = n_separated_spes * nrow(prop_cell_types), ncol = length(separated_prop_SAC_df_colnames)))
colnames(separated_prop_SAC_df) <- separated_prop_SAC_df_colnames

separated_prop_prevalence_df_colnames <- c("spe", "reference", "target", thresholds_colnames)
separated_prop_prevalence_df <- data.frame(matrix(nrow = n_separated_spes * nrow(prop_cell_types), ncol = length(separated_prop_prevalence_df_colnames)))
colnames(separated_prop_prevalence_df) <- separated_prop_prevalence_df_colnames


entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))

separated_entropy_SAC_df_colnames <- c("spe", "cell_types", "entropy")
separated_entropy_SAC_df <- data.frame(matrix(nrow = n_separated_spes * nrow(entropy_cell_types), ncol = length(separated_entropy_SAC_df_colnames)))
colnames(separated_entropy_SAC_df) <- separated_entropy_SAC_df_colnames

separated_entropy_prevalence_df_colnames <- c("spe", "cell_types", thresholds_colnames)
separated_entropy_prevalence_df <- data.frame(matrix(nrow = n_separated_spes * nrow(entropy_cell_types), ncol = length(separated_entropy_prevalence_df_colnames)))
colnames(separated_entropy_prevalence_df) <- separated_entropy_prevalence_df_colnames


# Loop through each separated spes and get SAC and prevalence
setwd("~/Objects/separated_spes")
for (i in seq(n_separated_spes)) {
  
  # Read in current separated spe
  separated_spe_name <- paste("separated_spe_", i, sep = "")
  separated_spe_file_name <- paste(separated_spe_name, ".rds", sep = "")
  separated_spe <- readRDS(separated_spe_file_name)
  
  # Get proportion grid metrics
  for (j in seq_len(nrow(prop_cell_types))) {
    proportion_grid_metrics <- determine_cell_proportion_grid_metrics3D(separated_spe, 
                                                                        n_splits,
                                                                        strsplit(prop_cell_types$ref[j], ",")[[1]], 
                                                                        strsplit(prop_cell_types$tar[j], ",")[[1]],
                                                                        plot_image = F)
    
    proportion_SAC <- determine_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                          "proportion",
                                                          "binary")
    
    proportion_prevalence_df <- determine_prevalence_gradient3D(proportion_grid_metrics,
                                                                "proportion",
                                                                show_AUC = T,
                                                                plot_image = F)
    
    index <- nrow(prop_cell_types) * (i - 1) + j
    separated_prop_SAC_df[index, c("spe", "reference", "target")] <- c(separated_spe_name, prop_cell_types$ref[j], prop_cell_types$tar[j])
    separated_prop_SAC_df[index, "proportion"] <- proportion_SAC
    
    separated_prop_prevalence_df[index, c("spe", "reference", "target")] <- c(separated_spe_name, prop_cell_types$ref[j], prop_cell_types$tar[j])
    separated_prop_prevalence_df[index, thresholds_colnames] <- proportion_prevalence_df$prevalence
  }
  
  # Get entropy grid metrics
  for (j in seq_len(nrow(entropy_cell_types))) {
    entropy_grid_metrics <- determine_entropy_grid_metrics3D(separated_spe, 
                                                             n_splits,
                                                             strsplit(entropy_cell_types$cell_types[j], ",")[[1]], 
                                                             plot_image = F)
    
    entropy_SAC <- determine_spatial_autocorrelation3D(entropy_grid_metrics, 
                                                       "entropy",
                                                       "binary")
    
    entropy_prevalence_df <- determine_prevalence_gradient3D(entropy_grid_metrics,
                                                             "entropy",
                                                             show_AUC = T,
                                                             plot_image = F)
    
    index <- nrow(entropy_cell_types) * (i - 1) + j
    separated_entropy_SAC_df[index, c("spe", "cell_types")] <- c(separated_spe_name, entropy_cell_types$cell_types[j])
    separated_entropy_SAC_df[index, "entropy"] <- entropy_SAC
    
    separated_entropy_prevalence_df[index, c("spe", "cell_types")] <- c(separated_spe_name, entropy_cell_types$cell_types[j])
    separated_entropy_prevalence_df[index, thresholds_colnames] <- entropy_prevalence_df$prevalence
  }
}

setwd("~/Objects/separated_spes/analysis_3D")
write.table(separated_prop_SAC_df, file = "separated_prop_SAC_df.csv")
write.table(separated_prop_prevalence_df, file = "separated_prop_prevalence_df.csv")
write.table(separated_entropy_SAC_df, file = "separated_entropy_SAC_df.csv")
write.table(separated_entropy_prevalence_df, file = "separated_entropy_prevalence_df.csv")


### 4.0. Noisy background - Generating simulation -----------------------------

# I will also generate this simulation in this script since it is a special, and easy case
setwd("~/Objects/spes_metadata")
bg_spes_metadata <- readRDS("bg_spes_metadata_100_100_100.rds")
bg_spe_metadata_A_B <- bg_spes_metadata$AB
bg_spe_A_B <- simulate_spe_metadata3D(bg_spe_metadata_A_B)

### 4.1. Noisy background - cc distance metrics -----------------------------

setwd("~/Objects/background_spe")
bg_spe_A_B <- readRDS("bg_spe_A_B.rds")

n_bg_spes <- 1
bg_spe_name <- "bg_spe_A_B"

# APD and AMD
cell_types <- c("A", "B")
APD_pairs <- c("A/A", "A/B", "B/B")
AMD_pairs <- c("A/A", "A/B", "B/A", "B/B")

bg_APD_df <- data.frame(matrix(nrow = n_bg_spes * length(APD_pairs), ncol = 3))
colnames(bg_APD_df) <- c("spe", "pair", "APD")

bg_AMD_df <- data.frame(matrix(nrow = n_bg_spes * length(AMD_pairs), ncol = 4))
colnames(bg_AMD_df) <- c("spe", "reference", "target", "AMD")


pairwise_distance_data <- calculate_pairwise_distances_between_cell_types3D(bg_spe_A_B,
                                                                            cell_types,
                                                                            show_summary = F,
                                                                            plot_image = F)
pairwise_distance_data_summary <- summarise_distances_between_cell_types3D(pairwise_distance_data)

bg_APD_df[1:3, "spe"] <- bg_spe_name
bg_APD_df[1:3, "pair"] <- pairwise_distance_data_summary$pair
bg_APD_df[1:3, "APD"] <- pairwise_distance_data_summary$mean


minimum_distance_data <- calculate_minimum_distances_between_cell_types3D(bg_spe_A_B,
                                                                          cell_types,
                                                                          show_summary = F,
                                                                          plot_image = F)

minimum_distance_data_summary <- summarise_distances_between_cell_types3D(minimum_distance_data)

bg_AMD_df[1:4, "spe"] <- bg_spe_name
bg_AMD_df[1:4, "reference"] <- minimum_distance_data_summary$reference
bg_AMD_df[1:4, "target"] <- minimum_distance_data_summary$target
bg_AMD_df[1:4, "AMD"] <- minimum_distance_data_summary$mean

setwd("~/Objects/background_spe")
# write.table(bg_APD_df, "bg_APD_df.csv")
# write.table(bg_AMD_df, "bg_AMD_df.csv")

### 4.2. Noisy background - cc gradient based metrics -----------------------------

setwd("~/Objects/background_spe")
bg_spe_A_B <- readRDS("bg_spe_A_B.rds")

n_bg_spes <- 1
bg_spe_name <- "bg_spe_A_B"


# MS, NMS, ACIN, ACINP, CKR and AE
radii <- 50
radii_colnames <- paste("r", seq(radii), sep = "")

bg_MS_df <- data.frame(matrix(nrow = n_bg_spes * length(cell_types), ncol = 3 + radii))
colnames(bg_MS_df) <- c("spe", "reference", "target", radii_colnames)

bg_NMS_df <- data.frame(matrix(nrow = n_bg_spes * length(cell_types), ncol = 3 + radii))
colnames(bg_NMS_df) <- c("spe", "reference", "target", radii_colnames)

bg_ACINP_df <- data.frame(matrix(nrow = n_bg_spes * length(cell_types), ncol = 2 + radii))
colnames(bg_ACINP_df) <- c("spe", "reference", radii_colnames)

bg_AE_df <- data.frame(matrix(nrow = n_bg_spes * length(cell_types), ncol = 2 + radii))
colnames(bg_AE_df) <- c("spe", "reference", radii_colnames)

bg_ACIN_df <- data.frame(matrix(nrow = n_bg_spes * length(cell_types)^2, ncol = 3 + radii))
colnames(bg_ACIN_df) <- c("spe", "reference", "target", radii_colnames)

bg_CKR_df <- data.frame(matrix(nrow = n_bg_spes * length(cell_types)^2, ncol = 3 + radii))
colnames(bg_CKR_df) <- c("spe", "reference", "target", radii_colnames)


index1 <- 1
index2 <- 1
for (reference_cell_type in cell_types) {
  ## Calculate MS, NMS, ACINP and AE first
  target_cell_type <- setdiff(cell_types, reference_cell_type)
  
  # MS and NMS
  MS_NMS_data <- calculate_mixing_scores_gradient3D(bg_spe_A_B,
                                                    reference_cell_type,
                                                    target_cell_type,
                                                    radii,
                                                    plot_image = F)
  bg_MS_df[index1, c("spe", "reference", "target")] <- c(bg_spe_name, reference_cell_type, target_cell_type)
  bg_MS_df[index1, radii_colnames] <- MS_NMS_data$mixing_score
  
  bg_NMS_df[index1, c("spe", "reference", "target")] <- c(bg_spe_name, reference_cell_type, target_cell_type)
  bg_NMS_df[index1, radii_colnames] <- MS_NMS_data$normalised_mixing_score
  
  
  # ACINP
  ACINP_data <- calculate_cells_in_neighbourhood_proportions_gradient3D(bg_spe_A_B,
                                                                        reference_cell_type,
                                                                        cell_types, # Use both cell types, but keep prop for A
                                                                        radii,
                                                                        plot_image = F)
  bg_ACINP_df[index1, c("spe", "reference")] <- c(bg_spe_name, reference_cell_type)
  bg_ACINP_df[index1, radii_colnames] <- ACINP_data[["A"]]
  
  
  # AE
  AE_data <- calculate_entropy_gradient3D(bg_spe_A_B,
                                          reference_cell_type,
                                          cell_types, # Use both A and B
                                          radii,
                                          plot_image = F)
  bg_AE_df[index1, c("spe", "reference")] <- c(bg_spe_name, reference_cell_type)
  bg_AE_df[index1, radii_colnames] <- AE_data$entropy
  
  index1 <- index1 + 1
  
  
  for (target_cell_type in cell_types) {
    ## Calculate ACIN and CKR as target cell type can also be the reference cell type
    
    # ACIN
    ACIN_data <- calculate_cells_in_neighbourhood_gradient3D(bg_spe_A_B,
                                                             reference_cell_type,
                                                             target_cell_type,
                                                             radii,
                                                             plot_image = F)
    bg_ACIN_df[index2, c("spe", "reference", "target")] <- c(bg_spe_name, reference_cell_type, target_cell_type)
    bg_ACIN_df[index2, radii_colnames] <- ACIN_data[[target_cell_type]]
    
    
    # CKR
    CK_data <- calculate_cross_K_gradient3D(bg_spe_A_B,
                                            reference_cell_type,
                                            target_cell_type,
                                            radii,
                                            plot_image = F)
    bg_CKR_df[index2, c("spe", "reference", "target")] <- c(bg_spe_name, reference_cell_type, target_cell_type)
    bg_CKR_df[index2, radii_colnames] <- CK_data$observed_cross_K / CK_data$expected_cross_K
    
    index2 <- index2 + 1
  }
}

setwd("~/Objects/background_spe")
# write.table(bg_MS_df, "bg_MS_df.csv")
# write.table(bg_NMS_df, "bg_NMS_df.csv")
# write.table(bg_ACIN_df, "bg_ACIN_df.csv")
# write.table(bg_ACINP_df, "bg_ACINP_df.csv")
# write.table(bg_CKR_df, "bg_CKR_df.csv")
# write.table(bg_AE_df, "bg_AE_df.csv")



### 4.3. Noisy background - heterogeneity metrics -----------------------------

setwd("~/Objects/background_spe")
bg_spe_A_B <- readRDS("bg_spe_A_B.rds")

n_bg_spes <- 1
bg_spe_name <- "bg_spe_A_B"

# SAC and prevalence
n_splits <- 10
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))

bg_prop_SAC_df_colnames <- c("spe", "reference", "target", "proportion")
bg_prop_SAC_df <- data.frame(matrix(nrow = n_bg_spes * nrow(prop_cell_types), ncol = length(bg_prop_SAC_df_colnames)))
colnames(bg_prop_SAC_df) <- bg_prop_SAC_df_colnames

bg_prop_prevalence_df_colnames <- c("spe", "reference", "target", thresholds_colnames)
bg_prop_prevalence_df <- data.frame(matrix(nrow = n_bg_spes * nrow(prop_cell_types), ncol = length(bg_prop_prevalence_df_colnames)))
colnames(bg_prop_prevalence_df) <- bg_prop_prevalence_df_colnames

entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))

bg_entropy_SAC_df_colnames <- c("spe", "cell_types", "entropy")
bg_entropy_SAC_df <- data.frame(matrix(nrow = n_bg_spes * nrow(entropy_cell_types), ncol = length(bg_entropy_SAC_df_colnames)))
colnames(bg_entropy_SAC_df) <- bg_entropy_SAC_df_colnames

bg_entropy_prevalence_df_colnames <- c("spe", "cell_types", thresholds_colnames)
bg_entropy_prevalence_df <- data.frame(matrix(nrow = n_bg_spes * nrow(entropy_cell_types), ncol = length(bg_entropy_prevalence_df_colnames)))
colnames(bg_entropy_prevalence_df) <- bg_entropy_prevalence_df_colnames


# Get proportion grid metrics
for (i in seq_len(nrow(prop_cell_types))) {
  proportion_grid_metrics <- determine_cell_proportion_grid_metrics3D(bg_spe_A_B, 
                                                                      n_splits,
                                                                      strsplit(prop_cell_types$ref[i], ",")[[1]], 
                                                                      strsplit(prop_cell_types$tar[i], ",")[[1]],
                                                                      plot_image = F)
  
  proportion_SAC <- determine_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                        "proportion",
                                                        "binary")
  
  proportion_prevalence_df <- determine_prevalence_gradient3D(proportion_grid_metrics,
                                                              "proportion",
                                                              show_AUC = T,
                                                              plot_image = F)
  
  bg_prop_SAC_df[i, c("spe", "reference", "target")] <- c(bg_spe_name, prop_cell_types$ref[i], prop_cell_types$tar[i])
  bg_prop_SAC_df[i, "proportion"] <- proportion_SAC
  
  bg_prop_prevalence_df[i, c("spe", "reference", "target")] <- c(bg_spe_name, prop_cell_types$ref[i], prop_cell_types$tar[i])
  bg_prop_prevalence_df[i, thresholds_colnames] <- proportion_prevalence_df$prevalence
}

# Get entropy grid metrics
for (i in seq_len(nrow(entropy_cell_types))) {
  entropy_grid_metrics <- determine_entropy_grid_metrics3D(bg_spe_A_B, 
                                                           n_splits,
                                                           strsplit(entropy_cell_types$cell_types[i], ",")[[1]], 
                                                           plot_image = F)
  
  entropy_SAC <- determine_spatial_autocorrelation3D(entropy_grid_metrics, 
                                                     "entropy",
                                                     "binary")
  
  entropy_prevalence_df <- determine_prevalence_gradient3D(entropy_grid_metrics,
                                                           "entropy",
                                                           show_AUC = T,
                                                           plot_image = F)
  
  bg_entropy_SAC_df[i, c("spe", "cell_types")] <- c(bg_spe_name, entropy_cell_types$cell_types[i])
  bg_entropy_SAC_df[i, "entropy"] <- entropy_SAC
  
  bg_entropy_prevalence_df[i, c("spe", "cell_types")] <- c(bg_spe_name, entropy_cell_types$cell_types[i])
  bg_entropy_prevalence_df[i, thresholds_colnames] <- entropy_prevalence_df$prevalence
}

setwd("~/Objects/background_spe")
# write.table(bg_prop_SAC_df, "bg_prop_SAC_df.csv")
# write.table(bg_prop_prevalence_df, "bg_prop_prevalence_df.csv")
# write.table(bg_entropy_SAC_df, "bg_entropy_SAC_df.csv")
# write.table(bg_entropy_prevalence_df, "bg_entropy_prevalence_df.csv")
