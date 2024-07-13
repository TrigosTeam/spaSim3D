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

  # # Read in current mixed spe
  # mixed_spe_name <- paste("mixed_spe_", i, sep = "")
  # mixed_spe_file_name <- paste(mixed_spe_name, ".rds", sep = "")
  # mixed_spe <- readRDS(mixed_spe_file_name)
  # 
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
cell_types <- c("A", "B")
n_splits <- 10

mixed_SAC_df <- data.frame(matrix(nrow = n_mixed_spes, ncol = 3))
colnames(mixed_SAC_df) <- c("spe", "proportion", "entropy")

thresholds <- seq(0, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")
mixed_prevalence_df <- data.frame(matrix(nrow = n_mixed_spes, ncol = 2 + length(thresholds)))
colnames(mixed_prevalence_df) <- c("spe", "metric", thresholds_colnames)


# Loop through each mixed spes and get SAC and prevalence
setwd("~/Objects/mixed_spes")
for (i in seq(n_mixed_spes)) {
  
  # Read in current mixed spe
  mixed_spe_name <- paste("mixed_spe_", i, sep = "")
  mixed_spe_file_name <- paste(mixed_spe_name, ".rds", sep = "")
  mixed_spe <- readRDS(mixed_spe_file_name)
  
  # Get grid metrics
  proportion_grid_metrics <- determine_cell_proportion_grid_metrics3D(mixed_spe, 
                                                                      n_splits,
                                                                      cell_types[2], 
                                                                      cell_types[1], # Assume A is target, but doesn't matter
                                                                      plot_image = F)
  entropy_grid_metrics <- determine_entropy_grid_metrics3D(mixed_spe,
                                                           n_splits,
                                                           cell_types,
                                                           plot_image = F)
  
  # Calculate SACs
  proportion_SAC <- determine_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                        "proportion",
                                                        "binary")
  entropy_SAC <- determine_spatial_autocorrelation3D(entropy_grid_metrics,
                                                     "entropy",
                                                     "binary")
  
  
  # Calculate prevalence gradients
  proportion_prevalence_data <- determine_prevalence_gradient3D(proportion_grid_metrics,
                                                                "proportion",
                                                                plot_image = F)
  entropy_prevalence_data <- determine_prevalence_gradient3D(entropy_grid_metrics,
                                                             "entropy",
                                                             plot_image = F)
  

  # Add SACs to mixed_SAC_df
  mixed_SAC_df[i, "spe"] <- mixed_spe_name
  mixed_SAC_df[i, "proportion"] <- proportion_SAC
  mixed_SAC_df[i, "entropy"] <- entropy_SAC
  
  
  ## Add prevalences to mixed_prevalence_df
  
  # Fill in 2 rows at a time (proportion and entropy)
  index <- 2 * (i - 1) + 1 # index is 1, 3, 5, 7 ...
  mixed_prevalence_df[index:(index + 1), "spe"] <- mixed_spe_name
  mixed_prevalence_df[index:(index + 1), "metric"] <- c("proportion", "entropy")
  
  mixed_prevalence_df[index, thresholds_colnames] <- proportion_prevalence_data$prevalence
  mixed_prevalence_df[index + 1, thresholds_colnames] <- entropy_prevalence_data$prevalence
}

setwd("~/Objects/mixed_spes/analysis_3D")
write.table(mixed_SAC_df, file = "mixed_SAC_df.csv")
write.table(mixed_prevalence_df, file = "mixed_prevalence_df.csv")


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
cell_types <- c("A", "B")
n_splits <- 10

ringed_SAC_df <- data.frame(matrix(nrow = n_ringed_spes, ncol = 3))
colnames(ringed_SAC_df) <- c("spe", "proportion", "entropy")

thresholds <- seq(0, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")
ringed_prevalence_df <- data.frame(matrix(nrow = n_ringed_spes, ncol = 2 + length(thresholds)))
colnames(ringed_prevalence_df) <- c("spe", "metric", thresholds_colnames)


# Loop through each ringed spes and get SAC and prevalence
setwd("~/Objects/ringed_spes")
for (i in seq(n_ringed_spes)) {
  
  # Read in current ringed spe
  ringed_spe_name <- paste("ringed_spe_", i, sep = "")
  ringed_spe_file_name <- paste(ringed_spe_name, ".rds", sep = "")
  ringed_spe <- readRDS(ringed_spe_file_name)
  
  # Get grid metrics
  proportion_grid_metrics <- determine_cell_proportion_grid_metrics3D(ringed_spe, 
                                                                      n_splits,
                                                                      cell_types[2], 
                                                                      cell_types[1], # Assume A is target, but doesn't matter
                                                                      plot_image = F)
  entropy_grid_metrics <- determine_entropy_grid_metrics3D(ringed_spe,
                                                           n_splits,
                                                           cell_types,
                                                           plot_image = F)
  
  # Calculate SACs
  proportion_SAC <- determine_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                        "proportion",
                                                        "binary")
  entropy_SAC <- determine_spatial_autocorrelation3D(entropy_grid_metrics,
                                                     "entropy",
                                                     "binary")
  
  
  # Calculate prevalence gradients
  proportion_prevalence_data <- determine_prevalence_gradient3D(proportion_grid_metrics,
                                                                "proportion",
                                                                plot_image = F)
  entropy_prevalence_data <- determine_prevalence_gradient3D(entropy_grid_metrics,
                                                             "entropy",
                                                             plot_image = F)
  
  
  # Add SACs to ringed_SAC_df
  ringed_SAC_df[i, "spe"] <- ringed_spe_name
  ringed_SAC_df[i, "proportion"] <- proportion_SAC
  ringed_SAC_df[i, "entropy"] <- entropy_SAC
  
  
  ## Add prevalences to ringed_prevalence_df
  
  # Fill in 2 rows at a time (proportion and entropy)
  index <- 2 * (i - 1) + 1 # index is 1, 3, 5, 7 ...
  ringed_prevalence_df[index:(index + 1), "spe"] <- ringed_spe_name
  ringed_prevalence_df[index:(index + 1), "metric"] <- c("proportion", "entropy")
  
  ringed_prevalence_df[index, thresholds_colnames] <- proportion_prevalence_data$prevalence
  ringed_prevalence_df[index + 1, thresholds_colnames] <- entropy_prevalence_data$prevalence
}

setwd("~/Objects/ringed_spes/analysis_3D")
write.table(ringed_SAC_df, file = "ringed_SAC_df.csv")
write.table(ringed_prevalence_df, file = "ringed_prevalence_df.csv")


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
  
  # Read in current separated spe
  separated_spe_name <- paste("separated_spe_", i, sep = "")
  separated_spe_file_name <- paste(separated_spe_name, ".rds", sep = "")
  separated_spe <- readRDS(separated_spe_file_name)
  
  pairwise_distance_data <- calculate_pairwise_distances_between_cell_types3D(separated_spe,
                                                                              cell_types,
                                                                              show_summary = F,
                                                                              plot_image = F)
  pairwise_distance_data_summary <- summarise_distances_between_cell_types3D(pairwise_distance_data)
  
  ## Fill in 3 rows at a time for APD df (as we have A/A, A/B, B/B)
  index <- 3 * (i - 1) + 1 # index is 1, 4, 7, 10 ...
  separated_APD_df[index:(index + 2), "spe"] <- separated_spe_name
  separated_APD_df[index:(index + 2), "pair"] <- pairwise_distance_data_summary$pair
  separated_APD_df[index:(index + 2), "APD"] <- pairwise_distance_data_summary$mean
  
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
cell_types <- c("A", "B")
n_splits <- 10

separated_SAC_df <- data.frame(matrix(nrow = n_separated_spes, ncol = 3))
colnames(separated_SAC_df) <- c("spe", "proportion", "entropy")

thresholds <- seq(0, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")
separated_prevalence_df <- data.frame(matrix(nrow = n_separated_spes, ncol = 2 + length(thresholds)))
colnames(separated_prevalence_df) <- c("spe", "metric", thresholds_colnames)


# Loop through each separated spes and get SAC and prevalence
setwd("~/Objects/separated_spes")
for (i in seq(n_separated_spes)) {
  
  # Read in current separated spe
  separated_spe_name <- paste("separated_spe_", i, sep = "")
  separated_spe_file_name <- paste(separated_spe_name, ".rds", sep = "")
  separated_spe <- readRDS(separated_spe_file_name)
  
  # Get grid metrics
  proportion_grid_metrics <- determine_cell_proportion_grid_metrics3D(separated_spe, 
                                                                      n_splits,
                                                                      cell_types[2], 
                                                                      cell_types[1], # Assume A is target, but doesn't matter
                                                                      plot_image = F)
  entropy_grid_metrics <- determine_entropy_grid_metrics3D(separated_spe,
                                                           n_splits,
                                                           cell_types,
                                                           plot_image = F)
  
  # Calculate SACs
  proportion_SAC <- determine_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                        "proportion",
                                                        "binary")
  entropy_SAC <- determine_spatial_autocorrelation3D(entropy_grid_metrics,
                                                     "entropy",
                                                     "binary")
  
  
  # Calculate prevalence gradients
  proportion_prevalence_data <- determine_prevalence_gradient3D(proportion_grid_metrics,
                                                                "proportion",
                                                                plot_image = F)
  entropy_prevalence_data <- determine_prevalence_gradient3D(entropy_grid_metrics,
                                                             "entropy",
                                                             plot_image = F)
  
  
  # Add SACs to separated_SAC_df
  separated_SAC_df[i, "spe"] <- separated_spe_name
  separated_SAC_df[i, "proportion"] <- proportion_SAC
  separated_SAC_df[i, "entropy"] <- entropy_SAC
  
  
  ## Add prevalences to separated_prevalence_df
  
  # Fill in 2 rows at a time (proportion and entropy)
  index <- 2 * (i - 1) + 1 # index is 1, 3, 5, 7 ...
  separated_prevalence_df[index:(index + 1), "spe"] <- separated_spe_name
  separated_prevalence_df[index:(index + 1), "metric"] <- c("proportion", "entropy")
  
  separated_prevalence_df[index, thresholds_colnames] <- proportion_prevalence_data$prevalence
  separated_prevalence_df[index + 1, thresholds_colnames] <- entropy_prevalence_data$prevalence
}

setwd("~/Objects/separated_spes/analysis_3D")
write.table(separated_SAC_df, file = "separated_SAC_df.csv")
write.table(separated_prevalence_df, file = "separated_prevalence_df.csv")


### Spacer --------------------------------

