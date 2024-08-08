### Mixed spes analysis --------------------------------------------

setwd("~/Objects/supervised/spes_metadata")
mixed_spes_metadata <- readRDS("mixed_spes_metadata_supervised.rds")

# Get number of mixed spes
n_mixed_spes <- length(mixed_spes_metadata)



# Define APD and AMD data frames as well as constants
cell_types <- c("A", "B")
APD_pairs <- c("A/A", "A/B", "B/B")
AMD_pairs <- c("A/A", "A/B", "B/A", "B/B")

mixed_APD_df <- data.frame(matrix(nrow = n_mixed_spes * length(APD_pairs), ncol = 3))
colnames(mixed_APD_df) <- c("spe", "pair", "APD")

mixed_AMD_df <- data.frame(matrix(nrow = n_mixed_spes * length(AMD_pairs), ncol = 4))
colnames(mixed_AMD_df) <- c("spe", "reference", "target", "AMD")


# Define MS, NMS, ACIN, ACINP, CKR, AE data frames as well as constants
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


for (i in seq_len(n_mixed_spes)) {
  
  mixed_spe <- simulate_spe_metadata3D(mixed_spes_metadata[[i]], plot_image = F)
  mixed_spe_name <- paste("mixed_spe_", i, sep = "")  
  
  
  
  
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
  
  
  
  
  
  index1 <- 2 * (i - 1) + 1 # index1 is 1, 3, 5, ...
  index2 <- 4 * (i - 1) + 1 # index2 is 1, 5, 9, 13...
  for (reference_cell_type in cell_types) {
    gradient_data <- calculate_all_gradient_cc_metrics3D(mixed_spe,
                                                         reference_cell_type,
                                                         cell_types,
                                                         radii,
                                                         plot_image = F)
    
    target_cell_type <- setdiff(cell_types, reference_cell_type)
    
    mixed_MS_df[index1, c("spe", "reference", "target")] <- c(mixed_spe_name, reference_cell_type, target_cell_type)
    mixed_MS_df[index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
    
    mixed_NMS_df[index1, c("spe", "reference", "target")] <- c(mixed_spe_name, reference_cell_type, target_cell_type)
    mixed_NMS_df[index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
    
    mixed_ACINP_df[index1, c("spe", "reference")] <- c(mixed_spe_name, reference_cell_type)
    mixed_ACINP_df[index1, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][["A"]]
    
    mixed_AE_df[index1, c("spe", "reference")] <- c(mixed_spe_name, reference_cell_type)
    mixed_AE_df[index1, radii_colnames] <- gradient_data[["entropy"]]$entropy
    
    index1 <- index1 + 1
    
    for (target_cell_type in cell_types) {
      ## Calculate ACIN and CKR as target cell type can also be the reference cell type
      
      # ACIN
      mixed_ACIN_df[index2, c("spe", "reference", "target")] <- c(mixed_spe_name, reference_cell_type, target_cell_type)
      mixed_ACIN_df[index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
      
      
      # CKR
      mixed_CKR_df[index2, c("spe", "reference", "target")] <- c(mixed_spe_name, reference_cell_type, target_cell_type)
      mixed_CKR_df[index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]]$cross_K_ratio
      
      index2 <- index2 + 1
    }
  }


  # Get proportion grid metrics
  for (j in seq_len(nrow(prop_cell_types))) {
    proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(mixed_spe, 
                                                                        n_splits,
                                                                        strsplit(prop_cell_types$ref[j], ",")[[1]], 
                                                                        strsplit(prop_cell_types$tar[j], ",")[[1]],
                                                                        plot_image = F)
    
    proportion_SAC <- calculate_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                          "proportion",
                                                          "binary")
    
    proportion_prevalence_df <- calculate_prevalence_gradient3D(proportion_grid_metrics,
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
    entropy_grid_metrics <- calculate_entropy_grid_metrics3D(mixed_spe, 
                                                             n_splits,
                                                             strsplit(entropy_cell_types$cell_types[j], ",")[[1]], 
                                                             plot_image = F)
    
    entropy_SAC <- calculate_spatial_autocorrelation3D(entropy_grid_metrics, 
                                                       "entropy",
                                                       "binary")
    
    entropy_prevalence_df <- calculate_prevalence_gradient3D(entropy_grid_metrics,
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


setwd("~/Objects/supervised/mixed_spes/analysis_3D")
write.table(mixed_APD_df, file = "mixed_APD_df.csv")
write.table(mixed_AMD_df, file = "mixed_AMD_df.csv")

write.table(mixed_MS_df, file = "mixed_MS_df.csv")
write.table(mixed_NMS_df, file = "mixed_NMS_df.csv")
write.table(mixed_ACINP_df, file = "mixed_ACINP_df.csv")
write.table(mixed_AE_df, file = "mixed_AE_df.csv")
write.table(mixed_ACIN_df, file = "mixed_ACIN_df.csv")
write.table(mixed_CKR_df, file = "mixed_CKR_df.csv")

write.table(mixed_prop_SAC_df, file = "mixed_prop_SAC_df.csv")
write.table(mixed_prop_prevalence_df, file = "mixed_prop_prevalence_df.csv")
write.table(mixed_entropy_SAC_df, file = "mixed_entropy_SAC_df.csv")
write.table(mixed_entropy_prevalence_df, file = "mixed_entropy_prevalence_df.csv")



### Ringed spes analysis --------------------------------------------

setwd("~/Objects/supervised/spes_metadata")
ringed_spes_metadata <- readRDS("ringed_spes_metadata_supervised.rds")

# Get number of ringed spes
n_ringed_spes <- length(ringed_spes_metadata)



# Define APD and AMD data frames as well as constants
cell_types <- c("A", "B")
APD_pairs <- c("A/A", "A/B", "B/B")
AMD_pairs <- c("A/A", "A/B", "B/A", "B/B")

ringed_APD_df <- data.frame(matrix(nrow = n_ringed_spes * length(APD_pairs), ncol = 3))
colnames(ringed_APD_df) <- c("spe", "pair", "APD")

ringed_AMD_df <- data.frame(matrix(nrow = n_ringed_spes * length(AMD_pairs), ncol = 4))
colnames(ringed_AMD_df) <- c("spe", "reference", "target", "AMD")


# Define MS, NMS, ACIN, ACINP, CKR, AE data frames as well as constants
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


for (i in seq_len(n_ringed_spes)) {
  
  ringed_spe <- simulate_spe_metadata3D(ringed_spes_metadata[[i]], plot_image = F)
  ringed_spe_name <- paste("ringed_spe_", i, sep = "")  
  
  
  
  
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
  
  
  
  
  
  index1 <- 2 * (i - 1) + 1 # index1 is 1, 3, 5, ...
  index2 <- 4 * (i - 1) + 1 # index2 is 1, 5, 9, 13...
  for (reference_cell_type in cell_types) {
    gradient_data <- calculate_all_gradient_cc_metrics3D(ringed_spe,
                                                         reference_cell_type,
                                                         cell_types,
                                                         radii,
                                                         plot_image = F)
    
    target_cell_type <- setdiff(cell_types, reference_cell_type)
    
    ringed_MS_df[index1, c("spe", "reference", "target")] <- c(ringed_spe_name, reference_cell_type, target_cell_type)
    ringed_MS_df[index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
    
    ringed_NMS_df[index1, c("spe", "reference", "target")] <- c(ringed_spe_name, reference_cell_type, target_cell_type)
    ringed_NMS_df[index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
    
    ringed_ACINP_df[index1, c("spe", "reference")] <- c(ringed_spe_name, reference_cell_type)
    ringed_ACINP_df[index1, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][["A"]]
    
    ringed_AE_df[index1, c("spe", "reference")] <- c(ringed_spe_name, reference_cell_type)
    ringed_AE_df[index1, radii_colnames] <- gradient_data[["entropy"]]$entropy
    
    index1 <- index1 + 1
    
    for (target_cell_type in cell_types) {
      ## Calculate ACIN and CKR as target cell type can also be the reference cell type
      
      # ACIN
      ringed_ACIN_df[index2, c("spe", "reference", "target")] <- c(ringed_spe_name, reference_cell_type, target_cell_type)
      ringed_ACIN_df[index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
      
      
      # CKR
      ringed_CKR_df[index2, c("spe", "reference", "target")] <- c(ringed_spe_name, reference_cell_type, target_cell_type)
      ringed_CKR_df[index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]]$cross_K_ratio
      
      index2 <- index2 + 1
    }
  }
  
  
  # Get proportion grid metrics
  for (j in seq_len(nrow(prop_cell_types))) {
    proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(ringed_spe, 
                                                                        n_splits,
                                                                        strsplit(prop_cell_types$ref[j], ",")[[1]], 
                                                                        strsplit(prop_cell_types$tar[j], ",")[[1]],
                                                                        plot_image = F)
    
    proportion_SAC <- calculate_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                          "proportion",
                                                          "binary")
    
    proportion_prevalence_df <- calculate_prevalence_gradient3D(proportion_grid_metrics,
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
    entropy_grid_metrics <- calculate_entropy_grid_metrics3D(ringed_spe, 
                                                             n_splits,
                                                             strsplit(entropy_cell_types$cell_types[j], ",")[[1]], 
                                                             plot_image = F)
    
    entropy_SAC <- calculate_spatial_autocorrelation3D(entropy_grid_metrics, 
                                                       "entropy",
                                                       "binary")
    
    entropy_prevalence_df <- calculate_prevalence_gradient3D(entropy_grid_metrics,
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


setwd("~/Objects/supervised/ringed_spes/analysis_3D")
write.table(ringed_APD_df, file = "ringed_APD_df.csv")
write.table(ringed_AMD_df, file = "ringed_AMD_df.csv")

write.table(ringed_MS_df, file = "ringed_MS_df.csv")
write.table(ringed_NMS_df, file = "ringed_NMS_df.csv")
write.table(ringed_ACINP_df, file = "ringed_ACINP_df.csv")
write.table(ringed_AE_df, file = "ringed_AE_df.csv")
write.table(ringed_ACIN_df, file = "ringed_ACIN_df.csv")
write.table(ringed_CKR_df, file = "ringed_CKR_df.csv")

write.table(ringed_prop_SAC_df, file = "ringed_prop_SAC_df.csv")
write.table(ringed_prop_prevalence_df, file = "ringed_prop_prevalence_df.csv")
write.table(ringed_entropy_SAC_df, file = "ringed_entropy_SAC_df.csv")
write.table(ringed_entropy_prevalence_df, file = "ringed_entropy_prevalence_df.csv")
### Separated spes analysis --------------------------------------------

setwd("~/Objects/supervised/spes_metadata")
separated_spes_metadata <- readRDS("separated_spes_metadata_supervised.rds")

# Get number of separated spes
n_separated_spes <- length(separated_spes_metadata)



# Define APD and AMD data frames as well as constants
cell_types <- c("A", "B")
APD_pairs <- c("A/A", "A/B", "B/B")
AMD_pairs <- c("A/A", "A/B", "B/A", "B/B")

separated_APD_df <- data.frame(matrix(nrow = n_separated_spes * length(APD_pairs), ncol = 3))
colnames(separated_APD_df) <- c("spe", "pair", "APD")

separated_AMD_df <- data.frame(matrix(nrow = n_separated_spes * length(AMD_pairs), ncol = 4))
colnames(separated_AMD_df) <- c("spe", "reference", "target", "AMD")


# Define MS, NMS, ACIN, ACINP, CKR, AE data frames as well as constants
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


for (i in seq_len(n_separated_spes)) {
  
  separated_spe <- simulate_spe_metadata3D(separated_spes_metadata[[i]], plot_image = F)
  separated_spe_name <- paste("separated_spe_", i, sep = "")  
  
  
  
  
  # pairwise_distance_data <- calculate_pairwise_distances_between_cell_types3D(separated_spe,
  #                                                                             cell_types,
  #                                                                             show_summary = F,
  #                                                                             plot_image = F)
  # pairwise_distance_data_summary <- summarise_distances_between_cell_types3D(pairwise_distance_data)
  # 
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
  
  
  
  
  
  index1 <- 2 * (i - 1) + 1 # index1 is 1, 3, 5, ...
  index2 <- 4 * (i - 1) + 1 # index2 is 1, 5, 9, 13...
  for (reference_cell_type in cell_types) {
    gradient_data <- calculate_all_gradient_cc_metrics3D(separated_spe,
                                                         reference_cell_type,
                                                         cell_types,
                                                         radii,
                                                         plot_image = F)
    
    target_cell_type <- setdiff(cell_types, reference_cell_type)
    
    separated_MS_df[index1, c("spe", "reference", "target")] <- c(separated_spe_name, reference_cell_type, target_cell_type)
    separated_MS_df[index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
    
    separated_NMS_df[index1, c("spe", "reference", "target")] <- c(separated_spe_name, reference_cell_type, target_cell_type)
    separated_NMS_df[index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
    
    separated_ACINP_df[index1, c("spe", "reference")] <- c(separated_spe_name, reference_cell_type)
    separated_ACINP_df[index1, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][["A"]]
    
    separated_AE_df[index1, c("spe", "reference")] <- c(separated_spe_name, reference_cell_type)
    separated_AE_df[index1, radii_colnames] <- gradient_data[["entropy"]]$entropy
    
    index1 <- index1 + 1
    
    for (target_cell_type in cell_types) {
      ## Calculate ACIN and CKR as target cell type can also be the reference cell type
      
      # ACIN
      separated_ACIN_df[index2, c("spe", "reference", "target")] <- c(separated_spe_name, reference_cell_type, target_cell_type)
      separated_ACIN_df[index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
      
      
      # CKR
      separated_CKR_df[index2, c("spe", "reference", "target")] <- c(separated_spe_name, reference_cell_type, target_cell_type)
      separated_CKR_df[index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]]$cross_K_ratio
      
      index2 <- index2 + 1
    }
  }
  
  
  # Get proportion grid metrics
  for (j in seq_len(nrow(prop_cell_types))) {
    proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(separated_spe, 
                                                                        n_splits,
                                                                        strsplit(prop_cell_types$ref[j], ",")[[1]], 
                                                                        strsplit(prop_cell_types$tar[j], ",")[[1]],
                                                                        plot_image = F)
    
    proportion_SAC <- calculate_spatial_autocorrelation3D(proportion_grid_metrics, 
                                                          "proportion",
                                                          "binary")
    
    proportion_prevalence_df <- calculate_prevalence_gradient3D(proportion_grid_metrics,
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
    entropy_grid_metrics <- calculate_entropy_grid_metrics3D(separated_spe, 
                                                             n_splits,
                                                             strsplit(entropy_cell_types$cell_types[j], ",")[[1]], 
                                                             plot_image = F)
    
    entropy_SAC <- calculate_spatial_autocorrelation3D(entropy_grid_metrics, 
                                                       "entropy",
                                                       "binary")
    
    entropy_prevalence_df <- calculate_prevalence_gradient3D(entropy_grid_metrics,
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


setwd("~/Objects/supervised/separated_spes/analysis_3D")
write.table(separated_APD_df, file = "separated_APD_df.csv")
write.table(separated_AMD_df, file = "separated_AMD_df.csv")

write.table(separated_MS_df, file = "separated_MS_df.csv")
write.table(separated_NMS_df, file = "separated_NMS_df.csv")
write.table(separated_ACINP_df, file = "separated_ACINP_df.csv")
write.table(separated_AE_df, file = "separated_AE_df.csv")
write.table(separated_ACIN_df, file = "separated_ACIN_df.csv")
write.table(separated_CKR_df, file = "separated_CKR_df.csv")

write.table(separated_prop_SAC_df, file = "separated_prop_SAC_df.csv")
write.table(separated_prop_prevalence_df, file = "separated_prop_prevalence_df.csv")
write.table(separated_entropy_SAC_df, file = "separated_entropy_SAC_df.csv")
write.table(ringed_entropy_prevalence_df, file = "ringed_entropy_prevalence_df.csv")