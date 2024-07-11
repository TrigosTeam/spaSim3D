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
  
  pairwise_distance_data <- calculate_pairwise_distances_between_cell_types3D(mixed_spe,
                                                                              cell_types,
                                                                              show_summary = F,
                                                                              plot_image = F)
  pairwise_distance_data_summary <- summarise_distances_between_cell_types3D(pairwise_distance_data)
  
  ## Fill in 3 rows at a time for APD df (as we have A/A, A/B, B/B)
  index <- 3 * (i - 1) + 1 # index is 1, 4, 7, 10 ...
  mixed_APD_df[index:(index + 2), "spe"] <- mixed_spe_name
  mixed_APD_df[index:(index + 2), "pair"] <- pairwise_distance_data_summary$pair
  mixed_APD_df[index:(index + 2), "APD"] <- pairwise_distance_data_summary$mean
  
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
  
  if (i == 5) break
  
}

setwd("~/Objects/mixed_spes/analysis_3D")
write.table(mixed_MS_df, file = "mixed_MS_df.csv")
write.table(mixed_NMS_df, file = "mixed_NMS_df.csv")
write.table(mixed_ACINP_df, file = "mixed_ACINP_df.csv")
write.table(mixed_AE_df, file = "mixed_AE_df.csv")
write.table(mixed_ACIN_df, file = "mixed_ACIN_df.csv")
write.table(mixed_CKR_df, file = "mixed_CKR_df.csv")


### 1.3. Mixed spes - heterogeneity metrics --------------------------------------



### Spacer --------------------------------

