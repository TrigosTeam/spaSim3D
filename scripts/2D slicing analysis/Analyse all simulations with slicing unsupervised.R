### Slicing function -------------------------------------------------

# Function to get slices from spe
get_spe_slices_list <- function(spe) {
  
  spe_slices_list <- list()
  
  slice_z_coords <- seq(150, 450, 60)
  
  for (i in seq(length(slice_z_coords) - 1)) {
    bottom_z_coord <- slice_z_coords[i]
    top_z_coord <- slice_z_coords[i + 1]
    spe_z_coords <- spatialCoords(spe)[ , "Cell.Z.Position"]
    spe_slice <- spe[, bottom_z_coord < spe_z_coords & spe_z_coords < top_z_coord]
    spatialCoords(spe_slice) <- spatialCoords(spe_slice)[ , c("Cell.X.Position", "Cell.Y.Position")]
    
    spe_slices_list[[i]] <- spe_slice
  }
  return(spe_slices_list)
}


### Mixed spes 3D setup --------------------------------------------

setwd("~/Objects/unsupervised/spes_metadata")
mixed_spes_metadata <- readRDS("mixed_spes_metadata_unsupervised.rds")

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
radii <- seq(10, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

mixed_MS_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types), ncol = 3 + length(radii)))
colnames(mixed_MS_df) <- c("spe", "reference", "target", radii_colnames)

mixed_NMS_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types), ncol = 3 + length(radii)))
colnames(mixed_NMS_df) <- c("spe", "reference", "target", radii_colnames)

# Target is always A and B together
# Only choose prop(A) as prop(B) = 1 - prop(A) always
mixed_ACINP_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types), ncol = 2 + length(radii)))
colnames(mixed_ACINP_df) <- c("spe", "reference", radii_colnames)

# Target is always A and B together
mixed_AE_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types), ncol = 2 + length(radii)))
colnames(mixed_AE_df) <- c("spe", "reference", radii_colnames)

## ACIN and CKR are twice as large
# (ref A and tar A or B) OR (ref B and tar B or A)
mixed_ACIN_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types)^2, ncol = 3 + length(radii)))
colnames(mixed_ACIN_df) <- c("spe", "reference", "target", radii_colnames)

# (ref A and tar A or B) OR (ref B and tar B or A)
mixed_CKR_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types)^2, ncol = 3 + length(radii)))
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






### Mixed spes 2D setup --------------------------------------------

n_slices <- 5

# Define APD and AMD data frames as well as constants
cell_types <- c("A", "B")
APD_pairs <- c("A/A", "A/B", "B/B")
AMD_pairs <- c("A/A", "A/B", "B/A", "B/B")

mixed_slices_APD_df_colnames <- c("spe", "slice", "pair", "APD")
mixed_slices_APD_df <- data.frame(matrix(nrow = n_mixed_spes * length(APD_pairs) * n_slices, ncol = length(mixed_slices_APD_df_colnames)))
colnames(mixed_slices_APD_df) <- mixed_slices_APD_df_colnames

mixed_slices_AMD_df_colnames <- c("spe", "slice", "reference", "target", "AMD")
mixed_slices_AMD_df <- data.frame(matrix(nrow = n_mixed_spes * length(AMD_pairs) * n_slices, ncol = length(mixed_slices_AMD_df_colnames)))
colnames(mixed_slices_AMD_df) <- mixed_slices_AMD_df_colnames

# Define MS, NMS, ACIN, ACINP, CKR, AE data frames as well as constants
radii <- seq(10, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

mixed_slices_MS_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
mixed_slices_MS_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types) * n_slices, ncol = length(mixed_slices_MS_df_colnames)))
colnames(mixed_slices_MS_df) <- mixed_slices_MS_df_colnames

mixed_slices_NMS_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
mixed_slices_NMS_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types) * n_slices, ncol = length(mixed_slices_MS_df_colnames)))
colnames(mixed_slices_NMS_df) <- mixed_slices_NMS_df_colnames

# Target is always A and B together
# Only choose prop(A) as prop(B) = 1 - prop(A) always
mixed_slices_ACINP_df_colnames <- c("spe", "slice", "reference", radii_colnames)
mixed_slices_ACINP_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types) * n_slices, ncol = length(mixed_slices_ACINP_df_colnames)))
colnames(mixed_slices_ACINP_df) <- mixed_slices_ACINP_df_colnames

# Target is always A and B together
mixed_slices_AE_df_colnames <- c("spe", "slice", "reference", radii_colnames)
mixed_slices_AE_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types) * n_slices, ncol = length(mixed_slices_AE_df_colnames)))
colnames(mixed_slices_AE_df) <- mixed_slices_AE_df_colnames

## ACIN and CKR are twice as large
# (ref A and tar A or B) OR (ref B and tar B or A)
mixed_slices_ACIN_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
mixed_slices_ACIN_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types)^2 * n_slices, ncol = length(mixed_slices_ACIN_df_colnames)))
colnames(mixed_slices_ACIN_df) <- mixed_slices_ACIN_df_colnames

# (ref A and tar A or B) OR (ref B and tar B or A)
mixed_slices_CKR_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
mixed_slices_CKR_df <- data.frame(matrix(nrow = n_mixed_spes * length(cell_types)^2 * n_slices, ncol = length(mixed_slices_CKR_df_colnames)))
colnames(mixed_slices_CKR_df) <- mixed_slices_CKR_df_colnames

# Define SAC and prevalence data frames as well as constants
n_splits <- 10
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))

mixed_slices_prop_SAC_df_colnames <- c("spe", "slice", "reference", "target", "proportion")
mixed_slices_prop_SAC_df <- data.frame(matrix(nrow = n_mixed_spes * nrow(prop_cell_types) * n_slices, ncol = length(mixed_slices_prop_SAC_df_colnames)))
colnames(mixed_slices_prop_SAC_df) <- mixed_slices_prop_SAC_df_colnames

mixed_slices_prop_prevalence_df_colnames <- c("spe", "slice", "reference", "target", thresholds_colnames)
mixed_slices_prop_prevalence_df <- data.frame(matrix(nrow = n_mixed_spes * nrow(prop_cell_types) * n_slices, ncol = length(mixed_slices_prop_prevalence_df_colnames)))
colnames(mixed_slices_prop_prevalence_df) <- mixed_slices_prop_prevalence_df_colnames


entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))

mixed_slices_entropy_SAC_df_colnames <- c("spe", "slice", "cell_types", "entropy")
mixed_slices_entropy_SAC_df <- data.frame(matrix(nrow = n_mixed_spes * nrow(entropy_cell_types) * n_slices, ncol = length(mixed_slices_entropy_SAC_df_colnames)))
colnames(mixed_slices_entropy_SAC_df) <- mixed_slices_entropy_SAC_df_colnames

mixed_slices_entropy_prevalence_df_colnames <- c("spe", "slice", "cell_types", thresholds_colnames)
mixed_slices_entropy_prevalence_df <- data.frame(matrix(nrow = n_mixed_spes * nrow(entropy_cell_types) * n_slices, ncol = length(mixed_slices_entropy_prevalence_df_colnames)))
colnames(mixed_slices_entropy_prevalence_df) <- mixed_slices_entropy_prevalence_df_colnames




### Mixed spe 3D and 2D analysis ------------------------------------------

for (i in seq_len(n_mixed_spes)) {
  print(i)
  ### 3D analysis -----------------------------
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
                                                          "rook")
    
    proportion_prevalence_df <- calculate_prevalence_gradient3D(proportion_grid_metrics,
                                                                "proportion",
                                                                show_AUC = F,
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
                                                       "rook")
    
    entropy_prevalence_df <- calculate_prevalence_gradient3D(entropy_grid_metrics,
                                                             "entropy",
                                                             show_AUC = F,
                                                             plot_image = F)
    
    index <- nrow(entropy_cell_types) * (i - 1) + j
    mixed_entropy_SAC_df[index, c("spe", "cell_types")] <- c(mixed_spe_name, entropy_cell_types$cell_types[j])
    mixed_entropy_SAC_df[index, "entropy"] <- entropy_SAC
    
    mixed_entropy_prevalence_df[index, c("spe", "cell_types")] <- c(mixed_spe_name, entropy_cell_types$cell_types[j])
    mixed_entropy_prevalence_df[index, thresholds_colnames] <- entropy_prevalence_df$prevalence
  }  

  
  ### 2D slicing analysis -------------------------
  mixed_spe_slices <- get_spe_slices_list(mixed_spe)

  for (slice_index in seq(n_slices)) {
    mixed_spe_slice <- mixed_spe_slices[[slice_index]]
    
    # pairwise_distance_data <- calculate_pairwise_distances_between_cell_types2D(mixed_spe_slice,
    #                                                                             cell_types,
    #                                                                             show_summary = F,
    #                                                                             plot_image = F)
    # pairwise_distance_data_summary <- summarise_distances_between_cell_types2D(pairwise_distance_data)
    # 
    # ## Fill in 3 rows at a time for APD df (as we have A/A, A/B, B/B)
    # index <- n_slices * 3 * (i - 1) + 3 * (slice_index - 1) + 1
    # mixed_slices_APD_df[index:(index + 2), "spe"] <- mixed_spe_name
    # mixed_slices_APD_df[index:(index + 2), "slice"] <- slice_index
    # mixed_slices_APD_df[index:(index + 2), "pair"] <- pairwise_distance_data_summary$pair
    # mixed_slices_APD_df[index:(index + 2), "APD"] <- pairwise_distance_data_summary$mean
    
    
    minimum_distance_data <- calculate_minimum_distances_between_cell_types2D(mixed_spe_slice,
                                                                              cell_types,
                                                                              show_summary = F,
                                                                              plot_image = F)
    
    minimum_distance_data_summary <- summarise_distances_between_cell_types2D(minimum_distance_data)
    ## Fill in 4 rows at a time for AMD df (as we have A/A, A/B, B/A, B/B)
    index <- n_slices * 4 * (i - 1) + 4 * (slice_index - 1) + 1
    mixed_slices_AMD_df[index:(index + 3), "spe"] <- mixed_spe_name
    mixed_slices_AMD_df[index:(index + 3), "slice"] <- slice_index
    mixed_slices_AMD_df[index:(index + 3), "reference"] <- minimum_distance_data_summary$reference
    mixed_slices_AMD_df[index:(index + 3), "target"] <- minimum_distance_data_summary$target
    mixed_slices_AMD_df[index:(index + 3), "AMD"] <- minimum_distance_data_summary$mean
    
    
    
    
    
    index1 <- n_slices * 2 * (i - 1) + 2 * (slice_index - 1) + 1
    index2 <- n_slices * 4 * (i - 1) + 4 * (slice_index - 1) + 1
    for (reference_cell_type in cell_types) {
      gradient_data <- calculate_all_gradient_cc_metrics2D(mixed_spe_slice,
                                                           reference_cell_type,
                                                           cell_types,
                                                           radii,
                                                           plot_image = F)
      
      target_cell_type <- setdiff(cell_types, reference_cell_type)
      
      mixed_slices_MS_df[index1, c("spe", "slice", "reference", "target")] <- c(mixed_spe_name, slice_index, reference_cell_type, target_cell_type)
      mixed_slices_NMS_df[index1, c("spe", "slice", "reference", "target")] <- c(mixed_spe_name, slice_index, reference_cell_type, target_cell_type)
      mixed_slices_ACINP_df[index1, c("spe", "slice","reference")] <- c(mixed_spe_name, slice_index, reference_cell_type)
      mixed_slices_AE_df[index1, c("spe", "slice","reference")] <- c(mixed_spe_name, slice_index, reference_cell_type)
      
      if (!is.null(gradient_data)) {
        mixed_slices_MS_df[index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
        mixed_slices_NMS_df[index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
        mixed_slices_ACINP_df[index1, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][["A"]]
        mixed_slices_AE_df[index1, radii_colnames] <- gradient_data[["entropy"]]$entropy        
      }
      else {
        mixed_slices_MS_df[index1, radii_colnames] <- NA
        mixed_slices_NMS_df[index1, radii_colnames] <- NA
        mixed_slices_ACINP_df[index1, radii_colnames] <- NA
        mixed_slices_AE_df[index1, radii_colnames] <- NA
      }

      
      index1 <- index1 + 1
      
      for (target_cell_type in cell_types) {
        ## Calculate ACIN and CKR as target cell type can also be the reference cell type
        
        # ACIN & CKR
        mixed_slices_ACIN_df[index2, c("spe", "slice", "reference", "target")] <- c(mixed_spe_name, slice_index, reference_cell_type, target_cell_type)
        mixed_slices_CKR_df[index2, c("spe", "slice", "reference", "target")] <- c(mixed_spe_name, slice_index, reference_cell_type, target_cell_type)
        
        if (!is.null(gradient_data)) {
          mixed_slices_ACIN_df[index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
          mixed_slices_CKR_df[index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]]$cross_K_ratio
        }
        else {
          mixed_slices_ACIN_df[index2, radii_colnames] <- NA
          mixed_slices_CKR_df[index2, radii_colnames] <- NA
        }
        
        index2 <- index2 + 1
      }
    }
    
    
    # Get proportion grid metrics
    for (j in seq_len(nrow(prop_cell_types))) {
      proportion_grid_metrics <- calculate_cell_proportion_grid_metrics2D(mixed_spe_slice, 
                                                                          n_splits,
                                                                          strsplit(prop_cell_types$ref[j], ",")[[1]], 
                                                                          strsplit(prop_cell_types$tar[j], ",")[[1]],
                                                                          plot_image = F)
      
      if (!is.null(proportion_grid_metrics)) {
        proportion_SAC <- calculate_spatial_autocorrelation2D(proportion_grid_metrics, 
                                                              "proportion",
                                                              "rook")
        
        proportion_prevalence_df <- calculate_prevalence_gradient2D(proportion_grid_metrics,
                                                                    "proportion",
                                                                    show_AUC = F,
                                                                    plot_image = F)
      }
      else {
        proportion_SAC <- NA
        proportion_prevalence_df <- data.frame(threshold = seq(0.01, 1, 0.01), prevalence = NA)
      }
      

      
      index <- n_slices * nrow(prop_cell_types) * (i - 1) + nrow(prop_cell_types) * (slice_index - 1) + j
      
      mixed_slices_prop_SAC_df[index, c("spe", "slice", "reference", "target")] <- c(mixed_spe_name, slice_index, prop_cell_types$ref[j], prop_cell_types$tar[j])
      mixed_slices_prop_SAC_df[index, "proportion"] <- proportion_SAC
      
      mixed_slices_prop_prevalence_df[index, c("spe", "slice", "reference", "target")] <- c(mixed_spe_name, slice_index, prop_cell_types$ref[j], prop_cell_types$tar[j])
      mixed_slices_prop_prevalence_df[index, thresholds_colnames] <- proportion_prevalence_df$prevalence
    }
    
    # Get entropy grid metrics
    for (j in seq_len(nrow(entropy_cell_types))) {
      entropy_grid_metrics <- calculate_entropy_grid_metrics2D(mixed_spe_slice, 
                                                               n_splits,
                                                               strsplit(entropy_cell_types$cell_types[j], ",")[[1]], 
                                                               plot_image = F)
      
      if (!is.null(entropy_grid_metrics)) {
        entropy_SAC <- calculate_spatial_autocorrelation2D(entropy_grid_metrics, 
                                                           "entropy",
                                                           "rook")
        
        entropy_prevalence_df <- calculate_prevalence_gradient2D(entropy_grid_metrics,
                                                                 "entropy",
                                                                 show_AUC = F,
                                                                 plot_image = F)
      }
      else {
        proportion_SAC <- NA
        proportion_prevalence_df <- data.frame(threshold = seq(0.01, 1, 0.01), prevalence = NA)
      }
      
      index <- n_slices * nrow(entropy_cell_types) * (i - 1) + nrow(entropy_cell_types) * (slice_index - 1) + j
      
      mixed_slices_entropy_SAC_df[index, c("spe", "slice", "cell_types")] <- c(mixed_spe_name, slice_index, entropy_cell_types$cell_types[j])
      mixed_slices_entropy_SAC_df[index, "entropy"] <- entropy_SAC
      
      mixed_slices_entropy_prevalence_df[index, c("spe", "slice", "cell_types")] <- c(mixed_spe_name, slice_index, entropy_cell_types$cell_types[j])
      mixed_slices_entropy_prevalence_df[index, thresholds_colnames] <- entropy_prevalence_df$prevalence
    }  
  }
}


setwd("~/Objects/unsupervised/mixed_spes/analysis_3D")
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


setwd("~/Objects/unsupervised/mixed_spes/analysis_2D")
write.table(mixed_slices_APD_df, file = "mixed_slices_APD_df.csv")
write.table(mixed_slices_AMD_df, file = "mixed_slices_AMD_df.csv")

write.table(mixed_slices_MS_df, file = "mixed_slices_MS_df.csv")
write.table(mixed_slices_NMS_df, file = "mixed_slices_NMS_df.csv")
write.table(mixed_slices_ACINP_df, file = "mixed_slices_ACINP_df.csv")
write.table(mixed_slices_AE_df, file = "mixed_slices_AE_df.csv")
write.table(mixed_slices_ACIN_df, file = "mixed_slices_ACIN_df.csv")
write.table(mixed_slices_CKR_df, file = "mixed_slices_CKR_df.csv")

write.table(mixed_slices_prop_SAC_df, file = "mixed_slices_prop_SAC_df.csv")
write.table(mixed_slices_prop_prevalence_df, file = "mixed_slices_prop_prevalence_df.csv")
write.table(mixed_slices_entropy_SAC_df, file = "mixed_slices_entropy_SAC_df.csv")
write.table(mixed_slices_entropy_prevalence_df, file = "mixed_slices_entropy_prevalence_df.csv")



### ringed spes 3D setup --------------------------------------------

setwd("~/Objects/unsupervised/spes_metadata")
ringed_spes_metadata <- readRDS("ringed_spes_metadata_unsupervised.rds")

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
radii <- seq(10, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

ringed_MS_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types), ncol = 3 + length(radii)))
colnames(ringed_MS_df) <- c("spe", "reference", "target", radii_colnames)

ringed_NMS_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types), ncol = 3 + length(radii)))
colnames(ringed_NMS_df) <- c("spe", "reference", "target", radii_colnames)

# Target is always A and B together
# Only choose prop(A) as prop(B) = 1 - prop(A) always
ringed_ACINP_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types), ncol = 2 + length(radii)))
colnames(ringed_ACINP_df) <- c("spe", "reference", radii_colnames)

# Target is always A and B together
ringed_AE_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types), ncol = 2 + length(radii)))
colnames(ringed_AE_df) <- c("spe", "reference", radii_colnames)

## ACIN and CKR are twice as large
# (ref A and tar A or B) OR (ref B and tar B or A)
ringed_ACIN_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types)^2, ncol = 3 + length(radii)))
colnames(ringed_ACIN_df) <- c("spe", "reference", "target", radii_colnames)

# (ref A and tar A or B) OR (ref B and tar B or A)
ringed_CKR_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types)^2, ncol = 3 + length(radii)))
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






### ringed spes 2D setup --------------------------------------------

n_slices <- 5

# Define APD and AMD data frames as well as constants
cell_types <- c("A", "B")
APD_pairs <- c("A/A", "A/B", "B/B")
AMD_pairs <- c("A/A", "A/B", "B/A", "B/B")

ringed_slices_APD_df_colnames <- c("spe", "slice", "pair", "APD")
ringed_slices_APD_df <- data.frame(matrix(nrow = n_ringed_spes * length(APD_pairs) * n_slices, ncol = length(ringed_slices_APD_df_colnames)))
colnames(ringed_slices_APD_df) <- ringed_slices_APD_df_colnames

ringed_slices_AMD_df_colnames <- c("spe", "slice", "reference", "target", "AMD")
ringed_slices_AMD_df <- data.frame(matrix(nrow = n_ringed_spes * length(AMD_pairs) * n_slices, ncol = length(ringed_slices_AMD_df_colnames)))
colnames(ringed_slices_AMD_df) <- ringed_slices_AMD_df_colnames

# Define MS, NMS, ACIN, ACINP, CKR, AE data frames as well as constants
radii <- seq(10, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

ringed_slices_MS_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
ringed_slices_MS_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types) * n_slices, ncol = length(ringed_slices_MS_df_colnames)))
colnames(ringed_slices_MS_df) <- ringed_slices_MS_df_colnames

ringed_slices_NMS_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
ringed_slices_NMS_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types) * n_slices, ncol = length(ringed_slices_MS_df_colnames)))
colnames(ringed_slices_NMS_df) <- ringed_slices_NMS_df_colnames

# Target is always A and B together
# Only choose prop(A) as prop(B) = 1 - prop(A) always
ringed_slices_ACINP_df_colnames <- c("spe", "slice", "reference", radii_colnames)
ringed_slices_ACINP_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types) * n_slices, ncol = length(ringed_slices_ACINP_df_colnames)))
colnames(ringed_slices_ACINP_df) <- ringed_slices_ACINP_df_colnames

# Target is always A and B together
ringed_slices_AE_df_colnames <- c("spe", "slice", "reference", radii_colnames)
ringed_slices_AE_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types) * n_slices, ncol = length(ringed_slices_AE_df_colnames)))
colnames(ringed_slices_AE_df) <- ringed_slices_AE_df_colnames

## ACIN and CKR are twice as large
# (ref A and tar A or B) OR (ref B and tar B or A)
ringed_slices_ACIN_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
ringed_slices_ACIN_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types)^2 * n_slices, ncol = length(ringed_slices_ACIN_df_colnames)))
colnames(ringed_slices_ACIN_df) <- ringed_slices_ACIN_df_colnames

# (ref A and tar A or B) OR (ref B and tar B or A)
ringed_slices_CKR_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
ringed_slices_CKR_df <- data.frame(matrix(nrow = n_ringed_spes * length(cell_types)^2 * n_slices, ncol = length(ringed_slices_CKR_df_colnames)))
colnames(ringed_slices_CKR_df) <- ringed_slices_CKR_df_colnames

# Define SAC and prevalence data frames as well as constants
n_splits <- 10
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))

ringed_slices_prop_SAC_df_colnames <- c("spe", "slice", "reference", "target", "proportion")
ringed_slices_prop_SAC_df <- data.frame(matrix(nrow = n_ringed_spes * nrow(prop_cell_types) * n_slices, ncol = length(ringed_slices_prop_SAC_df_colnames)))
colnames(ringed_slices_prop_SAC_df) <- ringed_slices_prop_SAC_df_colnames

ringed_slices_prop_prevalence_df_colnames <- c("spe", "slice", "reference", "target", thresholds_colnames)
ringed_slices_prop_prevalence_df <- data.frame(matrix(nrow = n_ringed_spes * nrow(prop_cell_types) * n_slices, ncol = length(ringed_slices_prop_prevalence_df_colnames)))
colnames(ringed_slices_prop_prevalence_df) <- ringed_slices_prop_prevalence_df_colnames


entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))

ringed_slices_entropy_SAC_df_colnames <- c("spe", "slice", "cell_types", "entropy")
ringed_slices_entropy_SAC_df <- data.frame(matrix(nrow = n_ringed_spes * nrow(entropy_cell_types) * n_slices, ncol = length(ringed_slices_entropy_SAC_df_colnames)))
colnames(ringed_slices_entropy_SAC_df) <- ringed_slices_entropy_SAC_df_colnames

ringed_slices_entropy_prevalence_df_colnames <- c("spe", "slice", "cell_types", thresholds_colnames)
ringed_slices_entropy_prevalence_df <- data.frame(matrix(nrow = n_ringed_spes * nrow(entropy_cell_types) * n_slices, ncol = length(ringed_slices_entropy_prevalence_df_colnames)))
colnames(ringed_slices_entropy_prevalence_df) <- ringed_slices_entropy_prevalence_df_colnames




### ringed spe 3D and 2D analysis ------------------------------------------

for (i in seq_len(n_ringed_spes)) {
  print(i)
  ### 3D analysis -----------------------------
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
                                                          "rook")
    
    proportion_prevalence_df <- calculate_prevalence_gradient3D(proportion_grid_metrics,
                                                                "proportion",
                                                                show_AUC = F,
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
                                                       "rook")
    
    entropy_prevalence_df <- calculate_prevalence_gradient3D(entropy_grid_metrics,
                                                             "entropy",
                                                             show_AUC = F,
                                                             plot_image = F)
    
    index <- nrow(entropy_cell_types) * (i - 1) + j
    ringed_entropy_SAC_df[index, c("spe", "cell_types")] <- c(ringed_spe_name, entropy_cell_types$cell_types[j])
    ringed_entropy_SAC_df[index, "entropy"] <- entropy_SAC
    
    ringed_entropy_prevalence_df[index, c("spe", "cell_types")] <- c(ringed_spe_name, entropy_cell_types$cell_types[j])
    ringed_entropy_prevalence_df[index, thresholds_colnames] <- entropy_prevalence_df$prevalence
  }  
  
  
  ### 2D slicing analysis -------------------------
  ringed_spe_slices <- get_spe_slices_list(ringed_spe)
  
  for (slice_index in seq(n_slices)) {
    ringed_spe_slice <- ringed_spe_slices[[slice_index]]
    
    # pairwise_distance_data <- calculate_pairwise_distances_between_cell_types2D(ringed_spe_slice,
    #                                                                             cell_types,
    #                                                                             show_summary = F,
    #                                                                             plot_image = F)
    # pairwise_distance_data_summary <- summarise_distances_between_cell_types2D(pairwise_distance_data)
    # 
    # ## Fill in 3 rows at a time for APD df (as we have A/A, A/B, B/B)
    # index <- n_slices * 3 * (i - 1) + 3 * (slice_index - 1) + 1
    # ringed_slices_APD_df[index:(index + 2), "spe"] <- ringed_spe_name
    # ringed_slices_APD_df[index:(index + 2), "slice"] <- slice_index
    # ringed_slices_APD_df[index:(index + 2), "pair"] <- pairwise_distance_data_summary$pair
    # ringed_slices_APD_df[index:(index + 2), "APD"] <- pairwise_distance_data_summary$mean
    
    
    minimum_distance_data <- calculate_minimum_distances_between_cell_types2D(ringed_spe_slice,
                                                                              cell_types,
                                                                              show_summary = F,
                                                                              plot_image = F)
    
    minimum_distance_data_summary <- summarise_distances_between_cell_types2D(minimum_distance_data)
    ## Fill in 4 rows at a time for AMD df (as we have A/A, A/B, B/A, B/B)
    index <- n_slices * 4 * (i - 1) + 4 * (slice_index - 1) + 1
    ringed_slices_AMD_df[index:(index + 3), "spe"] <- ringed_spe_name
    ringed_slices_AMD_df[index:(index + 3), "slice"] <- slice_index
    ringed_slices_AMD_df[index:(index + 3), "reference"] <- minimum_distance_data_summary$reference
    ringed_slices_AMD_df[index:(index + 3), "target"] <- minimum_distance_data_summary$target
    ringed_slices_AMD_df[index:(index + 3), "AMD"] <- minimum_distance_data_summary$mean
    
    
    
    
    
    index1 <- n_slices * 2 * (i - 1) + 2 * (slice_index - 1) + 1
    index2 <- n_slices * 4 * (i - 1) + 4 * (slice_index - 1) + 1
    for (reference_cell_type in cell_types) {
      gradient_data <- calculate_all_gradient_cc_metrics2D(ringed_spe_slice,
                                                           reference_cell_type,
                                                           cell_types,
                                                           radii,
                                                           plot_image = F)
      
      target_cell_type <- setdiff(cell_types, reference_cell_type)
      
      ringed_slices_MS_df[index1, c("spe", "slice", "reference", "target")] <- c(ringed_spe_name, slice_index, reference_cell_type, target_cell_type)
      ringed_slices_NMS_df[index1, c("spe", "slice", "reference", "target")] <- c(ringed_spe_name, slice_index, reference_cell_type, target_cell_type)
      ringed_slices_ACINP_df[index1, c("spe", "slice","reference")] <- c(ringed_spe_name, slice_index, reference_cell_type)
      ringed_slices_AE_df[index1, c("spe", "slice","reference")] <- c(ringed_spe_name, slice_index, reference_cell_type)
      
      if (!is.null(gradient_data)) {
        ringed_slices_MS_df[index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
        ringed_slices_NMS_df[index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
        ringed_slices_ACINP_df[index1, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][["A"]]
        ringed_slices_AE_df[index1, radii_colnames] <- gradient_data[["entropy"]]$entropy        
      }
      else {
        ringed_slices_MS_df[index1, radii_colnames] <- NA
        ringed_slices_NMS_df[index1, radii_colnames] <- NA
        ringed_slices_ACINP_df[index1, radii_colnames] <- NA
        ringed_slices_AE_df[index1, radii_colnames] <- NA
      }
      
      
      index1 <- index1 + 1
      
      for (target_cell_type in cell_types) {
        ## Calculate ACIN and CKR as target cell type can also be the reference cell type
        
        # ACIN & CKR
        ringed_slices_ACIN_df[index2, c("spe", "slice", "reference", "target")] <- c(ringed_spe_name, slice_index, reference_cell_type, target_cell_type)
        ringed_slices_CKR_df[index2, c("spe", "slice", "reference", "target")] <- c(ringed_spe_name, slice_index, reference_cell_type, target_cell_type)
        
        if (!is.null(gradient_data)) {
          ringed_slices_ACIN_df[index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
          ringed_slices_CKR_df[index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]]$cross_K_ratio
        }
        else {
          ringed_slices_ACIN_df[index2, radii_colnames] <- NA
          ringed_slices_CKR_df[index2, radii_colnames] <- NA
        }
        
        index2 <- index2 + 1
      }
    }
    
    
    # Get proportion grid metrics
    for (j in seq_len(nrow(prop_cell_types))) {
      proportion_grid_metrics <- calculate_cell_proportion_grid_metrics2D(ringed_spe_slice, 
                                                                          n_splits,
                                                                          strsplit(prop_cell_types$ref[j], ",")[[1]], 
                                                                          strsplit(prop_cell_types$tar[j], ",")[[1]],
                                                                          plot_image = F)
      
      if (!is.null(proportion_grid_metrics)) {
        proportion_SAC <- calculate_spatial_autocorrelation2D(proportion_grid_metrics, 
                                                              "proportion",
                                                              "rook")
        
        proportion_prevalence_df <- calculate_prevalence_gradient2D(proportion_grid_metrics,
                                                                    "proportion",
                                                                    show_AUC = F,
                                                                    plot_image = F)
      }
      else {
        proportion_SAC <- NA
        proportion_prevalence_df <- data.frame(threshold = seq(0.01, 1, 0.01), prevalence = NA)
      }
      
      
      
      index <- n_slices * nrow(prop_cell_types) * (i - 1) + nrow(prop_cell_types) * (slice_index - 1) + j
      
      ringed_slices_prop_SAC_df[index, c("spe", "slice", "reference", "target")] <- c(ringed_spe_name, slice_index, prop_cell_types$ref[j], prop_cell_types$tar[j])
      ringed_slices_prop_SAC_df[index, "proportion"] <- proportion_SAC
      
      ringed_slices_prop_prevalence_df[index, c("spe", "slice", "reference", "target")] <- c(ringed_spe_name, slice_index, prop_cell_types$ref[j], prop_cell_types$tar[j])
      ringed_slices_prop_prevalence_df[index, thresholds_colnames] <- proportion_prevalence_df$prevalence
    }
    
    # Get entropy grid metrics
    for (j in seq_len(nrow(entropy_cell_types))) {
      entropy_grid_metrics <- calculate_entropy_grid_metrics2D(ringed_spe_slice, 
                                                               n_splits,
                                                               strsplit(entropy_cell_types$cell_types[j], ",")[[1]], 
                                                               plot_image = F)
      
      if (!is.null(entropy_grid_metrics)) {
        entropy_SAC <- calculate_spatial_autocorrelation2D(entropy_grid_metrics, 
                                                           "entropy",
                                                           "rook")
        
        entropy_prevalence_df <- calculate_prevalence_gradient2D(entropy_grid_metrics,
                                                                 "entropy",
                                                                 show_AUC = F,
                                                                 plot_image = F)
      }
      else {
        proportion_SAC <- NA
        proportion_prevalence_df <- data.frame(threshold = seq(0.01, 1, 0.01), prevalence = NA)
      }
      
      index <- n_slices * nrow(entropy_cell_types) * (i - 1) + nrow(entropy_cell_types) * (slice_index - 1) + j
      
      ringed_slices_entropy_SAC_df[index, c("spe", "slice", "cell_types")] <- c(ringed_spe_name, slice_index, entropy_cell_types$cell_types[j])
      ringed_slices_entropy_SAC_df[index, "entropy"] <- entropy_SAC
      
      ringed_slices_entropy_prevalence_df[index, c("spe", "slice", "cell_types")] <- c(ringed_spe_name, slice_index, entropy_cell_types$cell_types[j])
      ringed_slices_entropy_prevalence_df[index, thresholds_colnames] <- entropy_prevalence_df$prevalence
    }  
  }
}


setwd("~/Objects/unsupervised/ringed_spes/analysis_3D")
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


setwd("~/Objects/unsupervised/ringed_spes/analysis_2D")
write.table(ringed_slices_APD_df, file = "ringed_slices_APD_df.csv")
write.table(ringed_slices_AMD_df, file = "ringed_slices_AMD_df.csv")

write.table(ringed_slices_MS_df, file = "ringed_slices_MS_df.csv")
write.table(ringed_slices_NMS_df, file = "ringed_slices_NMS_df.csv")
write.table(ringed_slices_ACINP_df, file = "ringed_slices_ACINP_df.csv")
write.table(ringed_slices_AE_df, file = "ringed_slices_AE_df.csv")
write.table(ringed_slices_ACIN_df, file = "ringed_slices_ACIN_df.csv")
write.table(ringed_slices_CKR_df, file = "ringed_slices_CKR_df.csv")

write.table(ringed_slices_prop_SAC_df, file = "ringed_slices_prop_SAC_df.csv")
write.table(ringed_slices_prop_prevalence_df, file = "ringed_slices_prop_prevalence_df.csv")
write.table(ringed_slices_entropy_SAC_df, file = "ringed_slices_entropy_SAC_df.csv")
write.table(ringed_slices_entropy_prevalence_df, file = "ringed_slices_entropy_prevalence_df.csv")



### separated spes 3D setup --------------------------------------------

setwd("~/Objects/unsupervised/spes_metadata")
separated_spes_metadata <- readRDS("separated_spes_metadata_unsupervised.rds")

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
radii <- seq(10, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

separated_MS_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types), ncol = 3 + length(radii)))
colnames(separated_MS_df) <- c("spe", "reference", "target", radii_colnames)

separated_NMS_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types), ncol = 3 + length(radii)))
colnames(separated_NMS_df) <- c("spe", "reference", "target", radii_colnames)

# Target is always A and B together
# Only choose prop(A) as prop(B) = 1 - prop(A) always
separated_ACINP_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types), ncol = 2 + length(radii)))
colnames(separated_ACINP_df) <- c("spe", "reference", radii_colnames)

# Target is always A and B together
separated_AE_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types), ncol = 2 + length(radii)))
colnames(separated_AE_df) <- c("spe", "reference", radii_colnames)

## ACIN and CKR are twice as large
# (ref A and tar A or B) OR (ref B and tar B or A)
separated_ACIN_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types)^2, ncol = 3 + length(radii)))
colnames(separated_ACIN_df) <- c("spe", "reference", "target", radii_colnames)

# (ref A and tar A or B) OR (ref B and tar B or A)
separated_CKR_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types)^2, ncol = 3 + length(radii)))
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






### separated spes 2D setup --------------------------------------------

n_slices <- 5

# Define APD and AMD data frames as well as constants
cell_types <- c("A", "B")
APD_pairs <- c("A/A", "A/B", "B/B")
AMD_pairs <- c("A/A", "A/B", "B/A", "B/B")

separated_slices_APD_df_colnames <- c("spe", "slice", "pair", "APD")
separated_slices_APD_df <- data.frame(matrix(nrow = n_separated_spes * length(APD_pairs) * n_slices, ncol = length(separated_slices_APD_df_colnames)))
colnames(separated_slices_APD_df) <- separated_slices_APD_df_colnames

separated_slices_AMD_df_colnames <- c("spe", "slice", "reference", "target", "AMD")
separated_slices_AMD_df <- data.frame(matrix(nrow = n_separated_spes * length(AMD_pairs) * n_slices, ncol = length(separated_slices_AMD_df_colnames)))
colnames(separated_slices_AMD_df) <- separated_slices_AMD_df_colnames

# Define MS, NMS, ACIN, ACINP, CKR, AE data frames as well as constants
radii <- seq(10, 100, 10)
radii_colnames <- paste("r", radii, sep = "")

separated_slices_MS_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
separated_slices_MS_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types) * n_slices, ncol = length(separated_slices_MS_df_colnames)))
colnames(separated_slices_MS_df) <- separated_slices_MS_df_colnames

separated_slices_NMS_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
separated_slices_NMS_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types) * n_slices, ncol = length(separated_slices_MS_df_colnames)))
colnames(separated_slices_NMS_df) <- separated_slices_NMS_df_colnames

# Target is always A and B together
# Only choose prop(A) as prop(B) = 1 - prop(A) always
separated_slices_ACINP_df_colnames <- c("spe", "slice", "reference", radii_colnames)
separated_slices_ACINP_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types) * n_slices, ncol = length(separated_slices_ACINP_df_colnames)))
colnames(separated_slices_ACINP_df) <- separated_slices_ACINP_df_colnames

# Target is always A and B together
separated_slices_AE_df_colnames <- c("spe", "slice", "reference", radii_colnames)
separated_slices_AE_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types) * n_slices, ncol = length(separated_slices_AE_df_colnames)))
colnames(separated_slices_AE_df) <- separated_slices_AE_df_colnames

## ACIN and CKR are twice as large
# (ref A and tar A or B) OR (ref B and tar B or A)
separated_slices_ACIN_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
separated_slices_ACIN_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types)^2 * n_slices, ncol = length(separated_slices_ACIN_df_colnames)))
colnames(separated_slices_ACIN_df) <- separated_slices_ACIN_df_colnames

# (ref A and tar A or B) OR (ref B and tar B or A)
separated_slices_CKR_df_colnames <- c("spe", "slice", "reference", "target", radii_colnames)
separated_slices_CKR_df <- data.frame(matrix(nrow = n_separated_spes * length(cell_types)^2 * n_slices, ncol = length(separated_slices_CKR_df_colnames)))
colnames(separated_slices_CKR_df) <- separated_slices_CKR_df_colnames

# Define SAC and prevalence data frames as well as constants
n_splits <- 10
thresholds <- seq(0.01, 1, 0.01)
thresholds_colnames <- paste("t", thresholds, sep = "")

prop_cell_types <- data.frame(ref = c("A", "O"), tar = c("B", "A,B"))

separated_slices_prop_SAC_df_colnames <- c("spe", "slice", "reference", "target", "proportion")
separated_slices_prop_SAC_df <- data.frame(matrix(nrow = n_separated_spes * nrow(prop_cell_types) * n_slices, ncol = length(separated_slices_prop_SAC_df_colnames)))
colnames(separated_slices_prop_SAC_df) <- separated_slices_prop_SAC_df_colnames

separated_slices_prop_prevalence_df_colnames <- c("spe", "slice", "reference", "target", thresholds_colnames)
separated_slices_prop_prevalence_df <- data.frame(matrix(nrow = n_separated_spes * nrow(prop_cell_types) * n_slices, ncol = length(separated_slices_prop_prevalence_df_colnames)))
colnames(separated_slices_prop_prevalence_df) <- separated_slices_prop_prevalence_df_colnames


entropy_cell_types <- data.frame(cell_types = c("A,B", "A,B,O"))

separated_slices_entropy_SAC_df_colnames <- c("spe", "slice", "cell_types", "entropy")
separated_slices_entropy_SAC_df <- data.frame(matrix(nrow = n_separated_spes * nrow(entropy_cell_types) * n_slices, ncol = length(separated_slices_entropy_SAC_df_colnames)))
colnames(separated_slices_entropy_SAC_df) <- separated_slices_entropy_SAC_df_colnames

separated_slices_entropy_prevalence_df_colnames <- c("spe", "slice", "cell_types", thresholds_colnames)
separated_slices_entropy_prevalence_df <- data.frame(matrix(nrow = n_separated_spes * nrow(entropy_cell_types) * n_slices, ncol = length(separated_slices_entropy_prevalence_df_colnames)))
colnames(separated_slices_entropy_prevalence_df) <- separated_slices_entropy_prevalence_df_colnames




### separated spe 3D and 2D analysis ------------------------------------------

for (i in seq_len(n_separated_spes)) {
  print(i)
  ### 3D analysis -----------------------------
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
                                                          "rook")
    
    proportion_prevalence_df <- calculate_prevalence_gradient3D(proportion_grid_metrics,
                                                                "proportion",
                                                                show_AUC = F,
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
                                                       "rook")
    
    entropy_prevalence_df <- calculate_prevalence_gradient3D(entropy_grid_metrics,
                                                             "entropy",
                                                             show_AUC = F,
                                                             plot_image = F)
    
    index <- nrow(entropy_cell_types) * (i - 1) + j
    separated_entropy_SAC_df[index, c("spe", "cell_types")] <- c(separated_spe_name, entropy_cell_types$cell_types[j])
    separated_entropy_SAC_df[index, "entropy"] <- entropy_SAC
    
    separated_entropy_prevalence_df[index, c("spe", "cell_types")] <- c(separated_spe_name, entropy_cell_types$cell_types[j])
    separated_entropy_prevalence_df[index, thresholds_colnames] <- entropy_prevalence_df$prevalence
  }  
  
  
  ### 2D slicing analysis -------------------------
  separated_spe_slices <- get_spe_slices_list(separated_spe)
  
  for (slice_index in seq(n_slices)) {
    separated_spe_slice <- separated_spe_slices[[slice_index]]
    
    # pairwise_distance_data <- calculate_pairwise_distances_between_cell_types2D(separated_spe_slice,
    #                                                                             cell_types,
    #                                                                             show_summary = F,
    #                                                                             plot_image = F)
    # pairwise_distance_data_summary <- summarise_distances_between_cell_types2D(pairwise_distance_data)
    # 
    # ## Fill in 3 rows at a time for APD df (as we have A/A, A/B, B/B)
    # index <- n_slices * 3 * (i - 1) + 3 * (slice_index - 1) + 1
    # separated_slices_APD_df[index:(index + 2), "spe"] <- separated_spe_name
    # separated_slices_APD_df[index:(index + 2), "slice"] <- slice_index
    # separated_slices_APD_df[index:(index + 2), "pair"] <- pairwise_distance_data_summary$pair
    # separated_slices_APD_df[index:(index + 2), "APD"] <- pairwise_distance_data_summary$mean
    
    
    minimum_distance_data <- calculate_minimum_distances_between_cell_types2D(separated_spe_slice,
                                                                              cell_types,
                                                                              show_summary = F,
                                                                              plot_image = F)
    
    minimum_distance_data_summary <- summarise_distances_between_cell_types2D(minimum_distance_data)
    ## Fill in 4 rows at a time for AMD df (as we have A/A, A/B, B/A, B/B)
    index <- n_slices * 4 * (i - 1) + 4 * (slice_index - 1) + 1
    separated_slices_AMD_df[index:(index + 3), "spe"] <- separated_spe_name
    separated_slices_AMD_df[index:(index + 3), "slice"] <- slice_index
    separated_slices_AMD_df[index:(index + 3), "reference"] <- minimum_distance_data_summary$reference
    separated_slices_AMD_df[index:(index + 3), "target"] <- minimum_distance_data_summary$target
    separated_slices_AMD_df[index:(index + 3), "AMD"] <- minimum_distance_data_summary$mean
    
    
    
    
    
    index1 <- n_slices * 2 * (i - 1) + 2 * (slice_index - 1) + 1
    index2 <- n_slices * 4 * (i - 1) + 4 * (slice_index - 1) + 1
    for (reference_cell_type in cell_types) {
      gradient_data <- calculate_all_gradient_cc_metrics2D(separated_spe_slice,
                                                           reference_cell_type,
                                                           cell_types,
                                                           radii,
                                                           plot_image = F)
      
      target_cell_type <- setdiff(cell_types, reference_cell_type)
      
      separated_slices_MS_df[index1, c("spe", "slice", "reference", "target")] <- c(separated_spe_name, slice_index, reference_cell_type, target_cell_type)
      separated_slices_NMS_df[index1, c("spe", "slice", "reference", "target")] <- c(separated_spe_name, slice_index, reference_cell_type, target_cell_type)
      separated_slices_ACINP_df[index1, c("spe", "slice","reference")] <- c(separated_spe_name, slice_index, reference_cell_type)
      separated_slices_AE_df[index1, c("spe", "slice","reference")] <- c(separated_spe_name, slice_index, reference_cell_type)
      
      if (!is.null(gradient_data)) {
        separated_slices_MS_df[index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$mixing_score
        separated_slices_NMS_df[index1, radii_colnames] <- gradient_data[["mixing_score"]][[target_cell_type]]$normalised_mixing_score
        separated_slices_ACINP_df[index1, radii_colnames] <- gradient_data[["cells_in_neighbourhood_proportion"]][["A"]]
        separated_slices_AE_df[index1, radii_colnames] <- gradient_data[["entropy"]]$entropy        
      }
      else {
        separated_slices_MS_df[index1, radii_colnames] <- NA
        separated_slices_NMS_df[index1, radii_colnames] <- NA
        separated_slices_ACINP_df[index1, radii_colnames] <- NA
        separated_slices_AE_df[index1, radii_colnames] <- NA
      }
      
      
      index1 <- index1 + 1
      
      for (target_cell_type in cell_types) {
        ## Calculate ACIN and CKR as target cell type can also be the reference cell type
        
        # ACIN & CKR
        separated_slices_ACIN_df[index2, c("spe", "slice", "reference", "target")] <- c(separated_spe_name, slice_index, reference_cell_type, target_cell_type)
        separated_slices_CKR_df[index2, c("spe", "slice", "reference", "target")] <- c(separated_spe_name, slice_index, reference_cell_type, target_cell_type)
        
        if (!is.null(gradient_data)) {
          separated_slices_ACIN_df[index2, radii_colnames] <- gradient_data[["cells_in_neighbourhood"]][[target_cell_type]]
          separated_slices_CKR_df[index2, radii_colnames] <- gradient_data[["cross_K"]][[target_cell_type]]$cross_K_ratio
        }
        else {
          separated_slices_ACIN_df[index2, radii_colnames] <- NA
          separated_slices_CKR_df[index2, radii_colnames] <- NA
        }
        
        index2 <- index2 + 1
      }
    }
    
    
    # Get proportion grid metrics
    for (j in seq_len(nrow(prop_cell_types))) {
      proportion_grid_metrics <- calculate_cell_proportion_grid_metrics2D(separated_spe_slice, 
                                                                          n_splits,
                                                                          strsplit(prop_cell_types$ref[j], ",")[[1]], 
                                                                          strsplit(prop_cell_types$tar[j], ",")[[1]],
                                                                          plot_image = F)
      
      if (!is.null(proportion_grid_metrics)) {
        proportion_SAC <- calculate_spatial_autocorrelation2D(proportion_grid_metrics, 
                                                              "proportion",
                                                              "rook")
        
        proportion_prevalence_df <- calculate_prevalence_gradient2D(proportion_grid_metrics,
                                                                    "proportion",
                                                                    show_AUC = F,
                                                                    plot_image = F)
      }
      else {
        proportion_SAC <- NA
        proportion_prevalence_df <- data.frame(threshold = seq(0.01, 1, 0.01), prevalence = NA)
      }
      
      
      
      index <- n_slices * nrow(prop_cell_types) * (i - 1) + nrow(prop_cell_types) * (slice_index - 1) + j
      
      separated_slices_prop_SAC_df[index, c("spe", "slice", "reference", "target")] <- c(separated_spe_name, slice_index, prop_cell_types$ref[j], prop_cell_types$tar[j])
      separated_slices_prop_SAC_df[index, "proportion"] <- proportion_SAC
      
      separated_slices_prop_prevalence_df[index, c("spe", "slice", "reference", "target")] <- c(separated_spe_name, slice_index, prop_cell_types$ref[j], prop_cell_types$tar[j])
      separated_slices_prop_prevalence_df[index, thresholds_colnames] <- proportion_prevalence_df$prevalence
    }
    
    # Get entropy grid metrics
    for (j in seq_len(nrow(entropy_cell_types))) {
      entropy_grid_metrics <- calculate_entropy_grid_metrics2D(separated_spe_slice, 
                                                               n_splits,
                                                               strsplit(entropy_cell_types$cell_types[j], ",")[[1]], 
                                                               plot_image = F)
      
      if (!is.null(entropy_grid_metrics)) {
        entropy_SAC <- calculate_spatial_autocorrelation2D(entropy_grid_metrics, 
                                                           "entropy",
                                                           "rook")
        
        entropy_prevalence_df <- calculate_prevalence_gradient2D(entropy_grid_metrics,
                                                                 "entropy",
                                                                 show_AUC = F,
                                                                 plot_image = F)
      }
      else {
        proportion_SAC <- NA
        proportion_prevalence_df <- data.frame(threshold = seq(0.01, 1, 0.01), prevalence = NA)
      }
      
      index <- n_slices * nrow(entropy_cell_types) * (i - 1) + nrow(entropy_cell_types) * (slice_index - 1) + j
      
      separated_slices_entropy_SAC_df[index, c("spe", "slice", "cell_types")] <- c(separated_spe_name, slice_index, entropy_cell_types$cell_types[j])
      separated_slices_entropy_SAC_df[index, "entropy"] <- entropy_SAC
      
      separated_slices_entropy_prevalence_df[index, c("spe", "slice", "cell_types")] <- c(separated_spe_name, slice_index, entropy_cell_types$cell_types[j])
      separated_slices_entropy_prevalence_df[index, thresholds_colnames] <- entropy_prevalence_df$prevalence
    }  
  }
}


setwd("~/Objects/unsupervised/separated_spes/analysis_3D")
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
write.table(separated_entropy_prevalence_df, file = "separated_entropy_prevalence_df.csv")


setwd("~/Objects/unsupervised/separated_spes/analysis_2D")
write.table(separated_slices_APD_df, file = "separated_slices_APD_df.csv")
write.table(separated_slices_AMD_df, file = "separated_slices_AMD_df.csv")

write.table(separated_slices_MS_df, file = "separated_slices_MS_df.csv")
write.table(separated_slices_NMS_df, file = "separated_slices_NMS_df.csv")
write.table(separated_slices_ACINP_df, file = "separated_slices_ACINP_df.csv")
write.table(separated_slices_AE_df, file = "separated_slices_AE_df.csv")
write.table(separated_slices_ACIN_df, file = "separated_slices_ACIN_df.csv")
write.table(separated_slices_CKR_df, file = "separated_slices_CKR_df.csv")

write.table(separated_slices_prop_SAC_df, file = "separated_slices_prop_SAC_df.csv")
write.table(separated_slices_prop_prevalence_df, file = "separated_slices_prop_prevalence_df.csv")
write.table(separated_slices_entropy_SAC_df, file = "separated_slices_entropy_SAC_df.csv")
write.table(separated_slices_entropy_prevalence_df, file = "separated_slices_entropy_prevalence_df.csv")


