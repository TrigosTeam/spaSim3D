### 1.1. Mixed spes - cc distance metrics ----------------------------------------

# Get mixed_spes_metadata and mixed_spes_table
setwd("~/Objects/spes_metadata")
mixed_spes_metadata <- readRDS("mixed_spes_metadata.rds")

setwd("~/Objects/spes_table")
mixed_spes_table <- readRDS("mixed_spes_table.rds")

# Get number of mixed spes
n_mixed_spes <- nrow(mixed_spes_table)

# Define APD and AMD data frames as well as constants
cell_types <- c("A", "B")
APD_pairs <- c("A/A", "A/B", "B/B")
AMD_pairs <- c("A/A", "A/B", "B/A", "B/B")

mixed_APD_df <- data.frame(matrix(nrow = n_mixed_spes * length(APD_pairs), ncol = 3))
colnames(mixed_APD_df) <- c("spe", "pair", "APD")

mixed_AMD_df <- data.frame(matrix(nrow = n_mixed_spes * length(AMD_pairs), ncol = 3))
colnames(mixed_AMD_df) <- c("spe", "pair", "AMD")


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
  mixed_AMD_df[index:(index + 3), "pair"] <- minimum_distance_data_summary$pair
  mixed_AMD_df[index:(index + 3), "AMD"] <- minimum_distance_data_summary$mean
  
}



### 1.2. Mixed spes - cc gradient based metrics ----------------------------------



### 1.3. Mixed spes - heterogeneity metrics --------------------------------------



### Spacer --------------------------------

