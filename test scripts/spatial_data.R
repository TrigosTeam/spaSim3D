library(SPIAT)
setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/Data from Anna/CELL_SEGMENTATION/TUMOUR/BRCA1/B1")

data <- read.delim(file = 'B1_[40986,17944]_cell_seg_data.txt',
                   header = T,
                   sep = '\t')

## Put phenotypes and x,y coords from data into spe object
spatial_data <- format_image_to_spe(phenotypes = data$Phenotype,
                                    coord_x = data$Cell.X.Position,
                                    coord_y = data$Cell.Y.Position)

## Remove empty phenotypes
spatial_data <- spatial_data[, spatial_data$Phenotype != ""]

## Get unique phenotypes from spatial_data
phenos <- unique(spatial_data$Phenotype)
phenos

## Determine cell types for each phenotype
cell_types <- c("Tumour",
                "Other",
                "CD4 T",
                "CD8 T",
                "T-reg",
                "Double negative T",
                "PDL1+ Immune")


formatted_spatial_data <- define_celltypes(spe_object = spatial_data, 
                                           categories = phenos, 
                                           category_colname = "Phenotype",
                                           names = cell_types,
                                           new_colname = "Cell.Type",
                                           print_names = FALSE)



avg_d <- calculate_pairwise_distances_between_celltypes(spe_object = formatted_spatial_data, 
                                                        cell_types_of_interest = c("Tumour"), 
                                                        feature_colname = "Cell.Type")

min_d <- calculate_minimum_distances_between_celltypes(spe_object = formatted_spatial_data, 
                                                       feature_colname = "Cell.Type",
                                                       cell_types_of_interest = c("CD4 T", "Tumour"))

###---------------------------------

data("simulated_image")
simulated_image$Phenotype
assay(simulated_image[1:5, 1:5])

predicted_image <- predict_phenotypes(spe_object = simulated_image,
                                      thresholds = NULL,
                                      tumour_marker = "Tumour_marker",
                                      baseline_markers = c("Immune_marker1", 
                                                           "Immune_marker2", 
                                                           "Immune_marker3", 
                                                           "Immune_marker4"),
                                      reference_phenotypes = FALSE)
