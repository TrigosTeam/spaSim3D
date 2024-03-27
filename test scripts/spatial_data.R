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
                "Helper T",
                "Cytotoxic T",
                "Regulatory T",
                "Double negative T",
                "PDL1+ Immune")


formatted_spatial_data <- define_celltypes(spe_object = spatial_data, 
                                           categories = phenos, 
                                           category_colname = "Phenotype",
                                           names = cell_types,
                                           new_colname = "Cell.Type",
                                           print_names = FALSE)



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
