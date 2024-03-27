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



###----------------------
# Construct a dummy marker intensity matrix
## rows are markers, columns are cells
intensity_matrix <- matrix(c(14.557, 0.169, 1.655, 0.054,
                             17.588, 0.229, 1.188, 2.074, 
                             21.262, 4.206,  5.924, 0.021), nrow = 4, ncol = 3)
# define marker names as rownames
rownames(intensity_matrix) <- c("DAPI", "CD3", "CD4", "AMACR")
# define cell IDs as colnames
colnames(intensity_matrix) <- c("Cell_1", "Cell_2", "Cell_3") 

# Construct a dummy metadata (phenotypes, x/y coordinates)
# the order of the elements in these vectors correspond to the cell order 
# in `intensity matrix`
phenotypes <- c("OTHER",  "AMACR", "CD3,CD4")
coord_x <- c(82, 171, 184)
coord_y <- c(30, 22, 38)

general_format_image <- format_image_to_spe(format = "general", 
                                            intensity_matrix = intensity_matrix,
                                            phenotypes = phenotypes, 
                                            coord_x = coord_x,
                                            coord_y = coord_y)
colData(general_format_image)
spatialCoords(general_format_image)
assay(general_format_image)

