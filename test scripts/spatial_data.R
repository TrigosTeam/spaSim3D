library(SPIAT)
setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/Data from Anna/CELL_SEGMENTATION/TUMOUR/BRCA1/B1")

data <- read.delim(file = 'B1_[40986,17944]_cell_seg_data.txt',
                   header = T,
                   sep = '\t')

spatial_data <- data[ , c('Phenotype', 'Cell.X.Position', 'Cell.Y.Position')]


data <- define_celltypes(spe_object = spatial_data, 
                         categories = NULL, 
                         category_colname = "Phenotype", 
                         names = NULL, 
                         new_colname = "Cell.Type", 
                         print_names = FALSE)


