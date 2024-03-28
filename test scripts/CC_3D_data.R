setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/Publicly available 3D spatial data/Colorectal Cancer")

## Read data
data <- read.csv(file = "WD-76845-002.csv",
                 header = TRUE)

### Use SPIAT to format the data into SPE object
library(SPIAT)

## Make intensity matrix
# Find start and end column indicies
start <- "CD3"
end <- "Collagen"
start <- which(col_names %in% start)
end <- which(col_names %in% end)

col_names <- colnames(data)

intensity_matrix <- data[ , c(start:end)]

# Need to remove fluorescent dyes still remaining
intensity_matrix <- subset(intensity_matrix, select = -c(AF555,
                                                         AF647,
                                                         A555,
                                                         A647))

intensity_matrix <- t(intensity_matrix)
colnames(intensity_matrix) <- paste("Cell_", 1:ncol(intensity_matrix), sep="")

intensity_matrix <- intensity_matrix[, 1:10000] # Choose first 100 to see if it works

## Get x and y coords
coord_x <- data$X[1:10000]
coord_y <- data$Y[1:10000]


## Define phenotypes as a vector of NAs
phenotypes <- rep(NA, length(coord_x))


general_format_image <- format_image_to_spe(format = "general",
                                            intensity_matrix = intensity_matrix,
                                            phenotypes = phenotypes,
                                            coord_x = coord_x,
                                            coord_y = coord_y)

## Examine SPE object
colData(general_format_image)[1:5, 1:3]
spatialCoords(general_format_image)[1:5, ]
assay(general_format_image)[1:5, 1:5]


## Predict phenotypes for cells
markers <- rownames(assay(general_format_image))
markers

predicted_image <- predict_phenotypes(spe_object = general_format_image,
                                      thresholds = NULL,
                                      tumour_marker = "Ki67",
                                      baseline_markers = c("CD3",
                                                           "CD4",
                                                           "CD8a",
                                                           "CD163",
                                                           "CD45",
                                                           "CD45RO",
                                                           "CD20",
                                                           "CD68"),
                                      reference_phenotypes = FALSE)


predicted_image$Phenotype[1:100]


## Specify cell types using predicted phenotypes
formatted_image <- define_celltypes(spe_object = predicted_image,
                                    category_colname = "Phenotype",
                                    new_colname = "Cell.Type")
# Need to define cell_type to each phenotype
formatted_image$Cell.Type[1:30]


## Plotting
plot_cell_marker_levels(predicted_image, "CD3")
plot_marker_level_heatmap(predicted_image, num_splits = 100, "Keratin")
