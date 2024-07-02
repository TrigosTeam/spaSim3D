## Set working directory
setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/3D public spatial data/Lin et al - human colorectal cancer/CRC1 data")

## Read data
slice_nums <- c("002", "007", "014", "020", "025", "029",
                "034", "039", "044", "049", "050", "051",
                "052", "054", "059", "064", "069", "074",
                "078", "084", "086", "091", "097", "102", "106")
file_names <- paste("WD-76845-", slice_nums, ".csv", sep = "")

dataset <- vector(mode = "list", length = length(file_names))

for (i in seq(length(file_names))) {
  dataset[[i]] <- read.csv(file_names[i])
}

## Get relevant cell markers
all_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                      "CD45RO", "CD4", "CD8a", "CD68", 
                      "CD163", "FOXP3", "PD1", "PDL1", 
                      "CD31", "aSMA", "Desmin", "CD45")

## Check cell markers are in data
all_cell_markers %in% colnames(data)


### 1. Cell type annotations with only tumour and immune cells ----------------

# Tumour: Keratin+
# Immune: CD45+
# Other: ??

data <- dataset[[1]]

## Get chosen cell markers
chosen_cell_markers <- c("Keratin", "CD45")
chosen_cell_markers %in% colnames(data)

intensity_matrix <- t(data[, c(chosen_cell_markers)])
colnames(intensity_matrix) <- paste("Cell_", seq(ncol(intensity_matrix)), sep="")

## Get x and y coords
coord_x <- data[ , "X"]
coord_y <- data[ , "Y"]

## Use SPIAT to get cell phenotypes
library(SPIAT)

# Format into spe object
general_format_image <- format_image_to_spe(format = "general",
                                            intensity_matrix = intensity_matrix,
                                            phenotypes = "unknown",
                                            coord_x = coord_x,
                                            coord_y = coord_y)

# Predict phenotypes
predicted_image <- predict_phenotypes(spe_object = general_format_image,
                                      thresholds = NULL,
                                      tumour_marker = "Keratin",
                                      baseline_markers = c("CD45"),
                                      reference_phenotype = FALSE)

colData(predicted_image)
table(predicted_image$Phenotype)

# Get cell types from phenotypes
formatted_image <- define_celltypes(
  predicted_image, 
  categories = c("None", "Keratin", "CD45", "Keratin,CD45"), 
  category_colname = "Phenotype", 
  names = c("Other", "Tumour", "Immune", "Other"),
  new_colname = "Cell.Type"
)


## Define cell types
cell_types <- c("Tumour", "Immune", "Other")

plot_data <- data.frame(spatialCoords(formatted_image))
plot_data$Cell.Type <- factor(formatted_image$Cell.Type, levels = cell_types)


## Plot data
library(ggplot2)
library(dplyr)
my_colours <- c("red", "blue", "lightgray")

ggplot(sample_n(plot_data, 5000), aes(Cell.X.Position, Cell.Y.Position, colour = Cell.Type)) +
  geom_point() +
  scale_colour_manual(values = my_colours)




### 2. Cell type annotations with only tumour, immune and stroma cells ------

# Tumour: Keratin+
# Immune: CD45+
# Stroma: CD31+, aSMA+
# Other: ??

data <- dataset[[1]]

## Get chosen cell markers
chosen_cell_markers <- c("Keratin", "CD45", "CD31", "aSMA")
chosen_cell_markers %in% colnames(data)

intensity_matrix <- t(data[, c(chosen_cell_markers)])
colnames(intensity_matrix) <- paste("Cell_", seq(ncol(intensity_matrix)), sep="")

## Get x and y coords
coord_x <- data[ , "X"]
coord_y <- data[ , "Y"]
  
## Use SPIAT to get cell phenotypes
library(SPIAT)

# Format into spe object
general_format_image <- format_image_to_spe(format = "general",
                                            intensity_matrix = intensity_matrix,
                                            phenotypes = "unknown",
                                            coord_x = coord_x,
                                            coord_y = coord_y)

# Predict phenotypes
predicted_image <- predict_phenotypes(spe_object = general_format_image,
                                      thresholds = NULL,
                                      tumour_marker = "Keratin",
                                      baseline_markers = c("CD45", "CD31", "aSMA"),
                                      reference_phenotype = FALSE)

colData(predicted_image)
phenotype_df <- data.frame(table(predicted_image$Phenotype))

unique(predicted_image$Phenotype)
cell_types <- c("Other", "Tumour", "Immune",
                "Stroma", "Stroma", "Other",
                "Immune", "Other", "Tumour",
                "Other", "Stroma", "Stroma",
                "Other", "Other", "Other")


# Get cell types from phenotypes
formatted_image <- define_celltypes(
  predicted_image, 
  categories = unique(predicted_image$Phenotype), 
  category_colname = "Phenotype", 
  names = cell_types,
  new_colname = "Cell.Type"
)

chosen_cell_types <- c("Other", "Tumour", "Immune", "Stroma")

plot_data <- data.frame(spatialCoords(formatted_image))
plot_data$Cell.Type <- factor(formatted_image$Cell.Type, levels = chosen_cell_types)


## Plot data
library(ggplot2)
library(dplyr)
my_colours <- c("lightgray", "red", "blue", "darkgreen")

ggplot(sample_n(plot_data, 20000), aes(Cell.X.Position, Cell.Y.Position, colour = Cell.Type)) +
  geom_point() +
  scale_colour_manual(values = my_colours)




### Use Cell type dictionary to help annotate cell types ----------------------
setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/3D public spatial data/Lin et al - human colorectal cancer")
cell_type_dict <- read.csv("cell_type_dictionary.csv")

# format column names of dictionary
all_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                      "CD45RO", "CD4", "CD8a", "CD68", 
                      "CD163", "FOXP3", "PD1", "PDL1", 
                      "CD31", "aSMA", "Desmin", "CD45")
colnames(cell_type_dict)[1:length(all_cell_markers)] <- all_cell_markers 

# Get phenotype for each row
get_phenotype <- function(num_vec) {

  temp_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                        "CD45RO", "CD4", "CD8a", "CD68", 
                        "CD163", "FOXP3", "PD1", "PDL1", 
                        "CD31", "aSMA", "Desmin", "CD45")
  
  temp_cell_markers <- temp_cell_markers[num_vec == 1]
  phenotype <- paste(temp_cell_markers, collapse = ',')
  
  return(phenotype)
}

cell_type_dict$phenotype <- apply(cell_type_dict[, 1:length(all_cell_markers)], 1, get_phenotype)

### 3. Use cell type dictionary to determine cell identities (Tumour, Immune or Stroma, one slice) ------------------

data <- dataset[[1]]

# Get chosen cell markers
all_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                      "CD45RO", "CD4", "CD8a", "CD68", 
                      "CD163", "FOXP3", "PD1", "PDL1", 
                      "CD31", "aSMA", "Desmin", "CD45")

intensity_matrix <- t(data[, all_cell_markers])
colnames(intensity_matrix) <- paste("Cell_", seq(ncol(intensity_matrix)), sep="")

## Get x and y coords
coord_x <- data[ , "X"]
coord_y <- data[ , "Y"]

## Use SPIAT to get cell phenotypes
library(SPIAT)

# Format into spe object
general_format_image <- format_image_to_spe(format = "general",
                                            intensity_matrix = intensity_matrix,
                                            phenotypes = "unknown",
                                            coord_x = coord_x,
                                            coord_y = coord_y)

# Predict phenotypes
predicted_image <- predict_phenotypes(spe_object = general_format_image,
                                      thresholds = NULL,
                                      tumour_marker = "Keratin",
                                      baseline_markers = all_cell_markers[all_cell_markers != "Keratin"],
                                      reference_phenotype = FALSE)

colData(predicted_image)
phenotype_df <- data.frame(table(predicted_image$Phenotype))

phenotypes <- unique(predicted_image$Phenotype)

names <- ifelse(phenotypes %in% cell_type_dict$phenotype, cell_type_dict$Category, "Other")

# Get cell types from phenotypes
formatted_image <- define_celltypes(
  predicted_image, 
  categories = phenotypes, 
  category_colname = "Phenotype", 
  names = names,
  new_colname = "Cell.Type"
)

chosen_cell_types <- c("Other", "Tumor", "Immune", "Stroma")

plot_data <- data.frame(spatialCoords(formatted_image))
plot_data$Cell.Type <- factor(formatted_image$Cell.Type, levels = chosen_cell_types)


## Plot data
library(ggplot2)
library(dplyr)
my_colours <- c("lightgray", "red", "blue", "darkgreen")

ggplot(sample_n(plot_data, 5000), aes(Cell.X.Position, Cell.Y.Position, colour = Cell.Type)) +
  geom_point() +
  scale_colour_manual(values = my_colours)



### 4. Use cell type dictionary to determine cell identities (all cell types, one slice) ------------------

data <- dataset[[1]]

# Get chosen cell markers
all_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                      "CD45RO", "CD4", "CD8a", "CD68", 
                      "CD163", "FOXP3", "PD1", "PDL1", 
                      "CD31", "aSMA", "Desmin", "CD45")

intensity_matrix <- t(data[, all_cell_markers])
colnames(intensity_matrix) <- paste("Cell_", seq(ncol(intensity_matrix)), sep="")

## Get x and y coords
coord_x <- data[ , "X"]
coord_y <- data[ , "Y"]

## Use SPIAT to get cell phenotypes
library(SPIAT)

# Format into spe object
general_format_image <- format_image_to_spe(format = "general",
                                            intensity_matrix = intensity_matrix,
                                            phenotypes = "unknown",
                                            coord_x = coord_x,
                                            coord_y = coord_y)

# Predict phenotypes
predicted_image <- predict_phenotypes(spe_object = general_format_image,
                                      thresholds = NULL,
                                      tumour_marker = "Keratin",
                                      baseline_markers = all_cell_markers[all_cell_markers != "Keratin"],
                                      reference_phenotype = FALSE)

colData(predicted_image)

phenotypes <- unique(predicted_image$Phenotype)

names <- ifelse(phenotypes %in% cell_type_dict$phenotype, cell_type_dict$Type_Name, "Other")

# Get cell types from phenotypes
formatted_image <- define_celltypes(
  predicted_image, 
  categories = phenotypes, 
  category_colname = "Phenotype", 
  names = names,
  new_colname = "Cell.Type"
)


plot_data <- data.frame(spatialCoords(formatted_image))
plot_data$Cell.Type <- formatted_image$Cell.Type


## Plot data
library(ggplot2)
library(dplyr)

ggplot(sample_n(plot_data, 2000), aes(Cell.X.Position, Cell.Y.Position, colour = Cell.Type)) +
  geom_point()



### 5. Use cell type dictionary to determine cell identities (Tumour, Immune or Stroma, all slices) ---------

library(SPIAT)

# Get chosen cell markers
all_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                      "CD45RO", "CD4", "CD8a", "CD68", 
                      "CD163", "FOXP3", "PD1", "PDL1", 
                      "CD31", "aSMA", "Desmin", "CD45")


spatial_data_3D <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(spatial_data_3D) <- c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position", "Cell.Type")

## Loop through each slice

for (i in seq(length(dataset))) {
  data <- dataset[[i]]
  
  if (is.null(data)) {
    next
  }
  
  ## Get intensity matrix
  intensity_matrix <- t(data[, all_cell_markers])
  colnames(intensity_matrix) <- paste("Cell_", seq(ncol(intensity_matrix)), sep="")
  
  ## Get x, y and z coords
  coord_x <- data[ , "X"]
  coord_y <- data[ , "Y"]
  coord_z <- as.integer(slice_nums[i])
  coord_z <- rep(5 * coord_z, nrow(data)) # 5 microns between slices
  
  # Format into spe object
  general_format_image <- format_image_to_spe(format = "general",
                                              intensity_matrix = intensity_matrix,
                                              phenotypes = "unknown",
                                              coord_x = coord_x,
                                              coord_y = coord_y)
  
  # Predict phenotypes
  predicted_image <- predict_phenotypes(spe_object = general_format_image,
                                        thresholds = NULL,
                                        tumour_marker = "Keratin",
                                        baseline_markers = all_cell_markers[all_cell_markers != "Keratin"],
                                        reference_phenotype = FALSE)
  
  # Get cell type from phenotypes
  phenotypes <- unique(predicted_image$Phenotype)
  cell_types <- ifelse(phenotypes %in% cell_type_dict$phenotype, cell_type_dict$Category, "Other")
  
  all_phenotypes <- predicted_image$Phenotype
  all_cell_types <- cell_types[match(all_phenotypes, phenotypes)]
  
  # Add data to spatial_data_3D
  df <- data.frame(Cell.X.Position = coord_x,
                   Cell.Y.Position = coord_y,
                   Cell.Z.Position = coord_z,
                   Cell.Type = all_cell_types)
  
  spatial_data_3D <- rbind(spatial_data_3D, df)
}


### 5.1. Plotting multiple slices (tumour, immune, stroma) --------------------

setwd("C:/Users/Me/OneDrive - The University of Melbourne/PeterMac/Honours 2024/3D public spatial data/Lin et al - human colorectal cancer")
spatial_data_3D <- readRDS("spatial_data_3D.rds")

library(plotly)
library(dplyr)

plot_cell_categories3D <- function(data,
                                   cell_types_of_interest = NULL,
                                   colour_vector = NULL,
                                   size = 2,
                                   include_cell_types_of_no_interest = FALSE,
                                   feature_colname = "Cell.Type") {
  
  if (is.null(cell_types_of_interest)) {
    cell_types_of_interest <- unique(data$Cell.Type)
  }
  
  if (is.null(colour_vector)) {
    colour_vector <- hcl.colors(length(cell_types_of_interest), "Batlow")
  }
  
  if (length(cell_types_of_interest) != length(colour_vector)) {
    stop("Length of cell_types_of_interest is not equal to length of colour_vector")
  }
  
  ## Including non-interest cell types
  ## Define cell.id of non-interest cell types as "No Interest"
  cell_types_of_non_interest <- c()
  if (include_cell_types_of_no_interest) {
    cell_types_of_non_interest <- setdiff(unique(data[[feature_colname]]), cell_types_of_interest)
    
    data[data[[feature_colname]] %in% cell_types_of_non_interest, feature_colname] <- "No Interest"
    
    ## Add "No Interest" as a cell type of interest
    cell_types_of_interest <- c(cell_types_of_interest, "No Interest")
    
    ## Use lightgray for "No Interest" cell types
    colour_vector <- c(colour_vector, "#F0F0F0")
  }
  ## Excluding non-interest cell types
  ## Subset data to only include cell types of interest
  else {
    data <- data[data[[feature_colname]] %in% cell_types_of_interest, ]
  }
  
  ## Factor for feature column
  data[, feature_colname] <- factor(data[, feature_colname],
                                    levels = cell_types_of_interest)
  
  ## Plot
  fig <- plot_ly(data,
                 type = "scatter3d",
                 mode = 'markers',
                 x = ~Cell.X.Position,
                 y = ~Cell.Y.Position,
                 z = ~Cell.Z.Position,
                 color = ~Cell.Type,
                 colors = colour_vector,
                 marker = list(size = size))
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                     yaxis = list(title = 'y'),
                                     zaxis = list(title = 'z')))
  
  # fig <- fig %>% layout(scene = list(xaxis = list(title = '', showgrid = F, showaxeslabels = F, showticklabels = F),
  #                                    yaxis = list(title = '', showgrid = F, showaxeslabels = F, showticklabels = F),
  #                                    zaxis = list(title = '', showgrid = F, showaxeslabels = F, showticklabels = F)))
  
  return(fig)
  
}

my_colours <- c("lightgray", "red", "blue", "darkgreen")
chosen_cell_types <- c("Other", "Tumor", "Immune", "Stroma")

plot_cell_categories3D(sample_n(spatial_data_3D, 4000),
                       chosen_cell_types,
                       my_colours)

### 6. Use cell type dictionary to determine cell identities (all cell types, all slices)

