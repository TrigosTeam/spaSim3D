### Get cell type dictionary to help annotate cell types ----------------------
setwd("~/Lin et al - human colorectal cancer/Other data")

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

### 1a. Combine each slice together, then determine cell types ---------

library(SPIAT)

### Get all slice data
setwd("~/Lin et al - human colorectal cancer/CRC1 data")

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



## Loop through each slice and combine all slice data first (also add the z coord)
combined_data <- data.frame()

for (i in seq(length(dataset))) {
  data <- dataset[[i]]
  data$Z <- as.integer(slice_nums[i]) * 5
  
  combined_data <- rbind(combined_data, data)
}

# Get chosen cell markers
all_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                      "CD45RO", "CD4", "CD8a", "CD68", 
                      "CD163", "FOXP3", "PD1", "PDL1", 
                      "CD31", "aSMA", "Desmin", "CD45")

## Get intensity matrix
intensity_matrix <- t(combined_data[, all_cell_markers])
colnames(intensity_matrix) <- paste("Cell_", seq(ncol(intensity_matrix)), sep="")

## Get x, y and z coords
coord_x <- combined_data[ , "X"]
coord_y <- combined_data[ , "Y"]
coord_z <- combined_data[ , "Z"]

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

all_phenotypes <- predicted_image$Phenotype
phenotypes <- unique(all_phenotypes)

# Get cell types from phenotypes
cell_types1 <- ifelse(phenotypes %in% cell_type_dict$phenotype, cell_type_dict$Category, "Other") ## Tumour, immune, stromal
cell_types2 <- ifelse(phenotypes %in% cell_type_dict$phenotype, cell_type_dict$Type_Name, "Other") ## All cell types

all_cell_types1 <- cell_types1[match(all_phenotypes, phenotypes)]
all_cell_types2 <- cell_types2[match(all_phenotypes, phenotypes)]

# Get final data frame
df <- data.frame(Cell.X.Position = coord_x,
                 Cell.Y.Position = coord_y,
                 Cell.Z.Position = coord_z,
                 Cell.Type1 = all_cell_types1,
                 Cell.Type2 = all_cell_types2)

### 1b. Plot (tumour, immune, stromal) --------------------

setwd("~/Lin et al - human colorectal cancer/Other data")
df <- readRDS("lin_et_al_3D_spatial_data_combine_slices.rds")

## Get all unique cell types
cell_types <- c("Tumor", "Immune", "Stroma", "Other")

## Assign colour to each cell type
cell_colours <- c("orange", "skyblue2", "brown", "lightgray")

names(cell_colours) <- cell_types
df$Cell.Colour <- cell_colours[df$Cell.Type1]


## Plot
library(dplyr)
library(rgl)
options(rgl.printRglwidget = T)

df_plot <- sample_n(df, 100000)
plot3d(df_plot$Cell.X.Position,
       df_plot$Cell.Y.Position,
       df_plot$Cell.Z.Position,
       col = df_plot$Cell.Colour,
       size = 8,
       xlab = 'x',
       ylab = 'y',
       zlab = 'z')

# legend3d("topright", legend = cell_types, pch = 16, col = cell_colours, cex=0.5, inset=0.02)



### 1c. Plot (all cell types) --------------------

setwd("~/Lin et al - human colorectal cancer/Other data")
df <- readRDS("lin_et_al_3D_spatial_data_combine_slices.rds")

## Get all unique cell types
cell_types <- c("Tumor/Epi.", "Ki67+ Tumor/Epi.", "PDL1+ Tumor/Epi.", 
                
                "Endothelial", "Muscle/Fibroblast",      
                
                "Macrophage(I)", "Macrophage(II)", "Macrophage(III)", "Macrophage(IV)", "PDL1+ Macrophage", 
                "PDL1+ lymphocyte",  "DN Lymphocyte", "DP Lymphocyte", "Lymphocyte(III)",    
                "T helper", "PD1+ T helper", "Tc cell", "PD1+ Tc", "Treg", 
                "B cells",
                
                "Other")

## Assign colour to each cell type
cell_colours <- c("orange", "orange2", "orange3",       # Tumour
                
                "brown1", "brown",             # Endothelial & fibroblast
                
                "green", "green1", "green2", "green3", "green4", # Macrophages
                "purple1",  "purple2", "purple3", "purple4",    # Lymphocytes
                "skyblue", "skyblue1", "skyblue2", "skyblue3", "skyblue4", # T cells
                "orchid",                 # B cells
                
                "lightgray")

names(cell_colours) <- cell_types
df$Cell.Colour <- cell_colours[df$Cell.Type2]


## Plot
library(dplyr)
library(rgl)
options(rgl.printRglwidget = T)

df_plot <- sample_n(df, 100000)
plot3d(df_plot$Cell.X.Position,
       df_plot$Cell.Y.Position,
       df_plot$Cell.Z.Position,
       col = df_plot$Cell.Colour,
       size = 8,
       xlab = 'x',
       ylab = 'y',
       zlab = 'z')

# legend3d("topright", legend = cell_types, pch = 16, col = cell_colours, cex=0.5, inset=0.02)








### 3a. Examine each slice separately to determine cell types


### 2a. Examine each slice separately and determine cell types for each slice ---------

library(SPIAT)

### Get all slice data
setwd("~/Lin et al - human colorectal cancer/CRC1 data")

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



# Get chosen cell markers
all_cell_markers <- c("Keratin", "Ki67", "CD3", "CD20", 
                      "CD45RO", "CD4", "CD8a", "CD68", 
                      "CD163", "FOXP3", "PD1", "PDL1", 
                      "CD31", "aSMA", "Desmin", "CD45")


## Define result data.frame
df <- data.frame()


## Loop through each slice
i <- 1
for (data in dataset) {
  
  ## Get intensity matrix
  intensity_matrix <- t(data[, all_cell_markers])
  colnames(intensity_matrix) <- paste("Cell_", seq(ncol(intensity_matrix)), sep="")
  
  ## Get x, y and z coords
  coord_x <- data[ , "X"]
  coord_y <- data[ , "Y"]
  coord_z <- rep(as.integer(slice_nums[i]) * 5, nrow(data))
  
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
  
  all_phenotypes <- predicted_image$Phenotype
  phenotypes <- unique(all_phenotypes)
  
  # Get cell types from phenotypes
  cell_types1 <- ifelse(phenotypes %in% cell_type_dict$phenotype, cell_type_dict$Category, "Other") ## Tumour, immune, stromal
  cell_types2 <- ifelse(phenotypes %in% cell_type_dict$phenotype, cell_type_dict$Type_Name, "Other") ## All cell types
  
  all_cell_types1 <- cell_types1[match(all_phenotypes, phenotypes)]
  all_cell_types2 <- cell_types2[match(all_phenotypes, phenotypes)]
  
  df_temp <- data.frame(Cell.X.Position = coord_x,
                        Cell.Y.Position = coord_y,
                        Cell.Z.Position = coord_z,
                        Cell.Type1 = all_cell_types1,
                        Cell.Type2 = all_cell_types2)
  
  df <- rbind(df, df_temp)
  
  i <- i + 1
}


### 2b. Plot (tumour, immune, stromal) --------------------

setwd("~/Lin et al - human colorectal cancer/Other data")
df <- readRDS("lin_et_al_3D_spatial_data_separate_slices.rds")

## Get all unique cell types
cell_types <- c("Tumor", "Immune", "Stroma", "Other")

## Assign colour to each cell type
cell_colours <- c("orange", "skyblue2", "brown", "lightgray")

names(cell_colours) <- cell_types
df$Cell.Colour <- cell_colours[df$Cell.Type1]


## Plot
library(dplyr)
library(rgl)
options(rgl.printRglwidget = T)

df_plot <- sample_n(df, 100000)
plot3d(df_plot$Cell.X.Position,
       df_plot$Cell.Y.Position,
       df_plot$Cell.Z.Position,
       col = df_plot$Cell.Colour,
       size = 8,
       xlab = 'x',
       ylab = 'y',
       zlab = 'z')

# legend3d("topright", legend = cell_types, pch = 16, col = cell_colours, cex=0.5, inset=0.02)



### 2c. Plot (all cell types) --------------------

setwd("~/Lin et al - human colorectal cancer/Other data")
df <- readRDS("lin_et_al_3D_spatial_data_separate_slices.rds")

## Get all unique cell types
cell_types <- c("Tumor/Epi.", "Ki67+ Tumor/Epi.", "PDL1+ Tumor/Epi.", 
                
                "Endothelial", "Muscle/Fibroblast",      
                
                "Macrophage(I)", "Macrophage(II)", "Macrophage(III)", "Macrophage(IV)", "PDL1+ Macrophage", 
                "PDL1+ lymphocyte",  "DN Lymphocyte", "DP Lymphocyte", "Lymphocyte(III)",    
                "T helper", "PD1+ T helper", "Tc cell", "PD1+ Tc", "Treg", 
                "B cells",
                
                "Other")

## Assign colour to each cell type
cell_colours <- c("orange", "orange2", "orange3",       # Tumour
                  
                  "brown1", "brown",             # Endothelial & fibroblast
                  
                  "green", "green1", "green2", "green3", "green4", # Macrophages
                  "purple1",  "purple2", "purple3", "purple4",    # Lymphocytes
                  "skyblue", "skyblue1", "skyblue2", "skyblue3", "skyblue4", # T cells
                  "orchid",                 # B cells
                  
                  "lightgray")

names(cell_colours) <- cell_types
df$Cell.Colour <- cell_colours[df$Cell.Type2]


## Plot
library(dplyr)
library(rgl)
options(rgl.printRglwidget = T)

df_plot <- sample_n(df, 100000)
plot3d(df_plot$Cell.X.Position,
       df_plot$Cell.Y.Position,
       df_plot$Cell.Z.Position,
       col = df_plot$Cell.Colour,
       size = 8,
       xlab = 'x',
       ylab = 'y',
       zlab = 'z')

# legend3d("topright", legend = cell_types, pch = 16, col = cell_colours, cex=0.5, inset=0.02)








### 3a. Examine each slice separately to determine cell types