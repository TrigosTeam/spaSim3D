### Read updated data --------------------------------------------------------

setwd("~/Lin et al - human colorectal cancer/CRC1_data_updated")

slice_nums <- c("002", "007", "014", "020", "025", "029",
                "034", "039", "044", "049", "050", "051",
                "052", "054", "059", "064", "069", "074",
                "078", "084", "086", "091", "097")
file_names <- paste("CRC01-", slice_nums, ".csv", sep = "")

# Get z coord of each slice in micrometers using the slice number
get_slice_z_coord <- function(slice_num) {
  
  slice_num <- as.integer(slice_num)
  
  if (slice_num <= 85) {
    return(slice_num * 5)
  }
  
  else {
    return(85 * 5 + (slice_num - 85) * 4)
  }
}

# Store combined result in data frame
df <- data.frame()

# Loop through each slice
for (i in seq(length(slice_nums))) {
  
  file_name <- file_names[i]
  
  CRC1_df <- read.csv(file_name, sep = ",")
  CRC1_df <- CRC1_df[ , c("NewType", "Xtt", "Ytt")]
  CRC1_df$Ztt <- get_slice_z_coord(slice_nums[i])
  
  df <- rbind(df, CRC1_df)

}

# Current data frame has a number to represent each cell type
# Use the cell type dictionary to determine the cell type from the number

setwd("~/Lin et al - human colorectal cancer/Other_data")

cell_type_dictionary <- read.csv("cell_type_dictionary.csv")

specific_cell_types <- unique(cell_type_dictionary$Type_Name) # E.g. Tumour/Epi, Ki67+ Tumour/Epi...
df$Cell.Type.Specific <- specific_cell_types[df$NewType]

generic_cell_types <- unique(cell_type_dictionary$Category) # E.g. Tumour, Stroma, Immune
df$Cell.Type.Generic <- generic_cell_types[ifelse(df$NewType <= 3, 1, ifelse(df$NewType <= 6, 2, 3))]

df$NewType <- NULL

# Rename coordinate columns
colnames(df)[1:3] <- c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position")

# Save
# setwd("~/Lin et al - human colorectal cancer/CRC1_data_updated")
# saveRDS(df, "CRC1_df.rds")




### Plot data ----------------------------------------------------------------

setwd("~/Lin et al - human colorectal cancer/CRC1_data_updated")
df <- readRDS("CRC1_df.rds")

library(dplyr)
library(rgl)
plot_cells_rgl_3D <- function(df, n_cells, plot_cell_types = NULL, feature_colname) {
  
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
                    
                    "green1", "green", "green3", "green4", "darkgreen", # Macrophages
                    "purple1",  "purple2", "purple3", "purple4",    # Lymphocytes
                    "turquoise1", "turquoise2", "steelblue1", "steelblue2", "steelblue3", # T cells
                    "orchid",                 # B cells
                    
                    "lightgray")
  
  names(cell_colours) <- cell_types
  
  if (!is.null(plot_cell_types)) df <- df[df[[feature_colname]] %in% plot_cell_types, ]
  
  df_plot <- sample_n(df, n_cells)  
  df_plot$Cell.Colour <- cell_colours[df_plot[[feature_colname]]]
  
  options(rgl.printRglwidget = T)
  
  open3d()
  plot3d(df_plot$Cell.X.Position,
         df_plot$Cell.Y.Position,
         df_plot$Cell.Z.Position,
         col = df_plot$Cell.Colour,
         size = 4,
         xlab = 'x',
         ylab = 'y',
         zlab = 'z',
         xlim = NULL,
         ylim = NULL,
         zlim = c(0, 500),
         forceClipregion = TRUE)
  aspect3d(5, 5, 1)
  highlevel()
  
}

# Tumour cells
tumour_cell_types <- c("Tumor/Epi.", "Ki67+ Tumor/Epi.", "PDL1+ Tumor/Epi.")

# Stromal cells
stromal_cell_types <- c("Endothelial", "Muscle/Fibroblast")

# Immune cells
immune_cell_types <- c("Macrophage(I)", "Macrophage(II)", "Macrophage(III)", "Macrophage(IV)", "PDL1+ Macrophage",
                       "PDL1+ lymphocyte",  "DN Lymphocyte", "DP Lymphocyte", "Lymphocyte(III)",
                       "T helper", "PD1+ T helper", "Tc cell", "PD1+ Tc", "Treg",
                       "B cells")


plot_cells_rgl_3D(df, 50000, NULL, "Cell.Type.Specific")
plot_cells_rgl_3D(df, 50000, tumour_cell_types, "Cell.Type.Specific")
plot_cells_rgl_3D(df, 50000, stromal_cell_types, "Cell.Type.Specific")
plot_cells_rgl_3D(df, 50000, immune_cell_types, "Cell.Type.Specific")


plot_cells_rgl_3D(df, 50000, "Tumor/Epi.", "Cell.Type.Specific")


