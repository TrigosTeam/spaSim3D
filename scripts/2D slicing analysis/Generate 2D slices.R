### Mixed spes --------------------------------------------------------------

# Get mixed_spes_table
setwd("~/Objects/spes_table")
mixed_spes_table <- read.table("mixed_spes_table.csv")

# Get number of mixed spes
n_mixed_spes <- nrow(mixed_spes_table)

# Define slice parameters
n_slices <- 5
slice_width <- 10
first_slice_bottom_z_coord <- 25

slice_z_coords <- data.frame(slice_number = seq(n_slices),
                             bottom = seq(first_slice_bottom_z_coord, first_slice_bottom_z_coord + (n_slices - 1) * slice_width, slice_width))
slice_z_coords$top <- slice_z_coords$bottom + slice_width


# Get slice data for each mixed_spe
for (i in seq(n_mixed_spes)) {
  setwd("~/Objects/mixed_spes")
  
  # Read in current mixed spe
  mixed_spe_name <- paste("mixed_spe_", i, sep = "")
  mixed_spe_file_name <- paste(mixed_spe_name, ".rds", sep = "")
  mixed_spe <- readRDS(mixed_spe_file_name)
  
  # Get z-coords from current mixed_spe
  mixed_spe_z_coords <- spatialCoords(mixed_spe)[ , "Cell.Z.Position"]
  
  # Define list to contain all 2D slice spes
  mixed_spe_slices <- list()
  
  for (j in slice_z_coords$slice_number) {
    
    # Get z-coord boundaries for current slice
    bottom_z_coord <- slice_z_coords[j, "bottom"]
    top_z_coord <- slice_z_coords[j, "top"]
    
    # Get spe slice
    mixed_spe_slice <- mixed_spe[ , bottom_z_coord < mixed_spe_z_coords & mixed_spe_z_coords < top_z_coord]
    
    # Remove the z-coords
    mixed_spe_slice@int_colData@listData$spatialCoords <- mixed_spe_slice@int_colData@listData$spatialCoords[ , c("Cell.X.Position", "Cell.Y.Position")]
    
    mixed_spe_slices[[j]] <- mixed_spe_slice
  }
  setwd("~/Objects/mixed_spes/2D_slices")
  mixed_spe_slices_file_name <- paste("mixed_spe_slices_", i, ".rds", sep = "")
  saveRDS(mixed_spe_slices, mixed_spe_slices_file_name)
}



### Ringed spes --------------------------------------------------------------

# Get ringed_spes_table
setwd("~/Objects/spes_table")
ringed_spes_table <- read.table("ringed_spes_table.csv")

# Get number of ringed spes
n_ringed_spes <- nrow(ringed_spes_table)

# Define slice parameters
n_slices <- 5
slice_width <- 10
first_slice_bottom_z_coord <- 25

slice_z_coords <- data.frame(slice_number = seq(n_slices),
                             bottom = seq(first_slice_bottom_z_coord, first_slice_bottom_z_coord + (n_slices - 1) * slice_width, slice_width))
slice_z_coords$top <- slice_z_coords$bottom + slice_width


# Get slice data for each ringed_spe
for (i in seq(n_ringed_spes)) {
  setwd("~/Objects/ringed_spes")
  
  # Read in current ringed spe
  ringed_spe_name <- paste("ringed_spe_", i, sep = "")
  ringed_spe_file_name <- paste(ringed_spe_name, ".rds", sep = "")
  ringed_spe <- readRDS(ringed_spe_file_name)
  
  # Get z-coords from current ringed_spe
  ringed_spe_z_coords <- spatialCoords(ringed_spe)[ , "Cell.Z.Position"]
  
  # Define list to contain all 2D slice spes
  ringed_spe_slices <- list()
  
  for (j in slice_z_coords$slice_number) {
    
    # Get z-coord boundaries for current slice
    bottom_z_coord <- slice_z_coords[j, "bottom"]
    top_z_coord <- slice_z_coords[j, "top"]
    
    # Get spe slice
    ringed_spe_slice <- ringed_spe[ , bottom_z_coord < ringed_spe_z_coords & ringed_spe_z_coords < top_z_coord]
    
    # Remove the z-coords
    ringed_spe_slice@int_colData@listData$spatialCoords <- ringed_spe_slice@int_colData@listData$spatialCoords[ , c("Cell.X.Position", "Cell.Y.Position")]
    
    ringed_spe_slices[[j]] <- ringed_spe_slice
  }
  setwd("~/Objects/ringed_spes/2D_slices")
  ringed_spe_slices_file_name <- paste("ringed_spe_slices_", i, ".rds", sep = "")
  saveRDS(ringed_spe_slices, ringed_spe_slices_file_name)
}


### Separated spes --------------------------------------------------------------

# Get separated_spes_table
setwd("~/Objects/spes_table")
separated_spes_table <- read.table("separated_spes_table.csv")

# Only focus on separated simulations which use the same shapes and sizes for simplicity
separated_spes_table <- separated_spes_table[separated_spes_table$shapeA == separated_spes_table$shapeB &
                                               separated_spes_table$sizeA == separated_spes_table$sizeB, ]

# Get separated spes simulation numbers
separated_spes_numbers <- as.integer(rownames(separated_spes_table))

# Define slice parameters
n_slices <- 5
slice_width <- 10
first_slice_bottom_z_coord <- 25

slice_z_coords <- data.frame(slice_number = seq(n_slices),
                             bottom = seq(first_slice_bottom_z_coord, first_slice_bottom_z_coord + (n_slices - 1) * slice_width, slice_width))
slice_z_coords$top <- slice_z_coords$bottom + slice_width


# Get slice data for each separated_spe
for (i in separated_spes_table) {
  setwd("~/Objects/separated_spes")
  
  # Read in current separated spe
  separated_spe_name <- paste("separated_spe_", i, sep = "")
  separated_spe_file_name <- paste(separated_spe_name, ".rds", sep = "")
  separated_spe <- readRDS(separated_spe_file_name)
  
  # Get z-coords from current separated_spe
  separated_spe_z_coords <- spatialCoords(separated_spe)[ , "Cell.Z.Position"]
  
  # Define list to contain all 2D slice spes
  separated_spe_slices <- list()
  
  for (j in slice_z_coords$slice_number) {
    
    # Get z-coord boundaries for current slice
    bottom_z_coord <- slice_z_coords[j, "bottom"]
    top_z_coord <- slice_z_coords[j, "top"]
    
    # Get spe slice
    separated_spe_slice <- separated_spe[ , bottom_z_coord < separated_spe_z_coords & separated_spe_z_coords < top_z_coord]
    
    # Remove the z-coords
    separated_spe_slice@int_colData@listData$spatialCoords <- separated_spe_slice@int_colData@listData$spatialCoords[ , c("Cell.X.Position", "Cell.Y.Position")]
    
    separated_spe_slices[[j]] <- separated_spe_slice
  }
  setwd("~/Objects/separated_spes/2D_slices")
  separated_spe_slices_file_name <- paste("separated_spe_slices_", i, ".rds", sep = "")
  saveRDS(separated_spe_slices, separated_spe_slices_file_name)
}




### -------------------------------------------------------------------------