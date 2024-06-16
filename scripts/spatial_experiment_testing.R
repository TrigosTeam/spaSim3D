library(SpatialExperiment)
library(DropletUtils)

### 1. Class structure -----------------------------------------------
# Get example spe object
example(read10xVisium, echo = FALSE)
spe

# Get spatial coords (stored in int_colData)
head(spatialCoords(spe))

# Get image data (stored in int_metadata's imgData field)
spi <- imgData(spe)

# Show image data
plot(imgRaster(spe))



### 2. Object construction ------------------------------
# spatialCoords: numeric matrix with spatial coords
# spatialCoordsNames: character vector

# Method 1. Supply spatialCoords via colData
n <- length(z <- letters)
y <- matrix(nrow = n, ncol = n)
cd <- DataFrame(x = seq(n), y = seq(n), z)

spe1 <- SpatialExperiment(
  assay = y,
  colData = cd,
  spatialCoordsNames = c("x", "y")
)


# Method 2. Supply spatialCoords with a matrix
xy <- as.matrix(cd[ , c("x", "y")])

spe2 <- SpatialExperiment(
  assay = y,
  colData = cd["z"],
  spatialCoords = xy
)

identical(spe1, spe2)


# Using 3D coords???
n <- length(reference <- letters)
y <- matrix(nrow = n, ncol = n)
cd <- DataFrame(x = seq(n), y = seq(n), z = seq(n), reference)

spe3D <- SpatialExperiment(
  assay = y,
  colData = cd,
  spatialCoordsNames = c("x", "y", "z")
)




# Spot-based ST data
dir <- system.file(
  file.path("extdata", "10xVisium", "section1", "outs"),
  package = "SpatialExperiment")

# read in counts
fnm <- file.path(dir, "raw_feature_bc_matrix")
sce <- DropletUtils::read10xCounts(fnm)

# read in image data
img <- readImgData(
  path = file.path(dir, "spatial"),
  sample_id = "foo")

# read in spatial coordinates
fnm <- file.path(dir, "spatial", "tissue_positions_list.csv")
xyz <- read.csv(fnm, header = FALSE,
                col.names = c(
                  "barcode", "in_tissue", "array_row", "array_col",
                  "pxl_row_in_fullres", "pxl_col_in_fullres"))

# construct observation & feature metadata
rd <- S4Vectors::DataFrame(
  symbol = rowData(sce)$Symbol)

# construct 'SpatialExperiment'
(spe <- SpatialExperiment(
  assays = list(counts = assay(sce)),
  rowData = rd, 
  colData = DataFrame(xyz), 
  spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
  imgData = img,
  sample_id = "foo"))




### 3. Common operations ----------------------------------------------------

# Subsetting by sample_id
sub <- spe[ , spe$sample_id == "fpo"]


# Subsetting to keep observations only in tissue
sub <- spe[ , colData(spe)$in_tissue == 1]





