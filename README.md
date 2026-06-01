
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spaSim3D

<!-- badges: start -->

<!-- badges: end -->

The goal of spaSim3D (**spa**tial **Sim**ulator 3D) is to generate fully
customisable 3D spatial tissue data. It includes a diverse range of
functions to meet your simulation needs. This includes generation and
customisation of background cells in the 3D tissue. This also includes
generation and customisation of various cell clusters, such as spheres,
ellipsoids, cylinders and network clusters with or without cellular
rings and double rings. Collectively, these tools can help to mimic a
realistic 3D biological tissue.

# Installation

You can install the development version of spaSim3D from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("TrigosTeam/spaSim3D")
```

# Vignette

The vignette with an overview of the package can be accessed from the
top Menu under About or clicking here.

# Example

This is a basic example which shows how to simulate background cells
with different cell clusters. I can’t show the plot because it is an
interactive 3D plot…

``` r
library(spaSim3D)

spe_metadata <- spe_metadata_background_template("random")
spe_metadata <- spe_metadata_cluster_template("regular", "sphere", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("regular", "network", spe_metadata)
spe_metadata <- spe_metadata_cluster_template("ring", "ellipsoid", spe_metadata)
spe_clusters <- simulate_spe_metadata3D(spe_metadata,
                                        plot_image = F)
# Plot
# plot_cells3D(spe_clusters,
#              plot_cell_types = c("Tumour", "Immune", "Immune1"),
#              plot_colours = c("orange", "skyblue", "blue"))
```
