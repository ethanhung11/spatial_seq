# general
library(dplyr)
library(here)

# single-cell general
library(BiocParallel)
library(scater)
library(Seurat)
library(convert2anndata)
library(terra)
options(Seurat.object.assay.version = "v3")

# others
library(scry)
library(scDblFinder)
library(DoubletFinder)
library("loupeR")

options(max.print = 300)

rds_filename <- "CU048_ConfPositioned_seurat_spatial_merged.rds"
seurat.obj <- readRDS(file = here("data", "processed", "external", rds_filename))

View(seurat.obj)

plot <- SpatialDimPlot(seurat.obj, label = TRUE)
cropped.coords <- Crop(vizgen.obj[["slice1"]], x = c(1750, 3000), y = c(3750, 5250), coords = "plot")
