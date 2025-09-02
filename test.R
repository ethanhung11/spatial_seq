# general
options(max.print = 300)
library(dplyr)
library(here)

# single-cell general
library(BiocParallel)
library(scRNAseq)
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
library(scSHC)
library(CHOIR)

cluster_names <- c('leiden_0.5','leiden_0.6',
                  'leiden_0.7','leiden_0.8',
                  'leiden_0.9','leiden_1.0',
                  'leiden_1.1','leiden_1.2',
                  'leiden_1.3','leiden_1.4',
                  'leiden_1.5','leiden')

savedir <- here("data", "processed", "single cell", "3_annotated", "manDoublet-seuratV3-harmony-annotated.rds")
seurat.obj <- readRDS(savedir)
# need UMAP/LocalMAP dim red
# need nearest neighbors

# CHOIR
inferTree(seurat.obj, )


# scSHC
tmp <- testClusters(
  seurat.obj[["RNA"]]@counts,
  cluster_ids=unlist(seurat.obj[[cluster_names[1]]]),
  cores=10
)
saveRDS(tmp, here("data", "processed", "single cell", "3_annotated", "scSHC-leiden0.5.rds"))
View(tmp[[1]])



# library(SeuratData)
# library(scSHC)
# AvailableData()
# data("pbmc3k")
# pbmc3k <- UpdateSeuratObject(pbmc3k)
# print(dim(pbmc3k[["RNA"]]$counts))
# print(length(pbmc3k$seurat_annotations))
# tmp <- testClusters(
#   pbmc3k[["RNA"]]$counts,
#   as.character(pbmc3k$seurat_annotations),
#   cores=10
# )