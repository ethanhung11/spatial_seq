# Rscript src/scripts/Seurat2h5ad.R &> ./outs/seurat_h5ad_conversion.out

library(here)
library(tictoc)
Sys.setenv(RETICULATE_PYTHON = here(".pixi/envs/main/bin/python"))
library(reticulate)
library(Seurat)
library(scCustomize)

DATADIR <- here("references")
obj_name <- "FluSobj"

tic("converting to h5ad object")
obj <- readRDS(here(DATADIR, paste0(obj_name,".rds")))
print(Assays(obj))
print(Layers(obj))
as.anndata(x = obj,
  main_layer = "counts", other_layers = NULL,
  file_path = DATADIR, file_name = paste0(obj_name,".h5ad")
)
toc()