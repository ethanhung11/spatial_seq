# Rscript src/scripts/Seurat2h5ad.R &> ./outs/seurat_h5ad_conversion.out

library(here)
Sys.setenv(RETICULATE_PYTHON = here(".pixi/envs/main/bin/python"))
library(tictoc)
library(reticulate)
library(sceasy)
library(Seurat)

DATADIR <- here("data", "processed", "external")
obj_name <- "analyzed"

tic("converting to h5ad object")
seurat_object <- sceasy::convertFormat(
  here(DATADIR, paste0(obj_name,".rds")), 
  from = "seurat", 
  to = "anndata",
  outFile=here(DATADIR, paste0(obj_name,".h5ad")))
toc()