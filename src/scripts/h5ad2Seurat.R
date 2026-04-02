# Rscript src/scripts/h5ad2Seurat.R &> ./outs/seurat_h5ad_conversion.out

library(here)
Sys.setenv(RETICULATE_PYTHON = here(".pixi/envs/main/bin/python"))
library(tictoc)
library(reticulate)
library(sceasy)
library(Seurat)

DATADIR <- here("data", "processed", "single_cell", "combined")
obj_name <- "eWAT_Male-fb_macs"

tic("converting to Seurat object")
seurat_object <- sceasy::convertFormat(
  here(DATADIR, paste0(obj_name,".h5ad")), 
  from = "anndata", 
  to = "seurat",
  outFile=here(DATADIR, paste0(obj_name,".rds"))
)
toc()