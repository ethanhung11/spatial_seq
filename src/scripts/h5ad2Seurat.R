# Rscript src/scripts/h5ad2Seurat.R &> ./outs/seurat_h5ad_conversion.out

library(here)
library(tictoc)
Sys.setenv(RETICULATE_PYTHON = here(".pixi/envs/main/bin/python"))
library(reticulate)
library(sceasy)
library(Seurat)

DATADIR <- here("data", "processed", "single_cell", "Emont2022")
obj_name <- "cleaned"

tic("converting to Seurat object")
seurat_object <- sceasy::convertFormat(
  here(DATADIR, paste0(obj_name,".h5ad")), 
  from = "anndata", 
  to = "seurat",
  outFile=here(DATADIR, paste0(obj_name,".rds"))
)
toc()