# Rscript src/scripts/h5ad2Seurat.R &> ./outs/seurat_h5ad_conversion.out

library(here)
library(tictoc)
Sys.setenv(RETICULATE_PYTHON = here(".pixi/envs/main/bin/python"))
library(reticulate)
library(sceasy)
library(Seurat)

DATADIR <- here("data", "processed", "spatial", "Xenium", "kennedi_flu")
obj_name <- "integrated_slim"
obj_name_new <- "integrated"

tic("converting to Seurat object")
seurat_object <- sceasy::convertFormat(
  here(DATADIR, paste0(obj_name,".h5ad")), 
  from = "anndata", 
  to = "seurat",
  outFile=here(DATADIR, paste0(obj_name_new,".rds"))
)
toc()