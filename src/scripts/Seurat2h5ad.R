library(sceasy)
library(reticulate)
library(Seurat)
library(here)
DATADIR <- here("data", "processed", "external")
obj_name <- "analyzed"

seurat_object <- sceasy::convertFormat(
  here(DATADIR, paste0(obj_name,".rds")), 
  from = "seurat", 
  to = "seurat",
  outFile=here(DATADIR, paste0(obj_name,".h5ad")))