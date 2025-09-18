# general
options(max.print = 300)
options(timeout = 300)
options(stringsAsFactors = FALSE)
library(tidyverse)
library(here)
library(Seurat)
library(nichenetr)
source(here("src", "single_cell", "nichenet.R"))

# Prepare inputs
print("Read in seurat object...")
system.time(seurat_object <- readRDS(here("data","processed", "single cell", "5_analysis", "annotated.rds")))
celltype_col <- "cell_type_hires"

# Preprocess inputs
print("Start preparing inputs...")
nichenet.inputs <- prepNicheNet(seurat_object, celltype_col)

print("Pulling celltypes of interest...")
# celltypes of interest
celltypes <- c(nichenet.inputs$cell.types[grepl("fibro", nichenet.inputs$cell.types)],
               nichenet.inputs$cell.types[grepl("adipo", nichenet.inputs$cell.types)],
               nichenet.inputs$cell.types[grepl("macro", nichenet.inputs$cell.types)])
celltypes <- celltypes[!grepl("SAT", celltypes)]
nichenet_results <- list(celltypes=celltypes)

# run NicheNet & plots
for (celltype in celltypes) {
  print(paste("Beginning cell type:",celltype))
  
  nichenet_results[[celltype]] <- list()
  
  nichenet_results[[celltype]]$nichenet.outs <- runNicheNet(celltype, nichenet.inputs)
  
  nichenet_results[[celltype]]$agnostic.plots <- plotNicheNet_Agnostic(
    nichenet.inputs,
    nichenet_results[[celltype]]$nichenet.outs,
    n_ligands = 40)
  
  nichenet_results[[celltype]]$sender.plots <- plotNicheNet_SenderSpecific(
    nichenet.inputs,
    nichenet_results[[celltype]]$nichenet.outs,
    n_ligands = 20)
}

saveRDS(nichenet_results, here("data", "processed", "single cell", "5_analysis", "NicheNet_hires_results.rds"))
