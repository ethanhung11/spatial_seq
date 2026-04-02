# Rscript src/notebooks/pailin/tools/cellchat_NEURON_nonspatial.R &> ./outs/spatial/tools/cellchat_neuron_nonspatial.out
options(max.print = 300)
options(stringsAsFactors = FALSE)
library(here)
library(dplyr)
library(tictoc)
source(here("src", "single_cell", "cellchat.R"))

library(here)
Sys.setenv(RETICULATE_PYTHON = here("~/miniconda3/envs/bulk_rnaseq/bin/python"))

tic("start")
DATADIR <- here("data", "processed", "external")
CELLCHAT_DIR <- here(DATADIR, "tools", "cellchat")
dir.create(CELLCHAT_DIR, showWarnings = FALSE)
savename <- "cellchat_NEURON_nonspatial"
has.pathway <- TRUE

tic("read seurat object")
seurat_object <- readRDS(here(DATADIR, "analyzed.rds"))
seurat_object <- NormalizeData(seurat_object)
Idents(seurat_object) <- seurat_object$celltype
groups <- levels(seurat_object$timepoint)
toc()

tic("get NeuronChat database")
NeuronChatDB <- get_NeuronChatDB()
toc()

tic("run CellChat (all)")
cellchat <- prepCellChat(
  seurat_object,
  "celltype",
  db=NeuronChatDB)
cellchat <- subsetData(cellchat)
cellchat <- runCellChat(cellchat, pathway=has.pathway)

# renaming for custom database
int_names <- cellchat@LR$LRsig$interaction_name
dimnames(cellchat@net$prob)[[3]] <- int_names
dimnames(cellchat@net$pval)[[3]] <- int_names
names(cellchat@net$centr) <- int_names

saveRDS(cellchat, file = here(CELLCHAT_DIR, paste0(savename,"ALL.rds")))
toc()
print("\n")

# run on timepoint data
for (g in groups) {
  tic(paste("run CellChat", g))
  
  subset_seurat <- subset(seurat_object, subset = timepoint == g)
  cellchat <- prepCellChat(
    subset_seurat,
    "celltype",
    db=NeuronChatDB)
  cellchat <- subsetData(cellchat)
  
  cellchat <- runCellChat(cellchat, pathway=has.pathway)
  
  # renaming for custom database
  int_names <- cellchat@LR$LRsig$interaction_name
  dimnames(cellchat@net$prob)[[3]] <- int_names
  dimnames(cellchat@net$pval)[[3]] <- int_names
  names(cellchat@net$centr) <- int_names
  
  saveRDS(cellchat, file = here(CELLCHAT_DIR, paste0(savename, "_", g, ".rds")))
  
  toc()
  print("\n")
}

tic("run Merged CellChat")
cellchat.list <- lapply(
  here(CELLCHAT_DIR,paste0(savename, "_", groups, ".rds")),
  function(x) readRDS(x)
)
names(cellchat.list) <- groups

# CUSTOM - Merge all neurons
grouped_celltypes <- levels(cellchat.list[[1]]@idents)
grouped_celltypes[grepl("neu", grouped_celltypes, ignore.case = TRUE)] <- "Neurons"
grouped_celltypes <- factor(grouped_celltypes)
cellchat.list <- lapply(cellchat.list, function(x) {mergeInteractions(x, grouped_celltypes)})
for (g in groups) saveRDS(cellchat.list[[g]], file = here(CELLCHAT_DIR, paste0(savename, "_", g, ".rds")))

# Merge
cellchat.MERGED <- mergeCellChat.custom(cellchat.list, pathway=has.pathway)
saveRDS(cellchat.MERGED, file = here(CELLCHAT_DIR, paste0(savename, "_MERGED.rds")))
toc()

toc()