# Rscript src/notebooks/pailin/tools/cellchat_NEURON.R &> ./outs/spatial/tools/cellchat_neuron.out
options(max.print = 300)
options(stringsAsFactors = FALSE)
library(here)
library(dplyr)
library(tictoc)
source(here("src", "single_cell", "cellchat.R"))

tic("start")
DATADIR <- here("data", "processed", "external")
CELLCHAT_DIR <- here(DATADIR, "tools", "cellchat")
dir.create(CELLCHAT_DIR, showWarnings = FALSE)
savename <- "cellchat_NEURON"
has.pathway <- FALSE

tic("read seurat object")
seurat_object <- readRDS(here(DATADIR, "analyzed.rds"))
seurat_object <- NormalizeData(seurat_object)
Idents(seurat_object) <- seurat_object$celltype
groups <- levels(seurat_object$timepoint)
toc()

tic("pull NeuronChat database")
load(url("https://github.com/Wei-BioMath/NeuronChat/raw/main/data/interactionDB_mouse.rda"))
neuron_df <- as.data.frame(do.call(rbind, lapply(interactionDB_mouse, function(x) {
  data.frame(
    interaction_name = x$interaction_name,
    ligand = x$lig_contributor,
    receptor = x$receptor_subunit,
    lig_contributor = x$lig_contributor,
    receptor_subunit = x$receptor_subunit,
    interaction_type = x$interaction_type,
    ligant_type = x$ligand_type,
    stringsAsFactors = FALSE
  )
})))
neuronChatDB <- customCellChatDB(neuron_df, evidence = "NeuronChat")
toc()
print("\n")

tic("run CellChat (all)")
cellchat <- prepCellChat(
  seurat_object,
  "celltype",
  db=neuronChatDB, # input for custom DB
  spatial=TRUE,
  spatial_key="spatial",
  spatial_factors=data.frame(ratio = 4/3, tol = 5)) # 6mm/6000um vs ~4500px tissue width 
cellchat <- subsetData(cellchat)
cellchat <- runCellChat(cellchat, spatial=TRUE, pathway=has.pathway)
saveRDS(cellchat, file = here(CELLCHAT_DIR, paste0(savename,"ALL.rds")))
toc()
print("\n")

# run on timepoint data
for (g in groups) {
  tic(paste("run CellChat", g))
  
  subset_seurat <- subset(seurat_object, subset = timepoint == g)
  cellchat <- prepCellChat(
    subset(seurat_object, subset = timepoint == g),
    "celltype",
    db=neuronChatDB, # input for custom DB
    spatial=TRUE,
    spatial_key="spatial",
    spatial_factors=data.frame(ratio = 1, tol = 5))
  cellchat <- subsetData(cellchat)
  
  cellchat <- runCellChat(cellchat, spatial=TRUE, pathway=has.pathway)
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