# Rscript src/notebooks/pailin/tools/cellchat_DEFAULT.R &> ./outs/spatial/tools/cellchat.out
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
savename <- "cellchat_DEFAULT"
has.pathway <- TRUE

tic("read seurat object")
seurat_object <- readRDS(here(DATADIR, "analyzed.rds"))
seurat_object <- NormalizeData(seurat_object)
Idents(seurat_object) <- seurat_object$celltype
groups <- levels(seurat_object$timepoint)
toc()

tic("run CellChat (all)")
cellchat <- prepCellChat(
  seurat_object,
  "celltype",
  "sample",
  spatial=TRUE,
  spatial_key="spatial",
  spatial_factors=data.frame(ratio = 1, tol = 5))
cellchat <- subsetData(cellchat)
cellchat <- runCellChat(cellchat, spatial=TRUE, pathway=has.pathway)
saveRDS(cellchat, file = here(CELLCHAT_DIR, paste0(savename,"ALL.rds")))
toc()

# run on timepoint data
for (g in groups) {
  tic(paste("run CellChat", g))

  subset_seurat <- subset(seurat_object, subset = timepoint == g)
  cellchat <- prepCellChat(
    subset(seurat_object, subset = timepoint == g),
    "celltype",
    spatial=TRUE,
    spatial_key="spatial",
    spatial_factors=data.frame(ratio = 1, tol = 5))
  cellchat <- subsetData(cellchat)

  cellchat <- runCellChat(cellchat, spatial=TRUE, pathway=has.pathway)
  saveRDS(cellchat, file = here(CELLCHAT_DIR, paste0(savename, "_", g, ".rds")))

  toc()
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