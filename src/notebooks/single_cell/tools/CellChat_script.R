# Rscript src/notebooks/single_cell/2.2-CellChat_script.R &> ./outs/spatial/tools/cellchat.out
options(max.print = 300)
options(stringsAsFactors = FALSE)
library(here)
library(dplyr)
library(tictoc)
source(here("src", "single_cell", "cellchat.R"))

tic("start")
DATADIR <- here("data", "processed", "single_cell", "combined")
CELLCHAT_DIR <- here(DATADIR, "tools", "cellchat")
dir.create(CELLCHAT_DIR, showWarnings = FALSE)

tic("read seurat object")
seurat_object <- readRDS(here(DATADIR, "eWAT_Male.rds"))
seurat_object <- NormalizeData(seurat_object)
Idents(seurat_object) <- seurat_object$celltype
groups <- levels(seurat_object$Diet)
toc()

tic("run cellchat")
cellchat <- prepCellChat(
  seurat_object,
  "celltype",
  spatial=TRUE,
  spatial_key="spatial",
  spatial_factors=data.frame(ratio = 1, tol = 5))
cellchat <- subsetData(cellchat)
cellchat <- runCellChat(cellchat, spatial=TRUE, pathway=TRUE)
saveRDS(cellchat, file = here(CELLCHAT_DIR, "cellchat_ALL.rds"))
toc()

# run on timepoint data
for (g in groups) {
  tic(paste("run cellchat", g))
  
  subset_seurat <- subset(seurat_object, subset = timepoint == g)
  cellchat <- prepCellChat(
    subset(seurat_object, subset = timepoint == g),
    "celltype",
    spatial=TRUE,
    spatial_key="spatial",
    spatial_factors=data.frame(ratio = 1, tol = 5))
  cellchat <- subsetData(cellchat)
  
  cellchat <- runCellChat(cellchat, spatial=TRUE, pathway=TRUE)
  saveRDS(cellchat, file = here(CELLCHAT_DIR, paste0("cellchat_", g, ".rds")))
  
  toc()
}

cellchat.TIMEPOINTS <- lapply(
  here(CELLCHAT_DIR, paste0("cellchat_", groups, ".rds")),
  function(x) readRDS(x)
)
cellchat.TIMEPOINTS_MERGED <- mergeCellChat(cellchat.TIMEPOINTS, add.names = groups)
for (slot in c("net", "netP")) {
  cellchat.TIMEPOINTS_MERGED <- computeNetSimilarityPairwise(cellchat.TIMEPOINTS_MERGED, type = "functional", slot.name = slot)
  cellchat.TIMEPOINTS_MERGED <- netEmbedding(cellchat.TIMEPOINTS_MERGED, type = "functional", slot.name = slot)
  cellchat.TIMEPOINTS_MERGED <- netClustering(cellchat.TIMEPOINTS_MERGED, type = "functional", slot.name = slot)
}
saveRDS(cellchat.TIMEPOINTS_MERGED, file = here(CELLCHAT_DIR, paste0("cellchat_merged_timepoints.rds")))

toc()