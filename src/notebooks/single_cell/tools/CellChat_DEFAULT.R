# Rscript src/notebooks/single_cell/tools/CellChat_DEFAULT.R &> ./outs/single_cell/tools/cellchat.out
options(max.print = 300)
options(stringsAsFactors = FALSE)
library(here)
library(dplyr)
library(tictoc)
source(here("src", "single_cell", "cellchat.R"))

library(here)
Sys.setenv(RETICULATE_PYTHON = here("~/miniconda3/envs/bulk_rnaseq/bin/python"))

tic("start")
DATA_DIR <- here("data", "processed", "single_cell", "combined")
CELLCHAT_DIR <- here(DATA_DIR, "tools", "cellchat")
obj_name <- "eWAT_Male-fb_macs.rds"
dir.create(CELLCHAT_DIR, showWarnings = FALSE)
save_name <- "cellchat_DEFAULT"
has.pathway <- TRUE
# TODO: Make sure to `groups` (see next block) downstream!
# TODO: Make sure to ctrl+f `celltype_hires` downstream!

tic("read seurat object")
seurat_object <- readRDS(here(DATA_DIR, obj_name))
seurat_object <- NormalizeData(seurat_object)
Idents(seurat_object) <- seurat_object$celltype_hires
groups <- levels(seurat_object$Diet)
print(groups)
toc()

# tic("run CellChat (all)")
# cellchat <- prepCellChat(
#   seurat_object,
#   "celltype_hires",
#   "Identifier")
# cellchat <- subsetData(cellchat)
# cellchat <- runCellChat(cellchat, pathway=has.pathway)
# saveRDS(cellchat, file = here(CELLCHAT_DIR, paste0(save_name,"ALL.rds")))
# toc()

# run on group data
for (g in groups) {
  tic(paste("run CellChat, group", g))

  subset_seurat <- subset(seurat_object, subset = Diet == g)
  cellchat <- prepCellChat(
    subset(seurat_object, subset = Diet == g),
    "celltype_hires",
    "Identifier")
  cellchat <- subsetData(cellchat)

  cellchat <- runCellChat(cellchat, pathway=has.pathway)
  saveRDS(cellchat, file = here(CELLCHAT_DIR, paste0(save_name, "_", g, ".rds")))

  toc()
}

tic("run Merged CellChat")
cellchat.list <- lapply(
  here(CELLCHAT_DIR,paste0(save_name, "_", groups, ".rds")),
  function(x) readRDS(x)
)
names(cellchat.list) <- groups

# # CUSTOM - Merge all neurons
# grouped_celltypes <- levels(cellchat.list[[1]]@idents)
# grouped_celltypes[grepl("neu", grouped_celltypes, ignore.case = TRUE)] <- "Neurons"
# grouped_celltypes <- factor(grouped_celltypes)
# cellchat.list <- lapply(cellchat.list, function(x) {mergeInteractions(x, grouped_celltypes)})
# for (g in groups) saveRDS(cellchat.list[[g]], file = here(CELLCHAT_DIR, paste0(save_name, "_", g, ".rds")))

# Merge
cellchat.MERGED <- mergeCellChat.custom(cellchat.list, pathway=has.pathway)
saveRDS(cellchat.MERGED, file = here(CELLCHAT_DIR, paste0(save_name, "_MERGED.rds")))
toc()

toc()