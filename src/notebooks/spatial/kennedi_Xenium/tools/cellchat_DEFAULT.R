# Use this command to run it. Relies on having a [project].Rproj file and repo root, for the here() reference.
# Rscript src/notebooks/spatial/kennedi_Xenium/tools/cellchat_DEFAULT.R &> ./outs/spatial/tools/cellchat.out

options(max.print = 300)
options(stringsAsFactors = FALSE)
library(here)
library(dplyr)
library(tictoc)

# TODO: Adjust inputs
tic("start")
source(here("src", "single_cell", "cellchat.R"))
DATADIR <- here("data", "processed", "spatial", "Xenium", "kennedi_flu")
CELLCHAT_DIR <- here(DATADIR, "tools", "cellchat")
FILENAME <- "integrated.rds"
savename <- "cellchat_DEFAULT"

CELLTYPE_KEY <- "celltype0518"
SAMPLE_KEY <- "Sample"
GROUP_KEY <- "Groups"
has.pathway <- TRUE

USE_SPATIAL <- FALSE
SPATIAL_KEY <- "SPATIAL"
SPOT_SIZE <- 1
SPATIAL_FACTOR <- data.frame(ratio = 1, tol = SPOT_SIZE/2)
# Adjust `ratio` only -- this should be spot size / pixels per spot --> um / voxel
# e.g. Visium has 65um spots / 200 voxel units per spot --> 0.325
# e.g. Xenium has 1um units / 1 voxel --> 1
# Reference: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat_analysis_of_multiple_spatial_transcriptomics_datasets.html#create-a-cellchat-object:~:text=CellChat)%0Alibrary(patchwork)-,Part%20I%3A%20Data%20input%20%26%20processing%20and%20initialization%20of%20CellChat%20object,-Load%20data


# input processing
dir.create(CELLCHAT_DIR, showWarnings = FALSE)
if (USE_SPATIAL == FALSE) {
  SPATIAL_KEY <- NULL
  SPATIAL_FACTOR <- NULL
}

tic("read seurat object")
seurat_object <- readRDS(here(DATADIR, FILENAME))
seurat_object <- NormalizeData(seurat_object)
Idents(seurat_object) <- seurat_object@meta.data[CELLTYPE_KEY][[CELLTYPE_KEY]]
groups <- levels(as.factor(seurat_object@meta.data[GROUP_KEY][[GROUP_KEY]]))

groups <- grep("Tgfbr", groups, value = TRUE)
print(groups)
toc()

# Assume for Xenium 1:1 voxel to um. Double check with technology.
# Run across all cells, no separation of condition
tic("run CellChat (all)")
cellchat <- prepCellChat(
  seurat_object,
  CELLTYPE_KEY,
  SAMPLE_KEY,
  spatial=USE_SPATIAL,
  spatial_key=SPATIAL_KEY,
  spatial_factors=SPATIAL_FACTOR)
cellchat <- subsetData(cellchat)
cellchat <- runCellChat(cellchat, spatial=USE_SPATIAL, pathway=has.pathway)
saveRDS(cellchat, file = here(CELLCHAT_DIR, paste0(savename, "_", "ALL.rds")))
toc()

# Run across cells for each condition separately
for (g in groups) {
  tic(paste("run CellChat", g))

  subset_seurat <- seurat_object[, seurat_object[[GROUP_KEY]] == g]
  cellchat <- prepCellChat(
    subset_seurat,
    CELLTYPE_KEY,
    SAMPLE_KEY,
    spatial=USE_SPATIAL,
    spatial_key=SPATIAL_KEY,
    spatial_factors=SPATIAL_FACTOR
  )
  cellchat <- subsetData(cellchat)
  cellchat <- runCellChat(cellchat, pathway=has.pathway, spatial=USE_SPATIAL)
  saveRDS(cellchat, file = here(CELLCHAT_DIR, paste0(savename, "_", g, ".rds")))

  toc()
}

# Merge conditions where applicable for valid comparisons
# per Kennedi: rest-d3-d7-d14 for ctrl only & ctrl-tgfbf per timepoint
tic("run Merged CellChat")
cellchat.list <- lapply(
  here(CELLCHAT_DIR,paste0(savename, "_", groups, ".rds")),
  function(x) readRDS(x)
)
names(cellchat.list) <- groups

# # CUSTOM - Merge all celltypes of the same kind, if useful. My commented example below was for Pailin's various neuron celltypes.
# grouped_celltypes <- levels(cellchat.list[[1]]@idents)
# grouped_celltypes[grepl("neu", grouped_celltypes, ignore.case = TRUE)] <- "Neurons"
# grouped_celltypes <- factor(grouped_celltypes)
# cellchat.list <- lapply(cellchat.list, function(x) {mergeInteractions(x, grouped_celltypes)})
# for (g in groups) saveRDS(cellchat.list[[g]], file = here(CELLCHAT_DIR, paste0(savename, "_", g, ".rds")))

# Merge
cellchat.MERGED <- mergeCellChat.custom(cellchat.list, pathway=has.pathway)
saveRDS(cellchat.MERGED, file = here(CELLCHAT_DIR, paste0(savename, "-ALL_MERGED.rds")))

cellchat.tp <- cellchat.list[c('Ctrl-Rest', 'Ctrl-D3', 'Ctrl-D7', 'Ctrl-D14')]
cellchat.MERGED <- mergeCellChat.custom(cellchat.tp, pathway=has.pathway)
saveRDS(cellchat.MERGED, file = here(CELLCHAT_DIR, paste0(savename, "-TIMEPOINTS.rds")))

for (tp in c('Rest', 'D3', 'D7', 'D14')){
  cellchat.geno <- cellchat.list[paste0(c('Ctrl-', 'Tgfbr2F-'), tp)]
  cellchat.MERGED <- mergeCellChat.custom(cellchat.geno, pathway=has.pathway)
  saveRDS(cellchat.MERGED, file = here(CELLCHAT_DIR, paste0(savename, "-", tp, "_GENOTYPES.rds")))
}
toc()

toc()