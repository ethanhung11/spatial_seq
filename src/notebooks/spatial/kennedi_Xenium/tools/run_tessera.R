# Rscript src/notebooks/spatial/kennedi_Xenium/tools/run_tessera.R &> ./outs/spatial/tools/tessera.out

library(here)
Sys.setenv(RETICULATE_PYTHON = here(".pixi/envs/main/bin/python"))
library(reticulate)
library(tictoc)
library(dplyr)
library(tessera)
library(Seurat)
library(scCustomize)
library(sf)

options(future.globals.maxSize= 80*1024^3) # memory allowance
future::plan(future::multicore) # parallelize

DATADIR <- here("data", "processed", "spatial", "Xenium", "kennedi_flu")
SAVEDIR <- here(DATADIR, "tools", "tessera")
dir.create(SAVEDIR, showWarnings = FALSE, recursive = TRUE)

tic("read objects")
seurat_object <- readRDS(here(DATADIR, "integrated.rds"))
seurat_object$idx <- Cells(seurat_object)
coords <- Embeddings(object = seurat_object, reduction = "SPATIAL")
toc()

tic("run Tessera")
res <- GetTiles(
  coords[,1],
  coords[,2],
  counts = seurat_object[["RNA"]]$counts,
  embeddings = Embeddings(object = seurat_object, reduction = "integrated"),
  meta_data = seurat_object@meta.data,
  meta_vars_include = c("idx", "Slide", "Sample", "Genotype", "Timepoint"),
  group.by = "Sample",
  prune_thresh_quantile = 0.99, prune_min_cells = 1,
  max_npts = 50, min_npts = 5, alpha = 1,
  verbose = TRUE,
)
toc()

tic("save objects")
# process usable data
seurat_tiles = Seurat::CreateSeuratObject(
  counts = res$aggs$counts, 
  meta.data = tibble::column_to_rownames(data.frame(dplyr::select(res$aggs$meta_data, -shape)), 'id')
)
shapes = res$aggs$meta_data$shape
tile_mapping <- res$dmt$pts[,c("idx", "agg_id")] %>% as.data.frame()

# save data
st_write(shapes, here(SAVEDIR,"tiles.geojson"))
write.csv(tile_mapping, here(SAVEDIR, "point-tile-mapping.csv"))
as.anndata(x = seurat_tiles, main_layer = "counts", other_layers = NULL,
           file_path = SAVEDIR, file_name = "tessera_tiles.h5ad")
toc()