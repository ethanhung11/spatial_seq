import hotspot
import anndata as ad
import scanpy as sc
import pickle
from utils import stopwatch
from pathlib import Path

DATADIR = Path("data") / "processed" / "external"
SCRIPT_START = stopwatch(mode=2)

###

task = "initialize adata & hotspot object"
start = stopwatch(task, SCRIPT_START)
adata = ad.read_zarr(DATADIR / "analyzed.zarr")
adata = adata[adata.obs["sample"] == adata.obs["sample"][0]]
sc.pp.filter_genes(adata, min_cells=3)

hs = hotspot.Hotspot(
    adata,
    layer_key="counts",
    model='danb',
    latent_obsm_key="spatial",
    umi_counts_obs_key="total_counts"
)
hs.create_knn_graph(
    weighted_graph=False, n_neighbors=20,
)
stopwatch(task, start, mode=1)

###

task = "compute autocorrelations"
start = stopwatch(task, SCRIPT_START)
hs_results = hs.compute_autocorrelations(jobs=20)
for cutoff in [5e-2, 1e-2, 1e-3, 1e-4]:
    print(f"cuttoff {cutoff}: {(hs_results.FDR < cutoff).sum()} genes")
hs_genes = hs_results.index[hs_results.FDR < 0.01]
stopwatch(task, start, mode=1)

###

task = "compute local correlations from top autocorrelated genes"
start = stopwatch(task, SCRIPT_START)
lcz = hs.compute_local_correlations(hs_genes, jobs=20)
modules = hs.create_modules(
    min_gene_threshold=50, core_only=True, fdr_threshold=0.01
)
stopwatch(task, start, mode=1)

###

task = "save final object"
start = stopwatch(task, SCRIPT_START)

SAVEDIR = DATADIR / "tools" / "hotspotsc"
SAVEDIR.mkdir(parents=True, exist_ok=True)
with open(SAVEDIR / "hotspot.pkl", 'wb') as f:
    pickle.dump(hs, f)
stopwatch(task, start, mode=1)

###

stopwatch("script", SCRIPT_START, mode=1)