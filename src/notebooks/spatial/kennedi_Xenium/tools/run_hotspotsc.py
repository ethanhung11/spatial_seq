# pixi run -e main python -u src/notebooks/spatial/kennedi_Xenium/tools/run_hotspotsc.py &> outs/spatial/tools/hotspot.out

import sys
from pathlib import Path
import pickle
import scanpy as sc
import hotspot

sys.path.insert(0, "src")
from stopwatch import stopwatch

DATADIR = Path("data") / "processed" / "spatial" / "Xenium" / "kennedi_flu"
FILENAME = "integrated.h5ad"
SCRIPT_START = stopwatch(mode=-1)
SPATIAL_KEY = 'SPATIAL'

if __name__ == "__main__":
    ###
    TASK = "initialize adata & hotspot object"
    START = stopwatch(TASK, SCRIPT_START)
    adata = sc.read_h5ad(DATADIR / "integrated.h5ad")
    adata = adata[adata.obs["Sample"] == adata.obs["Sample"][0]]
    sc.pp.filter_genes(adata, min_cells=3)

    hs = hotspot.Hotspot(
        adata,
        layer_key="counts",
        model='danb',
        latent_obsm_key=SPATIAL_KEY,
        umi_counts_obs_key="total_counts"
    )
    hs.create_knn_graph(
        weighted_graph=False, n_neighbors=20,
    )
    stopwatch(TASK, START, mode=1)

    ###

    TASK = "COMPUTE AUTOCORRELATIONS"
    START = stopwatch(TASK, SCRIPT_START)
    hs_results = hs.compute_autocorrelations(jobs=20)
    for cutoff in [5e-2, 1e-2, 1e-3, 1e-4]:
        print(f"cuttoff {cutoff}: {(hs_results.FDR < cutoff).sum()} genes")
    hs_genes = hs_results.index[hs_results.FDR < 0.01]
    stopwatch(TASK, START, mode=1)

    ###

    TASK = "COMPUTE LOCAL CORRELATION FROM TOP AUTOCORRELATED GENES"
    START = stopwatch(TASK, SCRIPT_START)
    lcz = hs.compute_local_correlations(hs_genes, jobs=20)
    modules = hs.create_modules(
        min_gene_threshold=50, core_only=True, fdr_threshold=0.01
    )
    stopwatch(TASK, START, mode=1)

    ###

    TASK = "SAVE RUN"
    START = stopwatch(TASK, SCRIPT_START)

    SAVEDIR = DATADIR / "tools" / "hotspotSC"
    SAVEDIR.mkdir(parents=True, exist_ok=True)
    with open(SAVEDIR / "hotspot.pkl", 'wb') as f:
        pickle.dump(hs, f)
    stopwatch(TASK, START, mode=1)

    ###

    stopwatch("COMPLETED SCRIPT", SCRIPT_START, mode=1)