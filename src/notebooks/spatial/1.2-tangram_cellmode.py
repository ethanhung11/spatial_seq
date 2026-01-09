# taskset -c 21-32  uv run python ./src/notebooks/spatial/1.2-tangram_cellmode.py > outs/tangram-segmented.out 2>&1

from pathlib import Path
from time import time, ctime
from datetime import timedelta

import pandas as pd
import tangram as tg
import anndata as ad
import scanpy as sc

DATADIR = Path("data")

if __name__ == "__main__":
    start = time()
    print(f"begin script at {ctime(time())}")

    task_start = time()
    print(
        f"begin reading in datasets at {ctime(task_start)}, (total: {timedelta(seconds=task_start-start)})"
    )
    adata_sp = ad.read_zarr(
        DATADIR / "processed" / "spatial" / "analysis" / "analyzed_zarr"
    )
    adata_sc = sc.read_h5ad(
        DATADIR
        / "processed"
        / "single cell"
        / "So2025"
        / "5_analysis"
        / "annotated-ccc.h5ad"
    )
    print(f"finished in {timedelta(seconds=time()-task_start)}\n")

    task_start = time()
    print(
        f"begin finding DEGs at {ctime(task_start)}, (total: {timedelta(seconds=task_start-start)})"
    )
    DEGlist = sc.get.rank_genes_groups_df(
        adata_sc,
        group=None,
        key="de_all",
        pval_cutoff=1e-50,
        log2fc_min=0.1,
    ).sort_values("pvals_adj")
    DEGlist = pd.concat([df.head(200) for g, df in DEGlist.groupby("group")])
    markers = DEGlist.names.unique()
    print(f"finished in {timedelta(seconds=time()-task_start)}\n")

    task_start = time()
    print(
        f"begin running tangram at {ctime(task_start)}, (total: {timedelta(seconds=task_start-start)})"
    )
    tg.pp_adatas(
        adata_sc,
        adata_sp,
        genes=markers,
    )
    ad_map = tg.map_cells_to_space(
        adata_sc,
        adata_sp,
        mode="constrained",
        target_count=adata_sp.shape[0],
        density_prior="uniform",
        cluster_label="cell_type",
        device="cpu",
    )
    print(f"finished in {timedelta(seconds=time()-task_start)}\n")

    task_start = time()
    print(
        f"begin saving tangram outputs at {ctime(task_start)}, (total: {timedelta(seconds=task_start-start)})"
    )
    ad_map.write_zarr(DATADIR / "processed" / "spatial" / "analysis" / "tangram")
    print(f"finished in {timedelta(seconds=time()-task_start)}\n")
