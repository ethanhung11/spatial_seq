# pixi shell -e scenic
# python src/notebooks/single_cell/2.0-SCENIC.py &> ./outs/single_cell/tools/scenic.out

import os
import subprocess
from pathlib import Path
import scanpy as sc
import numpy as np
import pandas as pd
import logging

logging.basicConfig(level=logging.WARNING)
import warnings

warnings.filterwarnings("ignore")

from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase

# from arboreto.algo import grnboost2
# from pyscenic.utils import modules_from_adjacencies
# from pyscenic.prune import prune2df, df2regulons
# from pyscenic.aucell import aucell
import loompy as lp

from time import time, ctime
from datetime import timedelta


def stopwatch(taskname: str = "", start=time(), mode=0):
    if mode == 0:
        print(
            (
                f"\n\nbegin {taskname} at {ctime(time())}, "
                f"(total: {timedelta(seconds=time()-start)})"
            )
        )
    elif mode == 1:
        print(f"finished in {timedelta(seconds=time()-start)}\n")
    elif mode == 2:
        pass
    else:
        raise ValueError("must be mode 0, 1, or 2")
    return time()


if __name__ == "__main__":
    SCRIPT_START = stopwatch(mode=2)

    annotation = "eWAT_Male_SLIM.h5ad"
    CORES = 20
    DATADIR = Path("data") / "processed" / "single_cell"
    OUTSDIR = Path("outs") / "single_cell" / "tools"
    REFDIR = Path("references")
    MAIN_DIR = DATADIR / "combined"
    TOOL_DIR = MAIN_DIR / "tools" / "scenic"
    TOOL_DIR.mkdir(parents=True, exist_ok=True)

    GRN_RANKINGS = [
        RankingDatabase(fname=fname, name=os.path.basename(fname))
        for fname in list(map(str, (REFDIR / "scenic").glob("*feather")))
    ]
    MOTIF_ANNOTS = str(REFDIR / "scenic" / "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")
    TF_LIST = str(REFDIR / "scenic" / "mm_mgi_tfs.txt")
    AUCELL_OUTPUT = str(TOOL_DIR / f"{annotation}-aucell_output.csv")
    sc.settings.njobs = CORES

    # get adata
    start = stopwatch("read anndata", SCRIPT_START)
    adata = sc.read_h5ad(MAIN_DIR / annotation)
    stopwatch(start, mode=1)

    ## COMMAND LINE TOOL VERSION
    path_loom_input = str(TOOL_DIR / "input.loom")
    path_adjacency = str(TOOL_DIR / "scenic_adj.csv")
    path_regulon_pred = str(TOOL_DIR / "scenic_reg.csv")
    path_aucell_out = str(TOOL_DIR / "aucell_mtx.csv")

    # Create Loompy input file
    start = stopwatch("create loompy", SCRIPT_START)
    lp.create(
        path_loom_input,
        adata.X.transpose(),
        {
            "Gene": np.array(adata.var.index),
        },
        {"CellID": np.array(adata.obs.index)},
    )
    stopwatch(start, mode=1)

    # Run Grnboost2 to generate predicted network adjacency
    start = stopwatch("run Grnboost2", SCRIPT_START)
    if not Path(path_adjacency).exists():
        subprocess.run(
            [
                "pyscenic",
                "grn",
                path_loom_input,
                TF_LIST,
                "--num_workers",
                str(CORES),
                "--output",
                path_adjacency,
            ],
            check=True,
        )

    results_adjacencies = pd.read_csv(path_adjacency, index_col=False, sep=",")
    print(f"Number of associations: {results_adjacencies.shape[0]}")
    print(results_adjacencies.head())
    stopwatch(start, mode=1)

    # Run CTX regulon predicton
    start = stopwatch("run CTX regulon predicton", SCRIPT_START)
    if not Path(path_regulon_pred).exists():
        subprocess.run(
            [
                "time" "pyscenic",
                "ctx",
                path_adjacency,
                str(GRN_RANKINGS),
                "--annotations_fname",
                MOTIF_ANNOTS,
                "--expression_mtx_fname",
                path_loom_input,
                "--output",
                path_regulon_pred,
                "--num_workers",
                str(CORES),
                "--mask_dropouts",
            ],
            check=True,
        )

    results_regulons = pd.read_csv(path_regulon_pred, index_col=False, sep=",")
    print(f"Number of predicted regulons: {results_regulons.shape[0]}")
    print(results_regulons.head())
    stopwatch(start, mode=1)

    # Run AUCell
    start = stopwatch("run AUCell", SCRIPT_START)
    outfile = OUTSDIR / "scenic-3_aucell.out"
    if not Path(path_aucell_out).exists():
        subprocess.run(
            [
                "pyscenic",
                "aucell",
                path_loom_input,
                path_regulon_pred,
                "--output",
                path_aucell_out,
                "--num_workers",
                str(CORES),
            ],
            check=True,
        )
    auc_mtx = pd.read_csv(path_aucell_out)
    stopwatch(start, mode=1)

    print(f"SCENIC analysis complete. AUCell matrix shape: {auc_mtx.shape}")

    ### HOW TO ACCESS!
    # AUCELL_OUTPUT = str(MAIN_DIR / "tools" / "scenic" / f"{annotation}-aucell_output.csv")
    # auc_mtx = pd.read_csv(AUCELL_OUTPUT, index_col=0)
    # adata.obsm["SCENIC"] = auc_mtx
    # PCA(adata, gene_mask=None, key="SCENIC_PCA", obsm="SCENIC")
    # Visualize(adata,"SCENIC",obsm="SCENIC_PCA")
    # dc.pp.get_obsm(adata=adata, key="SCENIC")

    # # Step 1: GRN inference using GRNBoost2
    # start = stopwatch("XGBoost", begin)
    # ex_matrix = pd.DataFrame(
    #     adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X,
    #     index=adata.obs_names,
    #     columns=adata.var_names
    # ).astype(np.float32)
    # adjacencies = grnboost2(ex_matrix, tf_names=TF_LIST, verbose=True)
    # adjacencies.to_csv(TOOL_DIR / f"{annotation}-adjacencies.csv", index=False)
    # stopwatch("", start, 1)

    # # Step 2: Regulon prediction (cisTarget)
    # start = stopwatch("cisTarget", begin)
    # df = prune2df(GRN_RANKINGS, adjacencies, str(MOTIF_ANNOTS), num_workers=N_WORKERS)
    # df.to_csv(TOOL_DIR / f"{annotation}-regulons.csv", index=False)
    # regulons = df2regulons(df)
    # stopwatch("", start, 1)

    # # Step 3: AUCell scoring
    # start = stopwatch("AUCell", begin)
    # auc_mtx = aucell(ex_matrix, regulons, num_workers=N_WORKERS)
    # auc_mtx.to_csv(AUCELL_OUTPUT)
    # stopwatch("", start, 1)

    # # Add AUCell scores to adata
    # adata.obsm['X_aucell'] = auc_mtx.loc[adata.obs_names, :].values
    # adata.uns['regulons'] = {reg.name: list(reg.gene2weight.keys()) for reg in regulons}

    # # Save results
    # adata.write(TOOL_DIR / f"{annotation}-scenic.h5ad")
