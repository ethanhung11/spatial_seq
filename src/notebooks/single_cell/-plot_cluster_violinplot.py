# pixi run -e main python -u src/notebooks/single_cell/-plot_cluster_violinplot.py &> ./outs/single_cell/cluster_violin_plots.out

# base
import sys
import os
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
warnings.simplefilter("ignore", FutureWarning)

sys.path.insert(0, 'src')
from single_cell.plot import plot_cluster_violinplot
from utils import stopwatch

import logging

logging.basicConfig(level=logging.ERROR)

CORES = 10
DATADIR = Path("data")
PLOTDIR = Path("plots") / "adipose atlas"
REFDIR = Path("references")
MAIN_DIR = DATADIR / "processed" / "single_cell" / "combined"
SUBSETS_DIR = MAIN_DIR / "subsets"

METADATA = ["Diet", "Age", "Depot", "Sex"]
DOUBLETMETHODS = ["scDblFinder", "DoubletFinder", "doubletdetection", "scrublet"]
INT_KEY = "INT_harmony-Identifier"

if __name__ == "__main__":
    START = stopwatch(mode=-1)

    FILE = SUBSETS_DIR / f"macro.h5ad"
    CLUSTER_KEY = "leiden_macro"
    SELECT_RES = [0.4]
    COMPARISON = "Diet"

    targeted_markers = [
        "Lyz2", "Fn1", "Maf", "H2-Ab1", "Retnla", "Ly6e", "Mrc1", "Cd74", "Mgl2", 
        "Gpnmb", "Mmp12", "Ctsd", "Lyve1", "Apoe", "Folr2", "Top2a", "Ccl8", 
        "Itgam", "Itgax", "Cd36", "Soat1", "C1qa", "Pparg", "Ccr1", "Ccr2", 
        "Ifnar1", "Ifnar2", "Ifngr2", "Tgfbr1", "Tgfbr2", "Smad3", "Csf1r", 
        "Alcam", "Cd44", "S100a4", "Gas6", "Tm4sf19", "Lgals3", "Spp1",
    ]

    ################################
    ################################
    ################################

    # load
    adata = sc.read_h5ad(FILE)

    # check markers
    unique_markers = pd.Series(
        list(set([i for j in targeted_markers.values() for i in j]))
    ) if isinstance(targeted_markers, dict) else pd.Series(list(set(targeted_markers)))
    assert np.all(unique_markers.isin(adata.var_names))

    # fix group names
    for col in adata.obs.columns:
        if CLUSTER_KEY in col:
            adata.obs[col] = adata.obs[col].astype(int).astype("category")

    # plots cluster violins
    for res in SELECT_RES:
        res_key = f"{CLUSTER_KEY}_{res:.1f}"
        task = f"plotting {res_key} for {COMPARISON}"
        stopwatch(task, START, mode=0)

        f = plot_cluster_violinplot(
            adata, COMPARISON, res_key, targeted_markers,
        )

        # save plots
        dir_plot = PLOTDIR / "cluster_violins"
        os.makedirs(dir_plot, exist_ok=True)
        f.savefig(dir_plot / f"{res_key}-{COMPARISON}.jpg", bbox_inches='tight')

        stopwatch(task, START, mode=1)

    # os.system("zip -r ./plots/cluster_violins.zip ./plots/cluster_violins")
