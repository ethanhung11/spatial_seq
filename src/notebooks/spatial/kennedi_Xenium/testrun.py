# pixi run -e main python -u src/notebooks/spatial/kennedi_Xenium/testrun.py &> outs/spatial/testrun.out

import gc
import os
import sys
import logging
import warnings
from pathlib import Path
import matplotlib as mpl
import pickle
import session_info

HOMEDIR = Path('.')
SRCDIR = HOMEDIR / "src"
PLOTDIR = HOMEDIR / "plots" / "kennedi_xenium"
DATADIR = HOMEDIR / "data" / "processed" / "spatial" / "Xenium" / "kennedi_flu"
if not str(SRCDIR) in sys.path:
    sys.path.insert(0, str(SRCDIR))

logging.basicConfig(level="WARNING")
warnings.simplefilter("ignore", FutureWarning)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)
warnings.simplefilter("ignore", DeprecationWarning)

from single_cell.preprocess import *
from single_cell.plot import *
from single_cell.analysis import *
from spatial_seq.plot import *
from utils import *

import squidpy as sq
import cellcharter as cc
# import spatialdata as sd
# import spatialdata_plot as sdp

CORES = 20
mpl.rcdefaults()
plt.rcParams["figure.figsize"] = (8, 8)
gc.collect()

session_info.show()

if __name__ == "__main__":
    # load
    adata = sc.read_h5ad(DATADIR / "integrated.h5ad")
    THRESHOLD = adata.uns["hotspot"]["MinModuleSize"]
    obs_to_obsm(adata, "HotspotModule", "obs")

    results = adata.uns["hotspot"]["gene_df"]
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'grays', ['#DDDDDD', '#000000'])

    # Plot Hotspot genes spatially
    for mod in tqdm(p.arange(15)+1):
        genes = results[results["Module"] == mod].sort_values("Z").index

        f,axs = plt.subplots(2,7, figsize=(20*1.5,5*1.5), layout="constrained", dpi=400)
        for i,samp in enumerate(["Ctrl_D7_2", "Tgfbr2F_D7_2"]):
            sc.pl.embedding(
                adata[adata.obs["Sample"] == samp], 
                "SPATIAL",
                color=f"HotspotModule-{mod}", 
                cmap="Reds", 
                frameon=False,
                ax=axs[i][0], show=False,
                title=samp
            )
            for j,gene in enumerate(genes[:5]):
                sc.pl.embedding(
                    adata[adata.obs["Sample"] == samp],
                    "SPATIAL",
                    color=gene, 
                    cmap=cmap, 
                    frameon=False,
                    vmin='p0',
                    vmax='p95',
                    ax=axs[i][j+1], show=False,
                )

        N = 7
        chunks = [genes[i:i + N] for i in range(0, len(genes), N)]
        gene_str = '\n'.join([', '.join(c) for c in chunks])
        gene_str = f"Module {mod} genes:\n" + gene_str
        axs[0][-1].text(0, 1, gene_str, va="top", size=15, linespacing=2)
        axs[0][-1].axis('off')
        axs[1][-1].remove()

        f.suptitle(f"HotspotModule{mod}", size=25)
        f.savefig(PLOTDIR / "hotspot" / f"Spatial - HotspotSC-Thres{THRESHOLD} - HotspotModule-{mod} Genes.jpg")

    # # Plot Hotspot modules spatially
    # for key in adata.obs.columns[adata.obs.columns.str.contains("HotspotModule")]:
    #     f,ax = plt.subplots(1, 1, figsize=(50,20), layout="constrained", dpi=400)
    #     sc.pl.embedding(
    #         adata, basis='SPATIAL', color=key, 
    #         size=5, cmap="Reds",
    #         ax=ax, show=False, frameon=False)
    #     spatial_label(adata, "SPATIAL", "Sample", ax, fs=20)

    #     f.savefig(PLOTDIR / "hotspot" / f"Spatial - HotspotSC-Thres{THRESHOLD} - {key}.jpg")

    # # Plot Hotspot expression by label
    # f = plot_violinplot(
    #     adata, "LabelTransfer_OT",
    #     adata.obs.columns[adata.obs.columns.str.contains("HotspotModule")],
    #     y_fontsize=8)
    # f.savefig(PLOTDIR / "hotspot" / f"HotspotSC_Violin-Thres{THRESHOLD} - Celltype")
    # for key in ["cellcharter_k6", "cellcharter_k12"]:
    #     f = plot_violinplot(
    #         adata, key,
    #         adata.obs.columns[adata.obs.columns.str.contains("HotspotModule")],
    #         y_fontsize=8)
    #     f.savefig(PLOTDIR / "hotspot" / f"HotspotSC_Violin-Thres{THRESHOLD} - Niche {key}")