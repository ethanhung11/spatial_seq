# pixi run -e main python -u ./src/notebooks/spatial/kennedi_Xenium/tools/run-tangram.py > ./outs/spatial/tools/tangram.out 2>&1
# https://github.com/broadinstitute/Tangram

import sys
from pathlib import Path
import session_info
sys.path.insert(0, "src")

import tangram as tg
import scanpy as sc
from stopwatch import stopwatch

DATADIR = Path("data") / "processed" / "spatial" / "Xenium" / "kennedi_flu"
REFDIR = Path("references")
SAVEDIR = DATADIR / "tools" / "tangram" 

session_info.show()

if __name__ == "__main__":
    START = stopwatch(mode=-1)

    DEG_key = "DEG_global"

    TASK = "READ DATA"
    TASKSTART = stopwatch(TASK, START)
    adata_sp = sc.read_h5ad(DATADIR / "integrated.h5ad")
    adata_sc = sc.read_h5ad(REFDIR / "FluSobj.h5ad")
    START = stopwatch(TASK, TASKSTART, mode=1)

    TASK = "GET DEGs"
    TASKSTART = stopwatch(TASK, START)
    DEGlist = sc.get.rank_genes_groups_df(
        adata_sc,
        group=None,
        key=DEG_key,
        pval_cutoff=1e-50,
        log2fc_min=0.1,
    ).sort_values("pvals_adj")
    # DEGlist = pd.concat([df.head(200) for g, df in DEGlist.groupby("group")])
    markers = DEGlist.names.unique()
    START = stopwatch(TASK, TASKSTART, mode=1)

    TASK = "RUN TANGRAM"
    TASKSTART = stopwatch(TASK, START)
    tg.pp_adatas(
        adata_sc,
        adata_sp,
        genes=markers,
    )
    ad_map = tg.map_cells_to_space(
        adata_sc,
        adata_sp,
        mode="clusters",
        cluster_label='cell type'
        # target_count=adata_sp.shape[0],
        # density_prior="uniform",
        # cluster_label="cell_type",
        # device="cpu",
    )
    START = stopwatch(TASK, TASKSTART, mode=1)

    TASK = "SAVE RESULTS"
    TASKSTART = stopwatch(TASK, START)
    SAVEDIR.parent.mkdir(parents=True, exist_ok=True)
    ad_map.write(SAVEDIR / "tangram.h5ad")
    START = stopwatch(TASK, TASKSTART, mode=1)

    stopwatch("COMPLETED SCRIPT", START, mode=1)
