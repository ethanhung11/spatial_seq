# pixi run -e tacco python -u ./src/notebooks/spatial/kennedi_Xenium/tools/run-tacco.py > ./outs/spatial/tools/tacco.out 2>&1
# https://simonwm.github.io/tacco/index.html

import sys
from pathlib import Path
import session_info
sys.path.insert(0, "src")

import tacco as tc
import scanpy as sc
from stopwatch import stopwatch

DATADIR = Path("data") / "processed" / "spatial" / "Xenium" / "kennedi_flu"
REFDIR = Path("references")
SAVEDIR = DATADIR / "tools" / "tacco" 

session_info.show()

if __name__ == "__main__":
    START = stopwatch(mode=-1)
    ref_celltype_key = "cell type"

    TASK = "READ DATA"
    TASKSTART = stopwatch(TASK, START)
    adata_sp = sc.read_h5ad(SAVEDIR / "integrated.h5ad")
    adata_sc = sc.read_h5ad(REFDIR / "FluSobj.h5ad")
    stopwatch(TASK, TASKSTART, mode=1)

    TASK = "RUN TACCO "
    for method in ["OT"]:
        TASKSTART = stopwatch(TASK + method, START)
        adata_sp = tc.tl.annotate(
            adata_sp, adata_sc,
            annotation_key=ref_celltype_key,
            method=method,
            # counts_location=["layers", "counts"],
            result_key=f"LT_{method}")
        stopwatch(TASK, TASKSTART, mode=1)

    TASK = "SAVE RESULTS"
    TASKSTART = stopwatch(TASK, START)
    SAVEDIR.parent.mkdir(parents=True, exist_ok=True)
    adata_sp.write(SAVEDIR / "tacco.h5ad")
    stopwatch(TASK, TASKSTART, mode=1)

    stopwatch("COMPLETED SCRIPT", START, mode=1)
