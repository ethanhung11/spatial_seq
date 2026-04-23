# pixi run -e main python -u ./src/notebooks/spatial/kennedi_Xenium/tools/run-rctdpy.py > ./outs/spatial/tools/rctd.out 2>&1
# https://github.com/p-gueguen/rctd-py
# https://p-gueguen.github.io/rctd-py/tutorial.html

import pickle
import sys
from pathlib import Path
sys.path.insert(0, "src")

from rctd import Reference, run_rctd, RCTDConfig
import scanpy as sc
from stopwatch import stopwatch

DATADIR = Path("data") / "processed" / "spatial" / "Xenium" / "kennedi_flu"
REFDIR = Path("references")
SAVEDIR = DATADIR / "tools" / "rtcd"

if __name__ == "__main__":
    START = stopwatch(mode=-1)

    DEG_key = "DEG_global"

    TASK = "READ DATA"
    TASKSTART = stopwatch(TASK, START)
    adata_sp = sc.read_h5ad(DATADIR / "integrated.h5ad")
    adata_sc = Reference(sc.read_h5ad(REFDIR / "FluSobj.h5ad"), cell_type_col="cell type")
    START = stopwatch(TASK, TASKSTART, mode=1)

    TASK = "RUN RTCD"
    TASKSTART = stopwatch(TASK, START)
    result = run_rctd(adata_sp, adata_sc, mode="doublet", config=RCTDConfig(UMI_min=1, UMI_min_sigma=1))
    START = stopwatch(TASK, TASKSTART, mode=1)

    TASK = "SAVE RESULTS"
    TASKSTART = stopwatch(TASK, START)
    SAVEDIR.parent.mkdir(parents=True, exist_ok=True)
    with open(SAVEDIR, 'wb') as f:
        pickle.dump(result, f)
    START = stopwatch(TASK, TASKSTART, mode=1)

    stopwatch("COMPLETED SCRIPT", START, mode=1)
