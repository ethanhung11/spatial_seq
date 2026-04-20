# pixi run -e main python -u src/notebooks/spatial/kennedi_Xenium/xenium_downsample.py &> outs/spatial/xenium_downsample.out

import os
import sys
from pathlib import Path
import session_info
import warnings
warnings.simplefilter("ignore", FutureWarning)

import scanpy as sc
import spatialdata as sd
sc.settings.n_jobs = 20

SRCDIR = Path('../../..')
HOMEDIR = SRCDIR / '..'
DATADIR = HOMEDIR / "data" / "processed" / "spatial" / "Xenium" / "kennedi_flu"
if not str(SRCDIR) in sys.path:
    sys.path.insert(0, str(SRCDIR))

from utils import Downsample
from stopwatch import stopwatch
session_info.show()


if __name__ == "__main__":
    START = stopwatch(mode=-1)
    N = 50000
    str_N = f"{int(N/1e6)}M" if N>=1e6 else f"{int(N/1e3)}K"

    sdata = sd.read_zarr(DATADIR / "integrated.zarr")
    TASKSTART = stopwatch("DOWNSAMPLING", START, 0)
    adata_small = Downsample(sdata['table'], key="integrated", N=N, method="scsampler")
    stopwatch("DOWNSAMPLING", TASKSTART, 1)

    adata_small.write(DATADIR / f"downsampled_{str_N}.h5ad")
    stopwatch("COMPLETED", START, 1)