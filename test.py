import os
from pathlib import Path
import pandas as pd
import tangram as tg
import anndata as ad
import scanpy as sc

DATADIR = Path("data")

if __name__ == "__main__":
    adata_sp = ad.read_zarr(
        DATADIR / "processed" / "spatial" / "analysis" / "analyzed_zarr"
    )
    print("spatial data read\n")
