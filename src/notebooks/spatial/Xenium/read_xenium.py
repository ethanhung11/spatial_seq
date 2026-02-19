# uv run python src/notebooks/spatial/Xenium/read_xenium.py

from pathlib import Path
import spatialdata as sd
import spatialdata_plot as sdp
from spatialdata_io import xenium

DATADIR = Path("data") / "processed" / "spatial"

sdata = xenium(Path("../xenium/0038439/outs"))
sdata.write(DATADIR / "Xenium" / "0038439.zarr")
