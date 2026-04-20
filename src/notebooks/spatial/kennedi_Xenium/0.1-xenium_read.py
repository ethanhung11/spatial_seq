# pixi run -e main python -u xenium_downsample.py &> ../../../../outs/spatial/xenium_downsample.out

from pathlib import Path
from tqdm import tqdm
import spatialdata as sd
import spatialdata_plot as sdp
from spatialdata_io import xenium

DATADIR = Path("../../../kpyper") / "Xenium_Flu"
SAVEDIR = Path("data") / "processed" / "spatial" / "Xenium" / "kennedi_flu"

# Find all folders with "Molofski_Pyper" in their name
data_folders = [folder for folder in DATADIR.rglob("*Molofski_Pyper*") if folder.is_dir()]

for folder in data_folders:
    for out_folder in tqdm([subfolder for subfolder in folder.iterdir() if subfolder.is_dir()], desc=f"writing Xenium output as `h5ad` for folder {folder.name}"):
        print(out_folder)
        sdata = xenium(out_folder)
        sdata.write(SAVEDIR/ out_folder.name)