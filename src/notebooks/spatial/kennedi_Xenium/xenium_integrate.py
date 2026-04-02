# pixi run -e main python -u xenium_integrate.py &> ../../../../outs/spatial/integrate_xenium.out

import os
import sys
from pathlib import Path
import session_info

SRCDIR = Path('../../..')
HOMEDIR = SRCDIR / '..'
DATADIR = HOMEDIR / "data" / "processed" / "spatial" / "Xenium" / "kennedi_flu"
if not str(SRCDIR) in sys.path:
    sys.path.insert(0, str(SRCDIR))

from single_cell.preprocess import *
import spatialdata as sd
session_info.show()

CLUSTER_KEY = "global_leiden"
DE_KEY = "global_DEG"
resolutions = np.arange(1, 16) / 10

if __name__ == "__main__":
    # read all xenium outs
    sdatas = {
        fname.split("__")[1]:
        sd.read_zarr(DATADIR / fname)
        for fname in tqdm(os.listdir(DATADIR), desc="Reading spatialdata objects")
    }

    # read in annotations
    REFDIR = Path("../../../../../../../kpyper") / "Xenium_Flu"
    data_folders = [folder for folder in REFDIR.rglob("*Molofski_Pyper*") if folder.is_dir()]
    annotations = {}
    for folder in data_folders:
        for out_folder in [subfolder for subfolder in folder.iterdir() if subfolder.is_dir()]:
            slide = out_folder.name.split('__')[1]
            annotations[slide] = {}
            for file in os.listdir(out_folder / "regions"):
                sample = file.split(".csv")[0]
                df = pd.read_csv(out_folder / "regions" / file, skiprows=2)
                annotations[slide][sample] = df.copy()

    # filter on annotations
    for slide in annotations.keys():
        # switch to cell segmentation
        sdatas[slide]['table'].obs['region'] = "cell_boundaries"
        sdatas[slide]['table'].obs['region'] = sdatas[slide]['table'].obs['region'].astype(str).astype("category")
        sdatas[slide]['table'].uns['spatialdata_attrs']['region'] = "cell_boundaries"

        # annotate by Sample by cell_id
        sdatas[slide]['table'].obs["Slide"] = slide
        sdatas[slide]['table'].obs["Sample"] = None
        for sample in annotations[slide]:
            idx = sdatas[slide]['table'].obs['cell_id'].isin(annotations[slide][sample]['Cell ID'])
            sdatas[slide]['table'].obs.loc[idx,"Sample"] = sample

        # find cells with no Sample annotation
        to_remove = sdatas[slide]['table'].obs["Sample"].isna()
        print(f"Slide {slide}: {sum(to_remove)} cells removed")

        # remove said cells
        sdatas[slide]['table'] = sdatas[slide]['table'][~to_remove]

    # combine samples
    sdata = sd.concatenate(sdatas, concatenate_tables=True)

    # standard RNAseq preprocessing
    adata = sdata.tables["table"]
    adata.layers["counts"] = adata.X.copy()
    adata = Filter_QC(adata, 10, 0, 0)
    Normalize(adata)
    PCA(adata, gene_mask=None, key="PCA", comp=20)

    Integrate(adata, "Slide", gene_mask=None, pca_key="PCA")
    Visualize(adata)
    Cluster(adata, neighbor_key="neighbors", cluster_key=CLUSTER_KEY, resolutions=resolutions)

    sdata.write(DATADIR / "integrated.zarr")

