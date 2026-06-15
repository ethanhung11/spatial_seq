# pixi run -e main python -u src/notebooks/spatial/joanna_VisiumHD/plots.py &> outs/spatial/plotting.out

# base
import sys
import gc
import warnings
import logging
import session_info
from tqdm import tqdm
from pathlib import Path

# data manipulation
import pandas as pd
import matplotlib.pyplot as plt

# single cell
import scanpy as sc
import spatialdata as sd
import spatialdata_plot as sdp
import spatialdata_io

warnings.simplefilter("ignore", FutureWarning)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)
warnings.simplefilter("ignore", DeprecationWarning)
warnings.simplefilter("ignore", pd.errors.DtypeWarning)
warnings.simplefilter("ignore", pd.errors.PerformanceWarning)
mlogger = logging.getLogger("matplotlib")
mlogger.setLevel(logging.ERROR)

# analysis & device specific
CORES = 10
CORES = 10
HOMEDIR = Path(".")
SRCDIR = HOMEDIR / "src"
DATADIR = HOMEDIR / "data" / "processed" / "spatial" / "VisiumHD" / "analysis"
REFDIR = HOMEDIR / "references"
PLOTDIR = HOMEDIR / "plots" / "adipose visium"

# custom
sys.path.insert(0, str(SRCDIR))
from utils import *

gc.collect()
session_info.show()


if __name__ == "__main__":
    sdata = sd.read_zarr(HOMEDIR / "data" / "processed" / "spatial" / "VisiumHD" / "raws" / "COMBINED")
    adata_segment = sc.read_h5ad(DATADIR / "30_10-clustered.h5ad")

    # clean adata_segment
    clear_adata(adata_segment, ["dendrogram", "colors", "X_pca"])
    sd.sanitize_table(adata_segment)
    sdata.tables["segmented_bins"].uns['spatialdata_attrs']['region'] = list(adata_segment.obs.region.cat.categories)
    adata_segment.uns['spatialdata_attrs'] = sdata.tables["segmented_bins"].uns['spatialdata_attrs']
    sdata.tables["segmented_bins"] = adata_segment.copy()

    # add into sdata
    cleaned_sdata_shapes, sdata.tables["segmented_bins"] = sd.match_element_to_table(
        sdata, list(adata_segment.obs.region.cat.categories), "segmented_bins"
    )
    for key in cleaned_sdata_shapes:
        sdata[key] = cleaned_sdata_shapes[key].copy()

    # rename neighborhoods
    sample = "LFD-CTR-male-670-VAT"
    old_key = 'labels_scaled_gaussian_pc30_nc0.20_r0.50'
    coloring = "Niches"
    sdata['segmented_bins'].obs["Niches"] = sdata['segmented_bins'].obs[old_key].map({
        0 : "Adipo1",
        1 : "Endothelial",
        2 : "VAM & monocyte",
        3 : "FALCs",
        4 : "Neutrophils",
        6 : "Adipo2",
        7 : "LAM remodeling",
        9 : "Mesothelial",
        11 : "Lymphatics",
        13 : "Proliferating",
        14 : "Adipo3",
        16 : "Mast cell",
        17 : "Adipo4",
        19 : "RBC/hemoglobin",
    }).astype("category")

    custom_palette = {
        "Adipo1"            : "#FFDBA1",
        "Adipo2"            : "#FFDBA1",
        "Adipo3"            : "#FFDBA1",
        "Adipo4"            : "#FFDBA1",
        "Mesothelial"       : "#0F9430",
        "Lymphatics"        : "#D17D00",
        "Endothelial"       : "#F7209D",
        "FALCs"             : "#4A0875",
        "LAM remodeling"    : "#0000FF",
        "VAM & monocyte"    : "#24F0FF",
        "Neutrophils"       : "#FFFF00",
        "Mast cell"         : "#36FFDD",
        "RBC/hemoglobin"    : "#FF3030", 
        "Proliferating"     : "#595959",
    }

    for sample in tqdm(["LFD-cKO-male-780-VAT", "HFD-cKO-male-666-VAT", "LFD-CTR-male-670-VAT"]): 
        f,ax = plt.subplots(1,1, figsize=(30,30), layout="constrained")

        sdata.pl.render_images(
            f'hires_image-{sample}',
        ).pl.render_shapes(
            f'cell_segmentation-{sample}',
            color=coloring,
            fill_alpha=0.7,
            # outline_width=0.001,
            # outline_color="black",
            palette=custom_palette,
            groups=list(sdata['segmented_bins'].obs[coloring].cat.categories),
        ).pl.show(
            coordinate_systems=sample, title="Neighborhood", ax=ax
        )

        # sdata.pl.render_images(
        #     f'hires_image-{sample}',
        # ).pl.show(
        #     coordinate_systems=sample, title="Reference", ax=axs[1]
        # )

        f.savefig(PLOTDIR / f"neighborhoods - {sample}.jpg", dpi=800)