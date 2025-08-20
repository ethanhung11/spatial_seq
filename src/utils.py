import os
import re
import pandas as pd
from tqdm import tqdm
from typing import Iterable

# single cell
from scanpy import AnnData
import spatialdata as sd
import spatialdata_plot as sdp # noqa: F401
from spatialdata_io import visium_hd

# R
import rpy2.robjects as ro
from single_cell.R import R_preload, get_converter

def clean_string(s):
    return re.sub(r"[^a-zA-Z0-9._-]", "_", s)

def create_cloupe(adata: AnnData):
    with ro.conversion.localconverter(get_converter()):
        R_preload()
        ro.globalenv["sce"] = adata
        ro.r(
            """
        library("loupeR")
        clust <- as.list(colData(sce)) %>% lapply(as.factor)
        proj <- lapply(as.list(reducedDims(sce)), function(df) df %>% select(1:2))

        create_loupe(
            assay(sce, "X"),
            clusters = clust, 
            projections = proj,
            output_name = "output",
        )
        """
        )

def GetSpatialData(
    directory: str,
    savedir: str,
    intermediate_path: str = "outs",
    exclude:Iterable=[],
    overwrite:bool=True,
    get_sdata:bool=True,
    get_custom_barcodes:bool = True,
    **kwargs,
):
    """
    Gets all spatial data from a directory.
    **kwargs are sent to `spatialdata_io.visium_hd()`
    """
    sdatas = {}
    barcodes = {}
    samples = os.listdir(directory)
    for s in exclude:
        samples.remove(s)

    for n, sample in tqdm(enumerate(samples)):
        print(sample)

        savefile = os.path.join(savedir, sample, "object.zarr")

        if get_sdata is True:
            if overwrite is True:
                sdata = visium_hd(
                    os.path.join(directory, sample, intermediate_path),
                    dataset_id=sample,
                    load_all_images=True,
                    var_names_make_unique=True,
                )
                sdata.write(savefile, overwrite=overwrite, **kwargs)

            #
            sdata = sd.read_zarr(savefile)
            for table in sdata.tables.values():
                table.obs["Identifier"] = sample
                table.layers["counts"] = table.X

            # rename shapes
            for img_name in ["cytassist_image", "hires_image", "lowres_image"]:
                sdata[img_name] = sdata[f"{sample}_{img_name}"]
                del sdata[f"{sample}_{img_name}"]

            # rename shapes & annotations
            for bin_size in ["002", "008", "016"]:
                sdata[f"SHAPE_square_{bin_size}um"] = sdata[
                    f"{sample}_square_{bin_size}um"
                ]
                del sdata[f"{sample}_square_{bin_size}um"]
                sdata[f"square_{bin_size}um"].uns["spatialdata_attrs"]["region"] = [
                    f"SHAPE_square_{bin_size}um"
                ]
                sdata[f"square_{bin_size}um"].obs[
                    "region"
                ] = f"SHAPE_square_{bin_size}um"

            sdatas[sample] = sdata

        if get_custom_barcodes is True:
            spatial_barcode_subsets = os.listdir(os.path.join(savedir, sample))
            barcodes[sample] = {}
            for csv in spatial_barcode_subsets:
                if ".csv" in csv:
                    barcodes[sample][csv] = pd.read_csv(os.path.join(savedir, sample, csv))
                    barcodes[sample][csv]["Barcode"] = (
                        barcodes[sample][csv]["Barcode"] + f"-{sample}"
                    )
                    barcodes[sample][csv].set_index("Barcode", inplace=True)

    sdatas = sd.concatenate(sdatas, concatenate_tables=True)

    return samples, sdatas, barcodes