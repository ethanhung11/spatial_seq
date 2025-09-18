import os
import re
import pandas as pd
from tqdm import tqdm
from typing import Iterable

# single cell
import scanpy as sc
import spatialdata as sd
import spatialdata_plot as sdp  # noqa: F401
from spatialdata_io import visium_hd

# R
import rpy2.robjects as ro
from single_cell.R import R_preload, get_converter


def clear_uns(adata: sc.AnnData, search: str):
    for i in pd.Series(adata.uns.keys()):
        if search in i:
            del adata.uns[i]


def clear_obsm(adata: sc.AnnData, search: str):
    for i in pd.Series(adata.obsm.keys()):
        if search in i:
            del adata.obsm[i]


def clean_string(s):
    return re.sub(r"[^a-zA-Z0-9._-]", "_", s)


def create_cloupe(adata: sc.AnnData):
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


def GetSingleCellDataDIRECT(directory: str, name: str, filetype: str):
    if filetype == ".h5":
        adata = sc.read_10x_h5(directory + filetype)
    elif filetype == ".h5ad":
        adata = sc.read_h5ad(directory + filetype)
    elif filetype == ".mtx":
        adata = sc.read_10x_mtx(directory)
    else:
        raise ValueError(
            f"{filetype} is not a valid filetype. Must be h5, h5ad, or mtx."
        )

    adata.obs["Identifier"] = name
    adata.obs.index = adata.obs.index + "_" + name
    adata.layers["counts"] = adata.X.copy()

    return adata


def GetSingleCellData(
    data_dir: str,
    experiment_dir: str,
    filename: str,
    filetype: str,
    convertR: bool = False,
):
    if convertR is True:
        basedir = os.path.join(data_dir, experiment_dir)
        with ro.conversion.localconverter(get_converter()):
            ro.globalenv["datadir"] = basedir
            ro.globalenv["filename"] = filename
            ro.globalenv["filetype"] = filetype

            print(basedir, filename, filetype)

            R_preload()
            ro.r(
                """
            datadir <- here(datadir)
            experiments <- list.files(datadir)
            adatas <- c()
            sces <- c()
            for (exp in experiments) {
                if (filetype==".h5") {
                    filedir <- paste(datadir, exp, paste0(filename,filetype), sep="/")
                    print(filedir)
                    seurat_data <- Read10X_h5(filedir, use.names = T, unique.features=T)
                } else if (filetype==".mtx") {
                    filedir <- paste(datadir, exp, filename, sep="/")
                    print(filedir)
                    seurat_data <- Read10X(filedir, unique.features=T)
                }
                
                seurat_obj <- CreateSeuratObject(counts = seurat_data, project = exp)
                sce <- convert_seurat_to_sce(seurat_obj)
                sces <- c(sces,sce)
            }
            """
            )
            adatas = list(ro.globalenv["sces"])

        for ad in adatas:
            ad.obs["Identifier"] = ad.obs["orig.ident"]
            ad.obs.index = ad.obs.index + "_" + ad.obs["orig.ident"].astype(str)
            ad.layers["counts"] = ad.X.copy()
            del ad.uns

    else:
        basedir = os.path.join(data_dir, experiment_dir)
        samples = os.listdir(basedir)
        filt_files = [os.path.join(basedir, s, filename) for s in samples]
        adatas = [
            GetSingleCellDataDIRECT(file, samples[n], filetype)
            for n, file in enumerate(filt_files)
        ]

    return adatas


def GetSpatialData(
    directory: str,
    savedir: str,
    intermediate_path: str = "outs",
    exclude: Iterable = [],
    overwrite: bool = True,
    get_sdata: bool = True,
    get_custom_barcodes: bool = True,
    **kwargs,
):
    """
    Gets all spatial data from a directory.
    **kwargs are sent to `spatialdata_io.visium_hd()`
    Run sdatas = sd.concatenate(sdatas, concatenate_tables=True) to concatenate
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
                    barcodes[sample][csv] = pd.read_csv(
                        os.path.join(savedir, sample, csv)
                    )
                    barcodes[sample][csv]["Barcode"] = (
                        barcodes[sample][csv]["Barcode"] + f"-{sample}"
                    )
                    barcodes[sample][csv].set_index("Barcode", inplace=True)

    # sdatas = sd.concatenate(sdatas, concatenate_tables=True)

    return samples, sdatas, barcodes
