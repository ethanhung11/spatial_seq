import re
import gzip
import random
import pandas as pd
from pathlib import Path
from typing import Iterable, Literal

# single cell
import scanpy as sc
from geosketch import gs
from scsampler import scsampler
import rpy2.robjects as ro
from single_cell.R import R_preload, get_converter


def rename_uns(adata: sc.AnnData, pat: str, prefix: str=None, suffix: str=None,
                replace: str=None, filter: str=None):
    if (prefix or suffix) and replace:
        raise ValueError("Use prefix/suffix OR replace.")
    if not (prefix or suffix or replace):
        raise ValueError("Nothing to do.")

    for k in list(adata.uns.keys()):
        if bool(re.search(pat, k)) and (filter is None or not bool(re.search(filter, k))):
            if replace:
                new = re.sub(pat, replace, k)
            else:
                new = f"{prefix or ''}{k}{suffix or ''}"

            if new != k and not(new in adata.uns):
                print(f"renamed in uns: {k} --> {new}")
                adata.uns[new] = adata.uns[k].copy()
                del adata.uns[k]

    return True


def rename_obsm(adata: sc.AnnData, pat: str, prefix: str=None, suffix: str=None,
                replace: str=None, filter: str=None):
    if (prefix or suffix) and replace:
        raise ValueError("Use prefix/suffix OR replace.")
    if not (prefix or suffix or replace):
        raise ValueError("Nothing to do.")

    for k in list(adata.obsm.keys()):
        if bool(re.search(pat, k)) and (filter is None or not bool(re.search(filter, k))):
            new = k
            if replace:
                new = re.sub(pat, replace, k)
            else:
                new = f"{prefix or ''}{k}{suffix or ''}"

            if new != k and not(new in adata.obsm):
                print(f"renamed in obsm: {k} --> {new}")
                adata.obsm[new] = adata.obsm[k].copy()
                del adata.obsm[k]

    return True


def rename_obsp(adata: sc.AnnData, pat: str, prefix: str=None, suffix: str=None,
                replace: str=None, filter: str=None):
    if (prefix or suffix) and replace:
        raise ValueError("Use prefix/suffix OR replace.")
    if not (prefix or suffix or replace):
        raise ValueError("Nothing to do.")

    for k in list(adata.obsp.keys()):
        if bool(re.search(pat, k)) and (filter is None or not bool(re.search(filter, k))):
            if replace:
                new = re.sub(pat, replace, k)
            else:
                new = f"{prefix or ''}{k}{suffix or ''}"

            if new != k and not(new in adata.obsp):
                print(f"renamed in obsp: {k} --> {new}")
                adata.obsp[new] = adata.obsp[k].copy()
                del adata.obsp[k]

    return True


def rename_varm(adata: sc.AnnData, pat: str, prefix: str=None, suffix: str=None,
                replace: str=None, filter: str=None):
    if (prefix or suffix) and replace:
        raise ValueError("Use prefix/suffix OR replace.")
    if not (prefix or suffix or replace):
        raise ValueError("Nothing to do.")
    
    for k in list(adata.varm.keys()):
        if bool(re.search(pat, k)) and (filter is None or not bool(re.search(filter, k))):
            if replace:
                new = re.sub(pat, replace, k)
            else:
                new = f"{prefix or ''}{k}{suffix or ''}"

            if new != k and not(new in adata.varm):
                print(f"renamed in varm: {k} --> {new}")
                adata.varm[new] = adata.varm[k].copy()
                del adata.varm[k]

    return True


def clear_obs(adata: sc.AnnData, pat: str):
    print(
        "deleted from `obs`: ",
        adata.obs.columns[adata.obs.columns.str.contains(pat, regex=True)],
    )
    adata.obs = adata.obs.loc[:, ~adata.obs.columns.str.contains(pat, regex=True)]


def clear_var(adata: sc.AnnData, pat: str):
    print(
        "deleted from `var`: ",
        adata.var.columns[adata.var.columns.str.contains(pat)],
    )
    adata.var = adata.var.loc[:, ~adata.var.columns.str.contains(pat)]


def clear_uns(adata: sc.AnnData, pat: str, filter: str=None):
    for key in pd.Series(adata.uns.keys()):
        if bool(re.search(pat, key)) and (filter is None or not bool(re.search(filter, key))):
            print("deleted from `uns`: ", key)
            del adata.uns[key]


def clear_obsm(adata: sc.AnnData, pat: str, filter: str=None):
    for key in pd.Series(adata.obsm.keys()):
        if bool(re.search(pat, key)) and (filter is None or not bool(re.search(filter, key))):
            print("deleted from `obsm`: ", key)
            del adata.obsm[key]


def clear_obsp(adata: sc.AnnData, pat: str, filter: str=None):
    for key in pd.Series(adata.obsp.keys()):
        if bool(re.search(pat, key)) and (filter is None or not bool(re.search(filter, key))):
            print("deleted from `obsp`: ", key)
            del adata.obsp[key]


def clear_varm(adata: sc.AnnData, pat: str, filter: str=None):
    for key in pd.Series(adata.varm.keys()):
        if bool(re.search(pat, key)) and (filter is None or not bool(re.search(filter, key))):
            print("deleted from `varm`: ", key)
            del adata.varm[key]


def clear_adata(adata: sc.AnnData, search: str, filter=None):
    if (
        isinstance(search, str)
        or isinstance(search, int)
        or isinstance(search, float)
        or isinstance(search, bool)
    ):
        search = [str(search)]
    else:
        search = list(search)
    print(search)

    for term in search:
        clear_obs(adata, term)
        clear_var(adata, term)
        clear_uns(adata, term, filter)
        clear_obsm(adata, term, filter)
        clear_obsp(adata, term, filter)
        clear_varm(adata, term, filter)
        print("\n")


def obs_to_obsm(adata: sc.AnnData, pat: str, mode: Literal[1,-1]=1):
    if mode == 1:
        adata.obsm["HotspotModule"] = adata.obs[adata.obs.columns[adata.obs.columns.str.contains("HotspotModule")]].copy()
        clear_obs(adata, "HotspotModule")
    elif mode == -1:
        adata.obs = pd.merge(adata.obs, adata.obsm[pat], left_index=True, right_index=True)
    else:
        raise ValueError("mode must be 1 (to obsm) or -1 (to obs)")


def clean_string(s):
    return re.sub(r"[^a-zA-Z0-9._]", "_", s)


def generate_barcodes(
    N,
    cellranger_installation=Path.home() / "apps" / "cellranger-9.0.1",
    version="3M-3pgex-may-2023_TRU.txt.gz",
):
    with gzip.open(
        Path(cellranger_installation)
        / "lib"
        / "python"
        / "cellranger"
        / "barcodes"
        / version,
        "rt",
    ) as f:
        valid_barcodes = [line.strip() for line in f]
    return random.sample(valid_barcodes, k=N)


def validate_barcodes(
    barcodes,
    cellranger_installation=Path.home() / "apps" / "cellranger-9.0.1",
    version="3M-3pgex-may-2023_TRU.txt.gz",
):
    barcodes = set(barcodes)
    with gzip.open(
        Path(cellranger_installation)
        / "lib"
        / "python"
        / "cellranger"
        / "barcodes"
        / version,
        "rt",
    ) as f:
        valid_barcodes = set([line.strip() for line in f])

    valid = barcodes.intersection(valid_barcodes)
    not_valid = barcodes.difference(valid_barcodes)
    print(f"valid: {len(valid)}\nnot valid: {len(not_valid)}")
    return valid, not_valid


def find_barcodes(strings):
    matches = [re.findall(r"[ATCG]{16}", s) for s in strings]
    return pd.Series([m for sublist in matches for m in sublist]).astype(str)


def create_cloupe(
    adata: sc.AnnData,
    layer: str,
    output: str,
    obs: Iterable[str],
    obsm: Iterable[str],
    new_barcodes=None,
):
    output = Path(output)
    transfer = sc.AnnData(
        X=adata.layers[layer],
        obs=adata.obs[obs],
        var=adata.var.drop(columns=adata.var.columns),
        obsm={dr: adata.obsm[dr] for dr in obsm},
    )
    if new_barcodes is not None:
        if isinstance(new_barcodes, str):
            transfer.obs_names = adata.obs[new_barcodes]
        else:
            transfer.obs_names = new_barcodes
    print(transfer.obs_names)
    with ro.conversion.localconverter(get_converter()):
        R_preload()
        ro.globalenv["sce"] = transfer
        ro.globalenv["output_dir"] = str(output.parent)
        ro.globalenv["output"] = str(output.name)
        ro.r(r"library(loupeR)")
        ro.r(r"clust <- as.list(colData(sce)) %>% lapply(as.factor)")
        ro.r(
            r"proj <- sapply(reducedDims(sce), function(x) x[, 1:2], simplify = FALSE) %>% as.list()"
        )

        ro.r(r"""
        validate_barcodes(colnames(sce))$success %>% stopifnot()
        create_loupe(
            assay(sce, "X"),
            clusters = clust, 
            projections = proj,
            output_dir = here(output_dir),
            output_name = output,
        )
        """)

    return True


def Downsample(adata: sc.AnnData, key: str, N: int = 10000, method: Literal["geosketch", "scsampler"]="geosketch", **kwargs):
    if method == "geosketch":
        sketch_index = gs(adata.obsm[key], N, replac = False, **kwargs)
        return adata[sketch_index]
    elif method == "scsampler":
        return scsampler(adata, n_obs = N, copy = True, obsm = key, **kwargs)
    else:
        raise ValueError(f"`{method}` is not a valid method type! Must be among `geosketch`, `scsampler`")


def concat_spatial(x, y, ref_x, ref_y, offset_x=0, offset_y=0, mode="perc"):
    """
    Offsets starting from bottom right corner (Xmax, Ymin).
    Default offset is 1 full width of reference to the right of ref_x.
    """

    if mode == "perc":
        offset_x *= ref_x.max() - ref_x.min()
        offset_y *= ref_y.max() - ref_y.min()

    return x - x.min() + ref_x.max() + offset_x, y - y.min() + ref_y.min() + offset_y
