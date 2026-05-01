import numpy as np
import spatialdata as sd
import spatialdata_plot as sdp
from tqdm import tqdm
from typing import Iterable


def rasterize_to_plot(
    sdata: sd.SpatialData,
    samples: Iterable,
    sample_col: str,
    table: str,
    shapes: Iterable[str],
):
    for n, samp in enumerate(samples):
        sdata["raster"] = sdata[table][sdata[table].obs[sample_col] == samp].copy()
        sdata["raster"].X = sdata["raster"].X.tocsc()
        rasterized = sd.rasterize_bins(
            sdata,
            bins=shapes[n],
            table_name="raster",
            col_key="array_col",
            row_key="array_row",
        )
        sdata[f"rasterized_{samp}"] = rasterized
    del sdata["raster"]
    return sdata


def spatial_batchFlatten(adata, batch, old_key="spatial", new_key="SPATIAL"):
    assert adata.obs[batch].dtype == "category"

    adata.obs["x"], adata.obs["y"] = np.hsplit(adata.obsm[old_key], 2)

    for n, df in adata.obs.groupby(batch):
        adata.obs.loc[adata.obs[batch] == n, "x"] -= (
            df["x"].min() - adata.obs["x"].min()
        )
        adata.obs.loc[adata.obs[batch] == n, "y"] -= (
            df["y"].min() - adata.obs["y"].min()
        )

    adata.obsm[new_key] = adata.obs[["x", "y"]].to_numpy()


def spatial_batchMove(
    adata, batch, x_offset_perc, y_offset_perc, old_key="spatial", new_key="SPATIAL"
):
    assert adata.obs[batch].dtype == "category"

    adata.obs["x"], adata.obs["y"] = np.hsplit(adata.obsm[old_key], 2)

    x_range = np.mean(
        [(df["x"].max() - df["x"].min()) for _, df in adata.obs.groupby(batch)]
    )
    for i, t in enumerate(adata.obs[batch].cat.categories):
        adata.obs.loc[adata.obs[batch] == t, "x"] += x_offset_perc * x_range * i + 1

    y_range = np.mean(
        [(df["y"].max() - df["y"].min()) for _, df in adata.obs.groupby(batch)]
    )
    for i, t in enumerate(adata.obs[batch].cat.categories):
        adata.obs.loc[adata.obs[batch] == t, "y"] += y_offset_perc * y_range * i + 1

    adata.obsm[new_key] = adata.obs[["x", "y"]].to_numpy()


def spatial_label(adata, spatial_key, sample_key, ax, fs=None):
    if fs is None:
        fs = (adata.obsm[spatial_key].max(axis=0)[0] - adata.obsm[spatial_key].min(axis=0)[0]) / 250
    ax.legend(fontsize=fs, markerscale=fs/10)
    ax.set_title(ax.get_title(), fontsize=fs*1.2)
    for sample in tqdm(adata.obs[sample_key].cat.categories):
        coords = adata.obsm[spatial_key][adata.obs[sample_key] == sample]
        xmax, ymax = coords.max(axis=0)
        xmin, ymin = coords.min(axis=0)
        x = (xmin + xmax)/2
        y = ymax + (ymax - ymin)*0.03
        ax.text(x, y, sample, fontsize=fs, ha="center")
    ax.set_ylim(top=ax.get_ylim()[1]*1.03)
    return