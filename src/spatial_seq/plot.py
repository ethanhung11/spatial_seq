import numpy as np
import spatialdata as sd
import spatialdata_plot as sdp
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
