import spatialdata as sd
import spatialdata_plot as sdp
from typing import Iterable

def rasterize_to_plot(sdata:sd.SpatialData, samples:Iterable, sample_col:str, table:str, shapes:Iterable[str]):
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