import rpy2.robjects as ro
from .R import R_preload, get_converter
from scanpy import AnnData

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