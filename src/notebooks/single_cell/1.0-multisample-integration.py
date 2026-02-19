# uv run python src/notebooks/single_cell/1.0-multisample-integration.py > ./outs/scRNA_integration.out 2>&1

from pathlib import Path

import pandas as pd
import anndata as ad
import scanpy as sc
from single_cell.preprocess import (
    Convert_Genes,
    Filter_QC,
    Filter_GeneGroup,
    Normalize,
    FindVariableGenes,
    PCA,
    Integrate,
    Visualize,
)
from utils import find_barcodes, stopwatch

if __name__ == "__main__":
    ### Load & Merge
    CHECKPOINT = True
    DATADIR = Path("data") / "processed" / "single_cell"
    START = stopwatch(mode=2)

    if CHECKPOINT is False:
        #############################
        #######   Load Data   #######
        #############################
        task = "Reading all datasets"
        task_start = stopwatch(task, START, mode=0)

        adatas = {}

        annotation = "cleaned"
        study = "Emont2022"
        adatas[study] = ad.read_zarr(DATADIR / study / annotation)

        annotation = "doublet_cleaned"
        study = "So2025"
        adatas[study] = sc.read_h5ad(
            DATADIR / study / "1_preprocessed" / f"{annotation}.h5ad"
        )

        for study in [
            "Wang2025",
            "Sarvari2021",
            "GonzalezHurtado2025",
            "Jaitlin2019",
            "Kohda2025",
        ]:
            annotation = "doublet_cleaned"
            adatas[study] = sc.read_h5ad(DATADIR / study / f"{annotation}.h5ad")

        stopwatch(task, task_start, mode=1)

        #############################
        ####### Preprocessing #######
        #############################
        task = "Preprocessing (gene conversion, barcode convention, QC, normalization, hvg identification, PCA)"
        task_start = stopwatch(task, START, mode=0)

        # Map & combine unnamed genes/gene aliases
        adatas = Convert_Genes(adatas)

        # Merge
        adata = sc.concat(adatas, join="outer", label="Dataset")
        for dr in [
            "X_pca",
            "X_umap",
            "LocalMAP",
            "UMAP_doublet",
            "LocalMAP_doublet",
            "integrated",
        ]:
            del adata.obsm[dr]

        # Fix barcodes by renaming
        barcodes = find_barcodes(adata.obs_names)
        adata.obs_names = (
            barcodes
            + "_"
            + pd.Series(list(adata.obs["Dataset"].values))
            + "_"
            + pd.Series(list(adata.obs["Identifier"].values))
        ).str.replace("_", "-")

        # Begin actual preprocessing
        adata.layers["counts"] = adata.X.copy()
        adata.obs = adata.obs.dropna(axis=1, how="any")
        Filter_QC(adata)
        Filter_GeneGroup(adata, "ribo", "^Rps|^Rpl", perc_threshold=10)
        Filter_GeneGroup(adata, "mito", "^mt", perc_threshold=10)
        Filter_GeneGroup(adata, "unknown", "^ENMUSG|^Gm[0-9]|Rik")
        Normalize(adata)
        FindVariableGenes(adata, "seurat_v3", n_features=2000)

        PCA(adata, None, "PCA")
        Visualize(adata, "INT_none", obsm="global_PCA", localmap=False)
        PCA(adata)
        Visualize(adata, "INT_none_hvg", obsm="global_PCA_hvg", localmap=False)

        stopwatch(task, task_start, mode=1)

        #############################
        #######  Integration  #######
        #############################
        for integration_method in ["harmony", "scanorama", "scvi"]:
            for pca_type in ["all", "hvg"]:
                for batch_column in ["Dataset", "Identifier"]:
                    task = f"Integration Method `{integration_method}` with `{pca_type}` genes over label `{batch_column}`"
                    task_start = stopwatch(task, START, mode=0)

                    # Sort by batch to be contiguous (scanorama is picky or whateva)
                    idx = adata.obs.sort_values(batch_column).index
                    adata = adata[idx]

                    # Run integration
                    pca_key = f"global_PCA_{pca_type}"
                    int_key = f"global_INT_{integration_method}-{pca_type}-{batch_column}"
                    Integrate(
                        adata,
                        batch_column,
                        pca_key=pca_key,
                        kind=integration_method,
                        integration_key=int_key,
                    )
                    Visualize(adata, key=int_key, obsm=int_key, localmap=True, show=True)
                    
                    stopwatch(task, task_start, mode=1)

            adata.write_zarr(DATADIR / "combined" / "integrated_Outer")
    else:
        pass
