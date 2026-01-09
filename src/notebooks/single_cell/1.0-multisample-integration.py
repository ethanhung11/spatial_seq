# taskset -c 30-39 uv run python src/notebooks/single_cell/1.0-multisample-integration.py > ./outs/scRNA_integration.out 2>&1

from pathlib import Path
from time import time, ctime
from datetime import timedelta
import re

import pandas as pd
import anndata as ad
import scanpy as sc
from single_cell.preprocess import (
    Filter_QC,
    Normalize,
    FindVariableGenes,
    PCA,
    Integrate,
    Visualize,
)

DATADIR = Path("data")


def find_barcodes(strings):
    pattern = r"[ATCG]{16}"
    matches = [re.findall(pattern, s) for s in strings]
    return pd.Series([m for sublist in matches for m in sublist]).astype(str)


if __name__ == "__main__":
    ### Load & Merge
    savedir = DATADIR / "processed" / "single cell" / "combined"
    start = time()

    ### Load Data
    print(
        f"begin Reading Data at {ctime(time())}, (total: {timedelta(seconds=time()-start)})"
    )
    task_start = time()
    adatas = {}

    annotation = "cleaned"
    study = "Emont2022"
    folder = DATADIR / "processed" / "single cell" / study
    adatas[study] = ad.read_zarr(folder / annotation)

    annotation = "doublet_cleaned"
    study = "Sarvari2021"
    folder = DATADIR / "processed" / "single cell" / study
    adatas[study] = sc.read_h5ad(folder / f"{annotation}.h5ad")

    annotation = "doublet_cleaned"
    study = "GonzalezHurtado2025"
    folder = DATADIR / "processed" / "single cell" / study
    adatas[study] = sc.read_h5ad(folder / f"{annotation}.h5ad")

    annotation = "doublet_cleaned"
    study = "So2025"
    folder = DATADIR / "processed" / "single cell" / study / "1_preprocessed"
    adatas[study] = sc.read_h5ad(folder / f"{annotation}.h5ad")

    print(f"finished in {timedelta(seconds=time()-task_start)}\n")

    print(
        f"begin Preprocessing at {ctime(time())}, (total: {timedelta(seconds=time()-start)})"
    )
    task_start = time()

    adata = sc.concat(adatas, join="outer", label="Dataset")
    adata.layers["counts"] = adata.X.copy()
    adata.obs = adata.obs.dropna(axis=1, how="any")

    barcodes = find_barcodes(adata.obs_names)
    adata.obs_names = (
        barcodes
        + "_"
        + pd.Series(list(adata.obs["Dataset"].values))
        + "_"
        + pd.Series(list(adata.obs["Identifier"].values))
    ).str.replace("_", "-")

    adata = Filter_QC(adata)
    adata = Normalize(adata)
    adata = FindVariableGenes(adata, "seurat_v3", n_features=2000)

    adata = PCA(adata, None, "PCA")
    Visualize(adata, "INT_none", input_key="PCA", localmap=False)

    PCA(adata)
    Visualize(adata, "INT_none_hvg", input_key="PCA_hvg", localmap=False)

    print(f"finished in {timedelta(seconds=time()-task_start)}\n")

    print(
        f"begin Integration at {ctime(time())}, (total: {timedelta(seconds=time()-start)})"
    )
    task_start = time()

    pca_key = "PCA"

    for integration_method in ["harmony", "scanorama"]:

        for batch_column in ["Dataset", "Identifier"]:
            pca_key = "PCA_hvg"
            int_key = f"INT_{integration_method}_hvg-{batch_column}"
            Integrate(
                adata,
                batch_column,
                pca_key=pca_key,
                kind=integration_method,
                integration_key=int_key,
            )
            Visualize(adata, key=int_key, input_key=int_key, localmap=True)

            pca_key = "PCA"
            int_key = f"INT_{integration_method}-{batch_column}"
            Integrate(
                adata,
                batch_column,
                pca_key=pca_key,
                kind=integration_method,
                integration_key=int_key,
                gene_mask=None,
            )
            Visualize(adata, key=int_key, input_key=int_key, localmap=True)

    print(f"finished in {timedelta(seconds=time()-task_start)}\n")

    adata.write_zarr(
        DATADIR / "processed" / "single cell" / "combined" / "integrated_Outer"
    )


#     SCVI_RUN = True
#     merge_column = "Dataset"
#     savedir = DATADIR / "processed" / "single cell" / "combined"
#     USE_VAR_GENES = "highly_variable"

#     start = time()
#     print(f"begin script at {ctime(time())}")


#     ### Load Data
#     print(f"begin reading in datasets at {ctime(time())}, (total: {timedelta(seconds=time()-start)})")
#     task_start = time()
#     adatas = {}

#     annotation = "cleaned"
#     study = "Emont2022"
#     folder = DATADIR / "processed" / "single cell" / study
#     adatas[study] = ad.read_zarr(folder / annotation)

#     annotation = "doublet_cleaned"
#     study = "Sarvari2021"
#     folder = DATADIR / "processed" / "single cell" / study
#     adatas[study] = sc.read_h5ad(folder / f"{annotation}.h5ad")

#     annotation = "doublet_cleaned"
#     study = "GonzalezHurtado2025"
#     folder = DATADIR / "processed" / "single cell" / study
#     adatas[study] = sc.read_h5ad(folder / f"{annotation}.h5ad")

#     annotation = "doublet_cleaned"
#     study = "So2025"
#     folder = DATADIR / "processed" / "single cell" / study / "1_preprocessed"
#     adatas[study] = sc.read_h5ad(folder / f"{annotation}.h5ad")

#     print(f"finished in {timedelta(seconds=time()-task_start)}\n")
#     print(f"begin initial preprocessing at {ctime(time())}, (total: {timedelta(seconds=time()-start)})")
#     task_start = time()

#     print(f"finished in {timedelta(seconds=time()-task_start)}\n")


#     ### Preprocessing
#     print(f"begin concatenation & preprocessing at {ctime(time())}, (total: {timedelta(seconds=time()-start)})")
#     task_start = time()

#     adata = sc.concat(adatas, join="outer", label="Dataset")
#     adata.obs = adata.obs.dropna(axis=1, how='any')

#     barcodes = find_barcodes(adata.obs_names)
#     adata.obs_names = (barcodes + "_" + pd.Series(list(adata.obs["Dataset"].values)) + "_" + pd.Series(list(adata.obs["Identifier"].values))).str.replace("_", "-")

#     Filter_QC(adata)
#     Normalize(adata)
#     FindVariableGenes(adata, "seurat_v3", n_features=2000)
#     PCA(adata, "highly_variable", "PCA_hvg")
#     PCA(adata, None, "PCA")

#     adata.write_zarr(savedir / "integrated_Dataset")
#     print(f"finished in {timedelta(seconds=time()-task_start)}\n")


#     ### Integration

#     # # Seurat
#     # integration_method = "seurat"
#     # print(f"begin {integration_method} integration at {ctime(time())}, (total: {timedelta(seconds=time()-start)})")
#     # Integrate(
#     #     adata,
#     #     "Identifier",
#     #     kind=integration_method,
#     #     integration_key=f"integrated_{integration_method}",
#     #     gene_mask=USE_VAR_GENES)
#     # Visualize(
#     #     adata,
#     #     integration_method,
#     #     input_key=f"integrated_{integration_method}")
#     # print(f"finished in {timedelta(seconds=time()-task_start)}\n")


#     ## BY DATASET ##
#     for batch_column in ["Dataset"]:

#         # Harmony Integration
#         integration_method = "harmony"
#         print(f"begin {integration_method} integration at {ctime(time())}, (total: {timedelta(seconds=time()-start)})")
#         Integrate(
#             adata,
#             batch_column,
#             kind=integration_method,
#             integration_key=f"integrated_{integration_method}_{batch_column}",
#             pca_key="PCA_hvg")
#         Visualize(
#             adata,
#             f"{integration_method}_{batch_column}",
#             input_key=f"integrated_{integration_method}_{batch_column}")
#         print(f"finished in {timedelta(seconds=time()-task_start)}\n")


#         ### SCVI Integration
#         integration_method = "scvi"
#         SCVI_RUN = False if batch_column == "Dataset" else True
#         print(f"begin {integration_method} integration at {ctime(time())}, (total: {timedelta(seconds=time()-start)})")
#         print(f"SCVI_RUN set to {SCVI_RUN}")

#         # Run SCVI
#         task_start = time()
#         if SCVI_RUN is True:
#             Integrate(
#                 adata,
#                 batch_column,
#                 kind=integration_method,
#                 integration_key=f"integrated_{integration_method}_{batch_column}",
#                 gene_mask="highly_variable",
#                 scvi_params={"save_model" : True, "model_directory" : savedir / f"SCVI_{batch_column}"})

#             Visualize(
#                 adata,
#                 f"{integration_method}_{batch_column}",
#                 input_key=f"integrated_{integration_method}_{batch_column}")
#             print("scvi model successfully trained!")
#         else:
#             input_adata = sc.read_h5ad(savedir / f"SCVI_{batch_column}" / "input.h5ad")
#             model = scvi.model.SCVI.load(savedir / f"SCVI_{batch_column}" / "model", input_adata)
#             adata.obsm[f"integrated_{integration_method}_{batch_column}"] = model.get_latent_representation()
#             print("scvi model successfully loaded!")

#             Visualize(
#                 adata,
#                 f"{integration_method}_{batch_column}",
#                 input_key=f"integrated_{integration_method}_{batch_column}")

#         del input_adata, model
#         print(f"finished in {timedelta(seconds=time()-task_start)}\n")

#         adata.write_zarr(savedir / "integrated_Dataset")


#     ## WITHOUT HVG, BY IDENTIFIER ##
#     batch_column = "Identifier"

#     # Harmony Integration
#     integration_method = "harmony"
#     print(f"begin {integration_method} integration at {ctime(time())}, (total: {timedelta(seconds=time()-start)})")
#     Integrate(
#         adata,
#         batch_column,
#         kind=integration_method,
#         integration_key=f"integrated_{integration_method}_{batch_column}_hvgOFF",
#         gene_mask=None)
#     Visualize(
#         adata,
#         integration_method,
#         input_key=f"integrated_{integration_method}_{batch_column}_hvgOFF")
#     print(f"finished in {timedelta(seconds=time()-task_start)}\n")

#     # SCVI Integration
#     integration_method = "scvi"
#     print(f"begin {integration_method} integration at {ctime(time())}, (total: {timedelta(seconds=time()-start)})")
#     print(f"SCVI_RUN set to {SCVI_RUN}")

#     # Run SCVI
#     task_start = time()
#     if SCVI_RUN is True:
#         Integrate(
#             adata,
#             batch_column,
#             kind=integration_method,
#             integration_key=f"integrated_{integration_method}_{batch_column}",
#             gene_mask=USE_VAR_GENES,
#             scvi_params={"save_model" : True, "model_directory" : savedir / f"SCVI_{batch_column}"})

#         Visualize(adata, integration_method, input_key=f"integrated_{integration_method}_{batch_column}_hvgOFF")
#         print("scvi model successfully trained!")
#     else:
#         input_adata = sc.read_h5ad(savedir / f"SCVI_{batch_column}" / "input.h5ad")
#         model = scvi.model.SCVI.load(savedir / f"SCVI_{batch_column}" / "model", input_adata)
#         adata.obsm[f"integrated_{integration_method}_{batch_column}_hvgOFF"] = model.get_latent_representation()
#         print("scvi model successfully loaded!")
#         del input_adata, model

#     Visualize(
#         adata,
#         integration_method,
#         input_key=f"integrated_{integration_method}")
#     print(f"finished in {timedelta(seconds=time()-task_start)}\n")

#     # ### Cluster
#     # resolutions = np.arange(1, 25) / 20 # 0.05, 0.1, ... 1.20
#     # for method in ["harmony", "scvi", "seurat"]:
#     #     Cluster(adata, f"global_leiden_{method}", resolutions=resolutions, neighbor_key=f"neighbors_{method}")

#     ### Saving
#     print(f"begin saving at {ctime(time())}, (total: {timedelta(seconds=time()-start)})")
#     task_start = time()
#     adata.write_zarr(savedir / "integrated_Dataset")
#     print(f"finished in {timedelta(seconds=time()-task_start)}\n")
