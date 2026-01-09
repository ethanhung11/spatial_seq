from typing import Literal, List, Iterable
from datetime import timedelta
import numpy as np
import pandas as pd
from tqdm import tqdm

from anndata import AnnData
import scanpy as sc
import scanpy.external as sce
import scvi
import lightning.pytorch as pl
from doubletdetection import BoostClassifier
from pacmap import LocalMAP
from sklearn.metrics import silhouette_samples, silhouette_score

import rpy2.robjects as ro
from .R import get_converter, R_preload


def Filter_QC(
    adata: AnnData,
    GenePerCell: int = 250,
    CountPerCell: int = 500,
    CellPerGene: int = 10,
    verbose: bool = True,
):
    start = adata.shape[0]
    sc.pp.filter_cells(adata, min_genes=GenePerCell)
    sc.pp.filter_cells(adata, min_counts=CountPerCell)
    sc.pp.filter_genes(adata, min_cells=CellPerGene)

    Filter_GeneGroup(adata, key=None, verbose=True)
    adata.var_names_make_unique()
    if verbose is True:
        print(
            f"Cells removed by gene/count filters: {start - adata.shape[0]} ({(start - adata.shape[0]) / start * 100:.2f}%)"
        )
    return adata


def Filter_GeneGroup(
    adata: AnnData,
    key: str = "mito",
    marker: str = "mt",
    verbose: bool = False,
    perc_threshold: float = None,
):
    if key is None:
        sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=[20])
    else:
        # assign gene group
        adata.var[key] = adata.var_names.str.startswith(marker)

        # calculate & save metrics
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=key, inplace=True, percent_top=[20], log1p=True
        )

    # filter
    if perc_threshold is not None:
        thresholded = adata.obs[f"pct_counts_{key}"] <= perc_threshold
        if verbose is True:
            print(
                f"Cells with >{perc_threshold}% {key} genes: {adata.shape[0] - np.sum(thresholded)} ({(adata.shape[0] - np.sum(thresholded)) / adata.shape[0] * 100:.2f}%)"
            )
        adata = adata[thresholded]

    return adata


def Filter_Doublet(
    adata: AnnData,
    method: Literal[
        "doubletdetection", "scrublet", "scDblFinder", "DoubletFinder", "SOLO"
    ] = "doubletdetection",
    remove: bool = True,
    multipletRate: float = 0.075,
):
    if "methods" not in adata.uns:
        adata.uns["methods"] = {}

    # https://github.com/mousepixels/sanbomics_scripts/blob/main/sc2024/preprocessing.ipynb
    if method == "doubletdetection":
        clf = BoostClassifier(
            n_iters=10,
            clustering_algorithm="louvain",
            standard_scaling=True,
            pseudocount=0.1,
            n_jobs=10,
        )

        clf.fit(adata.X)
        doublets = clf.predict(p_thresh=1e-3, voter_thresh=0.5)
        doublet_score = clf.doublet_score()

        adata.obs["predicted_doublet"] = doublets.astype(bool)
        adata.obs["doublet_score"] = doublet_score
    # https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html#doublet-detection
    elif method == "scrublet":
        sc.pp.scrublet(adata)
    elif method == "scDblFinder":
        R_preload()
        with ro.conversion.localconverter(get_converter()):
            tmp = adata.copy()
            del tmp.uns
            ro.globalenv["sce"] = tmp
            ro.r(
                """
            library(scater)
            library(scDblFinder)
            library(BiocParallel)
                 
            sce <- scDblFinder(sce)
            doublet_score <- sce$scDblFinder.score
            doublet_class <- sce$scDblFinder.class
            """
            )
            print(pd.Series(ro.globalenv["doublet_class"]))
            adata.obs["predicted_doublet"] = (
                pd.Series(ro.globalenv["doublet_class"])
                .replace({"singlet": 0, "doublet": 1})
                .to_list()
            )
            adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)
            adata.obs["doublet_score"] = ro.globalenv["doublet_score"]
    elif method == "DoubletFinder":
        R_preload()
        with ro.conversion.localconverter(get_converter()):
            tmp = adata.copy()
            del tmp.uns
            ro.globalenv["sce"] = tmp
            ro.globalenv["multipletRate"] = multipletRate
            ro.r(
                """
            library(scater)
            library(Seurat)
            library(dplyr)
            library(DoubletFinder)
            options(Seurat.object.assay.version = "v3")
                 
            data <- as.Seurat(sce,data = NULL)
            data <- RenameAssays(data,"originalexp","RNA")
            data <- NormalizeData(data) %>%
                FindVariableFeatures() %>%
                ScaleData()
                 
            data <- RunPCA(data, nfeatures.print = 10)
            stdv <- data[["pca"]]@stdev
            sum.stdv <- sum(data[["pca"]]@stdev)
            percent.stdv <- (stdv / sum.stdv) * 100
            cumulative <- cumsum(percent.stdv)
            co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
            co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                                percent.stdv[2:length(percent.stdv)]) > 0.1), 
                        decreasing = T)[1] + 1
            min.pc <- min(co1, co2)
            min.pc
            
            # Finish pre-processing 
            data <- RunUMAP(data, dims = 1:min.pc) %>%
                FindNeighbors(dims = 1:min.pc) %>%
                FindClusters(resolution = 0.8)

            ## 1. pK Identification (no ground-truth)
            #----------------------------------------------------------------------------
            sweep.res.list <- paramSweep(data, PCs = 1:min.pc)
            sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
            bcmvn <- find.pK(sweep.stats)
            bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
            optimal.pk <- bcmvn.max$pK
            optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
            
            # 2. Homotypic Doublet Proportion Estimate
            #----------------------------------------------------------------------------
            annotations <- data@meta.data$seurat_clusters
            homotypic.prop <- modelHomotypic(annotations)
            nExp_poi <- round(multipletRate * nrow(data@meta.data))
            nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
            
            ## 3. Run DoubletFinder
            #----------------------------------------------------------------------------
            data <- doubletFinder(data, PCs = 1:min.pc, pK = optimal.pk, nExp = nExp_poi.adj)
            colnames(data@meta.data)[grepl('DF.classifications.*', colnames(data@meta.data))] <- "doublet_finder"
            
            doublets <- data@meta.data$doublet_finder
            """
            )

            adata.obs["predicted_doublet"] = (
                pd.Series(ro.globalenv["doublets"])
                .replace({"Singlet": 0, "Doublet": 1})
                .to_list()
            )
            adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)
            adata.obs["doublet_score"] = np.nan
    elif method == "SOLO":
        scvi.model.SCVI.setup_anndata(adata)
        vae = scvi.model.SCVI(adata)
        vae.train()
        solo = scvi.external.SOLO.from_scvi_model(vae)
        solo.train(early_stopping=True)
        adata.obs["prediction"] = solo.predict(soft=True)
        adata.obs["prediction"] = (
            solo.predict(soft=False)
            .prediction.map({"singlet": 0, "doublet": 1})
            .astype(bool)
        )
    else:
        raise ValueError("Not a valid doublet detection methodtype!")

    if remove is True:
        adata.uns["methods"]["doublet_remover"] = method
        print("doublets removed: ", adata.obs["predicted_doublet"].sum())
        adata = adata[adata.obs["predicted_doublet"] == 0]

    else:
        if "doublet_remover" in adata.uns:
            adata.uns["methods"]["doublet_remover"].append(method)
        else:
            adata.uns["methods"]["doublet_remover"] = [method]
        adata.obs[f"predicted_doublet-{method}"] = adata.obs["predicted_doublet"]
        adata.obs[f"doublet_score-{method}"] = adata.obs["doublet_score"]
        del adata.obs["predicted_doublet"], adata.obs["doublet_score"]

    return adata


def Normalize(
    adata: AnnData | List[AnnData], kind: Literal["log1p", "mnn"] = "log1p", **kwargs
) -> AnnData:
    if type(adata) is list:
        adata = sc.concat(adata, join="outer")
    if "methods" not in adata.uns:
        adata.uns["methods"] = {}

    print("Saving normalized counts in layer 'normalized'.")
    adata.layers["normalized"] = adata.layers["counts"].copy()

    if kind == "log1p":
        sc.pp.normalize_total(adata, layer="normalized", target_sum=1e4)
        sc.pp.log1p(adata, layer="normalized")
        adata.uns["methods"]["normalization"] = "log1p"

    elif kind == "mnn":
        tmp = adata.copy()
        sce.pp.mnn_correct(tmp)
        adata.layers["normalized"] = tmp.X
        adata.uns["methods"]["normalization"] = "mnn"

    else:
        raise ValueError("Not a valid normalization type!")

    return adata


def FindVariableGenes(
    adata: AnnData | List[AnnData],
    kind: Literal[
        "seurat", "seurat_v3", "seurat_v3_paper", "cell_ranger", "pearson", "deviant"
    ] = "seurat",
    batch_column: str = None,
    n_features=2000,
) -> AnnData:
    print("Finding HVGs...")

    if type(adata) is list:
        adata = sc.concat(adata, join="outer")
    if "methods" not in adata.uns:
        adata.uns["methods"] = {}

    if kind == "seurat" or kind == "cell_ranger":
        sc.pp.highly_variable_genes(
            adata,
            batch_key=batch_column,
            flavor=kind,
            layer="normalized",
            n_top_genes=n_features,
        )

    elif kind == "seurat_v3" or kind == "seurat_v3_paper":
        sc.pp.highly_variable_genes(
            adata,
            layer="counts",
            batch_key=batch_column,
            flavor=kind,
            n_top_genes=n_features,
        )

    elif kind == "pearson":
        sc.experimental.pp.highly_variable_genes(
            adata,
            layer="counts",
            batch_key=batch_column,
            flavor=kind,
            n_top_genes=n_features,
        )

    elif kind == "deviant":
        R_preload()
        with ro.conversion.localconverter(get_converter()):
            tmp = adata.copy()
            del tmp.uns
            ro.globalenv["sce"] = tmp
            ro.r(
                """
            library(scry)
            devianceFeatureSelection(sce, assay='X')
            """
            )
            result = ro.r("rowData(sce)$binomial_deviance").T

        idx = result.var["binomial_deviance"].argsort()[-n_features, :]
        mask = np.zeros_like(result, dtype=bool)
        mask[idx] = True
        adata.var["binomial_deviance"] = result.var["binomial_deviance"]
        adata.var["highly_variable"] = mask

    else:
        raise ValueError("Not a valid gene selection type!")

    adata.uns["methods"]["hvg"] = kind

    return adata


def PCA(
    adata: AnnData,
    gene_mask: str = "highly_variable",
    key: str = "PCA_hvg",
    layer: str = "normalized",
    comp: int = 50,
) -> AnnData:
    print(
        f"Staring PCA with gene mask {gene_mask} at {comp} comps, saved at .obsm[{key}]"
    )

    sc.pp.pca(adata, n_comps=comp, mask_var=gene_mask, layer=layer, key_added=key)
    return adata


def Integrate(
    adata: AnnData,
    batch_column: str,
    gene_mask: str = "highly_variable",
    pca_key: str = "PCA_hvg",
    kind: Literal["harmony", "scvi", "scanorama", "seurat"] = "harmony",
    integration_key: str = "integrated",
    scvi_params={"save_model": True, "model_directory": None},
    **kwargs,
) -> AnnData:
    print(f"Integrating by Column {batch_column}: {adata.obs[batch_column].unique()}")
    if gene_mask is not None:
        assert "hvg" in adata.uns
    if "methods" not in adata.uns:
        adata.uns["methods"] = {}

    elif kind == "harmony":
        sce.pp.harmony_integrate(
            adata,
            key=batch_column,
            basis=pca_key,
            max_iter_harmony=50,
            adjusted_basis=integration_key,
            **kwargs,
        )

    elif kind == "scanorama":
        sce.pp.scanorama_integrate(
            adata,
            key=batch_column,
            basis=pca_key,
            adjusted_basis=integration_key,
            **kwargs,
        )

    elif kind == "scvi":
        if gene_mask is not None:
            input_adata = adata[:, adata.var[gene_mask]].copy()
        else:
            input_adata = adata.copy()

        scvi.model.SCVI.setup_anndata(
            input_adata, layer="counts", batch_key=batch_column
        )
        model = scvi.model.SCVI(input_adata, **kwargs)
        model.train(
            accelerator="cpu",
            strategy=pl.strategies.DDPStrategy(
                timeout=timedelta(days=2), find_unused_parameters=True
            ),
            early_stopping=True,
        )

        if scvi_params["save_model"] is True:
            model.save(scvi_params["model_directory"] / "model")
            input_adata.write(scvi_params["model_directory"] / "input.h5ad")

        adata.obsm[integration_key] = model.get_latent_representation()

    elif kind in "seurat":
        R_preload()
        with ro.conversion.localconverter(get_converter()):
            subset = adata.copy()
            del subset.uns

            ro.globalenv["sce"] = subset
            ro.r(
                """
            seurat <- as.Seurat(sce, data = NULL)
            seurat <- RenameAssays(seurat, "originalexp", "RNA")
            counts <- GetAssayData(seurat, slot = "counts")
            print(counts)
            seurat <- subset(seurat, cells = colnames(counts)[colSums(counts) > 0])
            seurat <- subset(seurat, features = rownames(counts)[rowSums(counts) > 0])
            seurat <- SCTransform(seurat)
            """
            )

            ro.globalenv["batch_key"] = batch_column
            ro.r(
                """
            batch_list <- SplitObject(seurat, split.by = batch_key)
            anchors <- FindIntegrationAnchors(batch_list, anchor.features = rownames(seurat))
            integrated <- IntegrateData(anchors)
            integrated_expr <- GetAssayData(integrated)
            integrated_expr <- integrated_expr[rownames(seurat), colnames(seurat)]
            integrated_expr <- t(integrated_expr)
            """
            )
            subset = ro.globalenv["integrated_expr"]
            adata.obsm[integration_key] = subset.X

    elif kind == "scanvi":
        raise NotImplementedError

    adata.uns["methods"]["integration"] = kind

    return adata


def Visualize(
    adata: AnnData,
    key: str = None,
    umap=True,
    input_key: str = "integrated",
    neighbor_method: Literal["umap_ann", "bbknn"] = "umap_ann",
    localmap=True,
    show=True,
    **kwargs,
):
    if "methods" not in adata.uns:
        adata.uns["methods"] = {}

    if key is None:
        key = ""
        neighbor_key = None
    else:
        key = "_" + key
        neighbor_key = "neighbors" + key

    if umap is True:
        print("Starting UMAP...")
        if neighbor_method == "umap_ann":
            sc.pp.neighbors(adata, use_rep=input_key, key_added=neighbor_key, **kwargs)
        elif neighbor_method == "bbknn":
            sc.external.pp.bbknn(
                adata, use_rep=input_key, key_added=neighbor_key, **kwargs
            )
        adata.uns["methods"]["neighbors"] = neighbor_method
        sc.tl.umap(adata, key_added=f"UMAP{key}", neighbors_key=neighbor_key)
        if show is True:
            sc.pl.embedding(adata, basis=f"UMAP{key}")

    if localmap is True:
        print("Starting LocalMAP...")
        lm = LocalMAP(**kwargs)
        if input_key is None:
            adata.obsm[f"LocalMAP{key}"] = lm.fit_transform(adata.X, init="pca")
        elif input_key:
            adata.obsm[f"LocalMAP{key}"] = lm.fit_transform(
                adata.obsm[input_key], init="pca"
            )
        if show is True:
            sc.pl.embedding(adata, basis=f"LocalMAP{key}")

    return adata


def Cluster(
    adata: AnnData,
    obs_key: str = "leiden",
    resolutions: Iterable = np.arange(5, 16) / 10,
    neighbor_key: str = "neighbors",
    random_state: int = 123,
):
    for res in tqdm(resolutions):
        try:
            sc.tl.leiden(
                adata,
                resolution=res,
                neighbors_key=neighbor_key,
                key_added=f"{obs_key}_{res}",
                random_state=random_state,
            )
        except Exception as e:
            print(f"{res} failed to run: {e}")

    return adata


def Silhouette(
    adata: AnnData,
    obsm_key="integrated",
    obs_key="leiden",
    uns_key="silhouette",
):
    X = adata.obsm[obsm_key]
    cluster_labels = adata.obs[obs_key].astype("category").cat.codes

    if uns_key not in adata.uns:
        adata.uns[f"{uns_key}_avg"] = {}
        adata.uns[f"{uns_key}_scores"] = {}

    adata.uns[f"{uns_key}_avg"][obs_key] = silhouette_score(X, cluster_labels)
    adata.uns[f"{uns_key}_scores"][obs_key] = silhouette_samples(X, cluster_labels)

    return adata
