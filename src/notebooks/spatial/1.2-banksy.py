# Environment
import pickle
import numpy as np
import pandas as pd
from IPython.display import display

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import umap
from harmony import harmonize
from banksy_utils.load_data import display_adata
from banksy_utils.plot_utils import plot_qc_hist
from banksy_utils.filter_utils import filter_cells, filter_hvg
from banksy_utils.umap_pca import pca_umap
from banksy.main import median_dist_to_nearest_neighbour, concatenate_all
from banksy.cluster_methods import run_Leiden_partition
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix

sc.logging.print_header()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 0  # errors (0), warnings (1), info (2), hints (3)
plt.rcParams["font.family"] = "Arial"
sns.set_style("white")
mpl.rcdefaults()


def concat_spatial(x, y, ref_x, ref_y, offset_x=1, offset_y=0, mode="perc"):
    """
    Offsets starting from bottom right corner (Xmax, Ymin).
    Default offset is 1 full width of reference to the right of ref_x.
    """

    if mode == "perc":
        offset_x *= ref_x.max() - ref_x.min()
        offset_y *= ref_y.max() - ref_y.min()

    return x - x.min() + ref_x.max() + offset_x, y - y.min() + ref_y.min() + offset_y


if __name__ == "__main__":
    batch_key = "Identifier"
    pca_dims = [30]
    resolutions = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]  # clustering resolutions for UMAP

    # Initialization params
    plot_graph_weights = True
    k_geom = 15  # only for fixed type
    max_m = 1  # azumithal transform up to kth order
    nbr_weight_decay = (
        "scaled_gaussian"  # can also be "reciprocal", "uniform" or "ranked"
    )
    lambda_list = [0.2]  # list of lambda parameters

    import random

    seed = 1234
    np.random.seed(seed)
    random.seed(seed)

    print("Read & Preprocess")
    from pathlib import Path

    DATADIR = Path("../../../data")

    adata = sc.read(
        DATADIR / "processed" / "spatial" / "analysis" / "30_10-clustered.h5ad"
    )
    samples = list(adata.obs["Identifier"].unique())
    adata.obs["x"], adata.obs["y"] = np.hsplit(adata.obsm["spatial"], 2)
    coord_keys = ("x", "y", "spatial")

    adata.obs = adata.obs.rename(columns={"pct_counts_mito": "pct_counts_mt"})
    adata.X = adata.layers["counts"].copy()

    # Filter cells with each respective filters
    adata = filter_cells(
        adata, min_count=40, max_count=1000, MT_filter=20, gene_filter=10
    )

    samples = adata.obs["Identifier"].unique()
    ref = samples[0]
    for n in np.arange(1, len(samples)):
        ref = samples[n - 1]
        samp = samples[n]
        ref_adata = adata[adata.obs["Identifier"] == ref].copy()
        new_adata = adata[adata.obs["Identifier"] == samp].copy()

        (
            adata.obs["x"].loc[new_adata.obs_names],
            adata.obs["y"].loc[new_adata.obs_names],
        ) = concat_spatial(
            new_adata.obs["x"],
            new_adata.obs["y"],
            ref_adata.obs["x"],
            ref_adata.obs["y"],
            offset_x=0.5,
        )

    adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy()
    sc.pl.spatial(adata, color="Identifier", spot_size=100)
    adata.X = adata.layers["counts"].copy()

    # Filter hvgs
    adata, adata_allgenes = filter_hvg(adata, n_top_genes=2000, flavor="seurat")

    print("Run Banksy")

    # Find median distance to closest neighbours, the median distance will be `sigma`
    nbrs = median_dist_to_nearest_neighbour(adata, key="spatial")
    banksy_dict = initialize_banksy(
        adata, coord_keys, k_geom, nbr_weight_decay=nbr_weight_decay
    )

    # define spatial information
    banksy_dict, banksy_matrix = generate_banksy_matrix(
        adata, banksy_dict, lambda_list, max_m
    )

    # define nonspatial information
    banksy_dict["nonspatial"] = {
        # Here we simply append the nonspatial matrix (adata.X) to obtain the nonspatial clustering results
        0.0: {
            "adata": concatenate_all([adata.X], 0, adata=adata),
        }
    }

    # run PCA
    pca_umap(banksy_dict, pca_dims=pca_dims, add_umap=True)

    # harmony of PCA results
    for pca_dim in pca_dims:
        Z = harmonize(
            banksy_dict[nbr_weight_decay][0.2]["adata"].obsm[f"reduced_pc_{pca_dim}"],
            banksy_dict[nbr_weight_decay][0.2]["adata"].obs,
            batch_key=batch_key,
        )

        print(
            f'Replacing adata.obsm["reduced_pc_{pca_dim}"] with harmony corrected embeddings.'
        )
        banksy_dict[nbr_weight_decay][0.2]["adata"].obsm[f"reduced_pc_{pca_dim}"] = Z

        # Run UMAP
        reducer = umap.UMAP(transform_seed=42)
        umap_embedding = reducer.fit_transform(Z)
        banksy_dict[nbr_weight_decay][0.2]["adata"].obsm[
            f"reduced_pc_{pca_dim}_umap"
        ] = umap_embedding

    # partition
    results_df, max_num_labels = run_Leiden_partition(
        banksy_dict,
        resolutions,
        num_nn=50,
        num_iterations=-1,
        partition_seed=seed,
        match_labels=True,
    )

    filename = DATADIR / "processed" / "spatial" / "analysis" / "banksy_results.pkl"
    with open(filename, "wb") as file:
        pickle.dump((banksy_dict, results_df, max_num_labels), file)
