# Environment
import pickle
import numpy as np
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import umap
from harmony import harmonize
from banksy_utils.filter_utils import filter_cells, filter_hvg
from banksy_utils.umap_pca import pca_umap
from banksy.main import median_dist_to_nearest_neighbour, concatenate_all
from banksy.cluster_methods import run_Leiden_partition
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy.plot_banksy import plot_results

sc.logging.print_header()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 0  # errors (0), warnings (1), info (2), hints (3)
plt.rcParams["font.family"] = "Arial"
sns.set_style("white")
mpl.rcdefaults()


if __name__ == "__main__":
    DATADIR = Path("data")
    savedir = DATADIR / "processed" / "spatial" / "analysis" / "banksy"
    run_name = "banksy_visiumhd_vat"

    adata = sc.read(
        DATADIR / "processed" / "spatial" / "analysis" / "30_10-clustered.h5ad"
    )
    samples = list(adata.obs["Identifier"].unique())
    adata.obs["x"], adata.obs["y"] = np.hsplit(adata.obsm["SPATIAL"], 2)
    coord_keys = ("x", "y", "SPATIAL")

    # BANKSY PARAMS
    batch_key = "Identifier"
    pca_dims = [30]
    resolutions = [
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
    ]  # clustering resolutions for UMAP
    lambda_list = [0.2, 0.5]  # weight for cell vs neighborhood
    # Initialization params
    plot_graph_weights = True
    k_geom = 15  # only for fixed type
    max_m = 1  # azumithal transform up to kth order
    nbr_weight_decay = (
        "scaled_gaussian"  # can also be "reciprocal", "uniform" or "ranked"
    )

    ### START ###

    import random

    seed = 1234
    np.random.seed(seed)
    random.seed(seed)

    print("Preprocess")
    adata.obs = adata.obs.rename(columns={"pct_counts_mito": "pct_counts_mt"})
    adata.X = adata.layers["counts"].copy()
    # Filter cells with each respective filters
    adata = filter_cells(
        adata, min_count=40, max_count=1000, MT_filter=20, gene_filter=10
    )
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

    # plot and save clusters
    c_map = "tab20"  # specify color map
    weights_graph = banksy_dict["scaled_gaussian"]["weights"][0]
    plot_results(
        results_df,
        weights_graph,
        c_map,
        match_labels=True,
        coord_keys=coord_keys,
        max_num_labels=max_num_labels,
        dataset_name=run_name,
        save_path=savedir,
        save_fig=True,  # save the spatial map of all clusters
        save_seperate_fig=True,  # save the figure of all clusters plotted seperately
    )

    with open(savedir / "banksy_results.pkl", "wb") as file:
        pickle.dump((banksy_dict, results_df, max_num_labels), file)
