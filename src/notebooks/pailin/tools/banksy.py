# Environment
import random
from pathlib import Path
import pickle
import numpy as np

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
from banksy.plot_banksy import plot_results

DATADIR = Path("../../../data")
seed = 1234
np.random.seed(seed)
random.seed(seed)

if __name__ == "__main__":
    batch_key = "sample"
    pca_dims = [30]
    resolutions = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]  # UMAP resolutions
    lambda_list = [0.2, 0.5]  # weight for cell vs neighborhood

    # Initialization params
    plot_graph_weights = True
    k_geom = 15  # only for fixed type
    max_m = 1  # azumithal transform up to kth order
    nbr_weight_decay = (
        "scaled_gaussian"  # can also be "reciprocal", "uniform" or "ranked"
    )

    print("Read & Preprocess")

    adata = sc.read_h5ad(DATADIR / "processed" / "external" / "initial_clean.h5ad")
    samples = list(adata.obs[batch_key].unique())
    coord_keys = ("x", "y", "SPATIAL")

    # Filter HVGs
    adata, adata_allgenes = filter_hvg(adata, n_top_genes=2000, flavor="seurat")

    print("Run Banksy")

    # Find median distance to closest neighbours, the median distance will be `sigma`
    nbrs = median_dist_to_nearest_neighbour(adata, key=coord_keys[2])
    banksy_dict = initialize_banksy(
        adata,
        coord_keys,
        k_geom,
        nbr_weight_decay=nbr_weight_decay,
        plt_edge_hist=False,
        plt_nbr_weights=False,
        plt_agf_angles=False,
        plt_theta=False,
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

    c_map = "tab20"  # specify color map
    weights_graph = banksy_dict["scaled_gaussian"]["weights"][0]

    plot_results(
        results_df,
        weights_graph,
        c_map,
        match_labels=True,
        coord_keys=coord_keys,
        max_num_labels=max_num_labels,
        save_path=DATADIR / "processed" / "external" / "banksy",
        save_fig=True,  # save the spatial map of all clusters
        save_seperate_fig=True,  # save the figure of all clusters plotted seperately
    )

    filename = DATADIR / "processed" / "external" / "banksy" / "banksy_results.pkl"
    with open(filename, "wb") as file:
        pickle.dump((banksy_dict, results_df, max_num_labels), file)
