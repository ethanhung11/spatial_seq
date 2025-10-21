import SpaGCN as spg
import scanpy as sc
from scipy.sparse import block_diag


def run_spagcn_combined(adata, batch_col="Identifier", n_clusters=7, p=0.5):
    adata = adata.copy()

    adata.obs["x_pixel"] = adata.obsm["spatial"][:, 0]
    adata.obs["y_pixel"] = adata.obsm["spatial"][:, 1]
    adata.obs["x_array"] = adata.obs["x_pixel"].copy()
    adata.obs["y_array"] = adata.obs["y_pixel"].copy()

    spg.prefilter_genes(adata, min_cells=3)
    spg.prefilter_specialgenes(adata)

    batch_labels = adata.obs[batch_col].values
    dataset_names = adata.obs[batch_col].unique()
    x_pixel = adata.obs["x_pixel"].values
    y_pixel = adata.obs["y_pixel"].values

    adj_matrices = []
    for dataset_name in dataset_names:
        mask = batch_labels == dataset_name
        x_batch = x_pixel[mask]
        y_batch = y_pixel[mask]
        adj_batch = spg.calculate_adj_matrix(x=x_batch, y=y_batch, histology=False)
        adj_matrices.append(adj_batch)

    adj = block_diag(adj_matrices, format="csr")

    l = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
    res = spg.search_res(
        adata,
        adj,
        l,
        n_clusters,
        start=0.7,
        step=0.1,
        tol=5e-3,
        lr=0.05,
        max_epochs=20,
        r_seed=100,
        t_seed=100,
        n_seed=100,
    )

    clf = spg.SpaGCN()
    clf.set_l(l)
    clf.train(
        adata,
        adj,
        init_spa=True,
        init="louvain",
        res=res,
        tol=5e-3,
        lr=0.05,
        max_epochs=200,
    )
    y_pred, prob = clf.predict()

    adata.obs["pred"] = y_pred
    adata.obs["pred"] = adata.obs["pred"].astype("category")

    refined_pred = spg.refine(
        sample_id=adata.obs.index.tolist(),
        pred=adata.obs["pred"].tolist(),
        dis=adj,
        shape="square",
    )
    adata.obs["refined_pred"] = refined_pred
    adata.obs["refined_pred"] = adata.obs["refined_pred"].astype("category")

    return adata


if __name__ == "main":
    from pathlib import Path
    import spatialdata as sd
    import spatialdata_plot as sdp

    DATADIR = Path("../../../data")

    # Load data
    sdata = sd.read_zarr(DATADIR / "processed" / "spatial" / "raws" / "COMBINED")
    adata = sdata["square_008um"].copy()

    # Preprocess
    spg.prefilter_genes(adata, min_cells=3)  # avoiding all genes are zeros
    spg.prefilter_specialgenes(adata)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    # Run SpaGCN
    combined_result = run_spagcn_combined(adata, n_clusters=10)
