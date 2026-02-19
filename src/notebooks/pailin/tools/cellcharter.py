# base
from pathlib import Path
import pickle
import numpy as np

# single cell
import scanpy as sc
import squidpy as sq
import cellcharter as cc
import scvi

# analysis & device specific
CORES = 20
DATADIR = Path("../../../data")
REFDIR = Path("../../../references")

import scvi
from lightning.pytorch import seed_everything

seed_everything(12345)
scvi.settings.seed = 12345

if __name__ == "__main__":
    print("load data")
    adata = sc.read_h5ad(DATADIR / "processed" / "external" / "initial_clean.h5ad")
    savedir = DATADIR / "processed" / "external" / "cellcharter"
    batch = "sample"
    SCVI_RUN = True

    # scale by sample
    adata.X = adata.layers["counts"].copy().toarray()
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    # SCVI correction
    print("start SCVI")
    if SCVI_RUN is False:
        scvi.model.SCVI.setup_anndata(
            adata,
            layer="counts",
            batch_key=batch,
        )
        model = scvi.model.SCVI(adata)
        model.train(early_stopping=True, enable_progress_bar=True)
        model.save(savedir / "SCVI")
    else:
        model = scvi.model.SCVI.load(savedir / "SCVI", adata)

    adata.obsm["X_scVI"] = model.get_latent_representation(adata).astype(np.float32)

    # start cellcharter
    print("start cellcharter")
    sq.gr.spatial_neighbors(
        adata,
        library_key=batch,
        coord_type="generic",
        delaunay=True,
        spatial_key="SPATIAL",
        percentile=99,
    )
    cc.gr.remove_long_links(adata)
    cc.gr.aggregate_neighbors(
        adata, n_layers=3, use_rep="X_scVI", out_key="X_cellcharter", sample_key=batch
    )
    adata.write(savedir / "cellcharter.h5ad")

    # unknown # of clusters, so check MFI stability
    print("check cluster stability")
    autok = cc.tl.ClusterAutoK(n_clusters=(2, 12), max_runs=10, convergence_tol=0.001)
    autok.fit(adata, use_rep="X_cellcharter")

    with open(savedir / "autok.pkl", "wb") as file:
        pickle.dump(autok, file)

    print("Completed run!")
