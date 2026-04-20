# pixi run -e cellcharter python -u cellcharter.py &> ../../../../../outs/spatial/tools/cellcharter.out

# base
import sys
from pathlib import Path
import pickle
import session_info

# single cell
import scanpy as sc
import squidpy as sq
import cellcharter as cc

# # integration
# import scvi
# from lightning.pytorch import seed_everything
# seed_everything(12345)
# scvi.settings.seed = 12345

# analysis & device specific
CORES = 20
HOMEDIR = Path("../../../../..")
DATADIR = HOMEDIR / "data" / "processed" / "spatial" / "Xenium" / "kennedi_flu"
sys.path.insert(0, str(HOMEDIR / "src"))
from stopwatch import stopwatch

session_info.show()

if __name__ == "__main__":
    START = stopwatch(mode=-1)
    batch = "Sample"
    spatial_key = "spatial"

    # Prepwork
    TASKSTART = stopwatch("READ DATA", START, 0)
    adata = sc.read_h5ad(DATADIR / "integrated.h5ad")
    SAVEDIR = DATADIR / "tools" / "cellcharter"
    SAVEDIR.mkdir(parents=True, exist_ok=True)
    stopwatch("READ DATA", TASKSTART, 1)

    ### IGNORE PREPROCESSING & INTEGRATION STEPS -- ALREADY DONE ###

    # Identify spatial neighbors
    TASKSTART = stopwatch("SPATIAL NEIGHBORS", START, 0)
    sq.gr.spatial_neighbors(
        adata,
        library_key=batch,
        coord_type="generic",
        delaunay=True,
        spatial_key=spatial_key,
        percentile=99,
    )
    cc.gr.remove_long_links(adata)
    stopwatch("SPATIAL NEIGHBORS", TASKSTART, 1)

    # Run cellcharter
    TASKSTART = stopwatch("RUN CELLCHARTER", START, 0)
    cc.gr.aggregate_neighbors(
        adata, n_layers=3, use_rep="integrated", out_key="X_cellcharter", sample_key=batch
    )
    adata.write(SAVEDIR / "cellcharter.h5ad")
    stopwatch("RUN CELLCHARTER", TASKSTART, 1)

    # check MFI stability
    TASKSTART = stopwatch("CLUSTER STABILITY", START, 0)
    print("check cluster stability")
    autok = cc.tl.ClusterAutoK(n_clusters=(2, 12), max_runs=10, convergence_tol=0.001)
    autok.fit(adata, use_rep="X_cellcharter")
    with open(SAVEDIR / "autok.pkl", "wb") as file:
        pickle.dump(autok, file)
    stopwatch("CLUSTER STABILITY", TASKSTART, 1)


    stopwatch("COMPLETED SCRIPT", START, 1)
