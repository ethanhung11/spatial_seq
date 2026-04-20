import anndata as ad
from pathlib import Path
import notebooks.spatial.pailin_CurioTrekker.tools.run_commot as ct
from single_cell.analysis import *

from time import ctime, time
from datetime import timedelta

DATADIR = Path("./data")

if __name__ == "__main__":
    start = time()
    print(f"began script at {ctime(time())}")

    adata = ad.read_zarr(DATADIR / "processed" / "external" / "neighborhood_analyzed")

    ccc_db = liana_mouse_resource()[["ligand", "receptor"]]
    # ccc_db = ct.pp.ligand_receptor_database(species="mouse", signaling_type=None, database="CellChat")
    ccc_db_filt = ct.pp.filter_lr_database(ccc_db, adata, min_cell_pct=0.05)

    print(
        f"begin Commot at {ctime(time())}, (total: {timedelta(seconds=time()-start)})"
    )
    task_start = time()

    ctdata = ct.tl.spatial_communication(
        adata,
        database_name="custom",
        df_ligrec=ccc_db_filt,
        dis_thr=200,
        heteromeric=True,
        pathway_sum=True,
        copy=True,
    )

    print(f"finished in {timedelta(seconds=time()-task_start)}\n")

    adata.write_zarr(DATADIR / "processed" / "external" / "Commot_custom")
