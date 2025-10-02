import os
import pandas as pd
import tangram as tg
import scanpy as sc
import numpy as np

if __name__ == "__main__":
    DATADIR = "data"

    annotation = "annotated-ccc"
    adata_sc = sc.read_h5ad(
        os.path.join(
            DATADIR, "processed", "single cell", "5_analysis", f"{annotation}.h5ad"
        )
    )
    print("sc data read\n")

    # Use DEGs are markers
    DEGlist = sc.get.rank_genes_groups_df(
        adata_sc,
        group=None,
        key="de_all",
        pval_cutoff=1e-50,
        log2fc_min=0.1,
    ).sort_values("pvals_adj")
    DEGlist = pd.concat([df.head(200) for g, df in DEGlist.groupby("group")])
    markers = DEGlist.names.unique()
    print(len(markers), "\n")

    annotation = "raw_processed-V5 WT_filter70"
    adata_sp = sc.read_h5ad(
        os.path.join(DATADIR, "processed", "spatial", "combined", f"{annotation}.h5ad")
    )
    print("spatial data read\n")

    tg.pp_adatas(
        adata_sc,
        adata_sp,
        genes=markers,
    )
    ad_map = tg.map_cells_to_space(
        adata_sc,
        adata_sp,
        mode="constrained",
        target_count=adata_sp.obs.cell_count.sum(),
        density_prior="uniform",
        cluster_label="cell_type",
        device="cpu",
    )

    annotation = "raw_processed-V5 WT_filter70 MAP FULL"
    ad_map.write(
        os.path.join(DATADIR, "processed", "spatial", "combined", f"{annotation}.h5ad")
    )
