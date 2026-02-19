# python src/notebooks/pailin/SpaGRN.py &> ./outs/spatial/tools/spagrn.out

import os
import sys
from pathlib import Path

import argparse
import pandas as pd
import scanpy as sc
from multiprocessing import cpu_count

from spagrn.regulatory_network import InferNetwork as irn

# path to unfiltered loom file (this will be created in the optional steps below)
CORES = 10
DATADIR = Path("data") / "processed" / "external"
OUTSDIR = Path("outs")
REFDIR = Path("references")

annotation = "analyzed_SLIM.h5ad"
data = irn.read_file(DATADIR / annotation)
sc.pp.filter_genes(data, min_cells=10)
grn = irn(data)
grn.add_params('spg',{'prune_auc_threshold': 0.05, 'rank_threshold': 9000, 'auc_threshold': 0.05})

GRN_RANKINGS = " ".join(map(str, (REFDIR / "scenic").glob("*feather")))
MOTIF_ANNOTS = str(REFDIR / "scenic" / "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")
TF_LIST = str(REFDIR / "scenic" / "mm_mgi_tfs.txt")
LR_pairs = pd.read_csv(REFDIR / "omnipath_all.csv", index_col = 0)[['source', 'target']].rename(columns={'source' : 'from', 'target' : 'to'})

grn.infer(
    GRN_RANKINGS,
    MOTIF_ANNOTS,
    TF_LIST,
    niche_df=LR_pairs,
    target_genes=None,
    num_workers=CORES,
    cache=False,
    save_tmp=True,
    c_threshold=0.2,
    layers=None,
    latent_obsm_key='spatial',
    model='danb',
    n_neighbors=30,
    weighted_graph=False,
    cluster_label='celltype',
    method='spg',
    prefix='project',
    noweights=False)