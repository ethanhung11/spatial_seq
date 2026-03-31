from typing import Iterable, Literal
import os
import numpy as np
from tqdm import tqdm

import scipy
import scanpy as sc
import pyucell as uc
import liana as li
import decoupler as dc
import gseapy as gp
import pandas as pd
from .ccc import liana_mouse_resource
from utils import clear_obs

def module_score(
    adata: sc.AnnData,
    genesets: dict,
    collection_name: str,
    method: str | Literal["Seurat", "UCell", "Decoupler"] = "Seurat",
    decoupler_method: str | Literal["ULM", "FGSEA", "GSVA", "AUCell", "ORA", "Consensus"] = None,
    zscore: bool = False,
    save_obs: bool = False,
):
    for name in genesets:
        genesets[name] = pd.Series(genesets[name]).unique()

    if method == "Seurat":
        for name in tqdm(genesets):
            sc.tl.score_genes(adata, genesets[name], score_name=f"{name}-Seurat")
        
        pat = f"-{method}|".join(list(genesets.keys())) + f"-{method}"
        adata.obsm[f"{collection_name}-{method}"] = adata.obs.loc[:,adata.obs.columns.str.contains(pat, regex=True)]
        clear_obs(adata, pat)

    elif method == "UCell":
        uc.compute_ucell_scores(adata, signatures=genesets, suffix="-UCell")
        pat = f"-{method}|".join(list(genesets.keys())) + f"-{method}"
        adata.obsm[f"{collection_name}-{method}"] = adata.obs.loc[:,adata.obs.columns.str.contains(pat, regex=True)]
        clear_obs(adata, pat)

    elif method == "Decoupler":
        assert decoupler_method is not None
        method = decoupler_method

        genesets = pd.concat(
            [
                pd.DataFrame({key: genesets[key]}).melt(
                    var_name="source", value_name="target"
                )
                for key in genesets
            ]
        )
        genesets = genesets.drop_duplicates().reset_index(drop=True)
        genesets['source'] = genesets['source'] + f"-{method}"

        if method == "ULM":
            dc.mt.ulm(data=adata, net=genesets, layer="normalized", verbose=True)
            adata.obsm[f"{collection_name}-{method}"] = adata.obsm["score_ulm"]
            adata.obsm[f"{collection_name}-{method}_padj"] = adata.obsm["padj_ulm"]
            del adata.obsm["score_ulm"], adata.obsm["padj_ulm"]
        elif method == "FGSEA":
            dc.mt.gsea(data=adata, net=genesets, layer="normalized", verbose=True)
            adata.obsm[f"{collection_name}-{method}"] = adata.obsm["score_gsea"]
            adata.obsm[f"{collection_name}-{method}_padj"] = adata.obsm["padj_gsea"]
            del adata.obsm["score_gsea"], adata.obsm["padj_gsea"]
        elif method == "GSVA":
            dc.mt.gsva(data=adata, net=genesets, layer="normalized", verbose=True)
            adata.obsm[f"{collection_name}-{method}"] = adata.obsm["score_gsva"]
            del adata.obsm["score_gsva"]
        elif method == "AUCell":
            dc.mt.aucell(data=adata, net=genesets, layer="normalized", verbose=True)
            adata.obsm[f"{collection_name}-{method}"] = adata.obsm["score_aucell"]
            del adata.obsm["score_aucell"]
        elif method == "ORA":
            dc.mt.ora(data=adata, net=genesets, layer="normalized", verbose=True)
            adata.obsm[f"{collection_name}-{method}"] = adata.obsm["score_ora"]
            del adata.obsm["score_ora"]
        elif method == "Consensus":
            dc.mt.consensus(data=adata, net=genesets, layer="normalized", verbose=True)
            # adata.obsm[f"{collection_name}-{method}"] = adata.obsm["score_consensus"]
            # del adata.obsm["score_consensus"]
        else:
            raise ValueError("Not a valid method!")          

    else:
        raise ValueError("Not a valid method!")

    if zscore is True:
        for name in genesets:
            col = f"{name}_{method}"
            adata.obs[col + "_z"] = scipy.stats.zscore(adata.obs[col])

    if save_obs is True:
        adata.obs = adata.obs.join(adata.obsm[f"{collection_name}-{method}"])
        # pat = f"-{method}|".join(list(genesets.keys())) + f"-{method}"
        # clear_obs(adata, pat)
        del adata.obsm[f"{collection_name}-{method}"]

    return adata


def cell2cell_interactions(
    adata: sc.AnnData,
    cell_group: str,
    resource: pd.DataFrame,
    key="ccc",
    methods: Iterable = None,
    filter_results: bool = True,
    cores: int = None,
):
    if cores is None:
        cores = os.cpu_count()

    # pull interaction databases
    if resource is None:
        resource = liana_mouse_resource()

    # run all methods
    if methods is None:
        rank_aggregate_custom = li.method.rank_aggregate
    else:
        rank_aggregate_custom = li.mt.AggregateClass(
            li.mt.aggregate_meta, methods=methods
        )

    rank_aggregate_custom(
        adata,
        groupby=cell_group,
        layer="normalized",
        use_raw=False,
        key_added=key,
        de_method="wilcoxon",
        return_all_lrs=True,
        verbose=True,
        n_jobs=cores,
        resource=resource[["ligand", "receptor"]],
    )

    # save reference database info
    adata.uns[key] = adata.uns[key].merge(
        resource[["ligand", "receptor", "db_sources"]],
        left_on=["ligand_complex", "receptor_complex"],
        right_on=["ligand", "receptor"],
        how="left",
    )

    if filter_results is True:
        # filter for quality interactions
        ccc_filters = []

        # Cell Specificity filters
        if "cellphone_pvals" in adata.uns[key].columns:  # CellphoneDB
            ccc_filters.append(adata.uns[key]["cellphone_pvals"] <= 0.05)
        if "gmean_pvals" in adata.uns[key].columns:  # CellphoneDB V2
            ccc_filters.append(adata.uns[key]["gmean_pvals"] <= 0.05)
        if "cellchat_pvals" in adata.uns[key].columns:  # CellChat
            ccc_filters.append(adata.uns[key]["cellchat_pvals"] <= 0.05)
        if "lr_logfc" in adata.uns[key].columns:  # log2FC
            ccc_filters.append(adata.uns[key]["lr_logfc"] > 0)
        if "scaled_weight" in adata.uns[key].columns:  # Connectome
            ccc_filters.append(
                adata.uns[key]["scaled_weight"]
                > adata.uns[key]["scaled_weight"].quantile(0.95)
            )
        if "spec_weight" in adata.uns[key].columns:  # NATMI
            ccc_filters.append(
                adata.uns[key]["spec_weight"]
                > adata.uns[key]["spec_weight"].quantile(0.95)
            )
        if "specificity_rank" in adata.uns[key].columns:  # Liana (aggregated score)
            ccc_filters.append(adata.uns[key]["specificity_rank"] <= 0.05)

        # Magnitue filters
        if "lr_probs" in adata.uns[key].columns:  # CellChat
            ccc_filters.append(adata.uns[key]["lr_probs"] <= 0.05)
        if "lrscore" in adata.uns[key].columns:  # SingleCellSignalR
            ccc_filters.append(adata.uns[key]["lrscore"] > 0.6)
        if "expr_prod" in adata.uns[key].columns:  # NATMI/Connectome
            ccc_filters.append(
                adata.uns[key]["expr_prod"] > adata.uns[key]["expr_prod"].quantile(0.95)
            )
        if "magnitude_rank" in adata.uns[key].columns:  # Liana (aggregated score)
            ccc_filters.append(adata.uns[key]["magnitude_rank"] <= 0.05)

        df = pd.concat(ccc_filters, axis=1)
        ccc_filter_all = df.all(axis=1)
        adata.uns[key + "_filtered"] = adata.uns[key][ccc_filter_all].reset_index(
            drop=True
        )

    return adata


def GO_Enrich(
    adata: sc.AnnData,
    groupby: str,
    key: str,
    sources: Literal["GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC"] = [
        "GO:MF",
        "GO:CC",
        "GO:BP",
        "KEGG",
        "REAC",
    ],
    pval_cutoff: float = 1e-4,
    log2fc_min: float = 2,
):
    # see here for other gprofilier args: https://biit.cs.ut.ee/gprofiler/page/apis

    GO_enrichments = {}

    for category in adata.obs[groupby].unique():
        df = sc.get.rank_genes_groups_df(
            adata,
            group=category,
            key=key,
            pval_cutoff=pval_cutoff,
            log2fc_min=log2fc_min,
        )
        for src in sources:
            if src not in GO_enrichments:
                GO_enrichments[src] = {}
            GO_enrichments[src][category] = sc.queries.enrich(
                df.names.to_list(),
                org="mmusculus",
                gprofiler_kwargs={"sources": [src], "all_results": True},
            )

    return GO_enrichments


def gseapy_gsea(
    adata: sc.AnnData,
    groupby: str,
    genesets: Iterable,
    method: Literal["ssgea", "gsea", "prerank", "gsva"] = "gsea",
    threads: int = 10,
    seed: int = 6,
    one_condition = None,
    **kwargs,
):
    if len(adata.obs[groupby].unique()) < 2:
        return ValueError(f"Group '{groupby}' does not have at least 2 groups!")

    elif one_condition is not None:
        first_idx = np.where(adata.obs[groupby] == one_condition)[0][0]
        order = np.arange(np.shape(adata)[0])
        order[first_idx], order[0] = order[0], order[first_idx]
        tmp = adata[order]

        out = gseapy_helper(
            tmp,
            genesets,
            method,
            groupby,
            group,
            threads,
            seed,
            **kwargs,
        )
        result = out.res2d
        names = result["Term"].str.split("__", expand=True)
        result["Collection"] = names[0]
        result["Term"] = names[1]
        gsea_results = result

    else:
        gsea_results = {}
        for group in tqdm(adata.obs[groupby].unique()):
            # reorder group to be first
            first_idx = np.where(adata.obs[groupby] == group)[0][0]
            order = np.arange(np.shape(adata)[0])
            order[first_idx], order[0] = order[0], order[first_idx]
            tmp = adata[order]

            out = gseapy_helper(
                tmp,
                genesets,
                method,
                groupby,
                group,
                threads,
                seed,
                **kwargs,
            )
            result = out.res2d
            names = result["Term"].str.split("__", expand=True)
            result["Collection"] = names[0]
            result["Term"] = names[1]
            gsea_results[group] = result

    return gsea_results


def gseapy_helper(
    adata=None,
    genesets=None,
    method=None,
    groupby=None,
    group=None,
    threads=None,
    seed=None,
    **kwargs,
):
    if method == "prerank":
        df = sc.get.rank_genes_groups_df(adata, group=group, **kwargs)
        ranking = df[["names", "scores"]]
        out = gp.prerank(
            rnk=ranking,
            gene_sets=genesets,
            seed=seed,
            threads=threads,
            permutation_num=100,
        )
    elif method == "gsea":
        class_vector = adata.obs[groupby] == group
        assert class_vector[0] == True
        out = gp.gsea(
            data=adata.to_df().T,
            gene_sets=genesets,
            cls=class_vector.to_list(),
            threads=threads,
            seed=seed,
            permutation_num=100,
        )
    elif method == "ssgea":
        out = gp.ssgsea(
            data=adata[adata.obs[groupby] == group].to_df().T,
            gene_sets=genesets,
            threads=threads,
            seed=seed,
            permutation_num=100,
        )
    elif method == "gsva":
        out = gp.gsva(
            data=adata[adata.obs[groupby] == group].to_df().T,
            gene_sets=genesets,
            threads=threads,
            seed=seed,
            permutation_num=100,
        )

    return out
