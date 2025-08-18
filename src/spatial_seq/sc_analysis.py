from typing import Iterable, Literal
from anndata import AnnData

import os
import scanpy as sc
import liana as li
import decoupler as dc
import gseapy as gp
import pandas as pd
from .sc_ccc import liana_mouse_resource


def cell2cell_interactions(
    adata: AnnData,
    cell_group="cell_type",
    key="ccc",
    resource_opts: Iterable[str] = [
        "baccin2019",
        "cellcall",
        "cellchatdb",
        "cellinker",
        "cellphonedbv5",
        "cellphonedb",
        "celltalkdb",
        "connectomedb2020",
        "consensus",
        "embrace",
        "guide2pharma",
        "hpmr",
        "icellnet",
        "italk",
        "kirouac2010",
        "lrdb",
        "ramilowski2015",
        "mouseconsensus",
    ],
    methods: Iterable = None,
    cores: int = None,
):
    if cores is None:
        cores = os.cpu_count()

    # pull interaction databases
    ccc_db = liana_mouse_resource(resource_opts)

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
        resource=ccc_db[["ligand", "receptor"]],
    )

    # save reference database info
    adata.uns[key] = adata.uns[key].merge(
        ccc_db[["ligand", "receptor", "db_sources"]],
        left_on=["ligand_complex", "receptor_complex"],
        right_on=["ligand", "receptor"],
        how="left",
    )

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
            adata.uns[key]["spec_weight"] > adata.uns[key]["spec_weight"].quantile(0.95)
        )

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
    adata.uns[key + "_filtered"] = adata.uns[key][ccc_filter_all].reset_index(drop=True)

    return adata


def GO_Enrich(
    adata,
    groupby,
    key,
    sources=["GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC"],
    pval_cutoff=1e-4,
    log2fc_min=2,
):
    # see here for other gprofilier args: https://biit.cs.ut.ee/gprofiler/page/apis

    GO_enrichments = {}
    for src in sources:
        GO_enrichments[src] = {}

    for category in adata.obs[groupby].unique():
        df = sc.get.rank_genes_groups_df(
            adata,
            group=category,
            key=key,
            pval_cutoff=pval_cutoff,
            log2fc_min=log2fc_min,
        )
        for src in sources:
            GO_enrichments[src][category] = sc.queries.enrich(
                df.names.to_list(),
                org="mmusculus",
                gprofiler_kwargs={"sources": [src], "all_results": True},
            )

    return GO_enrichments


def GSEA_decoupler(
    adata: AnnData,
    name: str,
    type: Literal["ULM", "GSEA", "GSVA", "AUCELL"] = "ULM",
    geneset_dir: str = None,
    geneset=None,
    remove_prefix=False,
):
    if geneset is None:
        assert geneset_dir is not None
        geneset = dc.pp.read_gmt(geneset_dir)
        if remove_prefix is True:
            prefixes = geneset.source.str.split("_", expand=True)[0].unique().tolist()
            for pre in prefixes:
                geneset.source = geneset.source.str.replace(f"{pre}_", "")

    if type == "ULM":
        dc.mt.ulm(data=adata, net=geneset, layer="normalized", verbose=True)
        adata.obsm[f"{name}_score_ulm"] = adata.obsm["score_ulm"]
        adata.obsm[f"{name}_padj_ulm"] = adata.obsm["padj_ulm"]
        del adata.obsm["score_ulm"], adata.obsm["padj_ulm"]

    elif type == "GSEA":
        dc.mt.gsea(data=adata, net=geneset, layer="normalized", verbose=True)
        adata.obsm[f"{name}_score_gsea"] = adata.obsm["score_gsea"]
        adata.obsm[f"{name}_padj_gsea"] = adata.obsm["padj_gsea"]
        del adata.obsm["score_gsea"], adata.obsm["padj_gsea"]

    elif type == "GSVA":
        dc.mt.gsea(data=adata, net=geneset, layer="normalized", verbose=True)
        adata.obsm[f"{name}_score_gsea"] = adata.obsm["score_gsva"]
        del adata.obsm["score_gsva"]

    elif type == "AUCell":
        dc.mt.aucell(data=adata, net=geneset, layer="normalized", verbose=True)
        adata.obsm[f"{name}_score_aucell"] = adata.obsm["score_aucell"]
        del adata.obsm["score_aucell"]

    else:
        raise ValueError("Not a valid method!")

    return adata


def GSEA_gseapy(
    adata: AnnData,
    genesets: list,
    method: Literal["ssgea", "gsea", "prerank", "gsva"] = "prerank",
    groupby: str = None,
    threads: int = 10,
    seed: int = 6,
    **kwargs,
):
    if len(adata.obs[groupby].unique()) < 2:
        return ValueError(f"Group '{groupby}' does not have at least 2 groups!")

    elif len(adata.obs[groupby].unique()) == 2:

        out = GSEA_gseapy_helper(
            adata,
            genesets,
            method,
            groupby,
            adata.obs[groupby].unique()[1],
            threads,
            seed,
            **kwargs,
        )
        gsea_results = out.res2d
        names = gsea_results["Term"].str.split("__", expand=True)
        gsea_results["Collection"] = names[0]
        gsea_results["Term"] = names[1]

    else:
        gsea_results = {}
        for group in adata.obs[groupby].unique():
            out = GSEA_gseapy_helper(
                adata,
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


def GSEA_gseapy_helper(
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
        out = gp.gsea(
            data=adata.to_df().T,
            gene_sets=genesets,
            cls=adata.obs[groupby],
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
