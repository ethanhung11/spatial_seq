from typing import Iterable
from anndata import AnnData

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import decoupler as dc
import scanpy as sc
from glasbey import create_palette


def clear_uns(adata: AnnData, search: str):
    for i in pd.Series(adata.uns.keys()):
        if search in i:
            del adata.uns[i]


def empty_axs(axs: np.ndarray):
    for rs in axs[:]:
        for ax in rs[:]:
            ax.remove()
    return


def order_obs(adata: AnnData, col: str, order: Iterable[str]):
    adata.obs[col] = pd.Categorical(adata.obs[col], categories=order, ordered=True)
    return


def color_gen(groups: pd.Series | Iterable | np.array, custom_index=None):
    cs = create_palette(palette_size=len(groups.unique()))
    if custom_index is not None:
        return pd.Series(cs, index=custom_index)
    else:
        return pd.Series(cs, index=groups.unique())


def check_integration(
    adata: AnnData,
    category: str,
    f,
    embeddings: Iterable[str] = ["X_umap", "LocalMAP"],
    nrow: int = None,
    ncol: int = None,
    mini=False,
):
    print(f"Category {category} has {len(adata.obs[category].unique())} groups!")

    sf = f.subfigures(1, len(embeddings))
    int_colors = color_gen(adata.obs[category], adata.obs[category].cat.categories)

    for e, obsm in enumerate(embeddings):
        if mini is True:
            ax = sf[e].subplots(1, 1)
            sc.pl.embedding(
                adata,
                basis=obsm,
                color=category,
                ax=ax,
                show=False,
                palette=int_colors.to_list(),
            )

            ax.annotate(
                f"n = {adata.shape[0]}",
                size=10,
                fontweight="bold",
                xy=(0.98, 0.02),
                xycoords="axes fraction",
                horizontalalignment="right",
                verticalalignment="bottom",
            )
            # legend_loc='none' if e<len(embeddings)-1 else "right margin")

        else:
            assert nrow is not None and ncol is not None
            axs = sf[e].subplots(nrow * 2, ncol)
            gs = axs[0, 0].get_gridspec()
            empty_axs(axs)
            ax = sf[e].add_subplot(gs[:nrow, :])
            sc.pl.embedding(
                adata,
                basis=obsm,
                color=category,
                ax=ax,
                show=False,
                palette=int_colors.to_list(),
            )
            ax.annotate(
                f"n = {adata.shape[0]}",
                size=10,
                fontweight="bold",
                xy=(0.98, 0.02),
                xycoords="axes fraction",
                horizontalalignment="right",
                verticalalignment="bottom",
            )

            for n, group in enumerate(adata.obs[category].unique()):
                ax = sf[e].add_subplot(gs[nrow + n // nrow, n % nrow])
                sc.pl.embedding(
                    adata[adata.obs[category] == group],
                    basis=obsm,
                    color=category,
                    ax=ax,
                    show=False,
                    palette=[int_colors[group]],
                    legend_loc="none",
                )
                ax.set_title(group)

    return


def checkDoublets(
    adata,
    embedding="X_umap",
    cluster_key="leiden",
    doubletMethods=["scDblFinder", "DoubletFinder", "doubletdetection", "scrublet"],
):
    f = plt.figure(figsize=(20, 7), layout="constrained")
    sf = f.subfigures(1, 3, width_ratios=(1, 2, 1.5))

    axs = sf[0].subplots(2, 1)
    sc.pl.embedding(
        adata, basis=embedding, color=cluster_key, vmax=7000, ax=axs[0], show=False
    )
    sc.pl.embedding(
        adata, basis=embedding, color="n_genes", vmax=7000, ax=axs[1], show=False
    )

    axs = sf[1].subplots(2, 2)
    for n, method in enumerate(doubletMethods):
        adata.obs[f"predicted_doublet-{method}"] = adata.obs[
            f"predicted_doublet-{method}"
        ].astype(bool)
        sc.pl.embedding(
            adata,
            basis=embedding,
            color=f"predicted_doublet-{method}",
            ax=axs[n // 2, n % 2],
            show=False,
            alpha=0.7,
            palette=create_palette(2),
        )

    ax = sf[2].subplots(1, 1)
    doublet_df = pd.DataFrame(
        [
            adata.obs.groupby(cluster_key)[f"predicted_doublet-{method}"]
            .mean()
            .sort_values()
            for method in doubletMethods
        ]
    ).T
    doublet_df.columns = doublet_df.columns.str.replace(
        "predicted_doublet-", "", regex=True
    )
    sns.heatmap(doublet_df, annot=True, vmax=1, ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    return


def plot_violinplot(
    adata,
    markers,
    group: str,
    f,
    layer: str = "normalized",
    useStripPlot=True,
    bracket_params=None,
    ylabel_size=12,
):
    axs = f.subplots(len(markers), 1)
    for n, m in enumerate(markers):
        sc.pl.violin(
            adata,
            m,
            groupby=group,
            use_raw=False,
            layer=layer,
            show=False,
            ax=axs[n],
            stripplot=useStripPlot,
        )
        if n < len(markers) - 1:
            axs[n].set_xlabel("")
            axs[n].set_xticklabels([""] * len(axs[n].get_xticklabels()))
        axs[n].set_ylabel(axs[n].get_ylabel(), size=ylabel_size)

    if bracket_params is not None:
        ratios = bracket_params["ratio"] / np.sum(bracket_params["ratio"])
        ends = np.append(0, np.cumsum(ratios))
        bar_label_locs = [ends[i] + ratios[i] / 2 for i in range(len(ratios))]
        bar_bracket_widths = ratios * f.get_size_inches()[0] * 3.1

        axs = f.get_axes()
        for (
            n,
            label,
        ) in enumerate(bracket_params["labels"]):
            axs[-1].annotate(
                label,
                xy=(bar_label_locs[n], -bracket_params["bracket_y"]),
                xytext=(bar_label_locs[n], -bracket_params["label_y"]),
                xycoords="axes fraction",
                ha="center",
                va="bottom",
                bbox=dict(boxstyle="square", fc="none", color="none"),
                arrowprops=dict(
                    arrowstyle=f"-[, widthB={bar_bracket_widths[n]}, lengthB=0.3",
                    lw=1.0,
                    color="k",
                ),
            )
        axs[-1].set_xlabel(axs[-1].get_xlabel(), labelpad=bracket_params["padding"])

    return


def plot_cluster_violinplot(
    adata,
    group: str,
    clusters: str,
    markers,
    f,
):

    clusts = adata.obs[clusters].cat.categories
    cols = color_gen(adata.obs[group], adata.obs[group].cat.categories)
    sf = f.subfigures(2, len(clusts), height_ratios=[4, 1])

    crosstab_counts = pd.crosstab(adata.obs[clusters], adata.obs[group])
    crosstab_pct = crosstab_counts.div(crosstab_counts.sum(axis=0), axis=1) * 100

    for n, cluster in enumerate(clusts):
        cdata = adata[adata.obs[clusters] == cluster]
        plot_violinplot(cdata, markers, group, sf[0, n], useStripPlot=False)

        ax = sf[1, n].subplots(1, 1)
        crosstab_pct.loc[cluster][crosstab_pct.loc[cluster] > 0].plot(
            kind="bar", color=cols[crosstab_pct.loc[cluster] > 0], ax=ax
        )
        ax.tick_params(axis="x", rotation=0)
        ax.bar_label(
            ax.containers[0],
            labels=[
                f"{c:.2f}%\n({crosstab_counts.loc[cluster][n]})"
                for n, c in enumerate(crosstab_pct.loc[cluster])
            ],
        )

        sf[0, n].suptitle(f"Cluster {cluster}")

    return


def plot_cluster_stackedbarplot(
    adata,
    groupby: str,
    clusters: str,
    pct: bool = False,
    colors: list = None,
    ax=None,
):
    if colors is None:
        colors = color_gen(adata.obs[clusters])
    if ax is None:
        f, ax = plt.subplots(1, 2, figsize=(5, 5), layout="constrained")

    crosstab_counts = pd.crosstab(adata.obs[groupby], adata.obs[clusters])

    if pct is False:
        # Counts
        crosstab_counts.plot(kind="bar", stacked=True, ax=ax, color=colors)
        ax.set_title(
            f"{groupby} composition\n({clusters}, counts)",
            fontsize=14,
            fontweight="bold",
        )
        ax.set_xlabel(groupby, fontsize=12)
        ax.set_ylabel("Cell Count", fontsize=12)
        ax.legend(title=clusters, bbox_to_anchor=(1.05, 1), loc="upper left")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        ax.bar_label(ax.containers[-1], padding=2)
        ax.set_ylim(top=ax.get_ylim()[1] * 1.05)

    else:
        # Percentages
        crosstab_pct = crosstab_counts.div(crosstab_counts.sum(axis=1), axis=0) * 100

        crosstab_pct.plot(kind="bar", stacked=True, ax=ax, color=colors)
        ax.set_title(
            f"{groupby} composition\n({clusters}, percentage)",
            fontsize=14,
            fontweight="bold",
        )
        ax.set_xlabel(groupby, fontsize=12)
        ax.set_ylabel("Percentage (%)", fontsize=12)
        ax.legend(title=clusters, bbox_to_anchor=(1.05, 1), loc="upper left")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        ax.set_ylim(0, 100)

    return


def plot_cluster_silhouette(
    adata,
    obs_key="leiden",
    figsize=(10, 6),
    uns_key={"avg" : "silhouette_avg", "score" : "silhouette_scores"},
):
    cluster_labels = adata.obs[obs_key].astype("category").cat.codes
    n_clusters = len(np.unique(cluster_labels))
    
    f, ax = plt.subplots(figsize=figsize)
    y_lower = 10

    silhouette_avg = adata.uns[uns_key["avg"]][obs_key]

    for i in range(n_clusters):
        cluster_silhouette_values = adata.uns[uns_key["score"]][obs_key][cluster_labels == i]
        cluster_silhouette_values.sort()

        size_cluster_i = cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = mpl.cm.nipy_spectral(float(i) / n_clusters)
        ax.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )

        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
        y_lower = y_upper + 10

    ax.set_xlabel("Silhouette coefficient values")
    ax.set_ylabel("Cluster label")
    ax.set_title(f"Silhouette Plot (Average Score: {silhouette_avg:.3f})")

    ax.axvline(
        x=silhouette_avg,
        color="red",
        linestyle="--",
        label=f"Average Score: {silhouette_avg:.3f}",
    )
    ax.legend()
    plt.tight_layout()


def plot_c2c(
    adata: AnnData,
    key: str,
    pval: float = 0.05,
    top_n: int = 20,
    sources: Iterable[str] = None,
    targets: Iterable[str] = None,
    figsize=(15, 8),
):
    """Lightweight cell-cell communication dotplot with controlled dot sizes."""

    df = adata.uns[key]
    df = df[df["magnitude_rank"] < pval]

    if sources:
        df = df[df["source"].isin(sources)]
    if targets:
        df = df[df["target"].isin(targets)]

    df = df.nsmallest(top_n, "magnitude_rank")

    # Use complex columns if available, fallback to simple
    lig_col = "ligand_complex" if "ligand_complex" in df.columns else "ligand"
    rec_col = "receptor_complex" if "receptor_complex" in df.columns else "receptor"
    df["lr"] = df[lig_col] + " -> " + df[rec_col]

    # Efficient size mapping: 5 quantile-based bins, min size 10
    mag_vals = -np.log10(np.maximum(df["magnitude_rank"], 1e-10))
    size_bins = np.array([10, 50, 100, 150, 200])
    # Use 1 as minimum, then map higher values to larger sizes
    plt.Normalize(vmin=1, vmax=mag_vals.max())
    df["dot_size"] = size_bins[
        np.searchsorted(
            np.percentile(mag_vals, [20, 40, 60, 80]), np.maximum(mag_vals, 1)
        )
    ]

    src_list = sorted(df["source"].unique())
    tgt_list = sorted(df["target"].unique())
    lr_list = df["lr"].drop_duplicates().tolist()

    fig, axes = plt.subplots(
        1, len(src_list), figsize=figsize, sharey=True, layout="constrained"
    )
    if len(src_list) == 1:
        axes = [axes]

    # Pre-compute color normalization
    spec_vals = -np.log10(np.maximum(df["specificity_rank"], 1e-10))
    norm = plt.Normalize(vmin=1, vmax=spec_vals.max())

    for i, src in enumerate(src_list):
        ax = axes[i]
        data = df[df["source"] == src]

        if len(data):
            x = [tgt_list.index(t) for t in data["target"]]
            y = [lr_list.index(lr) for lr in data["lr"]]
            colors = plt.cm.Greens(
                norm(-np.log10(np.maximum(data["specificity_rank"], 1e-10)))
            )
            ax.scatter(
                x, y, c=colors, s=data["dot_size"], alpha=0.8, edgecolors=colors, lw=1
            )

        ax.set(
            xlim=(-0.5, len(tgt_list) - 0.5),
            ylim=(-0.5, len(lr_list) - 0.5),
            xticks=range(len(tgt_list)),
            title=src,
        )
        ax.set_xticklabels(tgt_list, rotation=45, ha="right")
        ax.grid(alpha=0.3)

    axes[0].set(yticks=range(len(lr_list)), ylabel="Ligand -> Receptor")
    axes[0].set_yticklabels(lr_list)
    fig.suptitle("Source")
    fig.text(0.5, 0.02, "Target", ha="center")

    fig.subplots_adjust(right=0.75)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap="Greens", norm=norm)
    plt.colorbar(sm, ax=axes, label="Specificity (-log10 p-val)")

    # Create proper size legend using evenly spaced integer values
    if len(mag_vals) > 0:  # Only create legend if we have data
        # Create evenly spaced integer thresholds from 1 to max
        max_val = int(np.ceil(mag_vals.max()))
        if max_val <= 1:
            size_labels = ["1"]
        else:
            # Create 5 evenly spaced integers from 1 to max
            thresholds = np.linspace(1, max_val, 5).astype(int)
            # Remove duplicates while preserving order
            thresholds = np.array(sorted(set(thresholds)))
            size_labels = [str(t) for t in thresholds]

        # Create legend elements for matplotlib
        from matplotlib.lines import Line2D

        legend_elements = []
        # Use only the number of sizes we actually have
        for i, label in enumerate(size_labels):
            if i < len(size_bins):
                size = size_bins[i]
                legend_elements.append(
                    Line2D(
                        [0],
                        [0],
                        marker="o",
                        color="w",
                        markerfacecolor="gray",
                        markersize=np.sqrt(size / 10),
                        alpha=0.7,
                        markeredgecolor="k",
                        label=label,
                    )
                )

        # Add legend to the figure
        fig.legend(
            handles=legend_elements,
            title="Magnitude\n(-log10 p-val)",
            loc="upper left",
            frameon=True,
            bbox_to_anchor=(1, 0.95),
            fancybox=True,
            shadow=True,
        )

    return


def plot_c2c(
    adata: AnnData,
    key: str,
    pval: float = 0.05,
    top_n: int = 20,
    sources: Iterable[str] = None,
    targets: Iterable[str] = None,
    figsize=(15, 8),
):
    """Lightweight cell-cell communication dotplot with controlled dot sizes."""

    df = adata.uns[key]
    df = df[df["magnitude_rank"] < pval]

    if sources:
        df = df[df["source"].isin(sources)]
    if targets:
        df = df[df["target"].isin(targets)]

    df = df.nsmallest(top_n, "magnitude_rank")

    # Use complex columns if available, fallback to simple
    lig_col = "ligand_complex" if "ligand_complex" in df.columns else "ligand"
    rec_col = "receptor_complex" if "receptor_complex" in df.columns else "receptor"
    df["lr"] = df[lig_col] + " -> " + df[rec_col]

    # Efficient size mapping: 5 quantile-based bins, min size 10
    mag_vals = -np.log10(np.maximum(df["magnitude_rank"], 1e-10))
    size_bins = np.array([10, 50, 100, 150, 200])
    # Use 1 as minimum, then map higher values to larger sizes
    plt.Normalize(vmin=1, vmax=mag_vals.max())
    df["dot_size"] = size_bins[
        np.searchsorted(
            np.percentile(mag_vals, [20, 40, 60, 80]), np.maximum(mag_vals, 1)
        )
    ]

    src_list = sorted(df["source"].unique())
    tgt_list = sorted(df["target"].unique())
    lr_list = df["lr"].drop_duplicates().tolist()

    fig, axes = plt.subplots(
        1, len(src_list), figsize=figsize, sharey=True, layout="constrained"
    )
    if len(src_list) == 1:
        axes = [axes]

    # Pre-compute color normalization
    spec_vals = -np.log10(np.maximum(df["specificity_rank"], 1e-10))
    norm = plt.Normalize(vmin=1, vmax=spec_vals.max())

    for i, src in enumerate(src_list):
        ax = axes[i]
        data = df[df["source"] == src]

        if len(data):
            x = [tgt_list.index(t) for t in data["target"]]
            y = [lr_list.index(lr) for lr in data["lr"]]
            colors = plt.cm.Greens(
                norm(-np.log10(np.maximum(data["specificity_rank"], 1e-10)))
            )
            ax.scatter(
                x, y, c=colors, s=data["dot_size"], alpha=0.8, edgecolors=colors, lw=1
            )

        ax.set(
            xlim=(-0.5, len(tgt_list) - 0.5),
            ylim=(-0.5, len(lr_list) - 0.5),
            xticks=range(len(tgt_list)),
            title=src,
        )
        ax.set_xticklabels(tgt_list, rotation=45, ha="right")
        ax.grid(alpha=0.3)

    axes[0].set(yticks=range(len(lr_list)), ylabel="Ligand -> Receptor")
    axes[0].set_yticklabels(lr_list)
    fig.suptitle("Source")
    fig.text(0.5, 0.02, "Target", ha="center")

    fig.subplots_adjust(right=0.75)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap="Greens", norm=norm)
    plt.colorbar(sm, ax=axes, label="Specificity (-log10 p-val)")

    # Create proper size legend using evenly spaced integer values
    if len(mag_vals) > 0:  # Only create legend if we have data
        # Create evenly spaced integer thresholds from 1 to max
        max_val = int(np.ceil(mag_vals.max()))
        if max_val <= 1:
            size_labels = ["1"]
        else:
            # Create 5 evenly spaced integers from 1 to max
            thresholds = np.linspace(1, max_val, 5).astype(int)
            # Remove duplicates while preserving order
            thresholds = np.array(sorted(set(thresholds)))
            size_labels = [str(t) for t in thresholds]

        # Create legend elements for matplotlib
        from matplotlib.lines import Line2D

        legend_elements = []
        # Use only the number of sizes we actually have
        for i, label in enumerate(size_labels):
            if i < len(size_bins):
                size = size_bins[i]
                legend_elements.append(
                    Line2D(
                        [0],
                        [0],
                        marker="o",
                        color="w",
                        markerfacecolor="gray",
                        markersize=np.sqrt(size / 10),
                        alpha=0.7,
                        markeredgecolor="k",
                        label=label,
                    )
                )

        # Add legend to the figure
        fig.legend(
            handles=legend_elements,
            title="Magnitude\n(-log10 p-val)",
            loc="upper left",
            frameon=True,
            bbox_to_anchor=(1, 0.95),
            fancybox=True,
            shadow=True,
        )

    return


def plot_gsea_dc(
    adata, name, key="score_ulm", group="cell_type", n_markers=5, flip=True, f=None
):
    ax = f.get_axes() if f is not None else None

    score = dc.pp.get_obsm(adata=adata, key=key)
    df = dc.tl.rankby_group(
        adata=score, groupby=group, reference="rest", method="wilcoxon"
    )
    df = df[df["stat"] > 0]

    source_markers = (
        df.groupby("group")
        .head(n_markers)
        .drop_duplicates("name")
        .groupby("group")["name"]
        .apply(lambda x: list(x))
        .to_dict()
    )

    sc.pl.matrixplot(
        adata=score,
        var_names=source_markers,
        groupby=group,
        dendrogram=False,
        standard_scale="var",
        colorbar_title="Z-scaled scores",
        cmap="Reds",
        swap_axes=flip,
        title=name,
        ax=ax,
    )


def plot_go_enrichment(
    df_dict,
    pvalue_col,
    score_col,
    names_to_plot=None,
    pvalue_threshold=None,
    score_threshold=None,
    rank_by=None,
    top_n=None,
    use_log_pvalue=False,
    use_log_score=False,
    database=None,
    figsize=(8, 10),
    **kwargs,
):
    """
    Create a vertical dot plot from multiple dataframes with discrete size intervals.

    Parameters:
    -----------
    df_dict : dict
        Dictionary where keys are dataframe names and values are pandas DataFrames
    pvalue_col : str
        Column name containing p-values (used for color)
    score_col : str
        Column name containing scores (used for size intervals)
    names_to_plot : dict, optional
        Dictionary mapping dataframe keys to lists of names to plot
    pvalue_threshold : float, optional
        Only include names with p-values below this threshold
    score_threshold : float, optional
        Only include names with scores above this threshold
    rank_by : str, optional
        Column to rank by ('pvalue' or 'score'). Use with top_n
    top_n : int, optional
        Number of top-ranked names to include per dataframe
    use_log_pvalue : bool, default False
        Whether to apply -log10 transformation to p-values
    use_log_score : bool, default False
        Whether to apply -log10 transformation to scores
    figsize : tuple, default (8, 10)
        Figure size
    **kwargs : dict
        Additional arguments passed to plt.scatter

    Returns:
    --------
    fig, ax : matplotlib figure and axes objects
    """

    # Set default scatter parameters
    scatter_params = {
        "alpha": 0.7,
        "cmap": "viridis",
        "edgecolors": "black",
        "linewidth": 0.5,
    }
    scatter_params.update(kwargs)

    # Filter names based on criteria
    if names_to_plot is None:
        names_to_plot = {}
        for df_name, df in df_dict.items():
            df_filtered = df.copy()

            if pvalue_threshold is not None:
                df_filtered = df_filtered[df_filtered[pvalue_col] <= pvalue_threshold]
            if score_threshold is not None:
                df_filtered = df_filtered[df_filtered[score_col] >= score_threshold]

            if rank_by is not None and top_n is not None:
                if rank_by == "pvalue":
                    df_filtered = df_filtered.nsmallest(top_n, pvalue_col)
                elif rank_by == "score":
                    df_filtered = df_filtered.nlargest(top_n, score_col)

            names_to_plot[df_name] = df_filtered["name"].tolist()

    # Get all unique names and sort by best p-value across all dataframes
    all_names = set()
    for df_name, df in df_dict.items():
        if df_name in names_to_plot:
            all_names.update(names_to_plot[df_name])
        else:
            all_names.update(df["name"].values)

    # Create a mapping of names to their best (lowest) p-value across all dataframes
    name_to_best_pvalue = {}
    for name in all_names:
        best_pvalue = float("inf")
        for df_name, df in df_dict.items():
            if name in df["name"].values:
                pvalue = df[df["name"] == name][pvalue_col].iloc[0]
                if np.isfinite(pvalue) and pvalue > 0:
                    best_pvalue = min(best_pvalue, pvalue)
        name_to_best_pvalue[name] = best_pvalue

    # Sort names by best p-value (lowest first)
    all_names = sorted(all_names, key=lambda x: name_to_best_pvalue[x])
    name_to_y = {name: i for i, name in enumerate(all_names)}

    # Collect all valid data for normalization
    all_pvalues, all_scores = [], []
    for df_name, df in df_dict.items():
        df_filtered = df[df["name"].isin(names_to_plot.get(df_name, df["name"]))]

        pvalues = df_filtered[pvalue_col].values
        scores = df_filtered[score_col].values

        valid_mask = (
            np.isfinite(pvalues) & np.isfinite(scores) & (pvalues > 0) & (scores > 0)
        )
        all_pvalues.extend(pvalues[valid_mask])
        all_scores.extend(scores[valid_mask])

    # Transform data and create normalizations
    pvalues_transformed = -np.log10(all_pvalues) if use_log_pvalue else all_pvalues
    scores_transformed = -np.log10(all_scores) if use_log_score else all_scores

    pvalue_norm = mpl.colors.Normalize(
        vmin=np.min(pvalues_transformed), vmax=np.max(pvalues_transformed)
    )
    score_bins = np.percentile(scores_transformed, [0, 20, 40, 60, 80, 100])
    score_sizes = [40, 100, 160, 220, 280]

    # Create plot
    fig, ax = plt.subplots(figsize=figsize, layout="constrained")

    # Plot data for each dataframe (sort dataframes alphabetically)
    df_names_sorted = sorted(df_dict.keys())
    for col_idx, df_name in enumerate(df_names_sorted):
        df = df_dict[df_name]
        df_filtered = df[df["name"].isin(names_to_plot.get(df_name, df["name"]))]

        for _, row in df_filtered.iterrows():
            name, pvalue, score = row["name"], row[pvalue_col], row[score_col]

            if not (
                np.isfinite(pvalue) and np.isfinite(score) and pvalue > 0 and score > 0
            ):
                continue

            pvalue_plot = -np.log10(pvalue) if use_log_pvalue else pvalue
            score_plot = -np.log10(score) if use_log_score else score

            # Get size based on percentile intervals
            size_idx = np.searchsorted(score_bins[1:], score_plot)
            size_to_use = score_sizes[min(size_idx, len(score_sizes) - 1)]

            ax.scatter(
                col_idx,
                name_to_y[name],
                c=pvalue_plot,
                s=size_to_use,
                norm=pvalue_norm,
                **scatter_params,
            )

    # Setup axes
    ax.set_xlim(-0.5, len(df_dict) - 0.5)
    ax.set_ylim(-0.5, len(all_names) - 0.5)
    ax.set_xticks(range(len(df_dict)))
    ax.set_xticklabels(df_names_sorted, rotation=45, ha="right")
    ax.set_yticks(range(len(all_names)))
    ax.set_yticklabels(all_names)
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("Dataframe")
    ax.set_ylabel(f"GO Enrichment Term ({database})")

    # Add legends
    # Color legend
    dummy_scatter = ax.scatter(
        [], [], c=[], norm=pvalue_norm, cmap=scatter_params["cmap"]
    )
    cbar = plt.colorbar(dummy_scatter, ax=ax, shrink=0.8)
    cbar.set_label(
        f"-log10({pvalue_col})" if use_log_pvalue else pvalue_col,
        rotation=270,
        labelpad=20,
    )

    # Size legend with actual score values
    legend_elements = [
        plt.scatter([], [], s=s, c="gray", alpha=0.6, edgecolors="black", linewidth=0.5)
        for s in score_sizes
    ]

    # Create labels with actual score ranges
    labels = []
    for i in range(len(score_bins) - 1):
        min_val = score_bins[i]
        max_val = score_bins[i + 1]
        if use_log_score:
            min_val = 10 ** (-min_val) if min_val != 0 else 0
            max_val = 10 ** (-max_val) if max_val != 0 else 0
        labels.append(f"{min_val:.1f} - {max_val:.1f}")

    ax.legend(
        legend_elements,
        labels,
        scatterpoints=1,
        loc="upper left",
        bbox_to_anchor=(1.3, 1),
        title=f"{score_col} ranges",
    )

    # Add title
    ax.set_title(f"GO Enrichment via {database}")

    return fig
