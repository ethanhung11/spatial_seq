from typing import Iterable

from collections import Counter
from tqdm import tqdm
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import plotly.graph_objects as go
from glasbey import create_palette
from pyclustree import clustree

import scanpy as sc
import decoupler as dc


def empty_axs(axs: np.ndarray):
    for rs in axs[:]:
        for ax in rs[:]:
            ax.remove()
    return


def order_obs(adata: sc.AnnData, col: str, order: Iterable[str]):
    adata.obs[col] = pd.Categorical(adata.obs[col], categories=order, ordered=True)
    return


def color_gen(groups: pd.Series | Iterable | np.array, name_matched=False):
    n_colors = len(groups)
    if n_colors > 200:
        raise ValueError(
            "Cannot have more than 200 colors!\nCheck that you are submitting pd.Series.cat.categories, not the whole pd.Series"
        )
    colors = create_palette(palette_size=n_colors)
    if name_matched is True:
        return pd.Series(colors, index=groups)
    else:
        return colors


def color_gen_seq(
    groups: pd.Series | Iterable | np.array, name_matched=False, cmap_name="RdBu"
):
    n_colors = len(groups)
    if n_colors > n_colors:
        raise ValueError(
            "Cannot have more than 50 colors!\nCheck that you are submitting pd.Series.cat.categories, not the whole pd.Series"
        )
    cmap = plt.get_cmap(cmap_name)
    colors = [mpl.colors.to_hex(color) for color in cmap(np.linspace(0, 1, n_colors))]

    if name_matched is True:
        return pd.Series(colors, index=groups)
    else:
        return colors


def check_QCPlot(df, value, groupby):
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    g = sns.FacetGrid(
        df, row=groupby, hue=groupby, aspect=15, height=0.5, palette="tab20"
    )

    g.map(sns.kdeplot, value, clip_on=False, fill=True, alpha=1, linewidth=1.5)
    g.map(sns.kdeplot, value, clip_on=False, color="w", lw=2)

    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    def label(x, color, label):
        ax = plt.gca()
        ax.text(
            0,
            0.2,
            label,
            fontweight="bold",
            color=color,
            ha="left",
            va="center",
            transform=ax.transAxes,
        )

    g.map(label, value)

    g.figure.subplots_adjust(hspace=-0.6)

    g.set_titles("")
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)

    for ax in g.axes.flat:
        ax.axvline(x=df[value].median(), color="r", linestyle="-")

    return g.figure


def check_integration(
    adata: sc.AnnData,
    category: str,
    f = None,
    embeddings: Iterable[str] = ["X_umap", "LocalMAP"],
    nrow: int = None,
    palette: pd.Series = None,
    mini=False,
):
    if f is None:
        f = plt.figure(figsize=(10,10))
    
    if adata.obs[category].dtype != "category":
        print(f"Converting adata.obs[{category}] to `category` dtype.")
        adata.obs[category] = adata.obs[category].astype("category")
    print(f"found {len(adata.obs[category].cat.categories)} groups in {category}")
    if palette is None:
        if adata.uns.get(f"{category}_colors") is None:
            print("Palette not found! Generating palette.")
            palette = color_gen(adata.obs[category].cat.categories, True)
        else:
            palette = pd.Series(
                adata.uns[f"{category}_colors"],
                index=adata.obs[category].cat.categories,
            )

    sf = f.subfigures(1, len(embeddings))
    if len(embeddings) == 1:
        sf = [sf]

    for e, obsm in enumerate(embeddings):
        if mini is True:
            ax = sf[e].subplots(1, 1)
            sc.pl.embedding(
                adata,
                basis=obsm,
                color=category,
                ax=ax,
                show=False,
                palette=palette.to_list(),
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
            size = 100000 / adata.shape[0]
            assert nrow is not None
            axs = sf[e].subplots(nrow * 2, nrow)
            gs = axs[0, 0].get_gridspec()
            empty_axs(axs)
            ax = sf[e].add_subplot(gs[:nrow, :])
            sc.pl.embedding(
                adata,
                basis=obsm,
                color=category,
                size=size,
                ax=ax,
                show=False,
                palette=palette.to_list(),
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
            xlims = ax.get_xlim()
            ylims = ax.get_ylim()

            for n, group in enumerate(
                tqdm(
                    adata.obs[category].cat.categories,
                    desc=f"Plotting individual {category}",
                )
            ):
                tmp = adata.obs[category]
                adata.obs["temp"] = tmp.mask(tmp != group)
                ax = sf[e].add_subplot(gs[nrow + n // nrow, n % nrow])
                sc.pl.embedding(
                    adata,
                    basis=obsm,
                    color="temp",
                    size=size,
                    ax=ax,
                    show=False,
                    # palette=[palette[group]],
                    palette=palette.tolist(),
                    legend_loc="none"
                )
                ax.set_title(group)
                ax.set_xlim(*xlims)
                ax.set_ylim(*ylims)
                ax.set(xlabel=None, ylabel=None)

            adata.obs.drop(columns=["temp"], inplace=True)

    return


def check_doublets(
    adata: sc.AnnData,
    embedding: str = "X_umap",
    cluster_key: str = "leiden",
    doubletMethods=["scDblFinder", "DoubletFinder", "doubletdetection", "scrublet"],
):
    f = plt.figure(figsize=(20, 7), layout="constrained")
    sf = f.subfigures(1, 3, width_ratios=(1, 2, 1.5))

    axs = sf[0].subplots(2, 1)
    sc.pl.embedding(
        adata,
        basis=embedding,
        color=cluster_key,
        ax=axs[0],
        show=False,
        palette=create_palette(len(adata.obs[cluster_key].unique())),
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
    adata: sc.AnnData,
    group: str,
    markers: Iterable[str] | dict[str, Iterable[str]],
    f=None,
    layer: str = "normalized",
    useStripPlot: bool = False,
    palette: Iterable = None,
    title = None,
    inner = "box",
    title_fontsize: int = 15,
    y_fontsize: int = 12,
    y_bracket_fontsize: int = 12,
    y_bracket_offset: float = -0.05,
    y_bracket_text_offset: float = 0.05,
    y_blank: bool = False,
    x_fontsize: int = 10,
    x_params={"angle": 60, "align": "right"},
    x_bracket_params = None,
    downscale: float = 1.0
):

    title_fontsize /= downscale
    y_fontsize /= downscale
    x_fontsize /= downscale
    y_bracket_fontsize /= downscale

    if adata.obs[group].dtype != "category":
        print(f"Converting adata.obs[{group}] to `category` dtype.")
        adata.obs[group] = adata.obs[group].astype("category")
    if palette is None and adata.uns.get(f"{group}_colors") is None:
        print("Palette not found! Generating palette.")
        palette = color_gen(adata.obs[group].cat.categories)

    # Convert markers to flat list for plotting
    is_dict = isinstance(markers, dict)
    marker_list = list(markers.values()) if is_dict else [[m] for m in markers]
    flat_markers = [m for group_markers in marker_list for m in group_markers]
    n_groups = len(adata.obs[group].unique())

    # Create subplots
    if f:
        axs = f.subplots(len(flat_markers), 1)
    else:
        f, axs = plt.subplots(
            len(flat_markers),
            1,
            figsize=(n_groups / downscale, len(flat_markers) / downscale),
        )

    if type(axs) is not np.ndarray:
        axs = [axs]

    if title is not None:
        f.suptitle(title, size=title_fontsize)

    # Plot violin plots for each marker
    for n, m in enumerate(tqdm(flat_markers, desc=f"Plotting violins ({title})", smoothing=50/len(flat_markers))):
        sc.pl.violin(
            adata,
            m,
            groupby=group,
            use_raw=False,
            layer=layer,
            show=False,
            ax=axs[n],
            stripplot=useStripPlot,
            palette=palette,
            inner=inner,
        )
        if hasattr(axs[n], "legend_") and axs[n].legend_ is not None:
            axs[n].legend_.remove
        # Hide x-axis labels and ticks for all but bottom plot
        if n < len(flat_markers) - 1:
            axs[n].set_xlabel("")
            axs[n].set_xticks([])
            axs[n].set_xticklabels([])
        if y_blank is False:
            axs[n].set_ylabel(axs[n].get_ylabel(), size=y_fontsize)
        else:
            axs[n].set_ylabel(None)

    # Add left-side brackets for marker sets if markers is a dictionary
    if is_dict and y_blank is False:
        print("Generating gene set brackets!")
        idx = 0
        for set_name, set_markers in markers.items():
            a = idx
            b = idx + len(set_markers) - 1

            # Add brackets
            bracket = mpl.patches.ConnectionPatch(
                xyA=(y_bracket_offset, 1), coordsA=axs[a].transAxes,
                xyB=(y_bracket_offset, 0), coordsB=axs[b].transAxes,
                connectionstyle=f"bar,armA={20/downscale:.2f},armB={20/downscale:.2f},fraction=0",
                linewidth=2, color="black"
            )
            f.add_artist(bracket)

            # Add text
            f.canvas.draw()
            top = axs[a].get_position().y1
            bot = axs[b].get_position().y0
            f.text(y_bracket_text_offset, (top + bot) / 2, set_name, ha="center", va="center", rotation=90, fontsize=y_bracket_fontsize)

            idx += len(set_markers)

    # Format x-axis labels on bottom plot
    axs[-1].set_xticklabels(
        axs[-1].get_xticklabels(),
        size=x_fontsize,
        rotation=x_params.get("angle"),
        ha=x_params.get("align"),
        rotation_mode="anchor",
    )

    # Add bottom x-axis brackets if bracket_params provided
    if x_bracket_params is not None:
        ratios = x_bracket_params["ratio"] / np.sum(x_bracket_params["ratio"])
        ends = np.append(0, np.cumsum(ratios))
        bar_label_locs = [ends[i] + ratios[i] / 2 for i in range(len(ratios))]
        bar_bracket_widths = ratios * f.get_size_inches()[0] * 3.1

        axs = f.get_axes()
        for n, label in enumerate(x_bracket_params["labels"]):
            axs[-1].annotate(
                label,
                xy=(bar_label_locs[n], -x_bracket_params["bracket_y"]),
                xytext=(bar_label_locs[n], -x_bracket_params["label_y"]),
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
        axs[-1].set_xlabel(axs[-1].get_xlabel(), labelpad=x_bracket_params["padding"])

    return f


def plot_cluster_violinplot(
    adata, group: str, clusters: str, markers, f=None, downscale=1.0
):
    n_markers = (
        np.sum(len(m) for m in markers.values())
        if isinstance(markers, dict)
        else len(markers)
    )
    n_groups = len(adata.obs[clusters].unique()) * len(adata.obs[group].unique())
    if f is None:
        # print(n_markers)
        f = plt.figure(
            figsize=((n_groups * 1.2) / downscale, (n_markers * 1.2 + 3) / downscale), layout="constrained"
        )

    clusts = adata.obs[clusters].cat.categories
    cols = color_gen(adata.obs[group].cat.categories, True)
    sf = f.subfigures(
        2,
        len(clusts),
        height_ratios=([len(markers) * 1.2, 3]),
    )

    crosstab_counts = pd.crosstab(adata.obs[clusters], adata.obs[group])
    crosstab_pct = crosstab_counts.div(crosstab_counts.sum(axis=0), axis=1) * 100

    for n, cluster in enumerate(clusts):
        cdata = adata[adata.obs[clusters] == cluster]
        plot_violinplot(
            cdata,
            group,
            markers,
            sf[0, n],
            useStripPlot=False,
            palette=cols.to_list(),
            y_blank=True if n > 0 else False,
            y_bracket_offset=-0.2,
            y_bracket_text_offset=-0.1,
            title=cluster,
            downscale=downscale
        )

        ax = sf[1, n].subplots(1, 1)
        crosstab_pct.loc[cluster][crosstab_pct.loc[cluster] > 0].plot(
            kind="bar", color=cols[crosstab_pct.loc[cluster] > 0], ax=ax
        )
        ax.set_ylim(top=ax.set_ylim()[1] * 1.2)
        ax.tick_params(axis="x", rotation=0)
        ax.bar_label(
            ax.containers[0],
            labels=[
                f"{c:.2f}%\n({crosstab_counts.loc[cluster][n]})"
                for n, c in enumerate(crosstab_pct.loc[cluster])
                if c > 0
            ],
        )

        sf[0, n].suptitle(cluster, size=15)

    f.suptitle(f"{clusters} expression by {group}", size=25)

    return f


def plot_cluster_barplots(
    adata: sc.AnnData,
    group: str,
    clusters: str,
    f,
):
    crosstab_counts = pd.crosstab(adata.obs[clusters], adata.obs[group])
    crosstab_pct = crosstab_counts.div(crosstab_counts.sum(axis=0), axis=1) * 100

    ax = f.subplots(1, 1)
    crosstab_pct.plot(kind="bar", ax=ax, color=color_gen(adata.obs[group]))

    for g, group in enumerate(crosstab_pct.columns):
        ax.bar_label(
            ax.containers[g],
            labels=[
                f"{perc:.2f}%\n({crosstab_counts[group][n]})"
                for n, perc in enumerate(crosstab_pct[group])
            ],
            size=7,
        )

    ax.set_ylim(top=ax.set_ylim()[1] * 1.2)
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=60,
        ha="right",
        rotation_mode="anchor",
    )

    return crosstab_pct


def plot_cluster_counts(
    adata: sc.AnnData,
    clusters: str,
    colors: list = None,
    ax=None,
):
    if colors is None:
        colors = color_gen(adata.obs[clusters].cat.categories)
    if ax is None:
        _, ax = plt.subplots(1, 2, figsize=(5, 5), layout="constrained")

    counts = pd.Series(Counter(adata.obs[clusters])).sort_index()

    pd.Series(counts.plot(kind="bar", color=colors, rot=30, ax=ax))
    ax.set_ylabel("Counts")
    ax.set_title(f"{clusters} counts")
    ax.bar_label(
        ax.containers[0],
        labels=[f"{c/counts.sum()*100:.2f}%\n({c})" for c in counts if c > 0],
    )
    ax.set_ylim(top=ax.get_ylim()[1] * 1.05)

    return


def plot_cluster_stackedbarplot(
    adata: sc.AnnData,
    groupby: str,
    clusters: str,
    pct: bool = False,
    colors: list = None,
    ax=None,
):
    if colors is None:
        colors = color_gen(adata.obs[clusters].cat.categories)
    if ax is None:
        f, ax = plt.subplots(1, 1, figsize=(5, 5), layout="constrained")

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
        for x, y in enumerate(crosstab_counts.sum(axis=1)):
            ax.annotate(y, (x, y * 1.01), ha="center")
        ax.set_ylim(top=ax.get_ylim()[1] * 1.05)

    else:
        # Percentages
        crosstab_pct = crosstab_counts.copy()
        crosstab_pct.loc["all"] = crosstab_pct.sum(axis=0)
        crosstab_pct = crosstab_pct.div(crosstab_pct.sum(axis=1), axis=0) * 100

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


def plot_cluster_trees(
    adata: sc.AnnData,
    clusters: str,
    markers: Iterable[str] = None,
    threshold=0.05,
    title: str = None,
    layer="normalized",
    **kwargs,
):
    if markers is None:
        f = clustree(
            adata, clusters, title=title, edge_weight_threshold=threshold, **kwargs
        )
    else:
        print(f"setting X to layer '{layer}'")
        adata.X = adata.layers[layer].copy()

        if isinstance(markers, str):
            markers = [markers]
        for m in tqdm(markers):
            f = clustree(
                adata,
                clusters,
                title=title,
                edge_weight_threshold=threshold,
                node_color_gene=m,
                node_color_gene_use_raw=False,
                node_colormap="Reds",
                show_colorbar=True,
                **kwargs,
            )
    return f


def plot_cluster_riverplot(
    clusters_1, clusters_2, prefix_1="", prefix_2="", min_flow=50
):
    flow_counts = Counter(zip(clusters_1, clusters_2))
    flow_counts = {k: v for k, v in flow_counts.items() if v >= min_flow}

    clusters1 = sorted(set(clusters_1))
    clusters2 = sorted(set(clusters_2))

    source = [clusters1.index(c1) for c1, c2 in flow_counts]
    target = [len(clusters1) + clusters2.index(c2) for c1, c2 in flow_counts]
    value = list(flow_counts.values())

    counts1 = Counter(clusters_1)
    counts2 = Counter(clusters_2)
    labels = [f"{prefix_1}{c} ({counts1[c]})" for c in clusters1] + [
        f"{prefix_2}{c} ({counts2[c]})" for c in clusters2
    ]

    fig = go.Figure(
        go.Sankey(
            node=dict(label=labels, pad=15, thickness=20),
            link=dict(source=source, target=target, value=value),
        )
    )

    fig.update_layout(
        title=f"Riverplot Comparison of `{clusters_1.name}` and `{clusters_2.name}`",
        font_size=12,
    )
    fig.show()


def plot_cluster_silhouette(
    adata: sc.AnnData,
    obs_key: str = "leiden",
    figsize=(10, 6),
    uns_key={"avg": "silhouette_avg", "score": "silhouette_scores"},
):
    cluster_labels = adata.obs[obs_key].astype("category").cat.codes
    n_clusters = len(np.unique(cluster_labels))

    f, ax = plt.subplots(figsize=figsize)
    y_lower = 10

    silhouette_avg = adata.uns[uns_key["avg"]][obs_key]

    for i in range(n_clusters):
        cluster_silhouette_values = adata.uns[uns_key["score"]][obs_key][
            cluster_labels == i
        ]
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
    adata: sc.AnnData,
    key: str,
    pval: float = 0.05,
    top_n: int = 20,
    sources: Iterable[str] = None,
    targets: Iterable[str] = None,
    includes_genes: Iterable[str] = None,
    figsize=(15, 8),
):
    """Lightweight cell-cell communication dotplot with controlled dot sizes."""

    df = adata.uns[key]
    df = df[df["magnitude_rank"] < pval]

    if sources:
        df = df[df["source"].isin(sources)]
    if targets:
        df = df[df["target"].isin(targets)]
    if includes_genes:
        df = df[df["ligand_complex"].isin(includes_genes) | df["receptor_complex"].isin(includes_genes)]

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
    fig.text(0.5, -0.02, "Target", ha="center")

    fig.subplots_adjust(right=0.75)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap="Greens", norm=norm)
    plt.colorbar(sm, ax=axes, label="Specificity (-log10 p-val)")

    # Create proper size legend using evenly spaced integer values
    if len(mag_vals) > 0:  # Only create legend if we have data
        max_val = int(np.ceil(mag_vals.max()))
        if max_val <= 1:
            size_labels = ["1"]
        else:
            thresholds = np.linspace(1, max_val, 5).astype(int)
            thresholds = np.array(sorted(set(thresholds)))
            size_labels = [str(t) for t in thresholds]

        # Create legend elements
        from matplotlib.lines import Line2D

        legend_elements = []
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
