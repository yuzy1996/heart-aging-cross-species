"""
spatial_utils.py - Reusable functions for the Figure 3 senescence-niche analysis.

Faithfully refactored from 01_mouse_visium_spatial_niche.ipynb; analysis logic
is unchanged. Functions are pure/side-effect-light where possible and radii are
passed explicitly (the notebook used module-level radius variables).
"""
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
from tqdm import tqdm

GRADIENT_LABELS = ["Hspot", "dist100", "dist200", "dist300",
                   "dist400", "dist500", "rest"]


def set_seed(seed=42):
    np.random.seed(seed)


# ---------------------------------------------------------------------------
# Geometry / neighbourhood
# ---------------------------------------------------------------------------
def min_spot_distance(adata):
    """Minimum pairwise distance between Visium spots (coordinate units)."""
    dist = distance_matrix(adata.obsm["spatial"], adata.obsm["spatial"])
    np.fill_diagonal(dist, np.nan)
    return float(np.nanmin(dist))


def radii_from_min_distance(adata, multipliers=(1.2, 2.2, 3.2, 4.2, 5.2)):
    """Five concentric radii derived from the nearest-neighbour spot distance."""
    d = min_spot_distance(adata)
    return [d * m for m in multipliers]


def get_surrounding(adata, in_spot=None, bc_spot=None, radius=100.0,
                    get_bcs=True):
    """Indices (or barcodes) of spots within `radius` of a centre spot."""
    if in_spot is None and bc_spot is None:
        raise ValueError("Specify either in_spot or bc_spot")
    if in_spot is not None:
        centre = adata.obsm["spatial"][in_spot]
    else:
        centre = adata[adata.obs_names == bc_spot, :].obsm["spatial"][0]

    coords = adata.obsm["spatial"]
    d = np.sqrt((coords[:, 0] - centre[0]) ** 2 +
                (coords[:, 1] - centre[1]) ** 2)
    idx = np.where(d <= radius)[0]
    return list(adata.obs_names[idx]) if get_bcs else list(idx)


def add_senescence_gradient(adata, score_col, radii, percentile=95):
    """Label each spot as Hspot (top 5% score) or in concentric distance bins.

    Mirrors add_senescence_gradient_singlesample(): the top `percentile` spots
    are hotspots; surrounding spots fall into five non-overlapping rings;
    everything else is 'rest'. Writes obs['senescence_gradient'].
    """
    scores = adata.obs[score_col]
    threshold = np.percentile(scores, percentile)
    hotspots = scores[scores > threshold].index.tolist()

    annotation = {bc: [10] * 7 for bc in adata.obs_names}
    r100, r200, r300, r400, r500 = radii

    for bc in tqdm(hotspots, desc="senescence gradient"):
        annotation[bc][0] = 0
        ring1 = get_surrounding(adata, bc_spot=bc, radius=r100)
        ring2 = get_surrounding(adata, bc_spot=bc, radius=r200)
        ring3 = get_surrounding(adata, bc_spot=bc, radius=r300)
        ring4 = get_surrounding(adata, bc_spot=bc, radius=r400)
        ring5 = get_surrounding(adata, bc_spot=bc, radius=r500)

        # keep only the nearest (non-overlapping) ring for each spot
        ring2 = [v for v in ring2 if v not in ring1]
        ring3 = [v for v in ring3 if v not in ring1 + ring2]
        ring4 = [v for v in ring4 if v not in ring1 + ring2 + ring3]
        ring5 = [v for v in ring5 if v not in ring1 + ring2 + ring3 + ring4]
        if bc in ring1:
            ring1.remove(bc)  # hotspot is not its own neighbour

        for region, label in [(ring1, 1), (ring2, 2), (ring3, 3),
                              (ring4, 4), (ring5, 5)]:
            for inner in region:
                annotation[inner][label] = label

    annot = pd.DataFrame.from_dict(annotation).T
    annot["code"] = annot.min(axis=1)
    annot["gradient"] = annot["code"].replace({
        0: "Hspot", 1: "dist100", 2: "dist200", 3: "dist300",
        4: "dist400", 5: "dist500", 10: "rest"})

    if "senescence_gradient" in adata.obs.columns:
        del adata.obs["senescence_gradient"]
    adata.obs["senescence_gradient"] = annot["gradient"].values
    return adata


# ---------------------------------------------------------------------------
# Averaging / slide selection
# ---------------------------------------------------------------------------
def average_expression(adata, group_by, feature=None, out_format="long",
                       layer=None):
    """Mean (un-log) expression per group, returned log1p-transformed.

    Equivalent to the notebook AverageExpression(): the log transform is undone
    before averaging and reapplied afterwards.
    """
    import scipy.sparse
    data = adata[:, feature].copy() if feature is not None else adata.copy()
    if layer is not None:
        data.X = data.layers[layer].copy()
    X = data.X
    data.X = X.tocsr().expm1() if scipy.sparse.issparse(X) else np.expm1(X)

    assert out_format in ("wide", "long")
    if isinstance(group_by, str):
        group_by = [group_by]

    main = pd.DataFrame()
    for group_name, df in data.obs.groupby(group_by, as_index=False):
        tmp = pd.DataFrame(np.log1p(np.asarray(data[df.index].X.mean(axis=0)).T),
                           columns=["expr"])
        tmp["gene"] = data[df.index].var_names
        if isinstance(group_name, str):
            group_name = [group_name]
        for i, name in enumerate(group_name):
            tmp[f"group{i}"] = str(name).replace("-", "_")
        main = pd.concat([main, tmp], axis=0)
    main["expr"] = pd.to_numeric(main["expr"])
    main = main[[c for c in main.columns if c != "expr"] + ["expr"]]

    if out_format == "wide":
        main = pd.pivot_table(
            main, index="gene",
            columns=list(main.columns[main.columns.str.startswith("group")]),
            values="expr")
        if len(group_by) > 1:
            main.columns = main.columns.map("_".join)
    return main


def select_slide(adata, slide, slide_col="sample"):
    """Subset a multi-slide spatial AnnData to one experiment, cleaning uns."""
    slid = adata[adata.obs[slide_col].isin([slide]), :].copy()
    for key in list(slid.uns.get("spatial", {}).keys()):
        if key != slide:
            del slid.uns["spatial"][key]
    return slid


# ---------------------------------------------------------------------------
# Figure-3 summarisation / plotting
# ---------------------------------------------------------------------------
def gradient_median_props(adata, prop_key, order=GRADIENT_LABELS):
    """Median deconvolved cell-type abundance per gradient bin (genes x bins)."""
    df = adata.obsm[prop_key].copy()
    df["gradient"] = adata.obs["senescence_gradient"].values
    plot_df = df.groupby("gradient", observed=False).median().T
    cols = [c for c in order if c in plot_df.columns]
    plot_df = plot_df.reindex(columns=cols)
    plot_df = plot_df.dropna(how="all").fillna(0)
    plot_df = plot_df.replace([np.inf, -np.inf], 0)
    return plot_df.loc[plot_df.sum(axis=1) > 0]


def plot_prop_clustermap(plot_df, out_path, top_celltypes=None):
    """Row-z-scored abundance heat map across the gradient (Fig. 3G)."""
    mat = plot_df.loc[top_celltypes] if top_celltypes else plot_df
    sns.clustermap(mat, z_score=0, cmap="RdBu_r", center=0, col_cluster=False,
                   row_cluster=False, figsize=(5.2, 6.3),
                   cbar_kws={"orientation": "horizontal"}, robust=True,
                   linewidth=0.1)
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()


def gradient_cluster_composition(adata, cluster_index, order=GRADIENT_LABELS):
    """Fraction of each spatial cluster assigned to each gradient bin (Fig. 3D)."""
    df = adata.obs.value_counts(["senescence_gradient", "clusters"]).reset_index()
    total = df[["clusters", "count"]].groupby("clusters").sum()
    df["norm"] = [row["count"] / total.loc[row.clusters].values[0]
                  for _, row in df.iterrows()]
    df = df.pivot(index="clusters", columns="senescence_gradient",
                  values="norm")
    cols = [c for c in order if c in df.columns]
    idx = [i for i in cluster_index if i in df.index]
    return df.reindex(columns=cols, index=idx).fillna(0)


def plot_gradient_stackbar(df, out_path, colors):
    """Stacked bar of gradient-bin fractions per spatial cluster (Fig. 3D)."""
    ax = df.plot.bar(stacked=True, color=colors, figsize=(5, 6), width=0.9)
    sns.move_legend(ax, loc="center right", frameon=False,
                    title="Senescence\nGradient",
                    title_fontproperties={"weight": "bold"},
                    bbox_to_anchor=(1.3, 0.5), labelspacing=0.2, fontsize=15)
    ax.set_xlabel("")
    ax.set_xticklabels(ax.get_xticklabels(), fontweight="bold", fontsize=15,
                       va="top", ha="right", rotation=45)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=15)
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()


def scaled_gradient_gene_matrix(adata, genes, order=GRADIENT_LABELS):
    """Min-max scaled mean expression of marker genes across gradient bins (Fig. 3F)."""
    df = average_expression(adata, group_by=["senescence_gradient"],
                            feature=genes, out_format="wide")
    cols = [c for c in order if c in df.columns]
    df = df.reindex(index=genes, columns=cols)
    return df.sub(df.min(axis=1), axis=0).div(df.max(axis=1), axis=0)


def plot_gene_gradient_heatmap(ndf, out_path):
    """Red-scale heat map of scaled marker expression across the gradient (Fig. 3F)."""
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    hm = sns.heatmap(ndf, cmap="Reds", square=True, linewidths=0.1, ax=ax,
                     cbar=False, yticklabels=True, xticklabels=True)
    hm.set_xlabel(""); hm.set_ylabel("")
    hm.set_yticklabels(hm.get_yticklabels(), fontsize=15)
    hm.set_xticklabels([t.get_text().split("_")[0] for t in hm.get_xticklabels()],
                       fontsize=15, fontweight="bold", rotation=45,
                       ha="right", va="top")
    hm.spines[["top", "right", "left", "bottom"]].set_visible(True)
    cbar_ax = fig.add_axes([0.75, 0.2, 0.015, 0.15])
    sm = plt.cm.ScalarMappable(cmap="Reds", norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.ax.set_title("Expr", fontweight="bold", loc="left", fontsize=12)
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()
