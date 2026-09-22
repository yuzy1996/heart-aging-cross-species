"""
07_cross_species_niche_markers.py
Title : Conserved human-mouse senescence-niche marker genes (supplementary)
Figure: supports Fig. 3 / conserved-niche comparison
Description:
  On a freshly normalised copy of each labelled object, rank genes across the
  gradient bins (Wilcoxon), extract Hspot markers, and intersect the
  lower-cased human and mouse marker sets (logFC > 1.2, score > 2) to define
  conserved niche markers. Writes per-species marker tables, a Venn diagram and
  a four-sheet Excel workbook. The normalisation/log1p here is applied to a
  copy and only feeds the marker ranking (as in the original notebook).
Inputs : cfg.HUMAN_GRADIENT_H5, cfg.MOUSE_GRADIENT_H5
Outputs: df_human.xlsx, df_mouse.xlsx, marker_venn_diagram.svg,
         marker_genes_all.xlsx
"""
import os, sys, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
cfg = importlib.import_module("00_config")
su = importlib.import_module("spatial_utils")

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from venn import venn


def _ranked_hotspot_markers(adata):
    """Normalised copy -> Wilcoxon ranking -> Hspot marker table."""
    d = adata.copy()
    sc.pp.normalize_total(d)
    sc.pp.log1p(d)
    sc.tl.rank_genes_groups(d, groupby=cfg.GRADIENT_COL, method="wilcoxon")
    return d, sc.get.rank_genes_groups_df(d, group="Hspot")


def _marker_set(df, fc_cutoff, score_cutoff):
    d = df[(df["logfoldchanges"] > fc_cutoff) &
           (df["scores"] > score_cutoff)]
    return {g.lower() for g in d["names"].tolist()}


def main():
    su.set_seed(cfg.RANDOM_SEED)
    adata_human = sc.read_h5ad(cfg.HUMAN_GRADIENT_H5)
    old3 = sc.read_h5ad(cfg.MOUSE_GRADIENT_H5)

    adata_human, df_human = _ranked_hotspot_markers(adata_human)
    old3, df_mouse = _ranked_hotspot_markers(old3)

    df_human[df_human["logfoldchanges"] > 1].to_excel(
        os.path.join(cfg.FIG6_DIR, "df_human.xlsx"), index=False)
    df_mouse[df_mouse["logfoldchanges"] > 1].to_excel(
        os.path.join(cfg.FIG6_DIR, "df_mouse.xlsx"), index=False)

    markers_human = _marker_set(
        df_human, cfg.MARKER_LOGFC_CUTOFF, cfg.MARKER_SCORE_CUTOFF)
    markers_mouse = _marker_set(
        df_mouse, cfg.MARKER_LOGFC_CUTOFF, cfg.MARKER_SCORE_CUTOFF)

    plt.figure(figsize=(7, 5))
    venn({"Human dataset": markers_human, "Mouse Old_3": markers_mouse})
    plt.title("Marker genes overlap (lowercase)", fontsize=14)
    plt.savefig(os.path.join(cfg.FIG6_DIR, "marker_venn_diagram.svg"),
                bbox_inches="tight", dpi=300)
    plt.close()

    common = sorted(markers_human & markers_mouse)
    only_human = sorted(markers_human - markers_mouse)
    only_mouse = sorted(markers_mouse - markers_human)
    all_markers = sorted(markers_human | markers_mouse)
    print(f"common={len(common)} only_human={len(only_human)} "
          f"only_mouse={len(only_mouse)}")

    with pd.ExcelWriter(os.path.join(cfg.FIG6_DIR, "marker_genes_all.xlsx"),
                        engine="openpyxl") as writer:
        pd.DataFrame(common, columns=["common_markers"]).to_excel(
            writer, sheet_name="common", index=False)
        pd.DataFrame(only_human, columns=["only_human"]).to_excel(
            writer, sheet_name="only_human", index=False)
        pd.DataFrame(only_mouse, columns=["only_old3"]).to_excel(
            writer, sheet_name="only_old3", index=False)
        pd.DataFrame(all_markers, columns=["all_markers"]).to_excel(
            writer, sheet_name="all_markers", index=False)


if __name__ == "__main__":
    main()
