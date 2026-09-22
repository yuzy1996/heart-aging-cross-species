"""
05_fig3F_gradient_gene_heatmap.py
Title : Niche-marker expression along the senescence gradient (Fig. 3F)
Figure: Fig. 3F
Description:
  Mean expression of fibrosis/vascular/inflammation marker genes in each
  gradient bin (log1p mean, undone/reapplied as in AverageExpression), scaled
  per gene to [0,1] across bins, shown as a red heat map for human and mouse.
Inputs : cfg.HUMAN_GRADIENT_H5, cfg.MOUSE_GRADIENT_H5
Outputs: Heatmap_human.svg, Heatmap_mouse.svg
"""
import os, sys, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
cfg = importlib.import_module("00_config")
su = importlib.import_module("spatial_utils")

import scanpy as sc


def main():
    adata_human = sc.read_h5ad(cfg.HUMAN_GRADIENT_H5)
    old3 = sc.read_h5ad(cfg.MOUSE_GRADIENT_H5)

    human_mat = su.scaled_gradient_gene_matrix(
        adata_human, cfg.NICHE_GENES_HUMAN)
    mouse_mat = su.scaled_gradient_gene_matrix(
        old3, cfg.NICHE_GENES_MOUSE)

    su.plot_gene_gradient_heatmap(
        human_mat, os.path.join(cfg.FIG6_DIR, "Heatmap_human.svg"))
    su.plot_gene_gradient_heatmap(
        mouse_mat, os.path.join(cfg.FIG6_DIR, "Heatmap_mouse.svg"))


if __name__ == "__main__":
    main()
