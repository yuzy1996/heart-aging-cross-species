"""
06_fig3G_niche_deconvolution.py
Title : Cell-type abundance across the senescence gradient (Fig. 3G)
Figure: Fig. 3G
Description:
  Median deconvolved cell-type abundance (human cell2location 'prop';
  mouse 'c2l_prop', both computed upstream) in each gradient bin, shown as a
  row-z-scored heat map for the principal niche cell types.
Inputs : cfg.HUMAN_GRADIENT_H5, cfg.MOUSE_GRADIENT_H5
Outputs: Heatmap_human_prop.svg, Heatmap_mouse_prop.svg
"""
import os, sys, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
cfg = importlib.import_module("00_config")
su = importlib.import_module("spatial_utils")

import scanpy as sc


def main():
    adata_human = sc.read_h5ad(cfg.HUMAN_GRADIENT_H5)
    old3 = sc.read_h5ad(cfg.MOUSE_GRADIENT_H5)

    human_props = su.gradient_median_props(adata_human, cfg.HUMAN_PROP_KEY)
    mouse_props = su.gradient_median_props(old3, cfg.MOUSE_PROP_KEY)

    su.plot_prop_clustermap(
        human_props, os.path.join(cfg.FIG6_DIR, "Heatmap_human_prop.svg"),
        top_celltypes=cfg.HUMAN_TOP_CELLTYPES)
    su.plot_prop_clustermap(
        mouse_props, os.path.join(cfg.FIG6_DIR, "Heatmap_mouse_prop.svg"),
        top_celltypes=cfg.MOUSE_TOP_CELLTYPES)


if __name__ == "__main__":
    main()
