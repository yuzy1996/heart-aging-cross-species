"""
03_fig3D_gradient_bin_composition.py
Title : Fraction of each spatial cluster assigned to gradient bins (Fig. 3D)
Figure: Fig. 3D
Description:
  For every unsupervised spatial cluster, compute the proportion of spots
  assigned to each senescence-gradient bin (Hspot, dist100-dist500, rest) and
  render stacked bars for human and mouse.
Inputs : cfg.HUMAN_GRADIENT_H5, cfg.MOUSE_GRADIENT_H5
Outputs: human_Aged_Stack.svg, mouse_old3_Stack.svg
"""
import os, sys, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
cfg = importlib.import_module("00_config")
su = importlib.import_module("spatial_utils")

import scanpy as sc


def main():
    adata_human = sc.read_h5ad(cfg.HUMAN_GRADIENT_H5)
    old3 = sc.read_h5ad(cfg.MOUSE_GRADIENT_H5)

    human_df = su.gradient_cluster_composition(
        adata_human, cluster_index=[str(i) for i in range(7)])
    mouse_df = su.gradient_cluster_composition(
        old3, cluster_index=[str(i) for i in range(8)])

    su.plot_gradient_stackbar(
        human_df, os.path.join(cfg.FIG6_DIR, "human_Aged_Stack.svg"),
        cfg.GRADIENT_COLORS)
    su.plot_gradient_stackbar(
        mouse_df, os.path.join(cfg.FIG6_DIR, "mouse_old3_Stack.svg"),
        cfg.GRADIENT_COLORS)


if __name__ == "__main__":
    main()
