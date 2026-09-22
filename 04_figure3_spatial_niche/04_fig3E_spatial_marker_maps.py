"""
04_fig3E_spatial_marker_maps.py
Title : Spatial expression of niche markers THBS1 / MYH11 / IGFBP6 (Fig. 3E)
Figure: Fig. 3E (plus a supplementary multi-marker panel)
Description:
  Spatial maps of the fibroblast/vascular niche markers THBS1, MYH11 and
  IGFBP6 in the aged human section, and Thbs1/Myh11/Igfbp6 in the aged mouse
  section (Old_3). A broader marker panel (activated fibroblast, receptor,
  SASP and stress markers) is also rendered for human.
Inputs : cfg.HUMAN_GRADIENT_H5, cfg.MOUSE_GRADIENT_H5
Outputs: fig3E_human_markers.pdf, fig3E_mouse_markers.pdf,
         human_marker_panel.pdf (supplementary)
"""
import os, sys, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
cfg = importlib.import_module("00_config")
su = importlib.import_module("spatial_utils")

import matplotlib.pyplot as plt
import scanpy as sc


def _spatial(adata, genes, out, library_id=None, alpha_img=0.0, size=1.4):
    kw = dict(img_key="hires", color=genes, size=size, show=False,
              alpha_img=alpha_img)
    if library_id is not None:
        kw["library_id"] = library_id
    sc.pl.spatial(adata, **kw)
    cfg.savefig(plt.gcf(), out)
    plt.close()


def main():
    adata_human = sc.read_h5ad(cfg.HUMAN_GRADIENT_H5)
    old3 = sc.read_h5ad(cfg.MOUSE_GRADIENT_H5)

    # Fig. 3E
    _spatial(adata_human, ["THBS1", "MYH11", "IGFBP6"],
             "fig3E_human_markers.pdf", alpha_img=0.0, size=1.5)
    _spatial(old3, ["Thbs1", "Myh11", "Igfbp6"],
             "fig3E_mouse_markers.pdf", library_id=cfg.MOUSE_AGED_SLIDE,
             alpha_img=0.0, size=1.0)

    # Supplementary multi-marker panel (equivalent to the exploratory cell)
    _spatial(adata_human,
             ["THBS1", "FB4_activated", "CD36", "vCM3_stressed",
              "NPPB", "COL1A1", "IL6", "MMP3"],
             "human_marker_panel.pdf", alpha_img=1.0, size=1.5)


if __name__ == "__main__":
    main()
