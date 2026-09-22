"""
02_fig3BC_cluster_niche_maps.py
Title : Spatial clusters and senescence-niche score maps (Fig. 3B-C)
Figure: Fig. 3B (unsupervised spatial clusters), Fig. 3C (senescence score)
Description:
  Spatial maps of unsupervised Visium clusters and of the continuous
  senescence score for the aged human section and the representative aged
  mouse section (Old_3). Magnified insets are produced interactively from the
  same AnnData in the companion notebook.
Inputs : cfg.HUMAN_GRADIENT_H5, cfg.MOUSE_GRADIENT_H5
Outputs: human_clusters.pdf, human_niche_score.pdf,
         mouse_old3_clusters.pdf, mouse_old3_sasp.pdf
"""
import os, sys, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
cfg = importlib.import_module("00_config")
su = importlib.import_module("spatial_utils")

import matplotlib.pyplot as plt
import scanpy as sc


def _spatial(adata, color, out, cmap=None, **kw):
    sc.pl.spatial(adata, img_key="hires", color=color, size=1.4,
                  show=False, cmap=cmap, **kw)
    cfg.savefig(plt.gcf(), out)
    plt.close()


def main():
    adata_human = sc.read_h5ad(cfg.HUMAN_GRADIENT_H5)
    old3 = sc.read_h5ad(cfg.MOUSE_GRADIENT_H5)

    # Fig. 3B - unsupervised clusters
    _spatial(adata_human, ["clusters"], "human_clusters.pdf", cmap="BuPu", alpha_img=1)
    _spatial(old3, ["clusters"], "mouse_old3_clusters.pdf", cmap="BuPu", alpha_img=1)

    # Fig. 3C - continuous senescence/niche score
    _spatial(adata_human, [cfg.HUMAN_SCORE_COL], "human_niche_score.pdf", alpha_img=1)
    _spatial(old3, [cfg.MOUSE_SCORE_COL], "mouse_old3_sasp.pdf",
             legend_loc=None, alpha_img=1)


if __name__ == "__main__":
    main()
