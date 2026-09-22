"""
01_fig3A_senescence_gradient.py
Title : Define the senescence gradient and top-5% hotspots (Fig. 3A)
Figure: Fig. 3A (and provides the labelled objects for Fig. 3B-G)
Description:
  Load the aged human Visium object and the representative aged mouse slide
  (Old_3), then label each spot as Hspot (top 5% senescence score within the
  sample) or as one of five concentric distance bins (dist100-dist500); all
  remaining tissue is 'rest'. Human radii are derived from the minimum spot
  distance; mouse radii are fixed at 100-500 um. Mouse spots are also
  clustered (PCA/UMAP/Leiden) for Fig. 3B; human clusters come from upstream.
Inputs :
  cfg.HUMAN_INPUT  (FIG5/adata.h5ad, with SASP_gene score, clusters and prop)
  cfg.MOUSE_ST_INPUT (multi-slide Visium object, with senescence and c2l_prop)
Outputs:
  cfg.HUMAN_GRADIENT_H5, cfg.MOUSE_GRADIENT_H5 (labelled intermediate objects)
"""
import os, sys, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
cfg = importlib.import_module("00_config")
su = importlib.import_module("spatial_utils")

import scanpy as sc


def main():
    # ---- Human: radii derived from the minimum spot-to-spot distance ----
    adata_human = sc.read_h5ad(cfg.HUMAN_INPUT)
    human_radii = su.radii_from_min_distance(
        adata_human, cfg.HUMAN_RADIUS_MULTIPLIERS)
    print("Human ring radii:", [round(r, 1) for r in human_radii])
    adata_human = su.add_senescence_gradient(
        adata_human, cfg.HUMAN_SCORE_COL, human_radii, cfg.HOTSPOT_PERCENTILE)
    print(adata_human.obs[cfg.GRADIENT_COL].value_counts())
    adata_human.write_h5ad(cfg.HUMAN_GRADIENT_H5)

    # ---- Mouse: fixed 100-500 um radii on the representative aged slide ----
    adata_mouse = sc.read_h5ad(cfg.MOUSE_ST_INPUT)
    old3 = su.select_slide(adata_mouse, cfg.MOUSE_AGED_SLIDE)
    old3 = su.add_senescence_gradient(
        old3, cfg.MOUSE_SCORE_COL, cfg.MOUSE_RADII, cfg.HOTSPOT_PERCENTILE)

    # unsupervised spatial clustering (Fig. 3B, mouse); scanpy default settings
    # match the original notebook (pass random_state=cfg.RANDOM_SEED to pin)
    sc.pp.pca(old3)
    sc.pp.neighbors(old3)
    sc.tl.umap(old3)
    sc.tl.leiden(old3, key_added="clusters")
    print(old3.obs[cfg.GRADIENT_COL].value_counts())
    old3.write_h5ad(cfg.MOUSE_GRADIENT_H5)


if __name__ == "__main__":
    main()
