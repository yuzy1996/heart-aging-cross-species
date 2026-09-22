"""
08_export_hotspot_objects.py
Title : Export labelled hotspot objects for LARIS ligand-receptor analysis
Figure: provides inputs for Fig. 4G-J (notebooks 02/03)
Description:
  Write the gradient-labelled human and mouse (Old_3) objects consumed by the
  LARIS spatial ligand-receptor notebooks. The human object is exported in its
  upstream log-normalised state (matching the original export step); the mouse
  object is normalised and log1p-transformed before export, matching the point
  in the original notebook at which old3_Hspot.h5ad was written.
Inputs : cfg.HUMAN_GRADIENT_H5, cfg.MOUSE_GRADIENT_H5
Outputs: cfg.HUMAN_HSPOT_H5 (adata_human_Hspot.h5ad),
         cfg.MOUSE_HSPOT_H5 (old3_Hspot.h5ad)
"""
import os, sys, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
cfg = importlib.import_module("00_config")
su = importlib.import_module("spatial_utils")

import scanpy as sc


def main():
    adata_human = sc.read_h5ad(cfg.HUMAN_GRADIENT_H5)
    old3 = sc.read_h5ad(cfg.MOUSE_GRADIENT_H5)

    # human: exported before the marker-ranking normalisation
    adata_human.write_h5ad(cfg.HUMAN_HSPOT_H5)

    # mouse: exported after normalize_total + log1p in the original notebook
    sc.pp.normalize_total(old3)
    sc.pp.log1p(old3)
    old3.write_h5ad(cfg.MOUSE_HSPOT_H5)


if __name__ == "__main__":
    main()
