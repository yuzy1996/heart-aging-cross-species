# Figure 3 — spatial senescence niche

Spatial transcriptomic definition of the senescence hotspot (Hspot) and its
concentric distance bins, the niche composition, spatial markers and the
human–mouse conserved niche-marker comparison.

The analysis is provided both as **modular Python scripts** (this directory,
recommended for reproduction) and as the original interactive notebook
[`../../notebooks/01_mouse_visium_spatial_niche.ipynb`](../../notebooks/01_mouse_visium_spatial_niche.ipynb).
The two are the same analysis; the scripts remove duplicated cells, pass radii
explicitly, write outputs to fixed paths and add a command-line driver.

## Run

```bash
# set the data root (defaults to ~/project/cross_species/hu_mk_mu)
export HU_MK_MU_ROOT=/path/to/hu_mk_mu
python run_figure3.py          # all steps, in order
python run_figure3.py 03 06    # individual steps (after step 01 has run)
```

All figures/tables/intermediate objects are written to `$HU_MK_MU_ROOT/FIG6`
(override with `FIG3_OUT_DIR`).

## Scripts

| Step | Script | Figure | Content |
|---|---|---|---|
| 01 | `01_fig3A_senescence_gradient.py` | Fig. 3A | Top-5% Hspot per aged sample + five distance bins; mouse PCA/UMAP/Leiden clusters |
| 02 | `02_fig3BC_cluster_niche_maps.py` | Fig. 3B–C | Spatial cluster maps and continuous senescence-score maps |
| 03 | `03_fig3D_gradient_bin_composition.py` | Fig. 3D | Stacked bars: gradient-bin fraction within each spatial cluster |
| 04 | `04_fig3E_spatial_marker_maps.py` | Fig. 3E | Spatial THBS1/MYH11/IGFBP6 (human) and Thbs1/Myh11/Igfbp6 (mouse) |
| 05 | `05_fig3F_gradient_gene_heatmap.py` | Fig. 3F | Min–max scaled niche-marker expression across the gradient |
| 06 | `06_fig3G_niche_deconvolution.py` | Fig. 3G | Median deconvolved cell-type abundance across the gradient |
| 07 | `07_cross_species_niche_markers.py` | supp. | Hspot marker ranking, human–mouse Venn overlap and marker tables |
| 08 | `08_export_hotspot_objects.py` | feeds Fig. 4G–J | Writes `adata_human_Hspot.h5ad` / `old3_Hspot.h5ad` for LARIS |

Shared configuration and constants live in `00_config.py`; reusable functions
in `spatial_utils.py`.

## Inputs (produced upstream)

- `FIG5/adata.h5ad` — aged human Visium object already carrying the `SASP_gene`
  score, spatial `clusters` and the cell2location abundance frame `obsm['prop']`.
- `FIG6/Visium_YoungOld_AgingProject_28August_Cleaned.h5ad` — multi-slide mouse
  Visium object carrying the `senescence` score and `obsm['c2l_prop']`.
  `Old_3` is the representative aged section used throughout Fig. 3; the other
  aged/young slides are available in the object but were not used for the
  figure panels.
- `FIG6/snRNA_RefAging_Manuscript.h5ad` — single-nucleus reference used by the
  upstream cell2location run (read for traceability; deconvolution is not
  re-run here).

## Method notes (preserved from the notebook)

- Within each aged sample, spots above the 95th percentile of the senescence
  score are Hspots; surrounding spots fall into five non-overlapping rings.
  Human ring radii are the minimum spot distance × 1.2/2.2/3.2/4.2/5.2; mouse
  radii are fixed at 100–500 µm. Remaining tissue is `rest`.
- `AverageExpression` undoes log1p before the group mean and reapplies it;
  the marker ranking in step 07 uses a freshly normalised/log1p **copy** and
  does not alter the objects used by the other panels.

## Not covered here

- Fig. 3H is a schematic (no code).
- The senescence scores (`SASP_gene`/`senescence`) and cell2location
  abundances (`prop`/`c2l_prop`) are computed upstream; only the cell2location
  reference is read here.
- The LARIS ligand–receptor inference on the exported objects is in
  `notebooks/02_human_LARIS_spatial_LR.ipynb` and
  `notebooks/03_mouse_LARIS_spatial_LR.ipynb` (Fig. 4G–J).
