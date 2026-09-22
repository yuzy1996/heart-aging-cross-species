"""
00_config.py - Configuration and constants for the Figure 3 spatial niche analysis.

All paths are read from environment variables (defaults match the analysis
server). Set HU_MK_MU_ROOT / HCA_ROOT in your shell or .Renviron-equivalent before
running. Figures, tables and intermediate objects are written to FIG6_DIR,
mirroring the original notebook working directory.
"""
import os

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
HU_ROOT = os.environ.get(
    "HU_MK_MU_ROOT", os.path.expanduser("~/project/cross_species/hu_mk_mu"))
HCA_ROOT = os.environ.get(
    "HCA_ROOT", os.path.expanduser("~/project/heart_Cell_Atlas"))

FIG5_DIR = os.path.join(HU_ROOT, "FIG5")          # integrated human spatial object
FIG6_DIR = os.environ.get("FIG3_OUT_DIR", os.path.join(HU_ROOT, "FIG6"))
os.makedirs(FIG6_DIR, exist_ok=True)

# Inputs
HUMAN_INPUT = os.path.join(FIG5_DIR, "adata.h5ad")
MOUSE_ST_INPUT = os.path.join(
    FIG6_DIR, "Visium_YoungOld_AgingProject_28August_Cleaned.h5ad")
MOUSE_SC_REF = os.path.join(FIG6_DIR, "snRNA_RefAging_Manuscript.h5ad")  # cell2location reference (upstream)

# Intermediate objects produced by 01_* and reused by the panel scripts
HUMAN_GRADIENT_H5 = os.path.join(FIG6_DIR, "adata_human_gradient.h5ad")
MOUSE_GRADIENT_H5 = os.path.join(FIG6_DIR, "old3_gradient.h5ad")

# Hotspot objects exported for the LARIS ligand-receptor notebooks (Fig. 4G-J)
HUMAN_HSPOT_H5 = os.path.join(FIG6_DIR, "adata_human_Hspot.h5ad")
MOUSE_HSPOT_H5 = os.path.join(FIG6_DIR, "old3_Hspot.h5ad")

# ---------------------------------------------------------------------------
# Analysis constants
# ---------------------------------------------------------------------------
RANDOM_SEED = 42

# Senescence score columns already present on the upstream objects
HUMAN_SCORE_COL = "SASP_gene"
MOUSE_SCORE_COL = "senescence"
GRADIENT_COL = "senescence_gradient"

# Ordered distance bins (top-5% hotspot -> concentric 100-500 um rings -> rest)
GRADIENT_ORDER = ["Hspot", "dist100", "dist200", "dist300",
                  "dist400", "dist500", "rest"]
GRADIENT_COLORS = {
    "Hspot": "firebrick", "dist100": "tomato", "dist200": "lightsalmon",
    "dist300": "royalblue", "dist400": "cornflowerblue",
    "dist500": "lightsteelblue", "rest": "sandybrown",
}

# Human radii are derived from the minimum spot-to-spot distance;
# mouse radii are fixed in micrometres (as in the original notebook).
HUMAN_RADIUS_MULTIPLIERS = [1.2, 2.2, 3.2, 4.2, 5.2]
MOUSE_RADII = [100.0, 200.0, 300.0, 400.0, 500.0]

HOTSPOT_PERCENTILE = 95          # top 5% within each aged sample

# Mouse slides available in the multi-slide object; Old_3 is the representative
# aged section used throughout Fig. 3 (the others were explored interactively).
MOUSE_AGED_SLIDE = "Old_3"
MOUSE_AGED_SLIDES = [f"Old_{i}" for i in range(1, 6)]
MOUSE_YOUNG_SLIDES = [f"Young_{i}" for i in range(1, 6)]

# Deconvolution abundance keys produced upstream by cell2location
HUMAN_PROP_KEY = "prop"
MOUSE_PROP_KEY = "c2l_prop"
HUMAN_TOP_CELLTYPES = ["FB4_activated", "MoMP", "LYVE1+IGF1+MP",
                       "SMC2_art", "EC6_ven", "vCM1"]
MOUSE_TOP_CELLTYPES = ["Fibro_activ", "EndoEC", "MP", "Ccr2+MP", "SMC", "CM"]

# Fibrosis / vascular / inflammation markers averaged across the gradient (Fig. 3F)
NICHE_GENES_HUMAN = ["THBS1", "AEBP1", "CD248", "ELN", "MFAP4", "MYH11",
                     "TAGLN", "CCDC3", "PTGIS", "ISLR", "IL33", "IGFBP6", "PGF"]
NICHE_GENES_MOUSE = ["Thbs1", "Aebp1", "Cd248", "Eln", "Mfap4", "Myh11",
                     "Tagln", "Ccdc3", "Ptgis", "Islr", "Il33", "Igfbp6", "Pgf"]

# Curated SASP panel (available for supplementary spatial maps)
SASP_CORE_GENES_MOUSE = ["Il6", "Cxcl8", "Il1a", "Il1b", "Tnf", "Ccl2",
                         "Cxcl1", "Cxcl2", "Mmp2", "Mmp3", "Mmp9",
                         "Tgfb1", "Igfbp3", "Vegfa", "Icam1"]

# Marker-overlap cut-offs for the cross-species niche-marker comparison
MARKER_LOGFC_CUTOFF = 1.2
MARKER_SCORE_CUTOFF = 2.0


def savefig(fig, filename, **kwargs):
    """Save a figure into FIG6_DIR (svg/pdf default to tight bbox)."""
    path = os.path.join(FIG6_DIR, filename)
    fig.savefig(path, bbox_inches="tight", **kwargs)
    return path
