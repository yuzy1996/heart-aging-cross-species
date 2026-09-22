# =============================================================================
# 15_EC_hdWGCNA_modules.R
# Title : hdWGCNA endothelial modules
# Figure: Supplementary
# Module: 03_figure2_fibroblast_states
# Description:
#   hdWGCNA co-expression modules for endothelial cells.
# Inputs:
#   - (see script)
# Outputs:
#   - hu("hdWGCNA", "EC_PlotSoftPowers.pdf")
#   - PlotKMEs.pdf
#   - hdWGCNA_EC.rds
#   - EC_dot.pdf
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

# single-cell analysis package
library(Seurat)
# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# set random seed for reproducibility
set.seed(12345)
# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)
table(EC_cca_re$orig.ident)
table(EC_cca_re$orig.ident)
table(FB_label_re_cca$subtypes)
table(FB_label_re_cca$species)
p <- DimPlot(FB_label_re_cca, group.by='subtypes', label=TRUE) +
  umap_theme() + ggtitle('fibroblast') + NoLegend()

p
seurat_obj <- SetupForWGCNA(
  EC_cca_re,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "EC" # the name of the hdWGCNA experiment
)
# Construct metacells -----------------------------------------------------

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("subtypes", "species"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'CCA', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'subtypes' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "artery_ec", # the name of the group of interest in the group.by column
  group.by='subtypes', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)
# Select soft-power threshold ---------------------------------------------
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)
# assemble with patchwork
pdf(hu("hdWGCNA", "EC_PlotSoftPowers.pdf"),width = 10,height = 10)
wrap_plots(plot_list, ncol=2)
dev.off()
power_table <- GetPowerTable(seurat_obj)
head(power_table)
# setwd() removed: paths are resolved via R/00_config.R
# Construct co-expression network -----------------------------------------
# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'artery_ec' # name of the topoligical overlap matrix written to disk
)
PlotDendrogram(seurat_obj, main='artery_ec hdWGCNA Dendrogram')
TOM <- GetTOM(seurat_obj)

# Module Eigengenes and Connectivity --------------------------------------

# need to run ScaleData first or else harmony throws an error:
#seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# Compute harmonized module eigengenes ------------------------------------
# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="species"
)
# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# Compute module connectivity ---------------------------------------------

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'subtypes', group_name = 'artery_ec'
)
# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "artery_ec-M"
)
# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)
pdf('PlotKMEs.pdf',width = 10,height = 10)
p
dev.off()
# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
head(hub_df)
saveRDS(seurat_obj, file='hdWGCNA_EC.rds')
# Compute hub gene signature scores ---------------------------------------
# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)
# Basic Visualization -----------------------------------------------------

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)
# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

seurat_obj$cluster <- do.call(rbind, strsplit(as.character(seurat_obj$annotation), ' '))[,1]

ModuleRadarPlot(
  seurat_obj,
  group.by = 'cluster',
  barcodes = seurat_obj@meta.data %>% subset(cell_type == 'INH') %>% rownames(),
  axis.label.size=4,
  grid.label.size=4
)
# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'subtypes')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
pdf('EC_dot.pdf',width = 10,height = 10)
p
dev.off()