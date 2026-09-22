# =============================================================================
# 14_FB_hdWGCNA_modules.R
# Title : hdWGCNA fibroblast co-expression modules
# Figure: Fig. 2H
# Module: 03_figure2_fibroblast_states
# Description:
#   Soft-power selection, metacell/network construction, module eigengenes and module-feature plots for fibroblasts.
# Inputs:
#   - hdWGCNA_FB.rds
# Outputs:
#   - hu("hdWGCNA", "ModuleNetworks2", "FB_PlotSoftPowers.pdf")
#   - PlotKMEs_FB.pdf
#   - ModuleFeaturePlot1.pdf
#   - ModuleFeaturePlot2.pdf
#   - FB_dot.pdf
#   - ModuleUMAPPlot.pdf
#   - hdWGCNA_FB.rds
#   - ModuleUMAPPlot_species.pdf
#   - PlotDMEsLollipop_speices1.pdf
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
  FB_label_re_cca,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "fibroblast" # the name of the hdWGCNA experiment
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
  group_name = "FB_4", # the name of the group of interest in the group.by column
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

# setwd() removed: paths are resolved via R/00_config.R
pdf(hu("hdWGCNA", "ModuleNetworks2", "FB_PlotSoftPowers.pdf"),width = 10,height = 10)
wrap_plots(plot_list, ncol=2)
dev.off()
power_table <- GetPowerTable(seurat_obj)
head(power_table)
# Construct co-expression network -----------------------------------------
# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'FB_4' # name of the topoligical overlap matrix written to disk
)
PlotDendrogram(seurat_obj, main='FB_4 hdWGCNA Dendrogram')
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
  group.by = 'subtypes', group_name = 'FB_4'
)
# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "FB_4-M"
)
# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)
pdf('PlotKMEs_FB.pdf',width = 10,height = 10)
p
dev.off()
# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
# show the first 6 columns:
head(modules[,1:6])
# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
head(hub_df)

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
pdf('ModuleFeaturePlot1.pdf',width = 10,height = 10)
wrap_plots(plot_list, ncol=6)
dev.off()
# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
pdf('ModuleFeaturePlot2.pdf',width = 10,height = 10)
wrap_plots(plot_list, ncol=6)
dev.off()

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
pdf('FB_dot.pdf',width = 10,height = 10)
p
dev.off()

# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(magrittr)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# network analysis & visualization package:
library(igraph)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
options(future.globals.maxSize = 200 * 1024^3)
set.seed(12345)
ModuleNetworkPlot(
  seurat_obj, 
  outdir='ModuleNetworks3', # new folder name
  n_inner = 10, # number of genes in inner ring
  n_outer = 18, # number of genes in outer ring
  n_conns = Inf, # show all of the connections
  plot_size=c(7,7), # larger plotting area
  vertex.label.cex=1 # font size
)
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)
pdf('ModuleUMAPPlot.pdf',width = 10,height = 10)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=8 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)
dev.off()

write_rds(seurat_obj, file='hdWGCNA_FB.rds')
seurat_obj <- read_rds('hdWGCNA_FB.rds')


run_DME <- function(seurat_obj){
  group1 <- seurat_obj@meta.data %>% subset(subtypes == 'FB_4' & aged %in% 'Aged') %>% rownames
  group2 <- seurat_obj@meta.data %>% subset(subtypes == 'FB_4' & aged %in% 'Young') %>% rownames
  DMEs <- FindDMEs(
    seurat_obj,
    features = 'ModuleScores',
    barcodes1 = group1,
    barcodes2 = group2,
    test.use='wilcox',
    wgcna_name='fibroblast'
  )
  DMEs$module <- str_replace_all(DMEs$module,'FB-4','FB_4')
  return(DMEs)
}

seurat_obj_human <- subset(seurat_obj,species%in%'human')
seurat_obj_mouse <- subset(seurat_obj,species%in%'mouse')
seurat_obj_monkey <- subset(seurat_obj,species%in%'monkey')
DMEs_human <- run_DME(seurat_obj_human)
DMEs_mouse <- run_DME(seurat_obj_mouse)
DMEs_monkey <- run_DME(seurat_obj_monkey)
pdf('ModuleUMAPPlot_species.pdf',width = 10,height = 10)
ModuleUMAPPlot(
  seurat_obj_human,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=8 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)
ModuleUMAPPlot(
  seurat_obj_monkey,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=8 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)
ModuleUMAPPlot(
  seurat_obj_mouse,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=8 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)
dev.off()
pdf('PlotDMEsLollipop_speices1.pdf',width = 5,height = 5)
PlotDMEsLollipop(
  seurat_obj_human, 
  DMEs_human, 
  wgcna_name='fibroblast', 
  pvalue = "p_val"
)
PlotDMEsLollipop(
  seurat_obj_monkey, 
  DMEs_monkey, 
  wgcna_name='fibroblast', 
  pvalue = "p_val"
)
PlotDMEsLollipop(
  seurat_obj_mouse, 
  DMEs_mouse, 
  wgcna_name='fibroblast', 
  pvalue = "p_val"
)
dev.off()

# select TF of interest
cur_tfs <- c('RUNX2', 'RXRA', 'TCF4')

# plot with default settings
p1 <- TFNetworkPlot(
  seurat_obj_mouse, selected_tfs=cur_tfs, 
  target_type='positive', 
  label_TFs=0, depth=1
) + ggtitle("positive targets")
p2 <- TFNetworkPlot(
  seurat_obj, selected_tfs=cur_tfs, 
  target_type = 'both', 
  label_TFs=0, depth=1
) + ggtitle("pos & neg targets")

p3 <- TFNetworkPlot(
  seurat_obj, selected_tfs=cur_tfs, 
  target_type = 'negative', 
  label_TFs=0, depth=1
) + ggtitle("negative targets")

p1 | p2 | p3

