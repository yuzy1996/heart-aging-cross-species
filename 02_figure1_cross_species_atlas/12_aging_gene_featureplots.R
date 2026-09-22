# =============================================================================
# 12_aging_gene_featureplots.R
# Title : Feature plots of representative aging genes
# Figure: Supplementary
# Module: 02_figure1_cross_species_atlas
# Description:
#   FeaturePlot of ICAM1, NR4A1, POU2F2, EPHB1 and MYOCD across species and age.
# Inputs:
#   - (see script)
# Outputs:
#   - (see script)
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

combine_species_new_label_human
library(Seurat)
FeaturePlot(combine_species_new_label_human,
            split.by = 'aged',
            features = c('ICAM1','NR4A1','POU2F2','EPHB1','MYOCD'),order = T)

FeaturePlot(combine_species_new_label_monkey,
            split.by = 'aged',
            features = c('ICAM1','NR4A1','POU2F2','EPHB1','MYOCD'),order = T)

FeaturePlot(combine_species_new_label_mouse,
            split.by = 'aged',
            features = c('ICAM1','NR4A1','POU2F2','EPHB1','MYOCD'),order = T)
DotPlot(combine_species_new_label_human,
        features = c('EPHB1'),
        split.by = 'aged')


VlnPlot(combine_species_new_label_human)
DotPlot(combine_species_new_label_human,
        features = c('ICAM1','NR4A3','POU2F2'),
        split.by = 'aged',
        group.by = 'celltypes')
DotPlot(combine_species_new_label_human,
        features = c('ICAM1','NR4A3','POU2F2'),
        split.by = 'aged',
        group.by = 'celltypes')
DotPlot(combine_species_new_label_human,
        features = c('ICAM1','NR4A3','POU2F2'),
        split.by = 'aged',
        group.by = 'celltypes')


DotPlot(combine_species_new_label_human,
        features = c('ICAM1','NR4A3','POU2F2'),
        split.by = 'aged',
        group.by = 'celltypes')
DotPlot(combine_species_new_label_human,
        features = c('ICAM1','NR4A3','POU2F2'),
        split.by = 'aged')
DotPlot(combine_species_new_label_human,
        features = c('ICAM1','NR4A3','POU2F2'),
        split.by = 'aged',
        group.by = 'celltypes')
DotPlot(combine_species_new_label_monkey,
        features = c('VCAN','VCAM1'),
        split.by = 'aged',
        group.by = 'celltypes')
DotPlot(combine_species_new_label_mouse,
        features = c('VCAN','VCAM1','THBS1'),
        split.by = 'aged',
        group.by = 'celltypes')