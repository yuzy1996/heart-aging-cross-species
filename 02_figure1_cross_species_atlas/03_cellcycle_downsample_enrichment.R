# =============================================================================
# 03_cellcycle_downsample_enrichment.R
# Title : Cell-cycle scoring, down-sampling, enrichment networks
# Figure: Supplementary Fig. 3
# Module: 02_figure1_cross_species_atlas
# Description:
#   Cell-cycle scoring, balanced down-sampling of human nuclei, and enrichment map/word-cloud for age DEGs.
# Inputs:
#   - hu("res", "combine_species_new_label.rds")
#   - hu("res", "down_sub.rds")
#   - hu("res", "enrich_sub.rds")
#   - hu("res", "combine_species_new_label_S4_sub.rds")
# Outputs:
#   - hu("res", "down_sub.rds")
#   - hu("fig", "V1", "fig1_enrich_network_up.pdf")
#   - hu("fig", "V1", "fig1_enrich_map_up.pdf")
#   - hu("fig", "V1", "fig1_enrich_wordcloud_up.pdf")
#   - hu("fig", "V1", "fig1_enrich_network_down.pdf")
#   - hu("fig", "V1", "fig1_enrich_map_down.pdf")
#   - hu("fig", "V1", "fig1_enrich_wordcloud_fea_up.pdf")
#   - hu("res", "enrich_sub.rds")
#   - hu("fig", "V1", "fig1_enrich_map.pdf")
#   - hu("fig", "V1", "fig1_dot_marker.pdf")
#   - hu("res", "combine_species_new_label_S4_sub.rds")
#   - hu("fig", "V1", "fig1_feature.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

library(SCP)
library(tidyverse)
options(future.globals.maxSize = 10000* 1024^2)
# (R library search paths are managed by renv / .Renviron, not set per-script)
combine_species_new_label_S4 <- read_rds(hu("res", "combine_species_new_label.rds"))
combine_species_new_label_S4[['RNA']] <- as(combine_species_new_label_S4[['RNA']],Class = 'Assay')
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
combine_species_new_label_S4 <- CellCycleScoring(combine_species_new_label_S4, 
                                                 s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DefaultAssay(combine_species_new_label_S4) <- 'RNA'
Idents(combine_species_new_label_S4) <- 'annotation'
down_sub <- subset(combine_species_new_label_S4,downsample=2000)
write_rds(down_sub,file = hu("res", "down_sub.rds"))

down_sub <- readRDS(file = hu("res", "down_sub.rds"))

down_sub_human_Aged <- subset(down_sub,dataset%in%'human')
down_submonkey_Aged <- subset(down_sub,dataset%in%'monkey')
down_sub_mouse1_Aged <- subset(down_sub,dataset%in%'mouse1')

down_sub_human_Aged <- RunDEtest(down_sub_human_Aged,group_by = 'aged',fc.threshold = 0.8)
down_submonkey_Aged <- RunDEtest(down_submonkey_Aged,group_by = 'aged',fc.threshold = 0.8)
down_sub_mouse1_Aged <- RunDEtest(down_sub_mouse1_Aged,group_by = 'aged',fc.threshold = 0.8)
enrich_down_sub_human_Aged <- RunEnrichment(down_sub_human_Aged,group_by = 'aged',
                            db = 'GO_BP',species = 'Homo_sapiens',
                            DE_threshold = "avg_log2FC > 1 & p_val_adj < 0.05")
enrich_down_submonkey_Aged <- RunEnrichment(down_submonkey_Aged,group_by = 'aged',
                            db = 'GO_BP',species = 'Homo_sapiens',
                            DE_threshold = "avg_log2FC > 1 & p_val_adj < 0.05")
enrich_down_sub_mouse1_Aged <- RunEnrichment(down_sub_mouse1_Aged,group_by = 'aged',
                            db = 'GO_BP',species = 'Homo_sapiens',
                            DE_threshold = "avg_log2FC > 1 & p_val_adj < 0.05")
pdf(hu("fig", "V1", "fig1_enrich_network_up.pdf"),width = 12,height = 12)
EnrichmentPlot(srt = enrich_down_sub_human_Aged, group_by = "aged", group_use = "Aged",plot_type = "network")
EnrichmentPlot(srt = enrich_down_submonkey_Aged, group_by = "aged", group_use = "Aged",plot_type = "network")
EnrichmentPlot(srt = enrich_down_sub_mouse1_Aged, group_by = "aged", group_use = "Aged",plot_type = "network")
dev.off()
pdf(hu("fig", "V1", "fig1_enrich_map_up.pdf"),width = 12,height = 12)
EnrichmentPlot(srt = enrich_down_sub_human_Aged, group_by = "aged", group_use = "Aged",plot_type = "enrichmap")
EnrichmentPlot(srt = enrich_down_submonkey_Aged, group_by = "aged", group_use = "Aged",plot_type = "enrichmap")
EnrichmentPlot(srt = enrich_down_sub_mouse1_Aged, group_by = "aged", group_use = "Aged",plot_type = "enrichmap")
dev.off()
pdf(hu("fig", "V1", "fig1_enrich_wordcloud_up.pdf"),width = 12,height = 12)
EnrichmentPlot(srt = enrich_down_sub_human_Aged, group_by = "aged", group_use = "Aged",plot_type = "wordcloud")
EnrichmentPlot(srt = enrich_down_submonkey_Aged, group_by = "aged", group_use = "Aged",plot_type = "wordcloud")
EnrichmentPlot(srt = enrich_down_sub_mouse1_Aged, group_by = "aged", group_use = "Aged",plot_type = "wordcloud")
dev.off()
pdf(hu("fig", "V1", "fig1_enrich_network_down.pdf"),width = 12,height = 12)
EnrichmentPlot(srt = enrich_down_sub_human_Aged, group_by = "aged", group_use = "Young",plot_type = "network")
EnrichmentPlot(srt = enrich_down_submonkey_Aged, group_by = "aged", group_use = "Young",plot_type = "network")
EnrichmentPlot(srt = enrich_down_sub_mouse1_Aged, group_by = "aged", group_use = "Young",plot_type = "network")
dev.off()
pdf(hu("fig", "V1", "fig1_enrich_map_down.pdf"),width = 12,height = 12)
EnrichmentPlot(srt = enrich_down_sub_human_Aged, group_by = "aged", group_use = "Young",plot_type = "enrichmap")
EnrichmentPlot(srt = enrich_down_submonkey_Aged, group_by = "aged", group_use = "Young",plot_type = "enrichmap")
EnrichmentPlot(srt = enrich_down_sub_mouse1_Aged, group_by = "aged", group_use = "Young",plot_type = "enrichmap")
dev.off()
pdf(hu("fig", "V1", "fig1_enrich_wordcloud_up.pdf"),width = 12,height = 12)
EnrichmentPlot(srt = enrich_down_sub_human_Aged, group_by = "aged", group_use = "Young",plot_type = "wordcloud")
EnrichmentPlot(srt = enrich_down_submonkey_Aged, group_by = "aged", group_use = "Young",plot_type = "wordcloud")
EnrichmentPlot(srt = enrich_down_sub_mouse1_Aged, group_by = "aged", group_use = "Young",plot_type = "wordcloud")
dev.off()
pdf(hu("fig", "V1", "fig1_enrich_wordcloud_fea_up.pdf"),width = 12,height = 12)
EnrichmentPlot(srt = enrich_down_sub_human_Aged, group_by = "aged",
               group_use = "Aged",plot_type = "wordcloud",word_type = "feature")
EnrichmentPlot(srt = enrich_down_submonkey_Aged, group_by = "aged",
               group_use = "Aged",plot_type = "wordcloud",word_type = "feature")
EnrichmentPlot(srt = enrich_down_sub_mouse1_Aged, group_by = "aged", 
               group_use = "Aged",plot_type = "wordcloud",word_type = "feature")
dev.off()
enrich_sub <- RunEnrichment(down_sub,group_by = 'annotation',
                            db = 'GO_BP',species = 'Homo_sapiens',
                            DE_threshold = "avg_log2FC > 1 & p_val_adj < 0.05")
saveRDS(enrich_sub,file = hu("res", "enrich_sub.rds"))
enrich_sub <- read_rds(hu("res", "enrich_sub.rds"))
pdf(hu("fig", "V1", "fig1_enrich_map.pdf"),width = 12,height = 12)
EnrichmentPlot(
  srt = enrich_sub, group_by = "annotation", group_use = c("FB", "EC"),
  plot_type = "wordcloud"
)
EnrichmentPlot(
  srt = enrich_sub, group_by = "annotation", group_use = "FB",
  plot_type = "network"
)
EnrichmentPlot(
  srt = enrich_sub, group_by = "annotation", group_use = "FB",
  plot_type = "enrichmap"
)
dev.off()
pdf(hu("fig", "V1", "fig1_dot_marker.pdf"),width = 6,height = 9)
SCP::GroupHeatmap(srt = combine_species_new_label_S4,
                  features =  c("DCN","COL1A1","GSN","FBLN1","TCF21",#FB
                                "CDH5","PECAM1","RAMP2","EMCN","TEK",#EC
                                "CD163","CD68","MS4A6A","CSF1R",#macrophage
                                "CD3D","CD3E","CD4","CD8A","IL7R",#T
                                "IGKC","MS4A1","CD19","CD79A",#B
                                "RGS5","ABCC9","KCNJ8",#Pericytes
                                "MYH11","TAGLN","ACTA2",#SMC
                                "PLP1","NRXN1","NRXN3",#Neural
                                "ADIPOQ","PLIN1","GPD1","PNPLA3","PCK1"#Adipocyte
                                ),
                  cell_annotation = c('nFeature_RNA','Phase'),
                  group.by = 'annotation',heatmap_palette = "YlOrRd",
                  group_palcolor  = list(celltypes_col),
                  show_row_names = T,row_names_side = 'left',
                  features_label = F,
                  show_column_names = F)
dev.off()
celltypes_col <- c("#FF1493", "#FFB6C1", "#483D8B", "#00BFFF", "#228B22", "#7FFF00", "#FFFF00", "#FF4500","#8DD3C7")
names(celltypes_col) <- celltypes
species_col <- c("#FC8D62","#8DD3C7","#BC80BD")
aged_coll <- c("#483D8B","#00BFFF")
ht1 <- SCP::GroupHeatmap(srt = subset(combine_species_new_label_S4,dataset%in%'human'),
                  features = genes_conserved_up,
                  #cell_annotation = 'aged',
                  group.by = c('annotation'),
                  group_palcolor  = list(celltypes_col),
                  split.by = c('aged'),
                  heatmap_palette = "YlOrRd",
                  cell_split_palcolor = aged_coll,
                  show_row_names = T,row_names_side = 'left',
                  features_label = F,
                  show_column_names = F)
ht1$plot
ht2 <- SCP::GroupHeatmap(srt = subset(combine_species_new_label_S4,dataset%in%'monkey'),
                         features = genes_conserved_up,
                         #cell_annotation = 'aged',
                         group.by = c('annotation'),
                         group_palcolor  = list(celltypes_col),
                         split.by = c('aged'),
                         heatmap_palette = "YlOrRd",
                         cell_split_palcolor = aged_coll,
                         show_row_names = T,row_names_side = 'left',
                         features_label = F,
                         show_column_names = F)

ht2$plot
ht3 <- SCP::GroupHeatmap(srt = subset(combine_species_new_label_S4,dataset%in%'mouse1'),
                         features = genes_conserved_up,
                         #cell_annotation = 'aged',
                         group.by = c('annotation'),
                         group_palcolor  = list(celltypes_col),
                         split.by = c('aged'),
                         heatmap_palette = "YlOrRd",
                         cell_split_palcolor = aged_coll,
                         show_row_names = T,row_names_side = 'left',
                         features_label = F,
                         show_column_names = F)

ht3$plot
ht1$plot+ht2$plot+ht3$plot

combine_species_new_label_S4@meta.data$group_species <- paste(combine_species_new_label_S4$dataset,
                                                              combine_species_new_label_S4$aged,sep = '_')
Idents(combine_species_new_label_S4) <- 'group_species'
combine_species_new_label_S4_sub <- subset(combine_species_new_label_S4,downsample=3000)
combine_species_new_label_S4_sub <- subset(combine_species_new_label_S4_sub,dataset%in%c('human','monkey','mouse1'))
table(combine_species_new_label_S4_sub$group_species)
combine_species_new_label_S4_sub$group_species <- factor(combine_species_new_label_S4_sub$group_species,
                                                         levels = c('human_Young','human_Aged',
                                                                    'monkey_Young','monkey_Aged',
                                                                    'mouse1_Young','mouse1_Aged' ))
write_rds(combine_species_new_label_S4_sub,file = hu("res", "combine_species_new_label_S4_sub.rds"))
combine_species_new_label_S4_sub <- readRDS(hu("res", "combine_species_new_label_S4_sub.rds"))
pdf(hu("fig", "V1", "fig1_feature.pdf"),width = 4,height = 5)
SCP::FeatureStatPlot(combine_species_new_label_S4_sub, 
                     stat.by = "ICAM1", 
                     group.by = "group_species", plot_type = "col")
SCP::FeatureStatPlot(combine_species_new_label_S4_sub, 
                     stat.by = "VCAM1", 
                     group.by = "group_species", plot_type = "col")
SCP::FeatureStatPlot(combine_species_new_label_S4_sub, 
                     stat.by = "HSPA1A", 
                     group.by = "group_species", plot_type = "col")
SCP::FeatureStatPlot(combine_species_new_label_S4_sub, 
                     stat.by = "NR4A1", 
                     group.by = "group_species", plot_type = "col")
SCP::FeatureStatPlot(combine_species_new_label_S4_sub, 
                     stat.by = "EPHB1", 
                     group.by = "group_species", plot_type = "col")
dev.off()


target_gene <- c('',,,)













