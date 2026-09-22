# =============================================================================
# 14_fig1_conserved_DEG_heatmaps.R
# Title : Conserved age-DEG heat maps (aged/young)
# Figure: Fig. 1E-F
# Module: 02_figure1_cross_species_atlas
# Description:
#   Cross-species heat maps of age-DEGs within aged and young subsets.
# Inputs:
#   - hu("res", "combine_species_new_label.rds")
# Outputs:
#   - hu("fig", "V1", "fig1_celltype_heatmap_conserved_Aged.pdf")
#   - hu("fig", "V1", "fig1_celltype_heatmap_conserved_young.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

library(VennDiagram)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(Seurat)
combine_species_new_label <- read_rds(hu("res", "combine_species_new_label.rds"))
combine_species_new_label_Aged <- subset(combine_species_new_label,aged%in%'Aged')
combine_species_new_label_Aged <- JoinLayers(combine_species_new_label_Aged)
Idents(combine_species_new_label_Aged) <- 'annotation'
celltype_human_Aged <- subset(combine_species_new_label_Aged,dataset%in%'human')
celltype_monkey_Aged <- subset(combine_species_new_label_Aged,dataset%in%'monkey')
celltype_mouse1_Aged <- subset(combine_species_new_label_Aged,dataset%in%'mouse1')
celltype_human_markers_Aged <- FindAllMarkers(celltype_human_Aged,min.pct = 0.2)
celltype_monkey_markers_Aged <- FindAllMarkers(celltype_monkey_Aged,min.pct = 0.2)
celltype_mouse1_markers_Aged <- FindAllMarkers(celltype_mouse1_Aged,min.pct = 0.2)

celltypes <- c('FB','EC','Myeloid',
               'T','B','Pericytes',
               'SMC','Neural')
lapply(celltypes, function(cl){
  veen <- list('human'=filter(celltype_human_markers_Aged,cluster%in%cl&avg_log2FC>2) %>% .$gene,
               'monkey'=filter(celltype_monkey_markers_Aged,cluster%in%cl&avg_log2FC>2) %>% .$gene,
               'mouse'=filter(celltype_mouse1_markers_Aged,cluster%in%cl&avg_log2FC>2) %>% .$gene)
  return(veen)
})->celltypes_Aged_venn
names(celltypes_Aged_venn) <- celltypes

FB_Aged_veen <- get.venn.partitions(celltypes_Aged_venn$FB)
EC_Aged_veen <- get.venn.partitions(celltypes_Aged_venn$EC)
Myeloid_Aged_veen <- get.venn.partitions(celltypes_Aged_venn$Myeloid)
T_Aged_veen <- get.venn.partitions(celltypes_Aged_venn$T)
B_Aged_veen <- get.venn.partitions(celltypes_Aged_venn$B)
Pericytes_Aged_veen <- get.venn.partitions(celltypes_Aged_venn$Pericytes)
SMC_Aged_veen <- get.venn.partitions(celltypes_Aged_venn$SMC)
Neural_Aged_veen <- get.venn.partitions(celltypes_Aged_venn$Neural)

celltypes_genes_conserved_Aged <- c(FB_Aged_veen[1,5][[1]],
                                    EC_Aged_veen[1,5][[1]],
                                    Myeloid_Aged_veen[1,5][[1]],
                                    T_Aged_veen[1,5][[1]],
                                    B_Aged_veen[1,5][[1]],
                                    Pericytes_Aged_veen[1,5][[1]],
                                    SMC_Aged_veen[1,5][[1]],
                                    Neural_Aged_veen[1,5][[1]]) %>% unique()
celltypes_genes_species_huamn_Aged <- c(FB_Aged_veen[7,5][[1]],
                                        EC_Aged_veen[7,5][[1]],
                                        Myeloid_Aged_veen[7,5][[1]],
                                        T_Aged_veen[7,5][[1]],
                                        B_Aged_veen[7,5][[1]],
                                        Pericytes_Aged_veen[7,5][[1]],
                                        SMC_Aged_veen[7,5][[1]],
                                        Neural_Aged_veen[7,5][[1]]) %>% unique()
celltypes_genes_species_monkey_Aged <- c(FB_Aged_veen[6,5][[1]],
                                         EC_Aged_veen[6,5][[1]],
                                         Myeloid_Aged_veen[6,5][[1]],
                                         T_Aged_veen[6,5][[1]],
                                         B_Aged_veen[6,5][[1]],
                                         Pericytes_Aged_veen[6,5][[1]],
                                         SMC_Aged_veen[6,5][[1]],
                                         Neural_Aged_veen[6,5][[1]]) %>% unique()
celltypes_genes_species_mouse_Aged <- c(FB_Aged_veen[4,5][[1]],
                                        EC_Aged_veen[4,5][[1]],
                                        Myeloid_Aged_veen[4,5][[1]],
                                        T_Aged_veen[4,5][[1]],
                                        B_Aged_veen[4,5][[1]],
                                        Pericytes_Aged_veen[4,5][[1]],
                                        SMC_Aged_veen[4,5][[1]],
                                        Neural_Aged_veen[4,5][[1]]) %>% unique()
c(length(celltypes_genes_conserved_Aged),
  length(celltypes_genes_species_huamn_Aged),
  length(celltypes_genes_species_monkey_Aged),
  length(celltypes_genes_species_mouse_Aged))
# create exp --------------------------------------------------------------
combine_species_new_label_Aged@meta.data$new_group <- paste(combine_species_new_label_Aged$species,
                                                            combine_species_new_label_Aged$annotation,sep = '_')
Idents(combine_species_new_label_Aged) <- 'new_group'
mat_corss_Aged <- AverageExpression(combine_species_new_label_Aged)
mat_corss_Aged$RNA[1:10,1:10]
mat_corss_scale_Aged <- mat_corss_Aged$RNA %>% as.matrix() %>% t()
mat_corss_scale_Aged <- mat_corss_scale_Aged[,c(celltypes_genes_conserved_Aged,
                                                celltypes_genes_species_huamn_Aged,
                                                celltypes_genes_species_monkey_Aged,
                                                celltypes_genes_species_mouse_Aged)]
mat_corss_scale_Aged <- mat_corss_scale_Aged[c(paste('human',celltypes,sep = '-'),
                                               paste('monkey',celltypes,sep = '-'),
                                               paste('mouse',celltypes,sep = '-')),]
# label ann ---------------------------------------------------------------
library(ComplexHeatmap)
library(colorRamp2)
group = rep(c("Conserved", "Species_human","Species_monkey","Species_mouse"),
            times = c(length(c(celltypes_genes_conserved_Aged)),
                      length(celltypes_genes_species_huamn_Aged),
                      length(celltypes_genes_species_monkey_Aged),
                      length(celltypes_genes_species_mouse_Aged)))
col1 = colorRamp2(c(0, 5), c( "white", "#FC8D62"))
col2 = colorRamp2(c(0, 5), c("white", "#8DD3C7"))
col3 = colorRamp2(c(0, 5), c("white", "#BC80BD"))
col4 = colorRamp2(c(0, 5), c("white", "pink"))
clusters <- rep(celltypes, 3)
clusters_col <- rep(c("#FF1493", "#FFB6C1", "#483D8B", "#00BFFF", "#228B22", "#7FFF00", "#FFFF00", "#FF4500"), 3)
names(clusters_col) <- clusters
clusters <- factor(clusters,levels = celltypes)
# plot heatmap ------------------------------------------------------------
left_annotation= rowAnnotation(clusters = clusters,
                               col = list(clusters = clusters_col))
#write_rds(mat_corss_scale,hu("res", "mat_corss_scale.rds"))
#mat_corss_scale <- read_rds(hu("res", "mat_corss_scale.rds"))
ht1 = Heatmap(mat_corss_scale_Aged[, group == "Conserved"], col = col1, name = "Conserved",
              left_annotation = left_annotation,row_dend_reorder = celltypes,
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              #row_order = c('FB','EC','Myeloid','T','B','Pericytes','SMC','Neural'),
              show_row_names = F,show_column_names = F)
ht2 = Heatmap(mat_corss_scale_Aged[, group == "Species_human"], col = col2, name = "Species_human",
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht3 = Heatmap(mat_corss_scale_Aged[, group == "Species_monkey"], col = col3, name = "Species_monkey",
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht4 = Heatmap(mat_corss_scale_Aged[, group == "Species_mouse"], col = col4, name = "Species_mouse",
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
pdf(hu("fig", "V1", "fig1_celltype_heatmap_conserved_Aged.pdf"),width = 10,height = 4)
combined_heatmap <- ht1 + ht2 + ht3 + ht4
# 绘制合并后的热图，设置热图之间的间距
draw(combined_heatmap, padding = unit(c(20, 20, 20, 20), "mm"))
dev.off()

# young -------------------------------------------------------------------
combine_species_new_label_young <- subset(combine_species_new_label,aged%in%'Young')
combine_species_new_label_young <- JoinLayers(combine_species_new_label_young)
Idents(combine_species_new_label_young) <- 'annotation'
celltype_human_young <- subset(combine_species_new_label_young,dataset%in%'human')
celltype_monkey_young <- subset(combine_species_new_label_young,dataset%in%'monkey')
celltype_mouse1_young <- subset(combine_species_new_label_young,dataset%in%'mouse1')
celltype_human_markers_young <- FindAllMarkers(celltype_human_young,min.pct = 0.2)
celltype_monkey_markers_young <- FindAllMarkers(celltype_monkey_young,min.pct = 0.2)
celltype_mouse1_markers_young <- FindAllMarkers(celltype_mouse1_young,min.pct = 0.2)

celltypes <- c('FB','EC','Myeloid',
               'T','B','Pericytes',
               'SMC','Neural')
lapply(celltypes, function(cl){
  veen <- list('human'=filter(celltype_human_markers_young,cluster%in%cl&avg_log2FC>2) %>% .$gene,
               'monkey'=filter(celltype_monkey_markers_young,cluster%in%cl&avg_log2FC>2) %>% .$gene,
               'mouse'=filter(celltype_mouse1_markers_young,cluster%in%cl&avg_log2FC>2) %>% .$gene)
  return(veen)
})->celltypes_venn_young
names(celltypes_venn_young) <- celltypes
FB_veen_young <- get.venn.partitions(celltypes_venn_young$FB)
EC_veen_young <- get.venn.partitions(celltypes_venn_young$EC)
Myeloid_veen_young <- get.venn.partitions(celltypes_venn_young$Myeloid)
T_veen_young <- get.venn.partitions(celltypes_venn_young$T)
B_veen_young <- get.venn.partitions(celltypes_venn_young$B)
Pericytes_veen_young <- get.venn.partitions(celltypes_venn_young$Pericytes)
SMC_veen_young <- get.venn.partitions(celltypes_venn_young$SMC)
Neural_veen_young <- get.venn.partitions(celltypes_venn_young$Neural)

celltypes_genes_conserved_young <- c(FB_veen_young[1,5][[1]],
                               EC_veen_young[1,5][[1]],
                               Myeloid_veen_young[1,5][[1]],
                               T_veen_young[1,5][[1]],
                               B_veen_young[1,5][[1]],
                               Pericytes_veen_young[1,5][[1]],
                               SMC_veen_young[1,5][[1]],
                               Neural_veen_young[1,5][[1]]) %>% unique()
celltypes_genes_species_huamn_young <- c(FB_veen_young[7,5][[1]],
                                   EC_veen_young[7,5][[1]],
                                   Myeloid_veen_young[7,5][[1]],
                                   T_veen_young[7,5][[1]],
                                   B_veen_young[7,5][[1]],
                                   Pericytes_veen_young[7,5][[1]],
                                   SMC_veen_young[7,5][[1]],
                                   Neural_veen_young[7,5][[1]]) %>% unique()
celltypes_genes_species_monkey_young <- c(FB_veen_young[6,5][[1]],
                                    EC_veen_young[6,5][[1]],
                                    Myeloid_veen_young[6,5][[1]],
                                    T_veen_young[6,5][[1]],
                                    B_veen_young[6,5][[1]],
                                    Pericytes_veen_young[6,5][[1]],
                                    SMC_veen_young[6,5][[1]],
                                    Neural_veen_young[6,5][[1]]) %>% unique()
celltypes_genes_species_mouse_young <- c(FB_veen_young[4,5][[1]],
                                   EC_veen_young[4,5][[1]],
                                   Myeloid_veen_young[4,5][[1]],
                                   T_veen_young[4,5][[1]],
                                   B_veen_young[4,5][[1]],
                                   Pericytes_veen_young[4,5][[1]],
                                   SMC_veen_young[4,5][[1]],
                                   Neural_veen_young[4,5][[1]]) %>% unique()
c(length(celltypes_genes_conserved),
  length(celltypes_genes_species_huamn),
  length(celltypes_genes_species_monkey),
  length(celltypes_genes_species_mouse))
# create exp --------------------------------------------------------------
combine_species_new_label_young@meta.data$new_group <- paste(combine_species_new_label_young$species,
                                                             combine_species_new_label_young$annotation,sep = '_')
Idents(combine_species_new_label_young) <- 'new_group'
mat_corss_young <- AverageExpression(combine_species_new_label_young)
mat_corss_young$RNA[1:10,1:10]
mat_corss_scale_young <- mat_corss_young$RNA %>% as.matrix() %>% t()
mat_corss_scale_young <- mat_corss_scale_young[,c(celltypes_genes_conserved_young,
                                      celltypes_genes_species_huamn_young,
                                      celltypes_genes_species_monkey_young,
                                      celltypes_genes_species_mouse_young)]
celltypes1 <- c('FB','EC','Myeloid','T','B','Pericytes','SMC','Neural')
mat_corss_scale_young <- mat_corss_scale_young[c(paste('human',celltypes1,sep = '-'),
                                     paste('monkey',celltypes1,sep = '-'),
                                     paste('mouse',celltypes1,sep = '-')),]
# label ann ---------------------------------------------------------------
group = rep(c("Conserved", "Species_human","Species_monkey","Species_mouse"),
            times = c(length(c(celltypes_genes_conserved_young)),
                      length(celltypes_genes_species_huamn_young),
                      length(celltypes_genes_species_monkey_young),
                      length(celltypes_genes_species_mouse_young)))
col1 = colorRamp2(c(0, 5), c( "white", "#FC8D62"))
col2 = colorRamp2(c(0, 5), c("white", "#8DD3C7"))
col3 = colorRamp2(c(0, 5), c("white", "#BC80BD"))
col4 = colorRamp2(c(0, 5), c("white", "pink"))
clusters <- rep(celltypes, 3)
clusters_col <- rep(c("#FF1493", "#FFB6C1", "#483D8B", "#00BFFF", "#228B22", "#7FFF00", "#FFFF00", "#FF4500"), 3)
names(clusters_col) <- clusters
clusters <- factor(clusters,levels = celltypes)
# plot heatmap ------------------------------------------------------------
left_annotation= rowAnnotation(clusters = clusters,
                               col = list(clusters = clusters_col))
ht1 = Heatmap(mat_corss_scale_young[, group == "Conserved"], col = col1, name = "Conserved",
              left_annotation = left_annotation,
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              #row_order = c('FB','EC','Myeloid','T','B','Pericytes','SMC','Neural'),
              show_row_names = F,show_column_names = F)
ht2 = Heatmap(mat_corss_scale_young[, group == "Species_human"], col = col2, name = "Species_human",
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht3 = Heatmap(mat_corss_scale_young[, group == "Species_monkey"], col = col3, name = "Species_monkey",
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht4 = Heatmap(mat_corss_scale_young[, group == "Species_mouse"], col = col4, name = "Species_mouse",
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
pdf(hu("fig", "V1", "fig1_celltype_heatmap_conserved_young.pdf"),width = 10,height = 4)
combined_heatmap <- ht1 + ht2 + ht3 + ht4
# 绘制合并后的热图，设置热图之间的间距
draw(combined_heatmap, padding = unit(c(20, 20, 20, 20), "mm"))
dev.off()

