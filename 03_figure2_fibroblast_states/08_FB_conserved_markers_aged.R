# =============================================================================
# 08_FB_conserved_markers_aged.R
# Title : Conserved fibroblast markers (aged)
# Figure: Supplementary
# Module: 03_figure2_fibroblast_states
# Description:
#   Per-species marker intersection within aged fibroblasts.
# Inputs:
#   - hu("res", "FB_species_rpca_sub_label.rds")
# Outputs:
#   - hu("res", "heatmap_conserved_aged.pdf")
#   - hu("res", "heatmap_conserved_Young.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

library(VennDiagram)
library(RColorBrewer)
library(Seurat)
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
FB_species_rpca_sub_label <- read_rds(hu("res", "FB_species_rpca_sub_label.rds"))
FB_aged <- subset(FB_species_rpca_sub_label,aged%in%'Aged')
Idents(FB_aged) <- 'subtype'
FB_human_aged <- subset(FB_aged,species%in%'human')
FB_monkey_aged <- subset(FB_aged,species%in%'monkey')
FB_mouse1_aged <- subset(FB_aged,species%in%'mouse1')
FB_mouse2_aged <- subset(FB_aged,species%in%'mouse2')
FB_human_markers_aged <- FindAllMarkers(FB_human_aged,min.pct = 0.3)
FB_monkey_markers_aged <- FindAllMarkers(FB_monkey_aged,min.pct = 0.3)
FB_mouse1_markers_aged <- FindAllMarkers(FB_mouse1_aged,min.pct = 0.3)
FB_mouse2_markers_aged <- FindAllMarkers(FB_mouse2_aged,min.pct = 0.3)

clusters <- c('FB1','FB2','FB3','FB4','FB5')
lapply(clusters, function(cl){
  veen <- list('human'=filter(FB_human_markers_aged,cluster%in%cl&avg_log2FC>0.8) %>% .$gene,
               'monkey'=filter(FB_monkey_markers_aged,cluster%in%cl&avg_log2FC>0.8) %>% .$gene,
               'mouse'=filter(FB_mouse1_markers_aged,cluster%in%cl&avg_log2FC>0.8) %>% .$gene)
  return(veen)
})->clusters_venn_aged
names(clusters_venn_aged) <- clusters

FB1_veen_aged <- get.venn.partitions(clusters_venn_aged$FB1)
#FB2_veen_aged <- get.venn.partitions(clusters_venn_aged$FB2)
FB3_veen_aged <- get.venn.partitions(clusters_venn_aged$FB3)
FB4_veen_aged <- get.venn.partitions(clusters_venn_aged$FB4)
FB5_veen_aged <- get.venn.partitions(clusters_venn_aged$FB5)
FB_genes_conserved <- c(FB1_veen_aged[1,5][[1]],
                        #FB2_veen_aged[1,5][[1]],
                        FB3_veen_aged[1,5][[1]],
                        FB4_veen_aged[1,5][[1]],
                        FB5_veen_aged[1,5][[1]]) %>% unique()
FB_genes_species_huamn <- c(FB1_veen_aged[7,5][[1]],
                           # FB2_veen_aged[7,5][[1]],
                            FB3_veen_aged[7,5][[1]],
                            FB4_veen_aged[7,5][[1]],
                            FB5_veen_aged[7,5][[1]])%>% unique()
FB_genes_species_monkey <- c(FB1_veen_aged[6,5][[1]],
                            # FB2_veen_aged[6,5][[1]],
                             FB3_veen_aged[6,5][[1]],
                             FB4_veen_aged[6,5][[1]],
                             FB5_veen_aged[6,5][[1]])%>% unique()
FB_genes_species_mouse <- c(FB1_veen_aged[4,5][[1]],
                           # FB2_veen_aged[4,5][[1]],
                            FB3_veen_aged[4,5][[1]],
                            FB4_veen_aged[4,5][[1]],
                            FB5_veen_aged[4,5][[1]])%>% unique()
FB_aged@meta.data$new_group <- paste(FB_aged$species,FB_aged$subtype,sep = '_')
Idents(FB_aged) <- 'new_group'
mat_corss <- AverageExpression(FB_aged)
mat_corss$RNA[1:10,1:10]
# label ann ---------------------------------------------------------------
group = rep(c("Conserved", "Species_human","Species_monkey","Species_mouse"),
            times = c(length(c(FB_genes_conserved)),
                      length(FB_genes_species_huamn),
                      length(FB_genes_species_monkey),
                      length(FB_genes_species_mouse)))
col1 = colorRamp2(c(0, 10), c( "white", "#FC8D62"))
col2 = colorRamp2(c(0, 10), c("white", "#8DD3C7"))
col3 = colorRamp2(c(0, 10), c("white", "#BC80BD"))
col4 = colorRamp2(c(0, 10), c("white", "pink"))
clusters <- rep(c('FB1','FB2','FB3','FB4','FB5'), 3)
clusters_col <- rep(c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728",'#FFED6F'), 3)
names(clusters_col) <- clusters
left_annotation= rowAnnotation(clusters = clusters,
                               col = list(clusters = clusters_col))
dim(mat_corss_scale[, group == "Conserved"])
library(ComplexHeatmap)

mat_corss_scale <- mat_corss$RNA %>% as.matrix() %>% t()
mat_corss_scale <- mat_corss_scale[,c(FB_genes_conserved,
                                      FB_genes_species_huamn,
                                      FB_genes_species_monkey,
                                      FB_genes_species_mouse)]
dim(mat_corss_scale)
length(group)
mat_corss_scale <- mat_corss_scale[1:15,]
ht1 = Heatmap(mat_corss_scale[, group == "Conserved"], col = col1, name = "Conserved",
              left_annotation = left_annotation,
              row_split = c(rep('human',5),rep('monkey',5),rep('mouse',5)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht2 = Heatmap(mat_corss_scale[, group == "Species_human"], col = col2, name = "Species_human",
              row_split = c(rep('human',5),rep('monkey',5),rep('mouse',5)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht3 = Heatmap(mat_corss_scale[, group == "Species_monkey"], col = col3, name = "Species_monkey",
              row_split = c(rep('human',5),rep('monkey',5),rep('mouse',5)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht4 = Heatmap(mat_corss_scale[, group == "Species_mouse"], col = col4, name = "Species_mouse",
              row_split = c(rep('human',5),rep('monkey',5),rep('mouse',5)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
pdf(hu("res", "heatmap_conserved_aged.pdf"),width = 8,height = 4)
ht1 + ht2 + ht3 + ht4
dev.off()

FB_Young <- subset(FB_species_rpca_sub_label,aged%in%'Young')
Idents(FB_Young) <- 'subtype'
FB_human_Young <- subset(FB_Young,species%in%'human')
FB_monkey_Young <- subset(FB_Young,species%in%'monkey')
FB_mouse1_Young <- subset(FB_Young,species%in%'mouse1')
FB_mouse2_Young <- subset(FB_Young,species%in%'mouse2')
FB_human_markers_Young <- FindAllMarkers(FB_human_Young,min.pct = 0.3)
FB_monkey_markers_Young <- FindAllMarkers(FB_monkey_Young,min.pct = 0.3)
FB_mouse1_markers_Young <- FindAllMarkers(FB_mouse1_Young,min.pct = 0.3)
FB_mouse2_markers_Young <- FindAllMarkers(FB_mouse2_Young,min.pct = 0.3)

clusters <- c('FB1','FB2','FB3','FB4','FB5')
lapply(clusters, function(cl){
  veen <- list('human'=filter(FB_human_markers_Young,cluster%in%cl&avg_log2FC>0.8) %>% .$gene,
               'monkey'=filter(FB_monkey_markers_Young,cluster%in%cl&avg_log2FC>0.8) %>% .$gene,
               'mouse'=filter(FB_mouse1_markers_Young,cluster%in%cl&avg_log2FC>0.8) %>% .$gene)
  return(veen)
})->clusters_venn_Young
names(clusters_venn_Young) <- clusters

FB1_veen_Young <- get.venn.partitions(clusters_venn_Young$FB1)
#FB2_veen_Young <- get.venn.partitions(clusters_venn_Young$FB2)
FB3_veen_Young <- get.venn.partitions(clusters_venn_Young$FB3)
FB4_veen_Young <- get.venn.partitions(clusters_venn_Young$FB4)
FB5_veen_Young <- get.venn.partitions(clusters_venn_Young$FB5)
FB_genes_conserved <- c(FB1_veen_Young[1,5][[1]],
                        #FB2_veen_Young[1,5][[1]],
                        FB3_veen_Young[1,5][[1]],
                        FB4_veen_Young[1,5][[1]],
                        FB5_veen_Young[1,5][[1]]) %>% unique()
FB_genes_species_huamn <- c(FB1_veen_Young[7,5][[1]],
                            #FB2_veen_Young[7,5][[1]],
                            FB3_veen_Young[7,5][[1]],
                            FB4_veen_Young[7,5][[1]],
                            FB5_veen_Young[7,5][[1]])%>% unique()
FB_genes_species_monkey <- c(FB1_veen_Young[6,5][[1]],
                             #FB2_veen_Young[6,5][[1]],
                             FB3_veen_Young[6,5][[1]],
                             FB4_veen_Young[6,5][[1]],
                             FB5_veen_Young[6,5][[1]])%>% unique()
FB_genes_species_mouse <- c(FB1_veen_Young[4,5][[1]],
                            #FB2_veen_Young[4,5][[1]],
                            FB3_veen_Young[4,5][[1]],
                            FB4_veen_Young[4,5][[1]],
                            FB5_veen_Young[4,5][[1]])%>% unique()
FB_aged@meta.data$new_group <- paste(FB_aged$species,
                                     FB_aged$subtype,sep = '_')
Idents(FB_aged) <- 'new_group'
mat_corss <- AverageExpression(FB_aged)
mat_corss$RNA[1:10,1:10]
# label ann ---------------------------------------------------------------
group = rep(c("Conserved", "Species_human","Species_monkey","Species_mouse"),
            times = c(length(c(FB_genes_conserved)),
                      length(FB_genes_species_huamn),
                      length(FB_genes_species_monkey),
                      length(FB_genes_species_mouse)))
col1 = colorRamp2(c(0, 10), c( "white", "#FC8D62"))
col2 = colorRamp2(c(0, 10), c("white", "#8DD3C7"))
col3 = colorRamp2(c(0, 10), c("white", "#BC80BD"))
col4 = colorRamp2(c(0, 10), c("white", "pink"))
clusters <- rep(c('FB1','FB2','FB3','FB4','FB5'), 3)
clusters_col <- rep(c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728",'#FFED6F'), 3)
names(clusters_col) <- clusters
left_annotation= rowAnnotation(clusters = clusters,
                               col = list(clusters = clusters_col))
dim(mat_corss_scale[, group == "Conserved"])
library(ComplexHeatmap)

mat_corss_scale <- mat_corss$RNA %>% as.matrix() %>% t()
mat_corss_scale <- mat_corss_scale[,c(FB_genes_conserved,
                                      FB_genes_species_huamn,
                                      FB_genes_species_monkey,
                                      FB_genes_species_mouse)]
dim(mat_corss_scale)
length(group)
mat_corss_scale <- mat_corss_scale[1:15,]
ht1 = Heatmap(mat_corss_scale[, group == "Conserved"], col = col1, name = "Conserved",
              left_annotation = left_annotation,
              row_split = c(rep('human',5),rep('monkey',5),rep('mouse',5)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht2 = Heatmap(mat_corss_scale[, group == "Species_human"], col = col2, name = "Species_human",
              row_split = c(rep('human',5),rep('monkey',5),rep('mouse',5)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht3 = Heatmap(mat_corss_scale[, group == "Species_monkey"], col = col3, name = "Species_monkey",
              row_split = c(rep('human',5),rep('monkey',5),rep('mouse',5)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht4 = Heatmap(mat_corss_scale[, group == "Species_mouse"], col = col4, name = "Species_mouse",
              row_split = c(rep('human',5),rep('monkey',5),rep('mouse',5)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
pdf(hu("res", "heatmap_conserved_Young.pdf"),width = 8,height = 4)
ht1 + ht2 + ht3 + ht4
dev.off()

