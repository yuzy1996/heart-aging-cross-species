# =============================================================================
# 18_celltype_conserved_markers_young.R
# Title : Conserved cell-type markers (young)
# Figure: Supplementary
# Module: 02_figure1_cross_species_atlas
# Description:
#   Per-species marker detection and intersection within young hearts.
# Inputs:
#   - hu("res", "combine_species_new_label.rds")
# Outputs:
#   - hu("res", "celltype_heatmap_conserved.pdf")
#   - hu("res", "celltype_heatmap_conserved1.pdf")
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
combine_species_new_label_young <- subset(combine_species_new_label,aged%in%'Young')
combine_species_new_label_young <- JoinLayers(combine_species_new_label_young)
Idents(combine_species_new_label_young) <- 'annotation'
celltype_young_marker <- FindAllMarkers(combine_species_new_label_young,min.pct = 0.2)
celltype_young_marker_cut1 <- celltype_young_marker %>% dplyr::filter(avg_log2FC>1&pct.1>0.3)

Idents(combine_species_new_label_young) <- 'annotation'
celltype_human <- subset(combine_species_new_label_young,species%in%'human')
celltype_monkey <- subset(combine_species_new_label_young,species%in%'monkey')
celltype_mouse <- subset(combine_species_new_label_young,species%in%'mouse')
celltype_human_markers <- FindAllMarkers(celltype_human,min.pct = 0.2)
celltype_monkey_markers <- FindAllMarkers(celltype_monkey,min.pct = 0.2)
celltype_mouse_markers <- FindAllMarkers(celltype_mouse,min.pct = 0.2)
celltype_mouse3_markers <- FindAllMarkers(celltype_mouse3,min.pct = 0.2)

celltypes <- c('FB','EC','Myeloid',
               'T','B','Pericytes',
               'SMC','Neural')
lapply(celltypes, function(cl){
  veen <- list('human'=filter(celltype_human_markers,cluster%in%cl&avg_log2FC>2) %>% .$gene,
               'monkey'=filter(celltype_monkey_markers,cluster%in%cl&avg_log2FC>2) %>% .$gene,
               'mouse'=filter(celltype_mouse_markers,cluster%in%cl&avg_log2FC>2) %>% .$gene)
  return(veen)
})->celltypes_venn1
names(celltypes_venn1) <- celltypes
celltypes <- c('FB','EC','Myeloid',
               'T','B','Pericytes',
               'SMC','Neural')
lapply(celltypes, function(cl){
  veen <- list('human'=filter(celltype_human_markers,cluster%in%cl&avg_log2FC>2) %>% .$gene,
               'monkey'=filter(celltype_monkey_markers,cluster%in%cl&avg_log2FC>2) %>% .$gene,
               'mouse'=filter(celltype_mouse3_markers,cluster%in%cl&avg_log2FC>2) %>% .$gene)
  return(veen)
})->celltypes_venn2
names(celltypes_venn2) <- celltypes

celltypes_venn=celltypes_venn2
FB_veen <- get.venn.partitions(celltypes_venn$FB)
EC_veen <- get.venn.partitions(celltypes_venn$EC)
Myeloid_veen <- get.venn.partitions(celltypes_venn$Myeloid)
T_veen <- get.venn.partitions(celltypes_venn$T)
B_veen <- get.venn.partitions(celltypes_venn$B)
Pericytes_veen <- get.venn.partitions(celltypes_venn$Pericytes)
SMC_veen <- get.venn.partitions(celltypes_venn$SMC)
Neural_veen <- get.venn.partitions(celltypes_venn$Neural)

celltypes_genes_conserved <- c(FB_veen[1,5][[1]],
                               EC_veen[1,5][[1]],
                               Myeloid_veen[1,5][[1]],
                               T_veen[1,5][[1]],
                               B_veen[1,5][[1]],
                               Pericytes_veen[1,5][[1]],
                               SMC_veen[1,5][[1]],
                               Neural_veen[1,5][[1]]) %>% unique()
celltypes_genes_species_huamn <- c(FB_veen[7,5][[1]],
                                   EC_veen[7,5][[1]],
                                   Myeloid_veen[7,5][[1]],
                                   T_veen[7,5][[1]],
                                   B_veen[7,5][[1]],
                                   Pericytes_veen[7,5][[1]],
                                   SMC_veen[7,5][[1]],
                                   Neural_veen[7,5][[1]]) %>% unique()
celltypes_genes_species_monkey <- c(FB_veen[6,5][[1]],
                                    EC_veen[6,5][[1]],
                                    Myeloid_veen[6,5][[1]],
                                    T_veen[6,5][[1]],
                                    B_veen[6,5][[1]],
                                    Pericytes_veen[6,5][[1]],
                                    SMC_veen[6,5][[1]],
                                    Neural_veen[6,5][[1]]) %>% unique()
celltypes_genes_species_mouse <- c(FB_veen[4,5][[1]],
                                   EC_veen[4,5][[1]],
                                   Myeloid_veen[4,5][[1]],
                                   T_veen[4,5][[1]],
                                   B_veen[4,5][[1]],
                                   Pericytes_veen[4,5][[1]],
                                   SMC_veen[4,5][[1]],
                                   Neural_veen[4,5][[1]]) %>% unique()
c(length(celltypes_genes_conserved),
  length(celltypes_genes_species_huamn),
  length(celltypes_genes_species_monkey),
  length(celltypes_genes_species_mouse))
# create exp --------------------------------------------------------------
combine_species_new_label_young@meta.data$new_group <- paste(combine_species_new_label_young$species,
                                                             combine_species_new_label_young$annotation,sep = '_')
Idents(combine_species_new_label_young) <- 'new_group'
mat_corss <- AverageExpression(combine_species_new_label_young)
mat_corss$RNA[1:10,1:10]
mat_corss_scale <- mat_corss$RNA %>% as.matrix() %>% t()
mat_corss_scale <- mat_corss_scale[,c(celltypes_genes_conserved,
                                      celltypes_genes_species_huamn,
                                      celltypes_genes_species_monkey,
                                      celltypes_genes_species_mouse)]
celltypes1 <- c('FB','EC','Myeloid','T','B','Pericytes','SMC','Neural')
mat_corss_scale <- mat_corss_scale[c(paste('human',celltypes1,sep = '-'),
                                     paste('monkey',celltypes1,sep = '-'),
                                     paste('mouse1',celltypes1,sep = '-')),]
# label ann ---------------------------------------------------------------
group = rep(c("Conserved", "Species_human","Species_monkey","Species_mouse"),
            times = c(length(c(celltypes_genes_conserved)),
                      length(celltypes_genes_species_huamn),
                      length(celltypes_genes_species_monkey),
                      length(celltypes_genes_species_mouse)))
col1 = colorRamp2(c(0, 5), c( "white", "#FC8D62"))
col2 = colorRamp2(c(0, 5), c("white", "#8DD3C7"))
col3 = colorRamp2(c(0, 5), c("white", "#BC80BD"))
col4 = colorRamp2(c(0, 5), c("white", "pink"))
clusters <- rep(c('FB','EC','Myeloid',
                  'T','B','Pericytes',
                  'SMC','Neural'), 3)
clusters_col <- rep(c("#66A61E", "#E6AB02", "#A6761D", "#666666", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072"), 3)
names(clusters_col) <- clusters
library(ComplexHeatmap)
# plot heatmap ------------------------------------------------------------
left_annotation= rowAnnotation(clusters = clusters,
                               col = list(clusters = clusters_col))
dim(mat_corss_scale[, group == "Conserved"])


#write_rds(mat_corss_scale,hu("res", "mat_corss_scale.rds"))
#mat_corss_scale <- read_rds(hu("res", "mat_corss_scale.rds"))
ht1 = Heatmap(mat_corss_scale[, group == "Conserved"], col = col1, name = "Conserved",
              left_annotation = left_annotation,
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              #row_order = c('FB','EC','Myeloid','T','B','Pericytes','SMC','Neural'),
              show_row_names = F,show_column_names = F)
ht2 = Heatmap(mat_corss_scale[, group == "Species_human"], col = col2, name = "Species_human",
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht3 = Heatmap(mat_corss_scale[, group == "Species_monkey"], col = col3, name = "Species_monkey",
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht4 = Heatmap(mat_corss_scale[, group == "Species_mouse"], col = col4, name = "Species_mouse",
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
pdf(hu("res", "celltype_heatmap_conserved.pdf"),width = 10,height = 4)
combined_heatmap <- ht1 + ht2 + ht3 + ht4
# 绘制合并后的热图，设置热图之间的间距
draw(combined_heatmap, padding = unit(c(20, 20, 20, 20), "mm"))
dev.off()


pdf(hu("res", "celltype_heatmap_conserved1.pdf"),width = 22,height = 4)
Heatmap(mat_corss_scale[, group == "Conserved"], col = col1, name = "Conserved",
        left_annotation = left_annotation,
        row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
        cluster_rows = F,cluster_columns = F,
        #row_order = c('FB','EC','Myeloid','T','B','Pericytes','SMC','Neural'),
        show_row_names = F,show_column_names = T)
dev.off()

Myeloid_species_final@meta.data$new_group <- paste(Myeloid_species_final$species,
                                                   Myeloid_species_final$aged,
                                                   Myeloid_species_final$subtype,sep = '_')
Myeloid_species_final@meta.data$new_group  <- as.factor(Myeloid_species_final@meta.data$new_group)
Idents(Myeloid_species_final) <- 'new_group'
Myeloid_species_final_exp <- AverageExpression(Myeloid_species_final,assays = 'RNA')

mat = Myeloid_species_final_exp$RNA %>% as.matrix() %>% t()
mat <- mat[1:36,]
genes_conserved_up <- c(INF_MC_veen_up[1,5][[1]],
                        MHClow_MC_veen_up[1,5][[1]],
                        CCR2_MC_veen_up[1,5][[1]],
                        MHChi_MC_veen_up[1,5][[1]],
                        Mono_veen_up[1,5][[1]],
                        Timd4_resident_MC_veen_up[1,5][[1]])%>% unique()
genes_conserved_down <-c(INF_MC_veen_down[1,5][[1]],
                         MHClow_MC_veen_down[1,5][[1]],
                         CCR2_MC_veen_down[1,5][[1]],
                         MHChi_MC_veen_down[1,5][[1]],
                         Mono_veen_down[1,5][[1]],
                         Timd4_resident_MC_veen_down[1,5][[1]])%>% unique()
genes_species_huamn <- c(INF_MC_veen_up[7,5][[1]],
                         MHClow_MC_veen_up[7,5][[1]],
                         CCR2_MC_veen_up[7,5][[1]],
                         MHChi_MC_veen_up[7,5][[1]],
                         Mono_veen_up[7,5][[1]],
                         Timd4_resident_MC_veen_up[7,5][[1]],
                         INF_MC_veen_down[7,5][[1]],
                         MHClow_MC_veen_down[7,5][[1]],
                         CCR2_MC_veen_down[7,5][[1]],
                         MHChi_MC_veen_down[7,5][[1]],
                         Mono_veen_down[7,5][[1]],
                         Timd4_resident_MC_veen_down[7,5][[1]]) %>% unique()
genes_species_monkey <- c(INF_MC_veen_up[6,5][[1]],
                          MHClow_MC_veen_up[6,5][[1]],
                          CCR2_MC_veen_up[6,5][[1]],
                          MHChi_MC_veen_up[6,5][[1]],
                          Mono_veen_up[6,5][[1]],
                          Timd4_resident_MC_veen_up[6,5][[1]],
                          INF_MC_veen_down[6,5][[1]],
                          MHClow_MC_veen_down[6,5][[1]],
                          CCR2_MC_veen_down[6,5][[1]],
                          MHChi_MC_veen_down[6,5][[1]],
                          Mono_veen_down[6,5][[1]],
                          Timd4_resident_MC_veen_down[6,5][[1]])%>% unique()
genes_species_mouse <- c(INF_MC_veen_up[4,5][[1]],
                         MHClow_MC_veen_up[4,5][[1]],
                         CCR2_MC_veen_up[4,5][[1]],
                         MHChi_MC_veen_up[4,5][[1]],
                         Mono_veen_up[4,5][[1]],
                         Timd4_resident_MC_veen_up[4,5][[1]],
                         INF_MC_veen_down[4,5][[1]],
                         MHClow_MC_veen_down[4,5][[1]],
                         CCR2_MC_veen_down[4,5][[1]],
                         MHChi_MC_veen_down[4,5][[1]],
                         Mono_veen_down[4,5][[1]],
                         Timd4_resident_MC_veen_down[4,5][[1]]) %>% unique()
mat_corss <- mat[,unique(c(genes_conserved,genes_species_huamn,genes_species_monkey,genes_species_mouse))]
dim(mat_corss)
group = rep(c("Conserved", "Species_human","Species_monkey","Species_mouse"),
            times = c(length(c(genes_conserved_up,genes_conserved_down)),
                      length(genes_species_huamn),
                      length(genes_species_monkey),
                      length(genes_species_mouse)))
col1 = colorRamp2(c(-2, 2), c( "white", "red"))
col2 = colorRamp2(c(-2, 2), c("white", "blue"))
col3 = colorRamp2(c(-2, 2), c("white", "green"))
col4 = colorRamp2(c(-2, 2), c("white", "pink"))
clusters <- rep(c('CCR2_MC','INF_MC','MHChi_MC','MHClow_MC','Mono',"Timd4-resident-MC"), 6)
aged <- rep(c(rep('Young',6),rep('Aged',6)),3)
clusters_col <- rep(c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD","#66C2A5"), 6)
aged_col <- rep(c(rep('#7570B3',6),rep('#FFED6F',6)),3)
names(clusters_col) <- clusters
names(aged_col) <- aged
left_annotation= rowAnnotation(clusters = clusters,
                               aged=aged, 
                               col = list(clusters = clusters_col,
                                          aged=aged_col))

Heatmap(mat_corss_scale[, 1], col = col1, name = "Conserved",
        left_annotation = left_annotation,
        row_split = c(rep('human',12),rep('monkey',12),rep('mouse',12)),
        cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = F)
mat_corss_scale <- scale(mat_corss)
ht1 = Heatmap(mat_corss_scale[, group == "Conserved"], col = col1, name = "Conserved",
              left_annotation = left_annotation,
              row_order = row_order,
              row_split = c(rep('human',12),rep('monkey',12),rep('mouse',12)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht1
ht2 = Heatmap(mat_corss_scale[, group == "Species_human"], col = col2, name = "Species_human",
              row_split = c(rep('human',12),rep('monkey',12),rep('mouse',12)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht3 = Heatmap(mat_corss_scale[, group == "Species_monkey"], col = col3, name = "Species_monkey",
              row_split = c(rep('human',12),rep('monkey',12),rep('mouse',12)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht4 = Heatmap(mat_corss_scale[, group == "Species_mouse"], col = col4, name = "Species_mouse",
              row_split = c(rep('human',12),rep('monkey',12),rep('mouse',12)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht1 + ht2 + ht3 + ht4



combine_species_new_label_aged <- subset(combine_species_new_label,aged%in%'Aged')
