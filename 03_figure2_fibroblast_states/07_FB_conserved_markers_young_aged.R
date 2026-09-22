# =============================================================================
# 07_FB_conserved_markers_young_aged.R
# Title : Conserved fibroblast subtype markers
# Figure: Fig. 2B, 2E
# Module: 03_figure2_fibroblast_states
# Description:
#   FindAllMarkers within young/aged fibroblasts and intersect conserved subtype markers.
# Inputs:
#   - (see script)
# Outputs:
#   - hu("res", "heatmap_conserved.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

library(VennDiagram)
library(RColorBrewer)
library(circlize)
FB_young <- subset(FB_species_rpca_sub_label,aged%in%'Young')
FB_aged <- subset(FB_species_rpca_sub_label,aged%in%'Aged')

FB_young <- JoinLayers(FB_young)
Idents(FB_young) <- 'subtype'
FB_young_marker <- FindAllMarkers(FB_young,min.pct = 0.2)
FB_young_marker_cut1 <- FB_young_marker %>% dplyr::filter(avg_log2FC>1&pct.1>0.3)

Idents(FB_species_rpca_sub_label) <- 'subtype'
FB_human <- subset(FB_species_rpca_sub_label,species%in%'human')
FB_monkey <- subset(FB_species_rpca_sub_label,species%in%'monkey')
FB_mouse1 <- subset(FB_species_rpca_sub_label,species%in%'mouse1')
FB_mouse2 <- subset(FB_species_rpca_sub_label,species%in%'mouse2')
FB_human_markers <- FindAllMarkers(FB_human,min.pct = 0.2)
FB_monkey_markers <- FindAllMarkers(FB_monkey,min.pct = 0.2)
FB_mouse1_markers <- FindAllMarkers(FB_mouse1,min.pct = 0.2)
FB_mouse2_markers <- FindAllMarkers(FB_mouse2,min.pct = 0.2)

clusters <- c('FB1','FB2','FB3','FB4','FB5')
lapply(clusters, function(cl){
  veen <- list('human'=filter(FB_human_markers,cluster%in%cl&avg_log2FC>0.6) %>% .$gene,
               'monkey'=filter(FB_monkey_markers,cluster%in%cl&avg_log2FC>0.6) %>% .$gene,
               'mouse'=filter(FB_mouse1_markers,cluster%in%cl&avg_log2FC>0.6) %>% .$gene)
  return(veen)
})->clusters_venn
names(clusters_venn) <- clusters

FB1_veen <- get.venn.partitions(clusters_venn$FB1)
FB2_veen <- get.venn.partitions(clusters_venn$FB2)
FB3_veen <- get.venn.partitions(clusters_venn$FB3)
FB4_veen <- get.venn.partitions(clusters_venn$FB4)
FB5_veen <- get.venn.partitions(clusters_venn$FB5)
FB_genes_conserved <- c(FB1_veen[1,5][[1]],
                        FB2_veen[1,5][[1]],
                        FB3_veen[1,5][[1]],
                        FB4_veen[1,5][[1]],
                        FB5_veen[1,5][[1]]) %>% unique()
FB_genes_species_huamn <- c(FB1_veen[7,5][[1]],
                        FB2_veen[7,5][[1]],
                        FB3_veen[7,5][[1]],
                        FB4_veen[7,5][[1]],
                        FB5_veen[7,5][[1]])%>% unique()
FB_genes_species_monkey <- c(FB1_veen[6,5][[1]],
                        FB2_veen[6,5][[1]],
                        FB3_veen[6,5][[1]],
                        FB4_veen[6,5][[1]],
                        FB5_veen[6,5][[1]])%>% unique()
FB_genes_species_mouse <- c(FB1_veen[4,5][[1]],
                        FB2_veen[4,5][[1]],
                        FB3_veen[4,5][[1]],
                        FB4_veen[4,5][[1]],
                        FB5_veen[4,5][[1]])%>% unique()
FB_species_rpca_sub_label@meta.data$new_group <- paste(FB_species_rpca_sub_label$species,
                                                       FB_species_rpca_sub_label$subtype,sep = '_')
Idents(FB_species_rpca_sub_label) <- 'new_group'
mat_corss <- AverageExpression(FB_species_rpca_sub_label)
mat_corss$RNA[1:10,1:10]
# label ann ---------------------------------------------------------------
group = rep(c("Conserved", "Species_human","Species_monkey","Species_mouse"),
            times = c(length(c(FB_genes_conserved)),
                      length(FB_genes_species_huamn),
                      length(FB_genes_species_monkey),
                      length(FB_genes_species_mouse)))
col1 = colorRamp2(c(0, 3), c( "white", "#FC8D62"))
col2 = colorRamp2(c(0, 3), c("white", "#8DD3C7"))
col3 = colorRamp2(c(0, 3), c("white", "#BC80BD"))
col4 = colorRamp2(c(0, 3), c("white", "pink"))
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
pdf(hu("res", "heatmap_conserved.pdf"),width = 8,height = 4)
ht1 + ht2 + ht3 + ht4
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
