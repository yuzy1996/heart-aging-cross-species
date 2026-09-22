# =============================================================================
# 12_subcluster_DEG_venn_by_cluster.R
# Title : Subtype DEG Venn intersections by cluster
# Figure: Supplementary
# Module: 03_figure2_fibroblast_states
# Description:
#   Cluster-wise Venn intersection of age-DEGs across species.
# Inputs:
#   - (see script)
# Outputs:
#   - hu("res", "p_INF_MC_veen_up.pdf")
#   - hu("res", "p_INF_MC_veen_down.pdf")
#   - hu("res", "p_MHClow_MC_veen_up.pdf")
#   - hu("res", "p_CCR2_MC_veen_up.pdf")
#   - hu("res", "p_CCR2_MC_veen_down.pdf")
#   - hu("res", "p_MHChi_MC_veen_up.pdf")
#   - hu("res", "p_MHChi_MC_veen_down.pdf")
#   - hu("res", "p_Mono_veen_up.pdf")
#   - hu("res", "p_Mono_veen_down.pdf")
#   - hu("res", "p_Timd4_resident_MC_veen_up.pdf")
#   - hu("res", "p_Timd4_resident_MC_veen_down.pdf")
#   - hu("res", "p_FB4_veen_up1.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

run_cluster_diff <- function(obj){
  Idents(obj) <- 'aged'
  Clusters <- data.frame(table(obj$subtype)) %>% filter(Freq>10) %>% .$Var1
  Seurat::DefaultAssay(obj) <- 'RNA'
  marker_list <- lapply(Clusters, function(x){
    marker_c0 <- FindMarkers(subset(obj,subtype%in%c(x)),ident.1 = 'Aged',ident.2 = 'Young')
    marker_c1 <- marker_c0%>%dplyr::filter(p_val<0.05&abs(avg_log2FC)>0.8) %>% 
      mutate(diff = case_when(avg_log2FC>0.8~'up',
                              avg_log2FC<(-0.8)~'down'),
             cluster=rep(x,length(.$avg_log2FC)),
             gene=row.names(.))
    
    return(marker_c1)
  })
  return(marker_list)
}
Myeloid_human <- subset(Myeloid_species_final,species%in%'human')
Myeloid_monkey <- subset(Myeloid_species_final,species%in%'monkey')
Myeloid_mouse2 <- subset(Myeloid_species_final,species%in%'mouse2')
Myeloid_mouse3 <- subset(Myeloid_species_final,species%in%'mouse3')
Myeloid_human_diff_list <- run_cluster_diff(Myeloid_human)
Myeloid_monkey_diff_list <- run_cluster_diff(Myeloid_monkey)
Myeloid_mouse2_diff_list <- run_cluster_diff(Myeloid_mouse2)
Myeloid_mouse3_diff_list <- run_cluster_diff(Myeloid_mouse3)
Myeloid_human_diff <- Reduce(rbind,Myeloid_human_diff_list) 
Myeloid_monkey_diff <- Reduce(rbind,Myeloid_monkey_diff_list) 
Myeloid_mouse2_diff <- Reduce(rbind,Myeloid_mouse2_diff_list) 
Myeloid_mouse3_diff <- Reduce(rbind,Myeloid_mouse3_diff_list) 

FB_human <- subset(FB_species_rpca_sub_label,species%in%'human')
FB_monkey <- subset(FB_species_rpca_sub_label,species%in%'monkey')
FB_mouse2 <- subset(FB_species_rpca_sub_label,species%in%'mouse2')
FB_mouse3 <- subset(FB_species_rpca_sub_label,species%in%'mouse3')
FB_human_diff_list <- run_cluster_diff(FB_human)
FB_monkey_diff_list <- run_cluster_diff(FB_monkey)
FB_mouse2_diff_list <- run_cluster_diff(FB_mouse2)
FB_mouse3_diff_list <- run_cluster_diff(FB_mouse3)
FB_human_diff <- Reduce(rbind,FB_human_diff_list) 
FB_monkey_diff <- Reduce(rbind,FB_monkey_diff_list) 
FB_mouse2_diff <- Reduce(rbind,FB_mouse2_diff_list) 

data_veen <- function(data1,data2,data3,
                      clusters=clusters){
  lapply(clusters, function(cl){
    veen_up <- list('human'=filter(data1,cluster%in%cl&diff%in%'up') %>% .$gene,
                    'monkey'=filter(data2,cluster%in%cl&diff%in%'up') %>% .$gene,
                    'mouse'=filter(data3,cluster%in%cl&diff%in%'up') %>% .$gene)
    veen_down <- list('human'=filter(data1,cluster%in%cl&diff%in%'down') %>% .$gene,
                      'monkey'=filter(data2,cluster%in%cl&diff%in%'down') %>% .$gene,
                      'mouse'=filter(data3,cluster%in%cl&diff%in%'down') %>% .$gene)
    return(list('veen_up'=veen_up,
                'veen_down'=veen_down))
  })->clusters_venn
  names(clusters_venn) <- clusters
  return(clusters_venn)
}
data_veen_Myeloid <- data_veen(Myeloid_human_diff,
                               Myeloid_monkey_diff,
                               Myeloid_mouse2_diff,
                               clusters=c('INF_MC','MHClow_MC','CCR2_MC','MHChi_MC','Mono','Timd4_resident_MC'))
data_veen_FB <- data_veen(FB_human_diff,
                          FB_monkey_diff,
                          FB_mouse2_diff,
                         clusters=c('FB1','FB2','FB3','FB4'))


INF_MC_veen_up <- get.venn.partitions(data_veen_Myeloid$INF_MC$veen_up)
INF_MC_veen_down <- get.venn.partitions(data_veen_Myeloid$INF_MC$veen_down)

MHClow_MC_veen_up <- get.venn.partitions(data_veen_Myeloid$MHClow_MC$veen_up)
MHClow_MC_veen_down <- get.venn.partitions(data_veen_Myeloid$MHClow_MC$veen_down)

CCR2_MC_veen_up <- get.venn.partitions(data_veen_Myeloid$CCR2_MC$veen_up)
CCR2_MC_veen_down <- get.venn.partitions(data_veen_Myeloid$CCR2_MC$veen_down)

MHChi_MC_veen_up <- get.venn.partitions(data_veen_Myeloid$MHChi_MC$veen_up)
MHChi_MC_veen_down <- get.venn.partitions(data_veen_Myeloid$MHChi_MC$veen_down)

Mono_veen_up <- get.venn.partitions(data_veen_Myeloid$Mono$veen_up)
Mono_veen_down <- get.venn.partitions(data_veen_Myeloid$Mono$veen_down)

Timd4_resident_MC_veen_up <- get.venn.partitions(data_veen_Myeloid$Timd4_resident_MC$veen_up)
Timd4_resident_MC_veen_down <- get.venn.partitions(data_veen_Myeloid$Timd4_resident_MC$veen_down)


FB1_veen_up <- get.venn.partitions(data_veen_FB$FB1$veen_up)
FB1_veen_down <- get.venn.partitions(data_veen_FB$FB1$veen_down)

FB2_veen_up <- get.venn.partitions(data_veen_FB$FB2$veen_up)
FB2_veen_down <- get.venn.partitions(data_veen_FB$FB2$veen_down)

FB3_veen_up <- get.venn.partitions(data_veen_FB$FB3$veen_up)
FB3_veen_down <- get.venn.partitions(data_veen_FB$FB3$veen_down)

FB4_veen_up <- get.venn.partitions(data_veen_FB$FB4$veen_up)
FB4_veen_down <- get.venn.partitions(data_veen_FB$FB4$veen_down)
plot_veen <- function(data){
  venn.plot_up <- venn.diagram(
    x = data,filename = NULL,
    category.names = c("human", "monkey", "mouse"),
    fill = c("green", "yellow", "purple"),
    lty = "blank",cex = 1.5,fontface = "bold",
    fontfamily = "sans",cat.cex = 1.5,
    cat.fontface = "bold",cat.fontfamily = "sans",
    cat.pos = 0,cat.dist = 0.05,alpha = 0.3)
  return(venn.plot_up)
}
p_INF_MC_veen_up <- plot_veen(data_veen_Myeloid$INF_MC$veen_up)
p_INF_MC_veen_down <- plot_veen(data_veen_Myeloid$INF_MC$veen_down)

p_MHClow_MC_veen_up <- plot_veen(data_veen_Myeloid$MHClow_MC$veen_up)
p_MHClow_MC_veen_down <- plot_veen(data_veen_Myeloid$MHClow_MC$veen_down)

p_CCR2_MC_veen_up <- plot_veen(data_veen_Myeloid$CCR2_MC$veen_up)
p_CCR2_MC_veen_down <- plot_veen(data_veen_Myeloid$CCR2_MC$veen_down)

p_MHChi_MC_veen_up <- plot_veen(data_veen_Myeloid$MHChi_MC$veen_up)
p_MHChi_MC_veen_down <- plot_veen(data_veen_Myeloid$MHChi_MC$veen_down)

p_Mono_veen_up <- plot_veen(data_veen_Myeloid$Mono$veen_up)
p_Mono_veen_down <- plot_veen(data_veen_Myeloid$Mono$veen_down)

p_Timd4_resident_MC_veen_up <- plot_veen(data_veen_Myeloid$Timd4_resident_MC$veen_up)
p_Timd4_resident_MC_veen_down <- plot_veen(data_veen_Myeloid$Timd4_resident_MC$veen_down)
pdf(hu("res", "p_INF_MC_veen_up.pdf"),width = 6,height = 6)
grid.draw(p_INF_MC_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_INF_MC_veen_down.pdf"),width = 6,height = 6)
grid.draw(p_INF_MC_veen_down) %>% print()
dev.off()
pdf(hu("res", "p_MHClow_MC_veen_up.pdf"),width = 6,height = 6)
grid.draw(p_MHClow_MC_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_MHClow_MC_veen_up.pdf"),width = 6,height = 6)
grid.draw(p_MHClow_MC_veen_down) %>% print()
dev.off()
pdf(hu("res", "p_CCR2_MC_veen_up.pdf"),width = 6,height = 6)
grid.draw(p_CCR2_MC_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_CCR2_MC_veen_down.pdf"),width = 6,height = 6)
grid.draw(p_CCR2_MC_veen_down) %>% print()
dev.off()
pdf(hu("res", "p_MHChi_MC_veen_up.pdf"),width = 6,height = 6)
grid.draw(p_MHChi_MC_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_MHChi_MC_veen_down.pdf"),width = 6,height = 6)
grid.draw(p_MHChi_MC_veen_down) %>% print()
dev.off()
pdf(hu("res", "p_Mono_veen_up.pdf"),width = 6,height = 6)
grid.draw(p_Mono_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_Mono_veen_down.pdf"),width = 6,height = 6)
grid.draw(p_Mono_veen_down) %>% print()
dev.off()
pdf(hu("res", "p_Timd4_resident_MC_veen_up.pdf"),width = 6,height = 6)
grid.draw(p_Timd4_resident_MC_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_Timd4_resident_MC_veen_down.pdf"),width = 6,height = 6)
grid.draw(p_Timd4_resident_MC_veen_down) %>% print()
dev.off()


p_FB1_veen_up <- plot_veen(data_veen_FB$FB1$veen_up)
p_FB1_veen_down <- plot_veen(data_veen_FB$FB1$veen_down)

p_FB2_veen_up <- plot_veen(data_veen_FB$FB2$veen_up)
p_FB2_veen_down <- plot_veen(data_veen_FB$FB2$veen_down)

p_FB3_veen_up <- plot_veen(data_veen_FB$FB3$veen_up)
p_FB3_veen_down <- plot_veen(data_veen_FB$FB3$veen_down)

p_FB4_veen_up <- plot_veen(data_veen_FB$FB4$veen_up)
p_FB4_veen_down <- plot_veen(data_veen_FB$FB4$veen_down)

pdf(hu("res", "p_FB4_veen_up1.pdf"),width = 6,height = 6)
grid.draw(p_FB4_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_FB4_veen_down1.pdf"),width = 6,height = 6)
grid.draw(p_FB4_veen_down) %>% print()
dev.off()
pdf(hu("res", "p_FB3_veen_up1.pdf"),width = 6,height = 6)
grid.draw(p_FB3_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_FB3_veen_down1.pdf"),width = 6,height = 6)
grid.draw(p_FB3_veen_down) %>% print()
dev.off()
pdf(hu("res", "p_FB2_veen_up1.pdf"),width = 6,height = 6)
grid.draw(p_FB2_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_FB2_veen_down1.pdf"),width = 6,height = 6)
grid.draw(p_FB2_veen_down) %>% print()
dev.off()
pdf(hu("res", "p_FB1_veen_up1.pdf"),width = 6,height = 6)
grid.draw(p_FB1_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_FB1_veen_down1.pdf"),width = 6,height = 6)
grid.draw(p_FB1_veen_down) %>% print()
dev.off()

Myeloid_species_final@meta.data$new_group <- paste(Myeloid_species_final$species,
                                                   Myeloid_species_final$aged,
                                                   Myeloid_species_final$subtype,sep = '_')
Myeloid_species_final@meta.data$new_group  <- as.factor(Myeloid_species_final@meta.data$new_group)
Idents(Myeloid_species_final) <- 'new_group'
Myeloid_species_final_exp <- AverageExpression(Myeloid_species_final,assays = 'RNA')

mat_Myeloid = Myeloid_species_final_exp$RNA %>% as.matrix() %>% t()
Myeloid_genes_conserved_up <- c(INF_MC_veen_up[1,5][[1]],
                        MHClow_MC_veen_up[1,5][[1]],
                        CCR2_MC_veen_up[1,5][[1]],
                        MHChi_MC_veen_up[1,5][[1]],
                        Mono_veen_up[1,5][[1]],
                        Timd4_resident_MC_veen_up[1,5][[1]])%>% unique()
Myeloid_genes_conserved_down <-c(INF_MC_veen_down[1,5][[1]],
                         MHClow_MC_veen_down[1,5][[1]],
                         CCR2_MC_veen_down[1,5][[1]],
                         MHChi_MC_veen_down[1,5][[1]],
                         Mono_veen_down[1,5][[1]],
                         Timd4_resident_MC_veen_down[1,5][[1]])%>% unique()
Myeloid_genes_species_huamn <- c(INF_MC_veen_up[7,5][[1]],
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
Myeloid_genes_species_monkey <- c(INF_MC_veen_up[6,5][[1]],
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
Myeloid_genes_species_mouse <- c(INF_MC_veen_up[4,5][[1]],
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

FB_species_rpca_sub_label@meta.data$new_group <- paste(FB_species_rpca_sub_label$species,
                                                   FB_species_rpca_sub_label$aged,
                                                   FB_species_rpca_sub_label$subtype,sep = '_')
FB_species_rpca_sub_label@meta.data$new_group  <- as.factor(FB_species_rpca_sub_label@meta.data$new_group)
Idents(FB_species_rpca_sub_label) <- 'new_group'
FB_species_rpca_sub_label <- AverageExpression(FB_species_rpca_sub_label,assays = 'RNA')

mat_FB = FB_species_rpca_sub_label$RNA %>% as.matrix() %>% t()
head(mat_FB)
FB_genes_conserved_up <- c(FB1_veen_up[1,5][[1]],
                           FB2_veen_up[1,5][[1]],
                           FB3_veen_up[1,5][[1]],
                           FB4_veen_up[1,5][[1]])%>% unique()
FB_genes_conserved_down <-c(FB1_veen_up[1,5][[1]],
                            FB2_veen_up[1,5][[1]],
                            FB3_veen_up[1,5][[1]],
                            FB4_veen_up[1,5][[1]])%>% unique()
FB_genes_species_huamn <- c(FB1_veen_up[7,5][[1]],
                            FB2_veen_up[7,5][[1]],
                            FB3_veen_up[7,5][[1]],
                            FB4_veen_up[7,5][[1]],
                            FB5_veen_up[7,5][[1]],
                            FB1_veen_down[7,5][[1]],
                            FB2_veen_down[7,5][[1]],
                            FB3_veen_down[7,5][[1]],
                            FB4_veen_down[7,5][[1]],
                            FB5_veen_down[7,5][[1]])%>% unique()
FB_genes_species_monkey <- c(FB1_veen_up[6,5][[1]],
                             FB2_veen_up[6,5][[1]],
                             FB3_veen_up[6,5][[1]],
                             FB4_veen_up[6,5][[1]],
                             FB5_veen_up[6,5][[1]],
                             FB1_veen_down[6,5][[1]],
                             FB2_veen_down[6,5][[1]],
                             FB3_veen_down[6,5][[1]],
                             FB4_veen_down[6,5][[1]],
                             FB5_veen_down[6,5][[1]])%>% unique()
FB_genes_species_mouse <- c(FB1_veen_up[4,5][[1]],
                            FB2_veen_up[4,5][[1]],
                            FB3_veen_up[4,5][[1]],
                            FB4_veen_up[4,5][[1]],
                            FB5_veen_up[4,5][[1]],
                            FB1_veen_down[4,5][[1]],
                            FB2_veen_down[4,5][[1]],
                            FB3_veen_down[4,5][[1]],
                            FB4_veen_down[4,5][[1]],
                            FB5_veen_down[4,5][[1]])%>% unique()
mat_corss <- mat_FB[,unique(c(FB_genes_conserved_up,
                              FB_genes_conserved_down,
                              FB_genes_species_huamn,
                              FB_genes_species_monkey,
                              FB_genes_species_mouse))]
dim(mat_corss)
group = rep(c("Conserved", "Species_human","Species_monkey","Species_mouse"),
            times = c(length(c(genes_conserved_up,genes_conserved_down)),
                      length(genes_species_huamn),
                      length(FB_genes_species_monkey),
                      length(genes_species_mouse)))
length(group)
col1 = colorRamp2(c(-2, 2), c( "white", "red"))
col2 = colorRamp2(c(-2, 2), c("white", "blue"))
col3 = colorRamp2(c(-2, 2), c("white", "green"))
col4 = colorRamp2(c(-2, 2), c("white", "pink"))
clusters <- rep(c('FB1','FB2','FB3','FB4'), 6)
aged <- rep(c(rep('Young',4),rep('Aged',4)),3)
clusters_col <- rep(c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"), 6)
aged_col <- rep(c(rep('#7570B3',4),rep('#FFED6F',4)),3)
names(clusters_col) <- clusters
names(aged_col) <- aged
left_annotation= rowAnnotation(clusters = clusters,
                               aged=aged, 
                               col = list(clusters = clusters_col,
                                          aged=aged_col))
mat_corss_scale <- scale(mat_corss)
dim(mat_corss_scale)
mat_corss_scale <- mat_corss_scale[1:24,]
ht1 = Heatmap(mat_corss_scale[, group == "Conserved"], col = col1, name = "Conserved",
              left_annotation = left_annotation,
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht1
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
ht1 + ht2 + ht3 + ht4
