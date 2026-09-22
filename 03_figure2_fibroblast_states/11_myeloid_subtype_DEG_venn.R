# =============================================================================
# 11_myeloid_subtype_DEG_venn.R
# Title : Myeloid subtype DEG Venn intersections
# Figure: Supplementary
# Module: 03_figure2_fibroblast_states
# Description:
#   Venn diagrams of age-DEGs for IFN, MHC-low and CCR2+ myeloid subtypes.
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
#   - hu("res", "p_FB5_veen_up.pdf")
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
    marker_c1 <- marker_c0%>%dplyr::filter(p_val<0.05&abs(avg_log2FC)>1.2) %>% 
      mutate(diff = case_when(avg_log2FC>1~'up',
                              avg_log2FC<(-1)~'down'),
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
                         clusters=c('FB1','FB2','FB3','FB4','FB5'))


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

FB5_veen_up <- get.venn.partitions(data_veen_FB$FB5$veen_up)
FB5_veen_down <- get.venn.partitions(data_veen_FB$FB5$veen_down)
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

p_FB5_veen_up <- plot_veen(data_veen_FB$FB5$veen_up)
p_FB5_veen_down <- plot_veen(data_veen_FB$FB5$veen_down)
pdf(hu("res", "p_FB5_veen_up.pdf"),width = 6,height = 6)
grid.draw(p_FB5_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_FB5_veen_down.pdf"),width = 6,height = 6)
grid.draw(p_FB5_veen_down) %>% print()
dev.off()
pdf(hu("res", "p_FB4_veen_up.pdf"),width = 6,height = 6)
grid.draw(p_FB4_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_FB4_veen_down.pdf"),width = 6,height = 6)
grid.draw(p_FB4_veen_down) %>% print()
dev.off()
pdf(hu("res", "p_FB3_veen_up.pdf"),width = 6,height = 6)
grid.draw(p_FB3_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_FB3_veen_down.pdf"),width = 6,height = 6)
grid.draw(p_FB3_veen_down) %>% print()
dev.off()
pdf(hu("res", "p_FB2_veen_up.pdf"),width = 6,height = 6)
grid.draw(p_FB2_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_FB2_veen_down.pdf"),width = 6,height = 6)
grid.draw(p_FB2_veen_down) %>% print()
dev.off()
pdf(hu("res", "p_FB1_veen_up.pdf"),width = 6,height = 6)
grid.draw(p_FB1_veen_up) %>% print()
dev.off()
pdf(hu("res", "p_FB1_veen_down.pdf"),width = 6,height = 6)
grid.draw(p_FB1_veen_down) %>% print()
dev.off()
