# =============================================================================
# 02_fig1_umap_split_sankey_enrichment.R
# Title : Split UMAPs, composition Sankey, marker enrichment
# Figure: Fig. 1B; Supplementary
# Module: 02_figure1_cross_species_atlas
# Description:
#   Species/age-split UMAPs, cell-type composition Sankey and up-regulated-gene enrichment.
# Inputs:
#   - hu("res", "combine_species_new_label.rds")
# Outputs:
#   - hu("fig", "V1", "fig1_umap_all.pdf")
#   - hu("fig", "V1", "fig1_umap_species_split.pdf")
#   - hu("fig", "V1", "fig1_umap_aged_split.pdf")
#   - hu("fig", "V1", "fig1_celltype_sankey.pdf")
#   - hu("fig", "V1", "fig1_upgene_enrich.pdf")
#   - hu("fig", "V1", "fig1_downgene_enrich.pdf")
#   - hu("fig", "V1", "fig1_cross_alltypes_veen_up.pdf")
#   - hu("fig", "V1", "fig1_cross_alltypes_veen_down.pdf")
#   - hu("res", "fig1_cross_alltypes.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(Seurat)
packageVersion('Seurat')
# Example data (you need to replace this with your actual data)
celltypes_col <- c("#FF1493", "#FFB6C1", "#483D8B", "#00BFFF", "#228B22", "#7FFF00", "#FFFF00", "#FF4500","#8DD3C7")
species_col <- c("#FC8D62","#8DD3C7","#BC80BD")
aged_coll <- c("#483D8B","#00BFFF")
celltypes <- c('FB','EC','Myeloid','T','B','Pericytes','SMC','Neural','Adipocyte')
species <- c('human','monkey','mouse')
aged <- c('Aged','Young')
combine_species_new_label <- read_rds(hu("res", "combine_species_new_label.rds"))
combine_species_new_label$dataset <- combine_species_new_label$species
combine_species_new_label$species[combine_species_new_label$species%in%c('mouse1','mouse3')] <- 'mouse'
DimPlot(combine_species_new_label,split.by = 'species') 
pdf(hu("fig", "V1", "fig1_umap_all.pdf"),width = 8,height = 6)
DimPlot(combine_species_new_label,group.by = 'annotation',
        label = T,cols = celltypes_col) 
DimPlot(combine_species_new_label,group.by = 'species',
        label = T,cols = species_col) 
DimPlot(combine_species_new_label,group.by = 'aged',
        label = T,cols = aged_coll) 
dev.off()
pdf(hu("fig", "V1", "fig1_umap_species_split.pdf"),width = 12,height = 4)
Seurat::DimPlot(combine_species_new_label,group.by = 'annotation',split.by = "species",
        label = T,cols = celltypes_col) 
dev.off()
pdf(hu("fig", "V1", "fig1_umap_aged_split.pdf"),width = 8,height = 4)
DimPlot(combine_species_new_label,group.by = 'annotation',split.by = "aged",
        label = T,cols = celltypes_col) 
dev.off()
# plot sankey -------------------------------------------------------------
library(RColorBrewer)
library(ggalluvial)
library(tidyverse)
new_meta <- combine_species_new_label@meta.data
cluster_counts_celltype <- new_meta %>%
  dplyr::group_by(species, aged, annotation) %>%
  summarize(count = n(), .groups = 'drop')
# 确保数据格式适合绘制桑基图
cluster_counts_celltype$celltypes <- factor(cluster_counts_celltype$annotation,
                                            levels = celltypes)
# 创建一个颜色映射列表，用于标签颜色
label_sankey <- c(species,aged,celltypes)
label_col <- c(species_col,aged_coll,celltypes_col)
names(label_col) <- label_sankey
pdf(hu("fig", "V1", "fig1_celltype_sankey.pdf"),width = 10,height = 8)
ggplot(cluster_counts_celltype, aes(axis1 = species, axis2 = annotation,axis3 = aged, y = count)) +
  geom_alluvium(aes(fill = aged), width = 0.1) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.1, color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),color = 'black',
            size = 4, vjust = -0.5) +
  scale_x_discrete(limits = c("Species", "celltypes", "aged"), expand = c(0.15, 0.05)) +
  scale_fill_manual(values = label_col) +
  scale_color_manual(values = label_col) +
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Sankey Diagram of Clusters Across Species and Groups",
       fill = "Cluster",
       color = "Group")
dev.off()
# all cells in aging change -----------------------------------------------
#write_rds(combine_species_new_label,hu("res", "combine_species_new_label.rds"))
combine_species_new_label_human <- subset(combine_species_new_label,dataset%in%'human')
combine_species_new_label_monkey <- subset(combine_species_new_label,dataset%in%'monkey')
combine_species_new_label_mouse <- subset(combine_species_new_label,dataset%in%'mouse1')
# find all cells diff gene ------------------------------------------------
Idents(combine_species_new_label_human) <- 'aged'
human_aging_diff <- FindMarkers(combine_species_new_label_human,ident.1 = 'Aged',ident.2 = 'Young')
Idents(combine_species_new_label_monkey) <- 'aged'
monkey_aging_diff <- FindMarkers(combine_species_new_label_monkey,ident.1 = 'Aged',ident.2 = 'Young')
Idents(combine_species_new_label_mouse) <- 'aged'
mouse_aging_diff <- FindMarkers(combine_species_new_label_mouse,ident.1 = 'Aged',ident.2 = 'Young')

human_aging_diff_cut0.8 <- human_aging_diff%>%dplyr::filter(p_val<0.05&abs(avg_log2FC)>0.8) %>% 
                         mutate(diff = case_when(avg_log2FC>0.8~'up',
                                                 avg_log2FC<(-0.8)~'down'),
                                gene=row.names(.))
monkey_aging_diff_cut0.8 <- monkey_aging_diff%>%dplyr::filter(p_val<0.05&abs(avg_log2FC)>0.8) %>% 
  mutate(diff = case_when(avg_log2FC>0.8~'up',avg_log2FC<(-0.8)~'down'),
         gene=row.names(.))
mouse_aging_diff_cut0.8 <- mouse_aging_diff%>%dplyr::filter(p_val<0.05&abs(avg_log2FC)>0.8) %>% 
  mutate(diff = case_when(avg_log2FC>0.8~'up',avg_log2FC<(-0.8)~'down'),
         gene=row.names(.))
DotPlot(combine_species_new_label_human,features = c('MRVI1','ITM2B','IL1R2','MGP'))
DotPlot(combine_species_new_label_human,features = c('MRVI1','ITM2B','IL1R2','MGP'))



genes_human_up <- filter(human_aging_diff_cut0.8,diff=='up') %>% .$gene %>% unique()
genes_monkey_up <- filter(monkey_aging_diff_cut0.8,diff=='up') %>% .$gene %>% unique()
genes_mouse_up <- filter(mouse_aging_diff_cut0.8,diff=='up') %>% .$gene %>% unique()
genes_human_down <- filter(human_aging_diff_cut0.8,diff=='down') %>% .$gene %>% unique()
genes_monkey_down <- filter(monkey_aging_diff_cut0.8,diff=='down') %>% .$gene %>% unique()
genes_mouse_down <- filter(mouse_aging_diff_cut0.8,diff=='down') %>% .$gene %>% unique()

genes_human_up_enrich <- enrichGO(gene = genes_human_up,
                            OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                            ont = "BP",keyType = "SYMBOL",
                            minGSSize=1,pvalueCutoff = 0.05,
                            pAdjustMethod = "BH")
genes_human_down_enrich <- enrichGO(gene = genes_human_down,
                                  OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                                  ont = "BP",keyType = "SYMBOL",
                                  minGSSize=1,pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH")
aging_term <- c('cellular response to interleukin-1',  
                'response to tumor necrosis factor',
                'cytokine-mediated signaling pathway',
                'cell chemotaxis',
                'heart contraction',
                'response to oxidative stress',
                'extracellular matrix organization',
                'response to toxic substance')
dotplot(filter(genes_human_up_enrich,Description%in%aging_term))
ego_human_up <- pairwise_termsim(genes_human_up_enrich)
ego_human_down <- pairwise_termsim(genes_human_down_enrich)
# monkey go enrich --------------------------------------------------------
genes_monkey_up_enrich <- enrichGO(gene = genes_monkey_up,
                                  OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                                  ont = "BP",keyType = "SYMBOL",
                                  minGSSize=1,pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH")
genes_monkey_down_enrich <- enrichGO(gene = genes_monkey_down,
                                   OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                                   ont = "BP",keyType = "SYMBOL",
                                   minGSSize=1,pvalueCutoff = 0.05,
                                   pAdjustMethod = "BH")
genes_monkey_up_enrich@result$Description[1:50]
ego_monkey_up <- pairwise_termsim(genes_monkey_up_enrich)
ego_monkey_down <- pairwise_termsim(genes_monkey_down_enrich)
# mosue go enrich ---------------------------------------------------------
genes_mouse_up_enrich <- enrichGO(gene = genes_mouse_up,
                                  OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                                  ont = "BP",keyType = "SYMBOL",
                                  minGSSize=1,pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH")
genes_mouse_down_enrich <- enrichGO(gene = genes_mouse_down,
                                  OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                                  ont = "BP",keyType = "SYMBOL",
                                  minGSSize=1,pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH")
genes_mouse_up_enrich@result$Description[1:50]
ego_mouse_up <- pairwise_termsim(genes_mouse_up_enrich)
ego_mouse_down <- pairwise_termsim(genes_mouse_down_enrich)
pdf(hu("fig", "V1", "fig1_upgene_enrich.pdf"),width = 8,height = 8)
emapplot(ego_human_up)
emapplot(ego_monkey_up)
emapplot(ego_mouse_up)
dev.off()
pdf(hu("fig", "V1", "fig1_downgene_enrich.pdf"),width = 8,height = 8)
emapplot(ego_human_down)
emapplot(ego_monkey_down)
emapplot(ego_mouse_down)
dev.off()
# cross species gene ------------------------------------------------------
gene_list_up <- list(human = genes_human_up, 
                     monkey = genes_monkey_up, 
                     mouse = c(genes_mouse_up))
gene_list_down <- list(human = genes_human_down, 
                       monkey = genes_monkey_down, 
                       mouse = c(genes_mouse_down))
gene_data_up <- fromList(gene_list_up)
gene_data_down <- fromList(gene_list_down)
pdf(hu("fig", "V1", "fig1_cross_alltypes_veen_up.pdf"),width = 6,height = 6)
venn.diagram(x = gene_list_up,category.names = c("human", "monkey", "mouse"),
             filename = NULL, output = FALSE,fill = c("green", "yellow", "purple"),alpha = 0.3) %>% grid.draw()
dev.off()
pdf(hu("fig", "V1", "fig1_cross_alltypes_veen_down.pdf"),width = 6,height = 6)
venn.diagram(x = gene_list_down,category.names = c("human", "monkey", "mouse"),
             filename = NULL, output = FALSE,fill = c("#FC8D62","#8DD3C7","#BC80BD"),alpha = 0.7) %>% grid.draw()
dev.off()
# 绘制 UpSet 图
pdf(hu("res", "fig1_cross_alltypes.pdf"),width = 12,height = 8)
upset(gene_data_up, 
      sets = c("human", "monkey", "mouse"), 
      keep.order = TRUE, 
      order.by = "freq",
      main.bar.color = "dodgerblue",
      sets.bar.color = "deepskyblue3",
      text.scale = c(2, 2, 2, 2, 2, 2))
upset(gene_data_down, 
      sets = c("human", "monkey", "mouse"), 
      keep.order = TRUE, 
      order.by = "freq",
      main.bar.color = "dodgerblue",
      sets.bar.color = "deepskyblue3",
      text.scale = c(2, 2, 2, 2, 2, 2))
dev.off()

up_veen <- get.venn.partitions(gene_list_up)
down_veen <- get.venn.partitions(gene_list_down)
writexl::write_xlsx(list(up_veen=up_veen,
                         down_veen=down_veen),
                    path = hu("fig", "V1", "fig1_cross_gene.xlsx"))
up_cross_gene <- c(up_veen[1,5][[1]])
down_cross_gene <- c(down_veen[1,5][[1]])


FeaturePlot(combine_species_new_label,features = up_cross_gene,split.by = 'aged')
FeaturePlot(combine_species_new_label,features = down_cross_gene,split.by = 'aged')
DotPlot(combine_species_new_label,features = up_cross_gene,
        group.by = 'new_group_aged')
DotPlot(combine_species_new_label,features = down_cross_gene,
        group.by = 'annotation',split.by = 'aged')


up_cross_enrich <- enrichGO(gene = up_cross_gene,
                            OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                            ont = "BP",keyType = "SYMBOL",
                            minGSSize=1,pvalueCutoff = 0.05,
                            pAdjustMethod = "BH")
down_cross_enrich <- enrichGO(gene = down_cross_gene,
                              OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                              ont = "BP",keyType = "SYMBOL",
                              minGSSize=1,pvalueCutoff = 0.05,
                              pAdjustMethod = "BH")
barplot(up_cross_enrich)
dotplot(down_cross_enrich)




celltype_human_markers <- FindAllMarkers(combine_species_new_label_human,min.pct = 0.2)
Idents(combine_species_new_label_monkey) <- 'annotation'
celltype_monkey_markers <- FindAllMarkers(combine_species_new_label_monkey,min.pct = 0.2)
Idents(combine_species_new_label_mouse) <- 'annotation'
celltype_mouse1_markers <- FindAllMarkers(combine_species_new_label_mouse,min.pct = 0.2)

celltypes <- c('FB','EC','Myeloid',
               'T','B','Pericytes',
               'SMC','Neural')

celltype_marker <- function(sce,celltypes){
  lapply(celltypes, function(cl){
    sub_sce <- subset(sce,annotation%in%cl)
    Idents(sub_sce) <- 'aged'
    age_diff <- FindMarkers(sub_sce,ident.1 = 'Aged',ident.2 = 'Young')
    marker_c1 <- age_diff%>%dplyr::filter(p_val<0.05&abs(avg_log2FC)>1) %>% 
      mutate(diff = case_when(avg_log2FC>1~'up',
                              avg_log2FC<(-1)~'down'),
             cluster=rep(cl,length(.$avg_log2FC)),
             gene=row.names(.))
    return(marker_c1)
  })->age_diff_list
  return(age_diff_list)
}
human_marker <- celltype_marker(combine_species_new_label_human,celltypes)
monkey_marker <- celltype_marker(combine_species_new_label_monkey,celltypes)
mouse_marker <- celltype_marker(combine_species_new_label_mouse,celltypes)
  
human_marker_all <- Reduce(rbind,human_marker)
monkey_marker_all <- Reduce(rbind,monkey_marker)
mouse_marker_all <- Reduce(rbind,mouse_marker)


genes_human_up <- filter(human_marker_all,diff=='up') %>% .$gene %>% unique()
genes_monkey_up <- filter(monkey_marker_all,diff=='up') %>% .$gene %>% unique()
genes_mouse_up <- filter(mouse_marker_all,diff=='up') %>% .$gene %>% unique()
genes_human_down <- filter(human_marker_all,diff=='down') %>% .$gene %>% unique()
genes_monkey_down <- filter(monkey_marker_all,diff=='down') %>% .$gene %>% unique()
genes_mouse_down <- filter(mouse_marker_all,diff=='down') %>% .$gene %>% unique()
# 创建基因集列表
gene_list_up <- list(human = genes_human_up, 
                     monkey = genes_monkey_up, 
                     mouse = c(genes_mouse_up))
gene_list_down <- list(human = genes_human_down, 
                       monkey = genes_monkey_down, 
                       mouse = c(genes_mouse_down))


up_veen <- get.venn.partitions(gene_list_up)
down_veen <- get.venn.partitions(gene_list_down)
up_cross_gene <- c(up_veen[1,5][[1]])
down_cross_gene <- c(down_veen[1,5][[1]])
up_cross_enrich <- enrichGO(gene = up_cross_gene,
                 OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                 ont = "BP",keyType = "SYMBOL",
                 minGSSize=1,pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
down_cross_enrich <- enrichGO(gene = down_cross_gene,
                            OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                            ont = "BP",keyType = "SYMBOL",
                            minGSSize=1,pvalueCutoff = 0.05,
                            pAdjustMethod = "BH")
barplot(up_cross_enrich)
dotplot(down_cross_enrich)

library(VennDiagram)
library(UpSetR)
# 将基因集列表转换为数据框格式
gene_data_up <- fromList(gene_list_up)
gene_data_down <- fromList(gene_list_down)
# 绘制 UpSet 图
pdf(hu("res", "fig1_cross_alltypes.pdf"),width = 12,height = 8)
upset(gene_data_up, 
      sets = c("human", "monkey", "mouse"), 
      keep.order = TRUE, 
      order.by = "freq",
      main.bar.color = "dodgerblue",
      sets.bar.color = "deepskyblue3",
      text.scale = c(2, 2, 2, 2, 2, 2))
upset(gene_data_down, 
      sets = c("human", "monkey", "mouse"), 
      keep.order = TRUE, 
      order.by = "freq",
      main.bar.color = "dodgerblue",
      sets.bar.color = "deepskyblue3",
      text.scale = c(2, 2, 2, 2, 2, 2))
dev.off()
