# =============================================================================
# 01_FB_subcluster_rpca.R
# Title : Fibroblast cross-species sub-clustering
# Figure: Fig. 2A-B
# Module: 03_figure2_fibroblast_states
# Description:
#   Subset fibroblasts, down-sample human nuclei, RPCA-integrate, cluster and annotate FB_1-FB_5; marker and proportion plots.
# Inputs:
#   - hu("res", "fb_ann", "FB_cross.rds")
#   - hu("res", "FB_species_rpca_sub_label.rds")
#   - hu("data", "FB_species_rpca_sub.rds")
# Outputs:
#   - hu("res", "FB_species_rpca_sub_label.rds")
#   - hu("res", "FB_dim_split.pdf")
#   - hu("res", "FB_dim.pdf")
#   - hu("res", "FB_prop.pdf")
#   - hu("res", "fb_ann", "fea_FB.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

library(Seurat)
FB_species <- subset(combine_species_rpca,celltypes%in%'FB')
FB_species <- JoinLayers(FB_species)
FB_human_cells <- colnames(FB_species)[FB_species$species%in%'human']

FB_species_sub <- FB_species[,-sample(which(FB_species$species%in%'human'),30157)]
table(FB_species_sub$species)
DimPlot(FB_species_sub,split.by = 'species')
DimPlot(FB_species,split.by = 'species')

run_rpca <- function(sceobj,nfeatures,dims){
  sceobj <- JoinLayers(sceobj)
  sceobj_exp <- LayerData(sceobj, assay="RNA", layer='counts')
  sceobj_new <- CreateSeuratObject(sceobj_exp)
  identical(row.names(sceobj_new),row.names(sceobj))
  sceobj_new@meta.data <- sceobj@meta.data
  sceobj_new[["RNA"]] <- split(sceobj_new[["RNA"]], 
                               f = sceobj_new$orig.ident)
  sceobj_new <- NormalizeData(sceobj_new)
  sceobj_new <- FindVariableFeatures(sceobj_new,nfeatures = nfeatures)
  sceobj_new <- ScaleData(sceobj_new)
  sceobj_new <- RunPCA(sceobj_new)
  sceobj_new <- IntegrateLayers(
    object = sceobj_new, method = RPCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.rpca",
    verbose = T)
  sceobj_new <- FindNeighbors(sceobj_new,
                              reduction = "integrated.rpca", dims = dims)
  sceobj_new <- FindClusters(sceobj_new,resolution = c(0.2,0.4,0.6,0.8))
  sceobj_new <- RunUMAP(sceobj_new,
                        reduction = "integrated.rpca",
                        dims = dims, reduction.name = "umap.rpca")
  return(sceobj_new)
}
FB_species_rpca <- run_rpca(FB_species_sub,nfeatures = 1000,dims = 1:15)
FB_species_rpca_dim30 <- run_rpca(FB_species_sub,nfeatures = 1000,dims = 1:30)
FB_species_rpca_dim30_2000 <- run_rpca(FB_species_sub,nfeatures = 2000,dims = 1:30)

FB_species_rpca_dim30_2000 <- run_rpca(FB_cross_v2_new,nfeatures = 2000,dims = 1:30)

DimPlot(FB_species_rpca,split.by = 'species',group.by = 'RNA_snn_res.0.2',label = T)
DimPlot(FB_species_rpca_dim30,split.by = 'species',group.by = 'RNA_snn_res.0.2',label = T)

FB_species_rpca <- read_rds(hu("res", "fb_ann", "FB_cross.rds"))
FB_species_rpca <- FindClusters(FB_species_rpca,resolution = c(0.4,0.6,0.8,0.2))
DimPlot(FB_species_rpca,split.by = 'species',label = T)
FB_species_rpca <- JoinLayers(FB_species_rpca)

Idents(FB_species_rpca) <- 'RNA_snn_res.0.2'
FB_species_rpca <- JoinLayers(FB_species_rpca)
FB_marker <- FindAllMarkers(FB_species_rpca,min.pct = 0.2)
FB_marker_top10 <- FB_marker %>% dplyr::filter(pct.1>0.3)%>% dplyr::group_by(cluster) %>% dplyr::top_n(10,avg_log2FC)
FB_marker_top50 <- FB_marker %>% dplyr::filter(pct.1>0.3)%>% dplyr::group_by(cluster) %>% dplyr::top_n(50,avg_log2FC)
FB_marker_cut1 <- FB_marker %>% dplyr::filter(avg_log2FC>1&pct.1>0.3)
FB_marker_cut1_list <- split(FB_marker_cut1$gene,FB_marker_cut1$cluster)
FB_species_rpca_sub <- subset(FB_species_rpca,RNA_snn_res.0.2%in%c(0,1,2,3,7))
label <- c('0'="FB1",
           '1'="FB3",
           '2'="FB2",
           '3'="FB4",
           '7'='FB5')   
Idents(FB_species_rpca_sub) <- 'RNA_snn_res.0.2'
FB_species_rpca_sub_label <- Seurat::RenameIdents(FB_species_rpca_sub,label)
FB_species_rpca_sub_label$subtype <- Idents(FB_species_rpca_sub_label)
DimPlot(FB_species_rpca_sub_label,split.by = 'species',label = T)
DimPlot(FB_species_rpca_sub_label,label = T)
FB_species_rpca_sub_label <- JoinLayers(FB_species_rpca_sub_label)
FB_sub_marker <- FindAllMarkers(FB_species_rpca_sub_label,min.pct = 0.3)
FB_sub_marker_cut1 <- FB_sub_marker %>% dplyr::filter(avg_log2FC>1&pct.1>0.25)
FB_sub_marker_cut1_list <- split(FB_sub_marker_cut1$gene,FB_sub_marker_cut1$cluster)
FB_sub_marker_cut0.8 <- FB_sub_marker %>% dplyr::filter(avg_log2FC>0.8&pct.1>0.2)
FB_sub_marker_cut0.8_list <- split(FB_sub_marker_cut0.8$gene,
                                   FB_sub_marker_cut0.8$cluster)
FB_sub_marker_top10 <- FB_sub_marker %>% dplyr::filter(pct.1>0.1)%>% 
                       dplyr::group_by(cluster) %>% 
                       dplyr::top_n(10,avg_log2FC)
FB_sub_marker_top30 <- FB_sub_marker %>% dplyr::filter(pct.1>0.3)%>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(30,avg_log2FC)
write_rds(FB_species_rpca_sub_label,file = hu("res", "FB_species_rpca_sub_label.rds"))
FB_species_rpca_sub_label <- read_rds(hu("res", "FB_species_rpca_sub_label.rds"))
pdf(hu("res", "FB_dim_split.pdf"),width = 15,height = 5)
DimPlot(FB_species_rpca_sub_label,label = T,split.by = 'species')
dev.off()
pdf(hu("res", "FB_dim.pdf"),width = 15,height = 5)
DimPlot(FB_species_rpca_sub_label,label = T,group.by = 'species')
dev.off()
# prop --------------------------------------------------------------------
# prop --------------------------------------------------------------------
prop_aged <- data.frame(table(FB_species_rpca_sub_label$subtype,
                              FB_species_rpca_sub_label$aged)) %>%  
  dplyr::group_by(Var2) %>% dplyr::mutate(percent = Freq/sum(Freq)) %>% 
  dplyr::mutate(label = paste0(round(percent*100,digits = 1), "%"))
prop_species <- data.frame(table(FB_species_rpca_sub_label$subtype,
                                 FB_species_rpca_sub_label$species)) %>%  
  dplyr::group_by(Var2) %>% dplyr::mutate(percent = Freq/sum(Freq)) %>% 
  dplyr::mutate(label = paste0(round(percent*100,digits = 1), "%"))

ggplot(mapping = aes(x = prop_aged$Var1,y = prop_aged$percent,
                     fill = prop_aged$Var2)) + 
  geom_bar(stat = 'identity', position = "stack") + 
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  geom_text(aes(label = prop_aged$label),size=5,position = position_stack(.9)) +
  labs(x="group",y="Percentage",fill="cluster") +
  theme(panel.background = element_rect(fill = NA, color = NA)) +
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))

ggplot(mapping = aes(x = prop_species$Var1,y = prop_species$percent,
                     fill = prop_species$Var2)) + 
  geom_bar(stat = 'identity', position = "stack") + 
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  geom_text(aes(label = prop_species$label),size=4,position = position_stack(.9)) +
  labs(x="group",y="Percentage",fill="cluster") +
  theme(panel.background = element_rect(fill = NA, color = NA)) +
  theme_classic()+
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12))
cluster_species_counts <- FB_species_rpca_sub_label@meta.data %>%
  dplyr::group_by(species,subtype) %>%
  summarise(count = n()) %>%
  dplyr::ungroup()

# 计算每个物种中 cluster 的百分比
cluster_species_percentage <- cluster_species_counts %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(percentage = count / sum(count) * 100) %>%
  dplyr::ungroup()
# 绘制图形
pdf(hu("res", "FB_prop.pdf"),width = 8,height = 8)
ggplot(cluster_species_percentage, aes(x = species, y = percentage, fill = subtype)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Fraction(%)", fill = "FB cluster") +
  geom_text(aes(label = sprintf("%.1f%%", percentage), y = percentage), 
            position = position_stack(vjust = 0.5), size = 5, color = "black")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,size = 15, hjust = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"))
dev.off()

DotPlot(FB_species_rpca_sub_label,
        features = c())+
  scale_color_viridis_c()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
FeaturePlot(FB_species_rpca_sub_label,c('CD59','CD44','LTBP2','CCN2','NFATC1'),order = T)
FeaturePlot(FB_species_rpca_sub_label,c('ERBB4','DKK3','KCNMA1','ECRG4'))
FeaturePlot(FB_species_rpca_sub_label,c('INMT','HSD11B1','LPL','APOE'))
FeaturePlot(FB_species_rpca_sub_label,c('PI16','CD248','ACKR3','SEMA3C'))
FeaturePlot(FB_species_rpca_sub_label,c('IFIT1B','IFIT3','ISG15','USP18'))
FeaturePlot(FB_species_rpca_sub_label,c('IFIT1B','IFIT3','ISG15','USP18'))
FeaturePlot(FB_species_rpca_sub_label,c('LUM','DCN','GSN','COL1A1'))
# enrich plot --------------------------------------------------------------
library(enrichplot)
enrich_res_FB_1 <- clusterProfiler::compareCluster(FB_sub_marker_cut1_list,
                                                 fun='enrichGO',
                                                 OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                                                 ont = "BP",keyType = "SYMBOL",
                                                 minGSSize=1,pvalueCutoff = 0.05,
                                                 pAdjustMethod = "BH")
enrich_res_FB_0.8 <- clusterProfiler::compareCluster(FB_sub_marker_cut0.8_list,
                                                 fun='enrichGO',
                                                 OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                                                 ont = "BP",keyType = "SYMBOL",
                                                 minGSSize=1,pvalueCutoff = 0.05,
                                                 pAdjustMethod = "BH")
dotplot(enrich_res_FB_1)
dotplot(enrich_res_FB_0.8)
enrich_res_FB <- clusterProfiler::compareCluster(FB_enrich_res,
                                                 fun='enrichGO',OrgDb=org.Mm.eg.db,
                                                 ont = "BP",keyType = "SYMBOL",
                                                 minGSSize=1,pvalueCutoff = 0.05,
                                                 pAdjustMethod = "BH")
barplot(enrich_res_FB)

FeaturePlot(FB_species_rpca_sub_label,c('PTPRG','MALAT1','C7','NEAT1','ABCA8'))
VlnPlot(FB_species_rpca_sub_label,c('PTPRG','MALAT1','C7','NEAT1','ABCA8'))
pdf(hu("res", "fb_ann", "fea_FB.pdf"),width = 12,height = 12)
FeaturePlot(FB_species_rpca_sub_label,c('DCN','GSN','LAMA2','BICC1','COL6A3',
                              'LAMC1', 'MGP','POSTN','COL1A1','TNC'))
FeaturePlot(FB_species_rpca_sub_label,c('MEOX1','LTBP2','FN1','COL8A1','RUNX1',
                              'ASPN'))
FeaturePlot(FB_species_rpca_sub_label,c('ACKR3','SEMA3C','CD248','GFPT2','CD55',
                              'PCOLCE2'))
FeaturePlot(FB_species_rpca_sub_label,c('HSD11B1','EGR1','G0S2','JUND','FOS',
                              'JUNB'))
FeaturePlot(FB_species_rpca_sub_label,c('ERBB4','DKK3','KCNMA1','ECRG4','COMP',
                                        'POSTN','THBS1','THBS2'))
dev.off()
FeaturePlot(FB_species_rpca_sub_label,c('MEOX1','FN1','RUNX1','POSTN',
                                        'DCN','GSN'))
FeaturePlot(FB_species_rpca_sub_label,c('PTPRG','ABCA8','DCN','MGP'))
FeaturePlot(FB_species_rpca_sub_label,c('HSPA1B','FOSB','DNAJB1','CEBPD'))
FeaturePlot(FB_species_rpca_sub_label,c('PI16','CD248','ACKR3','CD55'))

DimPlot(FB_species_rpca_sub,reduction = 'pca',dims = 2:3,label = T)
DimPlot(FB_species_rpca_sub,reduction = 'pca',dims = 3:4,label = T)
write_rds(FB_species_rpca_sub,
          file = hu("data", "FB_species_rpca_sub.rds"))
FB_species_rpca_sub <- read_rds(hu("data", "FB_species_rpca_sub.rds"))
# plot marker -------------------------------------------------------------
DimPlot(FB_species_rpca,label = T)
FB_species_rpca_sub <- subset(FB_species_rpca,seurat_clusters%in%c(0:5))
DimPlot(FB_species_rpca_sub,label = T)
FeaturePlot(FB_species_rpca,features = c('IFIT3','ISG15','RTP4','IFIT1B','XAF1',
                                         'RNF213'))
FeaturePlot(FB_species_rpca_sub,features = c('DCN','GSN','LAMA2','BICC1','COL6A3',
                                             'LAMC1', 'MGP','POSTN','COL1A1','TNC'))
FeaturePlot(FB_species_rpca_sub,features = c('KCNMA1','EDIL3','CLU','KCNQ5','PRRX1','SOX6'))
FeaturePlot(FB_species_rpca_sub,features = c('IL6','ICAM1','TNFAIP6','NFKBIA','FOSB','KDM6B'))#immfibro
FeaturePlot(FB_species_rpca_sub,features = c('MEOX1','LTBP2','PTPRE','FN1','CCN2'),order = T)#activate
FeaturePlot(FB_species_rpca_sub,features = c('PI16','PDPN','NPR3','CX3CL1','PROX1'),order = T)
FeaturePlot(FB_species_rpca_sub,features = c('POSTN','MEOX1','SERPINE1','FN1'))
FeaturePlot(FB_species_rpca_sub,features = c('MFAP5','ACKR3','PI16','PLXDC2','KLF4'))#stromal fibroblast

DimPlot(FB_species_rpca_sub,group.by = 'rpca_clusters',label = T)
DimPlot(FB_species_rpca_sub,group.by = 'rpca_clusters',split.by = 'species',label = T)



