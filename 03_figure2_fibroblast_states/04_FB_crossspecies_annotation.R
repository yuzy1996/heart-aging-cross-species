# =============================================================================
# 04_FB_crossspecies_annotation.R
# Title : Fibroblast cross-species annotation
# Figure: Fig. 2A
# Module: 03_figure2_fibroblast_states
# Description:
#   Translate macaque IDs via BioMart and annotate integrated fibroblast subtypes from markers.
# Inputs:
#   - cr("mart_export_monkey.txt")
#   - cr("human", "crossspecies", "sceobj", "FB_cross_v2_new.rds")
# Outputs:
#   - hu("res", "fb_ann", "fea_FB.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

# load data ---------------------------------------------------------------
trans_gene <- read_delim(cr("mart_export_monkey.txt")) 
trans_gene$monkeygene <- trans_gene$`Crab-eating macaque gene stable ID`
trans_gene$mousegene <- trans_gene$`Mouse gene stable ID`
trans_gene$humangene <- trans_gene$`Gene name`

FB_cross_v2_new <- read_rds(file=cr("human", "crossspecies", "sceobj", "FB_cross_v2_new.rds"))
FeaturePlot(FB_cross_v2_new,features = c(''))

combine_species_batch_noCM_V2_sub1 <- combine_species_batch_noCM_V2_sub[,sample(1:length(colnames(combine_species_batch_noCM_V2_sub)),10000)]
combine_species_batch_noCM_V2_sub1 <- JoinLayers(combine_species_batch_noCM_V2_sub1)
marker_celltype <- FindAllMarkers(combine_species_batch_noCM_V2_sub1,min.pct = 0.5)
marker_celltype_top30 <- marker_celltype %>% 
  dplyr::group_by(cluster) %>% dplyr::top_n(30,wt = avg_log2FC)

FeaturePlot(combine_species_batch,features = )

FB_cross_v2_new <- FindNeighbors(FB_cross_v2_new,reduction = 'umap.rpca',dims = 1:20)
FB_cross_v2_new <- FindClusters(FB_cross_v2_new,resolution = 0.2)

FB_cross_v2_new <- FindNeighbors(FB_cross_v2_new,reduction = "integrated.rpca", dims = 1:20)
FB_cross_v2_new <- FindClusters(FB_cross_v2_new,resolution = 0.2)
FB_cross_v2_new <- RunUMAP(FB_cross_v2_new,
                        reduction = "integrated.rpca",
                        dims = 1:20)
DimPlot(FB_cross_v2_new,group.by = 'RNA_snn_res.0.2',label = T)

marker_0_1 <- FindMarkers(FB_cross_v2_new,ident.1 ='0' ,ident.2 = '1',min.pct = 0.5)


pdf(hu("res", "fb_ann", "fea_FB.pdf"),width = 12,height = 12)
FeaturePlot(FB_cross_v2_new,c('DCN','GSN','LAMA2','BICC1','COL6A3',
                              'LAMC1', 'MGP','POSTN','COL1A1','TNC'))
FeaturePlot(FB_cross_v2_new,c('MEOX1','LTBP2','FN1','COL8A1','RUNX1',
                              'ASPN'),order = T)
FeaturePlot(FB_cross_v2_new,c('IFIT3','ISG15','RTP4','IFIT1B','XAF1',
                              'RNF213'),split.by = 'species')
FeaturePlot(FB_cross_v2_new,c('MYH6','TNNC1','MYL3','TNNT2','ATP2A2',
                              'ACTC1'))
dev.off()

FeaturePlot(FB_cross_v2_new,features = c('IFIT3','ISG15','XAF1',
                              'RNF213'),split.by = 'species')

FB_species <- subset(combine_species_rpca,celltypes%in%'EC')
FB_species_rpca <- run_rpca(FB_species,nfeatures = 1000)
FB_marker <- FindAllMarkers(FB_species_rpca_sub,min.pct = 0.2)
FB_marker_top10 <- FB_marker %>% dplyr::filter(pct.1>0.3)%>% dplyr::group_by(cluster) %>% dplyr::top_n(10,avg_log2FC)
FB_marker_top50 <- FB_marker %>% dplyr::filter(pct.1>0.3)%>% dplyr::group_by(cluster) %>% dplyr::top_n(50,avg_log2FC)
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

p <- DimPlot(EC_species_rpca_sub,label = T)
cells <- Seurat::CellSelector(p)
EC_species_rpca_sub <- subset(EC_species_rpca_sub,cells = cells,invert=T)
DimPlot(EC_species_rpca_sub,label = T)

EC_species_rpca_sub <- run_rpca(EC_species_rpca_sub,nfeatures = 2000)
label <- c('0'="capillary_ec",
           '1'="capillary_ec",
           '2'="vein_ec",
           '3'="artery_ec",
           '4'="lymphatic_ec",
           '5'="capillary_ec_INF",
           '6'='other1',
           '7'='other2')   
Idents(EC_species_rpca_sub) <- 'seurat_clusters'
EC_species_rpca_sub_label <- Seurat::RenameIdents(EC_species_rpca_sub,label)
EC_species_rpca_sub_label$subtype <- 
  DimPlot(EC_species_rpca_sub_label,label = T)
DimPlot(EC_species_rpca,label = T,split.by = 'species')

pan_ec <- c('PECAM1','VWF','CDH5')
cluste_gene <- c('PKHD1L1','PRKG1','TMEM108','CDH11','PCDH7','CGNL1')
imm_gene <- c('ISG15','IFIT1B','RSAD2','RTP4','IFIT3','GBP4','MNDA','SP140')
vein_ec	<- c('ACKR1','PLVAP')#2,4
capillary_ec <- c('RGCC','KDR','JUN','FOS','ATF3')#0,1
artery_ec	<- c('SEMA3G','GJA5','DLL4','EFNB2')#3
lymphatic_ec <- c('PROX1','NPR3','PKHD1L1','RELN') 




table(EC_species_rpca$species,EC_species_rpca$seurat_clusters)
table(EC_species_rpca$aged,EC_species_rpca$seurat_clusters)
DimPlot(EC_species_rpca,label = T)
DimPlot(EC_species_rpca,label = T,split.by = 'species')
FeaturePlot(EC_species_rpca,features = c('XDH','ADH1A','ADH1C',
                                         'ADH1B','TIMP4','FMO1',
                                         'MNDA','TCF15'))
FeaturePlot(EC_species_rpca,features = c('XDH','ADH1A','ADH1C',
                                         'ADH1B','TIMP4','FMO1',
                                         'MNDA','TCF15'))

EC_7 <- enrichGO(gene = dplyr::filter(EC_marker,cluster%in%'7'&pct.1>0.5&pct.2<0.3&avg_log2FC>1)$gene,
                 OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                 ont = "BP",keyType = "SYMBOL",
                 minGSSize=1,pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
dotplot(EC_7)

FeaturePlot(EC_species_rpca,
            split.by = 'species',
            slot = 'scale.data',
            features = c('EGFL7','PTPRB',
                         'CDH5','KDR',
                         'VWF','FABP4'))

FeaturePlot(EC_species_rpca,split.by = 'species',features = c('TCF15','PECAM1','CDH5','KDR'))


