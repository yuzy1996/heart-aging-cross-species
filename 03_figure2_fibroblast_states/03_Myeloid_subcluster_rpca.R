# =============================================================================
# 03_Myeloid_subcluster_rpca.R
# Title : Myeloid sub-clustering
# Figure: Supplementary
# Module: 03_figure2_fibroblast_states
# Description:
#   Subset, integrate, cluster and annotate myeloid subtypes (DC, neutrophil, IFN, CCR2+ etc.).
# Inputs:
#   - hu("data", "Myeloid_species_rpca_sub_label.rds")
# Outputs:
#   - (see script)
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

Myeloid_species <- subset(combine_species_rpca,celltypes%in%'Myeloid')
Myeloid_species_rpca <- run_rpca(Myeloid_species,nfeatures = 00)
Myeloid_marker <- FindAllMarkers(Myeloid_species_rpca)
Myeloid_marker_top10 <- Myeloid_marker %>% filter(pct.1>0.5) %>% dplyr::group_by(cluster) %>% dplyr::top_n(10,wt = avg_log2FC)
Myeloid_marker_top30 <- Myeloid_marker %>% filter(pct.1>0.5) %>% dplyr::group_by(cluster) %>% dplyr::top_n(30,wt = avg_log2FC)
# plot marker -------------------------------------------------------------
DimPlot(Myeloid_species_rpca,label = T)
Myeloid_species_rpca_sub <- subset(Myeloid_species_rpca,seurat_clusters%in%c(0:5,7:8,11:13))
DC <- c('XCR1','CLEC9A','SEPTIN3')#
Neutrophil <- c('S100A8','CXCR2','CSF3R')#
INF <- c('IFIT3','ISG15','IFIT2','IRF7')#
Ccr2_marker <- c('Ccr2','Cd74','H2-Aa','H2-Ab1','H2-Eb1','H2-DMb1')
Timd4_marker <- c('TIMD4','LYVE1','CD163')
MHC_marker <- c("HLA-DQB2",'HLA-DMA','H3-3A','CXCL16')
Mono <- c('FABP4','CAV1','CD36','RGCC','ID3')
Cycle <- c('MKI67','TOP2A','SMC4','HMGB2')

FeaturePlot(Myeloid_species_rpca_sub,features = DC)
FeaturePlot(Myeloid_species_rpca_sub,features = INF)
FeaturePlot(Myeloid_species_rpca_sub,features = Neutrophil)
FeaturePlot(Myeloid_species_rpca_sub,features = MHC_marker)
FeaturePlot(Myeloid_species_rpca_sub,features = Timd4_marker)
FeaturePlot(Myeloid_species_rpca_sub,features = Mono)
FeaturePlot(Myeloid_species_rpca_sub,features = c('PECAM1','VWF','CDH5','RGCC','KDR'),order = T)
FeaturePlot(Myeloid_species_rpca_sub,features = 'CCR2',order = T)
FeaturePlot(Myeloid_species_rpca_sub,features = c('CD163','NEAT1','CTSB','MALAT1','CCL2'))
FeaturePlot(Myeloid_species_rpca_sub,features = c('VCAN','PTPRC','CTSB'),order = T)

DimPlot(Myeloid_species_rpca_sub,label = T)

DimPlot(Myeloid_species_rpca_sub,label = T)

EC_species_rpca_sub <- run_rpca(EC_species_rpca_sub,nfeatures = 2000)
label <- c('0'="Timd4_resident_MC",
           '1'="MHClow_MC",
           '2'="CCR2_MC",
           '3'="MHChi_MC",
           '4'="MHClow_MC",
           '5'="Mono",
           '7'='Cycle',
           '8'="MHClow_MC",
           '11'="INF_MC",
           '12'='Neutrophil',
           '13'='DC')   
Idents(Myeloid_species_rpca_sub) <- 'seurat_clusters'
Myeloid_species_rpca_sub_label <- Seurat::RenameIdents(Myeloid_species_rpca_sub,label)
write_rds(Myeloid_species_rpca_sub_label,
          file = hu("data", "Myeloid_species_rpca_sub_label.rds"))
Myeloid_species_rpca_sub_label$subtype <- Idents(Myeloid_species_rpca_sub_label)
DimPlot(Myeloid_species_rpca_sub_label,label = T)
DimPlot(Myeloid_species_rpca_sub_label,label = T,split.by = 'species')

pan_ec <- c('PECAM1','VWF','CDH5')
cluste_gene <- c('PKHD1L1','PRKG1','TMEM108','CDH11','PCDH7','CGNL1')
imm_gene <- c('ISG15','IFIT1B','RSAD2','RTP4','IFIT3','GBP4','MNDA','SP140')
vein_ec	<- c('ACKR1','PLVAP')#2,4
capillary_ec <- c('RGCC','KDR','JUN','FOS','ATF3')#0,1
artery_ec	<- c('SEMA3G','GJA5','DLL4','EFNB2')#3
lymphatic_ec <- c('PROX1','NPR3','PKHD1L1','RELN') 

FeaturePlot(EC_species_rpca_sub,features = c('TBX1','PDPN','NPR3','CX3CL1','PROX1'),order = T)
FeaturePlot(EC_species_rpca_sub,features = c('ADH1A','ADH1C','ADH1B'))
FeaturePlot(EC_species_rpca_sub,features = c('LAMC1','LHFPL6','PAM','PLXDC2'))


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
Myeloid_species_rpca_sub_label <- read_rds(hu("data", "Myeloid_species_rpca_sub_label.rds"))
