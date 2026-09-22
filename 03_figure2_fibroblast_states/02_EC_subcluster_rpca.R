# =============================================================================
# 02_EC_subcluster_rpca.R
# Title : Endothelial sub-clustering
# Figure: Supplementary
# Module: 03_figure2_fibroblast_states
# Description:
#   Subset, integrate, cluster and annotate endothelial subtypes.
# Inputs:
#   - hu("data", "EC_species_rpca_sub_label.rds")
# Outputs:
#   - (see script)
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)
EC_species <- subset(combine_species_rpca,celltypes%in%'EC')
EC_species_rpca <- run_rpca(EC_species,nfeatures = 2000)

EC_species_rpca <- FindClusters(EC_species_rpca,resolution = 0.2)
EC_marker <- FindAllMarkers(EC_species_rpca,min.pct = 0.2)
EC_marker_top10 <- EC_marker %>% filter(pct.1>0.5)%>% group_by(cluster) %>% top_n(10,avg_log2FC)
EC_marker_top50 <- EC_marker %>% dplyr::filter(pct.1>0.5)%>% dplyr::group_by(cluster) %>% dplyr::top_n(50,avg_log2FC)

# plot marker -------------------------------------------------------------
DimPlot(EC_species_rpca,label = T)
EC_species_rpca_sub <- subset(EC_species_rpca,seurat_clusters%in%c(0:7))

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
EC_species_rpca_sub_label$subtype <- Idents(EC_species_rpca_sub_label)
DimPlot(EC_species_rpca_sub_label,label = T)
DimPlot(EC_species_rpca,label = T,split.by = 'species')

pan_ec <- c('PECAM1','VWF','CDH5')
cluste_gene <- c('PKHD1L1','PRKG1','TMEM108','CDH11','PCDH7','CGNL1')
imm_gene <- c('ISG15','IFIT1B','RSAD2','RTP4','IFIT3','GBP4','MNDA','SP140')
vein_ec	<- c('ACKR1','PLVAP')#2,4
capillary_ec <- c('RGCC','KDR','JUN','FOS','ATF3')#0,1
artery_ec	<- c('SEMA3G','GJA5','DLL4','EFNB2')#3
lymphatic_ec <- c('PROX1','NPR3','PKHD1L1','RELN') 

FeaturePlot(EC_species_rpca_sub,features = capillary_ec)
FeaturePlot(EC_species_rpca_sub,features = pan_ec)
FeaturePlot(EC_species_rpca_sub,features = artery_ec)
FeaturePlot(EC_species_rpca_sub,features = vein_ec,order = T)
FeaturePlot(EC_species_rpca_sub,features = lymphatic_ec,order = T)
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
write_rds(EC_species_rpca_sub_label,
          file = hu("data", "EC_species_rpca_sub_label.rds"))
EC_species_rpca_sub_label <- read_rds(hu("data", "EC_species_rpca_sub_label.rds"))
DimPlot(EC_species_rpca_sub_label,label = T)
EC_species_final <- subset(EC_species_rpca_sub_label,
                           subtype%in%c( 'capillary_ec',
                                         'vein_ec',
                                         'artery_ec',
                                         'lymphatic_ec',
                                         'capillary_ec_INF'))

DimPlot(EC_species_final)

write_rds(EC_species_final,
          file = hu("data", "EC_species_final.rds"))


