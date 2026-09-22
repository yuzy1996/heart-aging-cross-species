# =============================================================================
# 03_cross_species_integration.R
# Title : Cross-species integration and cell-type annotation
# Figure: Fig. 1B
# Module: 01_preprocessing
# Description:
#   Integrate the three species, assign the nine major cell types and write the core object combine_species_new_label.rds; marker dot plots and composition Sankey.
# Inputs:
#   - cr("human", "crossspecies", "sceobj", "combine_species_batch_noCM_V2_sub.rds")
#   - hu("res", "combine_species_new_label.rds")
#   - hu("res", "sceobj", "combine_species_rpca.rds")
# Outputs:
#   - hu("fig", "V1", "fig1_dot_marker.pdf")
#   - hu("res", "combine_species_new_label.rds")
#   - hu("res", "celltype_sankey.pdf")
#   - hu("res", "sceobj", "combine_species_rpca.rds")
#   - hu("res", "sceobj", "combine_species.rds")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

options(future.globals.maxSize = 10000* 1024^2)
# (R library search paths are managed by renv / .Renviron, not set per-script)
library(Seurat)
library(glue)
library(tidyverse)
library(SCP)
packageVersion('Seurat')

combine_species_batch_noCM_V2_sub <- read_rds(cr("human", "crossspecies", "sceobj", "combine_species_batch_noCM_V2_sub.rds"))
table(combine_species_batch_noCM_V2_sub$species)
table(combine_species_batch_noCM_V2_sub$celltypes)
combine_species_sub <- subset(combine_species_batch_noCM_V2_sub,
                              species%in%c('human','monkey','mouse1','mouse3'))
combine_species_sub <- JoinLayers(combine_species_sub)
combine_species_sub_exp <- LayerData(combine_species_sub, assay="RNA", layer='counts')
combine_species_new <- CreateSeuratObject(combine_species_sub_exp)
identical(row.names(combine_species_new),row.names(combine_species_sub))
combine_species_new@meta.data <- combine_species_sub@meta.data
combine_species_new[["RNA"]] <- split(combine_species_new[["RNA"]], 
                             f = combine_species_new$orig.ident)
combine_species_new <- NormalizeData(combine_species_new)
combine_species_new <- FindVariableFeatures(combine_species_new,nfeatures = 2000)
combine_species_new <- ScaleData(combine_species_new)
combine_species_new <- RunPCA(combine_species_new)
combine_species_new <- IntegrateLayers(
  object = combine_species_new, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  verbose = T)
combine_species_new <- FindNeighbors(combine_species_new,
                            reduction = "integrated.rpca", dims = 1:30)
combine_species_new <- FindClusters(combine_species_new,resolution = c(0.4,0.6,0.8,0.2))
combine_species_new <- RunUMAP(combine_species_new,
                      reduction = "integrated.rpca",
                      dims = 1:30, reduction.name = "umap.rpca")

DimPlot(combine_species_new,label = T,split.by = 'species')
DimPlot(combine_species_new,label = T,group.by = 'species')
DimPlot(combine_species_new,label = T)

FeaturePlot(combine_species_new,features = c('CD19','CD79A','MS4A1'),order = T)

combine_species_new <- JoinLayers(combine_species_new)
marker_cluster <- FindAllMarkers(combine_species_new,min.pct = 0.4)
marker_cluster_top10 <- marker_cluster %>% 
                          dplyr::group_by(cluster) %>% 
                           dplyr::top_n(10,wt = avg_log2FC)
marker_cluster_top30 <- marker_cluster %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(30,wt = avg_log2FC)
DoHeatmap(subset(combine_species_new,downsample=500),slot = 'data',
          features = marker_celltypes_top10$gene)

DimPlot(combine_species_new,group.by = 'celltypes',label = T)|

DimPlot(combine_species_new,group.by = 'seurat_clusters',label = T)

label <- c('0'="FB",
           '1'="EC",
           '2'="Myeloid",
           '3'="Pericytes",
           '4'="T",
           '5'="EC",
           '6'='Neural',
           '7'='SMC',
           '8'="FB",
           '9'="Myeloid",
           '10'="Adipocyte",
           '11'="Myeloid",
           '12'="EC",
           '13'='B',
           '14'='Myeloid',
           '15'="Myeloid")   
Idents(combine_species_new) <- 'seurat_clusters'
combine_species_new_label <- Seurat::RenameIdents(combine_species_new,label)
combine_species_new_label$annotation <- Idents(combine_species_new_label)
combine_species_new_label$annotation <- factor(combine_species_new_label$annotation,
                                               levels = c('FB','EC','Myeloid',
                                                          'T','B','Pericytes',
                                                          'SMC','Neural','Adipocyte'))
pdf()
DimPlot(combine_species_new_label,
        group.by = 'annotation',
        split.by  = 'species',label = T)
DimPlot(combine_species_new_label,
        group.by = 'annotation',
        label = T)
DimPlot(combine_species_new_label,
        group.by = 'annotation',
        split.by = 'aged',
        label = T)
pdf(hu("fig", "V1", "fig1_dot_marker.pdf"),width = 12,height = 8)
DotPlot(combine_species_new_label,group.by = 'annotation',
        features = c("DCN","COL1A1","GSN","FBLN1","TCF21",#FB
                     "CDH5","PECAM1","RAMP2","EMCN","TEK",#EC
                     "CD163","CD68","MS4A6A","CSF1R",#macrophage
                     "CD3D","CD3E","CD4","CD8A","IL7R",#T
                     "IGKC","MS4A1","CD19","CD79A",#B
                     "RGS5","ABCC9","KCNJ8",#Pericytes
                     "MYH11","TAGLN2","ACTA2",#SMC
                     "PLP1","NRXN1","NRXN3",#
                     "ADIPOQ","PLIN1","GPD1","PNPLA3","PCK1"))+#Adipocyte
  scale_color_viridis_c(option = 'D')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()


write_rds(combine_species_new_label,file = hu("res", "combine_species_new_label.rds"))
combine_species_new_label <- read_rds(hu("res", "combine_species_new_label.rds"))
library(RColorBrewer)
library(ggalluvial)
library(tidyverse)
new_meta <- combine_species_new_label@meta.data
cluster_counts_celltype <- new_meta %>%
  group_by(species, aged, annotation) %>%
  summarize(count = n(), .groups = 'drop')
# 确保数据格式适合绘制桑基图
cluster_counts_celltype$celltypes <- factor(cluster_counts_celltype$annotation,
                                      levels = celltypes_new_level)
# 创建一个颜色映射列表，用于标签颜色
cell_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
label_colors <- c(human = "#E69F00", monkey = "#56B4E9", mouse1 = "#009E73",
                   mouse3 = "#CC79A7", 
                  Aged = "#F0E442", Young = "#0072B2",
                  FB = "#D55E00",
                  EC=cell_colors[1],
                  Myeloid=cell_colors[2],
                  T=cell_colors[3],
                  B=cell_colors[4],
                  Pericytes=cell_colors[5],
                  SMC=cell_colors[6],
                  Neural=cell_colors[7],
                  Adipocyte=cell_colors[8])
pdf(hu("res", "celltype_sankey.pdf"),width = 10,height = 8)
ggplot(cluster_counts_celltype, aes(axis1 = species, axis2 = annotation,axis3 = aged, y = count)) +
  geom_alluvium(aes(fill = aged), width = 0.1) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.1, color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),color = 'black',
            size = 4, vjust = -0.5) +
  scale_x_discrete(limits = c("Species", "celltypes", "aged"), expand = c(0.15, 0.05)) +
  scale_fill_manual(values = label_colors) +
  scale_color_manual(values = label_colors) +
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

DimPlot(combine_species_new_label,label = T)
ht <- FeatureHeatmap(
  srt = combine_species_new,
  group.by = "CellType", features = marker_celltypes_top10$gene, 
  feature_split = marker_celltypes_top10$gene,
  feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 5, width = 4
)

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

table(combine_species_batch_noCM_V2_sub$orig.ident)

combine_species <- subset(combine_species_batch_noCM_V2_sub,
                          species%in%c('human','monkey','mouse2','mouse3'))
combine_species_rpca <- run_rpca(combine_species,nfeatures = 3000)
write_rds(combine_species_rpca,file = hu("res", "sceobj", "combine_species_rpca.rds"))
write_rds(combine_species,file = hu("res", "sceobj", "combine_species.rds"))
combine_species_rpca <- read_rds(file = hu("res", "sceobj", "combine_species_rpca.rds"))
DimPlot(combine_species_rpca,label = T)
DimPlot(combine_species_rpca,group.by = 'celltypes',label = T)

Idents(combine_species_rpca) <- 'celltypes'
combine_species_rpca <- JoinLayers(combine_species_rpca)
marker_all <- FindAllMarkers(combine_species_rpca,min.pct = 0.4)
marker_all_top10 <- marker_all %>% group_by(cluster) %>% top_n(10,avg_log2FC)
marker_all_top30 <- marker_all %>% group_by(cluster) %>% top_n(30,avg_log2FC)

