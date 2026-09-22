# =============================================================================
# 02_per_species_rpca_integration.R
# Title : Within-species RPCA integration
# Figure: Preprocessing
# Module: 01_preprocessing
# Description:
#   Define run_rpca() and integrate multiple donors/samples within each species with Seurat v5 split-layer RPCA integration.
# Inputs:
#   - cr("human", "merge_monkey_sce_rpca.rds")
#   - cr("cross_FB", "res", "sce_mouse_injury_rpca.rds")
#   - cr("human", "crossspecies", "sceobj", "merge_mouse1_sce_rpca.rds")
#   - cr("human", "crossspecies", "sceobj", "merge_mouse3_sce_rpca")
#   - scatac("scRNA", "public", "CM_sce.rds")
# Outputs:
#   - (see script)
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

run_rpca <- function(sceobj){
  sceobj <- JoinLayers(sceobj)
  sceobj[["RNA"]] <- split(sceobj[["RNA"]], 
                           f = sceobj$orig.ident)
  sceobj <- NormalizeData(sceobj)
  sceobj <- FindVariableFeatures(sceobj,nfeatures = 2000)
  sceobj <- ScaleData(sceobj)
  sceobj <- RunPCA(sceobj)
  sceobj <- IntegrateLayers(
    object = sceobj, method = RPCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.rpca",
    verbose = FALSE)
  sceobj <- FindNeighbors(sceobj,
                          reduction = "integrated.rpca", dims = 1:30)
  sceobj <- FindClusters(sceobj,resolution = 0.4,
                         cluster.name = "rpca_clusters")
  sceobj <- RunUMAP(sceobj,
                    reduction = "integrated.rpca",
                    dims = 1:30, reduction.name = "umap.rpca")
  return(sceobj)
}

merge_monkey_sce_rpca <- read_rds(file = cr("human", "merge_monkey_sce_rpca.rds"))
sce_mouse_injury_rpca <- read_rds(file = cr("cross_FB", "res", "sce_mouse_injury_rpca.rds"))
merge_mouse1_sce_rpca <- read_rds(file = cr("human", "crossspecies", "sceobj", "merge_mouse1_sce_rpca.rds"))
merge_mouse3_sce_rpca <- read_rds(file = cr("human", "crossspecies", "sceobj", "merge_mouse3_sce_rpca"))
human_sce <- read_rds(file = scatac("scRNA", "public", "CM_sce.rds"))
human_sce_sub <- subset(human_sce,age%in%c('20-25','40-45','70-75'))
human_sce_sub[["RNA"]] <- as(object = human_sce_sub[["RNA"]], Class = "Assay5")

human_sce_FB <- subset(human_sce_sub,cell_type%in%'Fibroblast')
human_sce_FB$orig.ident <- human_sce_FB$donor
human_sce_FB[["RNA"]] <- as(object = human_sce_FB[["RNA"]], Class = "Assay5")

merge_monkey_sce_rpca_FB <- subset(merge_monkey_sce_rpca,seurat_clusters%in%c('3','13'))
sce_mouse_injury_rpca_FB <- subset(sce_mouse_injury_rpca,)
merge_mouse1_sce_rpca_FB <- subset(merge_mouse1_sce_rpca,seurat_clusters%in%c('0','2','16'))
merge_mouse3_sce_rpca_FB <- subset(merge_mouse3_sce_rpca,seurat_clusters%in%c('0','2','16'))

VlnPlot(merge_monkey_sce_rpca,features = c('DCN','GSN','LAMA2','BICC1','COL6A3'),pt.size = 0)
DimPlot(merge_monkey_sce_rpca,label = T)

VlnPlot(merge_mouse1_sce_rpca,features = c('Dcn','Gsn','Col6a3','Fn1','Meox1'),pt.size = 0)
DimPlot(merge_mouse1_sce_rpca,label = T)

VlnPlot(merge_mouse3_sce_rpca_FB,features = c('Dcn','Gsn','Col6a3','Fn1','Meox1'),pt.size = 0)
DimPlot(merge_mouse3_sce_rpca_FB,label = T)

human_sce_FB <- run_rpca(human_sce_FB)
merge_monkey_sce_rpca_FB <- run_rpca(merge_monkey_sce_rpca_FB)
merge_mouse1_sce_rpca_FB <- run_rpca(merge_mouse1_sce_rpca_FB)
merge_mouse3_sce_rpca_FB <- run_rpca(merge_mouse3_sce_rpca_FB)
