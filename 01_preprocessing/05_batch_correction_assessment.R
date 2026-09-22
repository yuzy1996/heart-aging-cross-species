# =============================================================================
# 05_batch_correction_assessment.R
# Title : Batch-correction assessment (scVI/scANVI, CCA/MNN/Harmony/RPCA)
# Figure: Supplementary (integration QC)
# Module: 01_preprocessing
# Description:
#   Convert the integrated object to h5ad, read externally computed scVI/scANVI UMAP coordinates and compare integration methods for batch mixing.
# Inputs:
#   - hu("FIG1", "X_scVI_umap.csv")
#   - hu("FIG1", "X_scANVI_umap.csv")
#   - hu("fig", "V2", "combine_species_cca.rds")
#   - hu("fig", "V2", "combine_species_mnn.rds")
#   - hu("fig", "V2", "combine_species_harmony.rds")
#   - hu("res", "combine_species_new_label.rds")
#   - hu("FIG1", "combine_species_batch_access.rds")
# Outputs:
#   - hu("FIG1", "batch_umap.pdf")
#   - hu("FIG1", "combine_species_batch_access.rds")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

SeuratDisk::SaveH5Seurat(combine_species_new_label_S4_, 
                         filename = hu("FIG1", "combine_species_S4.h5Seurat"))
SeuratDisk::Convert(hu("FIG1", "combine_species_S4.h5Seurat"), dest = "h5ad")
combine_species_new_label_S4_ <- SeuratDisk::LoadH5Seurat(hu("FIG1", "combine_species_S4.h5Seurat"))
umap_scvi_results <- read.csv(hu("FIG1", "X_scVI_umap.csv"), row.names = 1)
umap_scanvi_results <- read.csv(hu("FIG1", "X_scANVI_umap.csv"), row.names = 1)
library(Seurat)
library(tidyverse)
combine_species_cca <- read_rds(file = hu("fig", "V2", "combine_species_cca.rds"))
combine_species_mnn <- read_rds(file = hu("fig", "V2", "combine_species_mnn.rds"))
combine_species_harmony <- read_rds(file = hu("fig", "V2", "combine_species_harmony.rds"))
# 确保细胞顺序匹配
umap_scvi_results <- umap_scvi_results[match(Cells(combine_species_new_label_S4_),
                                             rownames(umap_scvi_results)), ]
umap_scanvi_results <- umap_scanvi_results[match(Cells(combine_species_new_label_S4_),
                                             rownames(umap_scanvi_results)), ]
# 将UMAP结果添加到Seurat对象中
combine_species_new_label_S4_[['SCVI_umap']] <- CreateDimReducObject(embeddings = as.matrix(umap_scvi_results),
                                                                     key = 'UMAP_')
combine_species_new_label_S4_[['SCANVI_umap']] <- CreateDimReducObject(embeddings = as.matrix(umap_scanvi_results), 
                                                                       key = 'UMAP_')
combine_species_new_label_S4_ <- read_rds(hu("res", "combine_species_new_label.rds"))
# 查看UMAP结果
head(combine_species_new_label_S4_@reductions$SCVANI_umap@cell.embeddings)
pdf(hu("FIG1", "batch_umap.pdf"),width = 5,height = 4)
DimPlot(combine_species_new_label_S4_,reduction = 'SCANVI_umap',group.by = 'species')
DimPlot(combine_species_new_label_S4_,reduction = 'SCVI_umap',group.by = 'species')
DimPlot(combine_species_new_label_S4_,reduction = 'umap.rpca',group.by = 'species')
DimPlot(combine_species_cca,reduction = 'umap.cca',group.by = 'species')
DimPlot(combine_species_mnn,reduction = 'umap.mnn',group.by = 'species')
DimPlot(combine_species_harmony,reduction = 'umap.harmony',group.by = 'species')
dev.off()
combine_species_batch_access <- combine_species_new_label_S4_
combine_species_batch_access@reductions$umap.cca <- combine_species_cca@reductions$umap.cca
combine_species_batch_access@reductions$umap.mnn <- combine_species_mnn@reductions$umap.mnn
combine_species_batch_access@reductions$umap.harmony <- combine_species_harmony@reductions$umap.harmony
combine_species_batch_access@reductions$umap.unintegrated <- combine_species_cca@reductions$umap.unintegrated

combine_species_batch_access@reductions$integrated.rpca <- combine_species_cca@reductions$integrated.rpca
combine_species_batch_access@reductions$integrated.cca <- combine_species_cca@reductions$integrated.cca
combine_species_batch_access@reductions$integrated.mnn <- combine_species_mnn@reductions$integrated.mnn
combine_species_batch_access@reductions$harmony <- combine_species_harmony@reductions$harmony
combine_species_batch_access@reductions$pca <- combine_species_cca@reductions$pca


pdf(hu("FIG1", "batch_umap.pdf"),width = 5,height = 4)
DimPlot(combine_species_batch_access,reduction = 'umap.unintegrated',group.by = 'species')
DimPlot(combine_species_batch_access,reduction = 'SCANVI_umap',group.by = 'species')
DimPlot(combine_species_batch_access,reduction = 'SCVI_umap',group.by = 'species')
DimPlot(combine_species_batch_access,reduction = 'umap.rpca',group.by = 'species')
DimPlot(combine_species_batch_access,reduction = 'umap.cca',group.by = 'species')
DimPlot(combine_species_batch_access,reduction = 'umap.mnn',group.by = 'species')
DimPlot(combine_species_batch_access,reduction = 'umap.harmony',group.by = 'species')
dev.off()
write_rds(combine_species_batch_access,file=hu("FIG1", "combine_species_batch_access.rds"))
combine_species_batch_access <- read_rds(hu("FIG1", "combine_species_batch_access.rds"))
combine_species_batch_access[['RNA']] <- as(combine_species_batch_access[['RNA']],Class='Assay')
combine_species_batch_access <- read_rds(file=hu("FIG1", "combine_species_batch_access.rds"))
SeuratDisk::SaveH5Seurat(combine_species_batch_access, 
                         filename = hu("FIG1", "combine_species_batch.h5Seurat"))
SeuratDisk::Convert(hu("FIG1", "combine_species_batch.h5Seurat"), dest = "h5ad")

