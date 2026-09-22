# =============================================================================
# 01_load_per_species_objects.R
# Title : Load per-species Seurat objects
# Figure: Preprocessing
# Module: 01_preprocessing
# Description:
#   Read the upstream human, macaque and mouse single-nucleus Seurat objects and subset age groups before cross-species integration.
# Inputs:
#   - scatac("scRNA", "public", "CM_sce.rds")
#   - cr("human", "merge_monkey_sce_rpca.rds")
#   - cr("human", "crossspecies", "sceobj", "merge_mouse1_sce_rpca.rds")
#   - cr("human", "crossspecies", "sceobj", "merge_mouse3_sce_rpca")
#   - cr("human", "crossspecies", "sceobj", "merge_mouse2_sce.rds")
# Outputs:
#   - (see script)
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

library(Seurat)
packageVersion('Seurat')

human_sce <- read_rds(file = scatac("scRNA", "public", "CM_sce.rds"))
human_sce_sub <- subset(human_sce,age%in%c('20-25','40-45','70-75'))
human_sce_sub[["RNA"]] <- as(object = human_sce_sub[["RNA"]], Class = "Assay5")


merge_monkey_sce_rpca <- read_rds(file = cr("human", "merge_monkey_sce_rpca.rds"))
#sce_mouse_injury_rpca <- read_rds(file = cr("cross_FB", "res", "sce_mouse_injury_rpca.rds"))
merge_mouse1_sce_rpca <- read_rds(file = cr("human", "crossspecies", "sceobj", "merge_mouse1_sce_rpca.rds"))
merge_mouse3_sce_rpca <- read_rds(file = cr("human", "crossspecies", "sceobj", "merge_mouse3_sce_rpca"))
merge_mouse2_sce_rpca <- read_rds(file = cr("human", "crossspecies", "sceobj", "merge_mouse2_sce.rds"))

merge_monkey_sce_rpca <- JoinLayers(merge_monkey_sce_rpca)
monkey_cluster <- FindAllMarkers(merge_monkey_sce_rpca,min.pct = 0.4)

merge_mouse1_sce_rpca <- JoinLayers(merge_mouse1_sce_rpca)
mouse1_cluster <- FindAllMarkers(merge_mouse1_sce_rpca,min.pct = 0.4)

merge_mouse3_sce_rpca <- JoinLayers(merge_mouse3_sce_rpca)
mouse3_cluster <- FindAllMarkers(merge_mouse3_sce_rpca,min.pct = 0.4)

row.names(merge_monkey_sce_rpca)


merge_mouse1_sce_rpca <- JoinLayers(merge_mouse1_sce_rpca)
merge_mouse1_sce_rpca[["RNA"]] <- as(object = merge_mouse1_sce_rpca[["RNA"]], Class = "Assay")

merge_mouse1_sce_rpca <- JoinLayers(merge_mouse1_sce_rpca)
merge_mouse1_sce_rpca[["RNA"]] <- as(object = merge_mouse1_sce_rpca[["RNA"]], Class = "Assay")

write_rds()