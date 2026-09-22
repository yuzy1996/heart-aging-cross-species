# =============================================================================
# 07_fig1C_geneset_group_means.R
# Title : Gene-set enrichment by species and age
# Figure: Fig. 1C
# Module: 02_figure1_cross_species_atlas
# Description:
#   Group-mean aging gene-set enrichment scores for young vs aged hearts of each species.
# Inputs:
#   - hu("res", "combine_species_new_label.rds")
#   - hu("Aging_hallmarker.rda")
# Outputs:
#   - (see script)
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

library(tidyverse)
library(SCP)
library(SeuratWrappers)
library(Seurat)
library(SeuratData)
celltypes_col <- c("#FF1493", "#FFB6C1", "#483D8B", "#00BFFF", "#228B22", "#7FFF00", "#FFFF00", "#FF4500","#8DD3C7")
species_col <- c("#FC8D62","#8DD3C7","#BC80BD")
aged_coll <- c("#483D8B","#00BFFF")
celltypes <- c('FB','EC','Myeloid','T','B','Pericytes','SMC','Neural','Adipocyte')
species <- c('human','monkey','mouse')
aged <- c('Aged','Young')
combine_species_new_label <- read_rds(hu("res", "combine_species_new_label.rds"))
combine_species_new_label_S4 <- combine_species_new_label
combine_species_new_label_S4[['RNA']] <- as(combine_species_new_label_S4[['RNA']],Class = 'Assay')
load(hu("Aging_hallmarker.rda"))
combine_species_new_label <- AddModuleScore(combine_species_new_label,
                                            features = Aging_hallmarker,
                                            name = names(Aging_hallmarker))
aging_term <- combine_species_new_label@meta.data[,c('orig.ident',
                                                     'aged','species',
                                                     'SASP gene set11',
                                                     'Inflammation-related genes12',
                                                     'Cardiac fibrosis-related genes14',
                                                     'ECM15',
                                                     'genomic_instability2',
                                                     'loss_of_proteostasis10')]
head(aging_term)
result <- plot_sample_module_scores(
  aging_term,
  sample_col  = "orig.ident",
  age_col     = "aged",
  species_col = "species"
)

# 分别显示三个物种
print(result$plots_by_species$Human)
print(result$plots_by_species$Macaque)
print(result$plots_by_species$Mouse)

# 总图：三行物种 × 六列模块
print(result$plot)

# 查看实际 P 值、校正后 q 值及有效样本数
result$statistics %>%
  dplyr::select(
    Species, Module,
    n_young, n_aged,
    mean_difference_aged_minus_young,
    P, BH_q, method, status
  )