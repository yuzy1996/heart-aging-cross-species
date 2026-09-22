# =============================================================================
# 04_cluster_to_celltype_label_map.R
# Title : Map Seurat clusters to cell-type labels
# Figure: Fig. 1B
# Module: 01_preprocessing
# Description:
#   Explicit mapping from unsupervised Seurat clusters to the nine major cell-type labels.
# Inputs:
#   - (see script)
# Outputs:
#   - (see script)
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

combine_species_ann <- combine_species_new_label
label <- c('0'="FB",
           '1'="EC",
           '2'="Myeloid",
           '3'="PericytesFB",
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
           '15'="Myeloid",
           '16'="Myeloid",
           '17'="Myeloid",
           '18'="Myeloid",
           '19'="Myeloid",
           '20'="Myeloid",
           '21'="Myeloid",
           '22'="Myeloid")   
Idents(combine_species_ann) <- 'RNA_snn_res.0.4'
combine_species_new_label <- Seurat::RenameIdents(combine_species_ann,label)
combine_species_new_label$annotation <- Idents(combine_species_new_label)
combine_species_new_label$annotation <- factor(combine_species_new_label$annotation,
                                               levels = c('FB','EC','Myeloid',
                                                          'T','B','Pericytes',
                                                          'SMC','Neural','Adipocyte'))

