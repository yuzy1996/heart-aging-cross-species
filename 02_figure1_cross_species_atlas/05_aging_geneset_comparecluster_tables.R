# =============================================================================
# 05_aging_geneset_comparecluster_tables.R
# Title : Per-species aging gene-set enrichment tables
# Figure: Fig. 1D
# Module: 02_figure1_cross_species_atlas
# Description:
#   Reshape compareCluster enrichment results for human/macaque/mouse into cell-type-ordered tables.
# Inputs:
#   - hu("data", "human.csv")
#   - hu("data", "human1.csv")
#   - hu("data", "human2.csv")
#   - hu("data", "human3.csv")
#   - hu("data", "human4.csv")
#   - hu("data", "human5.csv")
#   - hu("data", "mouse1.csv")
#   - hu("data", "mouse2.csv")
#   - hu("data", "mouse3.csv")
#   - hu("data", "mouse4.csv")
# Outputs:
#   - (see script)
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))


human_enrich_res <- data.frame(human_enrich_down@compareClusterResult)
monkey_enrich_res <- data.frame(monkey_enrich_down@compareClusterResult)
mouse_enrich_res <- data.frame(mouse_enrich_down@compareClusterResult)

cell_factor <- c('FB','EC','Myeloid','Pericytes','SMC','T','B')
human_enrich_res <- filter(human_enrich_res,Cluster%in%cell_factor)
monkey_enrich_res <- filter(monkey_enrich_res,Cluster%in%cell_factor)
mouse_enrich_res <- filter(mouse_enrich_res,Cluster%in%cell_factor)

human_enrich_res$Cluster <- factor(human_enrich_res$Cluster, cell_factor)
monkey_enrich_res$Cluster <- factor(monkey_enrich_res$Cluster, cell_factor)
mouse_enrich_res$Cluster <- factor(mouse_enrich_res$Cluster, cell_factor)
plot_dot(human_enrich_res,n = 3)/
  plot_dot(monkey_enrich_res,n = 3)/
  plot_dot(mouse_enrich_res,n = 3)

human_ageset <- read.csv(hu("data", "human.csv"))
human_ageset1 <- read.csv(hu("data", "human1.csv"))
human_ageset2 <- read.csv(hu("data", "human2.csv"))
human_ageset3 <- read.csv(hu("data", "human3.csv"))
human_ageset4 <- read.csv(hu("data", "human4.csv"))
human_ageset5 <- read.csv(hu("data", "human5.csv"))
human_ageset_all <- rbind(human_ageset,
                          human_ageset1,
                          human_ageset2,
                          human_ageset3,
                          human_ageset4,
                          human_ageset5)
mouse_ageset1 <- read.csv(hu("data", "mouse1.csv"))
mouse_ageset2 <- read.csv(hu("data", "mouse2.csv"))
mouse_ageset3 <- read.csv(hu("data", "mouse3.csv"))
mouse_ageset4 <- read.csv(hu("data", "mouse4.csv"))
mouse_ageset_all <- rbind(mouse_ageset1,
                          mouse_ageset2,
                          mouse_ageset3,
                          mouse_ageset4)
head(human_ageset_all)
human_ageset_all$X...Symbol
table(human_ageset_all$Gene_Set)
head(mouse_ageset_all)
mouse_ageset_all$X...Symbol

cellage <- read_tsv(lmx("GeneSet_aging", "cellAge", "cellage3.tsv"))
SAUL_SEN_MAYO <- read.gmt(lmx("GeneSet_aging", "SAUL_SEN_MAYO.v2024.1.Hs.gmt"))
REACTOME_SASP <- read.gmt(lmx("GeneSet_aging", "SASP.gmt"))
