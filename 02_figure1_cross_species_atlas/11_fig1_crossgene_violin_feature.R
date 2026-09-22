# =============================================================================
# 11_fig1_crossgene_violin_feature.R
# Title : Cross-species gene violin/feature plots
# Figure: Fig. 1F; Supplementary
# Module: 02_figure1_cross_species_atlas
# Description:
#   Violin/feature plots of conserved aging genes split by species and age.
# Inputs:
#   - (see script)
# Outputs:
#   - hu("fig", "V1", "fig1_pheatmap.pdf")
#   - hu("fig", "V1", "fig1_pheatmap1.pdf")
#   - hu("fig", "fig1_feature.pdf")
#   - hu("fig", "fig1_feature_human.pdf")
#   - hu("fig", "fig1_feature_monkey.pdf")
#   - hu("fig", "fig1_feature_mouse.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

Idents(combine_species_new_label) <- 'new_group_aged'
combine_species_new_label$aged <- factor(combine_species_new_label$aged,
                                         levels = c('Young','Aged'))
celltype_human <- subset(combine_species_new_label,dataset%in%'human')
celltypes <- c("FB","EC","Myeloid","T","B","Pericytes","SMC","Neural","Adipocyte")
human_celltype <- c("human_Young_Adipocyte","human_Aged_Adipocyte",
                    "human_Young_B","human_Aged_B" ,
                    "human_Young_EC","human_Aged_EC",
                    "human_Young_FB","human_Aged_FB",
                    "human_Young_Myeloid","human_Aged_Myeloid",
                    "human_Young_Neural","human_Aged_Neural",
                    "human_Young_Pericytes","human_Aged_Pericytes",
                    "human_Young_SMC","human_Aged_SMC",
                    "human_Young_T","human_Aged_T")
celltype_human$new_group_aged <- factor(as.character(celltype_human$new_group_aged),
                                        levels=human_celltype)
celltype_human$new_group_aged <- factor(as.character(celltype_human$new_group_aged),
                                        levels = c(glue('human_Young_{celltypes}'),
                                                   glue('human_Aged_{celltypes}')))
celltype_monkey <- subset(combine_species_new_label,dataset%in%'monkey')
celltype_monkey$new_group_aged <- factor(celltype_monkey$new_group_aged,
                                        levels = c(glue('monkey_Young_{celltypes}'),
                                                   glue('monkey_Aged_{celltypes}')))
celltype_mouse1 <- subset(combine_species_new_label,dataset%in%'mouse1')
celltype_mouse1$new_group_aged <- factor(celltype_mouse1$new_group_aged,
                                        levels = c(glue('mouse1_Young_{celltypes}'),
                                                   glue('mouse1_Aged_{celltypes}')))

VlnPlot(celltype_human,features = genes_conserved_up,pt.size = 0)
VlnPlot(celltype_monkey,features = genes_conserved_up,pt.size = 0)
VlnPlot(celltype_mouse1,features = genes_conserved_up,pt.size = 0)

VlnPlot(celltype_human,features = genes_conserved_down,pt.size = 0)
VlnPlot(celltype_monkey,features = genes_conserved_down,pt.size = 0)
VlnPlot(celltype_mouse1,features = genes_conserved_down,pt.size = 0)

combine_species_new_label@meta.data$new_group_aged <- paste(combine_species_new_label$dataset,
                                                            combine_species_new_label$aged,
                                                            combine_species_new_label$annotation,sep = '_')
combine_species_new_label@meta.data$new_group_aged  <- as.factor(combine_species_new_label@meta.data$new_group_aged)
Idents(combine_species_new_label) <- 'new_group_aged'
Idents(combine_species_new_label) <- 'annotation'
library(dittoSeq)


pdf(hu("fig", "V1", "fig1_pheatmap.pdf"),
    width = 12,height = 8)
dittoHeatmap(subset(celltype_human,downsample=1000),
             genes = c(genes_conserved_up,genes_conserved_down),
             annot.by = c("new_group_aged",'aged'))
dittoHeatmap(subset(celltype_monkey,downsample=1000),
             genes = c(genes_conserved_up,genes_conserved_down),
             annot.by = c("new_group_aged",'aged'))
dittoHeatmap(subset(celltype_mouse1,downsample=1000),
             genes = c(genes_conserved_up,genes_conserved_down),
             annot.by = c("new_group_aged",'aged'))
dev.off()
pdf(hu("fig", "V1", "fig1_pheatmap1.pdf"),
    width = 12,height = 8)
dittoHeatmap(celltype_human,
             genes = c(genes_conserved_up,genes_conserved_down),
             annot.by = c("new_group_aged",'aged'))
dittoHeatmap(celltype_monkey,
             genes = c(genes_conserved_up,genes_conserved_down),
             annot.by = c("new_group_aged",'aged'))
dittoHeatmap(celltype_mouse1,
             genes = c(genes_conserved_up,genes_conserved_down),
             annot.by = c("new_group_aged",'aged'))
dev.off()

celltype_human_S4 <- celltype_human
celltype_human_S4[['RNA']] <- as(celltype_human[['RNA']],Class = 'Assay')
celltype_monkey_S4 <- celltype_monkey
celltype_monkey_S4[['RNA']] <- as(celltype_monkey[['RNA']],Class = 'Assay')
celltype_mouse1_S4 <- celltype_mouse1
celltype_mouse1_S4[['RNA']] <- as(celltype_mouse1[['RNA']],Class = 'Assay')

SCP::FeatureStatPlot(celltype_monkey_S4, stat.by = genes_conserved_up, group.by = "new_group_aged", plot_type = "col")


combine_species_new_label_S4 <- combine_species_new_label
combine_species_new_label_S4[['RNA']] <- as(combine_species_new_label_S4[['RNA']],Class = 'Assay')
combine_species_new_label_S4@meta.data$group_species <- paste(combine_species_new_label_S4$species,
                                                    combine_species_new_label_S4$aged,sep = '_')

pdf(hu("fig", "fig1_feature.pdf"),width = 10,height = 12)
SCP::FeatureStatPlot(combine_species_new_label_S4, stat.by = genes_conserved_up, group.by = "group_species", plot_type = "col")
dev.off()

pdf(hu("fig", "fig1_feature_human.pdf"),width = 12,height = 8)
SCP::FeatureStatPlot(subset(celltype_human_S4,downsample=1000), 
                     stat.by = genes_conserved_up[1], 
                     group.by = "new_group_aged", 
                     plot_type = "col")
SCP::FeatureStatPlot(subset(celltype_human_S4,downsample=1000), 
                     stat.by = genes_conserved_up[2], 
                     group.by = "new_group_aged", 
                     plot_type = "col")
dev.off()
pdf(hu("fig", "fig1_feature_monkey.pdf"),width = 12,height = 8)
SCP::FeatureStatPlot(subset(celltype_monkey_S4,downsample=1000), 
                     stat.by = genes_conserved_up, 
                     group.by = "new_group_aged", 
                     plot_type = "col")
dev.off()
pdf(hu("fig", "fig1_feature_mouse.pdf"),width = 12,height = 8)
SCP::FeatureStatPlot(subset(celltype_monkey_S4,downsample=1000), 
                     stat.by = genes_conserved_up, 
                     group.by = "new_group_aged", 
                     plot_type = "col")
dev.off()


