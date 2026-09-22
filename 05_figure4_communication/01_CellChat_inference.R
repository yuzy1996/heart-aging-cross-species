# =============================================================================
# 01_CellChat_inference.R
# Title : CellChat ligand-receptor inference
# Figure: Fig. 4A-F
# Module: 05_figure4_communication
# Description:
#   Run CellChat per species and age; infer interaction counts and pathway-level communication networks.
# Inputs:
#   - hu("res", "combine_species_new_label1.rds")
#   - hu("FIG5", "sub.rds")
#   - hu("FIG2", "subcluster", "EC_recluster_label_S4.rds")
#   - hu("FIG2", "subcluster", "Myeloid_recluster_label.rds")
#   - hu("FIG2", "subcluster", "FB_label.rds")
#   - hu("FIG5", "allcelltypes.rds")
#   - hu("FIG5", "allcelltypes_monkey_Y1.rds")
#   - hu("FIG5", "allcelltypes_monkey_O1.rds")
#   - hu("FIG5", "allcelltypes_mouse_Y.rds")
#   - hu("FIG5", "allcelltypes_mouse_O.rds")
# Outputs:
#   - hu("FIG5", "sub.rds")
#   - hu("FIG5", "allcelltypes.rds")
#   - hu("FIG5", "allcelltypes_human_Y_obj1.rds")
#   - hu("FIG5", "allcelltypes_human_O_obj1.rds")
#   - hu("FIG5", "allcelltypes_monkey_Y.rds")
#   - hu("FIG5", "allcelltypes_monkey_O.rds")
#   - hu("FIG5", "allcelltypes_monkey_Y_obj1.rds")
#   - hu("FIG5", "allcelltypes_monkey_O_obj1.rds")
#   - hu("FIG5", "allcelltypes_mouse_Y.rds")
#   - hu("FIG5", "allcelltypes_mouse_O.rds")
#   - hu("FIG5", "allcelltypes_mouse_Y_obj3_1.rds")
#   - hu("FIG5", "allcelltypes_mouse_O_obj3_1.rds")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

# (R library search paths are managed by renv / .Renviron, not set per-script)
library(Seurat)
library(tidyverse)
library(CellChat)
CellChatDB <- CellChatDB.human # set CellChatDB <- CellChatDB.human if working on the human dataset
interaction_input <- CellChatDB$interaction
complex_input <- CellChatDB$complex
cofactor_input <- CellChatDB$cofactor
geneInfo <- CellChatDB$geneInfo

allcelltype <- readRDS(hu("res", "combine_species_new_label1.rds"))
allcelltype[['RNA']] <- as(allcelltype[['RNA']],Class = 'Assay')
sub <- subset(allcelltype,annotation%in%c('T','B','Pericytes','SMC','Neural','Adipocyte'))
sub[['RNA']] <- as(sub[['RNA']],Class = 'Assay')
write_rds(sub,hu("FIG5", "sub.rds"))

table(EC_sub$subtype)
table(Myeloid_sub$subtype)
table(FB_sub$subtype)
write_rds(allcelltypes,hu("FIG5", "allcelltypes.rds"))
sub <- readRDS(hu("FIG5", "sub.rds"))

EC_sub <- readRDS(file=hu("FIG2", "subcluster", "EC_recluster_label_S4.rds"))
Myeloid_sub <- readRDS(hu("FIG2", "subcluster", "Myeloid_recluster_label.rds"))
FB_sub <- readRDS(hu("FIG2", "subcluster", "FB_label.rds"))
sub@meta.data$subtype <- sub$annotation
allcelltypes <- merge(sub,list(EC_sub,
                               Myeloid_sub,
                               FB_sub))
table(allcelltypes$species,allcelltypes$dataset)

allcelltypes <- readRDS(hu("FIG5", "allcelltypes.rds"))
allcelltypes_human <- subset(allcelltypes,species%in%'human')
allcelltypes_human_Y <- subset(allcelltypes_human,aged%in%'Young')
allcelltypes_human_O <- subset(allcelltypes_human,aged%in%'Aged')
table(allcelltypes_human_Y$subtype,allcelltypes_human_Y$species)
allcelltypes_monkey <- subset(allcelltypes,species%in%'monkey')
allcelltypes_monkey <- subset(allcelltypes_monkey,subtype%in%c('artery_ec','capillary_ec','vein_ec','lymphatic_ec',
                                                               'FB1','FB2','FB3','FB4','FB5','FB6',
                                                               'MC_Resident','MC_CCR2','MC_INF','MC_MHChi','MC_MHClow',
                                                               'DC','Mono','Cycle','B','T',
                                                               'Pericytes','SMC','Neural'))
allcelltypes_monkey_Y <- subset(allcelltypes_monkey,aged%in%'Young')
allcelltypes_monkey_O <- subset(allcelltypes_monkey,aged%in%'Aged')

allcelltypes_mouse <- subset(allcelltypes,dataset%in%'mouse3')
allcelltypes_mouse_Y <- subset(allcelltypes_mouse,aged%in%'Young')
allcelltypes_mouse_O <- subset(allcelltypes_mouse,aged%in%'Aged')

idents_order <- c('artery_ec','capillary_ec','vein_ec','lymphatic_ec',
                  'FB1','FB2','FB3','FB4','FB5','FB6',
                  'MC_Resident','MC_CCR2','MC_INF','MC_MHChi','MC_MHClow',
                  'DC','Mono','Cycle','B','T',
                  'Pericytes','SMC','Neural','Adipocyte')
idents_order_monkey <- c('artery_ec','capillary_ec','vein_ec','lymphatic_ec',
  'FB1','FB2','FB3','FB4','FB5','FB6',
  'MC_Resident','MC_CCR2','MC_INF','MC_MHChi','MC_MHClow',
  'DC','Mono','Cycle','B','T',
  'Pericytes','SMC','Neural')
future::plan("multicore", workers = 12) # do parallel  这里似乎有一些bug，在Linux上居然不行。de了它。
run_cellchat <- function(sce_input){
  data.input  <- sce_input@assays$RNA@data
  meta <- sce_input@meta.data
  cellchat <- createCellChat(object = data.input,meta = meta,group.by = 'subtype')
  #idents_order <- idents_order
  #cellchat@idents <- factor(cellchat@idents,levels = idents_order)
  cellchat@DB <- CellChatDB.human# set the used database in the object
  cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  cellchat <- updateCellChat(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)  
  cellchat <- computeCommunProb(cellchat,type = "truncatedMean",trim = 0.05, population.size = TRUE)  
  cellchat <- filterCommunication(cellchat, min.cells = 3)
  #Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- computeCommunProbPathway(cellchat) %>% 
    netAnalysis_computeCentrality(slot.name = "netP")
  return(cellchat)
}
dplyr::glimpse(cellchat)
allcelltypes_human_Y_obj <- run_cellchat(allcelltypes_human_Y)
allcelltypes_human_O_obj <- run_cellchat(allcelltypes_human_O)
write_rds(allcelltypes_human_Y_obj,file = hu("FIG5", "allcelltypes_human_Y_obj1.rds"))
write_rds(allcelltypes_human_O_obj,file = hu("FIG5", "allcelltypes_human_O_obj1.rds"))


write_rds(allcelltypes_monkey_Y,file = hu("FIG5", "allcelltypes_monkey_Y.rds"))
write_rds(allcelltypes_monkey_O,file = hu("FIG5", "allcelltypes_monkey_O.rds"))
allcelltypes_monkey_Y <- read_rds(hu("FIG5", "allcelltypes_monkey_Y1.rds"))
allcelltypes_monkey_O <- read_rds(hu("FIG5", "allcelltypes_monkey_O1.rds"))
allcelltypes_monkey_Y <- subset(allcelltypes_monkey_Y,subtype%in%idents_order)
allcelltypes_monkey_O <- subset(allcelltypes_monkey_O,subtype%in%idents_order)
allcelltypes_monkey_Y_obj <- run_cellchat(allcelltypes_monkey_Y)
allcelltypes_monkey_O_obj <- run_cellchat(allcelltypes_monkey_O)
write_rds(allcelltypes_monkey_Y_obj,file = hu("FIG5", "allcelltypes_monkey_Y_obj1.rds"))
write_rds(allcelltypes_monkey_O_obj,file = hu("FIG5", "allcelltypes_monkey_O_obj1.rds"))

write_rds(allcelltypes_mouse_Y,file = hu("FIG5", "allcelltypes_mouse_Y.rds"))
write_rds(allcelltypes_mouse_O,file = hu("FIG5", "allcelltypes_mouse_O.rds"))
allcelltypes_mouse_Y <- read_rds(hu("FIG5", "allcelltypes_mouse_Y.rds"))
allcelltypes_mouse_O <- read_rds(hu("FIG5", "allcelltypes_mouse_O.rds"))
idents_order <- unique(allcelltypes_mouse_Y$subtype)
allcelltypes_mouse_Y_obj <- run_cellchat(allcelltypes_mouse_Y)
allcelltypes_mouse_O_obj <- run_cellchat(allcelltypes_mouse_O)
write_rds(allcelltypes_mouse_Y_obj,file = hu("FIG5", "allcelltypes_mouse_Y_obj3_1.rds"))
write_rds(allcelltypes_mouse_O_obj,file = hu("FIG5", "allcelltypes_mouse_O_obj3_1.rds"))

