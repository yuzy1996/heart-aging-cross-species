# =============================================================================
# 10_fig1_TF_GO_enrichment.R
# Title : Transcription-factor and GO enrichment
# Figure: Fig. 1G; Supplementary
# Module: 02_figure1_cross_species_atlas
# Description:
#   msigdbr/ClusterGVis GO and transcription-factor (TF/EpiFactor database) enrichment across cell types.
# Inputs:
#   - hu("res", "combine_species_new_label.rds")
#   - hu("res", "term.rda")
#   - GO_DATA.RData
#   - cr("human", "crossspecies", "res", "metabolic_pathway_all_gsva.rds")
#   - refdb("human", "DatabaseExtract_v_1.01.csv")
#   - refdb("human", "epiFactor", "EpiGenes_main.csv")
# Outputs:
#   - hu("fig", "V1", "fig1_GO_species.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

# (R library search paths are managed by renv / .Renviron, not set per-script)
library(msigdbr)
library(Seurat)
library(ClusterGVis)
library(tidyverse)
# load data ---------------------------------------------------------------
combine_species_new_label <- read_rds(hu("res", "combine_species_new_label.rds"))
combine_species_new_label@meta.data$new_group_aged <- paste(combine_species_new_label$dataset,
                                                            combine_species_new_label$aged,
                                                            combine_species_new_label$annotation,sep = '_')
combine_species_new_label@meta.data$new_group_aged  <- as.factor(combine_species_new_label@meta.data$new_group_aged)
Idents(combine_species_new_label) <- 'new_group_aged'

combine_species_new_label_group_aged <- AverageExpression(combine_species_new_label,assays = 'RNA')$RNA
head(combine_species_new_label_group_aged)
dim(combine_species_new_label_group_aged)


combine_species_new_label_sub <- subset(combine_species_new_label,
                                        dataset%in%c('human','monkey','mouse1'))
combine_species_new_label_sub@meta.data$species_aged <- paste(combine_species_new_label_sub$dataset,
                                                              combine_species_new_label_sub$aged,sep = '_')
combine_species_new_label_sub@meta.data$species_aged  <- factor(combine_species_new_label_sub@meta.data$species_aged,
                                                                levels = c('human_Young','human_Aged',
                                                                           'monkey_Young','monkey_Aged',
                                                                           'mouse1_Young','mouse1_Aged'))
Idents(combine_species_new_label_sub) <- 'species_aged'
exp_species_aged <- AverageExpression(combine_species_new_label_sub,assays = 'RNA')$RNA
pheatmap::pheatmap(exp_species_aged[1:10,])
#save(terms_up,terms_down,file = hu("res", "term.rda"))
load(file = hu("res", "term.rda"))

combine_species_new_label@meta.data$species_aged <- paste(combine_species_new_label$dataset,
                                                              combine_species_new_label$aged,sep = '_')
Idents(combine_species_new_label) <- 'species_aged'
table(combine_species_new_label$species_aged)
combine_species_new_label <- ScaleData(combine_species_new_label,assay = 'RNA',features = row.names(combine_species_new_label))
exp_species_aged <- AverageExpression(combine_species_new_label,assays = 'RNA',slot = 'scale.data')$RNA
dim(exp_species_aged)
res_gsva_exp_species_aged_up <- GSVA::gsva(expr = exp_species_aged,gset.idx.list = terms_up)
res_gsva_exp_species_aged_down <- GSVA::gsva(expr = exp_species_aged,gset.idx.list = terms_down)
mat_go <- rbind(res_gsva_exp_species_aged_up[,1:6],
                res_gsva_exp_species_aged_down[,1:6])
term_cm8_up <- clusterData(exp = res_gsva_exp_species_aged_up[,1:6],
                        cluster.method = "mfuzz",
                        cluster.num = 8)
term_cm12_down <- clusterData(exp = res_gsva_exp_species_aged_down[,1:6],
                        cluster.method = "mfuzz",
                        cluster.num = 12)
pdf(hu("fig", "V1", "fig1_GO_species.pdf"),width = 10,height = 10)
visCluster(object = term_cm8_up,
           plot.type = "heatmap",
           cluster.order = c(1:8)) 
visCluster(object = term_cm12_down,
           plot.type = "heatmap",
           cluster.order = c(1:12)) 
dev.off()

pheatmap::pheatmap(res_gsva_exp_species_aged_up[filter(term_cm8_up$wide.res,cluster=='2')$gene,1:6],
                   cluster_cols = F,scale = 'row')
pheatmap::pheatmap(res_gsva_exp_species_aged_down[filter(term_cm8_down$wide.res,cluster=='6')$gene,1:6],cluster_cols = F,scale = 'row')

# GO term -----------------------------------------------------------------
GO_DATA <- clusterProfiler:::get_GO_data("org.Hs.eg.db", "BP", "SYMBOL") 　
PATHID2EXTID <- GO_DATA$PATHID2EXTID 
PATHID2NAME <- GO_DATA$PATHID2NAME

names(PATHID2EXTID)
names(PATHID2NAME)

ID_up <- names(PATHID2NAME[PATHID2NAME%in%terms_cross_up])
terms_up <- PATHID2EXTID[ID_up]
names(terms_up) <- as.character(PATHID2NAME[PATHID2NAME%in%terms_cross_up])

ID_down <- names(PATHID2NAME[PATHID2NAME%in%terms_cross_down])
terms_down <- PATHID2EXTID[ID_down]
names(terms_down) <- as.character(PATHID2NAME[PATHID2NAME%in%terms_cross_down])



getGO <- function(ID){
  
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  allNAME = names(GO_DATA$PATHID2EXTID)
  if(ID %in% allNAME){
    geneSet = GO_DATA$PATHID2EXTID[ID]
    names(geneSet) = GO_DATA$PATHID2NAME[ID]
    return(geneSet)
  } else{
    cat("No results!\n")
  }
}


PATH_ID_NAME <- merge(PATHID2EXTID, PATHID2NAME, by="from")

hsa_kegg <- clusterProfiler::download_KEGG("hsa")

names(hsa_kegg)

head(hsa_kegg$KEGGPATHID2NAME)

head(hsa_kegg$KEGGPATHID2EXTID)

PATH2ID <- hsa_kegg$KEGGPATHID2EXTID
PATH2NAME <- hsa_kegg$KEGGPATHID2NAME
PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by="from")
colnames(PATH_ID_NAME) <- c("KEGGID", "ENTREZID", "DESCRPTION")


cellCycleGO <- names(GO_DATA$PATHID2NAME[grep("cell cycle|DNA replication|cell division|segregation", GO_DATA$PATHID2NAME)])

cellCycleGene <- unique(unlist(GO_DATA$PATHID2EXTID[cellCycleGO]))



GO_data <- unlist(GO_DATA$PATHID2EXTID)

terms_C5 <-  msigdbr::msigdbr(species = "human", category = c("C5")) %>%
  dplyr::distinct(gs_name, gene_symbol) %>% 
  dplyr::mutate(gs_name=tolower(.$gs_name))

terms_C5_aging_related <- dplyr::filter(terms_C5,str_detect(terms_C5$gs_name,pattern =c("^go")))

terms_cross_up <- c(genes_human_up_enrich@result$Description[1:50],
           genes_monkey_up_enrich@result$Description[1:50],
           genes_mouse_up_enrich@result$Description[1:50]) %>% unique()

terms_cross_up <- c("nuclear division",
                    "meiosis I cell cycle process",
                    "regulation of signaling receptor activity",
                    "defense response to virus",
                    "interleukin-27-mediated signaling pathway",
                    "regulation of ribonuclease activity",
                    )
terms_cross_down <- c(genes_human_down_enrich@result$Description[1:50],
                    genes_monkey_down_enrich@result$Description[1:50],
                    genes_mouse_down_enrich@result$Description[1:50]) %>% unique()



terms_aging <- dplyr::filter(terms_C5,str_detect(terms_C5$gs_name,pattern =c("^go"))) %>% 
  dplyr::filter(.,str_detect(.$gs_name,pattern =c("aging")))
terms_inflam <- dplyr::filter(terms_C5,str_detect(terms_C5$gs_name,pattern =c("^go"))) %>% 
  dplyr::filter(.,str_detect(.$gs_name,pattern =c("inflammatory")))
terms_inflam <- dplyr::filter(terms_C5,str_detect(terms_C5$gs_name,pattern =c("^go"))) %>% 
  dplyr::filter(.,str_detect(.$gs_name,pattern =c("inflammatory")))
terms_damage <- dplyr::filter(terms_C5,str_detect(terms_C5$gs_name,pattern =c("^go"))) %>% 
  dplyr::filter(.,str_detect(.$gs_name,pattern =c("dna_damage")))

terms_damage <- dplyr::filter(terms_C5,str_detect(terms_C5$gs_name,pattern =c("^go"))) %>% 
  dplyr::filter(.,str_detect(.$gs_name,pattern =c("dna_damage")))

terms_aging <- split(terms_aging$gene_symbol,terms_aging$gs_name)
terms_inflam <- split(terms_inflam$gene_symbol,terms_inflam$gs_name)
terms_damage <- split(terms_damage$gene_symbol,terms_damage$gs_name)
metabolic_pathway_all_gsva <- read_rds(cr("human", "crossspecies", "res", "metabolic_pathway_all_gsva.rds"))
TF <- read_csv(refdb("human", "DatabaseExtract_v_1.01.csv"))
EpiFactor <- read_csv(refdb("human", "epiFactor", "EpiGenes_main.csv"))
res_gsva_exp <- GSVA::gsva(expr = combine_species_new_label_group_aged,gset.idx.list = metabolic_pathway_all_gsva)

pheatmap::pheatmap(res_gsva_exp)
tf_up <- intersect(TF$`HGNC symbol`,cross_up_gene)
tf_down <- intersect(TF$`HGNC symbol`,cross_down_gene)

Epi_up <- intersect(EpiFactor$HGNC_symbol,cross_up_gene)
Epi_down <- intersect(EpiFactor$HGNC_symbol,cross_down_gene)

human_celltypes <- paste('human','Aged',celltypes,sep = '-')


