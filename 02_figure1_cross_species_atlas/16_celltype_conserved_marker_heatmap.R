# =============================================================================
# 16_celltype_conserved_marker_heatmap.R
# Title : Conserved marker heat map across species
# Figure: Supplementary
# Module: 02_figure1_cross_species_atlas
# Description:
#   Heat map of shared cell-type marker expression across human/macaque/mouse.
# Inputs:
#   - hu("res", "combine_species_new_label.rds")
# Outputs:
#   - hu("FIG1", "celltype_heatmap_conserved.pdf")
#   - hu("FIG1", "marker_species.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

library(VennDiagram)
library(RColorBrewer)
library(circlize)
library(Seurat)
library(tidyverse)
combine_species_new_label <- read_rds(hu("res", "combine_species_new_label.rds"))
Idents(combine_species_new_label) <- 'annotation'
DefaultAssay(combine_species_new_label) <- 'RNA'
celltype_human <- subset(combine_species_new_label,species%in%'human')
celltype_monkey <- subset(combine_species_new_label,species%in%'monkey')
celltype_mouse <- subset(combine_species_new_label,species%in%'mouse')
celltype_human_markers <- FindAllMarkers(celltype_human,min.pct = 0.2)
celltype_monkey_markers <- FindAllMarkers(celltype_monkey,min.pct = 0.2)
celltype_mouse_markers <- FindAllMarkers(celltype_mouse,min.pct = 0.2)

celltypes <- c('FB','EC','Myeloid',
              'T','B','Pericytes',
              'SMC','Neural','Adipocyte')
lapply(celltypes, function(cl){
  veen <- list('human'=filter(celltype_human_markers,cluster%in%cl&avg_log2FC>1) %>% .$gene,
               'monkey'=filter(celltype_monkey_markers,cluster%in%cl&avg_log2FC>1) %>% .$gene,
               'mouse'=filter(celltype_mouse_markers,cluster%in%cl&avg_log2FC>1) %>% .$gene)
  return(veen)
})->celltypes_venn_FC1
names(celltypes_venn_FC1) <- celltypes
lapply(celltypes, function(cl){
  veen <- list('human'=filter(celltype_human_markers,cluster%in%cl&avg_log2FC>2) %>% .$gene,
               'monkey'=filter(celltype_monkey_markers,cluster%in%cl&avg_log2FC>2) %>% .$gene,
               'mouse'=filter(celltype_mouse_markers,cluster%in%cl&avg_log2FC>2) %>% .$gene)
  return(veen)
})->celltypes_venn_FC2
names(celltypes_venn_FC2) <- celltypes
lapply(celltypes, function(cl){
  veen <- list('human'=filter(celltype_human_markers,cluster%in%cl&avg_log2FC>3) %>% .$gene,
               'monkey'=filter(celltype_monkey_markers,cluster%in%cl&avg_log2FC>3) %>% .$gene,
               'mouse'=filter(celltype_mouse_markers,cluster%in%cl&avg_log2FC>3) %>% .$gene)
  return(veen)
})->celltypes_venn_FC3
names(celltypes_venn_FC3) <- celltypes
conserved_veen <- function(celltypes_venn){
  FB_veen <- get.venn.partitions(celltypes_venn$FB)
  EC_veen <- get.venn.partitions(celltypes_venn$EC)
  Myeloid_veen <- get.venn.partitions(celltypes_venn$Myeloid)
  T_veen <- get.venn.partitions(celltypes_venn$T)
  B_veen <- get.venn.partitions(celltypes_venn$B)
  Pericytes_veen <- get.venn.partitions(celltypes_venn$Pericytes)
  SMC_veen <- get.venn.partitions(celltypes_venn$SMC)
  Neural_veen <- get.venn.partitions(celltypes_venn$Neural)
  Adipocyte_veen <- get.venn.partitions(celltypes_venn$Adipocyte)
  celltypes_genes_conserved <- c(FB_veen[1,5][[1]],
                                 EC_veen[1,5][[1]],
                                 Myeloid_veen[1,5][[1]],
                                 T_veen[1,5][[1]],
                                 B_veen[1,5][[1]],
                                 Pericytes_veen[1,5][[1]],
                                 SMC_veen[1,5][[1]],
                                 Neural_veen[1,5][[1]]) %>% unique()
  celltypes_genes_species_huamn <- c(FB_veen[7,5][[1]],
                                     EC_veen[7,5][[1]],
                                     Myeloid_veen[7,5][[1]],
                                     T_veen[7,5][[1]],
                                     B_veen[7,5][[1]],
                                     Pericytes_veen[7,5][[1]],
                                     SMC_veen[7,5][[1]],
                                     Neural_veen[7,5][[1]]) %>% unique()
  celltypes_genes_species_monkey <- c(FB_veen[6,5][[1]],
                                      EC_veen[6,5][[1]],
                                      Myeloid_veen[6,5][[1]],
                                      T_veen[6,5][[1]],
                                      B_veen[6,5][[1]],
                                      Pericytes_veen[6,5][[1]],
                                      SMC_veen[6,5][[1]],
                                      Neural_veen[6,5][[1]]) %>% unique()
  celltypes_genes_species_mouse <- c(FB_veen[4,5][[1]],
                                     EC_veen[4,5][[1]],
                                     Myeloid_veen[4,5][[1]],
                                     T_veen[4,5][[1]],
                                     B_veen[4,5][[1]],
                                     Pericytes_veen[4,5][[1]],
                                     SMC_veen[4,5][[1]],
                                     Neural_veen[4,5][[1]]) %>% unique()
  return(list('celltypes_genes_conserved'=celltypes_genes_conserved,
              'celltypes_genes_species_huamn'=celltypes_genes_species_huamn,
              'celltypes_genes_species_monkey'=celltypes_genes_species_monkey,
              'celltypes_genes_species_mouse'=celltypes_genes_species_mouse))
}
conserved_veen_FC3 <- conserved_veen(celltypes_venn_FC3)
conserved_veen_FC2 <- conserved_veen(celltypes_venn_FC2)
conserved_veen_FC1 <- conserved_veen(celltypes_venn_FC1)
celltypes_genes_conserved=conserved_veen_FC1$celltypes_genes_conserved
celltypes_genes_species_huamn=conserved_veen_FC2$celltypes_genes_species_huamn
celltypes_genes_species_monkey=conserved_veen_FC2$celltypes_genes_species_monkey
celltypes_genes_species_mouse=conserved_veen_FC3$celltypes_genes_species_mouse
combine_species_new_label@meta.data$new_group <- paste(combine_species_new_label$species,
                                                       combine_species_new_label$annotation,sep = '_')
combine_species_new_label@meta.data$new_group <- factor(combine_species_new_label@meta.data$new_group)
Idents(combine_species_new_label) <- 'new_group'
mat_corss <- AverageExpression(combine_species_new_label)
mat_corss$RNA[1:10,1:10]
# label ann ---------------------------------------------------------------
group = rep(c("Conserved", "Species_human","Species_monkey","Species_mouse"),
            times = c(length(c(celltypes_genes_conserved)),
                      length(celltypes_genes_species_huamn),
                      length(celltypes_genes_species_monkey),
                      length(celltypes_genes_species_mouse)))
col1 = colorRamp2(c(0, 5), c( "white", "#FC8D62"))
col2 = colorRamp2(c(0, 5), c("white", "#8DD3C7"))
col3 = colorRamp2(c(0, 5), c("white", "#BC80BD"))
col4 = colorRamp2(c(0, 5), c("white", "pink"))
clusters <- rep(c('FB','EC','Myeloid',
                  'T','B','Pericytes',
                  'SMC','Neural'), 3)
cell_colors <- c("#FF1493", "#FFB6C1", "#483D8B", "#00BFFF", "#228B22", "#7FFF00", "#FFFF00", "#FF4500","#8DD3C7")
clusters_col <- rep(cell_colors[1:8], 3)
names(clusters_col) <- clusters
library(ComplexHeatmap)
clusters <- factor(clusters,levels =celltypes1)
left_annotation= rowAnnotation(clusters = clusters,
                               col = list(clusters = clusters_col))
ha <- HeatmapAnnotation(treatment = treatments,
                        col = list(treatment = c("control" = "blue", "treated" = "red")),
                        annotation_name_gp = gpar(fontsize = 12), 
                        annotation_legend_param = list(title_gp = gpar(fontsize = 12)))

mat_corss_scale <- mat_corss$RNA %>% as.matrix() %>% t()
mat_corss_scale <- mat_corss_scale[,c(celltypes_genes_conserved,
                                      celltypes_genes_species_huamn,
                                      celltypes_genes_species_monkey,
                                      celltypes_genes_species_mouse)]
celltypes1 <- c('FB','EC','Myeloid','T','B','Pericytes','SMC','Neural')
mat_corss_scale <- mat_corss_scale[c(paste('human',celltypes1,sep = '-'),
                                     paste('monkey',celltypes1,sep = '-'),
                                     paste('mouse',celltypes1,sep = '-')),]
#write_rds(mat_corss_scale,hu("res", "mat_corss_scale.rds"))
#mat_corss_scale <- read_rds(hu("res", "mat_corss_scale.rds"))
ht1 = Heatmap(mat_corss_scale[, group == "Conserved"], 
              col = col1, name = "Conserved",
              left_annotation = left_annotation,
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht2 = Heatmap(mat_corss_scale[, group == "Species_human"], col = col2, name = "Species_human",
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht3 = Heatmap(mat_corss_scale[, group == "Species_monkey"], col = col3, name = "Species_monkey",
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht4 = Heatmap(mat_corss_scale[, group == "Species_mouse"], col = col4, name = "Species_mouse",
              row_split = c(rep('human',8),rep('monkey',8),rep('mouse',8)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
pdf(hu("FIG1", "celltype_heatmap_conserved.pdf"),width = 10,height = 4)
combined_heatmap <- ht1 + ht2 + ht3 + ht4
# 绘制合并后的热图，设置热图之间的间距
draw(combined_heatmap, padding = unit(c(20, 20, 20, 20), "mm"))
dev.off()



combine_species_new_label@meta.data$new_group1 <- paste(combine_species_new_label$annotation,
                                                       combine_species_new_label$species,sep = '_')
combine_species_new_label@meta.data$new_group1 <- factor(combine_species_new_label@meta.data$new_group1)

pdf(hu("FIG1", "marker_species.pdf"),width = 12,height = 12)
scop::FeatureStatPlot(
  combine_species_new_label,
  stat.by = c(
    "DCN", "COL1A1", "FBLN1", # Ductal
    'TCF21',
    'CDH5',
    'PECAM1',
    'RAMP2',
    'EMCN',
    'TEK',
    'CD163',
    'CD68',
    'MS4A6A',
    'CSF1R',
    'CD3D',
    'CD3E',
    'CD8A',
    'IL7R',
    'IGKC',
    'CD19',
    'CD79A',
    'RGS5',
    'ABCC9',
    'KCNJ8',
    'MYH11',
    'TAGLN',
    'ACTA2',
    'PLP1',
    'NRXN1',
    'NRXN3',
    'ADIPOQ',
    'PLIN1',
    'GPD1',
    'PNPLA3',
    'PCK1'
  ),
  legend.position = "top",
  legend.direction = "horizontal",
  group.by = "new_group1",
  bg.by = "annotation",
  stack = TRUE
)
dev.off()

scop::FeatureStatPlot(
  subset(combine_species_new_label,annotation%in%'FB'),
  stat.by = c("DCN", "COL1A1", "FBLN1", 'TCF21'),
  group.by = "new_group",
  plot_type = "bar"
)

species_col <- c("#FC8D62","#8DD3C7","#BC80BD")
pdf(hu("FIG1", "marker_species.pdf"),width = 5,height = 5)
scop::FeatureStatPlot(
  subset(combine_species_new_label,annotation%in%'EC'),
  stat.by = c('CDH5',
              'PECAM1',
              'TEK'),
  group.by = "species", palcolor = species_col,
  plot_type = "bar",stack = TRUE
)
scop::FeatureStatPlot(
  subset(combine_species_new_label,annotation%in%'FB'),
  stat.by = c("DCN", "COL1A1", "FBLN1"),
  group.by = "species", palcolor = species_col,
  plot_type = "bar",stack = TRUE
)
scop::FeatureStatPlot(
  subset(combine_species_new_label,annotation%in%'Myeloid'),
  stat.by = c( 'CD163',
               'MS4A6A',
               'CSF1R'),
  group.by = "species", palcolor = species_col,
  plot_type = "bar",stack = TRUE
)
dev.off()


scop::FeatureStatPlot(
  subset(combine_species_new_label,annotation%in%'Myeloid'),
  stat.by = c( 'CD163',
               'MS4A6A',
               'CSF1R'),
  group.by = "species", palcolor = species_col,
  plot_type = "bar",stack = TRUE
)