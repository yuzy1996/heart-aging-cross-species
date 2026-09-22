# =============================================================================
# 13_fig1_celltype_age_heatmaps.R
# Title : Cell-type x age conserved-marker heat maps / alluvial
# Figure: Fig. 1E; Supplementary
# Module: 02_figure1_cross_species_atlas
# Description:
#   Heat maps and alluvial/Sankey of conserved up/down markers across cell types and age.
# Inputs:
#   - hu("res", "combine_species_new_label.rds")
# Outputs:
#   - hu("downsample", "cross_celltype_up.pdf")
#   - hu("downsample", "cross_celltype_down.pdf")
#   - hu("downsample", "cross_up_gene.pdf")
#   - hu("fig", "V1", "fig1_cross.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

library(VennDiagram)
library(ComplexHeatmap)
library(colorRamp2)
library(clusterProfiler)
# plot sankey -------------------------------------------------------------
library(RColorBrewer)
library(ggalluvial)
library(tidyverse)
library(Seurat)
#library(clusterProfiler)
# Example data (you need to replace this with your actual data)
celltypes_col <- c("#FF1493", "#FFB6C1", "#483D8B", "#00BFFF", "#228B22", "#7FFF00", "#FFFF00", "#FF4500","#8DD3C7")
species_col <- c("#FC8D62","#8DD3C7","#BC80BD")
aged_coll <- c("#483D8B","#00BFFF")
celltypes <- c('FB','EC','Myeloid','T','B','Pericytes','SMC','Neural','Adipocyte')
species <- c('human','monkey','mouse')
aged <- c('Aged','Young')
combine_species_new_label <- read_rds(hu("res", "combine_species_new_label.rds"))

combine_species_new_label@meta.data$new_group_aged <- paste(combine_species_new_label$dataset,
                                                       combine_species_new_label$aged,
                                                       combine_species_new_label$annotation,sep = '_')
combine_species_new_label@meta.data$new_group_aged  <- as.factor(combine_species_new_label@meta.data$new_group_aged)
Idents(combine_species_new_label) <- 'new_group_aged'
combine_species_new_label <- ScaleData(combine_species_new_label,assay = 'RNA',features = row.names(combine_species_new_label))
combine_species_new_label_group_aged <- AverageExpression(combine_species_new_label,assays = 'RNA',slot = 'scale.data')
dim(combine_species_new_label_group_aged$RNA)
run_cluster_diff <- function(obj,Clusters){
  Idents(obj) <- 'aged'
  Seurat::DefaultAssay(obj) <- 'RNA'
  marker_list <- lapply(Clusters, function(x){
    marker_c0 <- FindMarkers(subset(obj,annotation%in%c(x)),ident.1 = 'Aged',ident.2 = 'Young')
    marker_c1 <- marker_c0%>%dplyr::filter(p_val<0.05&abs(avg_log2FC)>0.3) %>% 
      mutate(diff = case_when(avg_log2FC>0.3~'up',
                              avg_log2FC<(-0.3)~'down'),
             cluster=rep(x,length(.$avg_log2FC)),
             gene=row.names(.))
    
    return(marker_c1)
  })
  names(marker_list) <- Clusters
  return(marker_list)
}

celltypes <- c('FB','EC','Myeloid',
               'T','B','Pericytes',
               'SMC','Neural')
celltype_human <- subset(combine_species_new_label,dataset%in%'human')
celltype_monkey <- subset(combine_species_new_label,dataset%in%'monkey')
celltype_mouse1 <- subset(combine_species_new_label,dataset%in%'mouse1')

celltype_human_diff_list <- run_cluster_diff(celltype_human,Clusters=celltypes)
celltype_monkey_diff_list <- run_cluster_diff(celltype_monkey,Clusters=celltypes)
celltype_mouse1_diff_list <- run_cluster_diff(celltype_mouse1,Clusters=celltypes)

celltype_human_markers <- Reduce(rbind,celltype_human_diff_list)
celltype_monkey_markers <- Reduce(rbind,celltype_monkey_diff_list)
celltype_mouse_markers <- Reduce(rbind,celltype_mouse1_diff_list)

cross_up_gene <- intersect(filter(celltype_human_markers,avg_log2FC>1) %>% .$gene,
          filter(celltype_monkey_markers,avg_log2FC>1) %>% .$gene) %>% 
          intersect(filter(celltype_mouse1_markers,avg_log2FC>1) %>% .$gene)
cross_down_gene <- intersect(filter(celltype_human_markers,avg_log2FC<(-1)) %>% .$gene,
                           filter(celltype_monkey_markers,avg_log2FC<(-1)) %>% .$gene) %>% 
  intersect(filter(celltype_mouse1_markers,avg_log2FC<(-1)) %>% .$gene)
cross_up_gene_enrich <- enrichGO(gene = cross_up_gene,
                                  OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                                  ont = "BP",keyType = "SYMBOL",
                                  minGSSize=1,pvalueCutoff = 0.05,
                                  pAdjustMethod = "BH")
cross_down_gene_enrich <- enrichGO(gene = cross_down_gene,
                                    OrgDb=org.Hs.eg.db::org.Hs.eg.db,
                                    ont = "BP",keyType = "SYMBOL",
                                    minGSSize=1,pvalueCutoff = 0.05,
                                    pAdjustMethod = "BH")
dotplot(cross_up_gene_enrich)
dotplot(cross_down_gene_enrich)


plot_coord_polar <- function(obj,palette="YlGnBu"){
  gene_frequency <- obj %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(frequency = n())%>%
    dplyr::filter(frequency > 3) %>% arrange(desc(frequency))
  genes <- gene_frequency$gene[1:30]
  data <- obj %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(frequency = n())%>%
    dplyr::filter(frequency > 3) %>% arrange(desc(frequency))
  genes <- gene_frequency$gene[1:30]
  data <- obj %>%
    dplyr::filter(gene %in% genes)
  data$gene <- factor(data$gene,levels = genes)
  data$num <- 1
  data$cluster <- factor(data$cluster,levels = celltypes)
  label=table(data$gene) %>% data.frame()
  colnames(label) <- c('gene','freq')
  label$number <- 1:nrow(label)
  angle <- 90 - 360 * (label$number-0.5) /nrow(label)
  label$hjust <- ifelse(angle < -90, 1, 0)
  label$angle <- ifelse(angle < -90, angle+180, angle)
  data <- left_join(label,data[,c('cluster','num','gene')],by='gene')
  ggplot(data,aes(x = gene,y = num,fill=cluster))+
    geom_bar(position = "stack",stat = 'identity',alpha=0.8)+
    coord_polar()+
    geom_text(data = data,aes(x = gene,y=freq+0.5,label=gene,hjust=hjust,angle=angle),size=4)+
    scale_fill_brewer(palette = palette)+ylim(-3,12)+
    theme(panel.background = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank())
}


celltype_human_markers_up <- filter(celltype_human_markers,diff%in%'up')
celltype_monkey_markers_up <- filter(celltype_monkey_markers,diff%in%'up')
celltype_mouse_markers_up <- filter(celltype_mouse_markers,diff%in%'up')


celltype_human_markers_down <- filter(celltype_human_markers,diff%in%'down')
celltype_monkey_markers_down <- filter(celltype_monkey_markers,diff%in%'down')
celltype_mouse_markers_down <- filter(celltype_mouse_markers,diff%in%'down')

plot_coord_polar(celltype_human_markers_up,palette = 'YlOrRd')
plot_coord_polar(celltype_monkey_markers_up,palette = 'YlOrRd')
plot_coord_polar(celltype_mouse_markers_up,palette = 'YlOrRd')

pdf(hu("downsample", "cross_celltype_up.pdf"),width = 7,height = 7)
plot_coord_polar(celltype_human_markers_up,palette = 'YlOrRd')
plot_coord_polar(celltype_monkey_markers_up,palette = 'YlOrRd')
plot_coord_polar(celltype_mouse_markers_up,palette = 'YlOrRd')
dev.off()

pdf(hu("downsample", "cross_celltype_down.pdf"),width = 7,height = 7)
plot_coord_polar(celltype_human_markers_down)
plot_coord_polar(celltype_monkey_markers_down)
plot_coord_polar(celltype_mouse_markers_down)
dev.off()


pdf(hu("downsample", "cross_up_gene.pdf"),width = 8,height = 8)
for (gene in cross_down_gene) {
  FeaturePlot(combine_species_new_label_human,
              split.by = 'aged',ncol = 3,
              features = c(gene),order = T)
  FeaturePlot(combine_species_new_label_monkey,
              split.by = 'aged',ncol = 3,
              features = c(gene),)
  FeaturePlot(combine_species_new_label_mouse,
              split.by = 'aged',ncol = 3,
              features = c(gene))
  
}
dev.off()



lapply(celltypes, function(cl){
  veen <- list('human'=filter(celltype_human_markers,cluster%in%cl&avg_log2FC>1) %>% .$gene,
               'monkey'=filter(celltype_monkey_markers,cluster%in%cl&avg_log2FC>1) %>% .$gene,
               'mouse'=filter(celltype_mouse1_markers,cluster%in%cl&avg_log2FC>1) %>% .$gene)
  return(veen)
})->celltypes_venn_up
names(celltypes_venn_up) <- celltypes
lapply(celltypes, function(cl){
  veen <- list('human'=filter(celltype_human_markers,cluster%in%cl&avg_log2FC<(-1)) %>% .$gene,
               'monkey'=filter(celltype_monkey_markers,cluster%in%cl&avg_log2FC<(-1)) %>% .$gene,
               'mouse'=filter(celltype_mouse1_markers,cluster%in%cl&avg_log2FC<(-1)) %>% .$gene)
  return(veen)
})->celltypes_venn_down
names(celltypes_venn_down) <- celltypes

FB_Aged_veen_up <- get.venn.partitions(celltypes_venn_up$FB)
EC_Aged_veen_up <- get.venn.partitions(celltypes_venn_up$EC)
Myeloid_Aged_veen_up <- get.venn.partitions(celltypes_venn_up$Myeloid)
T_Aged_veen_up <- get.venn.partitions(celltypes_venn_up$T)
B_Aged_veen_up <- get.venn.partitions(celltypes_venn_up$B)
Pericytes_Aged_veen_up <- get.venn.partitions(celltypes_venn_up$Pericytes)
SMC_Aged_veen_up <- get.venn.partitions(celltypes_venn_up$SMC)
Neural_Aged_veen_up <- get.venn.partitions(celltypes_venn_up$Neural)

FB_Aged_veen_down <- get.venn.partitions(celltypes_venn_down$FB)
EC_Aged_veen_down <- get.venn.partitions(celltypes_venn_down$EC)
Myeloid_Aged_veen_down <- get.venn.partitions(celltypes_venn_down$Myeloid)
T_Aged_veen_down <- get.venn.partitions(celltypes_venn_down$T)
B_Aged_veen_down <- get.venn.partitions(celltypes_venn_down$B)
Pericytes_Aged_veen_down <- get.venn.partitions(celltypes_venn_down$Pericytes)
SMC_Aged_veen_down <- get.venn.partitions(celltypes_venn_down$SMC)
Neural_Aged_veen_down <- get.venn.partitions(celltypes_venn_down$Neural)



genes_conserved_up <- c(FB_Aged_veen_up[1,5][[1]],
                        EC_Aged_veen_up[1,5][[1]],
                        Myeloid_Aged_veen_up[1,5][[1]],
                        T_Aged_veen_up[1,5][[1]],
                        B_Aged_veen_up[1,5][[1]],
                        Pericytes_Aged_veen_up[1,5][[1]],
                        SMC_Aged_veen_up[1,5][[1]],
                        Neural_Aged_veen_up[1,5][[1]])%>% unique()
genes_conserved_down <- c(FB_Aged_veen_down[1,5][[1]],
                          EC_Aged_veen_down[1,5][[1]],
                          Myeloid_Aged_veen_down[1,5][[1]],
                          T_Aged_veen_down[1,5][[1]],
                          B_Aged_veen_down[1,5][[1]],
                          Pericytes_Aged_veen_down[1,5][[1]],
                          SMC_Aged_veen_down[1,5][[1]],
                          Neural_Aged_veen_down[1,5][[1]])%>% unique()
genes_conserved_human <- c(FB_Aged_veen_up[7,5][[1]],
                        EC_Aged_veen_up[7,5][[1]],
                        Myeloid_Aged_veen_up[7,5][[1]],
                        T_Aged_veen_up[7,5][[1]],
                        B_Aged_veen_up[7,5][[1]],
                        Pericytes_Aged_veen_up[7,5][[1]],
                        SMC_Aged_veen_up[7,5][[1]],
                        Neural_Aged_veen_up[7,5][[1]],
                        FB_Aged_veen_down[7,5][[1]],
                        EC_Aged_veen_down[7,5][[1]],
                        Myeloid_Aged_veen_down[7,5][[1]],
                        T_Aged_veen_down[7,5][[1]],
                        B_Aged_veen_down[7,5][[1]],
                        Pericytes_Aged_veen_down[7,5][[1]],
                        SMC_Aged_veen_down[7,5][[1]],
                        Neural_Aged_veen_down[7,5][[1]])%>% unique()
genes_conserved_monkey <- c(FB_Aged_veen_up[6,5][[1]],
                        EC_Aged_veen_up[6,5][[1]],
                        Myeloid_Aged_veen_up[6,5][[1]],
                        T_Aged_veen_up[6,5][[1]],
                        B_Aged_veen_up[6,5][[1]],
                        Pericytes_Aged_veen_up[6,5][[1]],
                        SMC_Aged_veen_up[6,5][[1]],
                        Neural_Aged_veen_up[6,5][[1]],
                        FB_Aged_veen_down[6,5][[1]],
                        EC_Aged_veen_down[6,5][[1]],
                        Myeloid_Aged_veen_down[6,5][[1]],
                        T_Aged_veen_down[6,5][[1]],
                        B_Aged_veen_down[6,5][[1]],
                        Pericytes_Aged_veen_down[6,5][[1]],
                        SMC_Aged_veen_down[6,5][[1]],
                        Neural_Aged_veen_down[6,5][[1]])%>% unique()
genes_conserved_mouse_up <- c(FB_Aged_veen_up[4,5][[1]],
                        EC_Aged_veen_up[4,5][[1]],
                        Myeloid_Aged_veen_up[4,5][[1]],
                        T_Aged_veen_up[4,5][[1]],
                        B_Aged_veen_up[4,5][[1]],
                        Pericytes_Aged_veen_up[4,5][[1]],
                        SMC_Aged_veen_up[4,5][[1]],
                        Neural_Aged_veen_up[4,5][[1]])%>% unique()
                        
genes_conserved_mouse_down <- c(FB_Aged_veen_down[4,5][[1]],
                        EC_Aged_veen_down[4,5][[1]],
                        Myeloid_Aged_veen_down[4,5][[1]],
                        T_Aged_veen_down[4,5][[1]],
                        B_Aged_veen_down[4,5][[1]],
                        Pericytes_Aged_veen_down[4,5][[1]],
                        SMC_Aged_veen_down[4,5][[1]],
                        Neural_Aged_veen_down[4,5][[1]])%>% unique()

mat <- combine_species_new_label_group_aged$RNA %>% as.matrix() %>% t()
mat_diff <- mat[c(glue('human-Young-{celltypes}'),
                  glue('human-Aged-{celltypes}'),
                  glue('monkey-Young-{celltypes}'),
                  glue('monkey-Aged-{celltypes}'),
                  glue('mouse1-Young-{celltypes}'),
                  glue('mouse1-Aged-{celltypes}')),]

mat_corss <- mat_diff[,unique(c(genes_conserved_up,genes_conserved_down,
                           genes_conserved_human,
                           genes_conserved_monkey,
                           genes_conserved_mouse))]

mat_cross1 <- cbind(mat_diff[,genes_conserved_up],
                    mat_diff[,genes_conserved_down],
                    mat_diff[,genes_conserved_human],
                    mat_diff[,genes_conserved_monkey],
                    mat_diff[,genes_conserved_mouse])
dim(mat_cross1)
group = rep(c("Conserved", "Species_human","Species_monkey","Species_mouse"),
            times = c(length(c(genes_conserved_up,genes_conserved_down)),
                      length(genes_conserved_human),
                      length(genes_conserved_monkey),
                      length(genes_conserved_mouse)))
length(group)
col1 = colorRamp2(c(0, 5), c( "white", "#FC8D62"))
col2 = colorRamp2(c(0, 5), c("white", "#8DD3C7"))
col3 = colorRamp2(c(0, 5), c("white", "#BC80BD"))
col4 = colorRamp2(c(0, 5), c("white", "pink"))
clusters <- rep(celltypes, 6)
aged <- rep(c(rep('Young',8),rep('Aged',8)),3)
names(aged_col) <- aged
clusters_col <- rep(c("#FF1493", "#FFB6C1", "#483D8B", "#00BFFF", "#228B22", "#7FFF00", "#FFFF00", "#FF4500"), 6)
names(clusters_col) <- clusters
clusters <- factor(clusters,levels = celltypes)
aged <- factor(aged,levels = c('Young','Aged'))
left_annotation= rowAnnotation(clusters = clusters,
                               aged=aged, 
                               col = list(clusters = clusters_col,
                                          aged=aged_col))

Heatmap(mat_corss_scale[, 1:6], col = col1, name = "Conserved",
        left_annotation = left_annotation,
        row_split = c(rep('human',16),rep('monkey',16),rep('mouse',16)),
        cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = T)
mat_corss_scale <- scale(mat_cross1)
dim(mat_corss_scale)
ht1 = Heatmap(mat_cross1[, group == "Conserved"], col = col1, name = "Conserved",
              left_annotation = left_annotation,
              row_split = c(rep('human',16),rep('monkey',16),rep('mouse',16)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht1
ht2 = Heatmap(mat_cross1[, group == "Species_human"], col = col2, name = "Species_human",
              row_split = c(rep('human',16),rep('monkey',16),rep('mouse',16)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht2
ht3 = Heatmap(mat_cross1[, group == "Species_monkey"], col = col3, name = "Species_monkey",
              row_split = c(rep('human',16),rep('monkey',16),rep('mouse',16)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht4 = Heatmap(mat_cross1[, group == "Species_mouse"], col = col4, name = "Species_mouse",
              row_split = c(rep('human',16),rep('monkey',16),rep('mouse',16)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht4
ht1
pdf(hu("fig", "V1", "fig1_cross.pdf"),width = 12,height = 10)
draw(ht1 + ht2 + ht3 + ht4, padding = unit(c(20, 20, 20, 20), "mm"))
dev.off()

Heatmap(mat_corss_scale[, group == "Conserved"], col = col1, name = "Conserved",
        left_annotation = left_annotation,
        row_split = c(rep('human',16),rep('monkey',16),rep('mouse',16)),
        cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = T)
