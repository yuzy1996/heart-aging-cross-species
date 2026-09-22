# =============================================================================
# 19_fig1D_senMayo_celltype_enrichment.R
# Title : SenMayo senescence score by cell type
# Figure: Fig. 1D
# Module: 02_figure1_cross_species_atlas
# Description:
#   AddModuleScore of the SenMayo gene set and cross-cell-type up/down enrichment (legacy file name carries Fig2).
# Inputs:
#   - hu("res", "combine_species_new_label_S4_sub.rds")
# Outputs:
#   - hu("FIG2", "Mayo.pdf")
#   - hu("FIG2", "cross_celltype_up.pdf")
#   - hu("FIG2", "cross_celltype_down.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

Mayo <- clusterProfiler::read.gmt(lmx("GeneSet_aging", "SAUL_SEN_MAYO.v2024.1.Hs.gmt"))
celltypes_col <- c("#FF1493", "#FFB6C1", "#483D8B", "#00BFFF", "#228B22", "#7FFF00", "#FFFF00", "#FF4500","#8DD3C7")
species_col <- c("#FC8D62","#8DD3C7","#BC80BD")
aged_coll <- c("#483D8B","#00BFFF")
celltypes <- c('FB','EC','Myeloid','T','B','Pericytes','SMC','Neural','Adipocyte')
species <- c('human','monkey','mouse')
aged <- c('Aged','Young')
combine_species_new_label <- read_rds(hu("res", "combine_species_new_label_S4_sub.rds"))


combine_species_new_label <- AddModuleScore(combine_species_new_label,
                                            features = list(Mayo$gene),name = 'Mayo')
combine_species_new_label$age <- factor(combine_species_new_label$age,levels = c('Young','Aged'))
combine_species_new_label$species_age_group <- paste(combine_species_new_label$species,combine_species_new_label$aged,sep = '_')
combine_species_new_label$species_age_group <- factor(combine_species_new_label$species_age_group,
                                     levels = c('human_Young','human_Aged',
                                                'monkey_Young','monkey_Aged',
                                                'mouse_Young','mouse_Aged'))
combine_species_new_label$species <- factor(combine_species_new_label$species,levels = c('human','monkey','mouse'))
pdf(hu("FIG2", "Mayo.pdf"),width = 5,height = 6)
SCP::FeatureStatPlot(combine_species_new_label,
                     group.by = 'species_age_group',
                     split.by = 'species',
                     stat.by = 'Mayo1',comparisons = list(c('human_Young','human_Aged'),
                                                          c('monkey_Young','monkey_Aged'),
                                                          c('mouse_Young','mouse_Aged')),
                     add_box = T,palcolor = c("#FC8D62","#8DD3C7","#BC80BD"),
                     add_trend = T)
dev.off()


library(VennDiagram)
library(ComplexHeatmap)
library(colorRamp2)
library(clusterProfiler)
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
    marker_c1 <- marker_c0%>%
      mutate(diff = case_when(avg_log2FC>0.25~'up',
                              avg_log2FC<(-0.25)~'down'),
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
celltype_human <- subset(combine_species_new_label,species%in%'human')
celltype_monkey <- subset(combine_species_new_label,species%in%'monkey')
celltype_mouse <- subset(combine_species_new_label,species%in%'mouse')

celltype_human_diff_list <- run_cluster_diff(celltype_human,Clusters=celltypes)
celltype_monkey_diff_list <- run_cluster_diff(celltype_monkey,Clusters=celltypes)
celltype_mouse_diff_list <- run_cluster_diff(celltype_mouse,Clusters=celltypes)

celltype_human_markers <- Reduce(rbind,celltype_human_diff_list)
celltype_monkey_markers <- Reduce(rbind,celltype_monkey_diff_list)
celltype_mouse_markers <- Reduce(rbind,celltype_mouse_diff_list)

celltype_human_markers$species='human'
celltype_monkey_markers$species='monkey'
celltype_mouse_markers$species='mouse'
celltype_human_markers_up <- filter(celltype_human_markers,diff%in%'up')
celltype_monkey_markers_up <- filter(celltype_monkey_markers,diff%in%'up')
celltype_mouse_markers_up <- filter(celltype_mouse_markers,diff%in%'up')

celltype_human_markers_down <- filter(celltype_human_markers,diff%in%'down')
celltype_monkey_markers_down <- filter(celltype_monkey_markers,diff%in%'down')
celltype_mouse_markers_down <- filter(celltype_mouse_markers,diff%in%'down')

# 使用dplyr的group_by和summarize函数统计每个基因的出现频率
gene_frequency_human_up <- celltype_human_markers_up %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(frequency = n())%>%
  dplyr::filter(frequency > 5) %>% arrange(desc(frequency))
gene_frequency_monkey_up <- celltype_monkey_markers_up %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(frequency = n())%>%
  dplyr::filter(frequency > 3) %>% arrange(desc(frequency))
gene_frequency_mouse_up <- celltype_mouse_markers_up %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(frequency = n())%>%
  dplyr::filter(frequency > 4) %>% arrange(desc(frequency))
# 将筛选后的基因与原始数据框进行匹配，保留相应的行
result_df_human_up <- celltype_human_markers_up %>%
  dplyr::filter(gene %in% gene_frequency_human_up$gene)
result_df_human_up$gene <- factor(result_df_human_up$gene,levels = gene_frequency_human_up$gene)
result_df_human_up$num <- 1

result_df_monkey_up <- celltype_monkey_markers_up %>%
  dplyr::filter(gene %in% gene_frequency_monkey_up$gene[1:30])
result_df_monkey_up$gene <- factor(result_df_monkey_up$gene,levels = gene_frequency_monkey_up$gene[1:30])
result_df_monkey_up$num <- 1

result_df_mouse_up <- celltype_mouse_markers_up %>%
  dplyr::filter(gene %in% gene_frequency_mouse_up$gene)
result_df_mouse_up$gene <- factor(result_df_mouse_up$gene,levels = gene_frequency_mouse_up$gene)
result_df_mouse_up$num <- 1


celltype_human_markers_down <- filter(celltype_human_markers,diff%in%'down')
celltype_monkey_markers_down <- filter(celltype_monkey_markers,diff%in%'down')
celltype_mouse_markers_down <- filter(celltype_mouse_markers,diff%in%'down')
# 使用dplyr的group_by和summarize函数统计每个基因的出现频率
gene_frequency_human_down <- celltype_human_markers_down %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(frequency = n())%>%
  dplyr::filter(frequency > 5) %>% arrange(desc(frequency))
gene_frequency_monkey_down <- celltype_monkey_markers_down %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(frequency = n())%>%
  dplyr::filter(frequency > 4) %>% arrange(desc(frequency))
gene_frequency_mouse_down <- celltype_mouse_markers_down %>%
  dplyr::group_by(gene) %>%
  dplyr::summarize(frequency = n())%>%
  dplyr::filter(frequency > 4) %>% arrange(desc(frequency))
# 将筛选后的基因与原始数据框进行匹配，保留相应的行
result_df_human_down <- celltype_human_markers_down %>%
  dplyr::filter(gene %in% gene_frequency_human_down$gene)
result_df_human_down$gene <- factor(result_df_human_down$gene,levels = gene_frequency_human_down$gene)
result_df_human_down$num <- 1

result_df_monkey_down <- celltype_monkey_markers_down %>%
  dplyr::filter(gene %in% gene_frequency_monkey_down$gene)
result_df_monkey_down$gene <- factor(result_df_monkey_down$gene,levels = gene_frequency_monkey_down$gene)
result_df_monkey_down$num <- 1

result_df_mouse_down <- celltype_mouse_markers_down %>%
  dplyr::filter(gene %in% gene_frequency_mouse_down$gene)
result_df_mouse_down$gene <- factor(result_df_mouse_down$gene,levels = gene_frequency_mouse_down$gene)
result_df_mouse_down$num <- 1

celltypes <- c('FB','EC','Myeloid',
               'T','B','Pericytes',
               'SMC','Neural','Adipocyte')
plot_coord_polar <- function(obj,palette="YlGnBu"){
  gene_frequency <- obj %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(frequency = n())%>%
    dplyr::filter(frequency > 3) %>% arrange(desc(frequency))
  genes <- gene_frequency$gene[1:40]
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

pdf(hu("FIG2", "cross_celltype_up.pdf"),width = 7,height = 7)
plot_coord_polar(celltype_human_markers_up,palette = 'YlOrRd')
plot_coord_polar(celltype_monkey_markers_up,palette = 'YlOrRd')
plot_coord_polar(celltype_mouse_markers_up,palette = 'YlOrRd')
dev.off()

celltype_human_markers_down_1 <- celltype_human_markers_down[!row.names(celltype_human_markers_down)%in%str_subset(row.names(celltype_human_markers_down),"\\."),]

celltype_monkey_markers_down_1 <- celltype_monkey_markers_down[!row.names(celltype_monkey_markers_down)%in%str_subset(row.names(celltype_monkey_markers_down),"\\."),]

celltype_mouse_markers_down_1 <- celltype_mouse_markers_down[!row.names(celltype_mouse_markers_down)%in%str_subset(row.names(celltype_mouse_markers_down),"\\."),]

pdf(hu("FIG2", "cross_celltype_down.pdf"),width = 6,height = 6)
plot_coord_polar(celltype_human_markers_down_1)
plot_coord_polar(celltype_monkey_markers_down_1)
plot_coord_polar(celltype_mouse_markers_down_1)
dev.off()
