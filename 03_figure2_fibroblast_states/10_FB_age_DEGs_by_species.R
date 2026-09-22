# =============================================================================
# 10_FB_age_DEGs_by_species.R
# Title : Fibroblast age DEGs by species
# Figure: Fig. 2E-G
# Module: 03_figure2_fibroblast_states
# Description:
#   run_cluster_diff wrapper producing per-species fibroblast age-DEG lists.
# Inputs:
#   - (see script)
# Outputs:
#   - cr("human", "crossspecies", "res", "FB_human_diff_list.pdf")
#   - cr("human", "crossspecies", "res", "FB_monkey_diff_list.pdf")
#   - cr("human", "crossspecies", "res", "FB_mouse1_diff_list.pdf")
#   - cr("human", "crossspecies", "res", "FB_mouse2_diff_list.pdf")
#   - cr("human", "crossspecies", "res", "FB_mouse3_diff_list.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

run_cluster_diff <- function(obj,top_n = 30){
  Idents(obj) <- 'aged'
  Clusters <- data.frame(table(obj$seurat_clusters)) %>% filter(Freq>10) %>% .$Var1
  Seurat::DefaultAssay(obj) <- 'RNA'
  marker_list <- lapply(Clusters, function(x){
    marker_c0 <- FindMarkers(subset(obj,seurat_clusters%in%c(x)),ident.1 = 'Aged',ident.2 = 'Young')
    marker_c1 <- marker_c0%>%dplyr::filter(p_val<0.05&abs(avg_log2FC)>0.5) %>% 
      mutate(diff = case_when(avg_log2FC>0.5~'up',
                              avg_log2FC<(-0.5)~'down'),
             cluster=rep(x,length(.$avg_log2FC)),
             gene=row.names(.))
    
    return(marker_c1)
  })
  return(marker_list)
}
library(scRNAtoolVis)

subFB_human <- subset(FB_cross_new,species%in%'human')
subFB_monkey <- subset(FB_cross_new,species%in%'monkey')
subFB_mouse1 <- subset(FB_cross_new,species%in%'mouse1')
subFB_mouse2 <- subset(FB_cross_new,species%in%'mouse2')
subFB_mouse3 <- subset(FB_cross_new,species%in%'mouse3')

# human cell types diff gene ----------------------------------------------
FB_human_diff_list <- run_cluster_diff(subFB_human)
FB_human_diff_list <- Reduce(rbind,FB_human_diff_list) 
pdf(cr("human", "crossspecies", "res", "FB_human_diff_list.pdf"),width = 12,height = 12)
jjVolcano(diffData = FB_human_diff_list,
          log2FC.cutoff = 2,
          topGeneN =5,
          tile.col = corrplot::COL2('RdBu', 15)[1:15],
          size  = 3,
          fontface = 'italic',
          polar = T)+ylim(-8,10)+labs(title="human")
dev.off()
# monkey ------------------------------------------------------------------
FB_monkey_diff_list <- run_cluster_diff(subFB_monkey)
FB_monkey_diff_list <- Reduce(rbind,FB_monkey_diff_list) 
pdf(cr("human", "crossspecies", "res", "FB_monkey_diff_list.pdf"),width = 12,height = 12)
jjVolcano(diffData = FB_monkey_diff_list %>% filter(.,abs(avg_log2FC)>2),
          log2FC.cutoff = 2,
          topGeneN =5,
          tile.col = corrplot::COL2('RdBu', 15)[1:15],
          size  = 3,
          fontface = 'italic',
          polar = T)+ylim(-8,10)
dev.off()
# mouse1 ------------------------------------------------------------------
FB_mouse1_diff_list <- run_cluster_diff(subFB_mouse1)
FB_mouse1_diff_list <- Reduce(rbind,FB_mouse1_diff_list) 
pdf(cr("human", "crossspecies", "res", "FB_mouse1_diff_list.pdf"),width = 12,height = 12)
jjVolcano(diffData = FB_mouse1_diff_list,log2FC.cutoff = 2,
          topGeneN =5,
          tile.col = corrplot::COL2('RdBu', 15)[1:15],
          size  = 3,
          fontface = 'italic',
          polar = T)+ylim(-8,10)
dev.off()
# mouse2 ------------------------------------------------------------------
FB_mouse2_diff_list <- run_cluster_diff(subFB_mouse2)
FB_mouse2_diff_list <- Reduce(rbind,FB_mouse2_diff_list) 
pdf(cr("human", "crossspecies", "res", "FB_mouse2_diff_list.pdf"),width = 12,height = 12)
jjVolcano(diffData = FB_mouse2_diff_list,log2FC.cutoff = 2,
          topGeneN =5,
          tile.col = corrplot::COL2('RdBu', 15)[1:15],
          size  = 3,
          fontface = 'italic',
          polar = T)+ylim(-8,10)
dev.off()
# mouse3 ------------------------------------------------------------------
FB_mouse3_diff_list <- run_cluster_diff(subFB_mouse3)
FB_mouse3_diff_list <- Reduce(rbind,FB_mouse3_diff_list)
pdf(cr("human", "crossspecies", "res", "FB_mouse3_diff_list.pdf"),width = 12,height = 12)
jjVolcano(diffData = FB_mouse3_diff_list,log2FC.cutoff = 2,
          topGeneN =5,
          tile.col = corrplot::COL2('RdBu', 15)[1:15],
          size  = 3,
          fontface = 'italic',
          polar = T)+ylim(-8,10)
dev.off()


# 安装并加载必要的包
# (package installation is handled by environment/00_setup_env.R)  if (!require(dplyr)) install.packages("dplyr")
# (package installation is handled by environment/00_setup_env.R)  if (!require(ggplot2)) install.packages("ggplot2")
# (package installation is handled by environment/00_setup_env.R)  if (!require(scales)) install.packages("scales")

library(dplyr)
library(ggplot2)
library(scales)

# 示例数据框
set.seed(123)
metadata <- data.frame(
  cell_id = paste0("cell", 1:2400),
  species = rep(c("Human", "Monkey", "Mouse"), each = 800),
  cluster = sample(paste0("Cluster", 1:24), 2400, replace = TRUE),
  group = sample(c("Control", "Disease"), 2400, replace = TRUE),
  sample_id = rep(paste0("Sample", 1:12), each = 200)
)

metadata <- FB_cross_new@meta.data
# 计算每个样本、每个物种、每个组中每个cluster的细胞数量
cluster_counts <- metadata %>%
  group_by(species, aged, seurat_clusters) %>%
  summarize(count = n(), .groups = 'drop')

# 计算每个样本的总细胞数量
total_counts <- metadata %>%
  group_by(species, aged) %>%
  summarize(total = n(), .groups = 'drop')

# 合并数据框以计算比例
proportion_data <- merge(cluster_counts, total_counts, by = c("species", "aged"))
proportion_data$proportion <- proportion_data$count / proportion_data$total * 100

# 自定义每个物种的颜色
species_colors <- c("human" = "#E69F00",
                    "monkey" = "#56B4E9",
                    "mouse1" = "#009E73",
                    "mouse2" = "#D55E00", 
                    "mouse3" = "#CC79A7")

# 使用 ggplot2 绘制分组条形图，显示百分比
ggplot(proportion_data, aes(x = seurat_clusters, y = proportion, 
                            fill = aged)) +
  geom_bar(stat = "identity", position = "dodge") +
  #geom_text(aes(label = sprintf("%.1f%%", proportion), y = proportion + 1),
  #  position = position_dodge(width = 0.9), vjust = 0, size = 3) +
  facet_wrap(~ species) +
  # scale_fill_manual(values = species_colors) +
  labs(title = "Percentage of Clusters in Different Species and Groups by aged",
       x = "Species",
       y = "Percentage",
       fill = "Species") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(), # 去掉主网格线
    panel.grid.minor = element_blank()  # 去掉次网格线
  )
