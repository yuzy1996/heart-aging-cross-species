# =============================================================================
# 02_fig4B_interaction_counts.R
# Title : Interaction counts and strengths
# Figure: Fig. 4B-D
# Module: 05_figure4_communication
# Description:
#   Plot numbers of significant L-R pairs and incoming/outgoing interaction strengths from saved CellChat objects.
# Inputs:
#   - hu("downsample", "crosstalk", "allcelltypes_human_Y_obj.rds")
#   - hu("downsample", "crosstalk", "allcelltypes_human_O_obj.rds")
#   - hu("downsample", "crosstalk", "allcelltypes_monkey_Y_obj.rds")
#   - hu("downsample", "crosstalk", "allcelltypes_monkey_O_obj.rds")
#   - hu("downsample", "crosstalk", "allcelltypes_mouse_Y_obj.rds")
#   - hu("downsample", "crosstalk", "allcelltypes_mouse_O_obj.rds")
# Outputs:
#   - (see script)
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

# 1. 安装并加载所需包（若未安装）
# (package installation is handled by environment/00_setup_env.R)  if (!require("ggplot2")) install.packages("ggplot2")
# (package installation is handled by environment/00_setup_env.R)  if (!require("dplyr")) install.packages("dplyr")
# (package installation is handled by environment/00_setup_env.R)  if (!require("tidyr")) install.packages("tidyr")

library(ggplot2)
library(dplyr)
library(tidyr)
allcelltypes_human_Y_obj <- read_rds(hu("downsample", "crosstalk", "allcelltypes_human_Y_obj.rds"))
allcelltypes_human_O_obj <- read_rds(hu("downsample", "crosstalk", "allcelltypes_human_O_obj.rds"))
allcelltypes_monkey_Y_obj <- read_rds(hu("downsample", "crosstalk", "allcelltypes_monkey_Y_obj.rds"))
allcelltypes_monkey_O_obj <- read_rds(hu("downsample", "crosstalk", "allcelltypes_monkey_O_obj.rds"))
allcelltypes_mouse_Y_obj <- read_rds(hu("downsample", "crosstalk", "allcelltypes_mouse_Y_obj.rds"))
allcelltypes_mouse_O_obj <- read_rds(hu("downsample", "crosstalk", "allcelltypes_mouse_O_obj.rds"))

CellChat::setIdent(allcelltypes_mouse_Y_obj)
CellChat::subsetCellChat(allcelltypes_mouse_Y_obj,group.by = )


object.list_human <- list(Young = allcelltypes_human_Y_obj, Aged = allcelltypes_human_O_obj)
cellchat_human <- mergeCellChat(object.list_human, add.names = names(object.list_human))
object.list_monkey <- list(Young = allcelltypes_monkey_Y_obj, Aged = allcelltypes_monkey_O_obj)
cellchat_monkey <- mergeCellChat(object.list_monkey, add.names = names(object.list_monkey))
object.list_mouse <- list(Young = allcelltypes_mouse_Y_obj, Aged = allcelltypes_mouse_O_obj)
cellchat_mouse <- mergeCellChat(object.list_mouse, add.names = names(object.list_mouse))


human_count_celltype <- lapply(object.list_human, count_celltype)
plot_dot(human_count_celltype$Young)+plot_dot(human_count_celltype$Aged)

compareInteractions(cellchat_human, show.legend = F, group = c(1,2))+
  compareInteractions(cellchat_monkey, show.legend = F, group = c(1,2))+
  compareInteractions(cellchat_mouse, show.legend = F, group = c(1,2))

compareInteractions(cellchat_human, show.legend = F, group = c(1,2), measure = "weight")
compareInteractions(cellchat_monkey, show.legend = F, group = c(1,2), measure = "weight")
compareInteractions(cellchat_mouse, show.legend = F, group = c(1,2), measure = "weight")

count_celltype <- function(obj){
  df <- as.data.frame(obj@net$count)
  # 3. 统计核心指标
  sender_total <- rowSums(df)    # 每个细胞作为sender的互作总数
  receiver_total <- colSums(df)  # 每个细胞作为receiver的互作总数
  cell_types <- c("Adipocyte", "artery_ec", "B", "capillary_ec", "Cycle", "DC",
                  "FB1", "FB2", "FB3", "FB4", "FB5", 
                  "lymphatic_ec", "MC_CCR2", "MC_INF", "MC_MHChi",
                  "MC_MHClow", "MC_Resident", "Mono", "Neural", "Neutrophil",
                  "Pericytes", "SMC", "T", "vein_ec")
  
  # 4. 整理结果表
  result <- data.frame(
    Cell_Type = cell_types,
    Sender_Total = sender_total,
    Receiver_Total = receiver_total,
    Total_Interactions = sender_total + receiver_total
  ) %>%
    arrange(desc(Total_Interactions))  # 按总互作数降序排列
  return(result)
}
human_count_celltype <- lapply(object.list_human, count_celltype)

# 5. 打印统计结果
print("=== 细胞互作受配体对统计结果 ===")
print(result, row.names = FALSE)

# 6. 可视化：sender vs receiver 对比图
result_long <- result %>%
  pivot_longer(cols = c(Sender_Total, Receiver_Total),
               names_to = "Role",
               values_to = "Count") %>%
  mutate(Role = factor(Role, levels = c("Sender_Total", "Receiver_Total"),
                       labels = c("Sender", "Receiver")))

plot_dot <- function(result){
  p <- ggplot(result, aes(x = Total_Interactions, y = reorder(Cell_Type, Total_Interactions))) +
    geom_point(aes(size = Total_Interactions, color = Total_Interactions), alpha = 0.7) +
    labs(x = "# of Interaction Pairs", y = "Cell Type",
         title = "Cell Interaction",
         fill = "Role") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.y = element_text(size = 9),
          legend.position = "top")
  return(p)
}
plot_dot(human_count_celltype$Young)+plot_dot(human_count_celltype$Aged)
human_count_celltype_Y <- human_count_celltype$Young
human_count_celltype_Y$group <- 'Young'
human_count_celltype_O <- human_count_celltype$Aged
human_count_celltype_O$group <- 'Aged'
human_count_celltype_O$Total_Interactions <- -(human_count_celltype_O$Total_Interactions)
combine_human_count_celltype <- rbind(human_count_celltype_Y,human_count_celltype_O)
ggplot(combine_human_count_celltype, aes(x = Total_Interactions, 
                                         y = reorder(Cell_Type, Total_Interactions))) +
  geom_bar(aes(size=Total_Interactions, color = group), alpha = 0.7) +
  labs(x = "# of Interaction Pairs", y = "Cell Type",
       title = "Cell Interaction",
       fill = "group") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 9),
        legend.position = "top")

ggplot(combine_human_count_celltype, aes(x = Total_Interactions, y = reorder(Cell_Type, Total_Interactions), fill = group)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Aged" = "#3182bd", "Young" = "#de2d26")) +
  labs(x = "# of Interaction Pairs", y = "Cell Type",
       title = "Cell Interaction Summary: Sender vs Receiver",
       fill = "Role") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 9),
        legend.position = "top")

# 4. 绘制图形
ggplot() +
  # 水平条形图
  geom_col(
    data = combine_human_count_celltype,
    aes(x = Total_Interactions, y = Cell_Type, fill = group),
    position = "identity",  # 条形不堆叠，直接按 x 值绘制
    width = 0.6             # 调整条形宽度
  ) +
  # 添加差值点（黑色实心点）
  geom_point(
    data = df,
    aes(x = Difference, y = Cell_Type),
    color = "black",
    size = 3,
    shape = 19  # 实心圆点
  ) +
  # 设置填充颜色（蓝色=Monocyte，红色=Macrophage）
  scale_fill_manual(
    values = c("Aged" = "#3182bd", "Young" = "#de2d26"),
    labels = c("Aged", "Young")
  ) +
  # 调整 x 轴范围和刻度（匹配原图的 47,31,16,0,16,31,47）
  scale_x_continuous(
    limits = c(-47, 47),
    breaks = seq(-47, 47, 16),
    labels = abs(seq(-47, 47, 16))  # 刻度标签显示绝对值，避免负数
  ) +
  # 坐标轴和标题设置
  labs(
    x = "# of interaction pairs",
    y = "Celltypes as Senders",
    title = "CellChat Interaction Comparison: Macrophage vs Monocyte",
    fill = "Cell Type"  # 图例标题
  ) +
  # 主题美化
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # 标题居中、加粗
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    panel.grid = element_blank(),  # 去掉网格线
    panel.border = element_rect(linewidth = 0.8)  # 边框粗细
  ) +
  # 调整图例符号大小
  guides(fill = guide_legend(override.aes = list(size = 3)))
