# =============================================================================
# 09_fig1G_deg_clusters_GO_network.R
# Title : Conserved/species-enriched DEG ordering and GO network
# Figure: Fig. 1E, 1G
# Module: 02_figure1_cross_species_atlas
# Description:
#   Order age-DEGs from conserved (three-species shared) to species-enriched; GO enrichment network/bar chart for clusters C1-C3.
# Inputs:
#   - (see script)
# Outputs:
#   - hu("FIG1", "Fig1_G_network_barchart_complete.png")
#   - hu("FIG1", "Fig1_G_network_barchart_complete.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

# ========================
# 1. 环境准备与包加载
# ========================

# 安装必要的包（如果尚未安装）
# (package installation is handled by environment/00_setup_env.R)  if (!require("tidyverse")) install.packages("tidyverse")
# (package installation is handled by environment/00_setup_env.R)  if (!require("tidygraph")) install.packages("tidygraph")
# (package installation is handled by environment/00_setup_env.R)  if (!require("ggraph")) install.packages("ggraph")
# (package installation is handled by environment/00_setup_env.R)  if (!require("grid")) install.packages("grid")
# (package installation is handled by environment/00_setup_env.R)  if (!require("glue")) install.packages("glue")

# 加载包
library(tidyverse)
library(tidygraph)
library(ggraph)
library(grid)
library(glue)

# ========================
# 2. 生成模拟数据
# ========================
# 创建详细的模拟数据集
data <- tibble(
  category = character(),
  set = character(),
  number = integer()
) %>%
  # Miscellaneous Metabolism
  add_row(
    category = categories[1],
    set = c("Carbohydrates", "Amino Acids", "Lipids", "Vitamins", "Cofactors"),
    number = c(25, 18, 12, 8, 6)
  ) %>%
  # Stress Response
  add_row(
    category = categories[2],
    set = c("Oxidative", "Heat", "DNA Damage", "Acid", "Osmotic"),
    number = c(22, 15, 10, 7, 5)
  ) %>%
  # Carbon Metabolism
  add_row(
    category = categories[3],
    set = c("Glycolysis", "TCA Cycle", "Pentose Phosphate", "Fermentation"),
    number = c(20, 16, 12, 9)
  ) %>%
  # Uncharacterized
  add_row(
    category = categories[4],
    set = c("Hypothetical Proteins", "Conserved Domains"),
    number = c(30, 18)
  ) %>%
  # Cellular Processes
  add_row(
    category = categories[5],
    set = c("Cell Division", "Motility", "Sporulation", "Biofilm Formation"),
    number = c(15, 12, 10, 8)
  ) %>%
  # Prophages
  add_row(
    category = categories[6],
    set = c("Phage Structural", "Phage Regulatory", "Phage Integration"),
    number = c(10, 8, 6)
  ) %>%
  # Other
  add_row(
    category = categories[7],
    set = c("Regulatory", "Transporters", "Membrane Proteins"),
    number = c(18, 14, 10)
  ) %>%
  # Single Genes
  add_row(
    category = categories[8],
    set = c("Unique Genes"),
    number = c(25)
  ) %>%
  # Virulence
  add_row(
    category = categories[9],
    set = c("Toxins", "Adhesins", "Invasion Factors", "Immune Evasion"),
    number = c(12, 10, 8, 6)
  ) %>%
  # AA/Nucleotide Metabolism
  add_row(
    category = categories[10],
    set = c("Purine Synthesis", "Pyrimidine Synthesis", "Amino Acid Degradation"),
    number = c(15, 12, 10)
  ) %>%
  # 创建category_set列
  mutate(category_set = as.character(glue("{category}_{set}")))

# 定义颜色映射
color_mapping <- tibble(
  category = categories,
  color = c("#f76502", "#f41237", "#026026", "#2f0d06", "#0372e9", 
            "#f00a0b", "#2e3426", "#dea34a", "#aa1f19", "#5728a0")
) %>% deframe()  # 转换为命名向量

# ========================
# 3. 构建网络图数据
# ========================
plot_network <- function(data){
  
  # 创建节点数据
  nodes <- tibble(
    node = c("root", unique(data$category), unique(data$category_set))
  ) %>%
    mutate(
      levels = case_when(
        node == "root" ~ 1,
        node %in% unique(data$category) ~ 2,
        node %in% unique(data$category_set) ~ 3,
        TRUE ~ 4
      )
    ) %>%
    left_join(data, by = c("node" = "category_set")) %>%
    mutate(category = factor(category, levels = categories))
  
  # 创建边数据
  edges_level_1 <- data %>%
    distinct(category) %>%
    mutate(from = "root") %>%
    rename(to = category)
  
  edges_level_2 <- data %>%
    distinct(category, category_set) %>%
    select(from = category, to = category_set)
  
  # 创建颜色映射数据框（用于边着色）
  color_edges <- tibble(
    category = categories,
    color = c("#f76502", "#f41237", "#026026", "#2f0d06", "#0372e9",
              "#f00a0b", "#2e3426", "#dea34a", "#aa1f19", "#5728a0")
  )
  
  # 整合边数据并分配颜色
  edges <- bind_rows(edges_level_1, edges_level_2) %>%
    left_join(color_edges, by = c("to" = "category")) %>%
    left_join(color_edges, by = c("from" = "category")) %>%
    mutate(color = coalesce(color.x, color.y)) %>%
    select(-color.x, -color.y)
  
  # ========================
  # 4. 绘制网络图
  # ========================
  
  # 构建图形对象
  graph_data <- tbl_graph(nodes, edges)
  
  # 创建网络图
  p_network <- ggraph(graph_data, layout = 'dendrogram', circular = TRUE) +
    # 绘制边
    geom_edge_diagonal(aes(color = color), alpha = 0.5, edge_width = 0.8) +
    # 添加节点标签（只显示第二层节点）
    geom_node_text(
      aes(label = node, filter = levels == 2, color = node),
      size = 5, fontface = "bold", family = "sans"
    ) +
    # 添加节点点
    geom_node_point(
      aes(filter = levels == 2, color = node, fill = node),
      size = 6, alpha = 0.7, shape = 21, stroke = 1.2
    ) +
    # 设置颜色
    scale_color_manual(values = color_mapping) +
    scale_fill_manual(values = color_mapping) +
    scale_edge_color_identity() +
    # 设置主题
    theme(
      plot.margin = margin(120, 120, 120, 120),
      panel.background = element_rect(fill = "white", color = "white"),
      legend.position = "none"
    )
  return(p_network)
}

# ========================
# 5. 绘制极坐标条形图
# ========================
plot_chart <- function(data){
  # 准备条形图数据
  nodes_bar <- nodes %>%
    filter(levels == 3) %>%
    group_by(category) %>%
    arrange(category, desc(number), node) %>%
    ungroup() %>%
    mutate(ID = seq(1, nrow(.)))
  
  # 计算角度和对齐方式
  angle <- 90 - 360 * (nodes_bar$ID - 0.5) / nrow(nodes_bar)
  nodes_bar$hjust <- ifelse(angle < -90, 1, 0)
  nodes_bar$angle <- ifelse(angle < -90, angle + 180, angle)
  
  # 绘制极坐标条形图
  p_barchart <- ggplot(nodes_bar, aes(x = ID, y = number)) +
    geom_bar(
      stat = "identity",
      aes(fill = category, color = category),
      alpha = 0.7, linewidth = 0.7, width = 0.8
    ) +
    geom_text(
      aes(y = number + 1, label = set, hjust = hjust, color = category),
      angle = nodes_bar$angle, size = 3.5, family = "sans"
    ) +
    coord_polar() +
    ylim(-30, max(nodes_bar$number) + 8) +
    scale_fill_manual(values = color_mapping) +
    scale_color_manual(values = color_mapping) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = unit(c(-1, -1, -1, -1), "cm")
    )
}
p_barchart+p_barchart+p_barchart
# ========================
# 6. 组合图形
# ========================

# 使用网格系统组合图形
grid.newpage()

# 创建主视图
pushViewport(viewport(width = 0.9, height = 0.9))

# 先绘制条形图（作为背景）
print(p_barchart, vp = viewport(width = 1, height = 1))

# 再绘制网络图（叠加在中心）
print(p_network, vp = viewport(
  width = 0.65, height = 0.65,
  x = 0.5, y = 0.5
))

# ========================
# 7. 保存图形
# ========================

# 保存为高清PNG
png(hu("FIG1", "Fig1_G_network_barchart_complete.png"), width = 3000, height = 2800, res = 300)
grid.newpage()
pushViewport(viewport(width = 0.9, height = 0.9))
print(p_barchart, vp = viewport(width = 1, height = 1))
print(p_network, vp = viewport(width = 0.65, height = 0.65, x = 0.5, y = 0.5))
dev.off()

# 保存为PDF（矢量图）
pdf(hu("FIG1", "Fig1_G_network_barchart_complete.pdf"), width = 12, height = 11)
grid.newpage()
pushViewport(viewport(width = 0.9, height = 0.9))
print(p_barchart, vp = viewport(width = 1, height = 1))
print(p_network, vp = viewport(width = 0.65, height = 0.65, x = 0.5, y = 0.5))
dev.off()

cat("图形已保存为 network_barchart_complete.png 和 network_barchart_complete.pdf\n")
cat("数据统计:\n")
cat("- 总类别数:", length(unique(data$category)), "\n")
cat("- 总子集数:", nrow(data), "\n")
cat("- 总基因数:", sum(data$number), "\n")


head(geneset)

geneset %>%pivot_longer(.,colnames(geneset)[4:8], names_to = "term") %>% 
  dplyr::group_by(species,annotation,aged,term) %>%
  dplyr::summarise(average = mean(value)) %>% 
  pivot_wider(names_from = aged, values_from = average) %>%
  mutate(difference = Aged - Young) -> data

data <- data[,c('species','annotation','term','difference')]
colnames(data) <- c('category','set','term','number')
data%>%
  # 创建category_set列
  mutate(category_set = as.character(glue::glue("{category}_{set}")))->data

data$number <- data$number*1000
data <- data[135:1,]

# 定义功能类别
categories <- c(
  "mouse",
  "monkey", 
  "human")

# 定义颜色映射
color_mapping <- tibble(
  category = categories,
  color =c("#BC80BD","#8DD3C7","#FC8D62")
) %>% deframe()  # 转换为命名向量

# ========================
# 3. 构建网络图数据
# ========================
plot_network <- function(data){
# 创建边数据
edges_level_1 <- data %>%
  distinct(category) %>%
  mutate(from = "root") %>%
  rename(to = category)

edges_level_2 <- data %>%
  distinct(category, category_set) %>%
  select(from = category, to = category_set)

# 创建颜色映射数据框（用于边着色）
color_edges <- tibble(
  category = categories,
  color = c("#BC80BD","#8DD3C7","#FC8D62")
)

# 整合边数据并分配颜色
edges <- bind_rows(edges_level_1, edges_level_2) %>%
  left_join(color_edges, by = c("to" = "category")) %>%
  left_join(color_edges, by = c("from" = "category")) %>%
  mutate(color = coalesce(color.x, color.y)) %>%
  select(-color.x, -color.y)

# ========================
# 4. 绘制网络图
# ========================

# 构建图形对象
graph_data <- tbl_graph(nodes, edges)

# 创建网络图
p_network <- ggraph(graph_data, layout = 'dendrogram', circular = TRUE) +
  # 绘制边
  geom_edge_diagonal(aes(color = color), alpha = 0.5, edge_width = 0.8) +
  # 添加节点标签（只显示第二层节点）
  geom_node_text(
    aes(label = node, filter = levels == 2, color = node),
    size = 5, fontface = "bold", family = "sans"
  ) +
  # 添加节点点
  geom_node_point(
    aes(filter = levels == 2, color = node, fill = node),
    size = 10, alpha = 0.7, shape = 21, stroke = 1.2
  ) +
  # 设置颜色
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  scale_edge_color_identity() +
  # 设置主题
  theme(
    plot.margin = margin(120, 120, 120, 120),
    panel.background = element_rect(fill = "white", color = "white"),
    legend.position = "none"
  )
p_network}
plot_chart <- function(data,title){
  nodes <- tibble(
    node = c("root", unique(data$category), unique(data$category_set))
  ) %>%
    mutate(
      levels = case_when(
        node == "root" ~ 1,
        node %in% unique(data$category) ~ 2,
        node %in% unique(data$category_set) ~ 3,
        TRUE ~ 4
      )
    ) %>%
    left_join(data, by = c("node" = "category_set")) %>%
    mutate(category = factor(category, levels = categories))
  # 准备条形图数据
  nodes_bar <- nodes %>%
    filter(levels == 3) %>%
    group_by(category) %>%
    arrange(category, desc(number), node) %>%
    ungroup() %>%
    mutate(ID = seq(1, nrow(.)))
  
  # 计算角度和对齐方式
  angle <- 90 - 360 * (nodes_bar$ID - 0.5) / nrow(nodes_bar)
  nodes_bar$hjust <- ifelse(angle < -90, 1, 0)
  nodes_bar$angle <- ifelse(angle < -90, angle + 180, angle)
  nodes_bar$text <- nodes_bar$number
  nodes_bar$text <- ifelse(nodes_bar$text < 0, 0, nodes_bar$text)
  
  # 绘制极坐标条形图
  p_barchart <- ggplot(nodes_bar, aes(x = ID, y = number)) +
    geom_bar(
      stat = "identity",
      aes(fill = category, color = category),
      alpha = 0.7, linewidth = 0.7, width = 0.8
    ) +
    geom_text(
      aes(y = text + 1, label = set, hjust = hjust, color = category),
      angle = nodes_bar$angle, size = 6, family = "sans"
    ) +
    coord_polar() +
    ylim(-30, max(nodes_bar$number) + 8) +
    scale_fill_manual(values = color_mapping) +
    scale_color_manual(values = color_mapping) +
    theme_void() +
    theme(
      legend.position = "none",
    #  plot.margin = unit(c(-1, -1, -1, -1), "cm")
    )
  return(p_barchart)
  }

p1 <- plot_chart(filter(data,term%in%'SASP'))
p2 <- plot_chart(filter(data,term%in%'Cardiac-fibrosis'))
p3 <- plot_chart(filter(data,term%in%'Inflammation'))
p4 <- plot_chart(filter(data,term%in%'ECM'))
p5 <- plot_chart(filter(data,term%in%'NF-kappaB'))
# 保存为PDF（矢量图）
pdf(hu("FIG1", "Fig1_G_network_barchart_complete.pdf"), width = 25, height = 5)
p1|p2|p3|p4|p5
dev.off()
