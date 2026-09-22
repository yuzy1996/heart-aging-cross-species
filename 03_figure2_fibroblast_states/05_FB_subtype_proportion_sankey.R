# =============================================================================
# 05_FB_subtype_proportion_sankey.R
# Title : Fibroblast subtype composition Sankey
# Figure: Fig. 2 (composition)
# Module: 03_figure2_fibroblast_states
# Description:
#   Sankey/alluvial of fibroblast subtype counts across species and age.
# Inputs:
#   - (see script)
# Outputs:
#   - hu("res", "FB_sankey.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

# (R library search paths are managed by renv / .Renviron, not set per-script)
# plot sankey diagram of cluster ------------------------------------------
library(RColorBrewer)
library(ggalluvial)
new_meta <- FB_species_rpca_sub_label@meta.data
cluster_counts_FB <- new_meta %>%
  group_by(species, aged, subtype) %>%
  summarize(count = n(), .groups = 'drop')
# 确保数据格式适合绘制桑基图
cluster_counts_FB$celltypes <- factor(cluster_counts_FB$celltypes,
                                      levels = celltypes_new_level)
# 创建一个颜色映射列表，用于标签颜色
cell_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
label_colors <- c(human = "#E69F00", monkey = "#56B4E9", mouse1 = "#009E73",
                  mouse2 = "#D55E00", mouse3 = "#CC79A7", 
                  Aged = "#F0E442", Young = "#0072B2",
                  FB1 = cell_colors[1], FB2 = cell_colors[2],
                  FB3 = cell_colors[3], FB4 = cell_colors[4],
                  FB5 = cell_colors[5])
pdf(hu("res", "FB_sankey.pdf"),width = 10,height = 8)
ggplot(cluster_counts_FB, aes(axis1 = species, axis2 = subtype,axis3 = aged, y = count)) +
  geom_alluvium(aes(fill = aged), width = 0.1) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.1, color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),color = 'black',
            size = 4, vjust = -0.5) +
  scale_x_discrete(limits = c("Species", "celltypes", "aged"), expand = c(0.15, 0.05)) +
  scale_fill_manual(values = label_colors) +
  scale_color_manual(values = label_colors) +
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Sankey Diagram of Clusters Across Species and Groups",
       fill = "Cluster",
       color = "Group")
dev.off()
