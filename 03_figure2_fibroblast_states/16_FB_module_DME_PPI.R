# =============================================================================
# 16_FB_module_DME_PPI.R
# Title : Fibroblast module enrichment and PPI
# Figure: Fig. 2H-I
# Module: 03_figure2_fibroblast_states
# Description:
#   Differential module eigengene (DME) plots and protein-protein interaction network for module genes.
# Inputs:
#   - (see script)
# Outputs:
#   - (see script)
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

# get the modules table:
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# make a copy of the DME table for plotting
plot_df <- DMEs

# set the factor level for the modules so they plot in the right order:
plot_df$module <- factor(as.character(plot_df$module), levels=mods)

# set a min/max threshold for plotting
maxval <- 0.5; minval <- -0.5
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)

# add significance levels
plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)

# change the text color to make it easier to see 
plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.2, 'black', 'white')

# make the heatmap with geom_tile
p <- plot_df %>% 
  ggplot(aes(y=cluster, x=module, fill=avg_log2FC)) +
  geom_tile() 

# add the significance levels
p <- p + 
  geom_text(label=plot_df$Significance, color=plot_df$textcolor) 

# customize the color and theme of the plot
p <- p + 
  scale_fill_gradient2(low='purple', mid='black', high='yellow') +
  RotatedAxis() +
  theme(
    panel.border = element_rect(fill=NA, color='black', size=1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    plot.margin=margin(0,0,0,0)
  ) + xlab('') + ylab('') +
  coord_equal()

p

# convert sex to factor
seurat_obj$msex <- as.factor(seurat_obj$msex)

# convert age_death to numeric
seurat_obj$age_death <- as.numeric(seurat_obj$age_death)

# list of traits to correlate
cur_traits <- c('aged', 'nCount_RNA', 'nFeature_RNA')

seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = cur_traits,
  group.by='subtypes'
)
PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = 'yellow',
  mid_color = 'black',
  low_color = 'purple',
  plot_max = 0.2,
  combine=TRUE
)