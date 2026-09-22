# =============================================================================
# 08_fig1F_heatmap.R
# Title : Aging-regulated gene heat map
# Figure: Fig. 1F
# Module: 02_figure1_cross_species_atlas
# Description:
#   Heat map of z-scored aging-regulated gene expression by species and age.
# Inputs:
#   - (see script)
# Outputs:
#   - hu("FIG1", "Fig1_F_v2.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

# 加载必要的R包
library(tidyverse)
library(viridis)
library(reshape2)
library(circlize)
library(ComplexHeatmap)


geneset %>%pivot_longer(.,colnames(geneset)[4:8], names_to = "term") %>% 
  dplyr::filter(species=='human'&aged%in%'Aged') %>% 
  dplyr::group_by(annotation,term) %>%
  dplyr::summarise(average = mean(value))%>% 
  pivot_wider(names_from = term, values_from = average) %>% 
  column_to_rownames(var = "annotation") %>%
  as.matrix() ->geneset_data_human

geneset %>%pivot_longer(.,colnames(geneset)[4:8], names_to = "term") %>% 
  dplyr::filter(species=='monkey'&aged%in%'Aged') %>% 
  dplyr::group_by(annotation,term) %>%
  dplyr::summarise(average = mean(value))%>% 
  pivot_wider(names_from = term, values_from = average) %>% 
  column_to_rownames(var = "annotation") %>%
  as.matrix() ->geneset_data_monkey

geneset %>%pivot_longer(.,colnames(geneset)[4:8], names_to = "term") %>% 
  dplyr::filter(species=='mouse'&aged%in%'Aged') %>% 
  dplyr::group_by(annotation,term) %>%
  dplyr::summarise(average = mean(value))%>% 
  pivot_wider(names_from = term, values_from = average) %>% 
  column_to_rownames(var = "annotation") %>%
  as.matrix() ->geneset_data_mouse

col1 = colorRamp2(c(min(geneset_data_human), 
                    max(geneset_data_human)), c( "white", "#FC8D62"))
col2 = colorRamp2(c(min(geneset_data_monkey),
                    max(geneset_data_monkey)), c("white", "#8DD3C7"))
col3 = colorRamp2(c(min(geneset_data_mouse),
                    max(geneset_data_mouse)), c("white", "#BC80BD"))
ht1 <- ComplexHeatmap::Heatmap(geneset_data_human[,5:1],cluster_rows = F,
                               column_names_rot = 45,border = TRUE,
                               cluster_columns = F,col=col1,name = 'human') 
ht2 <- ComplexHeatmap::Heatmap(geneset_data_monkey[,5:1],cluster_rows = F,
                               column_names_rot = 45,border = TRUE,
                               cluster_columns = F,col=col2,name = 'monkey')
ht3 <- ComplexHeatmap::Heatmap(geneset_data_mouse[,5:1],cluster_rows = F,
                               column_names_rot = 45,border = TRUE,
                               cluster_columns = F,col=col3,name = 'mouse')
vertical_combine <- ht1 %v% ht2 %v% ht3
pdf(hu("FIG1", "Fig1_F_v2.pdf"),width = 4.5,height = 7)
draw(vertical_combine,padding = unit(c(25, 25, 25, 25),"mm"))
dev.off()
