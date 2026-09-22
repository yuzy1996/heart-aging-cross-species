# =============================================================================
# 13_myeloid_conserved_programs.R
# Title : Conserved myeloid expression programs
# Figure: Supplementary
# Module: 03_figure2_fibroblast_states
# Description:
#   Average-expression matrix and conserved up/down programs for myeloid subtypes.
# Inputs:
#   - (see script)
# Outputs:
#   - (see script)
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

Myeloid_species_final@meta.data$new_group <- paste(Myeloid_species_final$species,
                                                   Myeloid_species_final$aged,
                                                   Myeloid_species_final$subtype,sep = '_')
Myeloid_species_final@meta.data$new_group  <- as.factor(Myeloid_species_final@meta.data$new_group)
Idents(Myeloid_species_final) <- 'new_group'
Myeloid_species_final_exp <- AverageExpression(Myeloid_species_final,assays = 'RNA')

mat = Myeloid_species_final_exp$RNA %>% as.matrix() %>% t()
mat <- mat[1:36,]
genes_conserved_up <- c(INF_MC_veen_up[1,5][[1]],
                     MHClow_MC_veen_up[1,5][[1]],
                     CCR2_MC_veen_up[1,5][[1]],
                     MHChi_MC_veen_up[1,5][[1]],
                     Mono_veen_up[1,5][[1]],
                     Timd4_resident_MC_veen_up[1,5][[1]])%>% unique()
genes_conserved_down <-c(INF_MC_veen_down[1,5][[1]],
                     MHClow_MC_veen_down[1,5][[1]],
                     CCR2_MC_veen_down[1,5][[1]],
                     MHChi_MC_veen_down[1,5][[1]],
                     Mono_veen_down[1,5][[1]],
                     Timd4_resident_MC_veen_down[1,5][[1]])%>% unique()
genes_species_huamn <- c(INF_MC_veen_up[7,5][[1]],
                         MHClow_MC_veen_up[7,5][[1]],
                         CCR2_MC_veen_up[7,5][[1]],
                         MHChi_MC_veen_up[7,5][[1]],
                         Mono_veen_up[7,5][[1]],
                         Timd4_resident_MC_veen_up[7,5][[1]],
                         INF_MC_veen_down[7,5][[1]],
                         MHClow_MC_veen_down[7,5][[1]],
                         CCR2_MC_veen_down[7,5][[1]],
                         MHChi_MC_veen_down[7,5][[1]],
                         Mono_veen_down[7,5][[1]],
                         Timd4_resident_MC_veen_down[7,5][[1]]) %>% unique()
genes_species_monkey <- c(INF_MC_veen_up[6,5][[1]],
                          MHClow_MC_veen_up[6,5][[1]],
                          CCR2_MC_veen_up[6,5][[1]],
                          MHChi_MC_veen_up[6,5][[1]],
                          Mono_veen_up[6,5][[1]],
                          Timd4_resident_MC_veen_up[6,5][[1]],
                          INF_MC_veen_down[6,5][[1]],
                          MHClow_MC_veen_down[6,5][[1]],
                          CCR2_MC_veen_down[6,5][[1]],
                          MHChi_MC_veen_down[6,5][[1]],
                          Mono_veen_down[6,5][[1]],
                          Timd4_resident_MC_veen_down[6,5][[1]])%>% unique()
genes_species_mouse <- c(INF_MC_veen_up[4,5][[1]],
                         MHClow_MC_veen_up[4,5][[1]],
                         CCR2_MC_veen_up[4,5][[1]],
                         MHChi_MC_veen_up[4,5][[1]],
                         Mono_veen_up[4,5][[1]],
                         Timd4_resident_MC_veen_up[4,5][[1]],
                         INF_MC_veen_down[4,5][[1]],
                         MHClow_MC_veen_down[4,5][[1]],
                         CCR2_MC_veen_down[4,5][[1]],
                         MHChi_MC_veen_down[4,5][[1]],
                         Mono_veen_down[4,5][[1]],
                         Timd4_resident_MC_veen_down[4,5][[1]]) %>% unique()
mat_corss <- mat[,unique(c(genes_conserved,genes_species_huamn,genes_species_monkey,genes_species_mouse))]
dim(mat_corss)
group = rep(c("Conserved", "Species_human","Species_monkey","Species_mouse"),
            times = c(length(c(genes_conserved_up,genes_conserved_down)),
                      length(genes_species_huamn),
                      length(genes_species_monkey),
                      length(genes_species_mouse)))
col1 = colorRamp2(c(-2, 2), c( "white", "red"))
col2 = colorRamp2(c(-2, 2), c("white", "blue"))
col3 = colorRamp2(c(-2, 2), c("white", "green"))
col4 = colorRamp2(c(-2, 2), c("white", "pink"))
clusters <- rep(c('CCR2_MC','INF_MC','MHChi_MC','MHClow_MC','Mono',"Timd4-resident-MC"), 6)
aged <- rep(c(rep('Young',6),rep('Aged',6)),3)
clusters_col <- rep(c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD","#66C2A5"), 6)
aged_col <- rep(c(rep('#7570B3',6),rep('#FFED6F',6)),3)
names(clusters_col) <- clusters
names(aged_col) <- aged
left_annotation= rowAnnotation(clusters = clusters,
                               aged=aged, 
                               col = list(clusters = clusters_col,
                                          aged=aged_col))

Heatmap(mat_corss_scale[, 1], col = col1, name = "Conserved",
        left_annotation = left_annotation,
        row_split = c(rep('human',12),rep('monkey',12),rep('mouse',12)),
        cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = F)
mat_corss_scale <- scale(mat_corss)
ht1 = Heatmap(mat_corss_scale[, group == "Conserved"], col = col1, name = "Conserved",
             left_annotation = left_annotation,
             row_split = c(rep('human',12),rep('monkey',12),rep('mouse',12)),
              cluster_rows = F,cluster_columns = F,
             show_row_names = F,show_column_names = F)
ht1
ht2 = Heatmap(mat_corss_scale[, group == "Species_human"], col = col2, name = "Species_human",
              row_split = c(rep('human',12),rep('monkey',12),rep('mouse',12)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht3 = Heatmap(mat_corss_scale[, group == "Species_monkey"], col = col3, name = "Species_monkey",
              row_split = c(rep('human',12),rep('monkey',12),rep('mouse',12)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht4 = Heatmap(mat_corss_scale[, group == "Species_mouse"], col = col4, name = "Species_mouse",
              row_split = c(rep('human',12),rep('monkey',12),rep('mouse',12)),
              cluster_rows = F,cluster_columns = F,
              show_row_names = F,show_column_names = F)
ht1 + ht2 + ht3 + ht4
