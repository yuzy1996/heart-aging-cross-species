# =============================================================================
# 17_PPI_network_biomart.R
# Title : Build PPI network from STRING/BioMart
# Figure: Fig. 2I
# Module: 03_figure2_fibroblast_states
# Description:
#   Map genes to proteins via BioMart and build/plot the igraph protein-protein interaction network.
# Inputs:
#   - hu("data", "ensembl_Grch38_gene_protein.txt")
# Outputs:
#   - dot.pdf
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

# remove the Taxonomy ID from the protein names
ppi <- igraph::set_vertex_attr(
  ppi, 'name', 
  value=stringr::str_replace(names(V(ppi)), '9606.', '')
)

# subset the table for genes that are in the Seurat obj & have a corresponding 
# protein entry in the PPI network
# load the table from BioMart
ens_df <- read.table(hu("data", "ensembl_Grch38_gene_protein.txt"), sep='\t', header=1)
head(ens_df)
ens_df <- subset(
  ens_df, 
  Gene.name %in% rownames(seurat_obj) & 
    Protein.stable.ID %in% names(V(ppi))
)

# subset the PPI network to just proteins found in our Ensembl data
ppi <- igraph::subgraph(
  ppi,
  vids = unique(ens_df$Protein.stable.ID)
)

# match the protein names w/ gene names, update the PPI network
ix <- match(names(V(ppi)), ens_df$Protein.stable.ID)
ppi <- igraph::set_vertex_attr(
  ppi, 'name', 
  value=ens_df$Gene.name[ix]
)
head(names(V(ppi)))
# keep genes that are in the co-expression network and in the PPI network
keep_genes <- GetWGCNAGenes(seurat_obj)
keep_genes <- keep_genes[keep_genes %in% names(V(ppi))]

# subset the PPI network
ppi <- igraph::subgraph(
  ppi,
  vids = keep_genes
)

# get the adjacency matrix representation of the PPI and
# binarize the matrix so 0 means no PPI, 1 means PPI 
ppi_mat <- igraph::as_adjacency_matrix(ppi, sparse=TRUE)
ppi_mat[ppi_mat > 0] <- 1

# get the coex network (TOM), 
TOM <- GetTOM(seurat_obj)
TOM <- TOM[rownames(ppi_mat), colnames(ppi_mat)]

seurat_obj_human <- subset(seurat_obj,species%in%'human')
seurat_obj_mouse <- subset(seurat_obj,species%in%'mouse')
seurat_obj_monkey <- subset(seurat_obj,species%in%'monkey')

hub_df <- GetHubGenes(seurat_obj, n_hubs = 300) %>% subset(module%in%'FB_4-M1')

Ribo_gene <- c('RPS8','RPS12','RPL23','TPT1','RPL19','RPLP1','RPL35A','RPL28','RPS29','RPL36')
celltype_human_FB <- filter(celltype_human_markers_down,cluster%in%'FB'&avg_log2FC<(-0.3))
Ribo_gene_s <- intersect(celltype_human_FB$gene,hub_df$gene_name) %>% intersect(ribosome_GS$gene_symbol)
  
seurat_obj$species_age <- paste(seurat_obj$species,seurat_obj$aged,sep = '_')
pdf('dot.pdf',width = 18,height = 2.5)
DotPlot(seurat_obj_human,
        group.by = 'aged',Ribo_gene_s,
        cols = c('#DC5450','#08C0F9'))|
  DotPlot(seurat_obj_monkey,
          group.by = 'aged',Ribo_gene_s,
          cols = c('#DC5450','#08C0F9'))|
  DotPlot(seurat_obj_mouse,
          group.by = 'aged',Ribo_gene_s,
          cols = c('#DC5450','#08C0F9'))
dev.off()
# get modules and from the seurat obj
modules <- GetModules(seurat_obj_human) %>% 
  subset(module != 'grey' & gene_name %in% keep_genes) %>%  
  mutate(module = droplevels(module))
mods <- levels(modules$module)

# get module colors for plotting 
mod_colors <- dplyr::select(modules, c(module, color)) %>% distinct()
mod_cp <- mod_colors$color; names(mod_cp) <- as.character(mod_colors$module)

# select the module to plot
cur_mod <- 'INH-M9'

# plot the original co-expression network 
tmp <- ModuleTopologyHeatmap(
  seurat_obj,
  mod = cur_mod,
  matrix = TOM,
  matrix_name = 'TOM',
  order_by = 'kME',
  high_color = mod_cp[cur_mod],
  type = 'unsigned',
  return_genes = TRUE # select this option to get the gene order
) 

# get the plot from the output list
p1 <- tmp[[1]] + ggtitle('Co-expression only') + NoLegend()

# get the gene order from the output list
gene_list <- tmp[[2]]

# plot the integrated network
p2 <- ModuleTopologyHeatmap(
  seurat_obj,
  mod = cur_mod,
  matrix = TOM_ppi,
  matrix_name = 'TOM',
  order_by = 'kME',
  high_color = mod_cp[cur_mod],
  type = 'unsigned',
  plot_max = tmp$plot_max,
  genes_order = gene_list
) + ggtitle('Co-expression + PPI')

# combine with patchwork
(p1 | p2) + plot_annotation(title =cur_mod) & theme(plot.title=element_text(hjust=0.5))