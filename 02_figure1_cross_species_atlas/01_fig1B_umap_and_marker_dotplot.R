# =============================================================================
# 01_fig1B_umap_and_marker_dotplot.R
# Title : Integrated UMAP and marker dot plot
# Figure: Fig. 1B
# Module: 02_figure1_cross_species_atlas
# Description:
#   UMAP of integrated nuclei coloured by cell type/species/age plus a marker-gene dot plot.
# Inputs:
#   - hu("res", "combine_species_new_label.rds")
#   - cr("human", "crossspecies", "sceobj", "combine_species_batch_noCM_V2_sub.rds")
#   - hu("fig", "V2", "combine_species.rds")
#   - hu("fig", "V2", "combine_species_harmony.rds")
# Outputs:
#   - hu("FIG1", "fig1_umap_all.pdf")
#   - hu("fig", "V2", "fig1_umap_aged_split.pdf")
#   - hu("FIG1", "fig1_dot_marker.pdf")
#   - hu("FIG1", "fig1_fea_marker.pdf")
#   - hu("fig", "V2", "fig1_umap_species_batch.pdf")
#   - hu("fig", "V2", "combine_species.rds")
#   - hu("fig", "V2", "combine_species_cca.rds")
#   - hu("fig", "V2", "combine_species_mnn.rds")
#   - hu("fig", "V2", "combine_species_harmony.rds")
#   - hu("fig", "V2", "fig1_umap_mnn.pdf")
#   - hu("fig", "V2", "fig1_umap_harmony.pdf")
#   - hu("fig", "V2", "fig1_umap_batch.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

library(tidyverse)
library(SCP)
library(SeuratWrappers)
library(Seurat)
library(SeuratData)
celltypes_col <- c("#FF1493", "#FFB6C1", "#483D8B", "#00BFFF", "#228B22", "#7FFF00", "#FFFF00", "#FF4500","#8DD3C7")
species_col <- c("#FC8D62","#8DD3C7","#BC80BD")
aged_coll <- c("#483D8B","#00BFFF")
celltypes <- c('FB','EC','Myeloid','T','B','Pericytes','SMC','Neural','Adipocyte')
species <- c('human','monkey','mouse')
aged <- c('Aged','Young')
combine_species_new_label <- read_rds(hu("res", "combine_species_new_label.rds"))
combine_species_new_label_S4 <- combine_species_new_label
combine_species_new_label_S4[['RNA']] <- as(combine_species_new_label_S4[['RNA']],Class = 'Assay')


pdf(hu("FIG1", "fig1_umap_all.pdf"),width = 5,height = 4)
DimPlot(combine_species_new_label,group.by = 'annotation',
        label = T,cols = celltypes_col) 
DimPlot(combine_species_new_label,group.by = 'species',
        label = T,cols = species_col) 
DimPlot(combine_species_new_label,group.by = 'aged',
        label = F,cols = aged_coll) 
dev.off()
pdf(hu("fig", "V2", "fig1_umap_aged_split.pdf"),width = 8,height = 4)
DimPlot(combine_species_new_label,group.by = 'annotation',split.by = "aged",cols = species_col,
        label = T,cols = celltypes_col) 
dev.off()

pdf(hu("fig", "V2", "fig1_umap_aged_split.pdf"),width = 8,height = 4)
DimPlot(combine_species_new_label,group.by = 'species',reduction = '',label = T) 
dev.off()
pdf(hu("fig", "V2", "fig1_umap_aged_split.pdf"),width = 8,height = 4)
DimPlot(combine_species_new_label,group.by = 'annotation',split.by = "aged",
        label = T,cols = celltypes_col) 
dev.off()

pdf(hu("FIG1", "fig1_dot_marker.pdf"),width = 6,height = 9)
SCP::GroupHeatmap(srt = combine_species_new_label_S4,
                  features =  c("DCN","COL1A1","GSN","FBLN1","TCF21",#FB
                                "CDH5","PECAM1","RAMP2","EMCN","TEK",#EC
                                "CD163","CD68","MS4A6A","CSF1R",#macrophage
                                "CD3D","CD3E","CD8A","IL7R",#T
                                "IGKC","MS4A1","CD19","CD79A",#B
                                "RGS5","ABCC9","KCNJ8",#Pericytes
                                "MYH11","TAGLN","ACTA2",#SMC
                                "PLP1","NRXN1","NRXN3",#Neural
                                "ADIPOQ","PLIN1","GPD1","PNPLA3","PCK1"#Adipocyte
                  ),
                  # cell_annotation = c('nFeature_RNA','Phase'),
                  group.by = 'annotation',heatmap_palette = "YlOrRd",
                  group_palcolor  = list(celltypes_col),
                  show_row_names = T,row_names_side = 'left',
                  features_label = F,
                  show_column_names = F)
dev.off()

pdf(hu("FIG1", "fig1_fea_marker.pdf"),width =25,height = 25)
FeatureDimPlot(
  srt = combine_species_new_label_S4, features = c("DCN","COL1A1","GSN","FBLN1","TCF21",#FB
                                                   "CDH5","PECAM1","RAMP2","EMCN","TEK",#EC
                                                   "CD163","CD68","MS4A6A","CSF1R",#macrophage
                                                   "CD3D","CD3E","CD8A","IL7R",#T
                                                   "IGKC","MS4A1","CD19","CD79A",#B
                                                   "RGS5","ABCC9","KCNJ8",#Pericytes
                                                   "MYH11","TAGLN","ACTA2",#SMC
                                                   "PLP1","NRXN1","NRXN3",#Neural
                                                   "ADIPOQ","PLIN1","GPD1","PNPLA3","PCK1"#Adipocyte
  ),
  reduction = "UMAP", theme_use = "theme_blank"
)%>% panel_fix(height = 2, raster = TRUE, dpi = 300)
dev.off()


combine_species_batch_noCM_V2_sub <- read_rds(cr("human", "crossspecies", "sceobj", "combine_species_batch_noCM_V2_sub.rds"))
pdf(hu("fig", "V2", "fig1_umap_species_batch.pdf"),width = 12,height = 4)
DimPlot(combine_species_batch_noCM_V2_sub,group.by = 'species',cols = species_col,reduction = '',) 
DimPlot(combine_species_batch_noCM_V2_sub,group.by = 'species',cols = species_col) 
DimPlot(combine_species_batch_noCM_V2_sub,group.by = 'species',cols = species_col) 
dev.off()

combine_species <- JoinLayers(combine_species_new_label)
write_rds(combine_species,file = hu("fig", "V2", "combine_species.rds"))

combine_species <- read_rds(hu("fig", "V2", "combine_species.rds"))
combine_species[["RNA"]] <- split(combine_species[["RNA"]], f = combine_species$orig.ident)
combine_species <- NormalizeData(combine_species)
combine_species <- FindVariableFeatures(combine_species)
combine_species <- ScaleData(combine_species)
combine_species <- RunPCA(combine_species)
combine_species <- FindNeighbors(combine_species, dims = 1:30, reduction = "pca")
combine_species <- FindClusters(combine_species, resolution = 2, cluster.name = "unintegrated_clusters")
combine_species <- RunUMAP(combine_species, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(combine_species, reduction = "umap.unintegrated", group.by = c("Method", "predicted.celltype.l2"))
combine_species_cca <- IntegrateLayers(
  object = combine_species, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
combine_species_harmony <- IntegrateLayers(
  object = combine_species, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
combine_species_mnn <- IntegrateLayers(
  object = combine_species, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)
combine_species <- IntegrateLayers(
  object = combine_species, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = SCVI_CONDA_ENV, verbose = FALSE
)
write_rds(combine_species,file = hu("fig", "V2", "combine_species_cca.rds"))

write_rds(combine_species_mnn,file = hu("fig", "V2", "combine_species_mnn.rds"))

write_rds(combine_species_harmony,file = hu("fig", "V2", "combine_species_harmony.rds"))
combine_species <- read_rds(hu("fig", "V2", "combine_species_harmony.rds"))
combine_species <- RunUMAP(combine_species_harmony, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
combine_species <- RunUMAP(combine_species, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
combine_species_mnn <- RunUMAP(combine_species_mnn, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")
pdf(hu("fig", "V2", "fig1_umap_mnn.pdf"),width = 5,height = 4)
DimPlot(combine_species_mnn,group.by = 'species',reduction = 'umap.mnn',label = T) 
dev.off()
pdf(hu("fig", "V2", "fig1_umap_harmony.pdf"),width = 5,height = 4)
DimPlot(combine_species_harmony,group.by = 'species',reduction = 'umap.harmony',label = T) 
dev.off()

pdf(hu("fig", "V2", "fig1_umap_batch.pdf"),width = 5,height = 4)
DimPlot(combine_species,group.by = 'species',reduction = 'umap.cca',label = T) 
DimPlot(combine_species,group.by = 'species',
        reduction = 'umap.unintegrated',
        cols = celltypes_col,label = T) 
dev.off()

combine_species <- JoinLayers(combine_species)
combine_species[["RNA"]] <- as(object = combine_species[["RNA"]], Class = "Assay")
SeuratDisk::SaveH5Seurat(combine_species, filename = hu("fig", "V2", "combine_species.h5Seurat"))
SeuratDisk::Convert(hu("fig", "V2", "combine_species.h5Seurat"), dest = "h5ad")

scANVIIntegration <- function(
    object,
    features = VariableFeatures(object),
    layers = "counts",
    conda_env = SCVI_CONDA_ENV,
    new.reduction = "integrated.scanvi",
    ndims = 30,
    nlayers = 2,
    gene_likelihood = "nb",
    max_epochs = NULL,
    ...) {
  
  # import python methods from specified conda env
  reticulate::use_condaenv(conda_env, required = TRUE)
  sc <- reticulate::import("scanpy", convert = FALSE)
  scvi <- reticulate::import("scvi", convert = FALSE)
  anndata <- reticulate::import("anndata", convert = FALSE)
  scipy <- reticulate::import("scipy", convert = FALSE)
  
  # if `max_epochs` is not set
  if (is.null(max_epochs)) {
    # convert `NULL` to python's `None`
    max_epochs <- reticulate::r_to_py(max_epochs)
  } else {
    # otherwise make sure it's an int
    max_epochs <- as.integer(max_epochs)
  }
  object <- read_rds(hu("fig", "V2", "combine_species.rds"))
  # build a meta.data-style data.frame indicating the batch for each cell
  batches <- .FindBatches(object, layers = layers)
  # scVI expects a single counts matrix so we'll join our layers together
  # it also expects the raw counts matrix
  # TODO: avoid hardcoding this - users can rename their layers arbitrarily
  # so there's no gauruntee that the usual naming conventions will be followed
  object <- JoinLayers(object = object, layers = "counts")
  # setup an `AnnData` python instance
  obs = object$orig.ident %>% as.data.frame()
  names(obs)='batch'
  adata <- sc$AnnData(
    X = scipy$sparse$csr_matrix(
      # TODO: avoid hardcoding per comment above
      Matrix::t(LayerData(object, layer = "counts")[features, ])
    ),
    obs = obs,
    var = features
  )
  adata$obs
  scvi$model$SCVI$setup_anndata(adata, batch_key = "batch")
  # initialize and train the model
  model <- scvi$model$SCVI(
    adata = adata,
    n_latent = as.integer(x = ndims),
    n_layers = as.integer(x = nlayers),
    gene_likelihood = gene_likelihood
  )
  scvi$model$SCANVI$setup_anndata(adata,batch_key="batch")
  # create the model
  model = scvi$model$SCANVI(adata)
  # train the model
  model$train()
  
  adata.obs["scanvi_prediction"] = model.predict()
  
  
  model$train(max_epochs = max_epochs)
  # extract the latent representation of the merged data
  latent <- model$get_latent_representation()
  latent <- adata$obsm['X_scANVI']
  latent <- as.matrix(latent)
  # pull the cell identifiers back out of the `AnnData` instance
  # in case anything was sorted under the hood
  rownames(latent) <- reticulate::py_to_r(adata$obs$index$values)
  # prepend the latent space dimensions with `new.reduction` to
  # give the features more readable names
  new.reduction='integrated.scanvi'
  colnames(latent) <- paste0(new.reduction, "_", 1:ncol(latent))
  # build a `DimReduc` instance
  suppressWarnings(
    object@reductions[new.reduction] <- CreateDimReducObject(
      embeddings = latent, 
      key = new.reduction
    )
  )
  return(object)
}

