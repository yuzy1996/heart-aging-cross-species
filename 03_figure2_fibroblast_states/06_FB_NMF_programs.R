# =============================================================================
# 06_FB_NMF_programs.R
# Title : NMF fibroblast programs
# Figure: Fig. 2E; Supplementary
# Module: 03_figure2_fibroblast_states
# Description:
#   GeneNMF factorisation and UCell program scores for fibroblast states.
# Inputs:
#   - hu("downsample", "FB_label_re_cca.rds")
#   - hu("FIG3", "FB_program.rda")
# Outputs:
#   - hu("downsample", "FB_label_marker.pdf")
#   - hu("downsample", "FB_program1.pdf")
#   - hu("downsample", "FB_program.pdf")
#   - hu("downsample", "FB_program_MP7.pdf")
#   - hu("downsample", "FB_program_MP_score_fea.pdf")
#   - hu("downsample", "FB_program_MP_score_MP6.pdf")
#   - hu("downsample", "geom_box.pdf")
#   - hu("downsample", "geom_eye.pdf")
#   - hu("downsample", "FB_program.rda")
#   - hu("FIG3", "FB_program_MP_score_fea.pdf")
#   - hu("FIG3", "FB_program.rda")
#   - hu("FIG3", "FB_program_MP_score.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(tidyverse)
library(SCP)
#install.packages("GeneNMF")
#remotes::install_github("carmonalab/GeneNMF")  #from GitHub
library(GeneNMF)
FB_label <- read_rds(file=hu("downsample", "FB_label_re_cca.rds"))
FB_label[['RNA']] <- as(FB_label[['RNA']],Class = 'assay')
marker <- c('SEMA3C',
            'PI16',
            'SCN7A','LPL',
            'APOE',
            'C4A','MEOX1','FN1','THBS1','POSTN',
            'DKK3',
            'KCNMA1','HDAC9','ERBB4')
ht1 <- GroupHeatmap(FB_label,
                    features = marker,
                    show_row_names = T,
                    group.by = c("subtypes")
)
pdf(hu("downsample", "FB_label_marker.pdf"),width = 5,height = 6)
ht1
dev.off()
FB_monkey <- subset(FB_label,species%in%c('monkey')) 
FB_human <- subset(FB_label,species%in%c('human'))
FB_mouse <- subset(FB_label,dataset%in%c('mouse'))

seu.list <- SplitObject(FB_label, split.by = "orig.ident")
geneNMF.programs <- multiNMF(seu.list, 
                             assay="RNA", slot="data", k=4:12,min.exp=0.05)
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        metric = "cosine",
                                        weight.explained = 0.8,
                                        nMP=8)
geneNMF.metaprograms$metaprograms.metrics
lapply(geneNMF.metaprograms$metaprograms.genes,head)
ph <- plotMetaPrograms(geneNMF.metaprograms1,downsample=500,
                       similarity.cutoff = c(0.2,1))
geneNMF.metaprograms$metaprograms.genes
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(FB_label), category = "C5", subcategory = "GO:BP")
})
set.seed(123)

pdf(hu("downsample", "FB_program1.pdf"),width = 14,height = 12)
ph
dev.off()
ph <- plotMetaPrograms(geneNMF.metaprograms,downsample=500,gaps=0,
                       similarity.cutoff = c(0.05,1))
pdf(hu("downsample", "FB_program.pdf"),width = 7,height = 6)
ph
dev.off()

library(msigdbr)
library(fgsea)
library(viridis)
mp.genes <- geneNMF.metaprograms$metaprograms.genes
FB_label <- AddModuleScore_UCell(FB_label, features = mp.genes, ncores=4, name = "")

pdf(hu("downsample", "FB_program_MP7.pdf"),width = 12,height = 4)
SCP::FeatureDimPlot(FB_label,features = 'MP7',
                    split.by =  "aged",cells.highlight = T,
                    reduction = 'UMAP',theme_use = 'theme_blank')
dev.off()

pdf(hu("downsample", "FB_program_MP_score_fea.pdf"),width = 9,height = 3)
SCP::FeatureDimPlot(FB_label,features = 'MP1',
                    split.by =  "aged",cells.highlight = T,
                    reduction = 'UMAP',theme_use = 'theme_blank')
SCP::FeatureDimPlot(FB_label,features = 'MP2',
                    split.by =  "aged",cells.highlight = T,
                    reduction = 'UMAP',theme_use = 'theme_blank')
SCP::FeatureDimPlot(FB_label,features = 'MP3',
                    split.by =  "aged",cells.highlight = T,
                    reduction = 'UMAP',theme_use = 'theme_blank')
SCP::FeatureDimPlot(FB_label,features = 'MP4',
                    split.by =  "aged",cells.highlight = T,
                    reduction = 'UMAP',theme_use = 'theme_blank')
SCP::FeatureDimPlot(FB_label,features = 'MP5',
                    split.by =  "aged",cells.highlight = T,
                    reduction = 'UMAP',theme_use = 'theme_blank')
SCP::FeatureDimPlot(FB_label,features = 'MP6',
                    split.by =  "aged",cells.highlight = T,
                    reduction = 'UMAP',theme_use = 'theme_blank')
SCP::FeatureDimPlot(FB_label,features = 'MP7',
                    split.by =  "aged",cells.highlight = T,
                    reduction = 'UMAP',theme_use = 'theme_blank')
SCP::FeatureDimPlot(FB_label,features = 'MP8',
                    split.by =  "aged",cells.highlight = T,
                    reduction = 'UMAP',theme_use = 'theme_blank')
dev.off()

FB_human$aged <- factor(FB_human$aged,levels = c('Young','Aged'))
FB_monkey$aged <- factor(FB_monkey$aged,levels = c('Young','Aged'))
FB_mouse$aged <- factor(FB_mouse$aged,levels = c('Young','Aged'))
pdf(hu("downsample", "FB_program_MP_score_MP6.pdf"),width = 8,height = 4)
SCP::FeatureDimPlot(FB_human,features = 'MP6',
                    split.by =  "aged",cells.highlight = T,
                    reduction = 'UMAP',theme_use = 'theme_blank')
SCP::FeatureDimPlot(FB_monkey,features = 'MP6',
                    split.by =  "aged",cells.highlight = T,
                    reduction = 'UMAP',theme_use = 'theme_blank')
SCP::FeatureDimPlot(FB_mouse,features = 'MP6',
                    split.by =  "aged",cells.highlight = T,
                    reduction = 'UMAP',theme_use = 'theme_blank')
dev.off()
FB4 <- subset(FB_label,subtypes%in%'FB_4')
FB_monkey <- subset(FB_label,species%in%c('monkey')) 
FB_human <- subset(FB_label,species%in%c('human'))
FB_mouse <- subset(FB_label,species%in%c('mouse'))
Data_human <- data.frame(Group='MP6',
                    Attribute=FB_human$aged,
                    value=FB_human$MP6)
Data_monkey <- data.frame(Group='MP6',
                         Attribute=FB_monkey$aged,
                         value=FB_monkey$MP6)
Data_mouse <- data.frame(Group='MP6',
                         Attribute=FB_mouse$aged,
                         value=FB_mouse$MP6)
Data_human$Attribute <- factor(Data_human$Attribute,levels = c('Young','Aged'))
Data_monkey$Attribute <- factor(Data_monkey$Attribute,levels = c('Young','Aged'))
Data_mouse$Attribute <- factor(Data_mouse$Attribute,levels = c('Young','Aged'))
pdf(hu("downsample", "geom_box.pdf"),width = 5,height = 5)
plot_eye(Data = Data_human,Attribute1 = 'Young',Attribute2 = 'Aged',label='M6')
plot_eye(Data = Data_monkey,Attribute1 = 'Young',Attribute2 = 'Aged',label='M6')
plot_eye(Data = Data_mouse,Attribute1 = 'Young',Attribute2 = 'Aged',label='M6') 
dev.off()

pdf(hu("downsample", "geom_eye.pdf"),width = 5,height = 2.5)
plot_eye(Data = Data1,Attribute1 = 'Young',Attribute2 = 'Aged')
plot_eye(Data = Data1,Attribute1 = 'Young',Attribute2 = 'Aged')
dev.off()
save(list=ls(),file=hu("downsample", "FB_program.rda"))
#依据分组对vale进行统计
library(ggdist)
plot_eye <- function(Data,Attribute1,Attribute2,label='Senescence Burden'){
  Data_summary <- summarySE(Data, measurevar="value", groupvars=c("Group","Attribute"))
  p <- ggplot(Data_summary, aes(x=Group, y=value,fill=Attribute)) + 
    stat_halfeye(data = Data[Data$Attribute == Attribute1,],aes(alpha = stat(f)), 
                 show_interval = FALSE, trim=FALSE,side = "left", width = 0.4) +
    stat_halfeye(data = Data[Data$Attribute == Attribute2,],aes(alpha = stat(f)), 
                 show_interval = FALSE, trim=FALSE,side = "right", width = 0.4) +
    geom_point(data = Data_summary,aes(x=Group, y=value),pch=19,position=position_dodge(0.5),size=1.5)+ #绘制均值为点图
    geom_errorbar(data = Data_summary,aes(ymin = value-sd, ymax=value+sd), #误差条表示均值±标准差
                  width=0.1, #误差条末端短横线的宽度
                  position=position_dodge(0.5), 
                  color="black",
                  alpha = 0.7,
                  size=0.5) +
    scale_fill_manual(values = c("#007B66", "#A52CA1"))+ #设置填充的颜色
    theme_classic()+ylab(label)+xlab('')
  return(p)
}
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

FB_label$age <- factor(FB_label$age,levels = c('Young','Aged'))
FB_label$species_age_group <- paste(FB_label$species,FB_label$aged,sep = '_')
FB_label$species_age_group <- factor(FB_label$species_age_group,
                                     levels = c('human_Young','human_Aged',
                                                'monkey_Young','monkey_Aged',
                                                'mouse_Young','mouse_Aged'))
FB_label$species <- factor(FB_label$species,levels = c('human','monkey','mouse'))

pdf(hu("FIG3", "FB_program_MP_score_fea.pdf"),width = 10,height = 6)
SCP::FeatureStatPlot(FB_label,
                     group.by = 'species_age_group',
                     split.by = 'species',
                     stat.by = c('MP7','MP6'),comparisons = list(c('human_Young','human_Aged'),
                                                                 c('monkey_Young','monkey_Aged'),
                                                                 c('mouse_Young','mouse_Aged')),
                     add_box = T,palcolor = c("#FC8D62","#8DD3C7","#BC80BD"),
                     add_trend = T)
dev.off()
save(list=ls(),file = hu("FIG3", "FB_program.rda"))
load(hu("FIG3", "FB_program.rda"))
matrix <- FB_label@meta.data[,names(mp.genes)]
#dimred <- scale(matrix)
dimred <- as.matrix(matrix)
colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
#New dim reduction
FB_label@reductions[["MPsignatures"]] <- new("DimReduc",
                                             cell.embeddings = dimred,
                                             assay.used = "RNA",
                                             key = "MP_",
                                             global = FALSE)
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(FB_label), category = "C5", subcategory = "GO:BP")
})
set.seed(123)
FB_label <- RunUMAP(FB_label, reduction="MPsignatures",
                    dims=1:length(FB_label@reductions[["MPsignatures"]]),
                    metric = "euclidean", reduction.name = "umap_MP")

pdf(hu("FIG3", "FB_program_MP_score.pdf"),width = 16,height = 16)
FeaturePlot(FB_label, features = names(mp.genes), reduction = "umap_MP", ncol=4) &
  scale_color_viridis(option="B") &
  theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank())
dev.off()