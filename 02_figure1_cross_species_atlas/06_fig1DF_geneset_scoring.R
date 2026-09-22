# =============================================================================
# 06_fig1DF_geneset_scoring.R
# Title : Senescence/ECM gene-set scoring
# Figure: Fig. 1D, 1F
# Module: 02_figure1_cross_species_atlas
# Description:
#   Read published aging gene sets and msigDB ECM genes, score nuclei (AddModuleScore/UCell) and plot gene-set enrichment.
# Inputs:
#   - hu("44161_2025_635_MOESM2_ESM.xlsx")
#   - hu("res", "combine_species_new_label.rds")
# Outputs:
#   - hu("res", "seuratObj_10pct.rds")
#   - hu("FIG1", "fig1_F.pdf")
#   - hu("FIG1", "Fig1_F_v2.pdf")
#   - hu("FIG1", "geneset_score.pdf")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

library(readxl)
gene_set <- readxl::read_xlsx(hu("44161_2025_635_MOESM2_ESM.xlsx"),sheet = 'Table 6',col_names = T) %>% as.data.frame()
gene_set <- gene_set[-1,]
gene_set_list <- split(gene_set[,1],gene_set[,2])
gene_set_list <- gene_set_list[c(5,3,4,2)]
terms_C5_ECM <- msigdbr::msigdbr(species = "human", category = c("C5")) %>%
  dplyr::distinct(gs_name, gene_symbol) %>%
  dplyr::mutate(gs_name=tolower(.$gs_name)) %>%
  dplyr::filter(.,str_detect(.$gs_name,pattern =c("gocc_collagen_containing_extracellular_matrix")))
gene_set_list[['ECM']] <- terms_C5_ECM$gene_symbol



mitochondrial_dysfunction <- c('POLG', 'TFAM', 'SIRT3', 'MT-RNR2' , 'MT-RNR4')
stem_cell_exhaustion <- c('POU5F1', 'SOX2', 'KLF4', 'MYC', 'FOXM1')
intercellular_communication <- c('CCL11', 'TGFB1', 'IL6', 'GDF11', 'VEGF', 
                                 'YAP1', 'WWTR1', 'TIMP2', 'IL37','TIMP1')
genomic_instability <- c('BUB1B', 'SIRT6', 'OGG1') 
telomere_attrition <- c()
epigenetic_alterations<- c('SIRT1', 'SIRT3', 'SIRT6', 'SIRT7', 'KAT7', 'MIR188', 'MIR455', 'LINE1')
loss_of_proteostasis <- c('HSPA1A' , 'LAMP2A', 'EIF2A', 'RPS23', 'RPS9')
disabled_macroautoph <- c('ATG5', 'BECN1', 'EP300')
deregulated_nutrient_sensing <- c('IGF1', 'IGF1R', 'PIK3CA', 'AKT1', 'MTOR', 'FOXO3', 'ALK')



combine_species <- read_rds(hu("res", "combine_species_new_label.rds"))
set.seed(1234)
seuratObj_10pct <- subset(combine_species, cells = sample(colnames(combine_species), 0.1*ncol(combine_species)))
write_rds(seuratObj_10pct,hu("res", "seuratObj_10pct.rds"))
human <- subset(combine_species,species%in%'human')
monkey<- subset(combine_species,species%in%'monkey')
mouse <- subset(combine_species,species%in%'mouse')

human <- AddModuleScore(human,features = Aging_hallmarker,name =names(Aging_hallmarker) )
monkey <- AddModuleScore(monkey,features = Aging_hallmarker,name =names(Aging_hallmarker) )
mouse <- AddModuleScore(mouse,features = Aging_hallmarker,name =names(Aging_hallmarker) )

geneset <- rbind(human@meta.data,monkey@meta.data,mouse@meta.data)
geneset <- geneset[,c(25,26,36,38:51)]
colnames(geneset)
colnames(geneset)[4:17] <- c("Telomere attrition",                
                      "Genomic instability",
                      "Mitochondrial dysfunction",
                      "Stem cell exhaustion",
                      "Altered intercellular communication",
                      "Dysbiosis related",
                      "Disabled macroautophagy",
                      "Epigenetic",
                      "Deregulated nutrient sensing",
                      "Loss of proteostasis",
                      "SASP",
                      "Inflammation","Cardiac fibrosis","ECM")
result <- geneset%>%pivot_longer(.,colnames(geneset)[4:17], names_to = "term")%>%
  dplyr::group_by(species,aged,term) %>% 
  dplyr::summarise(average = mean(value))
result$aged <- factor(result$aged,levels = c('Young','Aged'))
result_main <- dplyr::filter(result,term%in%c(
  'Cardiac fibrosis',
  'ECM',
  'Inflammation',
  'SASP',
  'Genomic instability',
  'Loss of proteostasis'
))
result_main$term <- factor(result_main$term,levels = 
                            c('SASP',
                              'Inflammation',
                              'Cardiac fibrosis',
                           'ECM',
                           'Loss of proteostasis',
                           'Genomic instability'
))
pdf(hu("FIG1", "fig1_F.pdf"),width=7,height = 3)
ggplot(data=result_main, aes(x=aged, y=average, group=species))+
  geom_line(aes(color=species))+
  geom_point(aes(color=species,shape= species),size=3)+
  scale_color_manual(values=c("#FC8D62","#8DD3C7","#BC80BD"))+
  theme_bw()+
  scale_x_discrete(expand = expansion(add = c(0.2, 0.2)))+
  #theme_minimal(base_size = 10)+
  ylab('Enrichment Score')+facet_grid(~term, scales = "free_y")
dev.off()



result_human <- geneset %>%pivot_longer(.,colnames(geneset)[c(4,12,13,14,15,16)], names_to = "term") %>% 
  dplyr::group_by(species,annotation,aged,term) %>%
  dplyr::summarise(average = mean(value)) %>% 
  pivot_wider(names_from = aged, values_from = average) %>%
  mutate(difference = Aged - Young) %>%   dplyr::filter(species=='human') %>% 
  .[,c('term','annotation','difference')] %>%
  tidyr::pivot_wider(names_from = term, values_from = difference) %>% 
  # 将SID列设为行名，形成真正的矩阵
  column_to_rownames(var = "annotation") %>%
  as.matrix()


result_human <- geneset %>%pivot_longer(.,colnames(geneset)[4:9], names_to = "term") %>% 
  dplyr::group_by(species,annotation,aged,term) %>%
  dplyr::summarise(average = mean(value)) %>% 
  pivot_wider(names_from = aged, values_from = average) %>%
  mutate(difference = Aged - Young) %>%   dplyr::filter(species=='human') %>% 
  .[,c('term','annotation','difference')] %>%
  tidyr::pivot_wider(names_from = term, values_from = difference) %>% 
  # 将SID列设为行名，形成真正的矩阵
  column_to_rownames(var = "annotation") %>%
  as.matrix()

result_monkey <- geneset %>%pivot_longer(.,colnames(geneset)[4:9], names_to = "term") %>% 
  dplyr::group_by(species,annotation,aged,term) %>%
  dplyr::summarise(average = mean(value)) %>% 
  pivot_wider(names_from = aged, values_from = average) %>%
  mutate(difference = Aged - Young) %>%   dplyr::filter(species=='monkey') %>% 
  .[,c('term','annotation','difference')] %>%
  tidyr::pivot_wider(names_from = term, values_from = difference) %>% 
  # 将SID列设为行名，形成真正的矩阵
  column_to_rownames(var = "annotation") %>%
  as.matrix()

result_mouse <- geneset %>%pivot_longer(.,colnames(geneset)[4:9], names_to = "term") %>% 
  dplyr::group_by(species,annotation,aged,term) %>%
  dplyr::summarise(average = mean(value)) %>% 
  pivot_wider(names_from = aged, values_from = average) %>%
  mutate(difference = Aged - Young) %>%   dplyr::filter(species=='mouse') %>% 
  .[,c('term','annotation','difference')] %>%
  tidyr::pivot_wider(names_from = term, values_from = difference) %>% 
  # 将SID列设为行名，形成真正的矩阵
  column_to_rownames(var = "annotation") %>%
  as.matrix()

geneset_data_human <- geneset %>%pivot_longer(.,colnames(geneset)[c(5,13,14,15,16,17)], names_to = "term") %>% 
  dplyr::filter(species=='human'&aged%in%'Aged') %>%
  dplyr::group_by(annotation,term) %>%
  dplyr::summarise(average = mean(value)) %>%  
  pivot_wider(names_from = term, values_from = average) %>%
  column_to_rownames(var = "annotation") %>%
  as.matrix()

geneset_data_monkey <- geneset %>%pivot_longer(.,colnames(geneset)[c(5,13,14,15,16,17)], names_to = "term") %>% 
  dplyr::filter(species=='monkey'&aged%in%'Aged') %>%
  dplyr::group_by(annotation,term) %>%
  dplyr::summarise(average = mean(value)) %>%  
  pivot_wider(names_from = term, values_from = average) %>%
  column_to_rownames(var = "annotation") %>%
  as.matrix()

geneset_data_mouse <- geneset %>%pivot_longer(.,colnames(geneset)[c(5,13,14,15,16,17)], names_to = "term") %>% 
  dplyr::filter(species=='mouse'&aged%in%'Aged') %>%
  dplyr::group_by(annotation,term) %>%
  dplyr::summarise(average = mean(value)) %>%  
  pivot_wider(names_from = term, values_from = average) %>%
  column_to_rownames(var = "annotation") %>%
  as.matrix()

col1 = colorRamp2(c(min(geneset_data_human), 
                    max(geneset_data_human)), c( "white", "#FC8D62"))
col2 = colorRamp2(c(min(geneset_data_monkey),
                    max(geneset_data_monkey)), c("white", "#8DD3C7"))
col3 = colorRamp2(c(min(geneset_data_mouse),
                    max(geneset_data_mouse)), c("white", "#BC80BD"))

ht1 <- ComplexHeatmap::Heatmap(geneset_data_human[,c(6,4,1,2,3,5)],cluster_rows = F,
                               column_names_rot = 45,border = TRUE,
                               cluster_columns = F,col=col1,name = 'human') 
ht2 <- ComplexHeatmap::Heatmap(geneset_data_monkey[,c(6,4,1,2,3,5)],cluster_rows = F,
                               column_names_rot = 45,border = TRUE,
                               cluster_columns = F,col=col2,name = 'monkey')
ht3 <- ComplexHeatmap::Heatmap(geneset_data_mouse[,c(6,4,1,2,3,5)],cluster_rows = F,
                               column_names_rot = 45,border = TRUE,
                               cluster_columns = F,col=col3,name = 'mouse')
vertical_combine <- ht1 %v% ht2 %v% ht3
pdf(hu("FIG1", "Fig1_F_v2.pdf"),width = 4.5,height = 8)
draw(vertical_combine,padding = unit(c(25, 25, 25, 25),"mm"))
dev.off()


seurat_integrated_sub <- AddModuleScore(seurat_integrated_sub,features = gene_set_list,name =names(gene_set_list) )

VlnPlot(seurat_integrated_sub,group.by='group',split.by='')

plot_eye <- function(Data,Attribute1,Attribute2,label,col){
  Data_summary <- summarySE(Data, measurevar="value", groupvars=c("Group","Attribute"))
  p <- ggplot(Data, aes(x=Group, y=value,fill=Attribute)) + 
    stat_halfeye(data = Data[Data$Attribute == Attribute1,],
                 #aes(alpha = stat(f)), 
                 show_interval = FALSE, trim=FALSE,side = "left", width = 0.4) +
    stat_halfeye(data = Data[Data$Attribute == Attribute2,],
                 #aes(alpha = stat(f)), 
                 show_interval = FALSE, trim=FALSE,side = "right", width = 0.4) +
    geom_point(data = Data_summary,aes(x=Group, y=value),pch=19,position=position_dodge(0.5),size=1.5)+ #绘制均值为点图
    geom_errorbar(data = Data_summary,aes(ymin = value-sd, ymax=value+sd), #误差条表示均值±标准差
                  width=0.1, #误差条末端短横线的宽度
                  position=position_dodge(0.5), 
                  color="black",
                  alpha = 0.7,
                  size=0.5) +
    scale_fill_manual(values = col)+ #设置填充的颜色
    theme_classic()+ylab(label)+xlab('')
  return(p)
}
df <- data.frame(
  CellType = combine_species_new_label$annotation, 
  aged = combine_species_new_label$aged,
  species = combine_species_new_label$species,
  Senescence = combine_species_new_label$`Aging-related gene sets1`
)
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
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
library(ggdist)
Data1_human <- data.frame(Group=human$annotation,
                    Attribute=factor(human$aged,levels = c('Young','Aged')),
                    species=human$species,
                    value=human$`SASP gene set5`)
Data2_human <- data.frame(Group=human$annotation,
                    Attribute=factor(human$aged,levels = c('Young','Aged')),
                    species=human$species,
                    value=human$`Cardiac fibrosis-related genes2`)
Data1_monkey <- data.frame(Group=monkey$annotation,
                          Attribute=factor(monkey$aged,levels = c('Young','Aged')),
                          species=monkey$species,
                          value=monkey$`SASP gene set5`)
Data2_monkey <- data.frame(Group=monkey$annotation,
                          Attribute=factor(monkey$aged,levels = c('Young','Aged')),
                          species=monkey$species,
                          value=monkey$`Cardiac fibrosis-related genes2`)
Data1_mouse <- data.frame(Group=mouse$annotation,
                          Attribute=factor(mouse$aged,levels = c('Young','Aged')),
                          species=mouse$species,
                          value=mouse$`SASP gene set5`)
Data2_mouse <- data.frame(Group=mouse$annotation,
                          Attribute=factor(mouse$aged,levels = c('Young','Aged')),
                          species=mouse$species,
                          value=mouse$`Cardiac fibrosis-related genes2`)
col1 <- c('#73AE8C','#C66FAD')
col2 <- c("#A1D5F2", "#EEC25F")
p_human_1 <- plot_eye(Data = subset(Data1,species%in%'human'),
               Attribute1 = 'Young',Attribute2 = 'Aged',label='SASP-gene score',col = col1)
p_human_2 <- plot_eye(Data = subset(Data2,species%in%'human'),
               Attribute1 = 'Young',Attribute2 = 'Aged',label='Cardiac-fibrosis-gene score',col = col2)

p_monkey_1 <- plot_eye(Data = subset(Data1,species%in%'monkey'),
                      Attribute1 = 'Young',Attribute2 = 'Aged',label='SASP-gene score',col = col1)
p_monkey_2 <- plot_eye(Data = subset(Data2,species%in%'monkey'),
                      Attribute1 = 'Young',Attribute2 = 'Aged',label='Cardiac-fibrosis-gene score',col = col2)
p_monkey_2
p_mouse_1 <- plot_eye(Data = subset(Data1,species%in%'mouse'),
                      Attribute1 = 'Young',Attribute2 = 'Aged',label='SASP-gene score',col = col1)
p_mouse_2 <- plot_eye(Data = subset(Data2,species%in%'mouse'),
                      Attribute1 = 'Young',Attribute2 = 'Aged',label='Cardiac-fibrosis-gene score',col = col2)
pdf(hu("FIG1", "geneset_score.pdf"),width = 10,height = 10)
p_human_1/p_monkey_1/p_mouse_1+
p_human_2/p_monkey_2/p_mouse_2
dev.off()
