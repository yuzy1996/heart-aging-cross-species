# =============================================================================
# 04_define_aging_hallmark_genesets.R
# Title : Curated aging-hallmark gene sets
# Figure: Fig. 1C-D
# Module: 02_figure1_cross_species_atlas
# Description:
#   Curated telomere/DDR/senescence and other aging-hallmark gene lists, saved as Aging_hallmarker.rda.
# Inputs:
#   - (see script)
# Outputs:
#   - hu("Aging_hallmarker.rda")
# =============================================================================
library(here)
source(here("R", "00_config.R"))
source(here("R", "utils.R"))

telomere_attrition_genes <- c(
  # 端粒保护复合体 (Shelterin Complex) 核心组分
  "TERF1", "TERF2", "POT1", "TINF2", "TPP1", "RAP1A",
  
  # 端粒酶核心组分及直接相关因子
  "TERT", "TERC", "DKC1", "NOLA1", "NOLA2", "NOLA3",
  
  # DNA损伤应答与修复相关因子
  "ATM", "ATR", "TP53", "CDKN1A", "CDKN2A", "WRN", "BLM", "PARP1", "BRCA1", "ERCC1", "XPF", "KU70", "KU80", "DNAPKcs", "53BP1",
  
  # 染色质与核膜相关因子
  "LMNA", "LMNB1", "LMNB2", "SIRT6", "SIRT1", "SIRT3",
  
  # 其他重要调控因子
  "HAPSTR1", "RTEL1", "PINX1", "TANKYRASE1", "TANKYRASE2"
)

genomic_instability_genes <- c(
  # DNA损伤应答核心激酶
  "ATM", "ATR", "PRKDC", # PRKDC即DNA-PKcs
  # 同源重组修复关键因子
  "BRCA1", "BRCA2", "RAD51", "RAD52", "PALB2", "BRIP1",
  # 非同源末端连接关键因子
  "XRCC5", "XRCC6", "LIG4", "XRCC4", "NHEJ1",
  # Fanconi贫血通路
  "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF",
  # 跨损伤合成与特殊DNA聚合酶
  "POLH", "POLK", "REV3L", "POLQ",
  # 转录-复制冲突与R-loop处理
  "SETX", "RNASEH2A", "RNASEH2B",
  # 染色质调控与基因组稳定性
  "SIRT1", "SIRT7", "KAT5",
  # 着丝粒/着丝粒功能与染色体分离
  "BUB1", "BUB1B", "MAD2L1", "CENPE", "CENPF"
)
mitochondrial_dysfunction_genes <- c(
  # 1. 线粒体电子传递链复合物亚基
  "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8", # 复合体 I (NADH脱氢酶)
  "SDHA", "SDHB", "SDHC", "SDHD", # 复合体 II (琥珀酸脱氢酶)
  "UQCRC1", "UQCRC2", "CYC1", "CYTB", # 复合体 III (bc1复合物)
  "COX1", "COX2", "COX3", "COX4I1", "COX5A", "COX6B1", # 复合体 IV (细胞色素c氧化酶)
  "ATP5A1", "ATP5B", "ATP5C1", "ATP5F1", "ATP5G1", "ATP5G2", "ATP5G3", # 复合体 V (ATP合酶)
  
  # 2. 线粒体DNA (mtDNA) 维护与复制
  "POLG", "POLG2", "TWNK", # mtDNA复制 (聚合酶γ, Twinkle解旋酶)
  "TFAM", "TFB1M", "TFB2M", # mtDNA转录
  "MGME1", "RNASEH1", "DNA2", # mtDNA加工与修复
  "TK2", "DGUOK", "SUCLA2", "SUCLG1", "RRM2B", "TYMP", # 线粒体核苷酸库维持
  
  # 3. 线粒体代谢与底物利用
  "PDHA1", "PDHB", "DLAT", "DLD", "PDP1", # 丙酮酸脱氢酶复合体
  "ACADM", "ACADS", "HADHA", "HADHB", # 脂肪酸β-氧化
  "OTC", "CPS1", "NAGS", # 尿素循环
  "GLUD1", "GOT2", "MDH2", # 氨基酸代谢与苹果酸-天冬氨酸穿梭
  
  # 4. 线粒体质量控制系统
  "PINK1", "PARK2", # PINK1-Parkin线粒体自噬
  "OPA1", "MFN1", "MFN2", "DNM1L", # 线粒体动力学 (融合/分裂)
  "HSPD1", "HSPE1", "TRAP1", "CLPX", # 线粒体分子伴侣与蛋白酶
  
  # 5. 线粒体膜转运与动力学
  "SLC25A1", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6", # 线粒体载体
  "SLC25A10", "SLC25A11", "SLC25A12", "SLC25A13",
  "SLC25A19", "SLC25A20", "SLC25A22", "SLC25A33", "SLC25A36", "SLC25A37", "SLC25A38", "SLC25A39", "SLC25A40", 
  "SLC25A42", "SLC25A44", "SLC25A46", "SLC25A47", "SLC25A48", "SLC25A51", 
  "VDAC1", "VDAC2", "VDAC3", # 电压依赖性阴离子通道
  "TIMM8A", "TIMM9", "TIMM10", "TIMM13", "TIMM17A", "TIMM22", "TIMM23", "TIMM44", "TIMM50", # 线粒体膜转位酶
  
  # 6. 氧化应激反应
  "SOD2", # 线粒体超氧化物歧化酶 (MnSOD)
  "GPX1", "GPX4", # 谷胱甘肽过氧化物酶
  "PRDX3", "PRDX5", # 过氧化物氧化还原蛋白
  "TXNRD1", "TXNRD2", "TXN2", # 硫氧还蛋白系统
  
  # 7. 线粒体疾病相关基因 (示例)
  "CMPK2", # 与线粒体核苷酸代谢和脑钙化相关[1,7](@ref)
  "PRODH", "MRPL40", "TANGO2", "ZDHHC8", "SLC25A1", "TXNRD2", "UFD1", "DGCR8" # 22q11.2缺失综合征中与线粒体功能相关的基因[3,5](@ref)
)

stem_cell_exhaustion_genes <- c(
  # 表观遗传调控因子
  "SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7", # Sirtuin家族去乙酰化酶[7,10](@ref)
  "EZH2",  # 组蛋白甲基转移酶，PRC2复合物核心组分[3](@ref)
  "BMI1",  # Polycomb组蛋白，PRC1复合物组分[3](@ref)
  
  # 细胞周期调控与衰老相关
  "CDKN1A", # p21，细胞周期蛋白依赖性激酶抑制剂[5](@ref)
  "CDKN2A", # p16INK4a，关键衰老标志物[5](@ref)
  "TP53",   # p53，细胞应激应答核心调控因子
  
  # 信号通路相关
  "WNT1", "WNT3A", "CTNNB1", # Wnt/β-catenin 信号通路[7](@ref)
  "YAP1",   # Hippo信号通路效应因子，调控组织生长和干细胞自我更新[8](@ref)
  "FOXD1",  # YAP的下游靶基因，与YAP共同构成“年轻通路”[8](@ref)
  
  # DNA损伤修复与基因组稳定性
  "WRN",    # DNA解旋酶，缺陷导致Werner综合征（早衰）[2](@ref)
  "ATM",    # DNA损伤应答核心激酶
  "BRCA1", "BRCA2", # DNA同源重组修复蛋白
  
  # 转录因子与特异性标志物
  "FOXC1", "NFATC1", # 毛囊干细胞维持的关键因子，缺失导致细胞“逃跑”[4](@ref)
  "PAPPA",  # 胎盘特异性蛋白，sirtuin缺失后异位表达驱动衰老[10](@ref)
  "PSG4",   # 妊娠特异性糖蛋白，衰老驱动因子[9](@ref)
  
  # 其他重要调控因子
  "TERF1", "TERF2", "POT1", # 端粒保护蛋白复合体成分
  "TERT", "TERC"    # 端粒酶核心组分
)

altered_intercellular_communication_genes <- c(
  # 间隙连接/直接接触通讯
  "GJA1", "GJB1", "GJB2", "GJC1", # Connexin 家族 (GJA1 即 Connexin 43)
  
  # 细胞黏附与接触
  "CDH1", "CDH5", "PECAM1", "JAM2", "ESAM",
  
  # 免疫调节与炎症因子（关键配体与受体）
  "IL6", "IL1B", "TNF", "TGFB1", "CCL2", "CXCL8", "CXCL12",
  "IL10", "IFNG", "CSF1",
  
  # 神经胶质-神经元通讯（特别在阿尔茨海默病等模型中凸显）
  "APP", "APOE", "AHNAK", "GRIN1", "GRIA1",
  
  # 细胞外囊泡 (EV) 介导的通讯
  "TSG101", "CD63", "CD9", "PDCD6IP", # EV 生物发生/标记
  "MIR9-3P", "MIR223", # EV 携带的关键 miRNA（其宿主基因或前体）
  
  # 衰老相关分泌表型 (SASP) 核心组分
  "CDKN1A", "CDKN2A", "TP53", "NFKB1", "RELA",
  
  # 神经递质与受体
  "SLC6A4", "HTR2A", "DRD2", "GRIN2B",
  
  # 血管生成与内皮通讯
  "VEGFA", "VEGFR2", "ANGPT1", "TEK",
  
  # Wnt 信号通路
  "WNT3A", "CTNNB1", "DVL1",
  
  # 其他重要调控因子
  "SIRT1", "MYC", "NOTCH1", "DLL4","THBS1"
)

dysbiosis_related_genes <- c(
  # 宿主模式识别受体
  "TLR2", "TLR4", "TLR5", "NOD1", "NOD2", "CD14",
  
  # 炎症信号通路
  "MYD88", "NFKB1", "RELA", "IL6", "IL1B", "IL10", "IL17A", "IL23", "TNF",
  "IFNG", "STAT3", "NLRP3", "NLRP6", "NLRP12", "CASP1", "IL18",
  
  # 抗菌肽
  "DEFAs", "DEFG1", "REG3G", "LYZ", "PLA2G2A",
  
  # 黏膜屏障功能
  "MUC2", "OCLN", "TJP1", "CDH1",
  
  # 胆汁酸代谢与法尼醇X受体信号
  "FXR", "SHP", "CYP7A1", "CYP8B1", "CYP7B1", "IBABP", "OSTalpha", "OSTbeta",
  
  # 代谢与运输
  "GPR43", "GPR41", "GPR109A", "PXR", "CAR",
  
  # 自噬与内质网应激
  "ATG16L1", "IRGM", "XBP1"
)
disabled_macroautophagy_genes <- c(
  # 自噬起始与成核
  "ULK1", "ULK2", "ATG13", "RB1CC1", 
  "BECN1", "AMBRA1", "PIK3C3", "PIK3R4", "UVRAG", "SH3GLB1",
  
  # 自噬体膜延伸与LC3脂化系统
  "ATG9A", "ATG9B", "WIPI1", "WIPI2",
  "ATG5", "ATG12", "ATG16L1", "ATG16L2",
  "ATG3", "ATG4A", "ATG4B", "ATG4C", "ATG4D", "ATG7", 
  "MAP1LC3A", "MAP1LC3B", "GABARAP", "GABARAPL1", "GABARAPL2",
  
  # 选择性自噬受体与衔接蛋白
  "SQSTM1", "NBR1", "OPTN", "CALCOCO2", "TAX1BP1", "BNIP3", "FUNDC1",
  
  # 溶酶体功能与融合
  "LAMP1", "LAMP2", "CTSB", "CTSD", "VAMP7", "STX17", "YKT6",
  
  # 关键调控激酶与信号通路
  "MTOR", "AKT1", "PRKAA1", "PRKAA2", "EEF2K", "MAPK14", "MAPK8", "DAPK1",
  
  # 转录调控与表观遗传因子
  "TFEB", "TFE3", "FOXO1", "FOXO3", "EZH2", "HDAC6",
  
  # 泛素化与蛋白质稳态
  "UBE2L3", "TRIMM21", "HERC2", "VCP",
  
  # 疾病相关基因
  "PINK1", "PARK2", "VPS35", "SNCA", "C9ORF72", "HTT", "TSC1", "TSC2"
)
epigenetic_genes <- c(
  # DNA甲基转移酶
  "DNMT1", "DNMT3A", "DNMT3B", "DNMT3L",
  
  # DNA去甲基化酶
  "TET1", "TET2", "TET3", "TDG",
  
  # 甲基化阅读器
  "MECP2", "MBD1", "MBD2", "MBD3", "MBD4"
)
deregulated_nutrient_sensing_genes <- c(
  # mTOR信号通路
  "MTOR", "RPTOR", "RICTOR", "MLST8", "AKT1", "AKT2", "AKT3",
  "TSC1", "TSC2", "RHEB", "EEF2K", "EIF4EBP1", "RPSSKB1", "RPS6KB2",
  
  # AMPK信号通路
  "PRKAA1", "PRKAA2", "PRKAB1", "PRKAB2", "PRKAG1", "PRKAG2", "PRKAG3",
  "STRADA", "STRADB", "CAB39", "CAB39L",
  
  # Sirtuins去乙酰化酶家族
  "SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7",
  
  # 胰岛素/IGF-1信号通路
  "INSR", "IGF1R", "IRS1", "IRS2", "IRS4", "PIK3CA", "PIK3CB", "PIK3CG",
  "PIK3R1", "PIK3R2", "PIK3R3", "PDPK1", "FOXO1", "FOXO3", "FOXO4", "FOXO6",
  
  # 氨基酸感应通路
  "GCN2", "EIF2AK4", "CASTOR1", "CASTOR2", "SESN2", "SESN3", "RRAGA", "RRAGB",
  "RRAGC", "RRAGD", "FLCN", "FNIP1", "FNIP2",
  
  # 葡萄糖感应
  "GCK", "GCKR", "SLC2A2", "SLC2A4", "SLC5A1", "SLC5A2",
  
  # 脂肪代谢感应
  "PPARA", "PPARD", "PPARG", "NR1H3", "NR1H2", "CREB1", "CRTC2",
  
  # 其他重要调控因子
  "LKB1", "SIK1", "CRTC2", "CRTC3", "CREB1", "CREB3", "CREB3L1", "CREB3L2",
  "CREB3L3", "CREB3L4", "CREB5", "ATF4", "ATF6", "XBP1", "ERN1", "EIF2A"
)
loss_of_proteostasis_genes <- c(
  # 分子伴侣与折叠相关
  "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA5", "HSPA8", "HSPA9", "HSPB1", 
  "DNAJA1", "DNAJB1", "DNAJC1", "HSPD1", "HSPE1", "HSPH1", 
  "BAG3", "BAG6", "CCT2", "CCT5", "CCT7", "PPIA", "PPIB", "PPID", 
  
  # 泛素-蛋白酶体系统
  "UBB", "UBC", "UBA1", "UBA52", "UBA80", "UBE2A", "UBE2B", "UBE2C", "UBE2D1", 
  "UBE2D2", "UBE2D3", "UBE2E1", "UBE2G1", "UBE2G2", "UBE2H", "UBE2I", 
  "UBE2J1", "UBE2J2", "UBE2K", "UBE2L3", "UBE2M", "UBE2N", "UBE2O", 
  "PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", 
  "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7", 
  "PSMB8", "PSMB9", "PSMB10", "PSMC1", "PSMC2", "PSMC3", "PSMC4", 
  "PSMC5", "PSMC6", "PSMD1", "PSMD2", "PSMD3", "PSMD4", "PSMD11", "PSMD12",
  
  # 自噬相关（巨自噬/伴侣介导自噬）
  "BECN1", "ATG3", "ATG5", "ATG7", "ATG9A", "ATG9B", "ATG10", "ATG12", 
  "ATG13", "ATG14", "ATG16L1", "ATG16L2", "MAP1LC3A", "MAP1LC3B", 
  "GABARAP", "GABARAPL1", "GABARAPL2", "ULK1", "ULK2", "ULK3", "ULK4", 
  "WIPI1", "WIPI2", "SQSTM1", "OPTN", "CALCOCO2", "NBR1", 
  "LAMP1", "LAMP2", "LAMP2A",
  
  # 内质网应激与未折叠蛋白反应
  "ATF4", "ATF6", "XBP1", "ERN1", "EIF2AK3", "EIF2A", "EIF2S1", "DDIT3",
  
  # 蛋白质聚集与清除相关
  "VCP", "TIA1", "TIAL1", "HNRNPA1", "HNRNPA2B1", "FUS", "TARDBP"
)
Aging_hallmarker <- list('telomere_attrition'=telomere_attrition_genes,
                         'genomic_instability'=genomic_instability_genes,
                         'mitochondrial_dysfunction'=mitochondrial_dysfunction_genes,
                         'stem_cell_exhaustion'=stem_cell_exhaustion_genes,
                         'altered_intercellular_communication'=altered_intercellular_communication_genes,
                         'dysbiosis_related'=dysbiosis_related_genes,
                         'disabled_macroautophagy'=disabled_macroautophagy_genes,
                         'epigenetic'=epigenetic_genes,
                         'deregulated_nutrient_sensing'=deregulated_nutrient_sensing_genes,
                         'loss_of_proteostasis'=loss_of_proteostasis_genes)

gene_set_list <- gene_set_list[-3]
Aging_hallmarker <- c(Aging_hallmarker,gene_set_list)
names(Aging_hallmarker)




save(Aging_hallmarker,file = hu("Aging_hallmarker.rda"))


