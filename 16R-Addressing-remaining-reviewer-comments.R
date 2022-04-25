#-----------------------------------------------------------------------------
# 16R-Addressing-remaining-reviewer-comments.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Address remaining reviewer comments
#-----------------------------------------------------------------------------

#----------------------------------------------------------#
# Load environments
#----------------------------------------------------------#

# Load dependencies
require(devtools)
require(doMC)
require(phyloseq)
require(microbiome)
require(vegan)
require(plyr)
require(dplyr)
require(tibble)
require(scales)
require(reshape2)
require(ggpubr)
require(ggsci)
require(ANCOMBC)
require(biomformat)
require(Rhdf5lib)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Breakdown how many fungal reads per sample and variance thereof (R2)
#----------------------------------------------------------#

# Using the same objects as those in Script 02
load("Interim_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_25Mar22.RData", verbose = T)

readSumRep200 <- rowSums(rep200Data_WGS_RNA_Matched)
readSumBacteria <- rowSums(rep200Data_WGS_RNA_Matched_Bacteria)
readSumFungi <- rowSums(rep200Data_WGS_RNA_Matched_Fungi)
# Sanity check
all(names(readSumRep200) == rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts)) # TRUE
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_rep200 <- unname(readSumRep200)
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_bacteria <- unname(readSumBacteria)
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_fungi <- unname(readSumFungi)

metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_rep200_total <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_rep200/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_total_reads
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_rep200_unmapped <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_rep200/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_unmapped_reads
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_bacteria_total <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_bacteria/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_total_reads
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_bacteria_unmapped <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_bacteria/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_unmapped
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_fungi_total <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_fungi/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_total_reads
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_fungi_unmapped <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_fungi/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_unmapped

metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$percent_rep200_total <- 100*metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_rep200_total
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$percent_rep200_unmapped <- 100*metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_rep200_unmapped
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$percent_bacteria_total <- 100*metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_bacteria_total
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$percent_bacteria_unmapped <- 100*metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_bacteria_unmapped
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$percent_fungi_total <- 100*metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_fungi_total
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$percent_fungi_unmapped <- 100*metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_fungi_unmapped

# Non-zero fungi samples
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_Nonzero <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>%
  filter(readSumFungi > 0) %>% droplevels()
rep200Data_WGS_RNA_Matched_Fungi_Nonzero <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_Nonzero),]

## Subset to WGS and RNA
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_WGS <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>%
  filter(experimental_strategy == "WGS") %>% droplevels()
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_RNA <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>%
  filter(experimental_strategy == "RNA-Seq") %>% droplevels()
rep200Data_WGS_RNA_Matched_Fungi_WGS <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_WGS),]
rep200Data_WGS_RNA_Matched_Fungi_RNA <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_RNA),]

require(Rmisc)
# Everything
summary(rowSums(rep200Data_WGS_RNA_Matched_Fungi))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0      10      96    7556    5076 2671215
CI(rowSums(rep200Data_WGS_RNA_Matched_Fungi))
# upper     mean    lower 
# 8275.395 7555.492 6835.589

# Non-zero
summary(rowSums(rep200Data_WGS_RNA_Matched_Fungi_Nonzero))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1      12     114    7780    5403 2671215
CI(rowSums(rep200Data_WGS_RNA_Matched_Fungi_Nonzero))
# upper     mean    lower 
# 8520.640 7779.674 7038.709


# WGS only
summary(rowSums(rep200Data_WGS_RNA_Matched_Fungi_WGS))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0    4744   10066   18782   17983 2671215
CI(rowSums(rep200Data_WGS_RNA_Matched_Fungi_WGS))
# upper     mean    lower 
# 20652.63 18781.54 16910.45

# RNA only
summary(rowSums(rep200Data_WGS_RNA_Matched_Fungi_RNA))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       5      23    2622     131 1718647
CI(rowSums(rep200Data_WGS_RNA_Matched_Fungi_RNA))
# upper     mean    lower 
# 3229.663 2621.698 2013.733

cols2Keep <- c("investigation","sample_type","data_submitting_center_label","experimental_strategy",
               "sum_rep200", "sum_bacteria", "sum_fungi")
bamCountWithFungi <- droplevels(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts[,cols2Keep])
bamCountWithFungi.melted <- bamCountWithFungi %>%
  rownames_to_column("sampleid") %>%
  mutate(investigation_short = gsub("TCGA-","",investigation)) %>%
  reshape2::melt(id.vars = c("sampleid","investigation_short","investigation","sample_type","data_submitting_center_label","experimental_strategy"))
bamCountWithFungi.melted$variable <- factor(bamCountWithFungi.melted$variable,
                                               levels = c("sum_rep200", "sum_bacteria", "sum_fungi"))

#-----------------------Plot all fungi read counts-----------------------#
require(rstatix)
bamCountWithFungi.melted %>%
  filter(variable %in% c("sum_fungi")) %>%
  # filter(sample_type == "Primary Tumor") %>%
  group_by(variable) %>%
  anova_test(value ~ investigation_short) %>% data.frame() -> anovaValAllSTs

filePath <- "Figures/Supplementary_Figures/"
bamCountWithFungi.melted %>%
  filter(variable %in% c("sum_fungi")) %>%
  # filter(sample_type == "Primary Tumor") %>%
  ggplot(aes(reorder(investigation_short, value, FUN=median),value, fill=investigation_short)) +
  geom_boxplot(position = "dodge") + theme_pubr() +
  ylab("Total number of reads") + xlab("TCGA Cancer Type") +
  ggtitle("Total number of reads in TCGA across all sample types mapped to fungal genomes in rep200") +
  rotate_x_text(90) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
  scale_fill_igv() + guides(fill="none") +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5)) +
  stat_summary(geom = "text", angle = 0,
               fun.data = function(x){c(y = -1, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0)) +
  annotate("text", x=4, y=1e5, label = paste0("ANOVA: F=",anovaValAllSTs$F,", p=",anovaValAllSTs$p), color = "red")
ggsave(filename = paste0(filePath,"total_num_reads_tcga_fungi_all_sample_types_AllSeqPlatforms_6Apr22.pdf"),
       dpi = "retina", units = "in", width = 12, height = 4)

#-----------------------Plot PT fungi read counts-----------------------#

require(rstatix)
bamCountWithFungi.melted %>%
  filter(variable %in% c("sum_fungi")) %>%
  filter(sample_type == "Primary Tumor") %>%
  group_by(variable) %>%
  anova_test(value ~ investigation_short) %>% data.frame() -> anovaValAllPTs

filePath <- "Figures/Supplementary_Figures/"
bamCountWithFungi.melted %>%
  filter(variable %in% c("sum_fungi")) %>%
  filter(sample_type == "Primary Tumor") %>%
  ggplot(aes(reorder(investigation_short, value, FUN=median),value, fill=investigation_short)) +
  geom_boxplot(position = "dodge") + theme_pubr() +
  ylab("Total number of reads") + xlab("TCGA Cancer Type") +
  ggtitle("Total number of reads in TCGA across primary tumors mapped to fungal genomes in rep200") +
  rotate_x_text(90) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
  scale_fill_igv() + guides(fill="none") +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5)) +
  stat_summary(geom = "text", angle = 0,
               fun.data = function(x){c(y = -1, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0)) +
  annotate("text", x=4, y=1e5, label = paste0("ANOVA: F=",anovaValAllPTs$F,", p=",anovaValAllPTs$p), color = "red")
ggsave(filename = paste0(filePath,"total_num_reads_tcga_fungi_PT_AllSeqPlatforms_6Apr22.pdf"),
       dpi = "retina", units = "in", width = 12, height = 4)

#----------------------------------------------------------#
# Rerun TCGA ML with other models to address reported concerns with boosting (R3)
#----------------------------------------------------------#

# Script S27R - random forest

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_RF <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_random_forest_6Apr22_ALL_6Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_RF$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_RF$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_RF <- mlPerfAll10k_Allcancer_RF[,!(colnames(mlPerfAll10k_Allcancer_RF) == "X")]
colnames(mlPerfAll10k_Allcancer_RF)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_Allcancer_RF$nullAUPR <- ifelse(mlPerfAll10k_Allcancer_RF$minorityClassName == "SolidTissueNormal",
                                                                       yes=mlPerfAll10k_Allcancer_RF$majorityClassSize/(mlPerfAll10k_Allcancer_RF$minorityClassSize+mlPerfAll10k_Allcancer_RF$majorityClassSize),
                                                                       no=mlPerfAll10k_Allcancer_RF$minorityClassSize/(mlPerfAll10k_Allcancer_RF$minorityClassSize+mlPerfAll10k_Allcancer_RF$majorityClassSize))
mlPerfAll10k_Allcancer_RF$nullAUROC <- 0.5

# Rename entries in the "datasetName" column

# VSNM
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "snmDataOGUFungiDecontamV2"] <- "RF: Fungi decontaminated"
# DecontamV2 per seq center
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "rep200_HiSeq_Fungi_DecontamV2_HMS"] <- "RF: HMS species decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "rep200_HiSeq_Fungi_DecontamV2_MDA"] <- "RF: MDA species decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "rep200_HiSeq_Fungi_DecontamV2_BCM"] <- "RF: BCM species decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "rep200_HiSeq_Fungi_DecontamV2_WashU"] <- "RF: WashU species decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "rep200_HiSeq_Fungi_DecontamV2_Broad_WGS"] <- "RF: Broad species decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "rep200_HiSeq_Fungi_DecontamV2_UNC"] <- "RF: UNC species decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "rep200_HiSeq_Fungi_DecontamV2_CMS"] <- "RF: CMS species decontaminated (RNA-Seq)"
## NOTE: Broad_RNA only included GBM tumors, so no ML comparisons were made

# HMS - decontam
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_phylum"] <- "RF: HMS phylum level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_class"] <- "RF: HMS class level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_order"] <- "RF: HMS order level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_family"] <- "RF: HMS family level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_genus"] <- "RF: HMS genus level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_species"] <- "RF: HMS species level decontaminated (WGS)"
# BCM - decontam
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_phylum"] <- "RF: BCM phylum level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_class"] <- "RF: BCM class level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_order"] <- "RF: BCM order level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_family"] <- "RF: BCM family level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_genus"] <- "RF: BCM genus level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_species"] <- "RF: BCM species level decontaminated (WGS)"
# MDA - decontam
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_phylum"] <- "RF: MDA phylum level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_class"] <- "RF: MDA class level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_order"] <- "RF: MDA order level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_family"] <- "RF: MDA family level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_genus"] <- "RF: MDA genus level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_species"] <- "RF: MDA species level decontaminated (WGS)"
# WashU - decontam
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_phylum"] <- "RF: WashU phylum level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_class"] <- "RF: WashU class level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_order"] <- "RF: WashU order level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_family"] <- "RF: WashU family level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_genus"] <- "RF: WashU genus level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_species"] <- "RF: WashU species level decontaminated (WGS)"
# Broad_WGS - decontam
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_phylum"] <- "RF: Broad phylum level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_class"] <- "RF: Broad class level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_order"] <- "RF: Broad order level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_family"] <- "RF: Broad family level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_genus"] <- "RF: Broad genus level decontaminated (WGS)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_species"] <- "RF: Broad species level decontaminated (WGS)"
# UNC - decontam
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_phylum"] <- "RF: UNC phylum level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_class"] <- "RF: UNC class level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_order"] <- "RF: UNC order level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_family"] <- "RF: UNC family level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_genus"] <- "RF: UNC genus level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_species"] <- "RF: UNC species level decontaminated (RNA-Seq)"
# CMS - decontam
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_phylum"] <- "RF: CMS phylum level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_class"] <- "RF: CMS class level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_order"] <- "RF: CMS order level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_family"] <- "RF: CMS family level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_genus"] <- "RF: CMS genus level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_RF$datasetName[mlPerfAll10k_Allcancer_RF$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_species"] <- "RF: CMS species level decontaminated (RNA-Seq)"
## NOTE: Broad_RNA only included GBM tumors, so no ML comparisons were made

mlPerfAll10k_Allcancer_RF$datasetName <- factor(mlPerfAll10k_Allcancer_RF$datasetName,
                                                                          levels = c(# VSNM
                                                                                     "RF: Fungi decontaminated",
                                                                                     # DecontamV2 per seq center
                                                                                     "RF: HMS species decontaminated (WGS)",
                                                                                     "RF: MDA species decontaminated (WGS)",
                                                                                     "RF: BCM species decontaminated (WGS)",
                                                                                     "RF: WashU species decontaminated (WGS)",
                                                                                     "RF: Broad species decontaminated (WGS)",
                                                                                     "RF: UNC species decontaminated (RNA-Seq)",
                                                                                     "RF: CMS species decontaminated (RNA-Seq)",
                                                                                     # HMS
                                                                                     "RF: HMS phylum level decontaminated (WGS)",
                                                                                     "RF: HMS class level decontaminated (WGS)",
                                                                                     "RF: HMS order level decontaminated (WGS)",
                                                                                     "RF: HMS family level decontaminated (WGS)",
                                                                                     "RF: HMS genus level decontaminated (WGS)",
                                                                                     "RF: HMS species level decontaminated (WGS)",
                                                                                     # BCM
                                                                                     "RF: BCM phylum level decontaminated (WGS)",
                                                                                     "RF: BCM class level decontaminated (WGS)",
                                                                                     "RF: BCM order level decontaminated (WGS)",
                                                                                     "RF: BCM family level decontaminated (WGS)",
                                                                                     "RF: BCM genus level decontaminated (WGS)",
                                                                                     "RF: BCM species level decontaminated (WGS)",
                                                                                     # MDA
                                                                                     "RF: MDA phylum level decontaminated (WGS)",
                                                                                     "RF: MDA class level decontaminated (WGS)",
                                                                                     "RF: MDA order level decontaminated (WGS)",
                                                                                     "RF: MDA family level decontaminated (WGS)",
                                                                                     "RF: MDA genus level decontaminated (WGS)",
                                                                                     "RF: MDA species level decontaminated (WGS)",
                                                                                     # WashU
                                                                                     "RF: WashU phylum level decontaminated (WGS)",
                                                                                     "RF: WashU class level decontaminated (WGS)",
                                                                                     "RF: WashU order level decontaminated (WGS)",
                                                                                     "RF: WashU family level decontaminated (WGS)",
                                                                                     "RF: WashU genus level decontaminated (WGS)",
                                                                                     "RF: WashU species level decontaminated (WGS)",
                                                                                     # Broad
                                                                                     "RF: Broad phylum level decontaminated (WGS)",
                                                                                     "RF: Broad class level decontaminated (WGS)",
                                                                                     "RF: Broad order level decontaminated (WGS)",
                                                                                     "RF: Broad family level decontaminated (WGS)",
                                                                                     "RF: Broad genus level decontaminated (WGS)",
                                                                                     "RF: Broad species level decontaminated (WGS)",
                                                                                     # UNC
                                                                                     "RF: UNC phylum level decontaminated (RNA-Seq)",
                                                                                     "RF: UNC class level decontaminated (RNA-Seq)",
                                                                                     "RF: UNC order level decontaminated (RNA-Seq)",
                                                                                     "RF: UNC family level decontaminated (RNA-Seq)",
                                                                                     "RF: UNC genus level decontaminated (RNA-Seq)",
                                                                                     "RF: UNC species level decontaminated (RNA-Seq)",
                                                                                     # CMS
                                                                                     "RF: CMS phylum level decontaminated (RNA-Seq)",
                                                                                     "RF: CMS class level decontaminated (RNA-Seq)",
                                                                                     "RF: CMS order level decontaminated (RNA-Seq)",
                                                                                     "RF: CMS family level decontaminated (RNA-Seq)",
                                                                                     "RF: CMS genus level decontaminated (RNA-Seq)",
                                                                                     "RF: CMS species level decontaminated (RNA-Seq)"))

save(mlPerfAll10k_Allcancer_RF, file = "Interim_data/mlPerfAll10k_Allcancer_RF_6Apr22.RData")

#--------------------------Summarized VSNM plot--------------------------#
load("Interim_data/mlPerfAll10k_Allcancer_2Apr22.RData", verbose = T)

mlPerfAll10k_Allcancer_RF_VSNM <- mlPerfAll10k_Allcancer_RF %>%
  filter(datasetName == "RF: Fungi decontaminated") %>% distinct() %>% 
  rename(AUROC_RF = AUROC, AUPR_RF = AUPR, rep_RF = rep) %>% 
  mutate(idx = paste(abbrev,rep_RF,gsub(" ","",sampleType),sep="_")) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  select(AUROC_RF, AUPR_RF, rep_RF, idx, abbrev, sampleType) %>% droplevels()
mlPerfAll10k_Allcancer_Species <- mlPerfAll10k_Allcancer %>%
  filter(datasetName == "Species (Decontaminated)") %>% distinct() %>% 
  rename(AUROC_GBM = AUROC, AUPR_GBM = AUPR, rep_GBM = rep) %>% 
  mutate(idx = paste(abbrev,rep_GBM,gsub(" ","",sampleType),sep="_")) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  select(AUROC_GBM, AUPR_GBM, rep_GBM, idx, abbrev, sampleType) %>% droplevels()

mlPerf_RF_GBM_VSNM_Joined <- left_join(mlPerfAll10k_Allcancer_RF_VSNM,
                                  mlPerfAll10k_Allcancer_Species,
                                  by = "idx")
require(ggpmisc)
mlPerf_RF_GBM_VSNM_Joined %>%
  rename(Investigation=abbrev.x) %>%
  ggplot(aes(x = AUROC_RF, y =  AUROC_GBM, color = Investigation)) +
  geom_point(alpha = 0.6) + coord_equal() + 
  theme_pubr() + scale_color_igv() + guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  facet_wrap(~sampleType.x) + geom_abline(linetype = 2) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  labs(x = "RF classifier per-fold performance",
       y = "GBM classifier per-fold performance",
       title = "AUROC | GBM vs. RF | Matched Per-Fold CV Performance | TCGA Batch-Corrected Data") +
  stat_cor(aes(x = AUROC_RF, y =  AUROC_GBM), method = "spearman", cor.coef.name = "rho", inherit.aes = FALSE) +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) -> roc.vsnm

mlPerf_RF_GBM_VSNM_Joined %>%
  rename(Investigation=abbrev.x) %>% 
  ggplot(aes(x = AUPR_RF, y =  AUPR_GBM, color = Investigation)) +
  geom_point(alpha = 0.6) + coord_equal() +
  theme_pubr() + scale_color_igv() + guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  facet_wrap(~sampleType.x) +
  geom_abline(linetype = 2) +
  labs(x = "RF classifier per-fold performance",
       y = "GBM classifier per-fold performance",
       title = "AUPR | GBM vs. RF | Matched Per-Fold CV Performance | TCGA Batch-Corrected Data") +
  stat_cor(aes(x = AUPR_RF, y =  AUPR_GBM), method = "spearman", cor.coef.name = "rho", inherit.aes = FALSE) +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) -> pr.vsnm

combinedPlotVSNM <- ggarrange(roc.vsnm, pr.vsnm, nrow = 2, legend = "top", common.legend = TRUE) 
ggsave(plot=combinedPlotVSNM,
       filename = "Figures/Supplementary_Figures/ml_gbm_vs_rf_vsnm_perf_combined_6Apr22.pdf",
       units = "in", width = 12, height = 10, dpi = "retina")

#--------------------------Summarized raw data plot--------------------------#
load("Interim_data/mlPerfAll10k_Allcancer_Raw_2Apr22.RData", verbose = T)

mlPerfAll10k_Allcancer_Raw_RF <- mlPerfAll10k_Allcancer_RF %>%
  filter(grepl("species decontaminated",datasetName)) %>% distinct() %>% 
  rename(AUROC_RF = AUROC, AUPR_RF = AUPR, rep_RF = rep) %>% 
  mutate(seqCenterTmp = gsub(" species decontaminated.+","",datasetName)) %>%
  mutate(seqCenter = gsub("^RF\\: ", "", seqCenterTmp)) %>%
  mutate(idx = paste(abbrev,rep_RF,gsub(" ","",sampleType),seqCenter,sep="_")) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  select(AUROC_RF, AUPR_RF, rep_RF, idx, abbrev, sampleType, seqCenter) %>% droplevels()
mlPerfAll10k_Allcancer_Raw_GBM <- mlPerfAll10k_Allcancer_Raw %>%
  filter(grepl("species decontaminated",datasetName)) %>% distinct() %>% 
  rename(AUROC_GBM = AUROC, AUPR_GBM = AUPR, rep_GBM = rep) %>% 
  mutate(seqCenter = gsub(" species decontaminated.+","",datasetName)) %>%
  mutate(idx = paste(abbrev,rep_GBM,gsub(" ","",sampleType),seqCenter,sep="_")) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  select(AUROC_GBM, AUPR_GBM, rep_GBM, idx, abbrev, sampleType, seqCenter) %>% droplevels()

mlPerf_RF_GBM_Raw_Joined <- left_join(mlPerfAll10k_Allcancer_Raw_RF,
                                      mlPerfAll10k_Allcancer_Raw_GBM, by = "idx")

mlPerf_RF_GBM_Raw_Joined %>%
  rename(Investigation=abbrev.x) %>%
  ggplot(aes(x = AUROC_RF, y =  AUROC_GBM, color = Investigation)) +
  geom_point(alpha = 0.6) + coord_equal() + 
  theme_pubr() + scale_color_igv() + guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  facet_wrap(~sampleType.x) + geom_abline(linetype = 2) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  labs(x = "RF classifier per-fold performance",
       y = "GBM classifier per-fold performance",
       title = "AUROC | GBM vs. RF | Matched Per-Fold CV Performance | TCGA Raw Decontaminated Data\n(All Sequencing Centers)") +
  stat_cor(aes(x = AUROC_RF, y =  AUROC_GBM), method = "spearman", cor.coef.name = "rho", inherit.aes = FALSE) +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) -> roc.raw

mlPerf_RF_GBM_Raw_Joined %>%
  rename(Investigation=abbrev.x) %>% 
  ggplot(aes(x = AUPR_RF, y =  AUPR_GBM, color = Investigation)) +
  geom_point(alpha = 0.6) + coord_equal() +
  theme_pubr() + scale_color_igv() + guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  facet_wrap(~sampleType.x) +
  geom_abline(linetype = 2) +
  labs(x = "RF classifier per-fold performance",
       y = "GBM classifier per-fold performance",
       title = "AUPR | GBM vs. RF | Matched Per-Fold CV Performance | TCGA Raw Decontaminated Data\n(All Sequencing Centers)") +
  stat_cor(aes(x = AUPR_RF, y =  AUPR_GBM), method = "spearman", cor.coef.name = "rho", inherit.aes = FALSE) +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) -> pr.raw

combinedPlotRaw <- ggarrange(roc.raw, pr.raw, nrow = 2, legend = "top", common.legend = TRUE) 
ggsave(plot=combinedPlotRaw,
       filename = "Figures/Supplementary_Figures/ml_gbm_vs_rf_raw_perf_combined_6Apr22.pdf",
       units = "in", width = 12, height = 10, dpi = "retina")

#--------------------------Summarized per taxa level plot--------------------------#
load("Interim_data/mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann_2Apr22.RData", verbose = T)

mlPerfAll10k_Allcancer_TaxaLevel_RF <- mlPerfAll10k_Allcancer_RF %>%
  filter(grepl("level",datasetName)) %>% distinct() %>% 
  rename(AUROC_RF = AUROC, AUPR_RF = AUPR, rep_RF = rep) %>% 
  mutate(seqCenterTmp1 = gsub("phylum|order|class|family|genus|species|[[:space:]]","",datasetName)) %>%
  mutate(seqCenterTmp2 = gsub("leveldecontaminated.+","",seqCenterTmp1)) %>%
  mutate(seqCenter = gsub("^RF\\:", "", seqCenterTmp2)) %>%
  mutate(taxaLevelTmp1 = gsub("^RF\\: ","",datasetName)) %>%
  mutate(taxaLevelTmp2 = gsub("HMS|BCM|MDA|WashU|Broad|UNC|CMS|[[:space:]]","",taxaLevelTmp1)) %>%
  mutate(taxaLevel = gsub("level.+","",taxaLevelTmp2)) %>%
  mutate(idx = paste(abbrev,rep_RF,gsub(" ","",sampleType),seqCenter,taxaLevel,sep="_")) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  select(AUROC_RF, AUPR_RF, rep_RF, idx, abbrev, sampleType, seqCenter, taxaLevel) %>% droplevels()
mlPerfAll10k_Allcancer_TaxaLevel_GBM <- mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(grepl("decontaminated",datasetName)) %>% distinct() %>% 
  rename(AUROC_GBM = AUROC, AUPR_GBM = AUPR, rep_GBM = rep) %>% 
  mutate(seqCenterTmp1 = gsub("phylum|order|class|family|genus|species|[[:space:]]","",datasetName)) %>%
  mutate(seqCenter = gsub("decontaminated.+","",seqCenterTmp1)) %>%
  mutate(taxaLevelTmp1 = gsub("HMS|BCM|MDA|WashU|Broad|UNC|CMS|[[:space:]]","",datasetName)) %>%
  mutate(taxaLevel = gsub("decontaminated.+","",taxaLevelTmp1)) %>%
  mutate(idx = paste(abbrev,rep_GBM,gsub(" ","",sampleType),seqCenter,taxaLevel,sep="_")) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  select(AUROC_GBM, AUPR_GBM, rep_GBM, idx, abbrev, sampleType, seqCenter,taxaLevel) %>% droplevels()

mlPerf_RF_GBM_TaxaLevel_Joined <- left_join(mlPerfAll10k_Allcancer_TaxaLevel_RF,
                                            mlPerfAll10k_Allcancer_TaxaLevel_GBM, by = "idx")

saveTaxaLevelPlots <- function(joinedTable,level){
  require(stringr)
  joinedTable %>%
    filter(taxaLevel.x == level) %>%
    rename(Investigation=abbrev.x) %>%
    ggplot(aes(x = AUROC_RF, y =  AUROC_GBM, color = Investigation)) +
    geom_point(alpha = 0.6) + coord_equal() + 
    theme_pubr() + scale_color_igv() + guides(color=guide_legend(nrow=3, byrow=TRUE)) +
    facet_wrap(~sampleType.x) + geom_abline(linetype = 2) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    labs(x = "RF classifier per-fold performance",
         y = "GBM classifier per-fold performance",
         title = paste0("AUROC | ",str_to_title(level)," | GBM vs. RF | Matched Per-Fold CV Performance | TCGA Raw Decontaminated Data\n(All Sequencing Centers)")) +
    stat_cor(aes(x = AUROC_RF, y =  AUROC_GBM), method = "spearman", cor.coef.name = "rho", inherit.aes = FALSE) +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) -> roc.plot
  
  joinedTable %>%
    filter(taxaLevel.x == level) %>%
    rename(Investigation=abbrev.x) %>% 
    ggplot(aes(x = AUPR_RF, y =  AUPR_GBM, color = Investigation)) +
    geom_point(alpha = 0.6) + coord_equal() +
    theme_pubr() + scale_color_igv() + guides(color=guide_legend(nrow=3, byrow=TRUE)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    facet_wrap(~sampleType.x) +
    geom_abline(linetype = 2) +
    labs(x = "RF classifier per-fold performance",
         y = "GBM classifier per-fold performance",
         title = paste0("AUPR | ",str_to_title(level)," | GBM vs. RF | Matched Per-Fold CV Performance | TCGA Raw Decontaminated Data\n(All Sequencing Centers)")) +
    stat_cor(aes(x = AUPR_RF, y =  AUPR_GBM), method = "spearman", cor.coef.name = "rho", inherit.aes = FALSE) +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) -> pr.plot
  
  combinedPlot <- ggarrange(roc.plot, pr.plot, nrow = 2, legend = "top", common.legend = TRUE) 
  print(combinedPlot)
  ggsave(plot=combinedPlot,
         filename = paste0("Figures/Supplementary_Figures/ml_gbm_vs_rf_",level,"_perf_combined_6Apr22.pdf"),
         units = "in", width = 12, height = 10, dpi = "retina")
}

saveTaxaLevelPlots(mlPerf_RF_GBM_TaxaLevel_Joined,
                   level = "phylum")
saveTaxaLevelPlots(mlPerf_RF_GBM_TaxaLevel_Joined,
                   level = "class")
saveTaxaLevelPlots(mlPerf_RF_GBM_TaxaLevel_Joined,
                   level = "order")
saveTaxaLevelPlots(mlPerf_RF_GBM_TaxaLevel_Joined,
                   level = "family")
saveTaxaLevelPlots(mlPerf_RF_GBM_TaxaLevel_Joined,
                   level = "genus")
saveTaxaLevelPlots(mlPerf_RF_GBM_TaxaLevel_Joined,
                   level = "species") # same as "Raw" plot above

#----------------------------------------------------------#
# Rerun TCGA ML with iterations of 90% train-10% test to address reported concerns with CV (R3)
#----------------------------------------------------------#

# Script S28R

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Train_Test <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_gbm_train_test_ALL_6Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_Train_Test$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_Train_Test$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_Train_Test <- mlPerfAll10k_Allcancer_Train_Test[,!(colnames(mlPerfAll10k_Allcancer_Train_Test) == "X")]
colnames(mlPerfAll10k_Allcancer_Train_Test)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_Allcancer_Train_Test$nullAUPR <- ifelse(mlPerfAll10k_Allcancer_Train_Test$minorityClassName == "SolidTissueNormal",
                                             yes=mlPerfAll10k_Allcancer_Train_Test$majorityClassSize/(mlPerfAll10k_Allcancer_Train_Test$minorityClassSize+mlPerfAll10k_Allcancer_Train_Test$majorityClassSize),
                                             no=mlPerfAll10k_Allcancer_Train_Test$minorityClassSize/(mlPerfAll10k_Allcancer_Train_Test$minorityClassSize+mlPerfAll10k_Allcancer_Train_Test$majorityClassSize))
mlPerfAll10k_Allcancer_Train_Test$nullAUROC <- 0.5

# Rename entries in the "datasetName" column

# VSNM
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "snmDataOGUFungiDecontamV2"] <- "Splits: Fungi decontaminated"
# DecontamV2 per seq center
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "rep200_HiSeq_Fungi_DecontamV2_HMS"] <- "Splits: HMS species decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "rep200_HiSeq_Fungi_DecontamV2_MDA"] <- "Splits: MDA species decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "rep200_HiSeq_Fungi_DecontamV2_BCM"] <- "Splits: BCM species decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "rep200_HiSeq_Fungi_DecontamV2_WashU"] <- "Splits: WashU species decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "rep200_HiSeq_Fungi_DecontamV2_Broad_WGS"] <- "Splits: Broad species decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "rep200_HiSeq_Fungi_DecontamV2_UNC"] <- "Splits: UNC species decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "rep200_HiSeq_Fungi_DecontamV2_CMS"] <- "Splits: CMS species decontaminated (RNA-Seq)"
## NOTE: Broad_RNA only included GBM tumors, so no ML comparisons were made

# HMS - decontam
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_phylum"] <- "Splits: HMS phylum level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_class"] <- "Splits: HMS class level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_order"] <- "Splits: HMS order level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_family"] <- "Splits: HMS family level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_genus"] <- "Splits: HMS genus level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_species"] <- "Splits: HMS species level decontaminated (WGS)"
# BCM - decontam
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_phylum"] <- "Splits: BCM phylum level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_class"] <- "Splits: BCM class level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_order"] <- "Splits: BCM order level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_family"] <- "Splits: BCM family level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_genus"] <- "Splits: BCM genus level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_species"] <- "Splits: BCM species level decontaminated (WGS)"
# MDA - decontam
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_phylum"] <- "Splits: MDA phylum level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_class"] <- "Splits: MDA class level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_order"] <- "Splits: MDA order level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_family"] <- "Splits: MDA family level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_genus"] <- "Splits: MDA genus level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_species"] <- "Splits: MDA species level decontaminated (WGS)"
# WashU - decontam
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_phylum"] <- "Splits: WashU phylum level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_class"] <- "Splits: WashU class level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_order"] <- "Splits: WashU order level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_family"] <- "Splits: WashU family level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_genus"] <- "Splits: WashU genus level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_species"] <- "Splits: WashU species level decontaminated (WGS)"
# Broad_WGS - decontam
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_phylum"] <- "Splits: Broad phylum level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_class"] <- "Splits: Broad class level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_order"] <- "Splits: Broad order level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_family"] <- "Splits: Broad family level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_genus"] <- "Splits: Broad genus level decontaminated (WGS)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_species"] <- "Splits: Broad species level decontaminated (WGS)"
# UNC - decontam
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_phylum"] <- "Splits: UNC phylum level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_class"] <- "Splits: UNC class level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_order"] <- "Splits: UNC order level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_family"] <- "Splits: UNC family level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_genus"] <- "Splits: UNC genus level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_species"] <- "Splits: UNC species level decontaminated (RNA-Seq)"
# CMS - decontam
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_phylum"] <- "Splits: CMS phylum level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_class"] <- "Splits: CMS class level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_order"] <- "Splits: CMS order level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_family"] <- "Splits: CMS family level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_genus"] <- "Splits: CMS genus level decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Train_Test$datasetName[mlPerfAll10k_Allcancer_Train_Test$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_species"] <- "Splits: CMS species level decontaminated (RNA-Seq)"
## NOTE: Broad_RNA only included GBM tumors, so no ML comparisons were made

mlPerfAll10k_Allcancer_Train_Test$datasetName <- factor(mlPerfAll10k_Allcancer_Train_Test$datasetName,
                                                levels = c(# VSNM
                                                  "Splits: Fungi decontaminated",
                                                  # DecontamV2 per seq center
                                                  "Splits: HMS species decontaminated (WGS)",
                                                  "Splits: MDA species decontaminated (WGS)",
                                                  "Splits: BCM species decontaminated (WGS)",
                                                  "Splits: WashU species decontaminated (WGS)",
                                                  "Splits: Broad species decontaminated (WGS)",
                                                  "Splits: UNC species decontaminated (RNA-Seq)",
                                                  "Splits: CMS species decontaminated (RNA-Seq)",
                                                  # HMS
                                                  "Splits: HMS phylum level decontaminated (WGS)",
                                                  "Splits: HMS class level decontaminated (WGS)",
                                                  "Splits: HMS order level decontaminated (WGS)",
                                                  "Splits: HMS family level decontaminated (WGS)",
                                                  "Splits: HMS genus level decontaminated (WGS)",
                                                  "Splits: HMS species level decontaminated (WGS)",
                                                  # BCM
                                                  "Splits: BCM phylum level decontaminated (WGS)",
                                                  "Splits: BCM class level decontaminated (WGS)",
                                                  "Splits: BCM order level decontaminated (WGS)",
                                                  "Splits: BCM family level decontaminated (WGS)",
                                                  "Splits: BCM genus level decontaminated (WGS)",
                                                  "Splits: BCM species level decontaminated (WGS)",
                                                  # MDA
                                                  "Splits: MDA phylum level decontaminated (WGS)",
                                                  "Splits: MDA class level decontaminated (WGS)",
                                                  "Splits: MDA order level decontaminated (WGS)",
                                                  "Splits: MDA family level decontaminated (WGS)",
                                                  "Splits: MDA genus level decontaminated (WGS)",
                                                  "Splits: MDA species level decontaminated (WGS)",
                                                  # WashU
                                                  "Splits: WashU phylum level decontaminated (WGS)",
                                                  "Splits: WashU class level decontaminated (WGS)",
                                                  "Splits: WashU order level decontaminated (WGS)",
                                                  "Splits: WashU family level decontaminated (WGS)",
                                                  "Splits: WashU genus level decontaminated (WGS)",
                                                  "Splits: WashU species level decontaminated (WGS)",
                                                  # Broad
                                                  "Splits: Broad phylum level decontaminated (WGS)",
                                                  "Splits: Broad class level decontaminated (WGS)",
                                                  "Splits: Broad order level decontaminated (WGS)",
                                                  "Splits: Broad family level decontaminated (WGS)",
                                                  "Splits: Broad genus level decontaminated (WGS)",
                                                  "Splits: Broad species level decontaminated (WGS)",
                                                  # UNC
                                                  "Splits: UNC phylum level decontaminated (RNA-Seq)",
                                                  "Splits: UNC class level decontaminated (RNA-Seq)",
                                                  "Splits: UNC order level decontaminated (RNA-Seq)",
                                                  "Splits: UNC family level decontaminated (RNA-Seq)",
                                                  "Splits: UNC genus level decontaminated (RNA-Seq)",
                                                  "Splits: UNC species level decontaminated (RNA-Seq)",
                                                  # CMS
                                                  "Splits: CMS phylum level decontaminated (RNA-Seq)",
                                                  "Splits: CMS class level decontaminated (RNA-Seq)",
                                                  "Splits: CMS order level decontaminated (RNA-Seq)",
                                                  "Splits: CMS family level decontaminated (RNA-Seq)",
                                                  "Splits: CMS genus level decontaminated (RNA-Seq)",
                                                  "Splits: CMS species level decontaminated (RNA-Seq)"))

save(mlPerfAll10k_Allcancer_Train_Test, file = "Interim_data/mlPerfAll10k_Allcancer_Train_Test_6Apr22.RData")

#--------------------------Summarized VSNM plot--------------------------#
load("Interim_data/mlPerfAll10k_Allcancer_2Apr22.RData", verbose = T)

mlPerfAll10k_Allcancer_Train_Test_VSNM <- mlPerfAll10k_Allcancer_Train_Test %>%
  filter(datasetName == "Splits: Fungi decontaminated") %>% distinct() %>% 
  rename(AUROC_Splits = AUROC, AUPR_Splits = AUPR, rep_Splits = rep) %>% 
  mutate(idx = paste(abbrev,rep_Splits,gsub(" ","",sampleType),sep="_")) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  select(AUROC_Splits, AUPR_Splits, rep_Splits, idx, abbrev, sampleType) %>% droplevels()
mlPerfAll10k_Allcancer_Species <- mlPerfAll10k_Allcancer %>%
  filter(datasetName == "Species (Decontaminated)") %>% distinct() %>% 
  rename(AUROC_CV = AUROC, AUPR_CV = AUPR, rep_CV = rep) %>% 
  mutate(idx = paste(abbrev,rep_CV,gsub(" ","",sampleType),sep="_")) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  select(AUROC_CV, AUPR_CV, rep_CV, idx, abbrev, sampleType) %>% droplevels()

#-----------------------Calculate 95% CIs-----------------------#
mlPerfAll10k_Allcancer_Train_Test_VSNM_CIs <- mlPerfAll10k_Allcancer_Train_Test_VSNM %>%
  reshape2::melt(id.vars = c("rep_Splits","abbrev","sampleType","idx")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","abbrev","sampleType")) %>% 
  mutate(idx = paste(gsub("_|Splits|CV","",variable),abbrev,gsub(" ","",sampleType),sep="_")) %>% 
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  rename(variable_Splits=variable, value_Splits = value, sd_Splits=sd, se_Splits=se, ci_Splits=ci, N_Splits=N)

mlPerfAll10k_Allcancer_Species_CIs <- mlPerfAll10k_Allcancer_Species %>%
  reshape2::melt(id.vars = c("rep_CV","abbrev","sampleType","idx")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","abbrev","sampleType")) %>% 
  mutate(idx = paste(gsub("_|Splits|CV","",variable),abbrev,gsub(" ","",sampleType),sep="_")) %>% 
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  rename(variable_CV=variable, value_CV=value, sd_CV=sd, se_CV=se, ci_CV=ci, N_CV=N)

mlPerf_Train_Test_VSNM_Joined_CIs <- left_join(mlPerfAll10k_Allcancer_Train_Test_VSNM_CIs,
                                           mlPerfAll10k_Allcancer_Species_CIs,
                                           by = "idx")

#-----------------------Plot using 95% CIs-----------------------#

mlPerf_Train_Test_VSNM_Joined_CIs %>%
  rename(Investigation=abbrev.x) %>%
  filter(grepl("AUROC",variable_Splits)) %>%
  filter(grepl("AUROC",variable_CV)) %>%
  ggplot(aes(x=value_Splits, y=value_CV, color=Investigation)) +
  facet_wrap(~sampleType.x) + geom_abline(linetype = 2) +
  geom_point(alpha=0.6) + 
  geom_errorbarh(aes(xmin=ifelse(value_Splits-ci_Splits<0,0,value_Splits-ci_Splits),
                     xmax=ifelse(value_Splits+ci_Splits>1,1,value_Splits+ci_Splits)), height=1) +
  geom_errorbar(aes(ymin=ifelse(value_CV-ci_CV<0,0,value_CV-ci_CV),
                    ymax=ifelse(value_CV+ci_CV>1,1,value_CV+ci_CV)), width=1) +
  coord_equal() + theme_pubr() +
  scale_color_igv() + guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  labs(x = "Train-test classifier performance (10 iterations)",
       y = "CV classifier performance (10 folds)",
       title = "AUROC | TCGA Batch-Corrected Data") +
  stat_cor(aes(x = value_Splits, y =  value_CV), method = "spearman", cor.coef.name = "rho", inherit.aes = FALSE) +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) -> roc.splits.vsnm

mlPerf_Train_Test_VSNM_Joined_CIs %>%
  rename(Investigation=abbrev.x) %>%
  filter(grepl("AUPR",variable_Splits)) %>%
  filter(grepl("AUPR",variable_CV)) %>%
  ggplot(aes(x=value_Splits, y=value_CV, color=Investigation)) +
  facet_wrap(~sampleType.x) + geom_abline(linetype = 2) +
  geom_point(alpha=0.6) + 
  geom_errorbarh(aes(xmin=ifelse(value_Splits-ci_Splits<0,0,value_Splits-ci_Splits),
                     xmax=ifelse(value_Splits+ci_Splits>1,1,value_Splits+ci_Splits)), height=1) +
  geom_errorbar(aes(ymin=ifelse(value_CV-ci_CV<0,0,value_CV-ci_CV),
                    ymax=ifelse(value_CV+ci_CV>1,1,value_CV+ci_CV)), width=1) +
  coord_equal() + theme_pubr() +
  scale_color_igv() + guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  labs(x = "Train-test classifier performance (10 iterations)",
       y = "CV classifier performance (10 folds)",
       title = "AUPR | TCGA Batch-Corrected Data") +
  stat_cor(aes(x = value_Splits, y =  value_CV), method = "spearman", cor.coef.name = "rho", inherit.aes = FALSE) +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) -> pr.splits.vsnm

combinedPlotSplitsVSNM <- ggarrange(roc.splits.vsnm, pr.splits.vsnm, nrow = 2, legend = "top", common.legend = TRUE) 
ggsave(plot=combinedPlotSplitsVSNM,
       filename = "Figures/Supplementary_Figures/ml_splits_vs_cv_vsnm_perf_combined_8Apr22.pdf",
       units = "in", width = 12, height = 10, dpi = "retina")

#--------------------------Summarized raw data plot--------------------------#
load("Interim_data/mlPerfAll10k_Allcancer_Raw_2Apr22.RData", verbose = T)

mlPerfAll10k_Allcancer_Train_Test_Raw <- mlPerfAll10k_Allcancer_Train_Test %>%
  filter(grepl("species decontaminated",datasetName)) %>% distinct() %>% 
  rename(AUROC_Splits = AUROC, AUPR_Splits = AUPR, rep_Splits = rep) %>% 
  mutate(seqCenterTmp = gsub(" species decontaminated.+","",datasetName)) %>%
  mutate(seqCenter = gsub("^Splits\\: ", "", seqCenterTmp)) %>%
  mutate(idx = paste(abbrev,rep_Splits,gsub(" ","",sampleType),seqCenter,sep="_")) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  select(AUROC_Splits, AUPR_Splits, rep_Splits, idx, abbrev, sampleType, seqCenter) %>% droplevels()
mlPerfAll10k_Allcancer_CV_Raw <- mlPerfAll10k_Allcancer_Raw %>%
  filter(grepl("species decontaminated",datasetName)) %>% distinct() %>% 
  rename(AUROC_CV = AUROC, AUPR_CV = AUPR, rep_CV = rep) %>% 
  mutate(seqCenter = gsub(" species decontaminated.+","",datasetName)) %>%
  mutate(idx = paste(abbrev,rep_CV,gsub(" ","",sampleType),seqCenter,sep="_")) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  select(AUROC_CV, AUPR_CV, rep_CV, idx, abbrev, sampleType, seqCenter) %>% droplevels()

#-----------------------Calculate 95% CIs-----------------------#
mlPerfAll10k_Allcancer_Train_Test_Raw_CIs <- mlPerfAll10k_Allcancer_Train_Test_Raw %>%
  reshape2::melt(id.vars = c("rep_Splits","abbrev","sampleType","idx","seqCenter")) %>% 
  summarySE(measurevar = "value", groupvars = c("variable","abbrev","sampleType","seqCenter")) %>%
  mutate(idx = paste(gsub("_|Splits|CV","",variable),abbrev,gsub(" ","",sampleType),seqCenter,sep="_")) %>% 
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  rename(variable_Splits=variable, value_Splits = value, sd_Splits=sd, se_Splits=se, ci_Splits=ci, N_Splits=N)

mlPerfAll10k_Allcancer_CV_Raw_CIs <- mlPerfAll10k_Allcancer_CV_Raw %>%
  reshape2::melt(id.vars = c("rep_CV","abbrev","sampleType","idx","seqCenter")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","abbrev","sampleType","seqCenter")) %>% 
  mutate(idx = paste(gsub("_|Splits|CV","",variable),abbrev,gsub(" ","",sampleType),seqCenter,sep="_")) %>% 
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  rename(variable_CV=variable, value_CV=value, sd_CV=sd, se_CV=se, ci_CV=ci, N_CV=N)

mlPerf_Train_Test_Raw_Joined_CIs <- left_join(mlPerfAll10k_Allcancer_Train_Test_Raw_CIs,
                                              mlPerfAll10k_Allcancer_CV_Raw_CIs,
                                               by = "idx")

#-----------------------Plot using 95% CIs-----------------------#

mlPerf_Train_Test_Raw_Joined_CIs %>%
  rename(Investigation=abbrev.x) %>%
  filter(grepl("AUROC",variable_Splits)) %>%
  filter(grepl("AUROC",variable_CV)) %>%
  ggplot(aes(x=value_Splits, y=value_CV, color=Investigation)) +
  facet_wrap(~sampleType.x) + geom_abline(linetype = 2) +
  geom_point(alpha=0.6) + 
  geom_errorbarh(aes(xmin=ifelse(value_Splits-ci_Splits<0,0,value_Splits-ci_Splits),
                     xmax=ifelse(value_Splits+ci_Splits>1,1,value_Splits+ci_Splits)), height=1) +
  geom_errorbar(aes(ymin=ifelse(value_CV-ci_CV<0,0,value_CV-ci_CV),
                    ymax=ifelse(value_CV+ci_CV>1,1,value_CV+ci_CV)), width=1) +
  coord_equal() + theme_pubr() +
  scale_color_igv() + guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  labs(x = "Train-test classifier performance (10 iterations)",
       y = "CV classifier performance (10 folds)",
       title = "AUROC | TCGA Raw Decontaminated Data | All Sequencing Centers") +
  stat_cor(aes(x = value_Splits, y =  value_CV), method = "spearman", cor.coef.name = "rho", inherit.aes = FALSE) +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) -> roc.splits.raw

mlPerf_Train_Test_Raw_Joined_CIs %>%
  rename(Investigation=abbrev.x) %>%
  filter(grepl("AUPR",variable_Splits)) %>%
  filter(grepl("AUPR",variable_CV)) %>%
  ggplot(aes(x=value_Splits, y=value_CV, color=Investigation)) +
  facet_wrap(~sampleType.x) + geom_abline(linetype = 2) +
  geom_point(alpha=0.6) + 
  geom_errorbarh(aes(xmin=ifelse(value_Splits-ci_Splits<0,0,value_Splits-ci_Splits),
                     xmax=ifelse(value_Splits+ci_Splits>1,1,value_Splits+ci_Splits)), height=1) +
  geom_errorbar(aes(ymin=ifelse(value_CV-ci_CV<0,0,value_CV-ci_CV),
                    ymax=ifelse(value_CV+ci_CV>1,1,value_CV+ci_CV)), width=1) +
  coord_equal() + theme_pubr() +
  scale_color_igv() + guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  labs(x = "Train-test classifier performance (10 iterations)",
       y = "CV classifier performance (10 folds)",
       title = "AUPR | TCGA Raw Decontaminated Data | All Sequencing Centers") +
  stat_cor(aes(x = value_Splits, y =  value_CV), method = "spearman", cor.coef.name = "rho", inherit.aes = FALSE) +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) -> pr.splits.raw

combinedPlotSplitsRaw <- ggarrange(roc.splits.raw, pr.splits.raw, nrow = 2, legend = "top", common.legend = TRUE) 
ggsave(plot=combinedPlotSplitsRaw,
       filename = "Figures/Supplementary_Figures/ml_splits_vs_cv_raw_perf_combined_8Apr22.pdf",
       units = "in", width = 12, height = 10, dpi = "retina")

#--------------------------Summarized per taxa level plot--------------------------#
load("Interim_data/mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann_2Apr22.RData", verbose = T)

mlPerfAll10k_Allcancer_Train_Test_TaxaLevel <- mlPerfAll10k_Allcancer_Train_Test %>%
  filter(grepl("level",datasetName)) %>% distinct() %>% 
  rename(AUROC_Splits = AUROC, AUPR_Splits = AUPR, rep_Splits = rep) %>% 
  mutate(seqCenterTmp1 = gsub("phylum|order|class|family|genus|species|[[:space:]]","",datasetName)) %>%
  mutate(seqCenterTmp2 = gsub("leveldecontaminated.+","",seqCenterTmp1)) %>%
  mutate(seqCenter = gsub("^Splits\\:", "", seqCenterTmp2)) %>%
  mutate(taxaLevelTmp1 = gsub("^Splits\\: ","",datasetName)) %>%
  mutate(taxaLevelTmp2 = gsub("HMS|BCM|MDA|WashU|Broad|UNC|CMS|[[:space:]]","",taxaLevelTmp1)) %>%
  mutate(taxaLevel = gsub("level.+","",taxaLevelTmp2)) %>%
  mutate(idx = paste(abbrev,rep_Splits,gsub(" ","",sampleType),seqCenter,taxaLevel,sep="_")) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  select(AUROC_Splits, AUPR_Splits, rep_Splits, idx, abbrev, sampleType, seqCenter, taxaLevel) %>% droplevels()
mlPerfAll10k_Allcancer_CV_TaxaLevel <- mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(grepl("decontaminated",datasetName)) %>% distinct() %>% 
  rename(AUROC_CV = AUROC, AUPR_CV = AUPR, rep_CV = rep) %>% 
  mutate(seqCenterTmp1 = gsub("phylum|order|class|family|genus|species|[[:space:]]","",datasetName)) %>%
  mutate(seqCenter = gsub("decontaminated.+","",seqCenterTmp1)) %>%
  mutate(taxaLevelTmp1 = gsub("HMS|BCM|MDA|WashU|Broad|UNC|CMS|[[:space:]]","",datasetName)) %>%
  mutate(taxaLevel = gsub("decontaminated.+","",taxaLevelTmp1)) %>%
  mutate(idx = paste(abbrev,rep_CV,gsub(" ","",sampleType),seqCenter,taxaLevel,sep="_")) %>%
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  select(AUROC_CV, AUPR_CV, rep_CV, idx, abbrev, sampleType, seqCenter,taxaLevel) %>% droplevels()

#-----------------------Calculate 95% CIs-----------------------#
mlPerfAll10k_Allcancer_Train_Test_TaxaLevel_CIs <- mlPerfAll10k_Allcancer_Train_Test_TaxaLevel %>%
  reshape2::melt(id.vars = c("rep_Splits","abbrev","sampleType","idx","seqCenter","taxaLevel")) %>% 
  summarySE(measurevar = "value", groupvars = c("variable","abbrev","sampleType","seqCenter","taxaLevel")) %>%
  mutate(idx = paste(gsub("_|Splits|CV","",variable),abbrev,gsub(" ","",sampleType),seqCenter,taxaLevel,sep="_")) %>% 
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  rename(variable_Splits=variable, value_Splits = value, sd_Splits=sd, se_Splits=se, ci_Splits=ci, N_Splits=N)

mlPerfAll10k_Allcancer_CV_TaxaLevel_CIs <- mlPerfAll10k_Allcancer_CV_TaxaLevel %>%
  reshape2::melt(id.vars = c("rep_CV","abbrev","sampleType","idx","seqCenter","taxaLevel")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","abbrev","sampleType","seqCenter","taxaLevel")) %>% 
  mutate(idx = paste(gsub("_|Splits|CV","",variable),abbrev,gsub(" ","",sampleType),seqCenter,taxaLevel,sep="_")) %>% 
  mutate(sampleType = factor(sampleType, levels = c("Primary Tumor", "Blood Derived Normal", "Primary Tumor vs Solid Tissue Normal"))) %>%
  rename(variable_CV=variable, value_CV=value, sd_CV=sd, se_CV=se, ci_CV=ci, N_CV=N)

mlPerf_Train_Test_TaxaLevel_Joined_CIs <- left_join(mlPerfAll10k_Allcancer_Train_Test_TaxaLevel_CIs,
                                              mlPerfAll10k_Allcancer_CV_TaxaLevel_CIs,
                                              by = "idx")

#-----------------------Plot using 95% CIs-----------------------#

saveTaxaLevelPlotsTrainTest <- function(joinedTable,level){
  
  require(stringr)
  joinedTable %>%
    filter(taxaLevel.x == level) %>%
    filter(grepl("AUROC",variable_Splits)) %>%
    filter(grepl("AUROC",variable_CV)) %>%
    rename(Investigation=abbrev.x) %>%
    ggplot(aes(x=value_Splits, y=value_CV, color=Investigation)) +
    facet_wrap(~sampleType.x) + geom_abline(linetype = 2) +
    geom_point(alpha=0.6) + 
    geom_errorbarh(aes(xmin=ifelse(value_Splits-ci_Splits<0,0,value_Splits-ci_Splits),
                       xmax=ifelse(value_Splits+ci_Splits>1,1,value_Splits+ci_Splits)), height=1) +
    geom_errorbar(aes(ymin=ifelse(value_CV-ci_CV<0,0,value_CV-ci_CV),
                      ymax=ifelse(value_CV+ci_CV>1,1,value_CV+ci_CV)), width=1) +
    coord_equal() + theme_pubr() +
    scale_color_igv() + guides(color=guide_legend(nrow=3, byrow=TRUE)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    labs(x = "Train-test classifier performance (10 iterations)",
         y = "CV classifier performance (10 folds)",
         title = paste0("AUROC | ",str_to_title(level)," Raw Decontaminated Data | All TCGA Sequencing Centers")) +
    stat_cor(aes(x = value_Splits, y =  value_CV), method = "spearman", cor.coef.name = "rho", inherit.aes = FALSE) +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) -> roc.plot
  
  joinedTable %>%
    filter(taxaLevel.x == level) %>%
    filter(grepl("AUPR",variable_Splits)) %>%
    filter(grepl("AUPR",variable_CV)) %>%
    rename(Investigation=abbrev.x) %>% 
    ggplot(aes(x=value_Splits, y=value_CV, color=Investigation)) +
    facet_wrap(~sampleType.x) + geom_abline(linetype = 2) +
    geom_point(alpha=0.6) + 
    geom_errorbarh(aes(xmin=ifelse(value_Splits-ci_Splits<0,0,value_Splits-ci_Splits),
                       xmax=ifelse(value_Splits+ci_Splits>1,1,value_Splits+ci_Splits)), height=1) +
    geom_errorbar(aes(ymin=ifelse(value_CV-ci_CV<0,0,value_CV-ci_CV),
                      ymax=ifelse(value_CV+ci_CV>1,1,value_CV+ci_CV)), width=1) +
    coord_equal() + theme_pubr() +
    scale_color_igv() + guides(color=guide_legend(nrow=3, byrow=TRUE)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    labs(x = "Train-test classifier performance (10 iterations)",
         y = "CV classifier performance (10 folds)",
         title = paste0("AUPR | ",str_to_title(level)," Raw Decontaminated Data | All TCGA Sequencing Centers")) +
    stat_cor(aes(x = value_Splits, y =  value_CV), method = "spearman", cor.coef.name = "rho", inherit.aes = FALSE) +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) -> pr.plot
  
  combinedPlot <- ggarrange(roc.plot, pr.plot, nrow = 2, legend = "top", common.legend = TRUE) 
  print(combinedPlot)
  ggsave(plot=combinedPlot,
         filename = paste0("Figures/Supplementary_Figures/ml_splits_vs_cv_",level,"_perf_combined_8Apr22.pdf"),
         units = "in", width = 12, height = 10, dpi = "retina")
}

saveTaxaLevelPlotsTrainTest(mlPerf_Train_Test_TaxaLevel_Joined_CIs,
                   level = "phylum")
saveTaxaLevelPlotsTrainTest(mlPerf_Train_Test_TaxaLevel_Joined_CIs,
                   level = "class")
saveTaxaLevelPlotsTrainTest(mlPerf_Train_Test_TaxaLevel_Joined_CIs,
                   level = "order")
saveTaxaLevelPlotsTrainTest(mlPerf_Train_Test_TaxaLevel_Joined_CIs,
                   level = "family")
saveTaxaLevelPlotsTrainTest(mlPerf_Train_Test_TaxaLevel_Joined_CIs,
                   level = "genus")
saveTaxaLevelPlotsTrainTest(mlPerf_Train_Test_TaxaLevel_Joined_CIs,
                   level = "species") # same as "Raw" plot above

tmp <- mlPerf_Train_Test_TaxaLevel_Joined_CIs %>%
  filter(taxaLevel.x == "genus", variable_Splits == "AUROC_Splits") %>% droplevels()
