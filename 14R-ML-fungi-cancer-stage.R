#-----------------------------------------------------------------------------
# 14R-ML-fungi-cancer-stage.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Test whether fungi can predict stage I vs stage IV cancers
# - Test whether blood-derived fungi still distinguish cancer types in stages I-II
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
require(reshape2)
require(ggpubr)
require(ggsci)
require(ANCOMBC)
require(biomformat)
require(Rhdf5lib)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Load TCGA fungi data
#----------------------------------------------------------#

load("Interim_data/data_for_pvca_tcga_taxa_levels_decontamV2_2Apr22.RData", verbose = T)

#----------------------------------------------------------#
# Add stage labels
#----------------------------------------------------------#
# Remove unclear or non-useful stages
metaQiitaCombined_Nonzero_DecontamV2_Path <- droplevels(metaQiitaCombined_Nonzero_DecontamV2[! (metaQiitaCombined_Nonzero_DecontamV2$pathologic_stage_label == "Not available" |
                                                                                                  metaQiitaCombined_Nonzero_DecontamV2$pathologic_stage_label == "I or II NOS" |
                                                                                                  metaQiitaCombined_Nonzero_DecontamV2$pathologic_stage_label == "Stage 0" |
                                                                                                  metaQiitaCombined_Nonzero_DecontamV2$pathologic_stage_label == "Stage IS" |
                                                                                                  metaQiitaCombined_Nonzero_DecontamV2$pathologic_stage_label == "Stage Tis" |
                                                                                                  metaQiitaCombined_Nonzero_DecontamV2$pathologic_stage_label == "Stage X"),])

tumorStageVector <- factor(metaQiitaCombined_Nonzero_DecontamV2_Path$pathologic_stage_label)
levels(tumorStageVector) <- list(StageI = c("Stage I", "Stage IA", "Stage IB", "Stage IC"),
                                 StageII = c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC"),
                                 StageIII = c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),
                                 StageIV = c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC"))
metaQiitaCombined_Nonzero_DecontamV2_Path$pathologic_stage_label_binned <- tumorStageVector
table(metaQiitaCombined_Nonzero_DecontamV2_Path$pathologic_stage_label_binned)

#---------------------------Subset VSNM data to primary tumors---------------------------#

metaQiitaCombined_Nonzero_DecontamV2_Path_PT <- metaQiitaCombined_Nonzero_DecontamV2_Path %>% 
  filter(sample_type == "Primary Tumor") %>% droplevels()

rep200FungiDecontamV2SpeciesVSNM <- rep200FungiDecontamV2SpeciesVSNM_Obj$snmData
rep200FungiDecontamV2SpeciesVSNM_Path_PT <- rep200FungiDecontamV2SpeciesVSNM[rownames(metaQiitaCombined_Nonzero_DecontamV2_Path_PT),]

save(metaQiitaCombined_Nonzero_DecontamV2_Path_PT,
     rep200FungiDecontamV2SpeciesVSNM_Path_PT,
     file = "Interim_data/data_for_ml_stage_5Apr22.RData")

# See scripts: S23R

#---------------------------Subset VSNM data blood derived normals---------------------------#

metaQiitaCombined_Nonzero_DecontamV2_Path_BDN <- metaQiitaCombined_Nonzero_DecontamV2_Path %>% 
  filter(sample_type == "Blood Derived Normal") %>%
  filter(pathologic_stage_label_binned %in% c("StageI","StageII")) %>% droplevels()
rep200FungiDecontamV2SpeciesVSNM_Path_BDN <- rep200FungiDecontamV2SpeciesVSNM[rownames(metaQiitaCombined_Nonzero_DecontamV2_Path_BDN),]

save(metaQiitaCombined_Nonzero_DecontamV2_Path_BDN,
     rep200FungiDecontamV2SpeciesVSNM_Path_BDN,
     file = "Interim_data/data_for_ml_vsnm_bdn_stageI_II_5Apr22.RData")

# See scripts: S24R

#---------------------------Subset raw data per seq center to blood derived normals---------------------------#

load("Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy_with_coverage_filter_decontamV2_2Apr22.RData")

metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_HMS <- metaQiitaCombined_Nonzero_DecontamV2_Path_BDN %>%
  filter(data_submitting_center_label == "Harvard Medical School") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_BCM <- metaQiitaCombined_Nonzero_DecontamV2_Path_BDN %>%
  filter(data_submitting_center_label == "Baylor College of Medicine") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_MDA <- metaQiitaCombined_Nonzero_DecontamV2_Path_BDN %>%
  filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_WashU <- metaQiitaCombined_Nonzero_DecontamV2_Path_BDN %>%
  filter(data_submitting_center_label == "Washington University School of Medicine") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_Broad_WGS <- metaQiitaCombined_Nonzero_DecontamV2_Path_BDN %>%
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>% droplevels()

rep200_HiSeq_Fungi_DecontamV2_Path_BDN_HMS <- rep200_HiSeq_Fungi_DecontamV2_HMS[rownames(metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_HMS),]
rep200_HiSeq_Fungi_DecontamV2_Path_BDN_BCM <- rep200_HiSeq_Fungi_DecontamV2_BCM[rownames(metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_BCM),]
rep200_HiSeq_Fungi_DecontamV2_Path_BDN_MDA <- rep200_HiSeq_Fungi_DecontamV2_MDA[rownames(metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_MDA),]
rep200_HiSeq_Fungi_DecontamV2_Path_BDN_WashU <- rep200_HiSeq_Fungi_DecontamV2_WashU[rownames(metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_WashU),]
rep200_HiSeq_Fungi_DecontamV2_Path_BDN_Broad_WGS <- rep200_HiSeq_Fungi_DecontamV2_Broad_WGS[rownames(metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_Broad_WGS),]

save(rep200_HiSeq_Fungi_DecontamV2_Path_BDN_HMS,
     rep200_HiSeq_Fungi_DecontamV2_Path_BDN_BCM,
     rep200_HiSeq_Fungi_DecontamV2_Path_BDN_MDA,
     rep200_HiSeq_Fungi_DecontamV2_Path_BDN_WashU,
     rep200_HiSeq_Fungi_DecontamV2_Path_BDN_Broad_WGS,
     
     metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_HMS,
     metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_BCM,
     metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_MDA,
     metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_WashU,
     metaQiitaCombined_Nonzero_DecontamV2_Path_BDN_Broad_WGS,
     file = "Interim_data/data_for_ml_raw_counts_bdn_stageI_II_5Apr22.RData")

# See scripts: S25R

#----------------------------------------------------------#
# Plot results for primary tumors
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_stage <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_stage_ALL_DecontamV2_5Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_stage$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_stage$diseaseType,"abbrev"]
mlPerfAll10k_stage <- mlPerfAll10k_stage[,!(colnames(mlPerfAll10k_stage) == "X")]
colnames(mlPerfAll10k_stage)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For stage, "StageIV" is the positive class and is used to calculate null AUPR
mlPerfAll10k_stage$nullAUPR <- ifelse(mlPerfAll10k_stage$minorityClassName == "StageIV",
                                      yes=mlPerfAll10k_stage$minorityClassSize/(mlPerfAll10k_stage$minorityClassSize+mlPerfAll10k_stage$majorityClassSize),
                                      no=mlPerfAll10k_stage$majorityClassSize/(mlPerfAll10k_stage$minorityClassSize+mlPerfAll10k_stage$majorityClassSize))
mlPerfAll10k_stage$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
table(mlPerfAll10k_stage$datasetName)
mlPerfAll10k_stage$datasetName[mlPerfAll10k_stage$datasetName == "rep200FungiDecontamV2SpeciesVSNM_Path_PT"] <- "VSNM species decontaminated"

mlPerfAll10k_stage$datasetName <- factor(mlPerfAll10k_stage$datasetName,
                                         levels = c("VSNM species decontaminated"))
#-------------------------Plot performance-------------------------#
mlPerfAll10k_stage %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Stage I vs Stage IV | Intratumoral | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(0) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_tcga_stage_decontamV2.pdf", dpi = "retina",
       width = 6, height = 4, units = "in")

#----------------------------------------------------------#
# Plot results for blood derived normals - VSNM data
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_BDN_EarlyStage_VSNM <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_vsnm_bdn_early_stage_ALL_5Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_BDN_EarlyStage_VSNM$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_BDN_EarlyStage_VSNM$diseaseType,"abbrev"]
mlPerfAll10k_BDN_EarlyStage_VSNM <- mlPerfAll10k_BDN_EarlyStage_VSNM[,!(colnames(mlPerfAll10k_BDN_EarlyStage_VSNM) == "X")]
colnames(mlPerfAll10k_BDN_EarlyStage_VSNM)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For stage, "StageIV" is the positive class and is used to calculate null AUPR
mlPerfAll10k_BDN_EarlyStage_VSNM$nullAUPR <- ifelse(mlPerfAll10k_BDN_EarlyStage_VSNM$minorityClassName == "SolidTissueNormal",
                                                    yes=mlPerfAll10k_BDN_EarlyStage_VSNM$majorityClassSize/(mlPerfAll10k_BDN_EarlyStage_VSNM$minorityClassSize+mlPerfAll10k_BDN_EarlyStage_VSNM$majorityClassSize),
                                                    no=mlPerfAll10k_BDN_EarlyStage_VSNM$minorityClassSize/(mlPerfAll10k_BDN_EarlyStage_VSNM$minorityClassSize+mlPerfAll10k_BDN_EarlyStage_VSNM$majorityClassSize))
mlPerfAll10k_BDN_EarlyStage_VSNM$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
table(mlPerfAll10k_BDN_EarlyStage_VSNM$datasetName)
mlPerfAll10k_BDN_EarlyStage_VSNM$datasetName[mlPerfAll10k_BDN_EarlyStage_VSNM$datasetName == "rep200FungiDecontamV2SpeciesVSNM_Path_BDN"] <- "VSNM species decontaminated"

mlPerfAll10k_BDN_EarlyStage_VSNM$datasetName <- factor(mlPerfAll10k_BDN_EarlyStage_VSNM$datasetName,
                                                       levels = c("VSNM species decontaminated"))
# Species level
mlPerfAll10k_BDN_EarlyStage_VSNM %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","minorityClassSize","majorityClassSize")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("TCGA Stage Ia-IIc only | VSNM Data | Blood Derived Normal | 1 Vs All | Fungi") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_VSNM_BDN_early_stage_fungi_decontamV2.svg", dpi = "retina",
       width = 8, height = 4, units = "in")

require(gmodels)
mlPerfAll10k_BDN_EarlyStage_VSNM %>%
  distinct() %>% droplevels() %>%
  pull(AUROC) %>% ci()
# Estimate    CI lower    CI upper  Std. Error 
# 0.912230419 0.895167944 0.929292893 0.008616977

mlPerfAll10k_BDN_EarlyStage_VSNM %>%
  distinct() %>% droplevels() %>%
  filter(abbrev == "BRCA") %>%
  pull(AUROC) %>% ci()
# Estimate    CI lower    CI upper  Std. Error 
# 0.995297271 0.992019708 0.998574834 0.001448866

#----------------------------------------------------------#
# Plot results for blood derived normals - Raw data
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_BDN_EarlyStage_Raw <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_bdn_early_stage_by_seq_center_ALL_DecontamV2_5Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_BDN_EarlyStage_Raw$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_BDN_EarlyStage_Raw$diseaseType,"abbrev"]
mlPerfAll10k_BDN_EarlyStage_Raw <- mlPerfAll10k_BDN_EarlyStage_Raw[,!(colnames(mlPerfAll10k_BDN_EarlyStage_Raw) == "X")]
colnames(mlPerfAll10k_BDN_EarlyStage_Raw)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For stage, "StageIV" is the positive class and is used to calculate null AUPR
mlPerfAll10k_BDN_EarlyStage_Raw$nullAUPR <- ifelse(mlPerfAll10k_BDN_EarlyStage_Raw$minorityClassName == "SolidTissueNormal",
                                                   yes=mlPerfAll10k_BDN_EarlyStage_Raw$majorityClassSize/(mlPerfAll10k_BDN_EarlyStage_Raw$minorityClassSize+mlPerfAll10k_BDN_EarlyStage_Raw$majorityClassSize),
                                                   no=mlPerfAll10k_BDN_EarlyStage_Raw$minorityClassSize/(mlPerfAll10k_BDN_EarlyStage_Raw$minorityClassSize+mlPerfAll10k_BDN_EarlyStage_Raw$majorityClassSize))
mlPerfAll10k_BDN_EarlyStage_Raw$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
table(mlPerfAll10k_BDN_EarlyStage_Raw$datasetName)
mlPerfAll10k_BDN_EarlyStage_Raw$datasetName[mlPerfAll10k_BDN_EarlyStage_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_Path_BDN_HMS"] <- "HMS species decontaminated (WGS)"
mlPerfAll10k_BDN_EarlyStage_Raw$datasetName[mlPerfAll10k_BDN_EarlyStage_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_Path_BDN_BCM"] <- "BCM species decontaminated (WGS)"
mlPerfAll10k_BDN_EarlyStage_Raw$datasetName[mlPerfAll10k_BDN_EarlyStage_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_Path_BDN_MDA"] <- "MDA species decontaminated (WGS)"
mlPerfAll10k_BDN_EarlyStage_Raw$datasetName[mlPerfAll10k_BDN_EarlyStage_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_Path_BDN_Broad_WGS"] <- "Broad species decontaminated (WGS)"

mlPerfAll10k_BDN_EarlyStage_Raw$datasetName <- factor(mlPerfAll10k_BDN_EarlyStage_Raw$datasetName,
                                                      levels = c("HMS species decontaminated (WGS)",
                                                                 "BCM species decontaminated (WGS)",
                                                                 "MDA species decontaminated (WGS)",
                                                                 "Broad species decontaminated (WGS)"))
# Species level
mlPerfAll10k_BDN_EarlyStage_Raw %>%
  # filter(grepl("MDA",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","minorityClassSize","majorityClassSize")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("TCGA Stage Ia-IIc only | Raw Data | Blood Derived Normal | 1 Vs All | Fungi") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Raw_BDN_early_stage_fungi_decontamV2.svg", dpi = "retina",
       width = 8, height = 4, units = "in")

require(gmodels)
mlPerfAll10k_BDN_EarlyStage_Raw %>%
  distinct() %>% droplevels() %>%
  pull(AUROC) %>% ci()
# Estimate   CI lower   CI upper Std. Error 
# 0.88372927 0.86255872 0.90489981 0.01069166

