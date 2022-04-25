#-----------------------------------------------------------------------------
# 10R-Control-validation-analyses-TCGA.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Split data into 2 stratified halves, normalize using VSNM, train ML models, and test on each other
# - Scramble metadata and shuffle samples and re-run ML as controls
# - Compare performance of ML controls to actual data for raw data (per seq center) and VSNM (all TCGA)
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

load("Interim_data/snmDataFungi_DecontamV2_25Mar22.RData", verbose = T)
rep200TaxSplit <- read.csv("Supporting_data/rep200_lineage_map_split.csv", stringsAsFactors = FALSE, row.names = 1)

#----------------------------------------------------------#
# Split into 2 stratified halves
#----------------------------------------------------------#
require(data.table)
require(splitstackshape)

splitVars <- c("data_submitting_center_label","sample_type","disease_type")
set.seed(42)
stratSamplingHalvesMetadata <- stratified(metaQiitaCombined_Nonzero_DecontamV2, 
                                      group = splitVars, 
                                      size = 0.5,
                                      keep.rownames = TRUE, 
                                      bothSets = TRUE, 
                                      replace = FALSE)
split1 <- as.data.frame(stratSamplingHalvesMetadata[[1]])
rownames(split1) <- split1$rn
split2 <- as.data.frame(stratSamplingHalvesMetadata[[2]])
rownames(split2) <- split2$rn

centerCompare <- data.frame(Split = c(rep("Split 1", dim(split1)[1]),
                                         rep("Split 2", dim(split2)[1])),
                               SeqCenter = c(as.character(split1$data_submitting_center_label),
                                             as.character(split2$data_submitting_center_label)),
                               SampleType = c(as.character(split1$sample_type),
                                              as.character(split2$sample_type)),
                               DiseaseType = c(as.character(split1$investigation),
                                               as.character(split2$investigation)),
                               PathStage = c(as.character(split1$pathologic_stage_label),
                                             as.character(split2$pathologic_stage_label)))
ggplot(centerCompare, aes(x = SeqCenter, fill = Split)) + geom_bar(position = "dodge",stat = "count") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) + scale_fill_nejm() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5, color="black", position = position_dodge(0.9)) + 
  ylim(c(0,floor(1.1*max(table(centerCompare$SeqCenter))/2))) + theme_bw() + rotate_x_text(15) +
  labs(x = "Sequencing Center", y = "Count", title = "Validation split - Sequence center distribution")
ggsave(filename = "Figures/Other_Figures/tcga_validation_splits_seq_center.pdf", dpi = "retina", width = 14, units = "in")

ggplot(centerCompare, aes(x = SampleType, fill = Split)) + geom_bar(position = "dodge",stat = "count") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) + scale_fill_nejm() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5, color="black", position = position_dodge(0.9)) + 
  ylim(c(0,floor(1.1*max(table(centerCompare$SampleType))/2))) + theme_bw() + rotate_x_text(0) +
  labs(x = "Sample Type", y = "Count", title = "Validation split - Sample type distribution")
ggsave(filename = "Figures/Other_Figures/tcga_validation_splits_sample_type.pdf", dpi = "retina", width = 12, units = "in")

ggplot(centerCompare, aes(x = DiseaseType, fill = Split)) + geom_bar(position = "dodge",stat = "count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_nejm() +
  geom_text(stat='count', aes(label=..count..), hjust=-0.25, color="black", position = position_dodge(0.9), angle=90) + 
  ylim(c(0,floor(1.15*max(table(centerCompare$DiseaseType))/2))) + theme_bw() + rotate_x_text(30) +
  labs(x = "Sample Type", y = "Count", title = "Validation split - Disease type distribution")
ggsave(filename = "Figures/Other_Figures/tcga_validation_splits_dz_type.pdf", dpi = "retina", width = 12, units = "in")

ggplot(centerCompare, aes(x = PathStage, fill = Split)) + geom_bar(position = "dodge",stat = "count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_nejm() +
  geom_text(stat='count', aes(label=..count..), hjust=-0.25, color="black", position = position_dodge(0.9), angle=90) + 
  ylim(c(0,floor(1.2*max(table(centerCompare$PathStage))/2))) + theme_bw() + rotate_x_text(30) +
  labs(x = "Sample Type", y = "Count", title = "Validation split - Pathologic Stage distribution")
ggsave(filename = "Figures/Other_Figures/tcga_validation_splits_stage.pdf", dpi = "retina", width = 8, units = "in")

# Subset vb data using splits

valSplit1Metadata <- split1
valSplit2Metadata <- split2

valSplit1FungiCountDecontamV2Data <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(valSplit1Metadata),]
valSplit2FungiCountDecontamV2Data <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(valSplit2Metadata),]

#----------------------------------------------------------#
# VSNM batch correct each split and full dataset
#----------------------------------------------------------#

source("00-Functions.R") # for vsnmFunctionTCGA() function
## TCGA full data - biological variable: sample_type | technical variables: data_submitting_center_label + experimental_strategy
# Split 1
valSplit1FungiCountDecontamV2Data_VSNM_Obj <- vsnmFunctionTCGA(qcData = valSplit1FungiCountDecontamV2Data, qcMetadata = valSplit1Metadata) # converges
valSplit1FungiCountDecontamV2Data_VSNM <- valSplit1FungiCountDecontamV2Data_VSNM_Obj$snmData
# Split 2
valSplit2FungiCountDecontamV2Data_VSNM_Obj <- vsnmFunctionTCGA(qcData = valSplit2FungiCountDecontamV2Data, qcMetadata = valSplit2Metadata) # converges
valSplit2FungiCountDecontamV2Data_VSNM <- valSplit2FungiCountDecontamV2Data_VSNM_Obj$snmData
# Full data
valFullFungiCountDecontamV2Data_VSNM_Obj <- vsnmFunctionTCGA(qcData = rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero, qcMetadata = metaQiitaCombined_Nonzero_DecontamV2) # converges
valFullFungiCountDecontamV2Data_VSNM <- valFullFungiCountDecontamV2Data_VSNM_Obj$snmData

save(valSplit1Metadata,
     valSplit1FungiCountDecontamV2Data,
     valSplit1FungiCountDecontamV2Data_VSNM,
     
     valSplit2Metadata,
     valSplit2FungiCountDecontamV2Data,
     valSplit2FungiCountDecontamV2Data_VSNM,
     
     metaQiitaCombined_Nonzero_DecontamV2,
     rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero,
     valFullFungiCountDecontamV2Data_VSNM,
     file = "Interim_data/validation_2halves_tcga_vsnm_data_4Apr22.RData")

# Scripts: S18R, S19R

#----------------------------------------------------------#
# Plot split performances
#----------------------------------------------------------#

source("00-Functions.R") # to load plotSplitPerfs() function

## Load data (note that the files each load a "perf" object, so it needs renaming)
rm(perf)
load("Interim_data/perf_raw_data_k10_4Apr22.RData", verbose = TRUE) # saved as perf object
perf_raw_data_k10 <- perf
rm(perf)
load("Interim_data/perf_vsnm_data_k10_4Apr22.RData", verbose = TRUE)
perf_val_halves_k10 <- perf
rm(perf)

## Plot
plotSplitPerfs(perfDf = perf_raw_data_k10, baseName = "val_halves_raw_data_k10")
plotSplitPerfs(perfDf = perf_val_halves_k10, baseName = "val_halves_vsnm_data_k10")

#----------------------------------------------------------#
# Scramble labels and data: raw data
#----------------------------------------------------------#

load("Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy_with_coverage_filter_decontamV2_2Apr22.RData")
#-----------------------Scramble metadata-----------------------#
# HMS
metaQiitaCombined_Nonzero_DecontamV2_HMS_scrambled <- metaQiitaCombined_Nonzero_DecontamV2_HMS
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_HMS_scrambled$disease_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_HMS_scrambled$disease_type)
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_HMS_scrambled$sample_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_HMS_scrambled$sample_type)
# BCM
metaQiitaCombined_Nonzero_DecontamV2_BCM_scrambled <- metaQiitaCombined_Nonzero_DecontamV2_BCM
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_BCM_scrambled$disease_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_BCM_scrambled$disease_type)
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_BCM_scrambled$sample_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_BCM_scrambled$sample_type)
# MDA
metaQiitaCombined_Nonzero_DecontamV2_MDA_scrambled <- metaQiitaCombined_Nonzero_DecontamV2_MDA
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_MDA_scrambled$disease_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_MDA_scrambled$disease_type)
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_MDA_scrambled$sample_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_MDA_scrambled$sample_type)
# WashU
metaQiitaCombined_Nonzero_DecontamV2_WashU_scrambled <- metaQiitaCombined_Nonzero_DecontamV2_WashU
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_WashU_scrambled$disease_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_WashU_scrambled$disease_type)
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_WashU_scrambled$sample_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_WashU_scrambled$sample_type)
# Broad_WGS
metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_scrambled <- metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_scrambled$disease_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_scrambled$disease_type)
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_scrambled$sample_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_scrambled$sample_type)
# UNC
metaQiitaCombined_Nonzero_DecontamV2_UNC_scrambled <- metaQiitaCombined_Nonzero_DecontamV2_UNC
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_UNC_scrambled$disease_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_UNC_scrambled$disease_type)
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_UNC_scrambled$sample_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_UNC_scrambled$sample_type)
# CMS
metaQiitaCombined_Nonzero_DecontamV2_CMS_scrambled <- metaQiitaCombined_Nonzero_DecontamV2_CMS
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_CMS_scrambled$disease_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_CMS_scrambled$disease_type)
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_CMS_scrambled$sample_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_CMS_scrambled$sample_type)
# Broad_RNA
metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_scrambled <- metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_scrambled$disease_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_scrambled$disease_type)
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_scrambled$sample_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_scrambled$sample_type)

#-----------------------Shuffle raw data-----------------------#

# HMS
set.seed(42)
rep200_HiSeq_Fungi_DecontamV2_HMS_shuffled <- rep200_HiSeq_Fungi_DecontamV2_HMS[sample(nrow(rep200_HiSeq_Fungi_DecontamV2_HMS)),]
rownames(rep200_HiSeq_Fungi_DecontamV2_HMS_shuffled) <- rownames(rep200_HiSeq_Fungi_DecontamV2_HMS)
# BCM
set.seed(42)
rep200_HiSeq_Fungi_DecontamV2_BCM_shuffled <- rep200_HiSeq_Fungi_DecontamV2_BCM[sample(nrow(rep200_HiSeq_Fungi_DecontamV2_BCM)),]
rownames(rep200_HiSeq_Fungi_DecontamV2_BCM_shuffled) <- rownames(rep200_HiSeq_Fungi_DecontamV2_BCM)
# MDA
set.seed(42)
rep200_HiSeq_Fungi_DecontamV2_MDA_shuffled <- rep200_HiSeq_Fungi_DecontamV2_MDA[sample(nrow(rep200_HiSeq_Fungi_DecontamV2_MDA)),]
rownames(rep200_HiSeq_Fungi_DecontamV2_MDA_shuffled) <- rownames(rep200_HiSeq_Fungi_DecontamV2_MDA)
# WashU
set.seed(42)
rep200_HiSeq_Fungi_DecontamV2_WashU_shuffled <- rep200_HiSeq_Fungi_DecontamV2_WashU[sample(nrow(rep200_HiSeq_Fungi_DecontamV2_WashU)),]
rownames(rep200_HiSeq_Fungi_DecontamV2_WashU_shuffled) <- rownames(rep200_HiSeq_Fungi_DecontamV2_WashU)
# Broad_WGS
set.seed(42)
rep200_HiSeq_Fungi_DecontamV2_Broad_WGS_shuffled <- rep200_HiSeq_Fungi_DecontamV2_Broad_WGS[sample(nrow(rep200_HiSeq_Fungi_DecontamV2_Broad_WGS)),]
rownames(rep200_HiSeq_Fungi_DecontamV2_Broad_WGS_shuffled) <- rownames(rep200_HiSeq_Fungi_DecontamV2_Broad_WGS)
# UNC
set.seed(42)
rep200_HiSeq_Fungi_DecontamV2_UNC_shuffled <- rep200_HiSeq_Fungi_DecontamV2_UNC[sample(nrow(rep200_HiSeq_Fungi_DecontamV2_UNC)),]
rownames(rep200_HiSeq_Fungi_DecontamV2_UNC_shuffled) <- rownames(rep200_HiSeq_Fungi_DecontamV2_UNC)
# CMS
set.seed(42)
rep200_HiSeq_Fungi_DecontamV2_CMS_shuffled <- rep200_HiSeq_Fungi_DecontamV2_CMS[sample(nrow(rep200_HiSeq_Fungi_DecontamV2_CMS)),]
rownames(rep200_HiSeq_Fungi_DecontamV2_CMS_shuffled) <- rownames(rep200_HiSeq_Fungi_DecontamV2_CMS)
# Broad_RNA
set.seed(42)
rep200_HiSeq_Fungi_DecontamV2_Broad_RNA_shuffled <- rep200_HiSeq_Fungi_DecontamV2_Broad_RNA[sample(nrow(rep200_HiSeq_Fungi_DecontamV2_Broad_RNA)),]
rownames(rep200_HiSeq_Fungi_DecontamV2_Broad_RNA_shuffled) <- rownames(rep200_HiSeq_Fungi_DecontamV2_Broad_RNA)

save(metaQiitaCombined_Nonzero_DecontamV2_HMS_scrambled,
     metaQiitaCombined_Nonzero_DecontamV2_BCM_scrambled,
     metaQiitaCombined_Nonzero_DecontamV2_MDA_scrambled,
     metaQiitaCombined_Nonzero_DecontamV2_WashU_scrambled,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_scrambled,
     metaQiitaCombined_Nonzero_DecontamV2_UNC_scrambled,
     metaQiitaCombined_Nonzero_DecontamV2_CMS_scrambled,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_scrambled,
     # Scrambled meta normal data
     rep200_HiSeq_Fungi_DecontamV2_HMS,
     rep200_HiSeq_Fungi_DecontamV2_BCM,
     rep200_HiSeq_Fungi_DecontamV2_MDA,
     rep200_HiSeq_Fungi_DecontamV2_WashU,
     rep200_HiSeq_Fungi_DecontamV2_Broad_WGS,
     rep200_HiSeq_Fungi_DecontamV2_UNC,
     rep200_HiSeq_Fungi_DecontamV2_CMS,
     rep200_HiSeq_Fungi_DecontamV2_Broad_RNA,
     # Shuffled data normal metadata
     metaQiitaCombined_Nonzero_DecontamV2_HMS,
     metaQiitaCombined_Nonzero_DecontamV2_BCM,
     metaQiitaCombined_Nonzero_DecontamV2_MDA,
     metaQiitaCombined_Nonzero_DecontamV2_WashU,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS,
     metaQiitaCombined_Nonzero_DecontamV2_UNC,
     metaQiitaCombined_Nonzero_DecontamV2_CMS,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA,

     rep200_HiSeq_Fungi_DecontamV2_HMS_shuffled,
     rep200_HiSeq_Fungi_DecontamV2_BCM_shuffled,
     rep200_HiSeq_Fungi_DecontamV2_MDA_shuffled,
     rep200_HiSeq_Fungi_DecontamV2_WashU_shuffled,
     rep200_HiSeq_Fungi_DecontamV2_Broad_WGS_shuffled,
     rep200_HiSeq_Fungi_DecontamV2_UNC_shuffled,
     rep200_HiSeq_Fungi_DecontamV2_CMS_shuffled,
     rep200_HiSeq_Fungi_DecontamV2_Broad_RNA_shuffled,
     file = "Interim_data/shuffled_controls_raw_data_TCGA_4Apr22.RData")
#----------------------------------------------------------#
# Scramble labels and data: VSNM data
#----------------------------------------------------------#

load("Interim_data/snmDataFungi_DecontamV2_25Mar22.RData")
load("Interim_data/validation_2halves_tcga_vsnm_data_4Apr22.RData")

# Scramble full metadata
metaQiitaCombined_Nonzero_DecontamV2_scrambled <- metaQiitaCombined_Nonzero_DecontamV2
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_scrambled$disease_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_scrambled$disease_type)
set.seed(42)
metaQiitaCombined_Nonzero_DecontamV2_scrambled$sample_type <- sample(metaQiitaCombined_Nonzero_DecontamV2_scrambled$sample_type)

# Shuffle VSNM data
set.seed(42)
valFullFungiCountDecontamV2Data_VSNM_Shuffled <-  valFullFungiCountDecontamV2Data_VSNM[sample(nrow(valFullFungiCountDecontamV2Data_VSNM)),]
rownames(valFullFungiCountDecontamV2Data_VSNM_Shuffled) <- rownames(valFullFungiCountDecontamV2Data_VSNM)

save(metaQiitaCombined_Nonzero_DecontamV2_scrambled,
     valFullFungiCountDecontamV2Data_VSNM,
     metaQiitaCombined_Nonzero_DecontamV2,
     valFullFungiCountDecontamV2Data_VSNM_Shuffled,
     file = "Interim_data/shuffled_controls_VSNM_data_TCGA_4Apr22.RData")

#----------------------------------------------------------#
# Run ML
#----------------------------------------------------------#

# Script: S20R

#----------------------------------------------------------#
# Plot scrambled ML perf
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_scrambled_and_shuffled_controls_ALL_DecontamV2_4Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM <- mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM[,!(colnames(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM) == "X")]
colnames(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$nullAUPR <- ifelse(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$minorityClassName == "SolidTissueNormal",
                                                             yes=mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$majorityClassSize/(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$minorityClassSize+mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$majorityClassSize),
                                                             no=mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$minorityClassSize/(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$minorityClassSize+mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$majorityClassSize))
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
table(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName)
table(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName)
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_HMS_scrambled"] <- "HMS metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_HMS_shuffled"] <- "HMS count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_MDA_scrambled"] <- "MDA metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_MDA_shuffled"] <- "MDA count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_BCM_scrambled"] <- "BCM metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_BCM_shuffled"] <- "BCM count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_WashU_scrambled"] <- "WashU metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_WashU_shuffled"] <- "WashU count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_scrambled"] <- "Broad metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_shuffled"] <- "Broad count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_UNC_scrambled"] <- "UNC metadata labels scrambled (RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_UNC_shuffled"] <- "UNC count data shuffled (RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_CMS_scrambled"] <- "CMS metadata labels scrambled (RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_CMS_shuffled"] <- "CMS count data shuffled (RNA-Seq)"
## NOTE: Broad_RNA only included GBM tumors, so no ML comparisons were made
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_scrambled"] <- "VSNM all metadata labels scrambled (WGS+RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_shuffled"] <- "VSNM all count data shuffled (WGS+RNA-Seq)"

mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName <- factor(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName,
                                                                levels = c("VSNM all metadata labels scrambled (WGS+RNA-Seq)",
                                                                           "VSNM all count data shuffled (WGS+RNA-Seq)",
                                                                           "HMS metadata labels scrambled (WGS)",
                                                                           "HMS count data shuffled (WGS)",
                                                                           "MDA metadata labels scrambled (WGS)",
                                                                           "MDA count data shuffled (WGS)",
                                                                           "BCM metadata labels scrambled (WGS)",
                                                                           "BCM count data shuffled (WGS)",
                                                                           "WashU metadata labels scrambled (WGS)",
                                                                           "WashU count data shuffled (WGS)",
                                                                           "Broad metadata labels scrambled (WGS)",
                                                                           "Broad count data shuffled (WGS)",
                                                                           "UNC metadata labels scrambled (RNA-Seq)",
                                                                           "UNC count data shuffled (RNA-Seq)",
                                                                           "CMS metadata labels scrambled (RNA-Seq)",
                                                                           "CMS count data shuffled (RNA-Seq)"))
table(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName)
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
# HMS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_HMS_PT_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")
# MDA
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("MDA",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Primary Tumor | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_MDA_PT_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")
# BCM
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("BCM",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_BCM_PT_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")
# WashU
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("WashU",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Washington University (WGS) | Raw Data | Primary Tumor | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_WashU_PT_species_decontamV2.pdf", dpi = "retina",
         width = 5, height = 4, units = "in")
# Broad WGS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Broad",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_Broad_WGS_PT_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")
# UNC
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("UNC",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University of North Carolina (RNA-Seq) | Raw Data | Primary Tumor | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_UNC_PT_species_decontamV2.pdf", dpi = "retina",
         width = 8, height = 4, units = "in")
# CMS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("CMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data | Primary Tumor | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_CMS_PT_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")
# Full dataset
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("VSNM",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Full dataset (WGS+RNA-Seq) | Voom-SNM | Primary Tumor | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_TCGA_VSNM_PT_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#
# HMS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_HMS_PT_vs_NAT_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")
# MDA does not have sufficient PT and NAT samples to compare
# BCM
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_BCM_PT_vs_NAT_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")
# WashU does not have sufficient PT and NAT samples to compare
# Broad WGS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")
# UNC
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("UNC",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University of North Carolina (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_UNC_PT_vs_NAT_species_decontamV2.pdf", dpi = "retina",
         width = 8, height = 4, units = "in")
# CMS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("CMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_CMS_PT_vs_NAT_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")

# Full dataset
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("VSNM",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Full dataset (WGS+RNA-Seq) | Voom-SNM | Primary Tumor vs Solid Tissue Normal | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_TCGA_VSNM_PT_vs_NAT_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")

#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
# HMS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Blood Derived Normal | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_HMS_BDN_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")
# MDA
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("MDA",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Blood Derived Normal | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_MDA_BDN_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")
# BCM
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Blood Derived Normal | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_BCM_BDN_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")
# WashU
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("WashU",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Washington University (WGS) | Raw Data | Blood Derived Normal | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_WashU_BDN_species_decontamV2.pdf", dpi = "retina",
         width = 5, height = 4, units = "in")
# Broad WGS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Blood Derived Normal | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_Broad_WGS_BDN_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")
# Full dataset
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("VSNM",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Full dataset (WGS+RNA-Seq) | Voom-SNM | Blood Derived Normal | 1 Vs All\nScrambled & Shuffled Controls") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_TCGA_VSNM_BDN_species_decontamV2.pdf", dpi = "retina",
         width = 6, height = 4, units = "in")

#----------------------------------------------------------#
# Overlay actual and (scrambled) control performance per seq center
#----------------------------------------------------------#

load("Interim_data/mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann_2Apr22.RData") # mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann
load("Interim_data/mlPerfAll10k_Allcancer_Raw_2Apr22.RData") # mlPerfAll10k_Allcancer_Raw

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Raw_Control_Overlay <- rbind(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM,
                                                    mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann,
                                                    mlPerfAll10k_Allcancer_Raw)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt <- mlPerfAll10k_Allcancer_Raw_Control_Overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt$datasetName <- factor(mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt$datasetName,
                                                              levels = rev(levels(mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt$datasetName)))

#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
# HMS - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_HMS_PT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# BCM - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("BCM",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_BCM_PT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# MDA - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("MDA",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_MDA_PT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# WashU - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("WashU",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_WashU_PT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# Broad_WGS - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Broad",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_Broad_WGS_PT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# UNC - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("UNC",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_UNC_PT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# CMS - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("CMS",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_CMS_PT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#
# HMS - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_HMS_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# BCM - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_BCM_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
## NOTE: Neither MDA nor WashU had sufficient samples to plot primary tumor vs. NAT performance

# Broad_WGS - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# UNC - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("UNC",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_UNC_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# CMS - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("CMS",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_CMS_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")

#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
# HMS - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_HMS_BDN_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# BCM - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_BCM_BDN_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# MDA - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("MDA",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_MDA_BDN_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# WashU - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("WashU",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_WashU_BDN_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# Broad_WGS - overlay
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/control_v2_scrambled_mlPerfAll10k_rep1_Broad_WGS_BDN_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")

#----------------------------------------------------------#
# Boxplot summaries of overlay per seq center
#----------------------------------------------------------#

mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short <- mlPerfAll10k_Allcancer_Raw_Control_Overlay_Filt %>%
  filter(grepl("species|shuffled|scrambled",datasetName)) %>% droplevels()
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- gsub(" \\(WGS\\)| \\(RNA-Seq\\)| \\(WGS\\+RNA-Seq\\)","",mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetName)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- gsub(" species","",mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- gsub("high coverage","high cov",mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- gsub(" count data","",mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- gsub(" metadata labels","",mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- gsub("decontaminated","decontam",mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- gsub("all scrambled","labels scrambled",mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- gsub("all shuffled","samples shuffled",mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort)
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort <- factor(mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort,
                                                                            levels = c("HMS decontam", "HMS high cov", "HMS ∩ WIS", "HMS scrambled", "HMS shuffled",
                                                                                       "BCM decontam", "BCM high cov", "BCM ∩ WIS", "BCM scrambled", "BCM shuffled",
                                                                                       "MDA decontam", "MDA high cov", "MDA ∩ WIS", "MDA scrambled", "MDA shuffled",
                                                                                       "WashU decontam", "WashU high cov", "WashU ∩ WIS", "WashU scrambled", "WashU shuffled",
                                                                                       "Broad decontam", "Broad high cov", "Broad ∩ WIS", "Broad scrambled", "Broad shuffled",
                                                                                       "UNC decontam", "UNC high cov", "UNC ∩ WIS", "UNC scrambled", "UNC shuffled",
                                                                                       "CMS decontam", "CMS high cov", "CMS ∩ WIS", "CMS scrambled", "CMS shuffled",
                                                                                       "VSNM labels scrambled", "VSNM samples shuffled"))
table(mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short$datasetNameShort)

#----------------------------------Plot overlays----------------------------------#

## Write plot function
source("00-Functions.R") # for plotControlsRaw() function
#-------------------------Plot primary tumor overlays-------------------------#
plotControlsRaw(seqCenter = "HMS", inputSampleType = "Primary Tumor", qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "BCM", inputSampleType = "Primary Tumor", qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "MDA", inputSampleType = "Primary Tumor", qvalSize = 2.5, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "WashU", inputSampleType = "Primary Tumor", qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "Broad", inputSampleType = "Primary Tumor", qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "UNC", inputSampleType = "Primary Tumor", qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "CMS", inputSampleType = "Primary Tumor", qvalSize = 3, qvalAsterisks = TRUE)
#-------------------------Plot primary tumor vs NAT overlays-------------------------#
plotControlsRaw(seqCenter = "HMS", inputSampleType = "Primary Tumor vs Solid Tissue Normal",
                statSpacingPR = 1.8, qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "BCM", inputSampleType = "Primary Tumor vs Solid Tissue Normal",
                qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "Broad", inputSampleType = "Primary Tumor vs Solid Tissue Normal",
                statSpacingPR = 1, qvalSize = 3, qvalAsterisks = TRUE)
## NOTE: Neither MDA nor WashU had sufficient samples to plot primary tumor vs. NAT performance
plotControlsRaw(seqCenter = "UNC", inputSampleType = "Primary Tumor vs Solid Tissue Normal",
                statSpacingPR = 6, statSpacingROC = 0.9, qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "CMS", inputSampleType = "Primary Tumor vs Solid Tissue Normal",
                statSpacingPR = 1.8, qvalSize = 3, qvalAsterisks = TRUE)
#-------------------------Plot blood derived normal overlays-------------------------#
plotControlsRaw(seqCenter = "HMS", inputSampleType = "Blood Derived Normal",
                qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "BCM", inputSampleType = "Blood Derived Normal", 
                statSpacingROC=1, qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "MDA", inputSampleType = "Blood Derived Normal", 
                qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "WashU", inputSampleType = "Blood Derived Normal",
                statSpacingPR = 0.8, statSpacingROC = 0.8, qvalSize = 3, qvalAsterisks = TRUE)
plotControlsRaw(seqCenter = "Broad", inputSampleType = "Blood Derived Normal",
                qvalSize = 3, qvalAsterisks = TRUE)

#----------------------------------------------------------#
# Boxplot summaries of overlay full data
#----------------------------------------------------------#

load("Interim_data/mlPerfAll10k_Allcancer_Shared_2Apr22.RData")
load("Interim_data/mlPerfAll10k_Allcancer_2Apr22.RData")
load("Interim_data/mlPerfAll10k_Allcancer_Cov_2Apr22.RData")

mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short_VSNM <- mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short %>%
  filter(grepl("VSNM",datasetNameShort)) %>% droplevels()
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short_VSNM$datasetName <- mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short_VSNM$datasetNameShort
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short_VSNM <- mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short_VSNM %>% select(-datasetNameShort)
mlPerfAll10k_Allcancer_Overlay_VSNM <- rbind(cbind(mlPerfAll10k_Allcancer, metadataName=NA),
                                        mlPerfAll10k_Allcancer_Shared,
                                        mlPerfAll10k_Allcancer_Cov,
                                        mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short_VSNM) %>%
  filter(grepl("Species|VSNM",datasetName) & !(grepl("\\(CT\\)",datasetName))) %>% droplevels()
mlPerfAll10k_Allcancer_Overlay_VSNM$statGroups <- ifelse(grepl("Species",mlPerfAll10k_Allcancer_Overlay_VSNM$datasetName),
                                                         yes = "Actual", no = "Control")

#-------------------------Plot VSNM overlays-------------------------#
source("00-Functions.R") # for plotControlsVSNM() function
plotControlsVSNM(inputSampleType = "Primary Tumor")
plotControlsVSNM(inputSampleType = "Primary Tumor vs Solid Tissue Normal")
plotControlsVSNM(inputSampleType = "Blood Derived Normal")
