#-----------------------------------------------------------------------------
# 13R-ML-fungi-vs-bacteria.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Subsample # of species (fungi vs bacteria) and compare multi-class performance differences
# - Subsample # of reads (fungi vs bacteria) and compare multi-class performance differences
# - NOTE: Due to non-uniform polyA tail availability in bacteria and uniform availability thereof in fungi, only WGS data is used
#-----------------------------------------------------------------------------

#----------------------------------------------------------#
# Load environments
#----------------------------------------------------------#

# Load dependencies
require(splitstackshape)
require(reshape2)
require(tidyr)
require(caret) # for model building
require(dplyr) # for data manipulation
require(doMC) # for parallel computing
require(gbm) # for machine learning
require(tibble) # for df operations
require(PRROC) # for precision-recall curves
require(MLmetrics) # for multi-class learning
require(cvAUC)
require(reshape2)
require(ggpubr)
require(ggsci)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Load TCGA fungi data
#----------------------------------------------------------#

load("Interim_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_25Mar22.RData", verbose = T)
load("Interim_data/snmDataFungi_DecontamV2_25Mar22.RData", verbose = T)
rep200TaxSplit <- read.csv("Supporting_data/rep200_lineage_map_split.csv", stringsAsFactors = FALSE, row.names = 1)

# - Develop multi-class ML code
# - Subsample # species for 1000 iterations. n=seq(25, 300, by = 25)
# - Rarefy # of reads for 10 iterations each time. n=100, 500, 1000, 5000, 10e3, 15e3, 20e3, 25e3, 30e3, 35e3

metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_WGS <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>% filter(experimental_strategy == "WGS") %>% droplevels()
rep200Data_WGS_RNA_Matched_Fungi_WGS <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_WGS),]
rep200Data_WGS_RNA_Matched_Bacteria_WGS <- rep200Data_WGS_RNA_Matched_Bacteria[rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_WGS),]

quantile(rowSums(rep200Data_WGS_RNA_Matched_Fungi_WGS), probs = seq(0.05, 1, by = 0.10))
# 5%      15%      25%      35%      45%      55%      65%      75%      85%      95% 
#   1048.00  2893.25  4744.00  6542.75  8725.50 11511.25 14581.75 17983.00 22953.75 37061.50
quantile(rowSums(rep200Data_WGS_RNA_Matched_Bacteria_WGS), probs = seq(0.05, 1, by = 0.10))
# 5%       15%       25%       35%       45%       55%       65%       75%       85%       95% 
#   7331.00  21926.50  34709.25  47844.75  63162.00  81523.00 105713.25 139862.75 212465.00 827940.50

#----------------------------------------------------------#
# Subset data
#----------------------------------------------------------#

# Subset to WGS and remove 2 fungi samples with 0 counts: 
# 13722.58cfa82ee4b0c9d6adf6a982, 13722.58cfa82fe4b0c9d6adf6b39c
# Note: None of the bacterial WGS data samples were nonzero
metaFungiVsBacteria <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>%
  filter(cgc_platform == "Illumina HiSeq", experimental_strategy == "WGS") %>% 
  filter(sample_name != "58cfa82ee4b0c9d6adf6a982") %>%
  filter(sample_name != "58cfa82fe4b0c9d6adf6b39c") %>% droplevels()
countsFungiWGS <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaFungiVsBacteria),]
countsBacteriaWGS <- rep200Data_WGS_RNA_Matched_Bacteria[rownames(metaFungiVsBacteria),]

save(metaFungiVsBacteria,
     countsFungiWGS,
     countsBacteriaWGS,
     file = "Interim_data/data_for_ml_fungi_vs_bacteria_numFeat_5Apr22.RData")

#----------------------------------------------------------#
# Identify feature numbers for WIS intersected fungi
#----------------------------------------------------------#

load("Interim_data/decontamResultsV2_25Mar22.RData")

ogusWISintersect <- decontamResultsV2 %>% filter(shared_with_WIS == "YES") %>% rownames()
indexWISfungi <- which(colnames(countsFungiWGS) %in% ogusWISintersect)
inverseIndexNonWISfungi <- which(!(1:319) %in% indexWISfungi)

save(indexWISfungi, inverseIndexNonWISfungi,
     file = "Interim_data/indices_WIS_and_nonWIS_fungi_for_simulations_5Apr22.RData")

#----------------------------------------------------------#
# Subset feature sets to WIS intersected fungi, bacteria, fungi+bacteria
#----------------------------------------------------------#
load("Interim_data/decontamResultsV2_25Mar22.RData")
load("Interim_data/shared_bacterial_features_at_each_taxa_level_29Mar22.RData")

ogusFungiWISintersect <- decontamResultsV2 %>% filter(shared_with_WIS == "YES") %>% rownames()
ogusBacteriaWISintersect <- unique(sharedSpeciesBacteria$intersectedOGUs)

countsFungiWGS_Shared <- countsFungiWGS[,colnames(countsFungiWGS) %in% ogusFungiWISintersect]
countsBacteriaWGS_Shared <- countsBacteriaWGS[,colnames(countsBacteriaWGS) %in% ogusBacteriaWISintersect]
dim(countsFungiWGS_Shared) # 4386   34
dim(countsBacteriaWGS_Shared) # 4386  267

# Sanity check
all(rownames(countsFungiWGS_Shared) == rownames(countsBacteriaWGS_Shared)) # TRUE
countsFungiBacteriaWGS_Shared <- cbind(countsFungiWGS_Shared,countsBacteriaWGS_Shared)
dim(countsFungiBacteriaWGS_Shared) # 4386  301

save(metaFungiVsBacteria,
     countsFungiWGS_Shared,
     countsBacteriaWGS_Shared,
     countsFungiBacteriaWGS_Shared,
     file = "Interim_data/data_for_ml_tcga_synergy_fungi_and_bacteria_5Apr22.RData")

#----------------------------------------------------------#
# Identify index for pan-seq-center WGS 70/30 split
#----------------------------------------------------------#
require(data.table)
require(splitstackshape)

splitVars <- c("data_submitting_center_label","sample_type","disease_type")
set.seed(42)
stratSamplingFungiVsBacteriaNumReads <- stratified(metaFungiVsBacteria, 
                                                   group = splitVars, 
                                                   size = 0.7,
                                                   keep.rownames = TRUE, 
                                                   bothSets = TRUE, 
                                                   replace = FALSE)
split1MetaFungiVsBacteria <- as.data.frame(stratSamplingFungiVsBacteriaNumReads[[1]])
rownames(split1MetaFungiVsBacteria) <- split1MetaFungiVsBacteria$rn
split2MetaFungiVsBacteria <- as.data.frame(stratSamplingFungiVsBacteriaNumReads[[2]])
rownames(split2MetaFungiVsBacteria) <- split2MetaFungiVsBacteria$rn

indexAllSeqCenterTrain <- which(rownames(metaFungiVsBacteria) %in% rownames(split1MetaFungiVsBacteria))
indexAllSeqCenterTest <- which(rownames(metaFungiVsBacteria) %in% rownames(split2MetaFungiVsBacteria))

# Load tax tables
rep200TaxSplit <- read.csv("Supporting_data/rep200_lineage_map_split.csv", stringsAsFactors = FALSE, row.names = 1)
rep200Kingdoms <- read.csv("Supporting_data/rep200_gOTU_kingdom_mapping.csv", stringsAsFactors = FALSE)
rep200Kingdoms_Fungi <- rep200Kingdoms[which(rep200Kingdoms$kingdom == "fungi"),]
rep200Kingdoms_Bacteria <- rep200Kingdoms[which(rep200Kingdoms$kingdom == "bacteria"),]

rep200TaxSplit_Fungi <- rep200TaxSplit[rownames(rep200TaxSplit) %in% rep200Kingdoms_Fungi$genomeID,,drop=FALSE]
rep200TaxSplit_Bacteria <- rep200TaxSplit[rownames(rep200TaxSplit) %in% rep200Kingdoms_Bacteria$genomeID,,drop=FALSE]

save(metaFungiVsBacteria,
     countsFungiWGS,
     countsBacteriaWGS,
     split1MetaFungiVsBacteria,
     split2MetaFungiVsBacteria,
     rep200TaxSplit_Fungi,
     rep200TaxSplit_Bacteria,
     file = "Interim_data/data_for_ml_fungi_vs_bacteria_numFeat_with_AllSeqCenter_index_5Apr22.RData")

#----------------------------------------------------------#
# Plot results of TCGA synergy (WIS fungi, WIS bacteria, WIS fungi+bacteria) simulations
#----------------------------------------------------------#

# NOTE: Script S21R was run to obtain the results for plotting
load("00-Functions.R") # for plotSynergyPerf() plotting function
# Primary Tumor
plotSynergyPerf(rDataFilePath = "Interim_data/tcga_CV_synergy_fungi_bacteria_PrimaryTumor_AllSeqCenters_numIter50_k10_modelType_xgbTree.RData", 
                modelType = "xgbTree", numIter = 50, inputSampleType = "Primary Tumor")
# Blood Derived Normal
plotSynergyPerf(rDataFilePath = "Interim_data/tcga_CV_synergy_fungi_bacteria_BloodDerivedNormal_AllSeqCenters_numIter50_k10_modelType_xgbTree.RData", 
                modelType = "xgbTree", numIter = 50, inputSampleType = "Blood Derived Normal")

# # Run on the cluster using the above function with the exact=TRUE for the wilcox
# plotSynergyPerf(rDataFilePath = "tcga_CV_synergy_fungi_bacteria_PrimaryTumor_AllSeqCenters_numIter50_k10_modelType_xgbTree.RData",
#                 modelType = "xgbTree", numIter = 50, inputSampleType = "Primary Tumor")
# plotSynergyPerf(rDataFilePath = "tcga_CV_synergy_fungi_bacteria_BloodDerivedNormal_AllSeqCenters_numIter50_k10_modelType_xgbTree.RData",
#                 modelType = "xgbTree", numIter = 50, inputSampleType = "Blood Derived Normal")

#----------------------------------------------------------#
# Plot results of WIS-intersected features vs non-intersected features
#----------------------------------------------------------#

# NOTE: Script S22R was run to obtain the results for plotting
load("00-Functions.R") # for plotWISVsNonWISFungiPerf() plotting function

# Primary Tumor
plotWISVsNonWISFungiPerf(rDataFilePath = "Interim_data/wis_vs_nonwis_fungi_CV_PrimaryTumor_AllSeqCenters_numIter50_k10_modelType_xgbTree.RData", 
                         seqCenter = "AllSeqCenters", allSeqCenterFlag = TRUE, modelType = "xgbTree", numIter = 50, inputSampleType = "Primary Tumor", statSize=3)
# Blood Derived Normal
plotWISVsNonWISFungiPerf(rDataFilePath = "Interim_data/wis_vs_nonwis_fungi_CV_BloodDerivedNormal_AllSeqCenters_numIter50_k10_modelType_xgbTree.RData", 
                         seqCenter = "AllSeqCenters", allSeqCenterFlag = TRUE, modelType = "xgbTree", numIter = 50, inputSampleType = "Blood Derived Normal", statSize=3)
