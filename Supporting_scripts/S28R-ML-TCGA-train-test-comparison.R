#-----------------------------------------------------------------------------
# S14-ML-fungi-10k-rep1-tcga-by-seqcenter-taxa-levels-and-WIS-intersect.R
# Copyright (c) 2021--, Greg Poore
# Purpose: Run machine learning on all available TCGA cancer types by seq center
#-----------------------------------------------------------------------------

#-------------------------------#
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

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
load("data_for_ml_tcga_wgs_vs_rna_decontamV2_2Apr22.RData", verbose=T)
load("data_for_ml_tcga_by_seq_center_and_experimental_strategy_taxa_levels_and_wz_intersect_decontamV2_2Apr22.RData", verbose=T)
load("data_for_ml_tcga_by_seq_center_and_experimental_strategy_with_coverage_filter_decontamV2_2Apr22.RData", verbose=T)

# GBM HYPERPARAMETER SEARCH GRID BELOW (DEFAULT PER CARET PACKAGE)
kfoldGBMGrid <- data.frame(n.trees=150, interaction.depth=3,
                                         shrinkage=0.1,
                                         n.minobsinnode=1)

# Model setup -- MODIFY AS NEEDED; ALL THESE COMPARISONS WILL BE TESTED
sampleTypeList <- c("Primary Tumor vs Solid Tissue Normal",
                    "Blood Derived Normal",
                    "Primary Tumor")
# MODIFY AS NEEDED
datasetList <- list( # VSNM
     snmDataOGUFungiDecontamV2,
     # HMS - by taxa level
     df_psRep200_HiSeq_Fungi_DecontamV2_HMS_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_HMS_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_HMS_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_HMS_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_HMS_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_HMS_species,
     # BCM - by taxa level
     df_psRep200_HiSeq_Fungi_DecontamV2_BCM_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_BCM_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_BCM_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_BCM_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_BCM_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_BCM_species,
     # MDA - by taxa level
     df_psRep200_HiSeq_Fungi_DecontamV2_MDA_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_MDA_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_MDA_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_MDA_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_MDA_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_MDA_species,
     # WashU - by taxa level
     df_psRep200_HiSeq_Fungi_DecontamV2_WashU_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_WashU_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_WashU_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_WashU_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_WashU_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_WashU_species,
     # Broad_WGS - by taxa level
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_species,
     # UNC - by taxa level
     df_psRep200_HiSeq_Fungi_DecontamV2_UNC_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_UNC_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_UNC_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_UNC_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_UNC_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_UNC_species,
     # CMS - by taxa level
     df_psRep200_HiSeq_Fungi_DecontamV2_CMS_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_CMS_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_CMS_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_CMS_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_CMS_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_CMS_species,
     # Broad_RNA - by taxa level
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_species,
     # Decontaminated fungi
     rep200_HiSeq_Fungi_DecontamV2_HMS,
     rep200_HiSeq_Fungi_DecontamV2_BCM,
     rep200_HiSeq_Fungi_DecontamV2_MDA,
     rep200_HiSeq_Fungi_DecontamV2_WashU,
     rep200_HiSeq_Fungi_DecontamV2_Broad_WGS,
     rep200_HiSeq_Fungi_DecontamV2_UNC,
     rep200_HiSeq_Fungi_DecontamV2_CMS,
     rep200_HiSeq_Fungi_DecontamV2_Broad_RNA)

datasetListNames <- c(# VSNM
     "snmDataOGUFungiDecontamV2",
     # HMS - by taxa level
     "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_phylum",
     "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_class",
     "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_order",
     "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_family",
     "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_genus",
     "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_species",
     # BCM - by taxa level
     "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_phylum",
     "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_class",
     "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_order",
     "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_family",
     "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_genus",
     "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_species",
     # MDA - by taxa level
     "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_phylum",
     "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_class",
     "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_order",
     "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_family",
     "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_genus",
     "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_species",
     # WashU - by taxa level
     "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_phylum",
     "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_class",
     "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_order",
     "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_family",
     "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_genus",
     "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_species",
     # Broad_WGS - by taxa level
     "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_phylum",
     "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_class",
     "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_order",
     "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_family",
     "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_genus",
     "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_species",
     # UNC - by taxa level
     "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_phylum",
     "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_class",
     "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_order",
     "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_family",
     "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_genus",
     "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_species",
     # CMS - by taxa level
     "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_phylum",
     "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_class",
     "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_order",
     "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_family",
     "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_genus",
     "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_species",
     # Broad_RNA - by taxa level
     "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_phylum",
     "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_class",
     "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_order",
     "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_family",
     "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_genus",
     "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_species",
     # Decontaminated fungi
     "rep200_HiSeq_Fungi_DecontamV2_HMS",
     "rep200_HiSeq_Fungi_DecontamV2_BCM",
     "rep200_HiSeq_Fungi_DecontamV2_MDA",
     "rep200_HiSeq_Fungi_DecontamV2_WashU",
     "rep200_HiSeq_Fungi_DecontamV2_Broad_WGS",
     "rep200_HiSeq_Fungi_DecontamV2_UNC",
     "rep200_HiSeq_Fungi_DecontamV2_CMS",
     "rep200_HiSeq_Fungi_DecontamV2_Broad_RNA")

metadataList <- list(# VSNM
     metaQiitaCombined_Nonzero_DecontamV2,
     # Metadata - HMS  - by taxa level
     metaQiitaCombined_Nonzero_DecontamV2_HMS,
     metaQiitaCombined_Nonzero_DecontamV2_HMS,
     metaQiitaCombined_Nonzero_DecontamV2_HMS,
     metaQiitaCombined_Nonzero_DecontamV2_HMS,
     metaQiitaCombined_Nonzero_DecontamV2_HMS,
     metaQiitaCombined_Nonzero_DecontamV2_HMS,
     # Metadata - BCM - by taxa level
     metaQiitaCombined_Nonzero_DecontamV2_BCM,
     metaQiitaCombined_Nonzero_DecontamV2_BCM,
     metaQiitaCombined_Nonzero_DecontamV2_BCM,
     metaQiitaCombined_Nonzero_DecontamV2_BCM,
     metaQiitaCombined_Nonzero_DecontamV2_BCM,
     metaQiitaCombined_Nonzero_DecontamV2_BCM,
     # Metadata - MDA - by taxa level
     metaQiitaCombined_Nonzero_DecontamV2_MDA,
     metaQiitaCombined_Nonzero_DecontamV2_MDA,
     metaQiitaCombined_Nonzero_DecontamV2_MDA,
     metaQiitaCombined_Nonzero_DecontamV2_MDA,
     metaQiitaCombined_Nonzero_DecontamV2_MDA,
     metaQiitaCombined_Nonzero_DecontamV2_MDA,
     # Metadata - WashU - by taxa level
     metaQiitaCombined_Nonzero_DecontamV2_WashU,
     metaQiitaCombined_Nonzero_DecontamV2_WashU,
     metaQiitaCombined_Nonzero_DecontamV2_WashU,
     metaQiitaCombined_Nonzero_DecontamV2_WashU,
     metaQiitaCombined_Nonzero_DecontamV2_WashU,
     metaQiitaCombined_Nonzero_DecontamV2_WashU,
     # Metadata - Broad_WGS - by taxa level
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS,
     # Metadata - UNC - by taxa level
     metaQiitaCombined_Nonzero_DecontamV2_UNC,
     metaQiitaCombined_Nonzero_DecontamV2_UNC,
     metaQiitaCombined_Nonzero_DecontamV2_UNC,
     metaQiitaCombined_Nonzero_DecontamV2_UNC,
     metaQiitaCombined_Nonzero_DecontamV2_UNC,
     metaQiitaCombined_Nonzero_DecontamV2_UNC,
     # Metadata - CMS - by taxa level
     metaQiitaCombined_Nonzero_DecontamV2_CMS,
     metaQiitaCombined_Nonzero_DecontamV2_CMS,
     metaQiitaCombined_Nonzero_DecontamV2_CMS,
     metaQiitaCombined_Nonzero_DecontamV2_CMS,
     metaQiitaCombined_Nonzero_DecontamV2_CMS,
     metaQiitaCombined_Nonzero_DecontamV2_CMS,
     # Metadata - Broad_RNA - by taxa level
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA,
     # Decontaminated fungi
     metaQiitaCombined_Nonzero_DecontamV2_HMS,
     metaQiitaCombined_Nonzero_DecontamV2_BCM,
     metaQiitaCombined_Nonzero_DecontamV2_MDA,
     metaQiitaCombined_Nonzero_DecontamV2_WashU,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS,
     metaQiitaCombined_Nonzero_DecontamV2_UNC,
     metaQiitaCombined_Nonzero_DecontamV2_CMS,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA)

metadataListNames <- list(# VSNM
     "metaQiitaCombined_Nonzero_DecontamV2",
     # Metadata - HMS - by taxa level
     "metaQiitaCombined_Nonzero_DecontamV2_HMS",
     "metaQiitaCombined_Nonzero_DecontamV2_HMS",
     "metaQiitaCombined_Nonzero_DecontamV2_HMS",
     "metaQiitaCombined_Nonzero_DecontamV2_HMS",
     "metaQiitaCombined_Nonzero_DecontamV2_HMS",
     "metaQiitaCombined_Nonzero_DecontamV2_HMS",
     # Metadata - BCM - by taxa level
     "metaQiitaCombined_Nonzero_DecontamV2_BCM",
     "metaQiitaCombined_Nonzero_DecontamV2_BCM",
     "metaQiitaCombined_Nonzero_DecontamV2_BCM",
     "metaQiitaCombined_Nonzero_DecontamV2_BCM",
     "metaQiitaCombined_Nonzero_DecontamV2_BCM",
     "metaQiitaCombined_Nonzero_DecontamV2_BCM",
     # Metadata - MDA - by taxa level
     "metaQiitaCombined_Nonzero_DecontamV2_MDA",
     "metaQiitaCombined_Nonzero_DecontamV2_MDA",
     "metaQiitaCombined_Nonzero_DecontamV2_MDA",
     "metaQiitaCombined_Nonzero_DecontamV2_MDA",
     "metaQiitaCombined_Nonzero_DecontamV2_MDA",
     "metaQiitaCombined_Nonzero_DecontamV2_MDA",
     # Metadata - WashU - by taxa level
     "metaQiitaCombined_Nonzero_DecontamV2_WashU",
     "metaQiitaCombined_Nonzero_DecontamV2_WashU",
     "metaQiitaCombined_Nonzero_DecontamV2_WashU",
     "metaQiitaCombined_Nonzero_DecontamV2_WashU",
     "metaQiitaCombined_Nonzero_DecontamV2_WashU",
     "metaQiitaCombined_Nonzero_DecontamV2_WashU",
     # Metadata - Broad_WGS - by taxa level
     "metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS",
     "metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS",
     "metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS",
     "metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS",
     "metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS",
     "metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS",
     # Metadata - UNC - by taxa level
     "metaQiitaCombined_Nonzero_DecontamV2_UNC",
     "metaQiitaCombined_Nonzero_DecontamV2_UNC",
     "metaQiitaCombined_Nonzero_DecontamV2_UNC",
     "metaQiitaCombined_Nonzero_DecontamV2_UNC",
     "metaQiitaCombined_Nonzero_DecontamV2_UNC",
     "metaQiitaCombined_Nonzero_DecontamV2_UNC",
     # Metadata - CMS - by taxa level
     "metaQiitaCombined_Nonzero_DecontamV2_CMS",
     "metaQiitaCombined_Nonzero_DecontamV2_CMS",
     "metaQiitaCombined_Nonzero_DecontamV2_CMS",
     "metaQiitaCombined_Nonzero_DecontamV2_CMS",
     "metaQiitaCombined_Nonzero_DecontamV2_CMS",
     "metaQiitaCombined_Nonzero_DecontamV2_CMS",
     # Metadata - Broad_RNA - by taxa level
     "metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA",
     "metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA",
     "metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA",
     "metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA",
     "metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA",
     "metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA",
     # Decontaminated fungi
     "metaQiitaCombined_Nonzero_DecontamV2_HMS",
     "metaQiitaCombined_Nonzero_DecontamV2_BCM",
     "metaQiitaCombined_Nonzero_DecontamV2_MDA",
     "metaQiitaCombined_Nonzero_DecontamV2_WashU",
     "metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS",
     "metaQiitaCombined_Nonzero_DecontamV2_UNC",
     "metaQiitaCombined_Nonzero_DecontamV2_CMS",
     "metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA")

baseNamePerDatasetResultsFile <- "perfFungi_10k_rep1_tcga_gbm_train_test_ALL_6Apr22"
baseNameAllResultsFile <- "rep_perfFungi_10k_rep1_tcga_gbm_train_test_ALL_6Apr22"

# MATCHED METADATA DATA FRAME FOR ALL COUNT DATA:
caretTuneGrid <- kfoldGBMGrid
numKFold <- 4
numResampleIter <- 1
prroc_roc <- list()
prroc_pr <- list()
perf <- list()
rep_perf <- list()
perfTmp <- list()
rep_perfTmp <- list()
perfTmp2 <- list()
rep_perfTmp2 <- list()

for(jj in seq_along(datasetList)){

  dataTmp <- datasetList[[jj]]
  datasetName <- datasetListNames[[jj]]

  metaTmpQC <- metadataList[[jj]]
  metaTmpQC$disease_type <- factor(metaTmpQC$disease_type)
  metadataName <- metadataListNames[[jj]]

  for(kk in seq_along(sampleTypeList)){
    st <- sampleTypeList[kk]
    print(st)
    
    if(st == "Primary Tumor vs Solid Tissue Normal"){
      metaTmp <- metaTmpQC
      metaTmp2 <- droplevels(metaTmp[metaTmp$sample_type %in% c("Primary Tumor",
                                                                "Solid Tissue Normal"),])
    } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
      metaTmp <- metaTmpQC
      metaTmp2 <- droplevels(metaTmp[metaTmp$sample_type %in% st,])
    } else if(st == "Stage I vs IV"){
      metaTmp <- metaTmpPath
      metaTmp2 <- droplevels(metaTmp[(metaTmp$sample_type %in% "Primary Tumor") &
                                       (metaTmp$pathologic_stage_label_binned %in% c("Stage1","Stage4")),])
    }
    
    print(seq_along(levels(metaTmp2$disease_type)))
    
    for(ii in seq_along(levels(metaTmp2$disease_type))){
      
      dt <- levels(metaTmp2$disease_type)[ii]
      print(dt)
      
      if(st == "Primary Tumor vs Solid Tissue Normal"){
        metaTmp3 <- droplevels(metaTmp2[metaTmp2$disease_type == dt,])
        metaTmp3$predY <- factor(gsub('([[:punct:]])|\\s+','',metaTmp3$sample_type))
        positiveClass <- "PrimaryTumor"
        negativeClass <- "SolidTissueNormal"
      } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
        metaTmp3 <- metaTmp2
        metaTmp3$predY <- factor(ifelse(metaTmp2$disease_type == dt, 
                                        yes = dt, 
                                        no = "OtherCancerType"),
                                 levels = c(dt, "OtherCancerType"))
        positiveClass <- gsub('([[:punct:]])|\\s+','',dt)
        negativeClass <- "OtherCancerType"
      } else if(st == "Stage I vs IV"){
        metaTmp3 <- droplevels(metaTmp2[metaTmp2$disease_type == dt,])
        metaTmp3$predY <- metaTmp3$pathologic_stage_label_binned
        positiveClass <- "Stage4"
        negativeClass <- "Stage1"
      }
      
      print(table(metaTmp3$predY))
      
      # SKIP CANCER TYPES THAT ONLY HAVE ONE CLASS OF A COMPARISON
      if(length(table(metaTmp3$predY)) < 2){next}
      
      # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF 20 SAMPLES IN EITHER CLASS
      if(any(table(metaTmp3$predY) < 20)){next}
      
      minorityClassSize <- min(table((metaTmp3$predY)))
      majorityClassSize <- max(table((metaTmp3$predY)))
      
      minorityClassName <- names(table(metaTmp3$predY)[which(table(metaTmp3$predY) == min(table(metaTmp3$predY)))])
      majorityClassName <- names(table(metaTmp3$predY)[which(table(metaTmp3$predY) == max(table(metaTmp3$predY)))])
      
      mlDataY <- metaTmp3
      mlDataX <- dataTmp[rownames(mlDataY),]

      repX_perf <- list()
      for(pp in 1:10){
        # USE 90% OF DATA FOR TRAINING AND 10% FOR TESTING
        set.seed(pp)
        index <- createDataPartition(mlDataY$predY, p = 0.9, list = FALSE)
        trainX <- mlDataX[index,]
        trainY <- mlDataY[index,"predY"]
        refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))

        testX <- mlDataX[-index,]
        testY <- mlDataY[-index,"predY"]
        refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
        
        set.seed(pp)
        ctrl <- trainControl(method = "repeatedcv",
                             number = numKFold,
                             repeats = numResampleIter,
                             sampling = "up",
                             summaryFunction = twoClassSummary,
                             classProbs = TRUE,
                             # verboseIter = TRUE,
                             savePredictions = TRUE,
                             allowParallel=TRUE)
        
        mlModel <- train(x = trainX,
                         y = refactoredTrainY,
                         method = "gbm",
                         # preProcess = c("scale","center"),
                         trControl = ctrl,
                         # verbose = TRUE,
                         metric = "ROC",
                         tuneGrid = caretTuneGrid)

        predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
        fg <- predProbs[refactoredTestY == positiveClass]
        bg <- predProbs[refactoredTestY == negativeClass]
        
        rep_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
        rep_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
        
        multiClass <- data.frame(obs = refactoredTestY,
                                 pred = predict(mlModel, newdata = testX),
                                 predict(mlModel, newdata = testX, type = "prob"))
        repX_perf[[pp]] <- data.frame(auroc=rep_roc$auc,
                               aupr=rep_pr$auc.integral,
                               rep=paste0("Fold",pp), 
                               diseaseType = dt,
                               sampleType = st,
                               datasetName = datasetName,
                               metadataName = metadataName,
                               minorityClassSize = minorityClassSize,
                               majorityClassSize = majorityClassSize,
                               minorityClassName = minorityClassName,
                               majorityClassName = majorityClassName)

        rm(mlModel)
      }
      rep_perf[[ii]] <- do.call(rbind, repX_perf)
      
      #--------------------------------------#
      # Save performance into relevant files #
      #--------------------------------------#
      
      filepathPerfStatsPerFold <- paste0("./perfDataPerFold__",datasetName)
      
      if(!( dir.exists( file.path(filepathPerfStatsPerFold)))){
        dir.create(file.path(filepathPerfStatsPerFold))
      }

      filenamePerFoldCSV <- paste0(filepathPerfStatsPerFold,"/",
                            dt,
                            " -- ",
                            st,
                            " -- PerfPerFold.csv")
      
      write.csv(rep_perf[[ii]], file = filenamePerFoldCSV)
      
    }
    # SUMMARIZE RESULTS
    rep_perfTmp[[kk]] <- do.call(rbind, rep_perf)
  }
  # SUMMARIZE RESULTS
  rep_perfTmp2[[jj]] <- do.call(rbind, rep_perfTmp)

  write.csv(rep_perfTmp2[[jj]], file = paste0("rep_",baseNamePerDatasetResultsFile,"_",datasetName,".csv"))
}

# SUMMARIZE RESULTS
rep_perf1VsAll <- do.call(rbind, rep_perfTmp2)

write.csv(rep_perf1VsAll, file = paste0("rep_",baseNameAllResultsFile,".csv"))

#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------