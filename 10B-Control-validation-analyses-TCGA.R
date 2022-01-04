#-----------------------------------------------------------------------------
# 10-Control-validation-analyses-TCGA.R
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

# snmDataOGUFungiDecontamV2,
# vdge_dataE_DecontamV2,
# rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero,
# metaQiitaCombined_Nonzero_DecontamV2,
load("Interim_data/snmDataFungi_DecontamV2_13Oct21.RData")
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
  ylim(c(0,floor(1.1*max(table(centerCompare$SeqCenter))/2))) +
  labs(x = "Sequencing Center", y = "Count", title = "Validation split - Sequence center distribution")
ggplot(centerCompare, aes(x = SampleType, fill = Split)) + geom_bar(position = "dodge",stat = "count") +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) + scale_fill_nejm() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5, color="black", position = position_dodge(0.9)) + 
  ylim(c(0,floor(1.1*max(table(centerCompare$SampleType))/2))) +
  labs(x = "Sample Type", y = "Count", title = "Validation split - Sample type distribution")
ggplot(centerCompare, aes(x = DiseaseType, fill = Split)) + geom_bar(position = "dodge",stat = "count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_nejm() +
  geom_text(stat='count', aes(label=..count..), hjust=-0.25, color="black", position = position_dodge(0.9), angle=90) + 
  ylim(c(0,floor(1.1*max(table(centerCompare$DiseaseType))/2))) +
  labs(x = "Sample Type", y = "Count", title = "Validation split - Disease type distribution")
ggplot(centerCompare, aes(x = PathStage, fill = Split)) + geom_bar(position = "dodge",stat = "count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_nejm() +
  geom_text(stat='count', aes(label=..count..), hjust=-0.25, color="black", position = position_dodge(0.9), angle=90) + 
  ylim(c(0,floor(1.1*max(table(centerCompare$PathStage))/2))) +
  labs(x = "Sample Type", y = "Count", title = "Validation split - Pathologic Stage distribution")

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

## TCGA full data - biological variables: disease_type + sample_type | technical variables: data_submitting_center_label + experimental_strategy
# Split 1
valSplit1FungiCountDecontamV2Data_VSNM_Obj_CT <- vsnmFunctionTCGA(qcData = valSplit1FungiCountDecontamV2Data, qcMetadata = valSplit1Metadata, cancerTypeFlag = TRUE) # converges
valSplit1FungiCountDecontamV2Data_VSNM_CT <- valSplit1FungiCountDecontamV2Data_VSNM_Obj$snmData
# Split 2
valSplit2FungiCountDecontamV2Data_VSNM_Obj_CT <- vsnmFunctionTCGA(qcData = valSplit2FungiCountDecontamV2Data, qcMetadata = valSplit2Metadata, cancerTypeFlag = TRUE) # converges
valSplit2FungiCountDecontamV2Data_VSNM_CT <- valSplit2FungiCountDecontamV2Data_VSNM_Obj$snmData
# Full data
valFullFungiCountDecontamV2Data_VSNM_Obj_CT <- vsnmFunctionTCGA(qcData = rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero, qcMetadata = metaQiitaCombined_Nonzero_DecontamV2, cancerTypeFlag = TRUE) # converges
valFullFungiCountDecontamV2Data_VSNM_CT <- valFullFungiCountDecontamV2Data_VSNM_Obj$snmData

save(valSplit1Metadata,
     valSplit1FungiCountDecontamV2Data,
     valSplit1FungiCountDecontamV2Data_VSNM,
     valSplit1FungiCountDecontamV2Data_VSNM_CT,
     
     valSplit2Metadata,
     valSplit2FungiCountDecontamV2Data,
     valSplit2FungiCountDecontamV2Data_VSNM,
     valSplit2FungiCountDecontamV2Data_VSNM_CT,
     
     metaQiitaCombined_Nonzero_DecontamV2,
     rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero,
     valFullFungiCountDecontamV2Data_VSNM,
     valFullFungiCountDecontamV2Data_VSNM_CT,
     file = "Interim_data/validation_2halves_tcga_vsnm_data_25Oct21.RData")

#----------------------------------------------------------#
# Build machine learning models on each split and half of full dataset
#----------------------------------------------------------#

valSplit1MetadataTest <- valSplit1Metadata %>% filter(disease_type %in% c("Breast Invasive Carcinoma", "Lung Adenocarcinoma", "Stomach Adenocarcinoma")) %>% droplevels()
valSplit2MetadataTest <- valSplit2Metadata %>% filter(disease_type %in% c("Breast Invasive Carcinoma", "Lung Adenocarcinoma", "Stomach Adenocarcinoma")) %>% droplevels()
valFullMetadataTest <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(disease_type %in% c("Breast Invasive Carcinoma", "Lung Adenocarcinoma", "Stomach Adenocarcinoma")) %>% droplevels()

valSplit1FungiCountDecontamV2Data_VSNM_Test <- valSplit1FungiCountDecontamV2Data_VSNM[rownames(valSplit1MetadataTest),]
valSplit2FungiCountDecontamV2Data_VSNM_Test <- valSplit2FungiCountDecontamV2Data_VSNM[rownames(valSplit2MetadataTest),]
valFullFungiCountDecontamV2Data_VSNM_Test <- valFullFungiCountDecontamV2Data_VSNM[rownames(valFullMetadataTest),]

valSplit1FungiCountDecontamV2Data_VSNM_CT_Test <- valSplit1FungiCountDecontamV2Data_VSNM_CT[rownames(valSplit1MetadataTest),]
valSplit2FungiCountDecontamV2Data_VSNM_CT_Test <- valSplit2FungiCountDecontamV2Data_VSNM_CT[rownames(valSplit2MetadataTest),]
valFullFungiCountDecontamV2Data_VSNM_CT_Test <- valFullFungiCountDecontamV2Data_VSNM_CT[rownames(valFullMetadataTest),]

## Write function
crossTestSplitML <- function(metadataSplit1,
                            metadataSplit2,
                            metadataFull,
                            vsnmDataSplit1,
                            vsnmDataSplit2,
                            vsnmDataFull,
                            sampleTypeComparison = "Primary Tumor",
                            caretModel = "gbm",
                            samplingStrategy = "up",
                            numResampleIter = 1,
                            numKFold = 4,
                            baseName = "val_halves_k4",
                            caretTuneGrid = data.frame(n.trees=150, interaction.depth=3,
                                                       shrinkage=0.1,
                                                       n.minobsinnode=1),
                            sampleTypeList = c("Primary Tumor")){
  require(caret) # for model building
  require(gbm) # for machine learning
  require(PRROC) # for precision-recall curves
  require(MLmetrics) # for multiclass ML
  require(e1071)
  require(ggpubr)
  require(ggrepel)
  
  perf <- list()
  perfTmp <- list()
  perfTmp2 <- list()
  
  for(kk in seq_along(sampleTypeList)){
    st <- sampleTypeList[kk]
    print(st)
    
    if(st == "Primary Tumor vs Solid Tissue Normal"){
      metadataSplit1_FiltSt <- droplevels(metadataSplit1[metadataSplit1$sample_type %in%
                                                         c("Primary Tumor", "Solid Tissue Normal"),])
      metadataSplit2_FiltSt <- droplevels(metadataSplit2[metadataSplit2$sample_type %in%
                                                         c("Primary Tumor", "Solid Tissue Normal"),])
      metadataFull_FiltSt <- droplevels(metadataFull[metadataFull$sample_type %in%
                                                         c("Primary Tumor", "Solid Tissue Normal"),])
    } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
      metadataSplit1_FiltSt <- droplevels(metadataSplit1[metadataSplit1$sample_type %in% st,])
      metadataSplit2_FiltSt <- droplevels(metadataSplit2[metadataSplit2$sample_type %in% st,])
      metadataFull_FiltSt <- droplevels(metadataFull[metadataFull$sample_type %in% st,])
    }
    
    # print(seq_along(levels(metaTmp2$disease_type)))
    
    for(ii in seq_along(levels(metadataSplit1_FiltSt$disease_type))){
      
      dt <- levels(metadataSplit1_FiltSt$disease_type)[ii]
      print(dt)
      
      if(st == "Primary Tumor vs Solid Tissue Normal"){
        metadataSplit1_FiltStDz <- droplevels(metadataSplit1_FiltSt[metadataSplit1_FiltSt$disease_type == dt,])
        metadataSplit2_FiltStDz <- droplevels(metadataSplit2_FiltSt[metadataSplit2_FiltSt$disease_type == dt,])
        metadataFull_FiltStDz <- droplevels(metadataFull_FiltSt[metadataFull_FiltSt$disease_type == dt,])
        
        metadataSplit1_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataSplit1_FiltStDz$sample_type))
        metadataSplit2_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataSplit2_FiltStDz$sample_type))
        metadataFull_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataFull_FiltStDz$sample_type))
        
        positiveClass <- "PrimaryTumor"
        negativeClass <- "SolidTissueNormal"
      } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
        metadataSplit1_FiltStDz <- metadataSplit1_FiltSt
        metadataSplit2_FiltStDz <- metadataSplit2_FiltSt
        metadataFull_FiltStDz <- metadataFull_FiltSt
        
        metadataSplit1_FiltStDz$predY <- ifelse(metadataSplit1_FiltStDz$disease_type == dt, yes = dt, no = "OtherCancerType")
        metadataSplit2_FiltStDz$predY <- ifelse(metadataSplit2_FiltStDz$disease_type == dt, yes = dt, no = "OtherCancerType")
        metadataFull_FiltStDz$predY <- ifelse(metadataFull_FiltStDz$disease_type == dt, yes = dt, no = "OtherCancerType")
        
        metadataSplit1_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataSplit1_FiltStDz$predY))
        metadataSplit2_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataSplit2_FiltStDz$predY))
        metadataFull_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataFull_FiltStDz$predY))
        
        positiveClass <- gsub('([[:punct:]])|\\s+','',dt)
        negativeClass <- "OtherCancerType"
      }
      
      # print(table(metadataSplit1_FiltStDz$predY))
      # print(table(metadataSplit2_FiltStDz$predY))
      # print(table(metadataFull_FiltStDz$predY))
      
      # SKIP CANCER TYPES THAT ONLY HAVE ONE CLASS OF A COMPARISON
      if(length(table(metadataSplit1_FiltStDz$predY)) < 2 | length(table(metadataSplit2_FiltStDz$predY)) < 2){next}
      
      # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF 20 SAMPLES IN EITHER CLASS
      if(any(table(metadataSplit1_FiltStDz$predY) < 20) | any(table(metadataSplit2_FiltStDz$predY) < 20)){next}
      
      mlDataY_Split1 <- metadataSplit1_FiltStDz
      mlDataY_Split2 <- metadataSplit2_FiltStDz
      mlDataY_Full <- metadataFull_FiltStDz
      
      mlDataX_Split1 <- vsnmDataSplit1[rownames(mlDataY_Split1),]
      mlDataX_Split2 <- vsnmDataSplit2[rownames(mlDataY_Split2),]
      mlDataX_Full <- vsnmDataFull[rownames(mlDataY_Full),]
      
      set.seed(42)
      ctrl <- trainControl(method = "repeatedcv",
                           number = numKFold,
                           repeats = numResampleIter,
                           sampling = "up",
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE,
                           verboseIter = FALSE,
                           savePredictions = TRUE,
                           allowParallel=TRUE)
      
      #--------------------------Train split 1--------------------------#
      set.seed(42)
      trainX_Split1 <- mlDataX_Split1
      trainY_Split1 <- mlDataY_Split1[,"predY"]
      testX_Split1 <- mlDataX_Split2
      testY_Split1 <- mlDataY_Split2[,"predY"]
      # refactoredTrainY_Split1 <- factor(gsub('([[:punct:]])|\\s+','',trainY_Split1))
      
      print("Now training model on split 1...")
      set.seed(42)
      mlModel_Split1 <- train(x = trainX_Split1,
                       y = trainY_Split1,
                       method = "gbm",
                       preProcess = c("scale","center"),
                       trControl = ctrl,
                       verbose = FALSE,
                       metric = "ROC",
                       tuneGrid = caretTuneGrid)
      
      #--------------------------Train split 2--------------------------#
      set.seed(42)
      trainX_Split2 <- mlDataX_Split2
      trainY_Split2 <- mlDataY_Split2[,"predY"]
      testX_Split2 <- mlDataX_Split1
      testY_Split2 <- mlDataY_Split1[,"predY"]
      # refactoredTrainY_Split2 <- factor(gsub('([[:punct:]])|\\s+','',trainY_Split2))
      
      print("Now training model on split 2...")
      set.seed(42)
      mlModel_Split2 <- train(x = trainX_Split2,
                              y = trainY_Split2,
                              method = "gbm",
                              preProcess = c("scale","center"),
                              trControl = ctrl,
                              verbose = FALSE,
                              metric = "ROC",
                              tuneGrid = caretTuneGrid)
      
      #--------------------------Train half of full data--------------------------#
      
      set.seed(42) # random number seed -- important if you want comparable results between data types
      index <- createDataPartition(mlDataY_Full[,"predY"], p = 0.5, list = FALSE)
      trainHalfFullX <- mlDataX_Full[index,]
      trainHalfFullY <- mlDataY_Full[index,"predY"]
      testHalfFullX <- mlDataX_Full[-index,]
      testHalfFullY <- mlDataY_Full[-index,"predY"]
      
      print("Now training model on half of full data...")
      set.seed(42)
      mlModel_HalfFull <- train(x = trainHalfFullX,
                              y = trainHalfFullY,
                              method = "gbm",
                              preProcess = c("scale","center"),
                              trControl = ctrl,
                              verbose = FALSE,
                              metric = "ROC",
                              tuneGrid = caretTuneGrid)
      
      #--------------------------Test on opposite split--------------------------#
      print("Now testing models on each other...")
      # Use model from split 1 to predict on split 2
      predClassS1S2 <- predict(mlModel_Split1, newdata = testX_Split1)
      predProbsS1S2 <- as.numeric(predict(mlModel_Split1, newdata = testX_Split1, type = "prob")[,positiveClass])
      # Use model from split 2 to predict on split 1
      predClassS2S1 <- predict(mlModel_Split2, newdata = testX_Split2)
      predProbsS2S1 <- as.numeric(predict(mlModel_Split2, newdata = testX_Split2, type = "prob")[,positiveClass])
      # Use model from half of full data to predict on holdout half
      predClassHalfFull <- predict(mlModel_HalfFull, newdata = testHalfFullX)
      predProbsHalfFull <- as.numeric(predict(mlModel_HalfFull, newdata = testHalfFullX, type = "prob")[,positiveClass])
      
      # Calculate perf of testing model from split 1 to predict on split 2
      fgS1S2 <- predProbsS1S2[testY_Split1 == positiveClass]
      bgS1S2 <- predProbsS1S2[testY_Split1 == negativeClass]
      prroc_rocS1S2 <- roc.curve(scores.class0 = fgS1S2, scores.class1 = bgS1S2, curve = T)
      prroc_prS1S2 <- pr.curve(scores.class0 = fgS1S2, scores.class1 = bgS1S2, curve = T, rand.compute=T)
      aurocS1S2 = prroc_rocS1S2$auc
      auprS1S2 = prroc_prS1S2$auc.integral
      # Calculate perf of testing model from split 2 to predict on split 1
      fgS2S1 <- predProbsS2S1[testY_Split2 == positiveClass]
      bgS2S1 <- predProbsS2S1[testY_Split2 == negativeClass]
      prroc_rocS2S1 <- roc.curve(scores.class0 = fgS2S1, scores.class1 = bgS2S1, curve = T)
      prroc_prS2S1 <- pr.curve(scores.class0 = fgS2S1, scores.class1 = bgS2S1, curve = T, rand.compute=T)
      aurocS2S1 = prroc_rocS2S1$auc
      auprS2S1 = prroc_prS2S1$auc.integral
      # Calculate perf of testing model from half of full data to predict on holdout half
      fgHalfFull <- predProbsHalfFull[testHalfFullY == positiveClass]
      bgHalfFull <- predProbsHalfFull[testHalfFullY == negativeClass]
      prroc_rocHalfFull <- roc.curve(scores.class0 = fgHalfFull, scores.class1 = bgHalfFull, curve = T)
      prroc_prHalfFull <- pr.curve(scores.class0 = fgHalfFull, scores.class1 = bgHalfFull, curve = T, rand.compute=T)
      aurocHalfFull = prroc_rocHalfFull$auc
      auprHalfFull = prroc_prHalfFull$auc.integral
      
      print("Perf for the model built on data from split 1 and tested on data from split 2...")
      # print(confusionMatrixSplit1)
      print(sprintf("AUROC & AUPR: %f.3 | %f.3", aurocS1S2, auprS1S2))
      
      print("Perf for the model built on data from split 2 and tested on data from split 1...")
      # print(confusionMatrixSplit2)
      print(sprintf("AUROC & AUPR: %f.3 | %f.3", aurocS2S1, auprS2S1))
      
      print("Perf for the model built on data from half of full data on holdout half...")
      # print(confusionMatrixSplit2)
      print(sprintf("AUROC & AUPR: %f.3 | %f.3", aurocHalfFull, auprHalfFull))
      cat('\n')
      
      perfTmp[[ii]] <- data.frame(disease_type = dt,
                            sample_type = st,
                            auroc_split12 = aurocS1S2,
                            aupr_split12 = auprS1S2,
                            auroc_split21 = aurocS2S1,
                            aupr_split21 = auprS2S1,
                            auroc_halffull = aurocHalfFull,
                            aupr_halffull = auprHalfFull)
      rm(mlModel_Split1)
      rm(mlModel_Split2)
      rm(mlModel_HalfFull)
    }
    perfTmp2[[kk]] <- do.call(rbind, perfTmp)
  }
  perf <- do.call(rbind, perfTmp2)
  perf$abbrev <- data.frame(abbrev=gsub("^TCGA-","",unique(metadataFull$investigation)), row.names=unique(metadataFull$disease_type))[perf$disease_type,"abbrev"]
  perf$sample_type <- factor(perf$sample_type, levels = c("Primary Tumor","Blood Derived Normal","Primary Tumor vs Solid Tissue Normal"))
  # Plot comparisons between split 1 and 2
  plot_auroc_split12_vs_split21 <- perf %>% distinct() %>%
    ggplot(aes(x = auroc_split21, y = auroc_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUROC split 2 on split 1", y = "AUROC split 1 on split 2", title = "Comparing AUROC on splits") + 
    ggsave(filename = paste0("auroc_",baseName,"_split12_vs_split21.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  plot_aupr_split12_vs_split21 <- perf %>% distinct() %>%
    ggplot(aes(x = aupr_split21, y = aupr_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUPR split 2 on split 1", y = "AUPR split 1 on split 2", title = "Comparing AUPR on splits") + 
    ggsave(filename = paste0("aupr_",baseName,"_split12_vs_split21.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  # Plot comparisons between split 1 and half of full
  plot_auroc_split12_vs_halffull <- perf %>% distinct() %>%
    ggplot(aes(x = auroc_halffull, y = auroc_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUROC half of full", y = "AUROC split 1", title = "Comparing AUROC on split 1 vs half of full") + 
    ggsave(filename = paste0("auroc_",baseName,"_split12_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  plot_aupr_split12_vs_halffull <- perf %>% distinct() %>%
    ggplot(aes(x = aupr_halffull, y = aupr_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUPR half of full", y = "AUPR split 1", title = "Comparing AUPR on split 1 vs half of full") + 
    ggsave(filename = paste0("aupr_",baseName,"_split12_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  # Plot comparisons between split 2 and half of full
  plot_auroc_split21_vs_halffull <- perf %>% distinct() %>%
    ggplot(aes(x = auroc_halffull, y = auroc_split21, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUROC half of full", y = "AUROC split 2", title = "Comparing AUROC on split 2 vs half of full") + 
    ggsave(filename = paste0("auroc_",baseName,"_split21_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  plot_aupr_split21_vs_halffull <- perf %>% distinct() %>%
    ggplot(aes(x = aupr_halffull, y = aupr_split21, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUPR half of full", y = "AUPR split 2", title = "Comparing AUPR on split 2 vs half of full") + 
    ggsave(filename = paste0("aupr_",baseName,"_split21_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  save(perf, file = paste0("perf_",baseName,".RData"))
  res <- list(
    perf=perf,
    plot_auroc_split12_vs_split21=plot_auroc_split12_vs_split21,
    plot_aupr_split12_vs_split21=plot_aupr_split12_vs_split21,
    plot_auroc_split12_vs_halffull=plot_auroc_split12_vs_halffull,
    plot_aupr_split12_vs_halffull=plot_aupr_split12_vs_halffull,
    plot_auroc_split21_vs_halffull=plot_auroc_split21_vs_halffull,
    plot_aupr_split21_vs_halffull=plot_aupr_split21_vs_halffull)
  return(res)
}

valResSplitsVsHalf <- crossTestSplitML(metadataSplit1 = valSplit1Metadata,
                 metadataSplit2 = valSplit2Metadata, 
                 metadataFull = metaQiitaCombined_Nonzero_DecontamV2,
                 vsnmDataSplit1 = valSplit1FungiCountDecontamV2Data_VSNM, 
                 vsnmDataSplit2 = valSplit2FungiCountDecontamV2Data_VSNM, 
                 vsnmDataFull = valFullFungiCountDecontamV2Data_VSNM,
                 sampleTypeList = c("Primary Tumor vs Solid Tissue Normal",
                                          "Blood Derived Normal",
                                          "Primary Tumor"),
                 numKFold = 4)
valResSplitsVsHalf %>% write.csv(file = "Interim_data/valResSplitsVsHalf_25Oct21.csv")
valResSplitsVsHalf$abbrev <- data.frame(abbrev=gsub("^TCGA-","",unique(metaQiitaCombined_Nonzero_DecontamV2$investigation)), row.names=unique(metaQiitaCombined_Nonzero_DecontamV2$disease_type))[valResSplitsVsHalf$disease_type,"abbrev"]
valResSplitsVsHalf$sample_type <- factor(valResSplitsVsHalf$sample_type, levels = c("Primary Tumor","Blood Derived Normal","Primary Tumor vs Solid Tissue Normal"))
save(valResSplitsVsHalf, file = "Interim_data/valResSplitsVsHalf_26Oct21.RData")

fit1 <- lm(aupr_split12 ~ aupr_split21, data = valResSplitsVsHalf)
summary(fit1)

require(ggrepel)
require(ggpmisc)
formula <- y ~ x  
valResSplitsVsHalf %>%
  ggplot(aes(x = aupr_split12, y = aupr_split21, label = abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
  facet_wrap(~sample_type) + 
  labs(x = "Perf half of full", y = "Perf split", title = "Comparing perf on split half vs half of full")

valResSplitsVsHalf_CT <- crossTestSplitML(metadataSplit1 = valSplit1Metadata,
                                       metadataSplit2 = valSplit2Metadata, 
                                       metadataFull = metaQiitaCombined_Nonzero_DecontamV2,
                                       vsnmDataSplit1 = valSplit1FungiCountDecontamV2Data_VSNM_CT, 
                                       vsnmDataSplit2 = valSplit2FungiCountDecontamV2Data_VSNM_CT, 
                                       vsnmDataFull = valFullFungiCountDecontamV2Data_VSNM_CT,
                                       sampleTypeList = c("Primary Tumor vs Solid Tissue Normal",
                                                          "Blood Derived Normal",
                                                          "Primary Tumor"),
                                       numKFold = 4)
valResSplitsVsHalf_CT %>% write.csv(file = "Interim_data/valResSplitsVsHalf_CT_25Oct21.csv")
valResSplitsVsHalf_CT$abbrev <- data.frame(abbrev=gsub("^TCGA-","",unique(metaQiitaCombined_Nonzero_DecontamV2$investigation)), row.names=unique(metaQiitaCombined_Nonzero_DecontamV2$disease_type))[valResSplitsVsHalf_CT$disease_type,"abbrev"]
valResSplitsVsHalf_CT$sample_type <- factor(valResSplitsVsHalf_CT$sample_type, levels = c("Primary Tumor","Blood Derived Normal","Primary Tumor vs Solid Tissue Normal"))
save(valResSplitsVsHalf_CT, file = "Interim_data/valResSplitsVsHalf_CT_26Oct21.RData")

valResSplitsVsHalf_CT %>%
  ggplot(aes(x = aupr_split12, y = aupr_split21, label = abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
  facet_wrap(~sample_type) + 
  labs(x = "Perf half of full", y = "Perf split", title = "Comparing perf on split half vs half of full")

save(valResSplitsVsHalf,
     valResSplitsVsHalf_CT,
     file = "Interim_data/valResSplitsVsHalf_and_CT_26Oct21.RData")

#--------------------------Barnacle runs--------------------------#

load("validation_2halves_tcga_vsnm_data_25Oct21.RData")

valResSplitsVsHalf_k10_v3 <- crossTestSplitML_Raw(metadataSplit1 = valSplit1Metadata,
                                       metadataSplit2 = valSplit2Metadata, 
                                       metadataFull = metaQiitaCombined_Nonzero_DecontamV2,
                                       vsnmDataSplit1 = valSplit1FungiCountDecontamV2Data_VSNM, 
                                       vsnmDataSplit2 = valSplit2FungiCountDecontamV2Data_VSNM, 
                                       vsnmDataFull = valFullFungiCountDecontamV2Data_VSNM,
                                       sampleTypeList = c("Primary Tumor vs Solid Tissue Normal",
                                                          "Blood Derived Normal",
                                                          "Primary Tumor"),
                                       numKFold = 10, baseName = "val_halves_k10_v3")

valResSplitsVsHalf_gridGBM <- crossTestSplitML(metadataSplit1 = valSplit1Metadata,
                                       metadataSplit2 = valSplit2Metadata, 
                                       metadataFull = metaQiitaCombined_Nonzero_DecontamV2,
                                       vsnmDataSplit1 = valSplit1FungiCountDecontamV2Data_VSNM, 
                                       vsnmDataSplit2 = valSplit2FungiCountDecontamV2Data_VSNM, 
                                       vsnmDataFull = valFullFungiCountDecontamV2Data_VSNM,
                                       sampleTypeList = c("Primary Tumor vs Solid Tissue Normal",
                                                          "Blood Derived Normal",
                                                          "Primary Tumor"),
                                       numKFold = 10, 
                                       caretTuneGrid = expand.grid(interaction.depth = seq(1,3),
                                                                   n.trees = floor((1:3) * 50),
                                                                   shrinkage = 0.1,
                                                                   n.minobsinnode = 5),
                                       baseName = "gridGBM_k10")

#-------------------Run on 2nd terminal-------------------#
crossTestSplitML_Raw <- function(metadataSplit1,
                             metadataSplit2,
                             metadataFull,
                             vsnmDataSplit1,
                             vsnmDataSplit2,
                             vsnmDataFull,
                             sampleTypeComparison = "Primary Tumor",
                             caretModel = "gbm",
                             samplingStrategy = "up",
                             numResampleIter = 1,
                             numKFold = 4,
                             baseName = "val_halves_k4",
                             caretTuneGrid = data.frame(n.trees=150, interaction.depth=3,
                                                        shrinkage=0.1,
                                                        n.minobsinnode=1),
                             sampleTypeList = c("Primary Tumor")){
  require(caret) # for model building
  require(gbm) # for machine learning
  require(PRROC) # for precision-recall curves
  require(MLmetrics) # for multiclass ML
  require(e1071)
  require(ggpubr)
  require(ggrepel)
  
  perf <- list()
  perfTmp <- list()
  perfTmp2 <- list()
  
  for(kk in seq_along(sampleTypeList)){
    st <- sampleTypeList[kk]
    print(st)
    
    if(st == "Primary Tumor vs Solid Tissue Normal"){
      metadataSplit1_FiltSt <- droplevels(metadataSplit1[metadataSplit1$sample_type %in%
                                                           c("Primary Tumor", "Solid Tissue Normal"),])
      metadataSplit2_FiltSt <- droplevels(metadataSplit2[metadataSplit2$sample_type %in%
                                                           c("Primary Tumor", "Solid Tissue Normal"),])
      metadataFull_FiltSt <- droplevels(metadataFull[metadataFull$sample_type %in%
                                                       c("Primary Tumor", "Solid Tissue Normal"),])
    } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
      metadataSplit1_FiltSt <- droplevels(metadataSplit1[metadataSplit1$sample_type %in% st,])
      metadataSplit2_FiltSt <- droplevels(metadataSplit2[metadataSplit2$sample_type %in% st,])
      metadataFull_FiltSt <- droplevels(metadataFull[metadataFull$sample_type %in% st,])
    }
    
    # print(seq_along(levels(metaTmp2$disease_type)))
    
    for(ii in seq_along(levels(metadataSplit1_FiltSt$disease_type))){
      
      dt <- levels(metadataSplit1_FiltSt$disease_type)[ii]
      print(dt)
      
      if(st == "Primary Tumor vs Solid Tissue Normal"){
        metadataSplit1_FiltStDz <- droplevels(metadataSplit1_FiltSt[metadataSplit1_FiltSt$disease_type == dt,])
        metadataSplit2_FiltStDz <- droplevels(metadataSplit2_FiltSt[metadataSplit2_FiltSt$disease_type == dt,])
        metadataFull_FiltStDz <- droplevels(metadataFull_FiltSt[metadataFull_FiltSt$disease_type == dt,])
        
        metadataSplit1_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataSplit1_FiltStDz$sample_type))
        metadataSplit2_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataSplit2_FiltStDz$sample_type))
        metadataFull_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataFull_FiltStDz$sample_type))
        
        positiveClass <- "PrimaryTumor"
        negativeClass <- "SolidTissueNormal"
      } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
        metadataSplit1_FiltStDz <- metadataSplit1_FiltSt
        metadataSplit2_FiltStDz <- metadataSplit2_FiltSt
        metadataFull_FiltStDz <- metadataFull_FiltSt
        
        metadataSplit1_FiltStDz$predY <- ifelse(metadataSplit1_FiltStDz$disease_type == dt, yes = dt, no = "OtherCancerType")
        metadataSplit2_FiltStDz$predY <- ifelse(metadataSplit2_FiltStDz$disease_type == dt, yes = dt, no = "OtherCancerType")
        metadataFull_FiltStDz$predY <- ifelse(metadataFull_FiltStDz$disease_type == dt, yes = dt, no = "OtherCancerType")
        
        metadataSplit1_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataSplit1_FiltStDz$predY))
        metadataSplit2_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataSplit2_FiltStDz$predY))
        metadataFull_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataFull_FiltStDz$predY))
        
        positiveClass <- gsub('([[:punct:]])|\\s+','',dt)
        negativeClass <- "OtherCancerType"
      }
      
      # print(table(metadataSplit1_FiltStDz$predY))
      # print(table(metadataSplit2_FiltStDz$predY))
      # print(table(metadataFull_FiltStDz$predY))
      
      # SKIP CANCER TYPES THAT ONLY HAVE ONE CLASS OF A COMPARISON
      if(length(table(metadataSplit1_FiltStDz$predY)) < 2 | length(table(metadataSplit2_FiltStDz$predY)) < 2){next}
      
      # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF 20 SAMPLES IN EITHER CLASS
      if(any(table(metadataSplit1_FiltStDz$predY) < 20) | any(table(metadataSplit2_FiltStDz$predY) < 20)){next}
      
      mlDataY_Split1 <- metadataSplit1_FiltStDz
      mlDataY_Split2 <- metadataSplit2_FiltStDz
      mlDataY_Full <- metadataFull_FiltStDz
      
      mlDataX_Split1 <- vsnmDataSplit1[rownames(mlDataY_Split1),]
      mlDataX_Split2 <- vsnmDataSplit2[rownames(mlDataY_Split2),]
      mlDataX_Full <- vsnmDataFull[rownames(mlDataY_Full),]
      
      set.seed(42)
      ctrl <- trainControl(method = "repeatedcv",
                           number = numKFold,
                           repeats = numResampleIter,
                           sampling = "up",
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE,
                           verboseIter = FALSE,
                           savePredictions = TRUE,
                           allowParallel=TRUE)
      
      #--------------------------Train split 1--------------------------#
      set.seed(42)
      trainX_Split1 <- mlDataX_Split1
      trainY_Split1 <- mlDataY_Split1[,"predY"]
      testX_Split1 <- mlDataX_Split2
      testY_Split1 <- mlDataY_Split2[,"predY"]
      # refactoredTrainY_Split1 <- factor(gsub('([[:punct:]])|\\s+','',trainY_Split1))
      
      print("Now training model on split 1...")
      set.seed(42)
      mlModel_Split1 <- train(x = trainX_Split1,
                              y = trainY_Split1,
                              method = "gbm",
                              # preProcess = c("scale","center"),
                              trControl = ctrl,
                              verbose = FALSE,
                              metric = "ROC",
                              tuneGrid = caretTuneGrid)
      
      #--------------------------Train split 2--------------------------#
      set.seed(42)
      trainX_Split2 <- mlDataX_Split2
      trainY_Split2 <- mlDataY_Split2[,"predY"]
      testX_Split2 <- mlDataX_Split1
      testY_Split2 <- mlDataY_Split1[,"predY"]
      # refactoredTrainY_Split2 <- factor(gsub('([[:punct:]])|\\s+','',trainY_Split2))
      
      print("Now training model on split 2...")
      set.seed(42)
      mlModel_Split2 <- train(x = trainX_Split2,
                              y = trainY_Split2,
                              method = "gbm",
                              # preProcess = c("scale","center"),
                              trControl = ctrl,
                              verbose = FALSE,
                              metric = "ROC",
                              tuneGrid = caretTuneGrid)
      
      #--------------------------Train half of full data--------------------------#
      
      set.seed(42) # random number seed -- important if you want comparable results between data types
      index <- createDataPartition(mlDataY_Full[,"predY"], p = 0.5, list = FALSE)
      trainHalfFullX <- mlDataX_Full[index,]
      trainHalfFullY <- mlDataY_Full[index,"predY"]
      testHalfFullX <- mlDataX_Full[-index,]
      testHalfFullY <- mlDataY_Full[-index,"predY"]
      
      print("Now training model on half of full data...")
      set.seed(42)
      mlModel_HalfFull <- train(x = trainHalfFullX,
                                y = trainHalfFullY,
                                method = "gbm",
                                # preProcess = c("scale","center"),
                                trControl = ctrl,
                                verbose = FALSE,
                                metric = "ROC",
                                tuneGrid = caretTuneGrid)
      
      #--------------------------Test on opposite split--------------------------#
      print("Now testing models on each other...")
      # Use model from split 1 to predict on split 2
      predClassS1S2 <- predict(mlModel_Split1, newdata = testX_Split1)
      predProbsS1S2 <- as.numeric(predict(mlModel_Split1, newdata = testX_Split1, type = "prob")[,positiveClass])
      # Use model from split 2 to predict on split 1
      predClassS2S1 <- predict(mlModel_Split2, newdata = testX_Split2)
      predProbsS2S1 <- as.numeric(predict(mlModel_Split2, newdata = testX_Split2, type = "prob")[,positiveClass])
      # Use model from half of full data to predict on holdout half
      predClassHalfFull <- predict(mlModel_HalfFull, newdata = testHalfFullX)
      predProbsHalfFull <- as.numeric(predict(mlModel_HalfFull, newdata = testHalfFullX, type = "prob")[,positiveClass])
      
      # Calculate perf of testing model from split 1 to predict on split 2
      fgS1S2 <- predProbsS1S2[testY_Split1 == positiveClass]
      bgS1S2 <- predProbsS1S2[testY_Split1 == negativeClass]
      prroc_rocS1S2 <- roc.curve(scores.class0 = fgS1S2, scores.class1 = bgS1S2, curve = T)
      prroc_prS1S2 <- pr.curve(scores.class0 = fgS1S2, scores.class1 = bgS1S2, curve = T, rand.compute=T)
      aurocS1S2 = prroc_rocS1S2$auc
      auprS1S2 = prroc_prS1S2$auc.integral
      # Calculate perf of testing model from split 2 to predict on split 1
      fgS2S1 <- predProbsS2S1[testY_Split2 == positiveClass]
      bgS2S1 <- predProbsS2S1[testY_Split2 == negativeClass]
      prroc_rocS2S1 <- roc.curve(scores.class0 = fgS2S1, scores.class1 = bgS2S1, curve = T)
      prroc_prS2S1 <- pr.curve(scores.class0 = fgS2S1, scores.class1 = bgS2S1, curve = T, rand.compute=T)
      aurocS2S1 = prroc_rocS2S1$auc
      auprS2S1 = prroc_prS2S1$auc.integral
      # Calculate perf of testing model from half of full data to predict on holdout half
      fgHalfFull <- predProbsHalfFull[testHalfFullY == positiveClass]
      bgHalfFull <- predProbsHalfFull[testHalfFullY == negativeClass]
      prroc_rocHalfFull <- roc.curve(scores.class0 = fgHalfFull, scores.class1 = bgHalfFull, curve = T)
      prroc_prHalfFull <- pr.curve(scores.class0 = fgHalfFull, scores.class1 = bgHalfFull, curve = T, rand.compute=T)
      aurocHalfFull = prroc_rocHalfFull$auc
      auprHalfFull = prroc_prHalfFull$auc.integral
      
      print("Perf for the model built on data from split 1 and tested on data from split 2...")
      # print(confusionMatrixSplit1)
      print(sprintf("AUROC & AUPR: %f.3 | %f.3", aurocS1S2, auprS1S2))
      
      print("Perf for the model built on data from split 2 and tested on data from split 1...")
      # print(confusionMatrixSplit2)
      print(sprintf("AUROC & AUPR: %f.3 | %f.3", aurocS2S1, auprS2S1))
      
      print("Perf for the model built on data from half of full data on holdout half...")
      # print(confusionMatrixSplit2)
      print(sprintf("AUROC & AUPR: %f.3 | %f.3", aurocHalfFull, auprHalfFull))
      cat('\n')
      
      perfTmp[[ii]] <- data.frame(disease_type = dt,
                                  sample_type = st,
                                  auroc_split12 = aurocS1S2,
                                  aupr_split12 = auprS1S2,
                                  auroc_split21 = aurocS2S1,
                                  aupr_split21 = auprS2S1,
                                  auroc_halffull = aurocHalfFull,
                                  aupr_halffull = auprHalfFull)
      rm(mlModel_Split1)
      rm(mlModel_Split2)
      rm(mlModel_HalfFull)
    }
    perfTmp2[[kk]] <- do.call(rbind, perfTmp)
  }
  perf <- do.call(rbind, perfTmp2)
  perf$abbrev <- data.frame(abbrev=gsub("^TCGA-","",unique(metadataFull$investigation)), row.names=unique(metadataFull$disease_type))[perf$disease_type,"abbrev"]
  perf$sample_type <- factor(perf$sample_type, levels = c("Primary Tumor","Blood Derived Normal","Primary Tumor vs Solid Tissue Normal"))
  # Plot comparisons between split 1 and 2
  plot_auroc_split12_vs_split21 <- perf %>% distinct() %>%
    ggplot(aes(x = auroc_split21, y = auroc_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUROC split 2 on split 1", y = "AUROC split 1 on split 2", title = "Comparing AUROC on splits") + 
    ggsave(filename = paste0("auroc_",baseName,"_split12_vs_split21.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  plot_aupr_split12_vs_split21 <- perf %>% distinct() %>%
    ggplot(aes(x = aupr_split21, y = aupr_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUPR split 2 on split 1", y = "AUPR split 1 on split 2", title = "Comparing AUPR on splits") + 
    ggsave(filename = paste0("aupr_",baseName,"_split12_vs_split21.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  # Plot comparisons between split 1 and half of full
  plot_auroc_split12_vs_halffull <- perf %>% distinct() %>%
    ggplot(aes(x = auroc_halffull, y = auroc_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUROC half of full", y = "AUROC split 1", title = "Comparing AUROC on split 1 vs half of full") + 
    ggsave(filename = paste0("auroc_",baseName,"_split12_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  plot_aupr_split12_vs_halffull <- perf %>% distinct() %>%
    ggplot(aes(x = aupr_halffull, y = aupr_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUPR half of full", y = "AUPR split 1", title = "Comparing AUPR on split 1 vs half of full") + 
    ggsave(filename = paste0("aupr_",baseName,"_split12_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  # Plot comparisons between split 2 and half of full
  plot_auroc_split21_vs_halffull <- perf %>% distinct() %>%
    ggplot(aes(x = auroc_halffull, y = auroc_split21, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUROC half of full", y = "AUROC split 2", title = "Comparing AUROC on split 2 vs half of full") + 
    ggsave(filename = paste0("auroc_",baseName,"_split21_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  plot_aupr_split21_vs_halffull <- perf %>% distinct() %>%
    ggplot(aes(x = aupr_halffull, y = aupr_split21, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUPR half of full", y = "AUPR split 2", title = "Comparing AUPR on split 2 vs half of full") + 
    ggsave(filename = paste0("aupr_",baseName,"_split21_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  save(perf, file = paste0("perf_",baseName,".RData"))
  res <- list(
    perf=perf,
    plot_auroc_split12_vs_split21=plot_auroc_split12_vs_split21,
    plot_aupr_split12_vs_split21=plot_aupr_split12_vs_split21,
    plot_auroc_split12_vs_halffull=plot_auroc_split12_vs_halffull,
    plot_aupr_split12_vs_halffull=plot_aupr_split12_vs_halffull,
    plot_auroc_split21_vs_halffull=plot_auroc_split21_vs_halffull,
    plot_aupr_split21_vs_halffull=plot_aupr_split21_vs_halffull)
  return(res)
}

valResSplitsVsHalf_rawData <- crossTestSplitML_Raw(metadataSplit1 = valSplit1Metadata,
                                          metadataSplit2 = valSplit2Metadata, 
                                          metadataFull = metaQiitaCombined_Nonzero_DecontamV2,
                                          vsnmDataSplit1 = valSplit1FungiCountDecontamV2Data, 
                                          vsnmDataSplit2 = valSplit2FungiCountDecontamV2Data, 
                                          vsnmDataFull = rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero,
                                          sampleTypeList = c("Primary Tumor vs Solid Tissue Normal",
                                                             "Blood Derived Normal",
                                                             "Primary Tumor"),
                                          numKFold = 10,
                                          baseName = "raw_data_k10")

valResSplitsVsHalf_CT <- crossTestSplitML(metadataSplit1 = valSplit1Metadata,
                                       metadataSplit2 = valSplit2Metadata, 
                                       metadataFull = metaQiitaCombined_Nonzero_DecontamV2,
                                       vsnmDataSplit1 = valSplit1FungiCountDecontamV2Data_VSNM_CT, 
                                       vsnmDataSplit2 = valSplit2FungiCountDecontamV2Data_VSNM_CT, 
                                       vsnmDataFull = valFullFungiCountDecontamV2Data_VSNM_CT,
                                       sampleTypeList = c("Primary Tumor vs Solid Tissue Normal",
                                                          "Blood Derived Normal",
                                                          "Primary Tumor"),
                                       numKFold = 10,
                                       baseName = "val_half_CT_data_k10")

valResSplitsVsHalf_k3 <- crossTestSplitML(metadataSplit1 = valSplit1Metadata,
                                               metadataSplit2 = valSplit2Metadata, 
                                               metadataFull = metaQiitaCombined_Nonzero_DecontamV2,
                                               vsnmDataSplit1 = valSplit1FungiCountDecontamV2Data_VSNM, 
                                               vsnmDataSplit2 = valSplit2FungiCountDecontamV2Data_VSNM, 
                                               vsnmDataFull = valFullFungiCountDecontamV2Data_VSNM,
                                               sampleTypeList = c("Primary Tumor vs Solid Tissue Normal",
                                                                  "Blood Derived Normal",
                                                                  "Primary Tumor"),
                                               numKFold = 3,
                                               baseName = "val_halves_k3")

#----------------------------------------------------------#
# Plot split performances
#----------------------------------------------------------#

plotSplitPerfs <- function(perfDf, baseName){
  # Plot comparisons between split 1 and 2
  plot_auroc_split12_vs_split21 <- perfDf %>% distinct() %>%
    ggplot(aes(x = auroc_split21, y = auroc_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUROC split 2 on split 1", y = "AUROC split 1 on split 2", title = "Comparing AUROC on splits") + 
    ggsave(filename = paste0("Figures/Supplementary_Figures/auroc_",baseName,"_split12_vs_split21.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  plot_aupr_split12_vs_split21 <- perfDf %>% distinct() %>%
    ggplot(aes(x = aupr_split21, y = aupr_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUPR split 2 on split 1", y = "AUPR split 1 on split 2", title = "Comparing AUPR on splits") + 
    ggsave(filename = paste0("Figures/Supplementary_Figures/aupr_",baseName,"_split12_vs_split21.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  # Plot comparisons between split 1 and half of full
  plot_auroc_split12_vs_halffull <- perfDf %>% distinct() %>%
    ggplot(aes(x = auroc_halffull, y = auroc_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUROC half of full", y = "AUROC split 1", title = "Comparing AUROC on split 1 vs half of full") + 
    ggsave(filename = paste0("Figures/Supplementary_Figures/auroc_",baseName,"_split12_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  plot_aupr_split12_vs_halffull <- perfDf %>% distinct() %>%
    ggplot(aes(x = aupr_halffull, y = aupr_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUPR half of full", y = "AUPR split 1", title = "Comparing AUPR on split 1 vs half of full") + 
    ggsave(filename = paste0("Figures/Supplementary_Figures/aupr_",baseName,"_split12_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  # Plot comparisons between split 2 and half of full
  plot_auroc_split21_vs_halffull <- perfDf %>% distinct() %>%
    ggplot(aes(x = auroc_halffull, y = auroc_split21, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUROC half of full", y = "AUROC split 2", title = "Comparing AUROC on split 2 vs half of full") + 
    ggsave(filename = paste0("Figures/Supplementary_Figures/auroc_",baseName,"_split21_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  plot_aupr_split21_vs_halffull <- perfDf %>% distinct() %>%
    ggplot(aes(x = aupr_halffull, y = aupr_split21, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUPR half of full", y = "AUPR split 2", title = "Comparing AUPR on split 2 vs half of full") + 
    ggsave(filename = paste0("Figures/Supplementary_Figures/aupr_",baseName,"_split21_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  
  res <- list(plot_auroc_split12_vs_split21=plot_auroc_split12_vs_split21,
              plot_aupr_split12_vs_split21=plot_aupr_split12_vs_split21,
              plot_auroc_split12_vs_halffull=plot_auroc_split12_vs_halffull,
              plot_aupr_split12_vs_halffull=plot_aupr_split12_vs_halffull,
              plot_auroc_split21_vs_halffull=plot_auroc_split21_vs_halffull,
              plot_aupr_split21_vs_halffull=plot_aupr_split21_vs_halffull)
  return(res)
}
rm(perf)
load("Interim_data/perf_raw_data_k10.RData", verbose = TRUE) # saved as perf object
perf_raw_data_k10 <- perf
rm(perf)
load("Interim_data/perf_val_halves_k10_v2.RData", verbose = TRUE)
perf_val_halves_k10 <- perf
rm(perf)
# Plot
plotSplitPerfs(perfDf = perf_raw_data_k10, baseName = "raw_data_k10")
plotSplitPerfs(perfDf = perf_val_halves_k10, baseName = "val_halves_k10")

#----------------------------------------------------------#
# Scramble labels and data: raw data
#----------------------------------------------------------#

load("Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy_with_coverage_filter_decontamV2_14Oct21.RData")
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

# save(metaQiitaCombined_Nonzero_DecontamV2_HMS_scrambled,
#      metaQiitaCombined_Nonzero_DecontamV2_BCM_scrambled,
#      metaQiitaCombined_Nonzero_DecontamV2_MDA_scrambled,
#      metaQiitaCombined_Nonzero_DecontamV2_WashU_scrambled,
#      metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_scrambled,
#      metaQiitaCombined_Nonzero_DecontamV2_UNC_scrambled,
#      metaQiitaCombined_Nonzero_DecontamV2_CMS_scrambled,
#      metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_scrambled,
#      # Scrambled meta normal data
#      rep200_HiSeq_Fungi_DecontamV2_HMS,
#      rep200_HiSeq_Fungi_DecontamV2_BCM,
#      rep200_HiSeq_Fungi_DecontamV2_MDA,
#      rep200_HiSeq_Fungi_DecontamV2_WashU,
#      rep200_HiSeq_Fungi_DecontamV2_Broad_WGS,
#      rep200_HiSeq_Fungi_DecontamV2_UNC,
#      rep200_HiSeq_Fungi_DecontamV2_CMS,
#      rep200_HiSeq_Fungi_DecontamV2_Broad_RNA,
#      # Shuffled data normal metadata
#      metaQiitaCombined_Nonzero_DecontamV2_HMS,
#      metaQiitaCombined_Nonzero_DecontamV2_BCM,
#      metaQiitaCombined_Nonzero_DecontamV2_MDA,
#      metaQiitaCombined_Nonzero_DecontamV2_WashU,
#      metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS,
#      metaQiitaCombined_Nonzero_DecontamV2_UNC,
#      metaQiitaCombined_Nonzero_DecontamV2_CMS,
#      metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA,
#      
#      rep200_HiSeq_Fungi_DecontamV2_HMS_shuffled,
#      rep200_HiSeq_Fungi_DecontamV2_BCM_shuffled,
#      rep200_HiSeq_Fungi_DecontamV2_MDA_shuffled,
#      rep200_HiSeq_Fungi_DecontamV2_WashU_shuffled,
#      rep200_HiSeq_Fungi_DecontamV2_Broad_WGS_shuffled,
#      rep200_HiSeq_Fungi_DecontamV2_UNC_shuffled,
#      rep200_HiSeq_Fungi_DecontamV2_CMS_shuffled,
#      rep200_HiSeq_Fungi_DecontamV2_Broad_RNA_shuffled,
#      file = "Interim_data/shuffled_controls_raw_data_TCGA_26Oct21.RData")
#----------------------------------------------------------#
# Scramble labels and data: VSNM data
#----------------------------------------------------------#

load("Interim_data/snmDataFungi_DecontamV2_13Oct21.RData")
load("Interim_data/validation_2halves_tcga_vsnm_data_25Oct21.RData")

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

# save(metaQiitaCombined_Nonzero_DecontamV2_scrambled,
#      valFullFungiCountDecontamV2Data_VSNM,
#      metaQiitaCombined_Nonzero_DecontamV2,
#      valFullFungiCountDecontamV2Data_VSNM_Shuffled,
#      file = "Interim_data/shuffled_controls_VSNM_data_TCGA_26Oct21.RData")

#----------------------------------------------------------#
# Plot scrambled ML perf - version 1 with global scrambling and shuffling
#----------------------------------------------------------#

source("Supporting_scripts/S05-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_scrambled_and_shuffled_controls_ALL_DecontamV2_26Oct21.csv", stringsAsFactors = FALSE)
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

mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "rep200_HiSeq_Fungi_DecontamV2_HMS"] <- "HMS metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "rep200_HiSeq_Fungi_DecontamV2_HMS_shuffled"] <- "HMS count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "rep200_HiSeq_Fungi_DecontamV2_MDA"] <- "MDA metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "rep200_HiSeq_Fungi_DecontamV2_MDA_shuffled"] <- "MDA count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "rep200_HiSeq_Fungi_DecontamV2_BCM"] <- "BCM metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "rep200_HiSeq_Fungi_DecontamV2_BCM_shuffled"] <- "BCM count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "rep200_HiSeq_Fungi_DecontamV2_WashU"] <- "WashU metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "rep200_HiSeq_Fungi_DecontamV2_WashU_shuffled"] <- "WashU count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "rep200_HiSeq_Fungi_DecontamV2_Broad_WGS"] <- "Broad metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "rep200_HiSeq_Fungi_DecontamV2_Broad_WGS_shuffled"] <- "Broad count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "rep200_HiSeq_Fungi_DecontamV2_UNC"] <- "UNC metadata labels scrambled (RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "rep200_HiSeq_Fungi_DecontamV2_UNC_shuffled"] <- "UNC count data shuffled (RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "rep200_HiSeq_Fungi_DecontamV2_CMS"] <- "CMS metadata labels scrambled (RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "rep200_HiSeq_Fungi_DecontamV2_CMS_shuffled"] <- "CMS count data shuffled (RNA-Seq)"
## NOTE: Broad_RNA only included GBM tumors, so no ML comparisons were made
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "valFullFungiCountDecontamV2Data_VSNM"] <- "VSNM all metadata labels scrambled (WGS+RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM$datasetName == "valFullFungiCountDecontamV2Data_VSNM_Shuffled"] <- "VSNM all count data shuffled (WGS+RNA-Seq)"

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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_HMS_PT_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_MDA_PT_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_BCM_PT_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_WashU_PT_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_Broad_WGS_PT_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_UNC_PT_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_CMS_PT_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_TCGA_VSNM_PT_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_HMS_PT_vs_NAT_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_BCM_PT_vs_NAT_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_UNC_PT_vs_NAT_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_CMS_PT_vs_NAT_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_TCGA_VSNM_PT_vs_NAT_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_HMS_BDN_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_MDA_BDN_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_BCM_BDN_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_WashU_BDN_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_Broad_WGS_BDN_species_decontamV2.jpeg", dpi = "retina",
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
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/control_scrambled_mlPerfAll10k_rep1_TCGA_VSNM_BDN_species_decontamV2.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")




#----------------------------------------------------------#
# Plot scrambled ML perf - version 2 with dynamic scrambling and shuffling just prior to ML
# NOTE: This is the preferred version
#----------------------------------------------------------#


source("Supporting_scripts/S05-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_scrambled_and_shuffled_controls_ALL_DecontamV2_17Nov21.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 <- mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2[,!(colnames(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2) == "X")]
colnames(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$nullAUPR <- ifelse(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$minorityClassName == "SolidTissueNormal",
                                                             yes=mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$majorityClassSize/(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$minorityClassSize+mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$majorityClassSize),
                                                             no=mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$minorityClassSize/(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$minorityClassSize+mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$majorityClassSize))
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
table(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName)
table(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName)
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_HMS_scrambled"] <- "HMS metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_HMS_shuffled"] <- "HMS count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_MDA_scrambled"] <- "MDA metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_MDA_shuffled"] <- "MDA count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_BCM_scrambled"] <- "BCM metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_BCM_shuffled"] <- "BCM count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_WashU_scrambled"] <- "WashU metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_WashU_shuffled"] <- "WashU count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_scrambled"] <- "Broad metadata labels scrambled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_shuffled"] <- "Broad count data shuffled (WGS)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_UNC_scrambled"] <- "UNC metadata labels scrambled (RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_UNC_shuffled"] <- "UNC count data shuffled (RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_CMS_scrambled"] <- "CMS metadata labels scrambled (RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_CMS_shuffled"] <- "CMS count data shuffled (RNA-Seq)"
## NOTE: Broad_RNA only included GBM tumors, so no ML comparisons were made
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_scrambled"] <- "VSNM all metadata labels scrambled (WGS+RNA-Seq)"
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName[mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$metadataName == "metaQiitaCombined_Nonzero_DecontamV2_shuffled"] <- "VSNM all count data shuffled (WGS+RNA-Seq)"

mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName <- factor(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName,
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
table(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2$datasetName)
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
# HMS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_HMS_PT_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")
# MDA
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_MDA_PT_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")
# BCM
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_BCM_PT_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")
# WashU
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_WashU_PT_species_decontamV2.svg", dpi = "retina",
         width = 5, height = 4, units = "in")
# Broad WGS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_Broad_WGS_PT_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")
# UNC
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_UNC_PT_species_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# CMS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_CMS_PT_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")
# Full dataset
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_TCGA_VSNM_PT_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#
# HMS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_HMS_PT_vs_NAT_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")
# MDA does not have sufficient PT and NAT samples to compare
# BCM
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_BCM_PT_vs_NAT_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")
# WashU does not have sufficient PT and NAT samples to compare
# Broad WGS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")
# UNC
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_UNC_PT_vs_NAT_species_decontamV2.svg", dpi = "retina",
         width = 8, height = 4, units = "in")
# CMS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_CMS_PT_vs_NAT_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")

# Full dataset
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_TCGA_VSNM_PT_vs_NAT_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")

#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
# HMS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_HMS_BDN_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")
# MDA
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_MDA_BDN_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")
# BCM
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_BCM_BDN_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")
# WashU
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_WashU_BDN_species_decontamV2.svg", dpi = "retina",
         width = 5, height = 4, units = "in")
# Broad WGS
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_Broad_WGS_BDN_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")
# Full dataset
mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2 %>%
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_TCGA_VSNM_BDN_species_decontamV2.svg", dpi = "retina",
         width = 6, height = 4, units = "in")








#----------------------------------------------------------#
# Overlay actual and (scrambled) control performance per seq center
#----------------------------------------------------------#

load("Interim_data/mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann.RData") # mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann
load("Interim_data/mlPerfAll10k_Allcancer_Raw.RData") # mlPerfAll10k_Allcancer_Raw

source("Supporting_scripts/S05-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Raw_Control_Overlay <- rbind(mlPerfAll10k_Allcancer_Scrambled_Raw_VSNM_v2,
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_HMS_PT_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_BCM_PT_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_MDA_PT_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_WashU_PT_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_Broad_WGS_PT_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_UNC_PT_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_CMS_PT_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_HMS_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_BCM_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_UNC_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_CMS_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_HMS_BDN_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_BCM_BDN_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_MDA_BDN_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_WashU_BDN_species_overlay_decontamV2.svg", dpi = "retina",
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
ggsave("Figures/Supplementary_Figures/control_v2_scrambled_mlPerfAll10k_rep1_Broad_WGS_BDN_species_overlay_decontamV2.svg", dpi = "retina",
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

#----------------------------------
# Testing rstatix

require(rstatix)
## AUROC
mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Broad",datasetNameShort)) %>%
  distinct() %>% droplevels() %>%
  wilcox_test(AUROC ~ datasetNameShort, exact = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(p.adj = signif(p.adj, digits=3)) %>%
  add_xy_position(x = "datasetNameShort", step.increase = 0.25) -> roc.stat.test
data.frame(roc.stat.test)

mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetNameShort)) %>%
  distinct() %>% droplevels() %>%
  ggboxplot(x = "datasetNameShort",
            y = "AUROC",
            notch = TRUE,
            add = "jitter",
            add.params = list(alpha=0.4),
            legend = "none",
            fill = "datasetNameShort",
            xlab = "Data type",
            palette = "nejm") +
  rotate_x_text(30) + 
  stat_pvalue_manual(roc.stat.test, 
                     # y.position = c(1.03, 1.12, 1.06),
                     # y.position = c(1.03, 1.10,1.17, 1.24, 1.31, 1.45),
                     label = "q = {p.adj}", 
                     size = 3) +
  geom_hline(yintercept = 0.5, linetype="dotted") + 
  labs(fill = "datasetNameShort") + theme(plot.title = element_text(hjust=0.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), 
                     limits = c(0,1.01*max(roc.stat.test$y.position))) #-> aggregatedTumorROC

mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetNameShort)) %>%
  distinct() %>% droplevels() %>%
  wilcox_test(AUPR ~ datasetNameShort, exact = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(p.adj = signif(p.adj, digits=3)) %>%
  add_xy_position(x = "datasetNameShort", step.increase = 0.25/2) -> pr.stat.test
data.frame(pr.stat.test)

mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetNameShort)) %>%
  distinct() %>% droplevels() %>%
  ggboxplot(x = "datasetNameShort",
            y = "AUPR",
            notch = TRUE,
            add = "jitter",
            add.params = list(alpha=0.4),
            legend = "none",
            fill = "datasetNameShort",
            xlab = "Data type",
            palette = "nejm") +
  rotate_x_text(30) + 
  stat_pvalue_manual(pr.stat.test, 
                     # y.position = c(1.03, 1.12, 1.06),
                     # y.position = c(1.03, 1.10,1.17, 1.24, 1.31, 1.45),
                     label = "q = {p.adj}", 
                     size = 3) +
  labs(fill = "datasetNameShort") + theme(plot.title = element_text(hjust=0.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), 
                     limits = c(0,1.01*max(pr.stat.test$y.position))) #-> aggregatedTumorROC
  


#----------------------------------

## Write plot function
plotControlsRaw <- function(seqCenter="HMS", 
                            inputSampleType="Primary Tumor",
                            qvalSize = 2.5,
                            statSpacingROC=0.5,
                            statSpacingPR=0.25,
                            alphaVal=0.2,
                            tipLength=0.01,
                            qvalAsterisks=FALSE){
  require(rstatix)
  if(inputSampleType=="Primary Tumor"){st = "PT"}
  if(inputSampleType=="Primary Tumor vs Solid Tissue Normal"){st = "PT_vs_NAT"}
  if(inputSampleType=="Blood Derived Normal"){st = "BDN"}
  ## AUROC
  mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short %>%
    filter(sampleType == inputSampleType) %>%
    filter(grepl(seqCenter,datasetNameShort)) %>%
    distinct() %>% droplevels() %>%
    wilcox_test(AUROC ~ datasetNameShort, exact = TRUE) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    mutate(p.adj = signif(p.adj, digits=3)) %>%
    add_xy_position(x = "datasetNameShort", step.increase = statSpacingROC) -> roc.stat.test
    
  if(all(roc.stat.test$y.position == 1)){
    roc.stat.test$y.position <- seq(from=1.05, by=0.2, length.out = length(roc.stat.test$y.position))
  }
  # print(data.frame(roc.stat.test))
  
  mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short %>%
    filter(sampleType == inputSampleType) %>%
    filter(grepl(seqCenter,datasetNameShort)) %>%
    distinct() %>% droplevels() %>%
    ggboxplot(x = "datasetNameShort",
             y = "AUROC", fill = "datasetNameShort",
             legend = "none",
             notch = TRUE,
             xlab = "",
             add = "jitter",
             add.params = list(alpha=alphaVal),
             palette = "nejm") +
    rotate_x_text(90) + 
    stat_pvalue_manual(roc.stat.test, 
                       # label = "q = {p.adj}", 
                       label = ifelse(qvalAsterisks, 
                                      yes = "{p.adj.signif}", 
                                      no = "q = {p.adj}"),
                       tip.length = tipLength,
                       size = qvalSize) +
    geom_hline(yintercept = 0.5, linetype="dotted") + 
    labs(fill = "datasetNameShort") + 
    theme(plot.title = element_text(hjust=0.5)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), 
                       limits = c(0,1.01*max(roc.stat.test$y.position))) -> plotROC
  
  ## AUPR
  mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short %>%
    filter(sampleType == inputSampleType) %>%
    filter(grepl(seqCenter,datasetNameShort)) %>%
    distinct() %>% droplevels() %>%
    wilcox_test(AUPR ~ datasetNameShort, exact = TRUE) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    mutate(p.adj = signif(p.adj, digits=3)) %>%
    add_xy_position(x = "datasetNameShort", step.increase = statSpacingPR) -> pr.stat.test
  
  if(all(pr.stat.test$y.position == 1)){
    pr.stat.test$y.position <- seq(from=1.05, by=0.2, length.out = length(pr.stat.test$y.position))
  }
  # print(data.frame(pr.stat.test))
  
  mlPerfAll10k_Allcancer_Raw_Control_Overlay_Short %>%
    filter(sampleType == inputSampleType) %>%
    filter(grepl(seqCenter,datasetNameShort)) %>%
    distinct() %>% droplevels() %>%
    ggboxplot(x = "datasetNameShort",
             y = "AUPR", fill = "datasetNameShort",
             legend = "none",
             notch = TRUE,
             add = "jitter",
             add.params = list(alpha=alphaVal),
             xlab = "",
             palette = "nejm") +
    rotate_x_text(90) + 
    stat_pvalue_manual(pr.stat.test, 
                       # label = "q = {p.adj}", 
                       label = "{p.adj.signif}", 
                       tip.length = tipLength,
                       size = qvalSize) +
    labs(fill = "datasetNameShort") + 
    theme(plot.title = element_text(hjust=0.5)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), 
                       limits = c(0,1.01*max(pr.stat.test$y.position))) -> plotPR
  
  combinedPlot <- ggarrange(plotROC, plotPR, ncol = 2) 
  combinedPlotAnnotated <- annotate_figure(combinedPlot, top = text_grob(paste0(seqCenter," | ",inputSampleType,"\nActual vs. Controls"), 
                                        color = "black", face = "bold", size = 14))
  print(combinedPlotAnnotated)
  baseName <- ifelse(qvalAsterisks, yes = "control_v2_signif_scrambled_overlay_raw_data",
                     no = "control_v2_scrambled_overlay_raw_data_")
  ggsave(filename = paste0("Figures/Supplementary_Figures/",baseName,"_",seqCenter,"_",st,".svg"),
           plot = combinedPlotAnnotated,
           dpi = "retina", units = "in", width = 6, height = 7)
}
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

load("Interim_data/mlPerfAll10k_Allcancer_Shared.RData")
load("Interim_data/mlPerfAll10k_Allcancer.RData")
load("Interim_data/mlPerfAll10k_Allcancer_Cov.RData")

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

## Write plot function
plotControlsVSNM <- function(inputSampleType="Primary Tumor"){
  if(inputSampleType=="Primary Tumor"){st = "PT"}
  if(inputSampleType=="Primary Tumor vs Solid Tissue Normal"){st = "PT_vs_NAT"}
  if(inputSampleType=="Blood Derived Normal"){st = "BDN"}
  mlPerfAll10k_Allcancer_Overlay_VSNM %>%
    filter(sampleType == inputSampleType) %>%
    distinct() %>% droplevels() %>%
    ggboxplot(x = "abbrev",
              y = "AUROC",
              fill = "datasetName",
              xlab = "Cancer type",
              palette = "nejm") +
    rotate_x_text(90) +
    geom_hline(yintercept = 0.5, linetype="dotted") + 
    labs(fill = "Dataset") + theme(plot.title = element_text(hjust=0.5)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    stat_compare_means(aes(group = statGroups), label = "p.signif") -> plotROC
  # NOTE: Default for grouped variables of stat_compare_means() is Wilcoxon
  
  mlPerfAll10k_Allcancer_Overlay_VSNM %>%
    filter(sampleType == inputSampleType) %>%
    distinct() %>% droplevels() %>%
    ggboxplot(x = "abbrev",
              y = "AUPR",
              fill = "datasetName",
              xlab = "Cancer type",
              legend = "none",
              palette = "nejm") +
    rotate_x_text(90) +
    labs(fill = "Dataset") + theme(plot.title = element_text(hjust=0.5)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    stat_compare_means(aes(group = statGroups), label = "p.signif") -> plotPR
  # NOTE: Default for grouped variables of stat_compare_means() is Wilcoxon
  
  combinedPlotTitle <- paste0("TCGA ML Performance: Actual vs. Control (scrambled metadata or shuffled samples)\n",
                              inputSampleType,"\n(per-cancer p-values derived from grouped actual vs. control performances)")
  combinedPlot <- ggarrange(plotROC, plotPR, nrow = 2) 
  combinedPlotAnnotated <- annotate_figure(combinedPlot, 
                                           top = text_grob(combinedPlotTitle, 
                                                           color = "black", face = "bold", size = 14))
  print(combinedPlotAnnotated)
  ggsave(filename = paste0("Figures/Supplementary_Figures/control_v2_scrambled_overlay_raw_data_TCGA_VSNM_",st,".svg"),
         plot = combinedPlotAnnotated,
         dpi = "retina", units = "in", width = 14, height = 10)
}
#-------------------------Plot VSNM overlays-------------------------#
plotControlsVSNM(inputSampleType = "Primary Tumor")
plotControlsVSNM(inputSampleType = "Primary Tumor vs Solid Tissue Normal")
plotControlsVSNM(inputSampleType = "Blood Derived Normal")
