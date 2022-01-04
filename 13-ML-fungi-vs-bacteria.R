#-----------------------------------------------------------------------------
# 13-ML-fungi-vs-bacteria.R
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

# save(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts, 
#      rep200Data_WGS_RNA_Matched,
#      rep200Data_WGS_RNA_Matched_Bacteria,
#      rep200Data_WGS_RNA_Matched_Fungi,
#      file = "Interim_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_29Sep21.RData")
load("Interim_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_29Sep21.RData")

# snmDataOGUFungiDecontamV2,
# vdge_dataE_DecontamV2,
# rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero,
# metaQiitaCombined_Nonzero_DecontamV2,
load("Interim_data/snmDataFungi_DecontamV2_13Oct21.RData")
rep200TaxSplit <- read.csv("Supporting_data/rep200_lineage_map_split.csv", stringsAsFactors = FALSE, row.names = 1)

# - Develop multi-class ML code
# - Subsample # species for 1000 iterations. n=seq(25, 300, by = 25)
# - Rarefy # of reads for 10 iterations each time. n=100, 500, 1000, 5000, 10e3, 15e3, 20e3, 25e3, 30e3, 35e3

metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_WGS <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>% filter(experimental_strategy == "WGS") %>% droplevels()
rep200Data_WGS_RNA_Matched_Fungi_WGS <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_WGS),]
rep200Data_WGS_RNA_Matched_Bacteria_WGS <- rep200Data_WGS_RNA_Matched_Bacteria[rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_WGS),]

quantile(rowSums(rep200Data_WGS_RNA_Matched_Fungi_WGS), probs = seq(0.05, 1, by = 0.10))
# 5%      15%      25%      35%      45%      55%      65%      75%      85%      95% 
#   1131.25  3047.25  4852.75  6726.25  9080.00 11851.25 14970.50 18318.50 23310.50 39521.00
quantile(rowSums(rep200Data_WGS_RNA_Matched_Bacteria_WGS), probs = seq(0.05, 1, by = 0.10))
# 5%       15%       25%       35%       45%       55%       65%       75%       85%       95% 
#   7652.25  22999.25  36441.75  49789.00  65063.00  84048.75 109104.50 146196.00 219860.25 851346.25

#----------------------------------------------------------#
# Subset data
#----------------------------------------------------------#

# Subset to WGS and remove 1 fungi sample with 0 counts: 13722.58cfa82ee4b0c9d6adf6a982
metaFungiVsBacteria <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>%
  filter(cgc_platform == "Illumina HiSeq", experimental_strategy == "WGS") %>% 
  filter(sample_name != "58cfa82ee4b0c9d6adf6a982") %>% droplevels()
countsFungiWGS <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaFungiVsBacteria),]
countsBacteriaWGS <- rep200Data_WGS_RNA_Matched_Bacteria[rownames(metaFungiVsBacteria),]

save(metaFungiVsBacteria,
     countsFungiWGS,
     countsBacteriaWGS,
     file = "Interim_data/data_for_ml_fungi_vs_bacteria_numFeat_28Oct21.RData")

load("Interim_data/data_for_ml_fungi_vs_bacteria_numFeat_28Oct21.RData")

#----------------------------------------------------------#
# Generate lists of randomly sampled OGU features
#----------------------------------------------------------#

# Fungi features
fungi10 <- replicate(1e3, expr = sample(1:ncol(countsFungiWGS), size = 10, replace = FALSE))
fungi25 <- replicate(1e3, expr = sample(1:ncol(countsFungiWGS), size = 25, replace = FALSE))
fungi50 <- replicate(1e3, expr = sample(1:ncol(countsFungiWGS), size = 50, replace = FALSE))
fungi75 <- replicate(1e3, expr = sample(1:ncol(countsFungiWGS), size = 75, replace = FALSE))
fungi100 <- replicate(1e3, expr = sample(1:ncol(countsFungiWGS), size = 100, replace = FALSE))
fungi125 <- replicate(1e3, expr = sample(1:ncol(countsFungiWGS), size = 125, replace = FALSE))
fungi150 <- replicate(1e3, expr = sample(1:ncol(countsFungiWGS), size = 150, replace = FALSE))
fungi175 <- replicate(1e3, expr = sample(1:ncol(countsFungiWGS), size = 175, replace = FALSE))
fungi200 <- replicate(1e3, expr = sample(1:ncol(countsFungiWGS), size = 200, replace = FALSE))
fungi225 <- replicate(1e3, expr = sample(1:ncol(countsFungiWGS), size = 225, replace = FALSE))
fungi250 <- replicate(1e3, expr = sample(1:ncol(countsFungiWGS), size = 250, replace = FALSE))
fungi275 <- replicate(1e3, expr = sample(1:ncol(countsFungiWGS), size = 275, replace = FALSE))
fungi300 <- replicate(1e3, expr = sample(1:ncol(countsFungiWGS), size = 300, replace = FALSE))
fungiXList <- list(fungi10, fungi25, fungi50, fungi75, fungi100, fungi125, fungi150,
                   fungi175, fungi200, fungi225, fungi250, fungi275, fungi300)
# Bacteria features
bacteria10 <- replicate(1e3, expr = sample(1:ncol(countsBacteriaWGS), size = 10, replace = FALSE))
bacteria25 <- replicate(1e3, expr = sample(1:ncol(countsBacteriaWGS), size = 25, replace = FALSE))
bacteria50 <- replicate(1e3, expr = sample(1:ncol(countsBacteriaWGS), size = 50, replace = FALSE))
bacteria75 <- replicate(1e3, expr = sample(1:ncol(countsBacteriaWGS), size = 75, replace = FALSE))
bacteria100 <- replicate(1e3, expr = sample(1:ncol(countsBacteriaWGS), size = 100, replace = FALSE))
bacteria125 <- replicate(1e3, expr = sample(1:ncol(countsBacteriaWGS), size = 125, replace = FALSE))
bacteria150 <- replicate(1e3, expr = sample(1:ncol(countsBacteriaWGS), size = 150, replace = FALSE))
bacteria175 <- replicate(1e3, expr = sample(1:ncol(countsBacteriaWGS), size = 175, replace = FALSE))
bacteria200 <- replicate(1e3, expr = sample(1:ncol(countsBacteriaWGS), size = 200, replace = FALSE))
bacteria225 <- replicate(1e3, expr = sample(1:ncol(countsBacteriaWGS), size = 225, replace = FALSE))
bacteria250 <- replicate(1e3, expr = sample(1:ncol(countsBacteriaWGS), size = 250, replace = FALSE))
bacteria275 <- replicate(1e3, expr = sample(1:ncol(countsBacteriaWGS), size = 275, replace = FALSE))
bacteria300 <- replicate(1e3, expr = sample(1:ncol(countsBacteriaWGS), size = 300, replace = FALSE))
bacteriaXList <- list(bacteria10, bacteria25, bacteria50, bacteria75, bacteria100, bacteria125, bacteria150,
                   bacteria175, bacteria200, bacteria225, bacteria250, bacteria275, bacteria300)

#----------------------------------------------------------#
# Identify feature numbers for WIS intersected fungi
#----------------------------------------------------------#

load("Interim_data/decontamResultsV2_13Oct21.RData")

ogusWISintersect <- decontamResultsV2 %>% filter(shared_with_WIS == "YES") %>% rownames()
indexWISfungi <- which(colnames(countsFungiWGS) %in% ogusWISintersect)
inverseIndexNonWISfungi <- which(!(1:319) %in% indexWISfungi)

save(indexWISfungi, inverseIndexNonWISfungi, 
     file = "Interim_data/indices_WIS_and_nonWIS_fungi_for_simulations.RData")

#----------------------------------------------------------#
# Subset feature sets to WIS intersected fungi, bacteria, fungi+bacteria
#----------------------------------------------------------#
load("Interim_data/decontamResultsV2_13Oct21.RData")
load("Interim_data/shared_bacterial_features_at_each_taxa_level_15Oct21.RData")

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
     file = "Interim_data/data_for_ml_tcga_synergy_fungi_and_bacteria_18Nov21.RData")


#------------------------------Plot results------------------------------#

plotSynergyPerf <- function(rDataFilePath, 
                             inputSampleType="Primary Tumor",
                             numIter = 100,
                             modelType = "xgbTree"){
  require(rstatix)
  if(inputSampleType=="Primary Tumor"){st = "PT"}
  if(inputSampleType=="Primary Tumor vs Solid Tissue Normal"){st = "PT_vs_NAT"}
  if(inputSampleType=="Blood Derived Normal"){st = "BDN"}
  if(modelType=="rf"){modelTypeFormatted = "Random Forest"}
  if(modelType %in% c("gbm","xgbTree")){modelTypeFormatted = "Gradient Boosting"}
  seqCenterFormatted <- "AllSeqCenters"
  seqCenterPlotTitle <- "All WGS Sequencing Centers"
  
  load(rDataFilePath) # loads perfCombinedAll object
  keepCols <- c("AUC","prAUC")
  perfCombinedAllFormatted <- reshape2::melt(perfCombinedAll, id.vars = c("iter")) %>%
    mutate(Dataset = factor(case_when(grepl("^f_",variable) ~ "Fungi",
                                      grepl("^b_",variable) ~ "Bacteria",
                                      grepl("^all_",variable) ~ "Fungi+Bacteria"),
                            levels = c("Fungi","Bacteria","Fungi+Bacteria"))) %>%
    mutate(variable = gsub("^f_|^b_|all_","",variable)) %>% filter(variable %in% keepCols)
  
  perfCombinedAllFormatted %>%
    filter(variable == "AUC") %>%
    wilcox_test(value ~ Dataset, exact = TRUE) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance() %>%
    add_xy_position() -> roc.stat

  perfCombinedAllFormatted %>%
    filter(variable == "AUC") %>%
    ggboxplot(x = "Dataset",
              y = "value",
              palette = "nejm",
              fill = "Dataset",
              legend = "none",
              add = "jitter",
              add.params = list(alpha=0.4),
              notch = TRUE,
              xlab = "Feature set ∩ WIS",
              ylab = "Average pan-cancer AUROC") +
    stat_pvalue_manual(roc.stat,
                       label = "q = {p.adj}") -> plotROC
  
  perfCombinedAllFormatted %>%
    filter(variable == "prAUC") %>%
    wilcox_test(value ~ Dataset, exact = TRUE) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance() %>%
    add_xy_position() -> pr.stat
  
  perfCombinedAllFormatted %>%
    filter(variable == "prAUC") %>%
    ggboxplot(x = "Dataset",
              y = "value",
              palette = "nejm",
              fill = "Dataset",
              legend = "none",
              add = "jitter",
              add.params = list(alpha=0.4),
              notch = TRUE,
              xlab = "Feature set ∩ WIS",
              ylab = "Average pan-cancer AUPR") +
    stat_pvalue_manual(pr.stat,
                       label = "q = {p.adj}") -> plotPR
  
  combinedPlotTitle <- paste0("TCGA ML: WIS-intersecting fungi vs. bacteria vs. fungi+bacteria\n",
                              inputSampleType," | ",seqCenterPlotTitle," | ",modelTypeFormatted,
                              "\n(",numIter," iterations of multi-class classification per box)")
  combinedPlot <- ggarrange(plotROC, plotPR, ncol = 2)
  combinedPlotAnnotated <- annotate_figure(combinedPlot,
                                           top = text_grob(combinedPlotTitle,
                                                           color = "black", face = "bold", size = 14))
  print(combinedPlotAnnotated)
  ggsave(filename = paste0("Figures/Supplementary_Figures/synergy_tcga_fungi_bacteria_",
                           seqCenterFormatted,"_",st,"_numIter",numIter,"_",modelType,".svg"),
         plot = combinedPlotAnnotated,
         dpi = "retina", units = "in", width = 8, height = 5)
  rm(perfCombinedAll)
  # return(perfCombinedAllFormatted)
}
# Primary Tumor
plotSynergyPerf(rDataFilePath = "Interim_data/tcga_synergy_fungi_bacteria_PrimaryTumor_AllSeqCenters_numIter100_k4_modelType_xgbTree.RData", 
                modelType = "xgbTree", numIter = 100, inputSampleType = "Primary Tumor")
# Blood Derived Normal
plotSynergyPerf(rDataFilePath = "Interim_data/tcga_synergy_fungi_bacteria_BloodDerivedNormal_AllSeqCenters_numIter100_k4_modelType_xgbTree.RData", 
                modelType = "xgbTree", numIter = 100, inputSampleType = "Blood Derived Normal")

#----------------------------------------------------------#
# Build multi-class ML function by sampling same number of 
# features as WIS intersect
#----------------------------------------------------------#
require(xgboost) # for machine learning

load("Interim_data/data_for_ml_fungi_vs_bacteria_numFeat_with_AllSeqCenter_index_29Oct21.RData")

fungi_WIS_vs_nonWIS_ML <- function(metaData = metaFungiVsBacteria,
                                         dataFungi = countsFungiWGS,
                                         sampleType = "Primary Tumor",
                                         seqCenter = "Harvard Medical School",
                                         modelType = "xgbTree",
                                         wisIndex = indexWISfungi,
                                         nonwisIndex = inverseIndexNonWISfungi,
                                         numFeatVec = 34,
                                         numIter = 1,
                                         numResampleIter = 1,
                                         numKFold = 4){
  st <- gsub('([[:punct:]])|\\s+','',sampleType)
  sc <- "AllSeqCenters"
  wisFungiSamplingList <- list()
  nonwisFungiSamplingList <- list()
  for(zz in 1:length(numFeatVec)){
    numSamp <- numFeatVec[zz]
    wisFungiSamplingList[[zz]] <- replicate(numIter, expr = sample(wisIndex, size = numSamp, replace = FALSE))
    nonwisFungiSamplingList[[zz]] <- replicate(numIter, expr = sample(nonwisIndex, size = numSamp, replace = FALSE))
  }
  
  metaDataFilt <- metaData %>% filter(sample_type == sampleType) %>% droplevels()
  metaDataFilt$predY <- factor(gsub("^TCGA-","",metaDataFilt$investigation))
  dataFungiFilt <- dataFungi[rownames(metaDataFilt),]
  
  #----------------------------
  indexAllSeqCenterTrainList <- list()
  indexAllSeqCenterTestList <- list()
  splitVars <- c("data_submitting_center_label","sample_type","disease_type")
  for(yy in 1:numIter){
    set.seed(yy)
    stratSampling <- stratified(metaDataFilt, group = splitVars, size = 0.7,
                                keep.rownames = TRUE, bothSets = TRUE, replace = FALSE)
    split1Meta <- as.data.frame(stratSampling[[1]])
    rownames(split1Meta) <- split1Meta$rn
    split2Meta <- as.data.frame(stratSampling[[2]])
    rownames(split2Meta) <- split2Meta$rn
    
    indexAllSeqCenterTrainList[[yy]] <- which(rownames(metaDataFilt) %in% rownames(split1Meta))
    indexAllSeqCenterTestList[[yy]] <- which(rownames(metaDataFilt) %in% rownames(split2Meta))
  }
  
  #---------------------------
  
  perfCombined <- list()
  perfCombinedTmp <- list()
  multiClass_wisFungi <- list()
  multiClass_nonwisFungi <- list()
  for(jj in 1:length(numFeatVec)){ # should be 1:length(numFeatVec)
    wisFungiListTmp <- wisFungiSamplingList[[jj]]
    nonwisFungiListTmp <- nonwisFungiSamplingList[[jj]]
    numFeatTmp <- dim(wisFungiListTmp)[1]
    print(sprintf("Number of features: %d", numFeatTmp))
    for(ii in 1:numIter){ # should be 1:numIter
      iterNum <- ii
      wisfungiFeat <- wisFungiListTmp[,iterNum]
      nonwisFungiFeat <- nonwisFungiListTmp[,iterNum]
      mlDataY <- metaDataFilt
      mlDataX_wisFungi <- dataFungiFilt[,wisfungiFeat]
      mlDataX_nonwisFungi <- dataFungiFilt[,nonwisFungiFeat]
      
      print(sprintf("Iteration %d",iterNum))
      
      indexAllSeqCenterTrain <- indexAllSeqCenterTrainList[[iterNum]]
      indexAllSeqCenterTest <- indexAllSeqCenterTestList[[iterNum]]
      
      trainX_wisFungi <- mlDataX_wisFungi[indexAllSeqCenterTrain,]
      trainX_nonwisFungi <- mlDataX_nonwisFungi[indexAllSeqCenterTrain,]
      trainY <- mlDataY[indexAllSeqCenterTrain,"predY"]
      testX_wisFungi <- mlDataX_wisFungi[indexAllSeqCenterTest,]
      testX_nonwisFungi <- mlDataX_nonwisFungi[indexAllSeqCenterTest,]
      testY <- mlDataY[indexAllSeqCenterTest,"predY"]
      
      # set.seed(42) # have to restate seed again, as ctrl defines the cross validation sampling during training
      # ctrl <- trainControl(method = "repeatedcv",
      #                      number = numKFold,
      #                      repeats = numResampleIter,
      #                      sampling = "up",
      #                      summaryFunction = multiClassSummary,
      #                      classProbs = TRUE,
      #                      verboseIter = FALSE,
      #                      savePredictions = TRUE,
      #                      allowParallel=TRUE)
      
      print("Now training WIS fungi model...")
      set.seed(42)
      mlModel_wisFungi <- train(x = trainX_wisFungi,
                             y = trainY,
                             method = modelType,
                             preProcess = c("zv"),
                             nthread = 1,
                             trControl = trainControl(method = "repeatedcv",
                                                      number = numKFold,
                                                      repeats = numResampleIter,
                                                      sampling = "up",
                                                      summaryFunction = multiClassSummary,
                                                      classProbs = TRUE,
                                                      verboseIter = FALSE,
                                                      savePredictions = TRUE,
                                                      allowParallel=TRUE),
                             # metric = "ROC",
                             # tuneGrid = kfoldGBMGrid,
                             verbose = FALSE)
      
      print("Now training non-WIS fungi model...")
      set.seed(42)
      mlModel_nonwisFungi <- train(x = trainX_nonwisFungi,
                                y = trainY,
                                method = modelType,
                                preProcess = c("zv"),
                                nthread = 1,
                                trControl = trainControl(method = "repeatedcv",
                                                         number = numKFold,
                                                         repeats = numResampleIter,
                                                         sampling = "up",
                                                         summaryFunction = multiClassSummary,
                                                         classProbs = TRUE,
                                                         verboseIter = FALSE,
                                                         savePredictions = TRUE,
                                                         allowParallel=TRUE),
                                # metric = "ROC",
                                # tuneGrid = kfoldGBMGrid,
                                verbose = FALSE)
      
      print("Obtaining performance values...")
      multiClass_wisFungi <- data.frame(obs = testY,
                                     pred = predict(mlModel_wisFungi, newdata = testX_wisFungi),
                                     predict(mlModel_wisFungi, newdata = testX_wisFungi, type = "prob"))
      multiClass_nonwisFungi <- data.frame(obs = testY,
                                        pred = predict(mlModel_nonwisFungi, newdata = testX_nonwisFungi),
                                        predict(mlModel_nonwisFungi, newdata = testX_nonwisFungi, type = "prob"))
      
      wisFungiPerf <- multiClassSummary(multiClass_wisFungi, lev = levels(multiClass_wisFungi$obs))
      nonwisFungiPerf <- multiClassSummary(multiClass_nonwisFungi, lev = levels(multiClass_nonwisFungi$obs))
      
      wisFungiPerfDf <- data.frame(as.list(wisFungiPerf))
      colnames(wisFungiPerfDf) <- paste0("wis_",colnames(wisFungiPerfDf))
      
      nonwisFungiPerfDf <- data.frame(as.list(nonwisFungiPerf))
      colnames(nonwisFungiPerfDf) <- paste0("nonwis_",colnames(nonwisFungiPerfDf))
      
      perfCombined[[ii]] <- cbind(wisFungiPerfDf, nonwisFungiPerfDf, iter = paste0("iter",iterNum), numFeat = numFeatTmp)
      # if((iterNum %% 10) == 0){
      #   checkpointPerfCombined <- perfCombined
      #   save(checkpointPerfCombined, file = paste0("checkpoints/checkpoint_numFeat",numFeatTmp,"_iterNum",iterNum,"_",st,"_",sc,"_k",numKFold,"_modelType_",modelType,".RData"))
      # }
    }
    perfCombinedTmp[[jj]] <- do.call(rbind,perfCombined)
    cat('\n')
  }
  
  perfCombinedAll <- do.call(rbind, perfCombinedTmp)
  # baseFilename <- paste0("fungi_vs_bacteria_numFeat_",st,"_",sc,"_numIter",numIter,"_k",numKFold,"_modelType_",modelType)
  # write.csv(perfCombinedAll, file = paste0(baseFilename,".csv"))
  # save(perfCombinedAll, fungiSamplingList, bacteriaSamplingList, file = paste0(baseFilename,".RData"))
  
  res <- list(perfCombinedAll=perfCombinedAll,
              multiClass_wisFungi=multiClass_wisFungi,
              multiClass_nonwisFungi=multiClass_nonwisFungi,
              indexAllSeqCenterTrainList=indexAllSeqCenterTrainList,
              indexAllSeqCenterTestList=indexAllSeqCenterTestList)
  return(res)
}
tmp2 <- fungi_WIS_vs_nonWIS_ML(sampleType = "Primary Tumor")

#----------------------------------------------------------#
# Plot results of WIS-intersected features vs non-intersected features
#----------------------------------------------------------#

plotWISVsNonWISFungiPerf <- function(rDataFilePath, 
                                       seqCenter = "Harvard Medical School",
                                       allSeqCenterFlag = FALSE,
                                       inputSampleType="Primary Tumor",
                                       numIter = 100,
                                       modelType = "xgbTree",
                                       statSize = 2.5){
  if(inputSampleType=="Primary Tumor"){st = "PT"}
  if(inputSampleType=="Primary Tumor vs Solid Tissue Normal"){st = "PT_vs_NAT"}
  if(inputSampleType=="Blood Derived Normal"){st = "BDN"}
  if(modelType=="rf"){modelTypeFormatted = "Random Forest"}
  if(modelType %in% c("gbm","xgbTree")){modelTypeFormatted = "Gradient Boosting"}
  if(allSeqCenterFlag){
    seqCenterFormatted <- "AllSeqCenters"
    seqCenterPlotTitle <- "All WGS Sequencing Centers"
  } else{
    seqCenterFormatted <- gsub('([[:punct:]])|\\s+','',seqCenter)
    seqCenterPlotTitle <- seqCenter
  }
  
  load(rDataFilePath) # loads perfCombinedAll object
  keepCols <- c("AUC","prAUC")
  perfCombinedAllFormatted <- reshape2::melt(perfCombinedAll, id.vars = c("iter","numFeat")) %>%
    mutate(Domain = factor(ifelse(grepl("^wis_",variable),yes = "Intersected",no = "Nonintersected"), levels = c("Intersected","Nonintersected"))) %>%
    mutate(variable = gsub("^wis_|^nonwis_","",variable)) %>% filter(variable %in% keepCols)
  perfCombinedAllFormatted %>%
    filter(variable == "AUC") %>%
    ggboxplot(x = "Domain",
           y = "value",
           palette = "nejm",
           fill = "Domain",
           legend = "none",
           add = "jitter",
           add.params = list(alpha=0.4),
           notch = TRUE,
           xlab = "TCGA fungi did or did not intersect with WIS data",
           ylab = "Average pan-cancer AUROC",
           numeric.x.axis = TRUE) +
    stat_compare_means(comparisons = list(c("Intersected","Nonintersected")), label = "p.format", size = statSize) -> plotROC

  perfCombinedAllFormatted %>%
    filter(variable == "prAUC") %>%
    ggboxplot(x = "Domain",
           y = "value",
           palette = "nejm",
           fill = "Domain",
           legend = "none",
           add = "jitter",
           add.params = list(alpha=0.4),
           notch = TRUE,
           xlab = "TCGA fungi did or did not intersect with WIS data",
           ylab = "Average pan-cancer AUPR",
           numeric.x.axis = TRUE) +
    stat_compare_means(comparisons = list(c("Intersected","Nonintersected")), label = "p.format", size = statSize) -> plotPR

  combinedPlotTitle <- paste0("TCGA Simulations: WIS-intersecting fungi vs. non-WIS-intersecting fungi species (equally fixed at 34 taxa)\n",
                              inputSampleType," | ",seqCenterPlotTitle," | ",modelTypeFormatted,
                              "\n(",numIter," iterations of multi-class classification per box)")
  combinedPlot <- ggarrange(plotROC, plotPR, ncol = 2, common.legend = TRUE)
  combinedPlotAnnotated <- annotate_figure(combinedPlot,
                                           top = text_grob(combinedPlotTitle,
                                                           color = "black", face = "bold", size = 14))
  print(combinedPlotAnnotated)
  ggsave(filename = paste0("Figures/Supplementary_Figures/wis_vs_nonwis_fungi",seqCenterFormatted,"_",st,"_numIter",numIter,"_",modelType,".jpeg"),
         plot = combinedPlotAnnotated,
         dpi = "retina", units = "in", width = 12, height = 5)
  rm(perfCombinedAll)
  # return(perfCombinedAllFormatted)
}
# Primary Tumor
plotWISVsNonWISFungiPerf(rDataFilePath = "Interim_data/wis_vs_nonwis_fungi_PrimaryTumor_AllSeqCenters_numIter100_k4_modelType_xgbTree.RData", 
                         seqCenter = "AllSeqCenters", allSeqCenterFlag = TRUE, modelType = "xgbTree", numIter = 100, inputSampleType = "Primary Tumor", statSize=3)
# Blood Derived Normal
plotWISVsNonWISFungiPerf(rDataFilePath = "Interim_data/wis_vs_nonwis_fungi_BloodDerivedNormal_AllSeqCenters_numIter100_k4_modelType_xgbTree.RData", 
                         seqCenter = "AllSeqCenters", allSeqCenterFlag = TRUE, modelType = "xgbTree", numIter = 100, inputSampleType = "Blood Derived Normal", statSize=3)

#----------------------------------------------------------#
# Build multi-class ML function by number of features
#----------------------------------------------------------#

fungi_vs_bacteria_numFeat_ML <- function(metaData = metaFungiVsBacteria,
                                         dataFungi = countsFungiWGS,
                                         dataBacteria = countsBacteriaWGS,
                                         sampleType = "Primary Tumor",
                                         seqCenter = "Harvard Medical School",
                                         modelType = "gbm",
                                         numFeatVec = c(10,25,50,75,100,125,150,175,200,225,250,275,300),
                                         numIter = 1000,
                                         numResampleIter = 1,
                                         numKFold = 4){
  st <- gsub('([[:punct:]])|\\s+','',sampleType)
  sc <- gsub('([[:punct:]])|\\s+','',seqCenter)
  fungiSamplingList <- list()
  bacteriaSamplingList <- list()
  for(zz in 1:length(numFeatVec)){
    numSamp <- numFeatVec[zz]
    fungiSamplingList[[zz]] <- replicate(numIter, expr = sample(1:ncol(dataFungi), size = numSamp, replace = FALSE))
    bacteriaSamplingList[[zz]] <- replicate(numIter, expr = sample(1:ncol(dataBacteria), size = numSamp, replace = FALSE))
  }
  
  metaDataFilt <- metaData %>% filter(sample_type == sampleType, data_submitting_center_label == seqCenter) %>% droplevels()
  metaDataFilt$predY <- factor(gsub("^TCGA-","",metaDataFilt$investigation))
  dataFungiFilt <- dataFungi[rownames(metaDataFilt),]
  dataBacteriaFilt <- dataBacteria[rownames(metaDataFilt),]
  
  perfCombined <- list()
  perfCombinedTmp <- list()
  for(jj in 1:length(numFeatVec)){ # should be length of numFeatVec
    fungiListTmp <- fungiSamplingList[[jj]]
    bacteriaListTmp <- bacteriaSamplingList[[jj]]
    numFeatTmp <- dim(fungiListTmp)[1]
    print(sprintf("Number of features: %d", numFeatTmp))
    for(ii in 1:2){ # should be numIter
      iterNum <- ii
      fungiFeat <- fungiListTmp[,iterNum]
      bacteriaFeat <- bacteriaListTmp[,iterNum]
      mlDataY <- metaDataFilt
      mlDataX_Fungi <- dataFungiFilt[,fungiFeat]
      mlDataX_Bacteria <- dataBacteriaFilt[,bacteriaFeat]
      
      print(sprintf("Iteration %d",iterNum))
      
      set.seed(42)
      index <- createDataPartition(mlDataY[,"predY"], p = 0.7, list = FALSE)
      trainX_Fungi <- mlDataX_Fungi[index,]
      trainX_Bacteria <- mlDataX_Bacteria[index,]
      trainY <- mlDataY[index,"predY"]
      testX_Fungi <- mlDataX_Fungi[-index,]
      testX_Bacteria <- mlDataX_Bacteria[-index,]
      testY <- mlDataY[-index,"predY"]
      
      set.seed(42) # have to restate seed again, as ctrl defines the cross validation sampling during training
      ctrl <- trainControl(method = "repeatedcv",
                           number = numKFold,
                           repeats = numResampleIter,
                           sampling = "up",
                           summaryFunction = multiClassSummary,
                           classProbs = TRUE,
                           verboseIter = FALSE,
                           savePredictions = TRUE,
                           allowParallel=TRUE)
      
      print("Now training fungi model...")
      set.seed(42)
      mlModel_Fungi <- train(x = trainX_Fungi,
                             y = trainY,
                             method = modelType,
                             preProcess = c("zv"),
                             trControl = ctrl,
                             # metric = "ROC",
                             # tuneGrid = kfoldGBMGrid,
                             verbose = FALSE)
      
      print("Now training bacteria model...")
      set.seed(42)
      mlModel_Bacteria <- train(x = trainX_Bacteria,
                                y = trainY,
                                method = modelType,
                                preProcess = c("zv"),
                                trControl = ctrl,
                                # metric = "ROC",
                                # tuneGrid = kfoldGBMGrid,
                                verbose = FALSE)
      
      print("Obtaining performance values...")
      multiClass_Fungi <- data.frame(obs = testY,
                                     pred = predict(mlModel_Fungi, newdata = testX_Fungi),
                                     predict(mlModel_Fungi, newdata = testX_Fungi, type = "prob"))
      multiClass_Bacteria <- data.frame(obs = testY,
                                        pred = predict(mlModel_Bacteria, newdata = testX_Bacteria),
                                        predict(mlModel_Bacteria, newdata = testX_Bacteria, type = "prob"))
      
      fungiPerf <- multiClassSummary(multiClass_Fungi, lev = levels(multiClass_Fungi$obs))
      bacteriaPerf <- multiClassSummary(multiClass_Bacteria, lev = levels(multiClass_Bacteria$obs))
      
      fungiPerfDf <- data.frame(as.list(fungiPerf))
      colnames(fungiPerfDf) <- paste0("f_",colnames(fungiPerfDf))
      
      bacteriaPerfDf <- data.frame(as.list(bacteriaPerf))
      colnames(bacteriaPerfDf) <- paste0("b_",colnames(bacteriaPerfDf))
      
      perfCombined[[ii]] <- cbind(fungiPerfDf, bacteriaPerfDf, iter = paste0("iter",iterNum), numFeat = numFeatTmp)
    }
    perfCombinedTmp[[jj]] <- do.call(rbind,perfCombined)
    cat('\n')
  }
  
  perfCombinedAll <- do.call(rbind, perfCombinedTmp)
  baseFilename <- paste0("fungi_vs_bacteria_numFeat_",st,"_",sc,"_numIter",numIter,"_k",numKFold)
  write.csv(perfCombinedAll, file = paste0(baseFilename,".csv"))
  save(perfCombinedAll, file = paste0(baseFilename,".RData"))
}
fungi_vs_bacteria_numFeat_ML(sampleType = "Primary Tumor", seqCenter = "Harvard Medical School")

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

save(metaFungiVsBacteria,
     countsFungiWGS,
     countsBacteriaWGS,
     split1MetaFungiVsBacteria,
     split2MetaFungiVsBacteria,
     file = "Interim_data/data_for_ml_fungi_vs_bacteria_numFeat_with_AllSeqCenter_index_29Oct21.RData")

save(metaFungiVsBacteria,
     countsFungiWGS,
     countsBacteriaWGS,
     split1MetaFungiVsBacteria,
     split2MetaFungiVsBacteria,
     rep200TaxSplit_Fungi,
     rep200TaxSplit_Bacteria,
     file = "Interim_data/data_for_ml_fungi_vs_bacteria_numFeat_with_AllSeqCenter_index_03Nov21.RData")

#----------------------------------------------------------#
# NumFeat: Plot results per seq center and WGS pan seq center
#----------------------------------------------------------#

plotFungiVsBacteriaNumFeat <- function(rDataFilePath, 
                                       seqCenter = "Harvard Medical School",
                                       allSeqCenterFlag = FALSE,
                                       inputSampleType="Primary Tumor",
                                       numIter = 100,
                                       modelType = "rf",
                                       outputPlotFlag = FALSE,
                                       statSize = 2.5){
  if(inputSampleType=="Primary Tumor"){st = "PT"}
  if(inputSampleType=="Primary Tumor vs Solid Tissue Normal"){st = "PT_vs_NAT"}
  if(inputSampleType=="Blood Derived Normal"){st = "BDN"}
  if(modelType=="rf"){modelTypeFormatted = "Random Forest"}
  if(modelType %in% c("gbm","xgbTree")){modelTypeFormatted = "Gradient Boosting"}
  if(allSeqCenterFlag){
    seqCenterFormatted <- "AllSeqCenters"
    seqCenterPlotTitle <- "All WGS Sequencing Centers"
    } else{
    seqCenterFormatted <- gsub('([[:punct:]])|\\s+','',seqCenter)
    seqCenterPlotTitle <- seqCenter
    }
  
  load(rDataFilePath) # loads perfCombinedAll object
  keepCols <- c("AUC","prAUC")
  perfCombinedAllFormatted <- reshape2::melt(perfCombinedAll, id.vars = c("iter","numFeat")) %>%
    mutate(Domain = factor(ifelse(grepl("f_",variable),yes = "Fungi",no = "Bacteria"), levels = c("Fungi","Bacteria"))) %>%
    mutate(variable = gsub("^f_|^b_","",variable)) %>% filter(variable %in% keepCols)
  perfCombinedAllFormatted %>%
    filter(variable == "AUC") %>%
    ggline(x = "numFeat",
           y = "value",
           palette = "nejm",
           color = "Domain",
           add = "mean_ci",
           xlab = "Number of microbial taxa in pan-cancer model",
           ylab = "Average AUROC",
           numeric.x.axis = TRUE) +
    scale_x_continuous(breaks = seq(0, 300, by = 25), limits = c(0,300)) -> plotROC
  plotROC + stat_compare_means(aes(group=Domain), label = "p.signif", size = statSize) -> sigPlotROC
  
  perfCombinedAllFormatted %>%
    filter(variable == "prAUC") %>%
    ggline(x = "numFeat",
           y = "value",
           palette = "nejm",
           color = "Domain",
           add = "mean_ci",
           xlab = "Number of microbial taxa in pan-cancer model",
           ylab = "Average AUPR",
           numeric.x.axis = TRUE) +
    scale_x_continuous(breaks = seq(0, 300, by = 25), limits = c(0,300)) -> plotPR
  plotPR + stat_compare_means(aes(group=Domain), label = "p.signif", size = statSize) -> sigPlotPR
  
  combinedPlotTitle <- paste0("TCGA Simulations: Fungi vs. Bacteria (varying number of taxa)\n",
                              inputSampleType," | ",seqCenterPlotTitle," | ",modelTypeFormatted,
                              "\n(",numIter," iterations of multi-class classification per point)")
  combinedPlot <- ggarrange(sigPlotROC, sigPlotPR, ncol = 2, common.legend = TRUE) 
  combinedPlotAnnotated <- annotate_figure(combinedPlot, 
                                           top = text_grob(combinedPlotTitle, 
                                                           color = "black", face = "bold", size = 14))
  print(combinedPlotAnnotated)
  ggsave(filename = paste0("Figures/Supplementary_Figures/fungi_vs_bacteria_numFeat",seqCenterFormatted,"_",st,"_numIter",numIter,"_",modelType,".jpeg"),
         plot = combinedPlotAnnotated,
         dpi = "retina", units = "in", width = 12, height = 5)
  rm(perfCombinedAll)
  if(outputPlotFlag){
    res <- list(plotROC=plotROC,
                plotPR=plotPR)
    return(res)
  }
}
# Primary Tumor
plotFungiVsBacteriaNumFeat(rDataFilePath = "Interim_data/fungi_vs_bacteria_numFeat_PrimaryTumor_HarvardMedicalSchool_numIter100_k4_modelType_rf.RData", 
                           seqCenter = "Harvard Medical School", statSize=3)
plotFungiVsBacteriaNumFeat(rDataFilePath = "Interim_data/fungi_vs_bacteria_numFeat_PrimaryTumor_BroadInstituteofMITandHarvard_numIter100_k4_modelType_rf.RData", 
                           seqCenter = "Broad Institute", statSize=3)
plotFungiVsBacteriaNumFeat(rDataFilePath = "Interim_data/fungi_vs_bacteria_numFeat_PrimaryTumor_BaylorCollegeofMedicine_numIter100_k4_modelType_rf.RData", 
                           seqCenter = "Baylor College of Medicine", statSize=3)
plotFungiVsBacteriaNumFeat(rDataFilePath = "Interim_data/fungi_vs_bacteria_numFeat_PrimaryTumor_MDAndersonInstituteforAppliedCancerScience_numIter100_k4_modelType_rf.RData", 
                           seqCenter = "MD Anderson", statSize=3)
plotFungiVsBacteriaNumFeat(rDataFilePath = "Interim_data/fungi_vs_bacteria_numFeat_PrimaryTumor_WashingtonUniversitySchoolofMedicine_numIter100_k4_modelType_rf.RData", 
                           seqCenter = "WashU", statSize=3)
# Blood Derived Normal
plotFungiVsBacteriaNumFeat(rDataFilePath = "Interim_data/fungi_vs_bacteria_numFeat_BloodDerivedNormal_HarvardMedicalSchool_numIter100_k4_modelType_rf.RData", 
                           seqCenter = "Harvard Medical School", inputSampleType = "Blood Derived Normal", statSize=3)
plotFungiVsBacteriaNumFeat(rDataFilePath = "Interim_data/fungi_vs_bacteria_numFeat_BloodDerivedNormal_BroadInstituteofMITandHarvard_numIter100_k4_modelType_rf.RData", 
                           seqCenter = "Broad Institute", inputSampleType = "Blood Derived Normal", statSize=3)

## Pan-seq-center: Primary Tumor
plotFungiVsBacteriaPT_WGS <- plotFungiVsBacteriaNumFeat(rDataFilePath = "Interim_data/fungi_vs_bacteria_numFeat_PrimaryTumor_AllSeqCenters_numIter100_k4_modelType_xgbTree.RData", 
                           seqCenter = "AllSeqCenters", allSeqCenterFlag = TRUE, modelType = "xgbTree", numIter = 100, inputSampleType = "Primary Tumor", statSize=3,
                           outputPlotFlag = TRUE)
plotFungiVsBacteriaNumFeat(rDataFilePath = "Interim_data/fungi_vs_bacteria_numFeat_PrimaryTumor_AllSeqCenters_numIter100_k4_modelType_rf.RData", 
                           seqCenter = "AllSeqCenters", allSeqCenterFlag = TRUE, modelType = "rf", numIter = 100, inputSampleType = "Primary Tumor", statSize=3)
plotFungiVsBacteriaNumFeat(rDataFilePath = "Interim_data/fungi_vs_bacteria_numFeat_PrimaryTumor_AllSeqCenters_numIter20_k4_modelType_xgbTree.RData", 
                           seqCenter = "AllSeqCenters", allSeqCenterFlag = TRUE, modelType = "xgbTree", numIter = 20, inputSampleType = "Primary Tumor", statSize=3)

## Pan-seq-center: Blood Derived Normal
plotFungiVsBacteriaBDN_WGS <- plotFungiVsBacteriaNumFeat(rDataFilePath = "Interim_data/fungi_vs_bacteria_numFeat_BloodDerivedNormal_AllSeqCenters_numIter100_k4_modelType_xgbTree.RData", 
                           seqCenter = "AllSeqCenters", allSeqCenterFlag = TRUE, modelType = "xgbTree", numIter = 100, inputSampleType = "Blood Derived Normal", statSize=3,
                           outputPlotFlag = TRUE)
plotFungiVsBacteriaNumFeat(rDataFilePath = "Interim_data/fungi_vs_bacteria_numFeat_BloodDerivedNormal_AllSeqCenters_numIter100_k4_modelType_rf.RData", 
                           seqCenter = "AllSeqCenters", allSeqCenterFlag = TRUE, modelType = "rf", numIter = 100, inputSampleType = "Blood Derived Normal", statSize=3)
plotFungiVsBacteriaNumFeat(rDataFilePath = "Interim_data/fungi_vs_bacteria_numFeat_BloodDerivedNormal_AllSeqCenters_numIter20_k4_modelType_xgbTree.RData", 
                           seqCenter = "AllSeqCenters", allSeqCenterFlag = TRUE, modelType = "xgbTree", numIter = 20, inputSampleType = "Blood Derived Normal", statSize=3)

#----------------------------------------------------------#
# Compare number feature performance using ranked pan-cancer features
#----------------------------------------------------------#

# See script S26 for determining ranked fungal and bacterial features and 
# calculating concomitant machine learning performance

plotOverlayFungiVsBacteriaNumFeatRankedPerf <- function(rDataFilePath, 
                                                        origPlotROC = plotFungiVsBacteriaPT_WGS$plotROC,
                                                        origPlotPR = plotFungiVsBacteriaPT_WGS$plotPR,
                                                        seqCenter = "AllSeqCenters",
                                                        allSeqCenterFlag = TRUE,
                                                        inputSampleType="Primary Tumor",
                                                        numIter = 100,
                                                        modelType = "xgbTree",
                                                        rankedFeatModelType = "rf",
                                                        trainingDataSetOnlyRanking = FALSE,
                                                        statSize = 3){
  if(inputSampleType=="Primary Tumor"){st = "PT"}
  if(inputSampleType=="Primary Tumor vs Solid Tissue Normal"){st = "PT_vs_NAT"}
  if(inputSampleType=="Blood Derived Normal"){st = "BDN"}
  if(rankedFeatModelType=="rf"){rankedFeatModelTypeFormatted = "Random Forest"}
  if(rankedFeatModelType=="rf" & trainingDataSetOnlyRanking){rankedFeatModelTypeFormatted = "Random Forest on training data only"}
  if(rankedFeatModelType %in% c("gbm","xgbTree")){rankedFeatModelTypeFormatted = "Gradient Boosting"}
  if(modelType=="rf"){modelTypeFormatted = "Random Forest"}
  if(modelType %in% c("gbm","xgbTree")){modelTypeFormatted = "Gradient Boosting"}
  if(allSeqCenterFlag){
    seqCenterFormatted <- "AllSeqCenters"
    seqCenterPlotTitle <- "All WGS Sequencing Centers"
  } else{
    seqCenterFormatted <- gsub('([[:punct:]])|\\s+','',seqCenter)
    seqCenterPlotTitle <- seqCenter
  }
  
  load(rDataFilePath)
  keepCols <- c("AUC","prAUC")
  perfCombinedAllFormatted <- reshape2::melt(perfCombinedAll, id.vars = c("iter","numFeat")) %>%
    mutate(Domain = factor(ifelse(grepl("f_",variable),yes = "Ranked Fungi",no = "Ranked Bacteria"), levels = c("Ranked Fungi","Ranked Bacteria"))) %>%
    mutate(variable = gsub("^f_|^b_","",variable)) %>% filter(variable %in% keepCols)
  perfCombinedAllFormattedAUROC <- perfCombinedAllFormatted %>% filter(variable == "AUC")
  perfCombinedAllFormattedAUPR <- perfCombinedAllFormatted %>% filter(variable == "prAUC")
  
  rocPlotWithRankedFeatPerf <- origPlotROC +
    geom_line(data = perfCombinedAllFormattedAUROC, aes(x = numFeat, y = value, color = Domain), linetype = "dashed") +
    geom_point(data = perfCombinedAllFormattedAUROC, aes(x = numFeat, y = value, color = Domain), shape=17)
  prPlotWithRankedFeatPerf <-  origPlotPR +
    geom_line(data = perfCombinedAllFormattedAUPR, aes(x = numFeat, y = value, color = Domain), linetype = "dashed") +
    geom_point(data = perfCombinedAllFormattedAUPR, aes(x = numFeat, y = value, color = Domain), shape=17)
  
  combinedPlotTitle <- paste0("TCGA Simulations: Fungi vs. Bacteria (varying number of taxa)\n",
                              "Overlaid with performance using ranked features (ranked by ",rankedFeatModelTypeFormatted,")\n",
                              inputSampleType," | ",seqCenterPlotTitle," | ",modelTypeFormatted)
  combinedRankedFeatPerfPlot <- ggarrange(rocPlotWithRankedFeatPerf, prPlotWithRankedFeatPerf, ncol = 2, common.legend = TRUE) 
  combinedRankedFeatPerfPlotAnnotated <- annotate_figure(combinedRankedFeatPerfPlot, 
                                                         top = text_grob(combinedPlotTitle, 
                                                                         color = "black", face = "bold", size = 14))
  print(combinedRankedFeatPerfPlotAnnotated)
  if(trainingDataSetOnlyRanking){
    ggFilename <- paste0("Figures/Supplementary_Figures/fungi_vs_bacteria_overlay_trainingDataOnly_rankedNumFeat",
                         seqCenterFormatted,"_",st,"_numIter",numIter,"_",modelType,"rankedFeatModelType_",
                         rankedFeatModelType,".jpeg")
  } else{
    ggFilename <- paste0("Figures/Supplementary_Figures/fungi_vs_bacteria_overlay_rankedNumFeat",
                         seqCenterFormatted,"_",st,"_numIter",numIter,"_",modelType,"rankedFeatModelType_",
                         rankedFeatModelType,".jpeg")
  }
  ggsave(filename = ggFilename,
         plot = combinedRankedFeatPerfPlotAnnotated,
         dpi = "retina", units = "in", width = 12, height = 5)
  rm(perfCombinedAll)
}
## Primary Tumor
plotOverlayFungiVsBacteriaNumFeatRankedPerf("Interim_data/fungi_vs_bacteria_rankedNumFeat_PrimaryTumor_AllSeqCenters_numIter1_k4_modelType_xgbTree_rankedFeatModelType_rf.RData")
plotOverlayFungiVsBacteriaNumFeatRankedPerf("Interim_data/fungi_vs_bacteria_rankedNumFeat_PrimaryTumor_AllSeqCenters_numIter1_k4_modelType_xgbTree_rankedFeatModelType_xgbTree.RData",
                                            rankedFeatModelType = "xgbTree")
# Ranking on training data set only
plotOverlayFungiVsBacteriaNumFeatRankedPerf("Interim_data/fungi_vs_bacteria_trainingSetOnly_rankedNumFeat_PrimaryTumor_AllSeqCenters_numIter1_k4_modelType_xgbTree_rankedFeatModelType_rf.RData",
                                            rankedFeatModelType = "rf", trainingDataSetOnlyRanking = TRUE)

## Blood Derived Normal
plotOverlayFungiVsBacteriaNumFeatRankedPerf("Interim_data/fungi_vs_bacteria_rankedNumFeat_BloodDerivedNormal_AllSeqCenters_numIter1_k4_modelType_xgbTree_rankedFeatModelType_rf.RData",
                                            inputSampleType = "Blood Derived Normal", origPlotROC = plotFungiVsBacteriaBDN_WGS$plotROC, origPlotPR = plotFungiVsBacteriaBDN_WGS$plotPR)
# Ranking on training data set only
plotOverlayFungiVsBacteriaNumFeatRankedPerf("Interim_data/fungi_vs_bacteria_trainingSetOnly_rankedNumFeat_BloodDerivedNormal_AllSeqCenters_numIter1_k4_modelType_xgbTree_rankedFeatModelType_rf.RData",
                                            inputSampleType = "Blood Derived Normal", origPlotROC = plotFungiVsBacteriaBDN_WGS$plotROC, origPlotPR = plotFungiVsBacteriaBDN_WGS$plotPR, trainingDataSetOnlyRanking = TRUE)

#----------------------------------------------------------#
# Plot results per cancer type when looking across WGS pan seq center data
#----------------------------------------------------------#

plotFungiVsBacteriaNumFeat_PerCT <- function(rDataFilePath, 
                                       allSeqCenterFlag = FALSE,
                                       inputSampleType="Primary Tumor",
                                       outputPlotFlag = FALSE,
                                       numIter = 100,
                                       numFeatLimit = 25,
                                       modelType = "gbm",
                                       ciVisibleFlag = FALSE,
                                       statSize = 2.5){
  # NOTE: To make CIs visible, set the following ggline plot options:
  # plot_type = "l"
  if(inputSampleType=="Primary Tumor"){st = "PT"}
  if(inputSampleType=="Blood Derived Normal"){st = "BDN"}
  if(modelType=="rf"){modelTypeFormatted = "Random Forest"}
  if(modelType %in% c("gbm","xgbTree")){modelTypeFormatted = "Gradient Boosting"}
  seqCenterFormatted <- "AllSeqCenters"
  seqCenterPlotTitle <- "All WGS Sequencing Centers"
  
  load(rDataFilePath) # loads perfCombinedAll object
  # keepCols <- c("auroc","aupr")
  perfCombinedAllFormatted <- reshape2::melt(perfCombinedAll, id.vars = c("iter","numFeat","diseaseType","sampleType","seqCenter")) %>%
    mutate(Domain = factor(ifelse(grepl("f_",variable),yes = "Fungi",no = "Bacteria"), levels = c("Fungi","Bacteria"))) %>%
    mutate(variable = gsub("^f_|^b_","",variable))
  perfCombinedAllFormatted %>%
    filter(variable == "auroc") %>%
    ggline(x = "numFeat",
           y = "value",
           palette = "nejm",
           color = "Domain",
           plot_type = ifelse(ciVisibleFlag, yes = "l", no = "b"),
           size = 0.1,
           # point.size = 0.01,
           add = "mean_ci",
           facet.by = "diseaseType",
           xlab = "Number of microbial taxa",
           ylab = "AUROC (1 vs All)",
           numeric.x.axis = TRUE) +
    rotate_x_text(90) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1.2)) +
    scale_x_continuous(breaks = seq(0, numFeatLimit, by = 25), limits = c(0,numFeatLimit)) -> plotROC
  plotROC + stat_compare_means(aes(group=Domain), label.y = 1, label = "p.signif", size = statSize) -> sigPlotROC

  perfCombinedAllFormatted %>%
    filter(variable == "aupr") %>%
    ggline(x = "numFeat",
           y = "value",
           plot_type = ifelse(ciVisibleFlag, yes = "l", no = "b"),
           size = 0.1,
           palette = "nejm",
           color = "Domain",
           add = "mean_ci",
           facet.by = "diseaseType",
           xlab = "Number of microbial taxa model",
           ylab = "AUPR (1 vs All)",
           numeric.x.axis = TRUE) +
    rotate_x_text(90) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0,1.2)) +
    scale_x_continuous(breaks = seq(0, numFeatLimit, by = 25), limits = c(0,numFeatLimit)) -> plotPR
  plotPR + stat_compare_means(aes(group=Domain), label.y = 1, label = "p.signif", size = statSize) -> sigPlotPR

  combinedPlotTitle <- paste0("TCGA Simulations per cancer type: Fungi vs. Bacteria (varying number of taxa)\n",
                              inputSampleType," | ",seqCenterPlotTitle," | ",modelTypeFormatted,
                              "\n(",numIter," iterations of one-versus-all-others classification per point)")
  combinedPlot <- ggarrange(sigPlotROC, sigPlotROC, nrow = 2, labels = "AUTO", font.label = list(size = 18, face = "bold"), common.legend = TRUE)
  combinedPlotAnnotated <- annotate_figure(combinedPlot,
                                           top = text_grob(combinedPlotTitle,
                                                           color = "black", face = "bold", size = 14))
  print(combinedPlotAnnotated)
  ggsave(filename = paste0("Figures/Supplementary_Figures/fungi_vs_bacteria_perCT_numFeat",seqCenterFormatted,"_",st,"_numIter",numIter,"_",modelType,".jpeg"),
         plot = combinedPlotAnnotated,
         dpi = "retina", units = "in", width = 14, height = 15)
  rm(perfCombinedAll)
  if(outputPlotFlag){
    res <- list(plotROC=plotROC,
                plotPR=plotPR)
    return(res)
  }
}
# Primary Tumor
plotFungiVsBacteriaNumFeat_PerCT_PT_WGS <- plotFungiVsBacteriaNumFeat_PerCT(rDataFilePath = "Interim_data/fungi_vs_bacteria_perCT_numFeat_10_25_50_75_100_125_150_175_200_225_250_275_300_PrimaryTumor_AllSeqCenters_numIter100_k4_modelType_gbm.RData", 
                                                                         outputPlotFlag = TRUE, modelType = "gbm", numFeatLimit = 300, numIter = 100, inputSampleType = "Primary Tumor", statSize=2)
# Blood Derived Normal
plotFungiVsBacteriaNumFeat_PerCT_BDN_WGS <- plotFungiVsBacteriaNumFeat_PerCT(rDataFilePath = "Interim_data/fungi_vs_bacteria_perCT_numFeat_10_25_50_75_100_125_150_175_200_225_250_275_300_BloodDerivedNormal_AllSeqCenters_numIter100_k4_modelType_gbm.RData", 
                                                                             outputPlotFlag = TRUE, modelType = "gbm", numFeatLimit = 300, numIter = 100, inputSampleType = "Blood Derived Normal", statSize=2)

#----------------------------------------------------------#
# Overlay per CT data with perf using ranked features
#----------------------------------------------------------#

plotOverlayFungiVsBacteriaNumFeatRankedPerf_PerCT <- function(rDataFilePath, 
                                                        origPlotROC = plotFungiVsBacteriaNumFeat_PerCT_PT_WGS$plotROC,
                                                        origPlotPR = plotFungiVsBacteriaNumFeat_PerCT_PT_WGS$plotPR,
                                                        seqCenter = "AllSeqCenters",
                                                        allSeqCenterFlag = TRUE,
                                                        rankedFeatLinetype = "solid",
                                                        inputSampleType="Primary Tumor",
                                                        numIter = 100,
                                                        modelType = "xgbTree",
                                                        rankedFeatModelType = "rf",
                                                        statSize = 3){
  if(inputSampleType=="Primary Tumor"){st = "PT"}
  if(inputSampleType=="Primary Tumor vs Solid Tissue Normal"){st = "PT_vs_NAT"}
  if(inputSampleType=="Blood Derived Normal"){st = "BDN"}
  if(rankedFeatModelType=="rf"){rankedFeatModelTypeFormatted = "Random Forest"}
  if(rankedFeatModelType %in% c("gbm","xgbTree")){rankedFeatModelTypeFormatted = "Gradient Boosting"}
  if(modelType=="rf"){modelTypeFormatted = "Random Forest"}
  if(modelType %in% c("gbm","xgbTree")){modelTypeFormatted = "Gradient Boosting"}
  if(allSeqCenterFlag){
    seqCenterFormatted <- "AllSeqCenters"
    seqCenterPlotTitle <- "All WGS Sequencing Centers"
  } else{
    seqCenterFormatted <- gsub('([[:punct:]])|\\s+','',seqCenter)
    seqCenterPlotTitle <- seqCenter
  }
  
  load(rDataFilePath)
  keepCols <- c("AUC","prAUC")
  perfCombinedAllFormatted <- reshape2::melt(perfCombinedAll, id.vars = c("iter","numFeat","diseaseType","sampleType","seqCenter","rankedFeatureModelType")) %>%
    mutate(Domain = factor(ifelse(grepl("f_",variable),yes = "Ranked Fungi",no = "Ranked Bacteria"), levels = c("Ranked Fungi","Ranked Bacteria"))) %>%
    mutate(variable = gsub("^f_|^b_","",variable))
  ciVisibleFlag <- FALSE
  numFeatLimit <- 300
  perfCombinedAllFormattedAUROC <- perfCombinedAllFormatted %>% filter(variable == "auroc")
  perfCombinedAllFormattedAUPR <- perfCombinedAllFormatted %>% filter(variable == "aupr")
  
  rocPlotWithRankedFeatPerf <- origPlotROC +
    geom_line(data = perfCombinedAllFormattedAUROC, aes(x = numFeat, y = value, color = Domain), linetype = rankedFeatLinetype, size = 0.5) +
    geom_point(data = perfCombinedAllFormattedAUROC, aes(x = numFeat, y = value, color = Domain), shape=17, size = 1.5)
  prPlotWithRankedFeatPerf <-  origPlotPR +
    geom_line(data = perfCombinedAllFormattedAUPR, aes(x = numFeat, y = value, color = Domain), linetype = rankedFeatLinetype, size = 0.5) +
    geom_point(data = perfCombinedAllFormattedAUPR, aes(x = numFeat, y = value, color = Domain), shape=17, size = 1.5)
  
  combinedPlotTitle <- paste0("TCGA Simulations: Fungi vs. Bacteria (varying number of taxa)\n",
                              "Overlaid with performance using ranked features (ranked by ",rankedFeatModelTypeFormatted,")\n",
                              inputSampleType," | ",seqCenterPlotTitle," | ",modelTypeFormatted)
  
  combinedPlot <- ggarrange(rocPlotWithRankedFeatPerf, prPlotWithRankedFeatPerf,
                            nrow = 2, labels = "AUTO", font.label = list(size = 18, face = "bold"), common.legend = TRUE)
  combinedPlotAnnotated <- annotate_figure(combinedPlot,
                                           top = text_grob(combinedPlotTitle,
                                                           color = "black", face = "bold", size = 14))
  print(combinedPlotAnnotated)
  ggsave(filename = paste0("Figures/Supplementary_Figures/fungi_vs_bacteria_perCT_rankedNumFeat",seqCenterFormatted,"_",st,"_numIter",numIter,"_",modelType,".jpeg"),
         plot = combinedPlotAnnotated,
         dpi = "retina", units = "in", width = 14, height = 15)
  rm(perfCombinedAll)
}
# Primary Tumor
plotOverlayFungiVsBacteriaNumFeatRankedPerf_PerCT("Interim_data/fungi_vs_bacteria_perCT_rankedNumFeat_10_25_50_75_100_125_150_175_200_225_250_275_300_PrimaryTumor_AllSeqCenters_numIter1_k4_modelType_gbm_rankedFeatModelType_rf.RData")
# Blood Derived Normal
plotOverlayFungiVsBacteriaNumFeatRankedPerf_PerCT("Interim_data/fungi_vs_bacteria_perCT_rankedNumFeat_10_25_50_75_100_125_150_175_200_225_250_275_300_BloodDerivedNormal_AllSeqCenters_numIter1_k4_modelType_gbm_rankedFeatModelType_rf.RData",
                                                  inputSampleType = "Blood Derived Normal", 
                                                  origPlotROC = plotFungiVsBacteriaNumFeat_PerCT_BDN_WGS$plotROC,
                                                  origPlotPR = plotFungiVsBacteriaNumFeat_PerCT_BDN_WGS$plotPR)

#----------------------------------------------------------#
# Subset data
#----------------------------------------------------------#

# Subset to WGS and remove 1 fungi sample with 0 counts: 13722.58cfa82ee4b0c9d6adf6a982
metaFungiVsBacteria <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>%
  filter(cgc_platform == "Illumina HiSeq", experimental_strategy == "WGS") %>% 
  filter(sample_name != "58cfa82ee4b0c9d6adf6a982") %>% droplevels()
countsFungiWGS <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaFungiVsBacteria),]
countsBacteriaWGS <- rep200Data_WGS_RNA_Matched_Bacteria[rownames(metaFungiVsBacteria),]

save(metaFungiVsBacteria,
     countsFungiWGS,
     countsBacteriaWGS,
     file = "Interim_data/data_for_ml_fungi_vs_bacteria_numFeat_28Oct21.RData")

load("Interim_data/data_for_ml_fungi_vs_bacteria_numFeat_28Oct21.RData")

#----------------------------------------------------------#
# Joined RF feature rankings
#----------------------------------------------------------#
load("Interim_data/decontamResultsV2_13Oct21.RData")
load("Interim_data/shared_fungi_features_at_each_taxa_level_13Sep21.RData")
rep200TaxSplit <- read.csv("Supporting_data/rep200_lineage_map_split.csv", stringsAsFactors = FALSE, row.names = 1)
rep200Kingdoms <- read.csv("Supporting_data/rep200_gOTU_kingdom_mapping.csv", stringsAsFactors = FALSE)
rep200Kingdoms_Fungi <- rep200Kingdoms[which(rep200Kingdoms$kingdom == "fungi"),]
rep200Kingdoms_Bacteria <- rep200Kingdoms[which(rep200Kingdoms$kingdom == "bacteria"),]

rep200TaxSplit_Fungi <- rep200TaxSplit[rownames(rep200TaxSplit) %in% rep200Kingdoms_Fungi$genomeID,,drop=FALSE]
dim(rep200TaxSplit_Fungi) # 320   7

rep200TaxSplit_Bacteria <- rep200TaxSplit[rownames(rep200TaxSplit) %in% rep200Kingdoms_Bacteria$genomeID,,drop=FALSE]
dim(rep200TaxSplit_Bacteria) # 11080     7

fungiOGUs <- rownames(rep200TaxSplit_Fungi)
bacteriaOGUs <- rownames(rep200TaxSplit_Bacteria)

rankedJoined_panWGS_PrimaryTumor_rf_10k <- read.csv("Interim_data/rankedJoined_panWGS_PrimaryTumor_rf_10k.csv", stringsAsFactors = FALSE)
which(rankedJoined_panWGS_PrimaryTumor_rf_10k$Feature %in% fungiOGUs)
rankedJoined_panWGS_PrimaryTumor_rf_10k$species <- rep200TaxSplit[rankedJoined_panWGS_PrimaryTumor_rf_10k$Feature,"Species"]
rankedJoined_panWGS_PrimaryTumor_rf_10k$species <- gsub("^s__","",rankedJoined_panWGS_PrimaryTumor_rf_10k$species)
rankedJoined_panWGS_PrimaryTumor_rf_10k$species <- gsub("\\s+","_",rankedJoined_panWGS_PrimaryTumor_rf_10k$species)
indexFungi <- which(rankedJoined_panWGS_PrimaryTumor_rf_10k$Feature %in% fungiOGUs)
tmp <- head(rankedJoined_panWGS_PrimaryTumor_rf_10k[indexFungi,],20)
tmp$decontamV2Decision <- decontamResultsV2[tmp$Feature, "decision"]
tmp$reason <- decontamResultsV2[tmp$Feature, "reason"]
# tmp$inWIS <- tmp$species %in% sharedSpecies
#----------------------------------------------------------#
# Subset data by reads (rarefy)
#----------------------------------------------------------#
require(phyloseq)
require(microbiome)
require(vegan)
# Rarefy # of reads for 10 iterations each time. n=100, 500, 1000, 5000, 10e3, 15e3, 20e3, 25e3, 30e3, 35e3
quantile(rowSums(countsFungiWGS), probs = seq(0.05, 1, by = 0.10))
# 5%      15%      25%      35%      45%      55%      65%      75%      85%      95% 
#   1426.25  3265.25  5075.75  6976.50  9332.50 11995.75 15081.00 18324.75 23191.25 37453.75
quantile(rowSums(countsBacteriaWGS), probs = seq(0.05, 1, by = 0.10))
# 5%       15%       25%       35%       45%       55%       65%       75%       85%       95% 
#   9094.00  23816.00  37344.25  50510.25  65467.75  83646.50 108026.00 143425.75 211967.00 854868.75
#----------------------------Load tax tables for phyloseq----------------------------#
rep200TaxSplit <- read.csv("Supporting_data/rep200_lineage_map_split.csv", stringsAsFactors = FALSE, row.names = 1)
rep200Kingdoms <- read.csv("Supporting_data/rep200_gOTU_kingdom_mapping.csv", stringsAsFactors = FALSE)
rep200Kingdoms_Fungi <- rep200Kingdoms[which(rep200Kingdoms$kingdom == "fungi"),]
rep200Kingdoms_Bacteria <- rep200Kingdoms[which(rep200Kingdoms$kingdom == "bacteria"),]

rep200TaxSplit_Fungi <- rep200TaxSplit[rownames(rep200TaxSplit) %in% rep200Kingdoms_Fungi$genomeID,,drop=FALSE]
dim(rep200TaxSplit_Fungi) # 320   7

rep200TaxSplit_Bacteria <- rep200TaxSplit[rownames(rep200TaxSplit) %in% rep200Kingdoms_Bacteria$genomeID,,drop=FALSE]
dim(rep200TaxSplit_Bacteria) # 11080     7

fungiOGUs <- rownames(rep200TaxSplit_Fungi)
bacteriaOGUs <- rownames(rep200TaxSplit_Bacteria)

save(metaFungiVsBacteria,
     countsFungiWGS,
     countsBacteriaWGS,
     rep200TaxSplit_Fungi,
     rep200TaxSplit_Bacteria,
     file = "Interim_data/data_for_ml_fungi_vs_bacteria_numReads_28Oct21.RData")

load("Interim_data/data_for_ml_fungi_vs_bacteria_numReads_28Oct21.RData")

#----------------------------------------------------------#
# Build function
#----------------------------------------------------------#

fungi_vs_bacteria_numReads_ML <- function(metaData = metaFungiVsBacteria,
                                          dataFungi = countsFungiWGS,
                                          dataBacteria = countsBacteriaWGS,
                                          sampleType = "Primary Tumor",
                                          seqCenter = "Harvard Medical School",
                                          taxTableFungi = rep200TaxSplit_Fungi,
                                          taxTableBacteria = rep200TaxSplit_Bacteria,
                                          modelType = "gbm",
                                          # readDepthVec = c(50, 100, 500, 1000, 5000, 10e3, 15e3, 20e3, 25e3, 30e3, 35e3),
                                          readDepthVec = c(25e3),
                                          numIter = 1,
                                          numResampleIter = 1,
                                          numKFold = 4){
  st <- gsub('([[:punct:]])|\\s+','',sampleType)
  sc <- gsub('([[:punct:]])|\\s+','',seqCenter)
  numberSeeds <- sample(1:1000, numIter, replace = FALSE)
  
  metaDataFilt <- metaData %>% filter(sample_type == sampleType, data_submitting_center_label == seqCenter) %>% droplevels()
  metaDataFilt$predY <- factor(gsub("^TCGA-","",metaDataFilt$investigation))
  dataFungiFilt <- dataFungi[rownames(metaDataFilt),]
  dataBacteriaFilt <- dataBacteria[rownames(metaDataFilt),]
  
  psFungi <- phyloseq(otu_table(dataFungiFilt, taxa_are_rows = FALSE), 
                      tax_table(as.matrix(taxTableFungi)), sample_data(metaDataFilt))
  psBacteria <- phyloseq(otu_table(dataBacteriaFilt, taxa_are_rows = FALSE), 
                         tax_table(as.matrix(taxTableBacteria)), sample_data(metaDataFilt))
  
  fungiDataListPerDepth <- list()
  bacteriaDataListPerDepth <- list()
  fungiDataListAll <- list()
  bacteriaDataListAll <- list()
  for(kk in 1:length(readDepthVec)){
    readDepth <- readDepthVec[kk]
    print(sprintf("Read depth: %d", readDepth))
    for(nn in 1:numIter){
      seed <- numberSeeds[nn]
      print(sprintf("Iteration %d",nn))
      print("Rarefying fungal data...")
      fungiRarefied <- rarefy_even_depth(psFungi, sample.size = readDepth, rngseed = seed, replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
      print("Rarefying bacterial data...")
      bacteriaRarefied <- rarefy_even_depth(psBacteria, sample.size = readDepth, rngseed = seed, replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
      fungiDataListPerDepth[[nn]] <- data.frame(otu_table(fungiRarefied))
      bacteriaDataListPerDepth[[nn]] <- data.frame(otu_table(bacteriaRarefied))
    }
    fungiDataListAll[[kk]] <- fungiDataListPerDepth
    bacteriaDataListAll[[kk]] <- bacteriaDataListPerDepth
  }
  
  perfCombined <- list()
  perfCombinedTmp <- list()
  for(jj in 1:length(readDepthVec)){ # should be length of readDepthVec
    readDepth <- readDepthVec[jj]
    print(sprintf("ML rarefied read depth: %d", readDepth))
    for(ii in 1:numIter){ # should be numIter
      iterNum <- ii
      fungiRarefiedData <- fungiDataListAll[[jj]][[ii]]
      bacteriaRarefiedData <- bacteriaDataListAll[[jj]][[ii]]
      # NOTE: Rarefying at various values will discard samples w/ less reads than the threshold.
      # Thus, to make an even comparison, we need to subset the data to the same samples
      intersectingSamples <- intersect(rownames(fungiRarefiedData),rownames(bacteriaRarefiedData))
      metaDataFiltSubset <- droplevels(metaDataFilt[intersectingSamples,])
      
      if(any(table(metaDataFiltSubset$predY) <= 1)){
        problemCTs <- names(table(metaDataFiltSubset$predY))[which(table(metaDataFiltSubset$predY) <= 1)]
        metaDataFiltSubset <- metaDataFiltSubset %>% filter(!(predY %in% problemCTs)) %>% droplevels()
      }
      
      mlDataY <- metaDataFiltSubset
      mlDataX_Fungi <- fungiRarefiedData[rownames(mlDataY),]
      mlDataX_Bacteria <- bacteriaRarefiedData[rownames(mlDataY),]
      
      # Rarefaction stats
      numIntersectingSamples <- length(intersectingSamples)
      numFeatFungi <- dim(mlDataX_Fungi)[2]
      numFeatBacteria <- dim(mlDataX_Bacteria)[2]
      
      print(sprintf("ML iteration %d",iterNum))
      print(sprintf("Number itersected samples: %d",numIntersectingSamples))
      print(sprintf("Number fungi: %d",numFeatFungi))
      print(sprintf("Number bacteria: %d",numFeatBacteria))
      
      set.seed(42)
      index <- createDataPartition(mlDataY[,"predY"], p = 0.7, list = FALSE)
      trainX_Fungi <- mlDataX_Fungi[index,]
      trainX_Bacteria <- mlDataX_Bacteria[index,]
      trainY <- mlDataY[index,"predY"]
      testX_Fungi <- mlDataX_Fungi[-index,]
      testX_Bacteria <- mlDataX_Bacteria[-index,]
      testY <- mlDataY[-index,"predY"]
      
      print(table(trainY))
      print(table(testY))
      
      set.seed(42) # have to restate seed again, as ctrl defines the cross validation sampling during training
      ctrl <- trainControl(method = "repeatedcv",
                           number = numKFold,
                           repeats = numResampleIter,
                           sampling = "up",
                           summaryFunction = multiClassSummary,
                           classProbs = TRUE,
                           verboseIter = FALSE,
                           savePredictions = TRUE,
                           allowParallel=TRUE)
      
      print("Now training fungi model...")
      set.seed(42)
      mlModel_Fungi <- train(x = trainX_Fungi,
                             y = trainY,
                             method = modelType,
                             preProcess = c("zv"),
                             trControl = ctrl,
                             # metric = "ROC",
                             # tuneGrid = kfoldGBMGrid,
                             verbose = FALSE)
      
      print("Now training bacteria model...")
      set.seed(42)
      mlModel_Bacteria <- train(x = trainX_Bacteria,
                                y = trainY,
                                method = modelType,
                                preProcess = c("zv"),
                                trControl = ctrl,
                                # metric = "ROC",
                                # tuneGrid = kfoldGBMGrid,
                                verbose = FALSE)
      
      print("Obtaining performance values...")
      multiClass_Fungi <- data.frame(obs = testY,
                                     pred = predict(mlModel_Fungi, newdata = testX_Fungi),
                                     predict(mlModel_Fungi, newdata = testX_Fungi, type = "prob"))
      multiClass_Bacteria <- data.frame(obs = testY,
                                        pred = predict(mlModel_Bacteria, newdata = testX_Bacteria),
                                        predict(mlModel_Bacteria, newdata = testX_Bacteria, type = "prob"))
      
      fungiPerf <- multiClassSummary(multiClass_Fungi, lev = levels(multiClass_Fungi$obs))
      bacteriaPerf <- multiClassSummary(multiClass_Bacteria, lev = levels(multiClass_Bacteria$obs))
      
      fungiPerfDf <- data.frame(as.list(fungiPerf))
      colnames(fungiPerfDf) <- paste0("f_",colnames(fungiPerfDf))
      
      bacteriaPerfDf <- data.frame(as.list(bacteriaPerf))
      colnames(bacteriaPerfDf) <- paste0("b_",colnames(bacteriaPerfDf))
      
      perfCombined[[ii]] <- cbind(fungiPerfDf, bacteriaPerfDf, 
                                  readDepth = readDepth,
                                  iter = paste0("iter",iterNum),
                                  numSamples = numIntersectingSamples, 
                                  numCT = length(table(mlDataY$predY)),
                                  numFungi = numFeatFungi, 
                                  numBacteria = numFeatBacteria,
                                  cancerTypes = paste(names(table(mlDataY$predY)), sep = "|"))
    }
    perfCombinedTmp[[jj]] <- do.call(rbind,perfCombined)
    cat('\n')
  }
  
  perfCombinedAll <- do.call(rbind, perfCombinedTmp)
  baseFilename <- paste0("fungi_vs_bacteria_numReads_",st,"_",sc,"_numIter",numIter,"_k",numKFold)
  write.csv(perfCombinedAll, file = paste0(baseFilename,".csv"))
  save(perfCombinedAll, file = paste0(baseFilename,".RData"))
  
}
# fungi_vs_bacteria_numReads_ML()

#----------------------------------------------------------#
# NumReads: Aggregate data across multiple runs
# (Had to be broken up due to memory requirements)
#----------------------------------------------------------#

load("Interim_data/fungi_vs_bacteria_numReads_10000_BloodDerivedNormal_AllSeqCenters_numIter10_k4_modelType_xgbTree.RData", verbose = T)

extractNumReadData <- function(numReads = "50", st = "PrimaryTumor"){
  filePath <- paste0("Interim_data/fungi_vs_bacteria_numReads_",numReads,"_",st,"_AllSeqCenters_numIter10_k4_modelType_xgbTree.RData")
  # Load perfCombinedAll object
  load(filePath)
  tmpPerf <- perfCombinedAll
  rm(perfCombinedAll)
  return(tmpPerf)
}
# PT
perfPT_numReads_50 <- extractNumReadData(numReads = 50)
perfPT_numReads_500 <- extractNumReadData(numReads = 500)
perfPT_numReads_1k <- extractNumReadData(numReads = 1e3)
perfPT_numReads_5k <- extractNumReadData(numReads = 5e3)
perfPT_numReads_10k <- extractNumReadData(numReads = 10e3)
perfPT_numReads_15k <- extractNumReadData(numReads = 15e3)
perfPT_numReads_20k <- extractNumReadData(numReads = 20e3)
perfPT_numReads_25k <- extractNumReadData(numReads = 25e3)
perfPT_numReads_30k <- extractNumReadData(numReads = 30e3)
perfPT_numReads_35k <- extractNumReadData(numReads = 35e3)
perfPT_numReads_All <- rbind(perfPT_numReads_50, perfPT_numReads_500, perfPT_numReads_1k,
                             perfPT_numReads_5k, perfPT_numReads_10k, perfPT_numReads_15k,
                             perfPT_numReads_20k, perfPT_numReads_25k, perfPT_numReads_30k, perfPT_numReads_35k)

keepCols <- c("AUC","prAUC")
perfPT_numReads_AllFormatted <- perfPT_numReads_All %>% select(-cancerTypes) %>% distinct() %>%
  reshape2::melt(id.vars = c("iter","readDepth","numSamples","numCT","numFungi","numBacteria")) %>%
  mutate(Domain = factor(ifelse(grepl("f_",variable),yes = "Fungi",no = "Bacteria"), levels = c("Fungi","Bacteria"))) %>%
  mutate(variable = gsub("^f_|^b_","",variable)) %>% filter(variable %in% keepCols) %>% filter(!is.na(value))
maxDepth <- max(perfPT_numReads_AllFormatted$readDepth)
perfPT_numReads_AllFormatted %>%
  filter(variable == "AUC") %>%
  ggline(x = "readDepth",
         y = "value",
         palette = "nejm",
         color = "Domain",
         add = "mean_ci",
         xlab = "Number of rarefied reads in pan-cancer model",
         ylab = "Average AUROC",
         numeric.x.axis = TRUE) +
  scale_x_continuous(breaks = seq(0, maxDepth, length.out = 5), limits = c(0,maxDepth)) +
  stat_compare_means(aes(group=Domain), label = "p.signif", size = 3)
# BDN
# perfBDN_numReads_50 <- extractNumReadData(numReads = 50, st = "BloodDerivedNormal")
# perfBDN_numReads_500 <- extractNumReadData(numReads = 500, st = "BloodDerivedNormal")
# perfBDN_numReads_1k <- extractNumReadData(numReads = 1e3, st = "BloodDerivedNormal")
# perfBDN_numReads_5k <- extractNumReadData(numReads = 5e3, st = "BloodDerivedNormal")
# perfBDN_numReads_10k <- extractNumReadData(numReads = 10e3, st = "BloodDerivedNormal")
# perfBDN_numReads_15k <- extractNumReadData(numReads = 15e3, st = "BloodDerivedNormal")
# perfBDN_numReads_20k <- extractNumReadData(numReads = 20e3, st = "BloodDerivedNormal")
# perfBDN_numReads_25k <- extractNumReadData(numReads = 25e3, st = "BloodDerivedNormal")
# perfBDN_numReads_30k <- extractNumReadData(numReads = 30e3, st = "BloodDerivedNormal")
# perfBDN_numReads_35k <- extractNumReadData(numReads = 35e3, st = "BloodDerivedNormal")
# perfBDN_numReads_All <- rbind(perfBDN_numReads_50, perfBDN_numReads_500, perfBDN_numReads_1k,
#                              perfBDN_numReads_5k, perfBDN_numReads_10k, perfBDN_numReads_15k,
#                              perfBDN_numReads_20k, perfBDN_numReads_25k, perfBDN_numReads_30k, perfBDN_numReads_35k)

#----------------------------------------------------------#
# NumReads: Plot results per seq center
#----------------------------------------------------------#

plotFungiVsBacteriaNumReads <- function(rDataFilePath, 
                                       seqCenter = "Harvard Medical School", 
                                       inputSampleType="Primary Tumor",
                                       statSize = 2.5){
  if(inputSampleType=="Primary Tumor"){st = "PT"}
  if(inputSampleType=="Primary Tumor vs Solid Tissue Normal"){st = "PT_vs_NAT"}
  if(inputSampleType=="Blood Derived Normal"){st = "BDN"}
  seqCenterFormatted <- gsub('([[:punct:]])|\\s+','',seqCenter)
  load(rDataFilePath) # loads perfCombinedAll object
  keepCols <- c("AUC","prAUC")
  perfCombinedAllFormatted <- perfCombinedAll %>% select(-cancerTypes) %>% distinct() %>%
    reshape2::melt(id.vars = c("iter","readDepth","numSamples","numCT","numFungi","numBacteria")) %>%
    mutate(Domain = factor(ifelse(grepl("f_",variable),yes = "Fungi",no = "Bacteria"), levels = c("Fungi","Bacteria"))) %>%
    mutate(variable = gsub("^f_|^b_","",variable)) %>% filter(variable %in% keepCols) %>% filter(!is.na(value))
  maxDepth <- max(perfCombinedAllFormatted$readDepth)
  perfCombinedAllFormatted %>%
    filter(variable == "AUC") %>%
    ggline(x = "readDepth",
           y = "value",
           palette = "nejm",
           color = "Domain",
           add = "mean_ci",
           xlab = "Number of rarefied reads in pan-cancer model",
           ylab = "Average AUROC",
           numeric.x.axis = TRUE) +
    scale_x_continuous(breaks = seq(0, maxDepth, length.out = 5), limits = c(0,maxDepth)) +
    stat_compare_means(aes(group=Domain), label = "p.signif", size = statSize) -> plotROC

  perfCombinedAllFormatted %>%
    filter(variable == "prAUC") %>%
    ggline(x = "readDepth",
           y = "value",
           palette = "nejm",
           color = "Domain",
           add = "mean_ci",
           xlab = "Number of rarefied reads in pan-cancer model",
           ylab = "Average AUPR",
           numeric.x.axis = TRUE) +
    scale_x_continuous(breaks = seq(0, maxDepth, length.out = 5), limits = c(0,maxDepth)) +
    stat_compare_means(aes(group=Domain), label = "p.signif", size = statSize) -> plotPR

  combinedPlotTitle <- paste0("TCGA Simulations: Fungi vs. Bacteria (varying number of rarefied reads)\n",
                              inputSampleType," | ",seqCenter,"\n(10 iterations of multi-class classification per point)")
  combinedPlot <- ggarrange(plotROC, plotPR, ncol = 2, common.legend = TRUE)
  combinedPlotAnnotated <- annotate_figure(combinedPlot,
                                           top = text_grob(combinedPlotTitle,
                                                           color = "black", face = "bold", size = 14))
  print(combinedPlotAnnotated)
  ggsave(filename = paste0("Figures/Supplementary_Figures/fungi_vs_bacteria_numReads",seqCenterFormatted,"_",st,".jpeg"),
         plot = combinedPlotAnnotated,
         dpi = "retina", units = "in", width = 12, height = 5)
  rm(perfCombinedAll)
  # return(perfCombinedAllFormatted)
}

plotFungiVsBacteriaNumReads(rDataFilePath = "Interim_data/fungi_vs_bacteria_numReads_PrimaryTumor_HarvardMedicalSchool_numIter5_k4_modelType_rf.RData", 
                           seqCenter = "Harvard Medical School", statSize=3)
plotFungiVsBacteriaNumReads(rDataFilePath = "Interim_data/fungi_vs_bacteria_numReads_PrimaryTumor_BroadInstituteofMITandHarvard_numIter5_k4_modelType_rf.RData", 
                           seqCenter = "Broad Institute", statSize=3)
plotFungiVsBacteriaNumReads(rDataFilePath = "Interim_data/fungi_vs_bacteria_numReads_PrimaryTumor_BaylorCollegeofMedicine_numIter5_k4_modelType_rf.RData",
                           seqCenter = "Baylor College of Medicine", statSize=3)
plotFungiVsBacteriaNumReads(rDataFilePath = "Interim_data/fungi_vs_bacteria_numReads_PrimaryTumor_MDAndersonInstituteforAppliedCancerScience_numIter5_k4_modelType_rf.RData",
                           seqCenter = "MD Anderson", statSize=3)
plotFungiVsBacteriaNumReads(rDataFilePath = "Interim_data/fungi_vs_bacteria_numReads_PrimaryTumor_WashingtonUniversitySchoolofMedicine_numIter5_k4_modelType_rf.RData",
                           seqCenter = "WashU", statSize=3)







