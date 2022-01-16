#-----------------------------------------------------------------------------
# S24-ML-WIS-fungi-vs-bacteria-number-features.R
# Copyright (c) 2021--, Greg Poore
# Purpose: Run machine learning on all available TCGA cancer types by seq center
#-----------------------------------------------------------------------------

#-------------------------------#
# Load dependencies
require(tidyr)
require(caret) # for model building
require(dplyr) # for data manipulation
require(doMC) # for parallel computing
require(gbm) # for machine learning
require(xgboost) # for machine learning
require(randomForest) # for machine learning
require(tibble) # for df operations
require(PRROC) # for precision-recall curves
require(MLmetrics) # for multi-class learning

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
load("data_for_WIS_fungi_vs_bacteria_04Nov21.RData")

wis_fungi_vs_bacteria_numFeat_ML <- function(metaData = metaWzFungiBacteriaFree_Species,
                                         dataFungi = wzFungiOnlyFreeCounts_Species,
                                         dataBacteria = wzBacteriaOnlyFreeCounts_Species,
                                         sampleType = "tumor",
                                         modelType = "xgbTree",
                                         numFeatVec = c(10,25,50,75,100,125,150,175,200,225,250,275,300),
                                         numIter = 100,
                                         numResampleIter = 1,
                                         numKFold = 4){
  st <- gsub('([[:punct:]])|\\s+','',sampleType)
  sc <- "WIS"
  fungiSamplingList <- list()
  bacteriaSamplingList <- list()
  for(zz in 1:length(numFeatVec)){
    numSamp <- numFeatVec[zz]
    fungiSamplingList[[zz]] <- replicate(numIter, expr = sample(1:ncol(dataFungi), size = numSamp, replace = FALSE))
    bacteriaSamplingList[[zz]] <- replicate(numIter, expr = sample(1:ncol(dataBacteria), size = numSamp, replace = FALSE))
  }
  
  metaDataFilt <- metaData %>% filter(sample_type == sampleType) %>% droplevels()
  metaDataFilt$predY <- factor(gsub('([[:punct:]])|\\s+','',metaDataFilt$tissue))
  dataFungiFilt <- dataFungi[rownames(metaDataFilt),]
  dataBacteriaFilt <- dataBacteria[rownames(metaDataFilt),]
  
  perfCombined <- list()
  perfCombinedTmp <- list()
  for(jj in 1:length(numFeatVec)){ # should be 1:length(numFeatVec)
    fungiListTmp <- fungiSamplingList[[jj]]
    bacteriaListTmp <- bacteriaSamplingList[[jj]]
    numFeatTmp <- dim(fungiListTmp)[1]
    print(sprintf("Number of features: %d", numFeatTmp))
    for(ii in 1:numIter){ # should be 1:numIter
      iterNum <- ii
      fungiFeat <- fungiListTmp[,iterNum]
      bacteriaFeat <- bacteriaListTmp[,iterNum]
      mlDataY <- metaDataFilt
      mlDataX_Fungi <- dataFungiFilt[,fungiFeat]
      mlDataX_Bacteria <- dataBacteriaFilt[,bacteriaFeat]
      
      print(sprintf("Iteration %d",iterNum))
      
      set.seed(42)
      index <- createDataPartition(mlDataY[,"predY"], p = 0.7, list = FALSE)

      # indexAllSeqCenterTrain <- which(rownames(mlDataY) %in% rownames(split1MetaFungiVsBacteria))
      # indexAllSeqCenterTest <- which(rownames(mlDataY) %in% rownames(split2MetaFungiVsBacteria))


      trainX_Fungi <- mlDataX_Fungi[index,]
      trainX_Bacteria <- mlDataX_Bacteria[index,]
      trainY <- mlDataY[index,"predY"]
      testX_Fungi <- mlDataX_Fungi[-index,]
      testX_Bacteria <- mlDataX_Bacteria[-index,]
      testY <- mlDataY[-index,"predY"]
      
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
      
      print("Now training fungi model...")
      set.seed(42)
      mlModel_Fungi <- train(x = trainX_Fungi,
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
      
      print("Now training bacteria model...")
      set.seed(42)
      mlModel_Bacteria <- train(x = trainX_Bacteria,
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
      if((iterNum %% 10) == 0){
        checkpointPerfCombined <- perfCombined
        save(checkpointPerfCombined, file = paste0("checkpoints/checkpoint_numFeat",numFeatTmp,"_iterNum",iterNum,"_",st,"_",sc,"_k",numKFold,"_modelType_",modelType,".RData"))
      }
    }
    perfCombinedTmp[[jj]] <- do.call(rbind,perfCombined)
    cat('\n')
  }
  
  perfCombinedAll <- do.call(rbind, perfCombinedTmp)
  baseFilename <- paste0("fungi_vs_bacteria_numFeat_",st,"_",sc,"_numIter",numIter,"_k",numKFold,"_modelType_",modelType)
  write.csv(perfCombinedAll, file = paste0(baseFilename,".csv"))
  save(perfCombinedAll, fungiSamplingList, bacteriaSamplingList, file = paste0(baseFilename,".RData"))
}
# Primary Tumor
wis_fungi_vs_bacteria_numFeat_ML(sampleType = "tumor")


