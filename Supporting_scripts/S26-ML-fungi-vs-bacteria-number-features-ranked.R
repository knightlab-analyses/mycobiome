#-----------------------------------------------------------------------------
# S26-ML-fungi-vs-bacteria-number-features.R
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
load("data_for_ml_fungi_vs_bacteria_numFeat_with_AllSeqCenter_index_29Oct21.RData")

#----------------------------------------------------------#
# Determine feature rankings
#----------------------------------------------------------#

calcRankedFeatures <- function(metaData = metaFungiVsBacteria,
                                dataFungi = countsFungiWGS,
                                dataBacteria = countsBacteriaWGS,
                                sampleType = "Primary Tumor",
                                modelType = "xgbTree",
                                numKFold = 10,
                                numResampleIter = 1){
  st <- gsub('([[:punct:]])|\\s+','',sampleType)
  metaDataFilt <- metaData %>% filter(sample_type == sampleType) %>% droplevels()
  metaDataFilt$predY <- factor(gsub("^TCGA-","",metaDataFilt$investigation))
  dataFungiFilt <- dataFungi[rownames(metaDataFilt),]
  dataBacteriaFilt <- dataBacteria[rownames(metaDataFilt),]

  trainX_Fungi <- dataFungiFilt
  trainX_Bacteria <- dataBacteriaFilt
  trainY <- metaDataFilt[,"predY"]

  if(modelType == "rf"){
    #--------------------------------RF model--------------------------------#

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
    print("Feature ranking fungi...")
    set.seed(42)
    mlModel_Fungi <- train(x = trainX_Fungi,
                           y = trainY,
                           method = modelType,
                           preProcess = c("zv"),
                           trControl = ctrl,
                           verbose = FALSE)
    # Extract feature importances for rf model
    varFungiImpBestModelDF <- as.data.frame(varImp( mlModel_Fungi$finalModel, scale = FALSE ))
    varFungiImpBestModelDF2 <- rownames_to_column(varFungiImpBestModelDF, "Feature")
    varFungiImpBestModelDF2Ordered <- varFungiImpBestModelDF2[order(-varFungiImpBestModelDF2$Overall),]
    colnames(varFungiImpBestModelDF2Ordered)[2] <- "varImp"
    varFungiImpBestModelDF2OrderedNonzero <- varFungiImpBestModelDF2Ordered[varFungiImpBestModelDF2Ordered$varImp != 0,]
    write.csv(varFungiImpBestModelDF2OrderedNonzero, file = paste0("rankedFungi_panWGS_",st,"_",modelType,"_10k.csv"), row.names = FALSE)
    print("Feature ranking bacteria...")
    set.seed(42)
    mlModel_Bacteria <- train(x = trainX_Bacteria,
                           y = trainY,
                           method = modelType,
                           preProcess = c("zv"),
                           trControl = ctrl,
                           verbose = FALSE)
    # Extract feature importances for rf model
    varBacteriaImpBestModelDF <- as.data.frame(varImp( mlModel_Bacteria$finalModel, scale = FALSE ))
    varBacteriaImpBestModelDF2 <- rownames_to_column(varBacteriaImpBestModelDF, "Feature")
    varBacteriaImpBestModelDF2Ordered <- varBacteriaImpBestModelDF2[order(-varBacteriaImpBestModelDF2$Overall),]
    colnames(varBacteriaImpBestModelDF2Ordered)[2] <- "varImp"
    varBacteriaImpBestModelDF2OrderedNonzero <- varBacteriaImpBestModelDF2Ordered[varBacteriaImpBestModelDF2Ordered$varImp != 0,]
    write.csv(varBacteriaImpBestModelDF2OrderedNonzero, file = paste0("rankedBacteria_panWGS_",st,"_",modelType,"_10k.csv"), row.names = FALSE)

    rm(mlModel_Fungi)
    rm(mlModel_Bacteria)
    } else if(modelType == "xgbTree"){
      #--------------------------------xgbTree model--------------------------------#
      print("Feature ranking fungi...")
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
                             verbose = FALSE)
      # Extract feature importances for xgbTree model
      varFungiImpBestModelDF <- data.frame(xgb.importance(model = mlModel_Fungi$finalModel))
      varFungiImpBestModelDFOrdered <- varFungiImpBestModelDF[order(-varFungiImpBestModelDF$Gain),]
      write.csv(varFungiImpBestModelDFOrdered, file = paste0("rankedFungi_panWGS_",st,"_",modelType,"_10k.csv"), row.names = FALSE)

      print("Feature ranking bacteria...")
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
                                verbose = FALSE)

      # Extract feature importances for xgbTree model
      varBacteriaImpBestModelDF <- data.frame(xgb.importance(model = mlModel_Bacteria$finalModel))
      varBacteriaImpBestModelDFOrdered <- varBacteriaImpBestModelDF[order(-varBacteriaImpBestModelDF$Gain),]
      write.csv(varBacteriaImpBestModelDFOrdered, file = paste0("rankedBacteria_panWGS_",st,"_",modelType,"_10k.csv"), row.names = FALSE)
      print("Feature ranking complete!")
      rm(mlModel_Fungi)
      rm(mlModel_Bacteria)
    }

}

#--------------------------------Function for running feature rankings--------------------------------#

fungi_vs_bacteria_rankedNumFeat_ML = function(metaData = metaFungiVsBacteria,
                                         dataFungi = countsFungiWGS,
                                         dataBacteria = countsBacteriaWGS,
                                         sampleType = "Primary Tumor",
                                         modelType = "xgbTree",
                                         rankedFeatModelType = "rf",
                                         basenameRankedFungiFeatures = "rankedFungi_panWGS_trainingSetOnly",
                                         basenameRankedBacteriaFeatures = "rankedBacteria_panWGS_trainingSetOnly",
                                         numFeatVec = c(10,25,50,75,100,125,150,175,200,225,250,275,300),
                                         # numFeatVec = c(10,25),
                                         numIter = 1,
                                         numResampleIter = 1,
                                         numKFold = 4){
  st <- gsub('([[:punct:]])|\\s+','',sampleType)
  sc <- "AllSeqCenters"

  fungiFeatRanked <- read.csv(file = paste0(basenameRankedFungiFeatures,"_",st,"_",rankedFeatModelType,"_10k.csv"), stringsAsFactors = FALSE, header = TRUE)
  bacteriaFeatRanked <- read.csv(file = paste0(basenameRankedBacteriaFeatures,"_",st,"_",rankedFeatModelType,"_10k.csv"), stringsAsFactors = FALSE, header = TRUE)

  fungiSamplingList <- list()
  bacteriaSamplingList <- list()
  for(zz in 1:length(numFeatVec)){
    numSamp <- numFeatVec[zz]
    fungiSamplingList[[zz]] <- fungiFeatRanked[1:numSamp,"Feature"]
    bacteriaSamplingList[[zz]] <- bacteriaFeatRanked[1:numSamp,"Feature"]
  }
  
  metaDataFilt <- metaData %>% filter(sample_type == sampleType) %>% droplevels()
  metaDataFilt$predY <- factor(gsub("^TCGA-","",metaDataFilt$investigation))
  dataFungiFilt <- dataFungi[rownames(metaDataFilt),]
  dataBacteriaFilt <- dataBacteria[rownames(metaDataFilt),]
  
  perfCombined <- list()
  perfCombinedTmp <- list()
  for(jj in 1:length(numFeatVec)){ # should be 1:length(numFeatVec)
    fungiListTmp <- fungiSamplingList[[jj]]
    bacteriaListTmp <- bacteriaSamplingList[[jj]]
    numFeatTmp <- numFeatVec[jj]
    print(sprintf("Number of features: %d", numFeatTmp))
    for(ii in 1:numIter){ # should be 1:numIter
      iterNum <- ii
      fungiFeat <- fungiListTmp
      bacteriaFeat <- bacteriaListTmp
      mlDataY <- metaDataFilt
      mlDataX_Fungi <- dataFungiFilt[,fungiFeat]
      mlDataX_Bacteria <- dataBacteriaFilt[,bacteriaFeat]
      
      print(sprintf("Iteration %d",iterNum))

      indexAllSeqCenterTrain <- which(rownames(mlDataY) %in% rownames(split1MetaFungiVsBacteria))
      indexAllSeqCenterTest <- which(rownames(mlDataY) %in% rownames(split2MetaFungiVsBacteria))


      trainX_Fungi <- mlDataX_Fungi[indexAllSeqCenterTrain,]
      trainX_Bacteria <- mlDataX_Bacteria[indexAllSeqCenterTrain,]
      trainY <- mlDataY[indexAllSeqCenterTrain,"predY"]
      testX_Fungi <- mlDataX_Fungi[indexAllSeqCenterTest,]
      testX_Bacteria <- mlDataX_Bacteria[indexAllSeqCenterTest,]
      testY <- mlDataY[indexAllSeqCenterTest,"predY"]

      print("Dimensions of training set (using fungi):")
      print(dim(trainX_Fungi))
      print("Dimensions of test set (using fungi):")
      print(dim(testX_Fungi))
      
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
        save(checkpointPerfCombined, file = paste0("checkpoints/checkpoint_trainingSetOnly_rankedNumFeat",numFeatTmp,
          "_iterNum",iterNum,"_",st,"_",sc,"_k",numKFold,"_modelType_",modelType,"_rankedFeatModelType_",rankedFeatModelType,".RData"))
      }

      rm(mlModel_Fungi)
      rm(mlModel_Bacteria)
    }
    perfCombinedTmp[[jj]] <- do.call(rbind,perfCombined)
    cat('\n')
  }
  
  perfCombinedAll <- do.call(rbind, perfCombinedTmp)
  baseFilename <- paste0("fungi_vs_bacteria_trainingSetOnly_rankedNumFeat_",st,"_",sc,"_numIter",
    numIter,"_k",numKFold,"_modelType_",modelType,"_rankedFeatModelType_",rankedFeatModelType)
  write.csv(perfCombinedAll, file = paste0(baseFilename,".csv"))
  save(perfCombinedAll, fungiSamplingList, bacteriaSamplingList, file = paste0(baseFilename,".RData"))
}

# Primary Tumor
# calcRankedFeatures(sampleType = "Primary Tumor", modelType = "rf")
# calcRankedFeatures(sampleType = "Primary Tumor", modelType = "xgbTree")
fungi_vs_bacteria_rankedNumFeat_ML(sampleType = "Primary Tumor", rankedFeatModelType = "rf")
fungi_vs_bacteria_rankedNumFeat_ML(sampleType = "Primary Tumor", rankedFeatModelType = "xgbTree", numFeatVec = c(10,25,50,75,100,125,150,175,200,225,250))


# Blood Derived Normal
# calcRankedFeatures(sampleType = "Blood Derived Normal", modelType = "rf")
# calcRankedFeatures(sampleType = "Blood Derived Normal", modelType = "xgbTree")
fungi_vs_bacteria_rankedNumFeat_ML(sampleType = "Blood Derived Normal", rankedFeatModelType = "rf")
fungi_vs_bacteria_rankedNumFeat_ML(sampleType = "Blood Derived Normal", rankedFeatModelType = "xgbTree", numFeatVec = c(10,25,50,75,100,125,150,175,200,225,250))


