#-----------------------------------------------------------------------------
# S29-calculate-ranked-features-fungi-with-bacteria-train-set-only-PT.R
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

calcRankedFeaturesJoinedTrainSetOnly <- function(metaData = metaFungiVsBacteria,
                                dataFungi = countsFungiWGS,
                                dataBacteria = countsBacteriaWGS,
                                sampleType = "Primary Tumor",
                                modelType = "rf",
                                numKFold = 10,
                                numResampleIter = 1){
  st <- gsub('([[:punct:]])|\\s+','',sampleType)
  sc <- "AllSeqCenters"

  metaDataFilt <- metaData %>% filter(sample_type == sampleType) %>% droplevels()
  metaDataFilt$predY <- factor(gsub("^TCGA-","",metaDataFilt$investigation))
  dataFungiFilt <- dataFungi[rownames(metaDataFilt),]
  dataBacteriaFilt <- dataBacteria[rownames(metaDataFilt),]

  mlDataY <- metaDataFilt
  mlDataX_Fungi <- dataFungiFilt
  mlDataX_Bacteria <- dataBacteriaFilt

  print("Dimensions of full fungal data set:")
  print(dim(mlDataX_Fungi))

  indexAllSeqCenterTrain <- which(rownames(mlDataY) %in% rownames(split1MetaFungiVsBacteria))
  indexAllSeqCenterTest <- which(rownames(mlDataY) %in% rownames(split2MetaFungiVsBacteria))

  trainX_Fungi <- mlDataX_Fungi[indexAllSeqCenterTrain,]
  trainX_Bacteria <- mlDataX_Bacteria[indexAllSeqCenterTrain,]
  trainY <- mlDataY[indexAllSeqCenterTrain,"predY"]
  trainX_Joined <- cbind(trainX_Fungi, trainX_Bacteria)

  print("Dimensions of joined training data and labels:")
  print(dim(trainX_Joined))
  print(length(trainY))

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
    print("Feature ranking with joined fungi and bacteria using RFs...")
    set.seed(42)
    mlModel_Joined <- train(x = trainX_Joined,
                           y = trainY,
                           method = modelType,
                           preProcess = c("zv"),
                           trControl = ctrl,
                           verbose = FALSE)
    # Extract feature importances for rf model
    varJoinedImpBestModelDF <- as.data.frame(varImp( mlModel_Joined$finalModel, scale = FALSE ))
    varJoinedImpBestModelDF2 <- rownames_to_column(varJoinedImpBestModelDF, "Feature")
    varJoinedImpBestModelDF2Ordered <- varJoinedImpBestModelDF2[order(-varJoinedImpBestModelDF2$Overall),]
    colnames(varJoinedImpBestModelDF2Ordered)[2] <- "varImp"
    varJoinedImpBestModelDF2OrderedNonzero <- varJoinedImpBestModelDF2Ordered[varJoinedImpBestModelDF2Ordered$varImp != 0,]
    write.csv(varJoinedImpBestModelDF2OrderedNonzero, file = paste0("rankedJoined_panWGS_trainingSetOnly_",st,"_",modelType,"_10k.csv"), row.names = FALSE)
    print("Feature ranking complete!")
    rm(mlModel_Joined)
    } else if(modelType == "xgbTree"){
      #--------------------------------xgbTree model--------------------------------#
      print("Feature ranking with joined fungi and bacteria using GBMs...")
      set.seed(42)
      mlModel_Joined <- train(x = trainX_Joined,
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
      varJoinedImpBestModelDF <- data.frame(xgb.importance(model = mlModel_Joined$finalModel))
      varJoinedImpBestModelDFOrdered <- varJoinedImpBestModelDF[order(-varJoinedImpBestModelDF$Gain),]
      write.csv(varJoinedImpBestModelDFOrdered, file = paste0("rankedJoined_panWGS_trainingSetOnly_",st,"_",modelType,"_10k.csv"), row.names = FALSE)
      print("Feature ranking complete!")
      rm(mlModel_Joined)
    }

}

calcRankedFeaturesJoinedTrainSetOnly(sampleType = "Primary Tumor", modelType = "rf")
calcRankedFeaturesJoinedTrainSetOnly(sampleType = "Blood Derived Normal", modelType = "rf")
#--------------------------------------------------------------------------------------------#

calcRankedFeaturesTrainSetOnly <- function(metaData = metaFungiVsBacteria,
                                dataFungi = countsFungiWGS,
                                dataBacteria = countsBacteriaWGS,
                                sampleType = "Primary Tumor",
                                modelType = "rf",
                                numKFold = 10,
                                numResampleIter = 1){
  st <- gsub('([[:punct:]])|\\s+','',sampleType)
  sc <- "AllSeqCenters"

  metaDataFilt <- metaData %>% filter(sample_type == sampleType) %>% droplevels()
  metaDataFilt$predY <- factor(gsub("^TCGA-","",metaDataFilt$investigation))
  dataFungiFilt <- dataFungi[rownames(metaDataFilt),]
  dataBacteriaFilt <- dataBacteria[rownames(metaDataFilt),]

  mlDataY <- metaDataFilt
  mlDataX_Fungi <- dataFungiFilt
  mlDataX_Bacteria <- dataBacteriaFilt

  print("Dimensions of full fungal data set:")
  print(dim(mlDataX_Fungi))

  indexAllSeqCenterTrain <- which(rownames(mlDataY) %in% rownames(split1MetaFungiVsBacteria))
  indexAllSeqCenterTest <- which(rownames(mlDataY) %in% rownames(split2MetaFungiVsBacteria))

  trainX_Fungi <- mlDataX_Fungi[indexAllSeqCenterTrain,]
  trainX_Bacteria <- mlDataX_Bacteria[indexAllSeqCenterTrain,]
  trainY <- mlDataY[indexAllSeqCenterTrain,"predY"]

  print("Dimensions of joined training data (fungi, bacteria, then labels):")
  print(dim(trainX_Fungi))
  print(dim(trainX_Bacteria))
  print(length(trainY))

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
    write.csv(varFungiImpBestModelDF2OrderedNonzero, file = paste0("rankedFungi_panWGS_trainingSetOnly_",st,"_",modelType,"_10k.csv"), row.names = FALSE)
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
    write.csv(varBacteriaImpBestModelDF2OrderedNonzero, file = paste0("rankedBacteria_panWGS_trainingSetOnly_",st,"_",modelType,"_10k.csv"), row.names = FALSE)

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
      write.csv(varFungiImpBestModelDFOrdered, file = paste0("rankedFungi_panWGS_trainingSetOnly_",st,"_",modelType,"_10k.csv"), row.names = FALSE)

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
      write.csv(varBacteriaImpBestModelDFOrdered, file = paste0("rankedBacteria_panWGS_trainingSetOnly_",st,"_",modelType,"_10k.csv"), row.names = FALSE)
      print("Feature ranking complete!")
      rm(mlModel_Fungi)
      rm(mlModel_Bacteria)
    }
}

calcRankedFeaturesTrainSetOnly(sampleType = "Primary Tumor", modelType = "rf")
calcRankedFeaturesTrainSetOnly(sampleType = "Blood Derived Normal", modelType = "rf")

