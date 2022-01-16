#-----------------------------------------------------------------------------
# S28-calculate-ranked-features-fungi-with-bacteria-PT.R
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
                                modelType = "rf",
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

  trainX_Joined <- cbind(trainX_Fungi, trainX_Bacteria)

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
    print("Feature ranking with joined fungi and bacteria...")
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
    write.csv(varJoinedImpBestModelDF2OrderedNonzero, file = paste0("rankedJoined_panWGS_",st,"_",modelType,"_10k.csv"), row.names = FALSE)
    print("Feature ranking complete!")
    rm(mlModel_Joined)
    } else if(modelType == "xgbTree"){
      #--------------------------------xgbTree model--------------------------------#
      print("Feature ranking with joined fungi and bacteria...")
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
      write.csv(varJoinedImpBestModelDFOrdered, file = paste0("rankedJoined_panWGS_",st,"_",modelType,"_10k.csv"), row.names = FALSE)
      print("Feature ranking complete!")
      rm(mlModel_Joined)
    }

}

calcRankedFeatures(sampleType = "Primary Tumor", modelType = "rf")