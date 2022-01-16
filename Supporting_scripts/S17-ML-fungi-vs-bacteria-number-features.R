#-----------------------------------------------------------------------------
# S17-ML-fungi-vs-bacteria-number-features.R
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
load("data_for_ml_fungi_vs_bacteria_numFeat_28Oct21.RData")

fungi_vs_bacteria_numFeat_ML = function(metaData = metaFungiVsBacteria,
                                         dataFungi = countsFungiWGS,
                                         dataBacteria = countsBacteriaWGS,
                                         sampleType = "Primary Tumor",
                                         seqCenter = "Harvard Medical School",
                                         modelType = "rf",
                                         numFeatVec = c(10,25,50,75,100,125,150,175,200,225,250,275,300),
                                         numIter = 100,
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
    for(ii in 1:numIter){ # should be numIter
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
      if((iterNum %% 250) == 0){
        checkpointPerfCombined <- perfCombined
        save(checkpointPerfCombined, file = paste0("checkpoint_numFeat",numFeatTmp,"_numIter",numIter,"_",st,"_",sc,"_k",numKFold))
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
fungi_vs_bacteria_numFeat_ML(sampleType = "Primary Tumor", seqCenter = "Harvard Medical School")
fungi_vs_bacteria_numFeat_ML(sampleType = "Primary Tumor", seqCenter = "Broad Institute of MIT and Harvard")
fungi_vs_bacteria_numFeat_ML(sampleType = "Primary Tumor", seqCenter = "Baylor College of Medicine")
fungi_vs_bacteria_numFeat_ML(sampleType = "Primary Tumor", seqCenter = "MD Anderson - Institute for Applied Cancer Science")
fungi_vs_bacteria_numFeat_ML(sampleType = "Primary Tumor", seqCenter = "Washington University School of Medicine")
# Blood Derived Normal
fungi_vs_bacteria_numFeat_ML(sampleType = "Blood Derived Normal", seqCenter = "Harvard Medical School")
fungi_vs_bacteria_numFeat_ML(sampleType = "Blood Derived Normal", seqCenter = "Broad Institute of MIT and Harvard")
fungi_vs_bacteria_numFeat_ML(sampleType = "Blood Derived Normal", seqCenter = "Baylor College of Medicine")
fungi_vs_bacteria_numFeat_ML(sampleType = "Blood Derived Normal", seqCenter = "MD Anderson - Institute for Applied Cancer Science")
fungi_vs_bacteria_numFeat_ML(sampleType = "Blood Derived Normal", seqCenter = "Washington University School of Medicine")




