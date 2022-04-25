#-----------------------------------------------------------------------------
# S34-ML-testing-synergy-fungi-bacteria.R
# Copyright (c) 2021--, Greg Poore
# Purpose: Using TCGA data, compare WIS-intersected fungi vs. bacteria vs. fungi+bacteria ML performance
#-----------------------------------------------------------------------------

# NOTE: Getting this function to run on the cluster took a bit of finagling:
# First start an interactive job: srun --nodes=1 --cpus-per-task=64 --mem=64gb --time=24:00:00 --pty bash -i
# Then open a r4 conda environment.
# Then interactively run the code in script S17B-ML-fungi-vs-bacteria-number-features.R for a few iterations
# Stop that function, load the data in this script, load the function, and test for a few iterations.
# If the first model iteration takes <1 min, it should run to completion

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

load("data_for_ml_fungi_vs_bacteria_numFeat_with_AllSeqCenter_index_5Apr22.RData")

fungi_vs_bacteria_numFeat_ML = function(metaData = metaFungiVsBacteria,
                                         dataFungi = countsFungiWGS,
                                         dataBacteria = countsBacteriaWGS,
                                         sampleType = "Primary Tumor",
                                         # seqCenter = "Harvard Medical School",
                                         modelType = "xgbTree",
                                         numFeatVec = c(10),
                                         numIter = 2,
                                         numResampleIter = 1,
                                         numKFold = 2){
  st <- gsub('([[:punct:]])|\\s+','',sampleType)
  sc <- "AllSeqCenters"
  fungiSamplingList <- list()
  bacteriaSamplingList <- list()
  for(zz in 1:length(numFeatVec)){
    numSamp <- numFeatVec[zz]
    fungiSamplingList[[zz]] <- replicate(numIter, expr = sample(1:ncol(dataFungi), size = numSamp, replace = FALSE))
    bacteriaSamplingList[[zz]] <- replicate(numIter, expr = sample(1:ncol(dataBacteria), size = numSamp, replace = FALSE))
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
      
      # set.seed(42)
      # index <- createDataPartition(mlDataY[,"predY"], p = 0.7, list = FALSE)

      indexAllSeqCenterTrain <- which(rownames(mlDataY) %in% rownames(split1MetaFungiVsBacteria))
      indexAllSeqCenterTest <- which(rownames(mlDataY) %in% rownames(split2MetaFungiVsBacteria))


      trainX_Fungi <- mlDataX_Fungi[indexAllSeqCenterTrain,]
      trainX_Bacteria <- mlDataX_Bacteria[indexAllSeqCenterTrain,]
      trainY <- mlDataY[indexAllSeqCenterTrain,"predY"]
      testX_Fungi <- mlDataX_Fungi[indexAllSeqCenterTest,]
      testX_Bacteria <- mlDataX_Bacteria[indexAllSeqCenterTest,]
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
        # save(checkpointPerfCombined, file = paste0("checkpoints/checkpoint_numFeat",numFeatTmp,"_iterNum",iterNum,"_",st,"_",sc,"_k",numKFold,"_modelType_",modelType,".RData"))
      }
    }
    perfCombinedTmp[[jj]] <- do.call(rbind,perfCombined)
    cat('\n')
  }
  
  perfCombinedAll <- do.call(rbind, perfCombinedTmp)
  baseFilename <- paste0("fungi_vs_bacteria_numFeat_",st,"_",sc,"_numIter",numIter,"_k",numKFold,"_modelType_",modelType)
  # write.csv(perfCombinedAll, file = paste0(baseFilename,".csv"))
  # save(perfCombinedAll, fungiSamplingList, bacteriaSamplingList, file = paste0(baseFilename,".RData"))
}





fungi_vs_bacteria_numFeat_ML(sampleType = "Blood Derived Normal")






## Import data
load("data_for_ml_tcga_synergy_fungi_and_bacteria_5Apr22.RData")

ml_synergy_fungiBacteria_CV <- function(metaData = metaFungiVsBacteria,
                                         dataFungi = countsFungiWGS_Shared,
                                         dataBacteria = countsBacteriaWGS_Shared,
                                         dataFungiBacteria = countsFungiBacteriaWGS_Shared,
                                         sampleType = "Primary Tumor",
                                         modelType = "xgbTree",
                                         numIter = 5,
                                         numResampleIter = 1,
                                         numKFold = 10){
  st <- gsub('([[:punct:]])|\\s+','',sampleType)
  sc <- "AllSeqCenters"
  
  metaDataFilt <- metaData %>% filter(sample_type == sampleType) %>% droplevels()
  metaDataFilt$predY <- factor(gsub("^TCGA-","",metaDataFilt$investigation))
  dataFungiFilt <- dataFungi[rownames(metaDataFilt),]
  dataBacteriaFilt <- dataBacteria[rownames(metaDataFilt),]
  dataFungiBacteriaFilt <- dataFungiBacteria[rownames(metaDataFilt),]

  xgbGrid <- data.frame(nrounds = 10,
                       max_depth = 4,
                       eta = .1,
                       gamma = 0,
                       colsample_bytree = .7,
                       min_child_weight = 1,
                       subsample = .8)
  
  #----------------------------
  
  perfCombined <- list()
  perfCombinedTmp <- list()
  for(ii in 1:numIter){ # should be 1:numIter
    iterNum <- ii
    mlDataY <- metaDataFilt

    mlDataX_fungi <- dataFungiFilt
    mlDataX_bacteria <- dataBacteriaFilt
    mlDataX_fungiBacteria <- dataFungiBacteriaFilt
    
    print(sprintf("Iteration %d",iterNum))

    # Define training sets
    trainX_fungi <- mlDataX_fungi
    trainX_bacteria <- mlDataX_bacteria
    trainX_fungiBacteria <- mlDataX_fungiBacteria
    trainY <- mlDataY[,"predY"]
    
    set.seed(ii) # have to restate seed again, as ctrl defines the cross validation sampling during training
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
  set.seed(ii)
  mlModel_fungi <- train(x = trainX_fungi,
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
                         tuneGrid = xgbGrid,
                         verbose = FALSE)
  
  print("Now training bacteria model...")
  set.seed(ii)
  mlModel_bacteria <- train(x = trainX_bacteria,
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
                            tuneGrid = xgbGrid,
                            verbose = FALSE)

  print("Now training fungi+bacteria model...")
  set.seed(ii)
  mlModel_fungiBacteria <- train(x = trainX_fungiBacteria,
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
                            tuneGrid = xgbGrid,
                            verbose = FALSE)

  resPredFun <- function(mlModel, type=NA){
    resPred <- mlModel$pred

    ## Split folds and calculate perf on each fold
    resPredSplit <- split(resPred, resPred$Resample)
    repX_perf <- list()
    for(zz in seq_along(resPredSplit)){
      resPredSingleRep <- resPredSplit[[zz]]
      predProbs <- resPredSingleRep
      multiClass <- resPredSingleRep
      rep_perfTmp <- multiClassSummary(multiClass, lev = levels(multiClass$obs))
      repX_perf[[zz]] <- data.frame(as.list(rep_perfTmp))

      if(type == "fungi"){
        colnames(repX_perf[[zz]]) <- paste0("f_",colnames(repX_perf[[zz]]))
      } else if(type == "bacteria"){
        colnames(repX_perf[[zz]]) <- paste0("b_",colnames(repX_perf[[zz]]))
      } else if(type == "fungi+bacteria"){
        colnames(repX_perf[[zz]]) <- paste0("all_",colnames(repX_perf[[zz]]))
      }

    }
    
    # SUMMARIZE MODEL PERFORMANCES
    rep_perf <- do.call(rbind, repX_perf)
    # print(rep_perf)
    return(rep_perf)
  }

  print("Obtaining performance values...")
  fungiPerfDf <- resPredFun(mlModel_fungi, type="fungi")
  bacteriaPerfDf <- resPredFun(mlModel_bacteria, type="bacteria")
  fungiBacteriaPerfDf <- resPredFun(mlModel_fungiBacteria, type="fungi+bacteria")
    
    perfCombined[[ii]] <- cbind(fungiPerfDf, bacteriaPerfDf, fungiBacteriaPerfDf, iter = paste0("iter",iterNum))
    if((iterNum %% 10) == 0){
      checkpointPerfCombined <- perfCombined
      save(checkpointPerfCombined, file = paste0("checkpoints/checkpoint_CV_synergy_iterNum",iterNum,"_",st,"_",sc,"_k",numKFold,"_modelType_",modelType,".RData"))
    }
    cat('\n')
  }
  
  perfCombinedAll <- do.call(rbind, perfCombined)
  baseFilename <- paste0("tcga_CV_synergy_fungi_bacteria_",st,"_",sc,"_numIter",numIter,"_k",numKFold,"_modelType_",modelType)
  write.csv(perfCombinedAll, file = paste0(baseFilename,".csv"))
  save(perfCombinedAll, file = paste0(baseFilename,".RData"))
  
  res <- list(perfCombinedAll=perfCombinedAll)
  return(res)
}

require(splitstackshape)
# Blood Derived Normal
tmp_BDN <- ml_synergy_fungiBacteria_CV(sampleType = "Blood Derived Normal", numIter = 50)


