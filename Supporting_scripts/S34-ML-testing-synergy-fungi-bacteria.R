#-----------------------------------------------------------------------------
# S34-ML-testing-synergy-fungi-bacteria.R
# Copyright (c) 2021--, Greg Poore
# Purpose: Using TCGA data, compare WIS-intersected fungi vs. bacteria vs. fungi+bacteria ML performance
#-----------------------------------------------------------------------------

# NOTE: Getting this function to run on the cluster took a bit of finagling:
# First start an interactive job: srun --nodes=1 --cpus-per-task=64 --mem=64gb --time=24:00:00 --pty bash -i
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

## Import data
load("data_for_ml_tcga_synergy_fungi_and_bacteria_18Nov21.RData")

ml_synergy_fungiBacteria <- function(metaData = metaFungiVsBacteria,
                                         dataFungi = countsFungiWGS_Shared,
                                         dataBacteria = countsBacteriaWGS_Shared,
                                         dataFungiBacteria = countsFungiBacteriaWGS_Shared,
                                         sampleType = "Primary Tumor",
                                         modelType = "xgbTree",
                                         # wisIndex = indexWISfungi,
                                         # nonwisIndex = inverseIndexNonWISfungi,
                                         # numFeatVec = 34,
                                         numIter = 5,
                                         numResampleIter = 1,
                                         numKFold = 4){
  st <- gsub('([[:punct:]])|\\s+','',sampleType)
  sc <- "AllSeqCenters"
  # wisFungiSamplingList <- list()
  # nonwisFungiSamplingList <- list()
  # for(zz in 1:length(numFeatVec)){
  #   numSamp <- numFeatVec[zz]
  #   wisFungiSamplingList[[zz]] <- replicate(numIter, expr = sample(wisIndex, size = numSamp, replace = FALSE))
  #   nonwisFungiSamplingList[[zz]] <- replicate(numIter, expr = sample(nonwisIndex, size = numSamp, replace = FALSE))
  # }
  
  metaDataFilt <- metaData %>% filter(sample_type == sampleType) %>% droplevels()
  metaDataFilt$predY <- factor(gsub("^TCGA-","",metaDataFilt$investigation))
  dataFungiFilt <- dataFungi[rownames(metaDataFilt),]
  dataBacteriaFilt <- dataBacteria[rownames(metaDataFilt),]
  dataFungiBacteriaFilt <- dataFungiBacteria[rownames(metaDataFilt),]
  
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
  for(ii in 1:numIter){ # should be 1:numIter
    iterNum <- ii
    mlDataY <- metaDataFilt

    mlDataX_fungi <- dataFungiFilt
    mlDataX_bacteria <- dataBacteriaFilt
    mlDataX_fungiBacteria <- dataFungiBacteriaFilt


    # mlDataX_wisFungi <- dataFungiFilt[,wisfungiFeat]
    # mlDataX_nonwisFungi <- dataFungiFilt[,nonwisFungiFeat]
    
    print(sprintf("Iteration %d",iterNum))
    
    indexAllSeqCenterTrain <- indexAllSeqCenterTrainList[[iterNum]]
    indexAllSeqCenterTest <- indexAllSeqCenterTestList[[iterNum]]

    # Define training sets
    trainX_fungi <- mlDataX_fungi[indexAllSeqCenterTrain,]
    trainX_bacteria <- mlDataX_bacteria[indexAllSeqCenterTrain,]
    trainX_fungiBacteria <- mlDataX_fungiBacteria[indexAllSeqCenterTrain,]
    trainY <- mlDataY[indexAllSeqCenterTrain,"predY"]

    # Define testing sets
    testX_fungi <- mlDataX_fungi[indexAllSeqCenterTest,]
    testX_bacteria <- mlDataX_bacteria[indexAllSeqCenterTest,]
    testX_fungiBacteria <- mlDataX_fungiBacteria[indexAllSeqCenterTest,]
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
                         # tuneGrid = kfoldGBMGrid,
                         verbose = FALSE)
  
  print("Now training bacteria model...")
  set.seed(42)
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
                            # tuneGrid = kfoldGBMGrid,
                            verbose = FALSE)

  print("Now training fungi+bacteria model...")
  set.seed(42)
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
                            # tuneGrid = kfoldGBMGrid,
                            verbose = FALSE)
    
    
    print("Obtaining performance values...")
    multiClass_fungi <- data.frame(obs = testY,
                                   pred = predict(mlModel_fungi, newdata = testX_fungi),
                                   predict(mlModel_fungi, newdata = testX_fungi, type = "prob"))
    multiClass_bacteria <- data.frame(obs = testY,
                                      pred = predict(mlModel_bacteria, newdata = testX_bacteria),
                                      predict(mlModel_bacteria, newdata = testX_bacteria, type = "prob"))
    multiClass_fungiBacteria <- data.frame(obs = testY,
                                      pred = predict(mlModel_fungiBacteria, newdata = testX_fungiBacteria),
                                      predict(mlModel_fungiBacteria, newdata = testX_fungiBacteria, type = "prob"))

    fungiPerf <- multiClassSummary(multiClass_fungi, lev = levels(multiClass_fungi$obs))
    bacteriaPerf <- multiClassSummary(multiClass_bacteria, lev = levels(multiClass_bacteria$obs))
    fungiBacteriaPerf <- multiClassSummary(multiClass_fungiBacteria, lev = levels(multiClass_fungiBacteria$obs))
    
    fungiPerfDf <- data.frame(as.list(fungiPerf))
    colnames(fungiPerfDf) <- paste0("f_",colnames(fungiPerfDf))

    bacteriaPerfDf <- data.frame(as.list(bacteriaPerf))
    colnames(bacteriaPerfDf) <- paste0("b_",colnames(bacteriaPerfDf))

    fungiBacteriaPerfDf <- data.frame(as.list(fungiBacteriaPerf))
    colnames(fungiBacteriaPerfDf) <- paste0("all_",colnames(fungiBacteriaPerfDf))
    
    perfCombined[[ii]] <- cbind(fungiPerfDf, bacteriaPerfDf, fungiBacteriaPerfDf, iter = paste0("iter",iterNum))
    if((iterNum %% 10) == 0){
      checkpointPerfCombined <- perfCombined
      save(checkpointPerfCombined, file = paste0("checkpoints/checkpoint_synergy_iterNum",iterNum,"_",st,"_",sc,"_k",numKFold,"_modelType_",modelType,".RData"))
    }
    cat('\n')
  }
  
  perfCombinedAll <- do.call(rbind, perfCombined)
  baseFilename <- paste0("tcga_synergy_fungi_bacteria_",st,"_",sc,"_numIter",numIter,"_k",numKFold,"_modelType_",modelType)
  write.csv(perfCombinedAll, file = paste0(baseFilename,".csv"))
  save(perfCombinedAll, file = paste0(baseFilename,".RData"))
  
  res <- list(perfCombinedAll=perfCombinedAll,
              indexAllSeqCenterTrainList=indexAllSeqCenterTrainList,
              indexAllSeqCenterTestList=indexAllSeqCenterTestList)
  return(res)
}

require(splitstackshape)

tmp_PT3 <- ml_synergy_fungiBacteria(sampleType = "Primary Tumor", numIter = 3)

require(splitstackshape)
# Primary Tumor
tmp_PT <- ml_synergy_fungiBacteria(sampleType = "Primary Tumor", numIter = 100)
# Blood Derived Normal
tmp_BDN <- ml_synergy_fungiBacteria(sampleType = "Blood Derived Normal", numIter = 100)


