#-----------------------------------------------------------------------------
# S27-ML-fungi-vs-bacteria-number-features-per-cancer-type-ranked-PT.R
# Copyright (c) 2021--, Greg Poore
# Purpose: Run machine learning on all available TCGA cancer types varying number of ranked taxa
#-----------------------------------------------------------------------------

#-------------------------------#
# Load dependencies
require(plyr)
require(caret) # for model building
require(dplyr) # for data manipulation
require(doMC) # for parallel computing
require(gbm) # for machine learning
# require(xgboost) # for machine learning
# require(randomForest) # for machine learning
require(tibble) # for df operations
require(PRROC) # for precision-recall curves
# require(MLmetrics) # for multi-class learning

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
load("data_for_ml_fungi_vs_bacteria_numFeat_with_AllSeqCenter_index_29Oct21.RData")

# c(10,25,50,75,100,125,150,175,200,225,250,275,300)

fungi_vs_bacteria_rankedNumFeat_perCT_ML = function(metaData = metaFungiVsBacteria,
                                         dataFungi = countsFungiWGS,
                                         dataBacteria = countsBacteriaWGS,
                                         sampleType = "Primary Tumor",
                                         modelType = "gbm",
                                         rankedFeatModelType = "rf",
                                         numFeatVec = c(10,25,50,75,100,125,150,175,200,225,250,275,300),
                                         numIter = 1,
                                         numResampleIter = 1,
                                         numKFold = 4){
  kfoldGBMGrid <- data.frame(n.trees=150, interaction.depth=3,
                                         shrinkage=0.1,
                                         n.minobsinnode=1)
  st <- gsub('([[:punct:]])|\\s+','',sampleType)
  sc <- "AllSeqCenters"

  fungiFeatRanked <- read.csv(file = paste0("rankedFungi_panWGS_",st,"_",rankedFeatModelType,"_10k.csv"), stringsAsFactors = FALSE, header = TRUE)
  bacteriaFeatRanked <- read.csv(file = paste0("rankedBacteria_panWGS_",st,"_",rankedFeatModelType,"_10k.csv"), stringsAsFactors = FALSE, header = TRUE)

  fungiSamplingList <- list()
  bacteriaSamplingList <- list()
  for(zz in 1:length(numFeatVec)){
    numSamp <- numFeatVec[zz]
    fungiSamplingList[[zz]] <- fungiFeatRanked[1:numSamp,"Feature"]
    bacteriaSamplingList[[zz]] <- bacteriaFeatRanked[1:numSamp,"Feature"]
  }

  metaDataFilt <- droplevels(metaData[metaData$sample_type == sampleType,])
  dataFungiFilt <- dataFungi[rownames(metaDataFilt),]
  dataBacteriaFilt <- dataBacteria[rownames(metaDataFilt),]

  metaDataFilt$investigationFormatted <- factor(gsub("^TCGA-","",metaDataFilt$investigation))
  split1MetaFungiVsBacteria$investigationFormatted <- factor(gsub("^TCGA-","",split1MetaFungiVsBacteria$investigation))
  split2MetaFungiVsBacteria$investigationFormatted <- factor(gsub("^TCGA-","",split2MetaFungiVsBacteria$investigation))

  # Extract disease types with 20 or more samples in training set
  dzVec <- names(table(split1MetaFungiVsBacteria$investigationFormatted))[unname(table(split1MetaFungiVsBacteria$investigationFormatted))>=20]
  
  perfCombined <- list()
  perfCombinedTmp2 <- list()
  perfCombinedTmp3 <- list()
  for(jj in 1:length(numFeatVec)){ # should be 1:length(numFeatVec)
    fungiListTmp <- fungiSamplingList[[jj]]
    bacteriaListTmp <- bacteriaSamplingList[[jj]]
    numFeatTmp <- numFeatVec[jj]
    print(sprintf("Number of features: %d", numFeatTmp))

    for(kk in 1:length(dzVec)){
      dz <- dzVec[kk]
      metaDataFilt$predY <- factor(ifelse(metaDataFilt$investigationFormatted == dz, yes = dz, no = "OtherCancerType"), levels = c(dz,"OtherCancerType"))
      positiveClass <- dz
      negativeClass <- "OtherCancerType"

      print(sprintf("Cancer type: %s (%d / %d)",dz,kk,length(dzVec)))

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
        
        set.seed(42) # have to restate seed again, as ctrl defines the cross validation sampling during training
        ctrl <- trainControl(method = "repeatedcv",
                             number = numKFold,
                             repeats = numResampleIter,
                             sampling = "up",
                             summaryFunction = twoClassSummary,
                             classProbs = TRUE,
                             verboseIter = FALSE,
                             savePredictions = TRUE,
                             allowParallel=TRUE)
        
        print("Now training fungi model...")
        set.seed(42)
        mlModel_Fungi <- train(x = trainX_Fungi,
                               y = trainY,
                               method = modelType,
                               trControl = ctrl,
                               preProcess = c("zv"),
                               metric = "ROC",
                               tuneGrid = kfoldGBMGrid,
                               verbose = FALSE)
        
        print("Now training bacteria model...")
        set.seed(42)
        mlModel_Bacteria <- train(x = trainX_Bacteria,
                                  y = trainY,
                                  method = modelType,
                                  trControl = ctrl,
                                  preProcess = c("zv"),
                                  metric = "ROC",
                                  tuneGrid = kfoldGBMGrid,
                                  verbose = FALSE)
        
        print("Obtaining performance values...")
        # Fungi
        predProbsFungi <- as.numeric(predict(mlModel_Fungi, newdata = testX_Fungi, type = "prob")[,positiveClass])
        fgFungi <- predProbsFungi[testY == positiveClass]
        bgFungi <- predProbsFungi[testY == negativeClass]
        prroc_roc_Fungi <- roc.curve(scores.class0 = fgFungi, scores.class1 = bgFungi, curve = T)
        prroc_pr_Fungi <- pr.curve(scores.class0 = fgFungi, scores.class1 = bgFungi, curve = T, rand.compute=T)
        # Bacteria
        predProbsBacteria <- as.numeric(predict(mlModel_Bacteria, newdata = testX_Bacteria, type = "prob")[,positiveClass])
        fgBacteria <- predProbsBacteria[testY == positiveClass]
        bgBacteria <- predProbsBacteria[testY == negativeClass]
        prroc_roc_Bacteria <- roc.curve(scores.class0 = fgBacteria, scores.class1 = bgBacteria, curve = T)
        prroc_pr_Bacteria <- pr.curve(scores.class0 = fgBacteria, scores.class1 = bgBacteria, curve = T, rand.compute=T)
        # Format into dataframe
        perfCombined[[ii]] <- data.frame(f_auroc = prroc_roc_Fungi$auc,
          f_aupr = prroc_pr_Fungi$auc.integral,
          b_auroc = prroc_roc_Bacteria$auc,
          b_aupr = prroc_pr_Bacteria$auc.integral,
          diseaseType = dz,
          sampleType = st,
          seqCenter = sc,
          iter = iterNum,
          rankedFeatureModelType = rankedFeatModelType,
          numFeat = numFeatTmp)

        # if((iterNum %% 50) == 0){
        #   checkpointPerfCombined <- perfCombined
        #   save(checkpointPerfCombined, file = paste0("checkpoints/checkpoint_perCT_rankedNumFeat",numFeatTmp,"_iterNum",
        #     iterNum,"_",st,"_",sc,"_k",numKFold,"_modelType_",modelType,"_rankedFeatModelType_",rankedFeatModelType,".RData"))
        # }
      }
      perfCombinedTmp2[[kk]] <- do.call(rbind,perfCombined)
    }
    
    perfCombinedTmp3[[jj]] <- do.call(rbind,perfCombinedTmp2)
    cat('\n')
  }
  
  perfCombinedAll <- do.call(rbind, perfCombinedTmp3)
  baseFilename <- paste0("fungi_vs_bacteria_perCT_rankedNumFeat_",paste(numFeatVec, collapse="_"),"_",st,"_",sc,
    "_numIter",numIter,"_k",numKFold,"_modelType_",modelType,"_rankedFeatModelType_",rankedFeatModelType)
  write.csv(perfCombinedAll, file = paste0(baseFilename,".csv"))
  save(perfCombinedAll, fungiSamplingList, bacteriaSamplingList, file = paste0(baseFilename,".RData"))
}
# Primary Tumor
fungi_vs_bacteria_rankedNumFeat_perCT_ML(sampleType = "Blood Derived Normal")


