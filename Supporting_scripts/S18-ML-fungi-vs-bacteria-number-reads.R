#-----------------------------------------------------------------------------
# S18-ML-fungi-vs-bacteria-number-reads.R
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
require(phyloseq) # for rarefaction

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
load("data_for_ml_fungi_vs_bacteria_numReads_28Oct21.RData")

fungi_vs_bacteria_numReads_ML <- function(metaData = metaFungiVsBacteria,
                                          dataFungi = countsFungiWGS,
                                          dataBacteria = countsBacteriaWGS,
                                          sampleType = "Primary Tumor",
                                          seqCenter = "Harvard Medical School",
                                          taxTableFungi = rep200TaxSplit_Fungi,
                                          taxTableBacteria = rep200TaxSplit_Bacteria,
                                          modelType = "rf",
                                          readDepthVec = c(50, 100, 500, 1000, 5000, 10e3, 15e3, 20e3, 25e3, 30e3, 35e3),
                                          numIter = 5,
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

      # Remove any cancer types with only 1 sample available after rarefaction (problematic for train/test splits)
      if(any(table(metaDataFiltSubset$predY) <= 1)){
        problemCTs <- names(table(metaDataFiltSubset$predY))[which(table(metaDataFiltSubset$predY) <= 1)]
        metaDataFiltSubset <- metaDataFiltSubset %>% filter(!(predY %in% problemCTs)) %>% droplevels()
      }
      # Skip if fewer than 2 cancer types available after rarefaction
      if(length(table(metaDataFiltSubset$predY)) < 2){next}
      
      # Skip if all cancer types have fewer than 15 samples after rarefaction
      if(all(table(metaDataFiltSubset$predY) < 15)){next}
      
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
  baseFilename <- paste0("fungi_vs_bacteria_numReads_",st,"_",sc,"_numIter",numIter,"_k",numKFold,"_modelType_",modelType)
  write.csv(perfCombinedAll, file = paste0(baseFilename,".csv"))
  save(perfCombinedAll, file = paste0(baseFilename,".RData")) 
}
# Primary Tumor
# fungi_vs_bacteria_numReads_ML(sampleType = "Primary Tumor", seqCenter = "Harvard Medical School")
fungi_vs_bacteria_numReads_ML(sampleType = "Primary Tumor", seqCenter = "Broad Institute of MIT and Harvard")
fungi_vs_bacteria_numReads_ML(sampleType = "Primary Tumor", seqCenter = "Baylor College of Medicine")
fungi_vs_bacteria_numReads_ML(sampleType = "Primary Tumor", seqCenter = "MD Anderson - Institute for Applied Cancer Science")
fungi_vs_bacteria_numReads_ML(sampleType = "Primary Tumor", seqCenter = "Washington University School of Medicine")
# Blood Derived Normal
fungi_vs_bacteria_numReads_ML(sampleType = "Blood Derived Normal", seqCenter = "Harvard Medical School")
fungi_vs_bacteria_numReads_ML(sampleType = "Blood Derived Normal", seqCenter = "Broad Institute of MIT and Harvard")
fungi_vs_bacteria_numReads_ML(sampleType = "Blood Derived Normal", seqCenter = "Baylor College of Medicine")
fungi_vs_bacteria_numReads_ML(sampleType = "Blood Derived Normal", seqCenter = "MD Anderson - Institute for Applied Cancer Science")
fungi_vs_bacteria_numReads_ML(sampleType = "Blood Derived Normal", seqCenter = "Washington University School of Medicine")




