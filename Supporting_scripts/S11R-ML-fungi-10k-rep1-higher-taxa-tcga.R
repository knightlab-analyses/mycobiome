#-----------------------------------------------------------------------------
# S04B-ML-fungi-10k-rep1-higher-taxa-tcga.R
# Copyright (c) 2021--, Greg Poore
# Purpose: Run machine learning on all available TCGA cancer types
#-----------------------------------------------------------------------------

#-------------------------------#
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

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
load("data_for_ml_tcga_decontamV2_2Apr22.RData")

# GBM HYPERPARAMETER SEARCH GRID BELOW (DEFAULT PER CARET PACKAGE)
kfoldGBMGrid <- data.frame(n.trees=150, interaction.depth=3,
                                         shrinkage=0.1,
                                         n.minobsinnode=1)

# Model setup -- MODIFY AS NEEDED; ALL THESE COMPARISONS WILL BE TESTED
sampleTypeList <- c("Primary Tumor vs Solid Tissue Normal",
                    "Blood Derived Normal",
                    "Primary Tumor")
# MODIFY AS NEEDED
datasetList <- list(rep200FungiDecontamV2OrderVSNM,
     rep200FungiDecontamV2FamilyVSNM,
     rep200FungiDecontamV2GenusVSNM,
     rep200FungiDecontamV2SpeciesVSNM,
     snmDataOGUFungiDecontamV2)

# MODIFY AS NEEDED
datasetListNames <- c("rep200FungiDecontamV2OrderVSNM",
     "rep200FungiDecontamV2FamilyVSNM",
     "rep200FungiDecontamV2GenusVSNM",
     "rep200FungiDecontamV2SpeciesVSNM",
     "snmDataOGUFungiDecontamV2")

# MATCHED METADATA DATA FRAME FOR ALL COUNT DATA:
metaTmpQC <- metaQiitaCombined_Nonzero_DecontamV2 # MODIFY AS NEEDED
metaTmpPath <- metaQiitaCombined_Nonzero_DecontamV2 # MODIFY AS NEEDED / PATHOLOGY METADATA DATA FRAME
metaTmpQC$disease_type <- factor(metaTmpQC$disease_type)
caretTuneGrid <- kfoldGBMGrid
numKFold <- 10
numResampleIter <- 1
prroc_roc <- list()
prroc_pr <- list()
perf <- list()
rep_perf <- list()
perfTmp <- list()
rep_perfTmp <- list()
perfTmp2 <- list()
rep_perfTmp2 <- list()

for(jj in seq_along(datasetList)){

  dataTmp <- datasetList[[jj]]
  datasetName <- datasetListNames[[jj]]

  for(kk in seq_along(sampleTypeList)){
    st <- sampleTypeList[kk]
    print(st)
    
    if(st == "Primary Tumor vs Solid Tissue Normal"){
      metaTmp <- metaTmpQC
      metaTmp2 <- droplevels(metaTmp[metaTmp$sample_type %in% c("Primary Tumor",
                                                                "Solid Tissue Normal"),])
    } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
      metaTmp <- metaTmpQC
      metaTmp2 <- droplevels(metaTmp[metaTmp$sample_type %in% st,])
    } else if(st == "Stage I vs IV"){
      metaTmp <- metaTmpPath
      metaTmp2 <- droplevels(metaTmp[(metaTmp$sample_type %in% "Primary Tumor") &
                                       (metaTmp$pathologic_stage_label_binned %in% c("Stage1","Stage4")),])
    }
    
    print(seq_along(levels(metaTmp2$disease_type)))
    
    for(ii in seq_along(levels(metaTmp2$disease_type))){
      
      dt <- levels(metaTmp2$disease_type)[ii]
      print(dt)
      
      if(st == "Primary Tumor vs Solid Tissue Normal"){
        metaTmp3 <- droplevels(metaTmp2[metaTmp2$disease_type == dt,])
        metaTmp3$predY <- factor(gsub('([[:punct:]])|\\s+','',metaTmp3$sample_type))
        positiveClass <- "PrimaryTumor"
        negativeClass <- "SolidTissueNormal"
      } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
        metaTmp3 <- metaTmp2
        metaTmp3$predY <- factor(ifelse(metaTmp2$disease_type == dt, 
                                        yes = dt, 
                                        no = "OtherCancerType"),
                                 levels = c(dt, "OtherCancerType"))
        positiveClass <- gsub('([[:punct:]])|\\s+','',dt)
        negativeClass <- "OtherCancerType"
      } else if(st == "Stage I vs IV"){
        metaTmp3 <- droplevels(metaTmp2[metaTmp2$disease_type == dt,])
        metaTmp3$predY <- metaTmp3$pathologic_stage_label_binned
        positiveClass <- "Stage4"
        negativeClass <- "Stage1"
      }
      
      print(table(metaTmp3$predY))
      
      # SKIP CANCER TYPES THAT ONLY HAVE ONE CLASS OF A COMPARISON
      if(length(table(metaTmp3$predY)) < 2){next}
      
      # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF 20 SAMPLES IN EITHER CLASS
      if(any(table(metaTmp3$predY) < 20)){next}
      
      minorityClassSize <- min(table((metaTmp3$predY)))
      majorityClassSize <- max(table((metaTmp3$predY)))
      
      minorityClassName <- names(table(metaTmp3$predY)[which(table(metaTmp3$predY) == min(table(metaTmp3$predY)))])
      majorityClassName <- names(table(metaTmp3$predY)[which(table(metaTmp3$predY) == max(table(metaTmp3$predY)))])
      
      mlDataY <- metaTmp3
      mlDataX <- dataTmp[rownames(mlDataY),]
      
      # USE 70% OF DATA FOR TRAINING AND 30% FOR TESTING
      set.seed(42)
      # index <- createDataPartition(mlDataY$predY, p = 0.7, list = FALSE)
      trainX <- mlDataX
      trainY <- mlDataY[,"predY"]
      refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
      
      set.seed(42)
      ctrl <- trainControl(method = "repeatedcv",
                           number = numKFold,
                           repeats = numResampleIter,
                           sampling = "up",
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE,
                           verboseIter = TRUE,
                           savePredictions = TRUE,
                           allowParallel=TRUE)
      
      mlModel <- train(x = trainX,
                       y = refactoredTrainY,
                       method = "gbm",
                       preProcess = c("scale","center"),
                       trControl = ctrl,
                       verbose = TRUE,
                       metric = "ROC",
                       tuneGrid = caretTuneGrid)

      resPred <- mlModel$pred
      
      ## Calculate performance on concatenated fold predictions
      predProbs <- resPred
      multiClass <- resPred
      multiClass$pred <- relevel(multiClass$pred, positiveClass)
      multiClass$obs <- relevel(multiClass$obs, positiveClass)
      fg <- predProbs[resPred$obs == positiveClass,positiveClass]
      bg <- predProbs[resPred$obs == negativeClass,positiveClass]

      confusionMatrix <- confusionMatrix(multiClass$pred, multiClass$obs)
      prroc_roc[[ii]] <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
      prroc_pr[[ii]] <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
      rocCurveData <- cbind(as.data.frame(prroc_roc[[ii]]$curve), disease_type = dt, sample_type = st)
      prCurveData <- cbind(as.data.frame(prroc_pr[[ii]]$curve), disease_type = dt, sample_type = st)

      ## Estimate AUROC CIs using cvAUC on concatenated fold predictions
      require(cvAUC)
      out95 <- ci.cvAUC(predictions = resPred[,positiveClass], labels = resPred$obs, label.ordering = c(negativeClass,positiveClass), folds = resPred$Resample, confidence = 0.95)
      out99 <- ci.cvAUC(predictions = resPred[,positiveClass], labels = resPred$obs, label.ordering = c(negativeClass,positiveClass), folds = resPred$Resample, confidence = 0.99)
      # resCvAUC <- data.frame(estimate = c(out95$cvAUC, out99$cvAUC), se = c(out95$se, out99$se), lowerCI = c(out95$ci[1], out99$ci[1]), upperCI = c(out95$ci[1], out99$ci[2]), levelCI = c(0.95,0.99))

      ## Split folds and calculate perf on each fold
      resPredSplit <- split(resPred, resPred$Resample)
      repX_perf <- list()
      for(zz in seq_along(resPredSplit)){
        resPredSingleRep <- resPredSplit[[zz]]
        predProbs <- resPredSingleRep
        multiClass <- resPredSingleRep
        multiClass$pred <- relevel(multiClass$pred, positiveClass)
        multiClass$obs <- relevel(multiClass$obs, positiveClass)
        fg <- predProbs[resPredSingleRep$obs == positiveClass,positiveClass]
        bg <- predProbs[resPredSingleRep$obs == negativeClass,positiveClass]
        
        rep_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
        rep_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
        
        repX_perf[[zz]] <- data.frame(auroc=rep_roc$auc,
                               aupr=rep_pr$auc.integral,
                               rep=paste0("Fold",zz), 
                               diseaseType = dt,
                               sampleType = st,
                               datasetName = datasetName,
                               minorityClassSize = minorityClassSize,
                               majorityClassSize = majorityClassSize,
                               minorityClassName = minorityClassName,
                               majorityClassName = majorityClassName)
      }
      
      # SUMMARIZE MODEL PERFORMANCES
      rep_perf[[ii]] <- do.call(rbind, repX_perf)
      perf[[ii]] <- data.frame(auroc = prroc_roc[[ii]]$auc,
                               aupr = prroc_pr[[ii]]$auc.integral,
                               aucEstimate = out95$cvAUC, # either out95 or out99 work (same result)
                               aucSE95 = out95$se,
                               lowerCI95 = out95$ci[1],
                               upperCI95 = out95$ci[2],
                               aucSE99 = out99$se,
                               lowerCI99 = out99$ci[1],
                               upperCI99 = out99$ci[2],
                               diseaseType = dt,
                               sampleType = st,
                               datasetName = datasetName,
                               minorityClassSize = minorityClassSize,
                               majorityClassSize = majorityClassSize,
                               minorityClassName = minorityClassName,
                               majorityClassName = majorityClassName)
      
      print(perf[[ii]])
      print(confusionMatrix)
      
      #--------------------------------------#
      # Save performance into relevant files #
      #--------------------------------------#
      
      filepathPerfPlots <- paste0("./perfPlots__",datasetName)
      filepathPerfPlotsDataPR <- paste0("./dataPR__",datasetName)
      filepathPerfPlotsDataROC <- paste0("./dataROC__",datasetName)
      filepathFeatures <- paste0("./features__",datasetName)
      filepathPerfStats <- paste0("./perfData__",datasetName)
      filepathPerfStatsPerFold <- paste0("./perfDataPerFold__",datasetName)
      filepathConfusionMatrix <- paste0("./confusionMatrix__",datasetName)
      
      if(!( dir.exists( file.path(filepathPerfPlots)))){
        dir.create(file.path(filepathPerfPlots))
      }
      if(!( dir.exists( file.path(filepathPerfPlotsDataPR)))){
        dir.create(file.path(filepathPerfPlotsDataPR))
      }
      if(!( dir.exists( file.path(filepathPerfPlotsDataROC)))){
        dir.create(file.path(filepathPerfPlotsDataROC))
      }
      if(!( dir.exists( file.path(filepathFeatures)))){
        dir.create(file.path(filepathFeatures))
      }
      if(!( dir.exists( file.path(filepathPerfStats)))){
        dir.create(file.path(filepathPerfStats))
      }
      if(!( dir.exists( file.path(filepathPerfStatsPerFold)))){
        dir.create(file.path(filepathPerfStatsPerFold))
      }
      if(!( dir.exists( file.path(filepathConfusionMatrix)))){
        dir.create(file.path(filepathConfusionMatrix))
      }
      
      filenameROC <- paste0(filepathPerfPlots,"/",
                            dt,
                            " -- ",
                            st,
                            " -- ROC.png")
      
      filenamePR <- paste0(filepathPerfPlots,"/",
                           dt,
                           " -- ",
                           st,
                           " -- PR.png")
      
      filenameROCData <- paste0(filepathPerfPlotsDataROC,"/",
                                dt,
                                " -- ",
                                st,
                                " -- ROC.csv")
      
      filenamePRData <- paste0(filepathPerfPlotsDataPR,"/",
                               dt,
                               " -- ",
                               st,
                               " -- PR.csv")
      
      filenameCM <- paste0(filepathConfusionMatrix,"/",
                           dt,
                           " -- ",
                           st,
                           " -- CM.txt")
      
      filenameFeatures <- paste0(filepathFeatures,"/",
                                 dt,
                                 " -- ",
                                 st,
                                 " -- Features.csv")
      
      filenameCSV <- paste0(filepathPerfStats,"/",
                            dt,
                            " -- ",
                            st,
                            " -- Perf.csv")

      filenamePerFoldCSV <- paste0(filepathPerfStatsPerFold,"/",
                            dt,
                            " -- ",
                            st,
                            " -- PerfPerFold.csv")
      
      write.csv(perf[[ii]], file = filenameCSV)
      write.csv(rep_perf[[ii]], file = filenamePerFoldCSV)
      
      png(filename=filenameROC, width = 6, height = 4, units = 'in', res = 300)
      plot(prroc_roc[[ii]])
      dev.off()
      
      png(filename=filenamePR, width = 6, height = 4, units = 'in', res = 300)
      plot(prroc_pr[[ii]])
      dev.off()
      
      write.table(prCurveData, sep=",", file = filenamePRData, col.names = FALSE)
      write.table(rocCurveData, sep=",", file = filenameROCData, col.names = FALSE)
      
      # EXTRACT AND SAVE RANKED FEATURE IMPORTANCE INFORMATION
      varImpBestModelDF <- as.data.frame(varImp( mlModel$finalModel, scale = FALSE ))
      varImpBestModelDF2 <- rownames_to_column(varImpBestModelDF, "Taxa")
      varImpBestModelDF2Ordered <- varImpBestModelDF2[order(-varImpBestModelDF2$Overall),]
      colnames(varImpBestModelDF2Ordered)[2] <- "varImp"
      varImpBestModelDF2OrderedNonzero <- varImpBestModelDF2Ordered[varImpBestModelDF2Ordered$varImp != 0,]
      write.csv(varImpBestModelDF2OrderedNonzero, file = filenameFeatures, row.names = FALSE)
      
      sink(file = filenameCM)
      print(dt)
      print(confusionMatrix)
      sink()
      
      rm(mlModel)
      
    }
    # SUMMARIZE RESULTS
    perfTmp[[kk]] <- do.call(rbind, perf)
    rep_perfTmp[[kk]] <- do.call(rbind, rep_perf)
  }
  # SUMMARIZE RESULTS
  perfTmp2[[jj]] <- do.call(rbind, perfTmp)
  rep_perfTmp2[[jj]] <- do.call(rbind, rep_perfTmp)

  write.csv(perfTmp2[[jj]], file = paste0("perfFungi_10k_rep1_tcga_decontamV2_14Oct21_",datasetName,".csv"))
  write.csv(rep_perfTmp2[[jj]], file = paste0("rep_perfFungi_10k_rep1_tcga_decontamV2_14Oct21_",datasetName,".csv"))
}

# SUMMARIZE RESULTS
perf1VsAll <- do.call(rbind, perfTmp2)
rep_perf1VsAll <- do.call(rbind, rep_perfTmp2)

write.csv(perf1VsAll, file = paste0("perfFungi_10k_rep1_tcga_higher_taxa_and_intersected_ALL_decontamV2_14Oct21.csv"))
write.csv(rep_perf1VsAll, file = paste0("rep_perfFungi_10k_rep1_tcga_higher_taxa_and_intersected_ALL_decontamV2_14Oct21.csv"))

#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------