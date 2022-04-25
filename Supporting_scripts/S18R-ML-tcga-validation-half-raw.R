#-----------------------------------------------------------------------------
# S18R-ML-tcga-validation-half-raw.R
# Copyright (c) 2021--, Greg Poore
# Purpose: Compare ML on TCGA validation halves using raw data
#-----------------------------------------------------------------------------
# Run using interactive slurm job using srun --nodes=1 --cpus-per-task=64 --mem=64gb --time=24:00:00 --pty bash -i
# Run using tcgaAnalysisPythonR conda environment

#-------------------------------#
# Load dependencies
require(splitstackshape)
require(reshape2)
require(plyr)
require(dplyr)
require(caret) # for model building
require(doMC) # for parallel computing
require(gbm) # for machine learning
require(tibble) # for df operations
require(PRROC) # for precision-recall curves
require(MLmetrics) # for multi-class learning
require(cvAUC) # for AUC confidence intervals
require(ggpubr) # for plotting
require(ggsci) # for plotting

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import data
load("validation_2halves_tcga_vsnm_data_4Apr22.RData")

crossTestSplitML_Raw <- function(metadataSplit1,
                             metadataSplit2,
                             metadataFull,
                             vsnmDataSplit1,
                             vsnmDataSplit2,
                             vsnmDataFull,
                             sampleTypeComparison = "Primary Tumor",
                             caretModel = "gbm",
                             samplingStrategy = "up",
                             numResampleIter = 1,
                             numKFold = 4,
                             baseName = "val_halves_k4",
                             caretTuneGrid = data.frame(n.trees=150, interaction.depth=3,
                                                        shrinkage=0.1,
                                                        n.minobsinnode=1),
                             sampleTypeList = c("Primary Tumor")){
  require(caret) # for model building
  require(gbm) # for machine learning
  require(PRROC) # for precision-recall curves
  require(MLmetrics) # for multiclass ML
  require(e1071)
  require(ggpubr)
  require(ggrepel)
  
  perf <- list()
  perfTmp <- list()
  perfTmp2 <- list()
  
  for(kk in seq_along(sampleTypeList)){
    st <- sampleTypeList[kk]
    print(st)
    
    if(st == "Primary Tumor vs Solid Tissue Normal"){
      metadataSplit1_FiltSt <- droplevels(metadataSplit1[metadataSplit1$sample_type %in%
                                                           c("Primary Tumor", "Solid Tissue Normal"),])
      metadataSplit2_FiltSt <- droplevels(metadataSplit2[metadataSplit2$sample_type %in%
                                                           c("Primary Tumor", "Solid Tissue Normal"),])
      metadataFull_FiltSt <- droplevels(metadataFull[metadataFull$sample_type %in%
                                                       c("Primary Tumor", "Solid Tissue Normal"),])
    } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
      metadataSplit1_FiltSt <- droplevels(metadataSplit1[metadataSplit1$sample_type %in% st,])
      metadataSplit2_FiltSt <- droplevels(metadataSplit2[metadataSplit2$sample_type %in% st,])
      metadataFull_FiltSt <- droplevels(metadataFull[metadataFull$sample_type %in% st,])
    }
    
    # print(seq_along(levels(metaTmp2$disease_type)))
    
    for(ii in seq_along(levels(metadataSplit1_FiltSt$disease_type))){
      
      dt <- levels(metadataSplit1_FiltSt$disease_type)[ii]
      print(dt)
      
      if(st == "Primary Tumor vs Solid Tissue Normal"){
        metadataSplit1_FiltStDz <- droplevels(metadataSplit1_FiltSt[metadataSplit1_FiltSt$disease_type == dt,])
        metadataSplit2_FiltStDz <- droplevels(metadataSplit2_FiltSt[metadataSplit2_FiltSt$disease_type == dt,])
        metadataFull_FiltStDz <- droplevels(metadataFull_FiltSt[metadataFull_FiltSt$disease_type == dt,])
        
        metadataSplit1_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataSplit1_FiltStDz$sample_type))
        metadataSplit2_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataSplit2_FiltStDz$sample_type))
        metadataFull_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataFull_FiltStDz$sample_type))
        
        positiveClass <- "PrimaryTumor"
        negativeClass <- "SolidTissueNormal"
      } else if(st %in% c("Blood Derived Normal", "Primary Tumor")){
        metadataSplit1_FiltStDz <- metadataSplit1_FiltSt
        metadataSplit2_FiltStDz <- metadataSplit2_FiltSt
        metadataFull_FiltStDz <- metadataFull_FiltSt
        
        metadataSplit1_FiltStDz$predY <- ifelse(metadataSplit1_FiltStDz$disease_type == dt, yes = dt, no = "OtherCancerType")
        metadataSplit2_FiltStDz$predY <- ifelse(metadataSplit2_FiltStDz$disease_type == dt, yes = dt, no = "OtherCancerType")
        metadataFull_FiltStDz$predY <- ifelse(metadataFull_FiltStDz$disease_type == dt, yes = dt, no = "OtherCancerType")
        
        metadataSplit1_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataSplit1_FiltStDz$predY))
        metadataSplit2_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataSplit2_FiltStDz$predY))
        metadataFull_FiltStDz$predY <- factor(gsub('([[:punct:]])|\\s+','',metadataFull_FiltStDz$predY))
        
        positiveClass <- gsub('([[:punct:]])|\\s+','',dt)
        negativeClass <- "OtherCancerType"
      }
      
      # print(table(metadataSplit1_FiltStDz$predY))
      # print(table(metadataSplit2_FiltStDz$predY))
      # print(table(metadataFull_FiltStDz$predY))
      
      # SKIP CANCER TYPES THAT ONLY HAVE ONE CLASS OF A COMPARISON
      if(length(table(metadataSplit1_FiltStDz$predY)) < 2 | length(table(metadataSplit2_FiltStDz$predY)) < 2){next}
      
      # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF 20 SAMPLES IN EITHER CLASS
      if(any(table(metadataSplit1_FiltStDz$predY) < 20) | any(table(metadataSplit2_FiltStDz$predY) < 20)){next}
      
      mlDataY_Split1 <- metadataSplit1_FiltStDz
      mlDataY_Split2 <- metadataSplit2_FiltStDz
      mlDataY_Full <- metadataFull_FiltStDz
      
      mlDataX_Split1 <- vsnmDataSplit1[rownames(mlDataY_Split1),]
      mlDataX_Split2 <- vsnmDataSplit2[rownames(mlDataY_Split2),]
      mlDataX_Full <- vsnmDataFull[rownames(mlDataY_Full),]
      
      set.seed(42)
      ctrl <- trainControl(method = "repeatedcv",
                           number = numKFold,
                           repeats = numResampleIter,
                           sampling = "up",
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE,
                           verboseIter = FALSE,
                           savePredictions = TRUE,
                           allowParallel=TRUE)
      
      #--------------------------Train split 1--------------------------#
      set.seed(42)
      trainX_Split1 <- mlDataX_Split1
      trainY_Split1 <- mlDataY_Split1[,"predY"]
      testX_Split1 <- mlDataX_Split2
      testY_Split1 <- mlDataY_Split2[,"predY"]
      # refactoredTrainY_Split1 <- factor(gsub('([[:punct:]])|\\s+','',trainY_Split1))
      
      print("Now training model on split 1...")
      set.seed(42)
      mlModel_Split1 <- train(x = trainX_Split1,
                              y = trainY_Split1,
                              method = "gbm",
                              # preProcess = c("scale","center"),
                              trControl = ctrl,
                              verbose = FALSE,
                              metric = "ROC",
                              tuneGrid = caretTuneGrid)
      
      #--------------------------Train split 2--------------------------#
      set.seed(42)
      trainX_Split2 <- mlDataX_Split2
      trainY_Split2 <- mlDataY_Split2[,"predY"]
      testX_Split2 <- mlDataX_Split1
      testY_Split2 <- mlDataY_Split1[,"predY"]
      # refactoredTrainY_Split2 <- factor(gsub('([[:punct:]])|\\s+','',trainY_Split2))
      
      print("Now training model on split 2...")
      set.seed(42)
      mlModel_Split2 <- train(x = trainX_Split2,
                              y = trainY_Split2,
                              method = "gbm",
                              # preProcess = c("scale","center"),
                              trControl = ctrl,
                              verbose = FALSE,
                              metric = "ROC",
                              tuneGrid = caretTuneGrid)
      
      #--------------------------Train half of full data--------------------------#
      
      set.seed(42) # random number seed -- important if you want comparable results between data types
      index <- createDataPartition(mlDataY_Full[,"predY"], p = 0.5, list = FALSE)
      trainHalfFullX <- mlDataX_Full[index,]
      trainHalfFullY <- mlDataY_Full[index,"predY"]
      testHalfFullX <- mlDataX_Full[-index,]
      testHalfFullY <- mlDataY_Full[-index,"predY"]
      
      print("Now training model on half of full data...")
      set.seed(42)
      mlModel_HalfFull <- train(x = trainHalfFullX,
                                y = trainHalfFullY,
                                method = "gbm",
                                # preProcess = c("scale","center"),
                                trControl = ctrl,
                                verbose = FALSE,
                                metric = "ROC",
                                tuneGrid = caretTuneGrid)
      
      #--------------------------Test on opposite split--------------------------#
      print("Now testing models on each other...")
      # Use model from split 1 to predict on split 2
      predClassS1S2 <- predict(mlModel_Split1, newdata = testX_Split1)
      predProbsS1S2 <- as.numeric(predict(mlModel_Split1, newdata = testX_Split1, type = "prob")[,positiveClass])
      # Use model from split 2 to predict on split 1
      predClassS2S1 <- predict(mlModel_Split2, newdata = testX_Split2)
      predProbsS2S1 <- as.numeric(predict(mlModel_Split2, newdata = testX_Split2, type = "prob")[,positiveClass])
      # Use model from half of full data to predict on holdout half
      predClassHalfFull <- predict(mlModel_HalfFull, newdata = testHalfFullX)
      predProbsHalfFull <- as.numeric(predict(mlModel_HalfFull, newdata = testHalfFullX, type = "prob")[,positiveClass])
      
      # Calculate perf of testing model from split 1 to predict on split 2
      fgS1S2 <- predProbsS1S2[testY_Split1 == positiveClass]
      bgS1S2 <- predProbsS1S2[testY_Split1 == negativeClass]
      prroc_rocS1S2 <- roc.curve(scores.class0 = fgS1S2, scores.class1 = bgS1S2, curve = T)
      prroc_prS1S2 <- pr.curve(scores.class0 = fgS1S2, scores.class1 = bgS1S2, curve = T, rand.compute=T)
      aurocS1S2 = prroc_rocS1S2$auc
      auprS1S2 = prroc_prS1S2$auc.integral
      # Calculate perf of testing model from split 2 to predict on split 1
      fgS2S1 <- predProbsS2S1[testY_Split2 == positiveClass]
      bgS2S1 <- predProbsS2S1[testY_Split2 == negativeClass]
      prroc_rocS2S1 <- roc.curve(scores.class0 = fgS2S1, scores.class1 = bgS2S1, curve = T)
      prroc_prS2S1 <- pr.curve(scores.class0 = fgS2S1, scores.class1 = bgS2S1, curve = T, rand.compute=T)
      aurocS2S1 = prroc_rocS2S1$auc
      auprS2S1 = prroc_prS2S1$auc.integral
      # Calculate perf of testing model from half of full data to predict on holdout half
      fgHalfFull <- predProbsHalfFull[testHalfFullY == positiveClass]
      bgHalfFull <- predProbsHalfFull[testHalfFullY == negativeClass]
      prroc_rocHalfFull <- roc.curve(scores.class0 = fgHalfFull, scores.class1 = bgHalfFull, curve = T)
      prroc_prHalfFull <- pr.curve(scores.class0 = fgHalfFull, scores.class1 = bgHalfFull, curve = T, rand.compute=T)
      aurocHalfFull = prroc_rocHalfFull$auc
      auprHalfFull = prroc_prHalfFull$auc.integral
      
      print("Perf for the model built on data from split 1 and tested on data from split 2...")
      # print(confusionMatrixSplit1)
      print(sprintf("AUROC & AUPR: %f.3 | %f.3", aurocS1S2, auprS1S2))
      
      print("Perf for the model built on data from split 2 and tested on data from split 1...")
      # print(confusionMatrixSplit2)
      print(sprintf("AUROC & AUPR: %f.3 | %f.3", aurocS2S1, auprS2S1))
      
      print("Perf for the model built on data from half of full data on holdout half...")
      # print(confusionMatrixSplit2)
      print(sprintf("AUROC & AUPR: %f.3 | %f.3", aurocHalfFull, auprHalfFull))
      cat('\n')
      
      perfTmp[[ii]] <- data.frame(disease_type = dt,
                                  sample_type = st,
                                  auroc_split12 = aurocS1S2,
                                  aupr_split12 = auprS1S2,
                                  auroc_split21 = aurocS2S1,
                                  aupr_split21 = auprS2S1,
                                  auroc_halffull = aurocHalfFull,
                                  aupr_halffull = auprHalfFull)
      rm(mlModel_Split1)
      rm(mlModel_Split2)
      rm(mlModel_HalfFull)
    }
    perfTmp2[[kk]] <- do.call(rbind, perfTmp)
  }
  perf <- do.call(rbind, perfTmp2)
  perf$abbrev <- data.frame(abbrev=gsub("^TCGA-","",unique(metadataFull$investigation)), row.names=unique(metadataFull$disease_type))[perf$disease_type,"abbrev"]
  perf$sample_type <- factor(perf$sample_type, levels = c("Primary Tumor","Blood Derived Normal","Primary Tumor vs Solid Tissue Normal"))
  # Plot comparisons between split 1 and 2
  plot_auroc_split12_vs_split21 <- perf %>% distinct() %>%
    ggplot(aes(x = auroc_split21, y = auroc_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUROC split 2 on split 1", y = "AUROC split 1 on split 2", title = "Comparing AUROC on splits") 
    ggsave(plot=plot_auroc_split12_vs_split21,
           filename = paste0("auroc_",baseName,"_split12_vs_split21.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  plot_aupr_split12_vs_split21 <- perf %>% distinct() %>%
    ggplot(aes(x = aupr_split21, y = aupr_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUPR split 2 on split 1", y = "AUPR split 1 on split 2", title = "Comparing AUPR on splits") 
    ggsave(plot=plot_aupr_split12_vs_split21,
           filename = paste0("aupr_",baseName,"_split12_vs_split21.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  # Plot comparisons between split 1 and half of full
  plot_auroc_split12_vs_halffull <- perf %>% distinct() %>%
    ggplot(aes(x = auroc_halffull, y = auroc_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUROC half of full", y = "AUROC split 1", title = "Comparing AUROC on split 1 vs half of full") 
    ggsave(plot=plot_auroc_split12_vs_halffull,
           filename = paste0("auroc_",baseName,"_split12_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  plot_aupr_split12_vs_halffull <- perf %>% distinct() %>%
    ggplot(aes(x = aupr_halffull, y = aupr_split12, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUPR half of full", y = "AUPR split 1", title = "Comparing AUPR on split 1 vs half of full") 
    ggsave(plot=plot_aupr_split12_vs_halffull,
           filename = paste0("aupr_",baseName,"_split12_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  # Plot comparisons between split 2 and half of full
  plot_auroc_split21_vs_halffull <- perf %>% distinct() %>%
    ggplot(aes(x = auroc_halffull, y = auroc_split21, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUROC half of full", y = "AUROC split 2", title = "Comparing AUROC on split 2 vs half of full") 
    ggsave(plot=plot_auroc_split21_vs_halffull,
           filename = paste0("auroc_",baseName,"_split21_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  plot_aupr_split21_vs_halffull <- perf %>% distinct() %>%
    ggplot(aes(x = aupr_halffull, y = aupr_split21, label = abbrev)) +
    geom_point(alpha = 0.4) + geom_smooth(method='lm', fullrange = TRUE) + 
    stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
    theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
    geom_text_repel() + geom_abline(linetype = 2) + xlim(c(0,1)) + ylim(c(0,1)) +
    facet_wrap(~sample_type) + 
    labs(x = "AUPR half of full", y = "AUPR split 2", title = "Comparing AUPR on split 2 vs half of full") 
    ggsave(plot=plot_aupr_split21_vs_halffull,
           filename = paste0("aupr_",baseName,"_split21_vs_halffull.jpeg"), dpi = "retina", units = "in",
           width = 8, height = 5)
  save(perf, file = paste0("perf_",baseName,".RData"))
  res <- list(
    perf=perf,
    plot_auroc_split12_vs_split21=plot_auroc_split12_vs_split21,
    plot_aupr_split12_vs_split21=plot_aupr_split12_vs_split21,
    plot_auroc_split12_vs_halffull=plot_auroc_split12_vs_halffull,
    plot_aupr_split12_vs_halffull=plot_aupr_split12_vs_halffull,
    plot_auroc_split21_vs_halffull=plot_auroc_split21_vs_halffull,
    plot_aupr_split21_vs_halffull=plot_aupr_split21_vs_halffull)
  return(res)
}

valResSplitsVsHalf_rawData <- crossTestSplitML_Raw(metadataSplit1 = valSplit1Metadata,
                                          metadataSplit2 = valSplit2Metadata, 
                                          metadataFull = metaQiitaCombined_Nonzero_DecontamV2,
                                          vsnmDataSplit1 = valSplit1FungiCountDecontamV2Data, 
                                          vsnmDataSplit2 = valSplit2FungiCountDecontamV2Data, 
                                          vsnmDataFull = rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero,
                                          sampleTypeList = c("Primary Tumor vs Solid Tissue Normal",
                                                             "Blood Derived Normal",
                                                             "Primary Tumor"),
                                          numKFold = 10,
                                          baseName = "raw_data_k10")


