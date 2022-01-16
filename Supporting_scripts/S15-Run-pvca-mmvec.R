#-----------------------------------------------------------------------------
# S15-Run-pvca-mmvec.R
# Copyright (c) 2021--, Greg Poore
# Purpose: Run PVCA on fungi and bacterial data for MMvec data
#-----------------------------------------------------------------------------

#-------------------------------#
# Load dependencies
require(tibble) # for data manipulation
require(dplyr) # for data manipulation
require(doMC) # for parallel computing
require(lme4) # for linear modelling calculations
require(ggplot2) # for plotting
require(ggpubr) # for plotting
require(ggsci) # for plotting

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import and run pvca script locally
source("S02-pvca-function.R")

## Import data
load("data_for_pvca_mmvec_fungi_and_bacteria_WGS_19Oct21.RData")

## Write function for running PVCA, saving the results, and plotting the results

runPVCA <- function(metaData, rawData, pvcaThreshold = 0.8, fileBaseName){
    metaFiltered <- metaData[,c("disease_type",
                              "data_submitting_center_label")]

    pvcaObj <- PVCA(counts = t(rawData),
                    meta = metaFiltered,
                    threshold = pvcaThreshold,
                    inter = FALSE)
    # Save data
    fileNameData <- paste0("pvca_",fileBaseName,"_19Oct21.RData")
    save(pvcaObj, file = fileNameData)
    return(pvcaObj)
}

plotPVCA <- function(pvcaObj1,pvcaObj2,pvcaObj3, description1, description2, description3, fileBaseName){
  # Round data
    pvcaObj1Round <- round(pvcaObj1,3)
    pvcaObj2Round <- round(pvcaObj2,3)
    pvcaObj3Round <- round(pvcaObj3,3)

    # Summarize results into data frame
    pvcaRes <- data.frame('Disease Type' = c(pvcaObj1Round[1], pvcaObj2Round[1], pvcaObj3Round[1]),
                          'Sequencing Center' = c(pvcaObj1Round[2], pvcaObj2Round[2], pvcaObj3Round[2]),
                          'Residual\n(not explained by\ntechnical variation)' = c(pvcaObj1Round[3], pvcaObj2Round[3], pvcaObj3Round[3]),
                          data_type = factor(c(description1,description2,description3),
                                             levels = c(description1,description2,description3)),
                          check.names = FALSE)

    # Melt data using reshape2, plot using ggpubr, and save plot
    pvcaRes.melted <- reshape2::melt(pvcaRes, id.vars = "data_type")
    fileNamePlot <- paste0("pvca_plot_",fileBaseName,"_19Oct21.jpeg")
    pvcaRes.melted %>%
      ggbarplot(x = "variable",
                y = "value",
                fill = "data_type",
                palette = "nejm",
                legend = "top",
                ylim = c(0,1),
                xlab = "Biological Effects & Technical Effects",
                ylab = "Weighted average proportion variance",
                label = TRUE,
                position = position_dodge(0.9)) +
      labs(fill = "Data type") +
      ggsave(fileNamePlot,
             dpi = "retina",
             width = 14,
             height = 3,
             units = "in")
}

# Decontam fungi data
decontamV2Fungi <- runPVCA(metaData=metaQiitaCombined_Nonzero_DecontamV2_PT_WGS, rawData=rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_PT_WGS, fileBaseName="mmvec_decontamV2_fungi_PT_WGS")
# WIS shared fungi data
fungiSharedWIS <- runPVCA(metaData=metaImmunePT_WGS, rawData=rep200Data_Matched2ImmunePT_WGS_Fungi_species, fileBaseName="mmvec_WIS_shared_fungi_PT_WGS")
# WIS shared bacterial data
bacteriaSharedWIS <- runPVCA(metaData=metaImmunePT_WGS_Bacteria_WIS_Nonzero, rawData=rep200Data_Matched2ImmunePT_WGS_Bacteria_WIS_Nonzero, fileBaseName="mmvec_WIS_shared_bacteria_PT_WGS")

# Plot PVCA results
plotPVCA(pvcaObj1=decontamV2Fungi,pvcaObj2=fungiSharedWIS,pvcaObj3=bacteriaSharedWIS, 
  description1="TCGA decontaminated fungi", description2="TCGA fungi shared with WIS", description3="TCGA bacteria shared with WIS", fileBaseName="pvca_plot_mmvec")


