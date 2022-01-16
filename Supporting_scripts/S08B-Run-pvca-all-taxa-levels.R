#-----------------------------------------------------------------------------
# S08B-Run-pvca-all-taxa-levels.R
# Copyright (c) 2021--, Greg Poore
# Purpose: Run PVCA on fungi data pre and post batch correction at various taxa levels
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
load("data_for_pvca_tcga_taxa_levels_decontamV2_14Oct21.RData")

## Write function for running PVCA, saving the results, and plotting the results

runPVCA <- function(metaData=metaQiitaCombined_Nonzero_DecontamV2, rawData, voomData, snmData, pvcaThreshold = 0.8, fileBaseName){
    metaFiltered <- metaData[,c("sample_type",
                              "disease_type",
                              "data_submitting_center_label",
                              "experimental_strategy")]

    pvcaRaw <- PVCA(counts = t(rawData),
                    meta = metaFiltered,
                    threshold = pvcaThreshold,
                    inter = FALSE)

    pvcaVoom <- PVCA(counts = t(voomData),
                     meta = metaFiltered,
                     threshold = pvcaThreshold,
                     inter = FALSE)

    pvcaVSNM <- PVCA(counts = t(snmData),
                     meta = metaFiltered,
                     threshold = pvcaThreshold,
                     inter = FALSE)
    # Save data
    fileNameData <- paste0("pvca_",fileBaseName,"_13Sep21.RData")
    save(pvcaRaw, pvcaVoom, pvcaVSNM, file = fileNameData)

    # Round data
    pvcaRawRound <- round(pvcaRaw,3)
    pvcaVoomRound <- round(pvcaVoom,3)
    pvcaVSNMRound <- round(pvcaVSNM,3)

    # Summarize results into data frame
    pvcaRes <- data.frame('Sample Type' = c(pvcaRawRound[1], pvcaVoomRound[1], pvcaVSNMRound[1]),
                          'Disease Type' = c(pvcaRawRound[2], pvcaVoomRound[2], pvcaVSNMRound[2]),
                          'Sequencing Center' = c(pvcaRawRound[3], pvcaVoomRound[3], pvcaVSNMRound[3]),
                          'Experimental Strategy' = c(pvcaRawRound[4], pvcaVoomRound[4], pvcaVSNMRound[4]),
                          'Residual\n(not explained by\ntechnical variation)' = c(pvcaRawRound[5], pvcaVoomRound[5], pvcaVSNMRound[5]),
                          data_type = factor(c("Raw count data","Voom Normalized Data","Voom Normalized & SNM Corrected Data"),
                                             levels = c("Raw count data","Voom Normalized Data","Voom Normalized & SNM Corrected Data")),
                          check.names = FALSE)

    # Melt data using reshape2, plot using ggpubr, and save plot
    pvcaRes.melted <- reshape2::melt(pvcaRes, id.vars = "data_type")
    fileNamePlot <- paste0("pvca_plot_",fileBaseName,"_13Sep21.jpeg")
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

# Order level
runPVCA(rawData=rep200FungiDecontamV2OrderVSNM_Obj$qcData, voomData=rep200FungiDecontamV2OrderVSNM_Obj$vdge_dataE, snmData=rep200FungiDecontamV2OrderVSNM_Obj$snmData, fileBaseName = "tcga_vsnm_order")
# Family level
runPVCA(rawData=rep200FungiDecontamV2FamilyVSNM_Obj$qcData, voomData=rep200FungiDecontamV2FamilyVSNM_Obj$vdge_dataE, snmData=rep200FungiDecontamV2FamilyVSNM_Obj$snmData, fileBaseName = "tcga_vsnm_family")
# Genus level
runPVCA(rawData=rep200FungiDecontamV2GenusVSNM_Obj$qcData, voomData=rep200FungiDecontamV2GenusVSNM_Obj$vdge_dataE, snmData=rep200FungiDecontamV2GenusVSNM_Obj$snmData, fileBaseName = "tcga_vsnm_genus")
# Species level
runPVCA(rawData=rep200FungiDecontamV2SpeciesVSNM_Obj$qcData, voomData=rep200FungiDecontamV2SpeciesVSNM_Obj$vdge_dataE, snmData=rep200FungiDecontamV2SpeciesVSNM_Obj$snmData, fileBaseName = "tcga_vsnm_species")
# Genus level CT
runPVCA(rawData=rep200FungiDecontamV2GenusVSNM_Obj_CT$qcData, voomData=rep200FungiDecontamV2GenusVSNM_Obj_CT$vdge_dataE, snmData=rep200FungiDecontamV2GenusVSNM_Obj_CT$snmData, fileBaseName = "tcga_vsnm_genus_CT")
# Species level CT
runPVCA(rawData=rep200FungiDecontamV2SpeciesVSNM_Obj_CT$qcData, voomData=rep200FungiDecontamV2SpeciesVSNM_Obj_CT$vdge_dataE, snmData=rep200FungiDecontamV2SpeciesVSNM_Obj_CT$snmData, fileBaseName = "tcga_vsnm_species_CT")

## 8 cancers matched to Weizmann
# Order level
runPVCA(rawData=rep200FungiDecontamV2Order_8cancer_VSNM_Obj$qcData, metaData=metaQiitaCombined_Nonzero_DecontamV2_8cancer,
  voomData=rep200FungiDecontamV2Order_8cancer_VSNM_Obj$vdge_dataE, snmData=rep200FungiDecontamV2Order_8cancer_VSNM_Obj$snmData, fileBaseName = "tcga_vsnm_order_8cancers")
# Family level
runPVCA(rawData=rep200FungiDecontamV2Family_8cancer_VSNM_Obj$qcData, metaData=metaQiitaCombined_Nonzero_DecontamV2_8cancer,
  voomData=rep200FungiDecontamV2Family_8cancer_VSNM_Obj$vdge_dataE, snmData=rep200FungiDecontamV2Family_8cancer_VSNM_Obj$snmData, fileBaseName = "tcga_vsnm_family_8cancers")
# Genus level
runPVCA(rawData=rep200FungiDecontamV2Genus_8cancer_VSNM_Obj$qcData, metaData=metaQiitaCombined_Nonzero_DecontamV2_8cancer,
  voomData=rep200FungiDecontamV2Genus_8cancer_VSNM_Obj$vdge_dataE, snmData=rep200FungiDecontamV2Genus_8cancer_VSNM_Obj$snmData, fileBaseName = "tcga_vsnm_genus_8cancers")
# Species level
runPVCA(rawData=rep200FungiDecontamV2Species_8cancer_VSNM_Obj$qcData, metaData=metaQiitaCombined_Nonzero_DecontamV2_8cancer,
  voomData=rep200FungiDecontamV2Species_8cancer_VSNM_Obj$vdge_dataE, snmData=rep200FungiDecontamV2Species_8cancer_VSNM_Obj$snmData, fileBaseName = "tcga_vsnm_species_8cancers")
# Genus level CT
runPVCA(rawData=rep200FungiDecontamV2Genus_8cancer_VSNM_Obj_CT$qcData, metaData=metaQiitaCombined_Nonzero_DecontamV2_8cancer,
  voomData=rep200FungiDecontamV2Genus_8cancer_VSNM_Obj_CT$vdge_dataE, snmData=rep200FungiDecontamV2Genus_8cancer_VSNM_Obj_CT$snmData, fileBaseName = "tcga_vsnm_genus_8cancers_CT")
# Species level CT
runPVCA(rawData=rep200FungiDecontamV2Species_8cancer_VSNM_Obj_CT$qcData, metaData=metaQiitaCombined_Nonzero_DecontamV2_8cancer,
  voomData=rep200FungiDecontamV2Species_8cancer_VSNM_Obj_CT$vdge_dataE, snmData=rep200FungiDecontamV2Species_8cancer_VSNM_Obj_CT$snmData, fileBaseName = "tcga_vsnm_species_8cancers_CT")


