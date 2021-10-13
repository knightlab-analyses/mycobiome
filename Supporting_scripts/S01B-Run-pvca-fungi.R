#-----------------------------------------------------------------------------
# S01-Run-pvca-fungi.R
# Copyright (c) 2021--, Greg Poore
# Purpose: Run PVCA on fungi data pre and post batch correction
#-----------------------------------------------------------------------------

#-------------------------------#
# Load dependencies
require(tibble)
require(dplyr) # for data manipulation
require(doMC) # for parallel computing
require(lme4)

numCores <- detectCores()
registerDoMC(cores=numCores)

## Import and run pvca script locally
source("S02-pvca-function.R")

## Import data
load("snmDataFungi_DecontamV2_13Oct21.RData")

metaFiltered <- metaQiitaCombined_Nonzero_DecontamV2[,c("sample_type",
						                              "disease_type",
						                              "data_submitting_center_label",
						                              "experimental_strategy")]

pvcaThreshold <- 0.8
pvcaRaw <- PVCA(counts = t(rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero),
                meta = metaFiltered,
                threshold = pvcaThreshold,
                inter = FALSE)

pvcaVoom <- PVCA(counts = t(vdge_dataE_DecontamV2),
                 meta = metaFiltered,
                 threshold = pvcaThreshold,
                 inter = FALSE)

pvcaVSNM <- PVCA(counts = t(snmDataOGUFungiDecontamV2),
                 meta = metaFiltered,
                 threshold = pvcaThreshold,
                 inter = FALSE)

save(pvcaRaw, 
	pvcaVoom,
	pvcaVSNM,
	file = "pvca_decontamV2_fungi_results_raw_Voom_VSNM_13Oct21.RData")

## Plot data
require(ggpubr)
require(ggsci)
require(reshape2)

pvcaRawRound <- round(pvcaRaw,3)
pvcaVoomRound <- round(pvcaVoom,3)
pvcaVSNMRound <- round(pvcaVSNM,3)

pvcaRes <- data.frame('Sample Type' = c(pvcaRawRound[1], pvcaVoomRound[1], pvcaVSNMRound[1]),
                      'Disease Type' = c(pvcaRawRound[2], pvcaVoomRound[2], pvcaVSNMRound[2]),
                      'Sequencing Center' = c(pvcaRawRound[3], pvcaVoomRound[3], pvcaVSNMRound[3]),
                      'Experimental Strategy' = c(pvcaRawRound[4], pvcaVoomRound[4], pvcaVSNMRound[4]),
                      'Residual\n(not explained by\ntechnical variation)' = c(pvcaRawRound[5], pvcaVoomRound[5], pvcaVSNMRound[5]),
                      data_type = factor(c("Raw count data","Voom Normalized Data","Voom Normalized & SNM Corrected Data"),
                                         levels = c("Raw count data","Voom Normalized Data","Voom Normalized & SNM Corrected Data")),
                      check.names = FALSE)

pvcaRes.melted <- reshape2::melt(pvcaRes, id.vars = "data_type")
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
  ggsave("pvca_plot_OGUs_DecontamV2_13Oct21.jpeg",
         dpi = "retina",
         width = 12,
         height = 3,
         units = "in")