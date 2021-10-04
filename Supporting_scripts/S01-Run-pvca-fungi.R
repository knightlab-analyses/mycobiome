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
load("snmDataFungi_13Sep21.RData")

metaFiltered <- metaQiitaCombined_Nonzero[,c("sample_type",
						                              "disease_type",
						                              "data_submitting_center_label",
						                              "experimental_strategy")]

pvcaThreshold <- 0.8
pvcaRaw <- PVCA(counts = t(rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero),
                meta = metaFiltered,
                threshold = pvcaThreshold,
                inter = FALSE)

pvcaVoom <- PVCA(counts = t(vdge_dataE),
                 meta = metaFiltered,
                 threshold = pvcaThreshold,
                 inter = FALSE)

pvcaVSNM <- PVCA(counts = t(snmDataOGUFungi),
                 meta = metaFiltered,
                 threshold = pvcaThreshold,
                 inter = FALSE)

save(pvcaRaw, 
	pvcaVoom,
	pvcaVSNM,
	file = "pvca_fungi_results_raw_Voom_VSNM_13Sep21.RData")