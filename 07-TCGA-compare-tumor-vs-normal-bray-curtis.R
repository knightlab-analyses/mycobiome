#-----------------------------------------------------------------------------
# 07-TCGA-compare-tumor-vs-normal-bray-curtis.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Calculate average relative abundances on TCGA data across groups (tumor vs. normal)
# - Use bray curtis to calculate distances between groups (tumor vs. normal)
# - Plot results in 2D (using ggplot) and 3D (using plotly)
#-----------------------------------------------------------------------------

#----------------------------------------------------------#
# Load environments
#----------------------------------------------------------#

# Load dependencies
require(devtools)
require(doMC)
require(phyloseq)
require(microbiome)
require(vegan)
require(dplyr)
require(reshape2)
require(ggpubr)
require(ggsci)
require(ggforce)
require(concaveman)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Import data: cancers matched to Weizmann
#----------------------------------------------------------#

# save(metaQiitaCombined_Nonzero_8cancer,
#      rep200FungiPhylum_8cancer,
#      rep200FungiClass_8cancer,
#      rep200FungiOrder_8cancer,
#      rep200FungiFamily_8cancer,
#      rep200FungiGenus_8cancer,
#      rep200FungiSpecies_8cancer,
#      file = "Interim_data/raw_data_tcga_8_cancers_16Sep21.RData")

load("Interim_data/raw_data_tcga_8_cancers_16Sep21.RData")

# save(metaQiitaCombined_Nonzero_8cancer,
#      rep200FungiOrder_8cancer_VSNM,
#      rep200FungiFamily_8cancer_VSNM,
#      rep200FungiGenus_8cancer_VSNM,
#      rep200FungiSpecies_8cancer_VSNM,
#      rep200FungiFamily_8cancer_VSNM_CT,
#      rep200FungiGenus_8cancer_VSNM_CT,
#      rep200FungiSpecies_8cancer_VSNM_CT,
#      file = "Interim_data/data_for_ml_8_cancers_13Sep21.RData")

load("Interim_data/data_for_ml_8_cancers_13Sep21.RData")

# save(metaQiitaCombined_Nonzero,
#      rep200FungiPhylum,
#      rep200FungiClass,
#      rep200FungiOrder,
#      rep200FungiFamily,
#      rep200FungiGenus,
#      rep200FungiSpecies,
#      file = "Interim_data/raw_data_tcga_all_cancers_16Sep21.RData")

load("Interim_data/raw_data_tcga_all_cancers_16Sep21.RData")

# save(metaQiitaCombined_Nonzero,
#      rep200FungiOrderVSNM,
#      rep200FungiFamilyVSNM,
#      rep200FungiGenusVSNM,
#      rep200FungiSpeciesVSNM,
#      rep200FungiGenusVSNM_CT,
#      rep200FungiSpeciesVSNM_CT,
#      snmDataOGUFungi,
#      file = "Interim_data/data_for_ml_tcga_13Sep21.RData")

load("Interim_data/data_for_ml_tcga_13Sep21.RData")

# # snmDataOGUFungi,
# # vdge_dataE,
# # rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero,
# # metaQiitaCombined_Nonzero,
load("Interim_data/snmDataFungi_13Sep21.RData") # To load the metadata and raw data objects
# load("Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy_16Sep21.RData")

#----------------------------------------------------------#
# Call Bray Curtis function
#----------------------------------------------------------#

all(rownames(rep200FungiSpecies_8cancer_VSNM) == rownames(metaQiitaCombined_Nonzero_8cancer)) # TRUE

source("00-Functions.R") # for compareTvsNATbrayCurtis() function

bcPlot_SpeciesVSNM <- compareTvsNATbrayCurtis(rep200FungiSpeciesVSNM, metaQiitaCombined_Nonzero, axes=c(1,2),
                                grepFlag = TRUE, 
                                greplString = "Lung|Breast|Colon|Rectum|Ovarian", 
                                scalar = 10^4)
bcPlot_SpeciesVSNM$brayPlot2D +
  ggsave("Figures/Figure_4/figure_4_I__tcga_2D_bray_curtis_species_vsnm_scaled.jpeg",
         dpi = "retina", units = "in", width = 8, height = 6)
bcPlot_SpeciesVSNM$brayData %>%
  write.csv("Figures_data/Figure_4/figure_4_I__tcga_2D_bray_curtis_species_vsnm_scaled.csv")


