#-----------------------------------------------------------------------------
# 07R-TCGA-compare-tumor-vs-normal-bray-curtis.R
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
# Import data
#----------------------------------------------------------#

load("Interim_data/data_for_ml_tcga_decontamV2_2Apr22.RData", verbose = T)

#----------------------------------------------------------#
# Call Bray Curtis function
#----------------------------------------------------------#

source("00-Functions.R") # for compareTvsNATbrayCurtis() function
bcPlot_SpeciesVSNM <- compareTvsNATbrayCurtis(rep200FungiDecontamV2SpeciesVSNM, 
                                              metaQiitaCombined_Nonzero_DecontamV2, 
                                              axes=c(1,2),
                                              grepFlag = TRUE, 
                                              greplString = "Lung|Breast|Colon|Rectum|Ovarian", 
                                              scalar = 10^4)
bcPlot_SpeciesVSNM$brayPlot2D
ggsave("Figures/Supplementary_Figures/tcga_2D_bray_curtis_species_vsnm_scaled_2Apr22.pdf",
         dpi = "retina", units = "in", width = 8, height = 6)
bcPlot_SpeciesVSNM$brayData %>%
  write.csv("Figures_data/Supplementary_Figures/tcga_2D_bray_curtis_species_vsnm_scaled_2Apr22.csv")

bcPlot_SpeciesVSNM$brayPlot3D

# Save data
write.csv(bcPlot_SpeciesVSNM$countDataRAfilt,
          file = "Interim_data/t_vs_nat_scaled_rel_abundances_TCGA_2Apr22.csv",
          row.names = TRUE)
write.csv(bcPlot_SpeciesVSNM$metaDataFilt,
          file = "Interim_data/t_vs_nat_scaled_metadata_TCGA_2Apr22.csv",
          row.names = TRUE)



