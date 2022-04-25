#-----------------------------------------------------------------------------
# 08R-TCGA-alpha-diversity-correlation-fungi-vs-bacteria.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Rarefy the original, full rep200 OTU table (output by Qiita)
# - Separate rarefied data into bacteria and fungi tables
# - Calculate alpha diversity for each table by cancer type and in aggregate
# - Correlate alpha diversities between bacteria and fungi
# - Repeat the above for each sequencing center (WGS only)
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
require(tibble)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Import data: bacteria and fungi data
#----------------------------------------------------------#

# psRep200All, fungiOGUs, bacteriaOGUs,
load("Interim_data/phyloseq_tcga_rep200_all_OGUs_25Mar25.RData")

#----------------------------------------------------------#
# Load function script and run
#----------------------------------------------------------#

source("00-Functions.R") # for calcCorrAlpha() function
# NOTE: If "rarefyNumber" is not filled, then the function becomes *interactive* 
# and will ask you to provide a rarefaction value after showing the sample read distribution.
# NOTE: This rarefaction value reflects *all* microbial reads mapped to the rep200 database.
# 115,000 reads/sample is approximately the 1st quartile of the sample read distribution
# NOTE: This step may take a few minutes to run
alphaCorr_All_WGS_PT <- calcCorrAlpha(psAll = psRep200All, 
                                      dataType = "WGS", 
                                      sampleType = "Primary Tumor",
                                      rarefyNumber = 115000,
                                      seqCenter = "All Seq Centers",
                                      statCorMethod = "Spearman",
                                      fungi_gOTUs = fungiOGUs, 
                                      bacteria_gOTUs = bacteriaOGUs,
                                      libraryCorrFlag = TRUE)
alphaCorr_All_WGS_PT$richnessCorrPlotPerCT

#-------------------------------Subset WGS analysis to Weizmann cancer types-------------------------------#
alphaCorr_All_WGS_PT$fungi_and_bacteria_alpha_with_meta %>%
  filter(investigation %in% c("TCGA-BRCA","TCGA-LUAD","TCGA-LUSC","TCGA-SKCM",
                              "TCGA-COAD","TCGA-READ","TCGA-GBM","TCGA-OV","TCGA-SARC","TCGA-PAAD")) %>%
  droplevels() -> alphaCorr_All_WGS_PT_Filt_Weizmann

alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol <- alphaCorr_All_WGS_PT_Filt_Weizmann$investigation
alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol[alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol == "TCGA-COAD"] <- "Colorectal"
alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol[alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol == "TCGA-READ"] <- "Colorectal"
alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol[alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol == "TCGA-LUAD"] <- "Lung"
alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol[alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol == "TCGA-LUSC"] <- "Lung"
alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol[alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol == "TCGA-SKCM"] <- "Melanoma"
alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol[alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol == "TCGA-GBM"] <- "Glioblastoma"
alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol[alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol == "TCGA-OV"] <- "Ovary"
alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol[alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol == "TCGA-SARC"] <- "Bone"
alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol[alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol == "TCGA-PAAD"] <- "Pancreas"
alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol[alphaCorr_All_WGS_PT_Filt_Weizmann$investigation_consol == "TCGA-BRCA"] <- "Breast"

maxFungiObserved <- max(alphaCorr_All_WGS_PT_Filt_Weizmann$fungi_Observed)
rarefactionTitle <- sprintf("Rarefied %d reads",115000)
plotTitlePerCT <- paste(rarefactionTitle, "All Seq Centers", "Primary Tumor", "WGS", "Pearson", sep = " | ")
alphaCorr_All_WGS_PT_Filt_Weizmann %>%
  ggplot(aes(bacteria_Observed,fungi_Observed)) +
  geom_point() +
  xlab("Bacterial richness") + ylab("Fungal richness") + theme_pubr() +
  facet_wrap(facets = vars(investigation_consol),  nrow = 2) +
  scale_y_continuous(limits = c(0, maxFungiObserved+20)) +
  stat_smooth(method = "lm", formula = y~x) +
  # geom_smooth() +
  ggtitle(plotTitlePerCT) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_igv() +
  stat_cor(method = "spearman", cor.coef.name = "rho") -> richnessCorrPlotPerCT_Weizmann
richnessCorrPlotPerCT_Weizmann

filePath <- "Figures/Main_Figures/"
baseName <- paste(paste0("rarefied",115000),
                  gsub('([[:punct:]])|\\s+','',"All Seq Centers"), 
                  gsub('([[:punct:]])|\\s+','',"Primary Tumor"),
                  gsub('([[:punct:]])|\\s+','',"WGS"), 
                  "Pearson", sep = "_")
ggsave(plot=richnessCorrPlotPerCT_Weizmann,
       filename = paste0("weizmann_matched_alpha_div_per_CT_",baseName,".pdf"), 
       path = filePath, dpi = "retina", units = "in", width = 8)

alphaCorr_All_WGS_PT_Filt_Weizmann %>% 
  write.csv(paste0("Figures_data/Main_Figures/weizmann_matched_alpha_div_per_CT_",baseName,".csv"))
