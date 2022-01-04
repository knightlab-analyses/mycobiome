#-----------------------------------------------------------------------------
# 08-TCGA-alpha-diversity-correlation-fungi-vs-bacteria.R
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
load("Interim_data/phyloseq_tcga_rep200_all_OGUs_16Sep21.RData")

#----------------------------------------------------------#
# Load function script and run
#----------------------------------------------------------#

source("00-Functions.R") # for calcCorrAlpha() function
# NOTE: If "rarefyNumber" is not filled, then the function becomes *interactive* 
# and will ask you to provide a rarefaction value after showing the sample read distribution.
# NOTE: This rarefaction value reflects *all* microbial reads mapped to the rep200 database.
# 130,000 reads/sample is approximately the 1st quartile of the sample read distribution
alphaCorr_All_WGS_PT <- calcCorrAlpha(psAll = psRep200All, 
                                      dataType = "WGS", 
                                      sampleType = "Primary Tumor",
                                      rarefyNumber = 130000,
                                      seqCenter = "All Seq Centers",
                                      statCorMethod = "Pearson",
                                      fungi_gOTUs = fungiOGUs, 
                                      bacteria_gOTUs = bacteriaOGUs,
                                      libraryCorrFlag = TRUE)

alphaCorr_All_RNA_PT <- calcCorrAlpha(psAll = psRep200All, 
                                      dataType="RNA-Seq",
                                      sampleType = "Primary Tumor",
                                      rarefyNumber = 150,
                                      seqCenter = "All Seq Centers", 
                                      statCorMethod = "Pearson",
                                      fungi_gOTUs = fungiOGUs, 
                                      bacteria_gOTUs = bacteriaOGUs,
                                      libraryCorrFlag = TRUE)

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
rarefactionTitle <- sprintf("Rarefied %d reads",130000)
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
  stat_cor(method = "pearson", cor.coef.name = "R") -> richnessCorrPlotPerCT_Weizmann
richnessCorrPlotPerCT_Weizmann

filePath <- "Figures/Figure_3/"
baseName <- paste(paste0("rarefied",130000),
                  gsub('([[:punct:]])|\\s+','',"All Seq Centers"), 
                  gsub('([[:punct:]])|\\s+','',"Primary Tumor"),
                  gsub('([[:punct:]])|\\s+','',"WGS"), 
                  "Pearson", sep = "_")
richnessCorrPlotPerCT_Weizmann + ggsave(filename = paste0("weizmann_matched_alpha_div_per_CT_",baseName,".jpeg"), 
                                        path = filePath,
                                        dpi = "retina", units = "in", width = 8)

alphaCorr_All_WGS_PT_Filt_Weizmann %>% 
  write.csv(paste0("Figures_data/Figure_3/weizmann_matched_alpha_div_per_CT_",baseName,".csv"))

#-------------------------------PT with Pearson correlation-------------------------------#
alphaCorr_HMS_WGS_PT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=150000,
                                        seqCenter = "Harvard Medical School", statCorMethod = "Pearson")
alphaCorr_BCM_WGS_PT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=120000,
                                        seqCenter = "Baylor College of Medicine", statCorMethod = "Pearson")
alphaCorr_WashU_WGS_PT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=120000,
                                          seqCenter = "Washington University School of Medicine", statCorMethod = "Pearson")
alphaCorr_Broad_WGS_PT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=60000,
                                          seqCenter = "Broad Institute of MIT and Harvard", statCorMethod = "Pearson")
alphaCorr_MDA_WGS_PT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=130000,
                                        seqCenter = "MD Anderson - Institute for Applied Cancer Science", statCorMethod = "Pearson")
alphaCorr_UNC_RNA_PT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=240,
                                        seqCenter = "University of North Carolina", dataType = "RNA-Seq", statCorMethod = "Pearson")
alphaCorr_CMS_RNA_PT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=40000,
                                        seqCenter = "Canada's Michael Smith Genome Sciences Centre", dataType = "RNA-Seq", statCorMethod = "Pearson")
alphaCorr_Broad_RNA_PT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=290000,
                                          seqCenter = "Broad Institute of MIT and Harvard", dataType = "RNA-Seq", statCorMethod = "Pearson")

#-------------------------------PT with Spearman correlation-------------------------------#
alphaCorr_HMS_WGS_PT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=150000,
                                        seqCenter = "Harvard Medical School", statCorMethod = "Spearman")
alphaCorr_BCM_WGS_PT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=120000,
                                        seqCenter = "Baylor College of Medicine", statCorMethod = "Spearman")
alphaCorr_WashU_WGS_PT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=120000,
                                          seqCenter = "Washington University School of Medicine", statCorMethod = "Spearman")
alphaCorr_Broad_WGS_PT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=60000,
                                          seqCenter = "Broad Institute of MIT and Harvard", statCorMethod = "Spearman")
alphaCorr_MDA_WGS_PT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=130000,
                                        seqCenter = "MD Anderson - Institute for Applied Cancer Science", statCorMethod = "Spearman")
alphaCorr_UNC_RNA_PT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=240,
                                        seqCenter = "University of North Carolina", dataType = "RNA-Seq", statCorMethod = "Spearman")
alphaCorr_CMS_RNA_PT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=40000,
                                        seqCenter = "Canada's Michael Smith Genome Sciences Centre", dataType = "RNA-Seq", statCorMethod = "Spearman")
alphaCorr_Broad_RNA_PT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=290000,
                                          seqCenter = "Broad Institute of MIT and Harvard", dataType = "RNA-Seq", statCorMethod = "Spearman")

#-------------------------------NAT with Pearson correlation-------------------------------#
alphaCorr_HMS_WGS_NAT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=290000,
                                         seqCenter = "Harvard Medical School", sampleType = "Solid Tissue Normal", statCorMethod = "Pearson")
alphaCorr_BCM_WGS_NAT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=90000,
                                         seqCenter = "Baylor College of Medicine", sampleType = "Solid Tissue Normal", statCorMethod = "Pearson")
alphaCorr_WashU_WGS_NAT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=80000,
                                           seqCenter = "Washington University School of Medicine", sampleType = "Solid Tissue Normal", statCorMethod = "Pearson")
alphaCorr_Broad_WGS_NAT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=70000,
                                           seqCenter = "Broad Institute of MIT and Harvard", sampleType = "Solid Tissue Normal", statCorMethod = "Pearson")
alphaCorr_MDA_WGS_NAT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=190000,
                                         seqCenter = "MD Anderson - Institute for Applied Cancer Science", sampleType = "Solid Tissue Normal", statCorMethod = "Pearson")
alphaCorr_UNC_RNA_NAT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=150,
                                         seqCenter = "University of North Carolina", dataType = "RNA-Seq", sampleType = "Solid Tissue Normal", statCorMethod = "Pearson")
alphaCorr_CMS_RNA_NAT_Pearson <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=40000,
                                         seqCenter = "Canada's Michael Smith Genome Sciences Centre", dataType = "RNA-Seq", sampleType = "Solid Tissue Normal", statCorMethod = "Pearson")

#-------------------------------NAT with Spearman correlation-------------------------------#
alphaCorr_HMS_WGS_NAT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=290000,
                                         seqCenter = "Harvard Medical School", sampleType = "Solid Tissue Normal", statCorMethod = "Spearman")
alphaCorr_BCM_WGS_NAT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=90000,
                                         seqCenter = "Baylor College of Medicine", sampleType = "Solid Tissue Normal", statCorMethod = "Spearman")
alphaCorr_WashU_WGS_NAT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=80000,
                                           seqCenter = "Washington University School of Medicine", sampleType = "Solid Tissue Normal", statCorMethod = "Spearman")
alphaCorr_Broad_WGS_NAT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=70000,
                                           seqCenter = "Broad Institute of MIT and Harvard", sampleType = "Solid Tissue Normal", statCorMethod = "Spearman")
alphaCorr_MDA_WGS_NAT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=190000,
                                         seqCenter = "MD Anderson - Institute for Applied Cancer Science", sampleType = "Solid Tissue Normal", statCorMethod = "Spearman")
alphaCorr_UNC_RNA_NAT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=150,
                                         seqCenter = "University of North Carolina", dataType = "RNA-Seq", sampleType = "Solid Tissue Normal", statCorMethod = "Spearman")
alphaCorr_CMS_RNA_NAT_Spearman <- calcCorrAlpha(psAll = psRep200All, rarefyNumber=40000,
                                         seqCenter = "Canada's Michael Smith Genome Sciences Centre", dataType = "RNA-Seq", sampleType = "Solid Tissue Normal", statCorMethod = "Spearman")


