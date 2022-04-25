#-----------------------------------------------------------------------------
# 15R-ML-WIS-samples.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Perform machine learning to discriminate between and within cancer types using WIS data
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
require(plyr)
require(dplyr)
require(reshape2)
require(ggpubr)
require(ggsci)
require(ANCOMBC)
require(biomformat)
require(Rhdf5lib)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Load WIS data
#----------------------------------------------------------#

#--------------------------WIS fungal data--------------------------#
# NOTE: "fungal_norm_unique_hit.rds" contains a list with 7 phyloseq objects
# containing fungal data summarized at all_rank, phylum, class, order, family,
# genus, and species levels
wzFungiRDS <- readRDS("Input_data/Weizmann_data/fungal_norm_unique_hit.rds")
wzFungiBacteriaRDS <- readRDS("Input_data/Weizmann_data/bacterial_fungal_all_rank_count_hit_list.rds")
psFungiWz_allRank <- wzFungiRDS$all_rank
psFungiBacteriaWz_allRank <- wzFungiBacteriaRDS$all_rank

# Subset to non-control samples to compare overlap
psFungiWz_allRank_Bio <- subset_samples(psFungiWz_allRank, type.detail %in% c("normal", "nat", "tumor"))
psFungiBacteriaWz_allRank_Bio <- subset_samples(psFungiBacteriaWz_allRank, type.detail %in% c("normal", "nat", "tumor"))
# Fungi
wzFungiFreeCounts <- data.frame(t(otu_table(psFungiWz_allRank_Bio)))
wzFungiFreeTax <- data.frame(tax_table(psFungiWz_allRank_Bio))
metaWzFungiFree <- data.frame(sample_data(psFungiWz_allRank_Bio))
# Fungi & Bacteria
wzFungiBacteriaFreeCounts <- data.frame(t(otu_table(psFungiBacteriaWz_allRank_Bio)))
wzFungiBacteriaFreeTax <- data.frame(tax_table(psFungiBacteriaWz_allRank_Bio))
metaWzFungiBacteriaFree <- data.frame(sample_data(psFungiBacteriaWz_allRank_Bio))
# Bacteria
wzBacteriaFreeCounts <- wzFungiBacteriaFreeCounts[,rownames(wzFungiBacteriaFreeTax)[wzFungiBacteriaFreeTax$kingdom == "Bacteria"]]

## Intersecting samples (note that metaWzFungiBacteriaFree is limiting w/ 840 samples)
intersectingSamples <- intersect(rownames(metaWzFungiFree), rownames(metaWzFungiBacteriaFree))
metaWzFungiBacteriaCombinedFreeFilt <- droplevels(metaWzFungiBacteriaFree[intersectingSamples,])

wzFungiFreeCountsFilt <- wzFungiFreeCounts[rownames(metaWzFungiBacteriaCombinedFreeFilt),]
wzFungiBacteriaFreeCountsFilt <- wzFungiBacteriaFreeCounts[rownames(metaWzFungiBacteriaCombinedFreeFilt),]
wzBacteriaFreeCountsFilt <- wzBacteriaFreeCounts[rownames(metaWzFungiBacteriaCombinedFreeFilt),]

#--------------------------Calculate RA tables and add diversity column--------------------------#

ra_wzFungiFreeCountsFilt <- wzFungiFreeCountsFilt/rowSums(wzFungiFreeCountsFilt)
ra_wzFungiBacteriaFreeCountsFilt <- wzFungiBacteriaFreeCountsFilt/rowSums(wzFungiBacteriaFreeCountsFilt)
ra_wzBacteriaFreeCountsFilt <- wzBacteriaFreeCountsFilt/rowSums(wzBacteriaFreeCountsFilt)

ra_wzFungiFreeCountsFilt$diversity <- rowSums(ra_wzFungiFreeCountsFilt>0)
ra_wzFungiBacteriaFreeCountsFilt$diversity <- rowSums(ra_wzFungiBacteriaFreeCountsFilt>0)
ra_wzBacteriaFreeCountsFilt$diversity <- rowSums(ra_wzBacteriaFreeCountsFilt>0)

#--------------------------Calculate presence/absence tables and add diversity column--------------------------#

## Binary matrices based on presence/absence
binary_wzFungiFreeCountsFilt <- data.frame((wzFungiFreeCountsFilt>0)*1)
binary_wzFungiBacteriaFreeCountsFilt <- data.frame((wzFungiBacteriaFreeCountsFilt>0)*1)
binary_wzBacteriaFreeCountsFilt <- data.frame((wzBacteriaFreeCountsFilt>0)*1)

binary_wzFungiFreeCountsFilt$diversity <- rowSums(binary_wzFungiFreeCountsFilt)
binary_wzFungiBacteriaFreeCountsFilt$diversity <- rowSums(binary_wzFungiBacteriaFreeCountsFilt)
binary_wzBacteriaFreeCountsFilt$diversity <- rowSums(binary_wzBacteriaFreeCountsFilt)

#--------------------------Create join count and presence/absence tables--------------------------#

count_binary_wzFungiFreeCountsFilt <- cbind(wzFungiFreeCountsFilt, binary_wzFungiFreeCountsFilt)
count_binary_wzFungiBacteriaFreeCountsFilt <- cbind(wzFungiBacteriaFreeCountsFilt, binary_wzFungiBacteriaFreeCountsFilt)
count_binary_wzBacteriaFreeCountsFilt <- cbind(wzBacteriaFreeCountsFilt, binary_wzBacteriaFreeCountsFilt)

#--------------------------Lastly add diversity column to main data (cannot do this before above)--------------------------#
wzFungiFreeCountsFilt$diversity <- rowSums(wzFungiFreeCountsFilt>0)
wzFungiBacteriaFreeCountsFilt$diversity <- rowSums(wzFungiBacteriaFreeCountsFilt>0)
wzBacteriaFreeCountsFilt$diversity <- rowSums(wzBacteriaFreeCountsFilt>0)

#--------------------------Save data--------------------------#

save(metaWzFungiBacteriaCombinedFreeFilt,
     wzFungiFreeCountsFilt,
     ra_wzFungiFreeCountsFilt,
     binary_wzFungiFreeCountsFilt,
     count_binary_wzFungiFreeCountsFilt,
     
     wzBacteriaFreeCountsFilt,
     ra_wzBacteriaFreeCountsFilt,
     binary_wzBacteriaFreeCountsFilt,
     count_binary_wzFungiBacteriaFreeCountsFilt,
     
     wzFungiBacteriaFreeCountsFilt,
     ra_wzFungiBacteriaFreeCountsFilt,
     binary_wzFungiBacteriaFreeCountsFilt,
     count_binary_wzBacteriaFreeCountsFilt,
     file = "Interim_data/data_for_ml_WIS_samples_fungi_bacteria_5Apr22.RData")

# Scripts: S26R

#----------------------------------------------------------#
# Save genus level fungi and bacterial data for MMvec
#----------------------------------------------------------#
# BRIEF NOTE: There are taxa differences between the files, so use the bacteria from
# bacterial_fungal_all_rank_count_hit_list.rds and fungi from fungal_norm_unique_hit.rds

#----------------------Get bacterial data----------------------#
wzFungiBacteriaRDS <- readRDS("Input_data/Weizmann_data/bacterial_fungal_all_rank_count_hit_list.rds")
psFungiBacteriaWz_genus <- wzFungiBacteriaRDS$genus_phy
wzFungiBacteriaGenusAll <- data.frame(t(otu_table(psFungiBacteriaWz_genus)))
wzMetaFungiBacteriaGenusAll <- data.frame(sample_data(psFungiBacteriaWz_genus))
wzFungiBacteriaTaxaTableAll <- data.frame(tax_table(psFungiBacteriaWz_genus))
wzBacteriaOnlyTaxaTableAll <- wzFungiBacteriaTaxaTableAll[wzFungiBacteriaTaxaTableAll$kingdom == "Bacteria",]
wzFungiOnlyTaxaTableAll <- wzFungiBacteriaTaxaTableAll[wzFungiBacteriaTaxaTableAll$kingdom == "Fungi",]
# Subset to fungi and bacteria only
wzFungiOnlyGenus <- wzFungiBacteriaGenusAll[,rownames(wzFungiBacteriaTaxaTableAll[wzFungiBacteriaTaxaTableAll$kingdom == "Fungi",])]
wzBacteriaOnlyGenus <- wzFungiBacteriaGenusAll[,rownames(wzFungiBacteriaTaxaTableAll[wzFungiBacteriaTaxaTableAll$kingdom == "Bacteria",])]
# Save files
wzFungiOnlyGenus %>% write.csv(file = "Interim_data/genus_data_WIS_full_fungi_count_data_genus_all_summarized_WIS.csv")
wzBacteriaOnlyGenus %>% write.csv(file = "Interim_data/genus_data_WIS_full_bacteria_count_data_genus_all_summarized_WIS.csv")
wzMetaFungiBacteriaGenusAll %>% write.csv(file = "Interim_data/genus_data_WIS_full_metadata_for_fungi_and_bacteria_genus_all_summarized_WIS.csv")
wzFungiBacteriaTaxaTableAll %>% write.csv(file = "Interim_data/genus_data_WIS_full_taxa_table_for_fungi_and_bacteria_genus_all_summarized_WIS.csv")

#----------------------Get fungal data----------------------#
wzRDS <- readRDS("Input_data/Weizmann_data/fungal_norm_unique_hit.rds")
psWz_genus <- wzRDS$genus_phy
wzFungiGenusAll <- data.frame(t(otu_table(psWz_genus)))
wzMetaFungiGenusAll <- data.frame(sample_data(psWz_genus))
wzFungiTaxaTableAll <- data.frame(tax_table(psWz_genus))

intersectedSamples <- intersect(rownames(wzMetaFungiBacteriaGenusAll), rownames(wzMetaFungiGenusAll))
wzMetaGenusIntersected <- droplevels(wzMetaFungiBacteriaGenusAll[intersectedSamples,])
wzBacteriaOnlyGenusIntersected <- wzBacteriaOnlyGenus[intersectedSamples,]
wzFungiGenusIntersected <- wzFungiGenusAll[intersectedSamples,]

wzFungiGenusIntersected %>% write.csv(file = "Interim_data/genus_data_WIS_full_revised_fungi_count_data_genus_all_summarized_WIS.csv")
wzBacteriaOnlyGenusIntersected %>% write.csv(file = "Interim_data/genus_data_WIS_full_revised_bacteria_count_data_genus_all_summarized_WIS.csv")
wzMetaGenusIntersected %>% write.csv(file = "Interim_data/genus_data_WIS_full_revised_metadata_for_fungi_and_bacteria_genus_all_summarized_WIS.csv")
wzBacteriaOnlyTaxaTableAll %>% write.csv(file = "Interim_data/genus_data_WIS_full_revised_taxa_table_for_bacteria_only_genus_all_summarized_WIS.csv")
wzFungiTaxaTableAll %>% write.csv(file = "Interim_data/genus_data_WIS_full_revised_taxa_table_for_fungi_only_genus_all_summarized_WIS.csv")

wzFungiGenusIntersected <- read.csv("Interim_data/genus_data_WIS_full_revised_fungi_count_data_genus_all_summarized_WIS.csv", row.names = 1)
wzMetaGenusIntersected <- read.csv("Interim_data/genus_data_WIS_full_revised_metadata_for_fungi_and_bacteria_genus_all_summarized_WIS.csv", row.names = 1)
wzFungiTaxaTableAll <- read.csv("Interim_data/genus_data_WIS_full_revised_taxa_table_for_fungi_only_genus_all_summarized_WIS.csv", row.names = 1)

#----------------------------------------------------------#
# Testing machine learning
#----------------------------------------------------------#

source("00-Functions.R") # for wzML1VsAll10k() function

table(metaWzFungiFree$tissue)
# breast     lung melanoma    colon      gbm    ovary     bone pancreas 
# 387      258       68       34       34       48       36       32 
wzFungiML_BRCA <- wzML1VsAll10k(metaData = metaWzFungiBacteriaFree_Species, 
                             snmData = wzBacteriaOnlyFreeCounts_Species,
                            modelType = "gbm", trainTestFlag = FALSE, numKFold = 10, varImpFlag = FALSE,
                            dzOfInterest = "breast")
wzFungiML_LC <- wzML1VsAll10k(metaData = metaWzFungiBacteriaFree, 
                                snmData = wzFungiBacteriaFreeCounts,
                                modelType = "gbm", trainTestFlag = FALSE, numKFold = 10, varImpFlag = TRUE,
                                dzOfInterest = "lung")
wzFungiML_SKCM <- wzML1VsAll10k(metaData = metaWzFungiBacteriaFree, 
                              snmData = wzFungiBacteriaFreeCounts,
                              modelType = "gbm", trainTestFlag = FALSE, numKFold = 10, varImpFlag = TRUE,
                              dzOfInterest = "melanoma")
wzFungiML_CRC <- wzML1VsAll10k(metaData = metaWzFungiBacteriaFree, 
                              snmData = wzFungiBacteriaFreeCounts,
                              modelType = "gbm", trainTestFlag = FALSE, numKFold = 10, varImpFlag = TRUE,
                              dzOfInterest = "colon")
wzFungiML_GBM <- wzML1VsAll10k(metaData = metaWzFungiBacteriaFree, 
                              snmData = wzFungiBacteriaFreeCounts,
                              modelType = "gbm", trainTestFlag = FALSE, numKFold = 10, varImpFlag = TRUE,
                              dzOfInterest = "gbm")
wzFungiML_OV <- wzML1VsAll10k(metaData = metaWzFungiBacteriaFree, 
                              snmData = wzFungiBacteriaFreeCounts,
                              modelType = "gbm", trainTestFlag = FALSE, numKFold = 10, varImpFlag = TRUE,
                              dzOfInterest = "ovary")
wzFungiML_SARC <- wzML1VsAll10k(metaData = metaWzFungiBacteriaFree, 
                              snmData = wzFungiBacteriaFreeCounts,
                              modelType = "gbm", trainTestFlag = FALSE, numKFold = 10, varImpFlag = TRUE,
                              dzOfInterest = "bone")
wzFungiML_PAAD <- wzML1VsAll10k(metaData = metaWzFungiBacteriaFree, 
                              snmData = wzFungiBacteriaFreeCounts,
                              modelType = "gbm", trainTestFlag = FALSE, numKFold = 10, varImpFlag = TRUE,
                              dzOfInterest = "pancreas")

#----------------------------------------------------------#
# Plot machine learning performance
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_WIS <- read.csv("Interim_data/rep_perfWIS_10k_rep1_fungi_bacteria_ALL_5Apr22.csv", stringsAsFactors = FALSE)
abbreviationsWIS <- read.csv("Supporting_data/wis_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_WIS$abbrev <- abbreviationsWIS[mlPerfAll10k_WIS$diseaseType,"abbrev"]
mlPerfAll10k_WIS <- mlPerfAll10k_WIS[,!(colnames(mlPerfAll10k_WIS) == "X")]
colnames(mlPerfAll10k_WIS)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_WIS$nullAUPR <- ifelse(mlPerfAll10k_WIS$minorityClassName == "nat" | mlPerfAll10k_WIS$minorityClassName == "normal",
                                        yes=mlPerfAll10k_WIS$majorityClassSize/(mlPerfAll10k_WIS$minorityClassSize+mlPerfAll10k_WIS$majorityClassSize),
                                        no=mlPerfAll10k_WIS$minorityClassSize/(mlPerfAll10k_WIS$minorityClassSize+mlPerfAll10k_WIS$majorityClassSize))
mlPerfAll10k_WIS$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
mlPerfAll10k_WIS$datasetName[mlPerfAll10k_WIS$datasetName == "wzFungiFreeCountsFilt"] <- "Fungi counts"
mlPerfAll10k_WIS$datasetName[mlPerfAll10k_WIS$datasetName == "ra_wzFungiFreeCountsFilt"] <- "Fungi relative abundances"
mlPerfAll10k_WIS$datasetName[mlPerfAll10k_WIS$datasetName == "binary_wzFungiFreeCountsFilt"] <- "Fungi presence-absence"
mlPerfAll10k_WIS$datasetName[mlPerfAll10k_WIS$datasetName == "count_binary_wzFungiFreeCountsFilt"] <- "Fungi counts & presence-absence"

mlPerfAll10k_WIS$datasetName[mlPerfAll10k_WIS$datasetName == "wzBacteriaFreeCountsFilt"] <- "Bacteria counts"
mlPerfAll10k_WIS$datasetName[mlPerfAll10k_WIS$datasetName == "ra_wzBacteriaFreeCountsFilt"] <- "Bacteria relative abundances"
mlPerfAll10k_WIS$datasetName[mlPerfAll10k_WIS$datasetName == "binary_wzBacteriaFreeCountsFilt"] <- "Bacteria presence-absence"
mlPerfAll10k_WIS$datasetName[mlPerfAll10k_WIS$datasetName == "count_binary_wzBacteriaFreeCountsFilt"] <- "Bacteria counts & presence-absence"

mlPerfAll10k_WIS$datasetName[mlPerfAll10k_WIS$datasetName == "wzFungiBacteriaFreeCountsFilt"] <- "Fungi+Bacteria counts"
mlPerfAll10k_WIS$datasetName[mlPerfAll10k_WIS$datasetName == "ra_wzFungiBacteriaFreeCountsFilt"] <- "Fungi+Bacteria relative abundances"
mlPerfAll10k_WIS$datasetName[mlPerfAll10k_WIS$datasetName == "binary_wzFungiBacteriaFreeCountsFilt"] <- "Fungi+Bacteria presence-absence"
mlPerfAll10k_WIS$datasetName[mlPerfAll10k_WIS$datasetName == "count_binary_wzFungiBacteriaFreeCountsFilt"] <- "Fungi+Bacteria counts & presence-absence"

mlPerfAll10k_WIS$datasetName <- factor(mlPerfAll10k_WIS$datasetName,
                                           levels = c("Fungi counts",
                                                      "Fungi relative abundances",
                                                      "Fungi presence-absence",
                                                      "Fungi counts & presence-absence",
                                                      "Bacteria counts",
                                                      "Bacteria relative abundances",
                                                      "Bacteria presence-absence",
                                                      "Bacteria counts & presence-absence",
                                                      "Fungi+Bacteria counts",
                                                      "Fungi+Bacteria relative abundances",
                                                      "Fungi+Bacteria presence-absence",
                                                      "Fungi+Bacteria counts & presence-absence"))
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
# All taxa levels and their intersections
mlPerfAll10k_WIS %>%
  filter(sampleType == "tumor") %>%
  filter(grepl("counts$",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>% 
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted") +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WIS | Primary Tumor | 1 Vs All | Free Rank") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Supplementary_Figures/ml_WIS_tumor_freerank_counts.svg", dpi = "retina",
         width = 10, height = 6, units = "in")

mlPerfAll10k_WIS %>%
  filter(sampleType == "tumor") %>%
  filter(!grepl("\\&",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>% 
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted") +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WIS | Primary Tumor | 1 Vs All | Free Rank") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_igv(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Supplementary_Figures/ml_WIS_tumor_freerank_everything.svg", dpi = "retina",
         width = 14, height = 6, units = "in")

#-------------------------Plot primary tumor vs nat/normal-------------------------#

# T vs NAT
mlPerfAll10k_WIS %>%
  filter(sampleType == "tumor vs nat") %>%
  # filter(grepl("counts$",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>% 
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted") +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WIS | Tumor vs NAT | Free Rank") + theme(plot.title = element_text(hjust = 0.5), legend.position="right") +
  rotate_x_text(90) + scale_color_igv(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Supplementary_Figures/ml_WIS_tumor_vs_nat_freerank_everything.svg", dpi = "retina",
         width = 10, height = 6, units = "in")

# T vs NAT
mlPerfAll10k_WIS %>%
  filter(sampleType == "tumor vs normal") %>%
  filter(grepl("counts$",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>% 
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="darkgray",lty="dotted") +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WIS | Tumor vs Normal | Free Rank") + theme(plot.title = element_text(hjust = 0.5), legend.position="right") +
  rotate_x_text(0) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Supplementary_Figures/ml_WIS_tumor_vs_normal_freerank_counts.svg", dpi = "retina",
         width = 6, height = 3, units = "in")

#-------------------------Plot aggregated boxplots-------------------------#
require(rstatix)
mlPerfAll10k_WIS %>%
  filter(sampleType == "tumor") %>%
  filter(grepl("counts$",datasetName)) %>%
  distinct() %>% droplevels() %>%
  group_by(datasetName) %>% data.frame() %>%
  wilcox_test(AUROC ~ datasetName, exact = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "datasetName") -> roc.tumor.stat.test

mlPerfAll10k_WIS %>%
  filter(sampleType == "tumor") %>%
  filter(grepl("counts$",datasetName)) %>%
  # filter(grepl("counts$",datasetName)) %>%
  distinct() %>% droplevels() %>%
  ggboxplot(x = "datasetName",
            y = "AUROC",
            fill = "datasetName",
            legend = "none",
            xlab = "WIS features",
            # title = "Aggregated primary\ntumor 1-Vs-All\nperformance",
            add = "jitter",
            add.params = list(alpha = 0.4),
            palette = "nejm",
            notch = TRUE) + 
  stat_pvalue_manual(roc.tumor.stat.test, label = "q = {p.adj}", y.position = c(1.03, 1.12, 1.06), size = 3) +
  # scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(breaks = seq(0, 1.1, by = 0.1), limits = c(0,1.3)) +
  rotate_x_text(30) + 
  theme(plot.title = element_text(hjust = 0.5)) -> aggregatedTumorROC

mlPerfAll10k_WIS %>%
  filter(sampleType == "tumor") %>%
  filter(grepl("counts$",datasetName)) %>%
  distinct() %>% droplevels() %>%
  group_by(datasetName) %>% data.frame() %>%
  wilcox_test(AUPR ~ datasetName, exact = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "datasetName") -> pr.tumor.stat.test

mlPerfAll10k_WIS %>%
  filter(sampleType == "tumor") %>%
  filter(grepl("counts$",datasetName)) %>%
  # filter(grepl("counts$",datasetName)) %>%
  distinct() %>% droplevels() %>%
  ggboxplot(x = "datasetName",
            y = "AUPR",
            fill = "datasetName",
            legend = "none",
            xlab = "WIS features",
            # title = "Aggregated primary\ntumor 1-Vs-All\nperformance",
            add = "jitter",
            add.params = list(alpha = 0.4),
            palette = "nejm",
            notch = TRUE) + 
  stat_pvalue_manual(pr.tumor.stat.test, label = "q = {p.adj}", y.position = c(1.03, 1.12, 1.06), size = 3) +
  # scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(breaks = seq(0, 1.1, by = 0.1), limits = c(0,1.3)) +
  rotate_x_text(30) + 
  theme(plot.title = element_text(hjust = 0.5)) -> aggregatedTumorPR

combinedAggregatedTumorPlot <- ggarrange(aggregatedTumorROC, aggregatedTumorPR, ncol = 2) 
combinedAggregatedTumorPlotAnnotated <- annotate_figure(combinedAggregatedTumorPlot, 
                                         top = text_grob("Aggregated WIS primary tumor 1-Vs-All performance", 
                                                         color = "black", face = "bold", size = 14))
print(combinedAggregatedTumorPlotAnnotated)
ggsave(filename = paste0("Figures/Main_Figures/ml_WIS_aggregated_tumor_1VsAll_counts_auroc_aupr.svg"),
       plot = combinedAggregatedTumorPlotAnnotated,
       dpi = "retina", units = "in", width = 6, height = 5)

#-------------------------Regress ML performance vs minority class sizes-------------------------#

ptWISPerfSummarizedMinClassSize <- mlPerfAll10k_WIS %>%
  filter(sampleType == "tumor") %>%
  # filter(grepl("counts$",datasetName)) %>%
  filter(grepl("Fungi counts$",datasetName)) %>%
  distinct() %>% droplevels() %>%
  group_by(sampleType, diseaseType) %>% 
  mutate(avgROC = mean(AUROC), avgPR = mean(AUPR)) %>%
  select(sampleType, diseaseType, abbrev, avgROC, avgPR, minorityClassSize, majorityClassSize) %>% 
  mutate(logAvgMinClass = log10(minorityClassSize), ratioMinMajClass=minorityClassSize/majorityClassSize) %>%
  unique() %>% data.frame()

## Regress using minority class size
summary(lm(avgPR ~ minorityClassSize, ptWISPerfSummarizedMinClassSize))
summary(lm(avgROC ~ minorityClassSize, ptWISPerfSummarizedMinClassSize))

require(ggrepel)
ptWISPerfSummarizedMinClassSize %>%
  ggplot(aes(x=minorityClassSize, y=avgPR, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red") + 
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  geom_text_repel() +
  labs(x = "Minority class size (per cancer type)", y = "Average AUPR", 
       title = "Correlating AUPR and minority class size\nfor primary tumor WIS models")
ggsave(filename = "Figures/Supplementary_Figures/aupr_vs_minority_class_size_ml_WIS_tumor.svg",
         dpi = "retina", units = "in", width = 4, height = 4)

ptWISPerfSummarizedMinClassSize %>%
  ggplot(aes(x=minorityClassSize, y=avgROC, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.y = 0.92) + 
  geom_text_repel() +
  theme_bw()+ theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  labs(x = "Minority class size (per cancer type)", y = "Average AUROC", 
       title = "Correlating AUROC and minority class size\nfor primary tumor WIS models")
ggsave(filename = "Figures/Supplementary_Figures/auroc_vs_minority_class_size_ml_WIS_tumor.svg",
         dpi = "retina", units = "in", width = 4, height = 4)

## Regress using minority class *ratio*
summary(lm(avgPR ~ ratioMinMajClass, ptWISPerfSummarizedMinClassSize))
summary(lm(avgROC ~ ratioMinMajClass, ptWISPerfSummarizedMinClassSize))

ptWISPerfSummarizedMinClassSize %>%
  ggplot(aes(x=ratioMinMajClass, y=avgPR, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.y = 1.1) + 
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  geom_text_repel() +
  labs(x = "Ratio of minority to majority class size\n(per cancer type)", y = "Average AUPR", 
       title = "Correlating AUPR and ratio of minority to majority class size\nfor primary tumor WIS models")
ggsave(filename = "Figures/Supplementary_Figures/aupr_vs_minority_to_majority_class_size_ml_WIS_tumor.svg",
         dpi = "retina", units = "in", width = 4, height = 4)

ptWISPerfSummarizedMinClassSize %>%
  ggplot(aes(x=ratioMinMajClass, y=avgROC, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.y = 1.05) + 
  geom_text_repel() +
  theme_bw()+ theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  labs(x = "Ratio of minority to majority class size\n(per cancer type)", y = "Average AUROC", 
       title = "Correlating AUROC and ratio of minority to majority class size\nfor primary tumor WIS models")
ggsave(filename = "Figures/Supplementary_Figures/auroc_vs_minority_to_majority_class_size_ml_WIS_tumor.svg",
         dpi = "retina", units = "in", width = 4, height = 4)

