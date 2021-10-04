#-----------------------------------------------------------------------------
# 12-Plasma-validation-cohort-2-Hopkins.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Import rep200 data for plasma validation cohort #2 (Cristiano et al. 2019. Nature)
# - Remove contaminants identified in TCGA and plasma validation cohort #1
# - Perform machine learning between cancer vs. healthy groups and between stages
#-----------------------------------------------------------------------------

#----------------------------------------------------------#
# Load environments
#----------------------------------------------------------#

# Load dependencies
require(devtools)
require(doMC)
require(plyr)
require(dplyr)
require(reshape2)
require(ggpubr)
require(ggsci)
require(tibble)
require(phyloseq)
require(microbiome)
require(vegan)
require(biomformat)
require(rhdf5)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Rep200 fungal species identification
#----------------------------------------------------------#

rep200TaxSplit <- read.csv("Supporting_data/rep200_lineage_map_split.csv", stringsAsFactors = FALSE, row.names = 1)
rep200Kingdoms <- read.csv("Supporting_data/rep200_gOTU_kingdom_mapping.csv", stringsAsFactors = FALSE)
rep200Kingdoms_Fungi <- rep200Kingdoms[which(rep200Kingdoms$kingdom == "fungi"),]
rep200Kingdoms_Bacteria <- rep200Kingdoms[which(rep200Kingdoms$kingdom == "bacteria"),]

rep200TaxSplit_Fungi <- rep200TaxSplit[rownames(rep200TaxSplit) %in% rep200Kingdoms_Fungi$genomeID,,drop=FALSE]
dim(rep200TaxSplit_Fungi) # 320   7

rep200TaxSplit_Bacteria <- rep200TaxSplit[rownames(rep200TaxSplit) %in% rep200Kingdoms_Bacteria$genomeID,,drop=FALSE]
dim(rep200TaxSplit_Bacteria) # 11080     7

fungiOGUs <- rownames(rep200TaxSplit_Fungi)
bacteriaOGUs <- rownames(rep200TaxSplit_Bacteria)

#----------------------------------------------------------#
# Microbial data import
#----------------------------------------------------------#
# The biom table and metadata were loaded from the following 
# Qiita analysis (to be made public upon publication)
# https://qiita.ucsd.edu/analysis/description/47106/

## Import Qiita metadata
metaCristiano <- read.csv("Input_data/qiita_metadata_plasma_validation_cohort_2.csv", stringsAsFactors = FALSE, row.names = 1)
## Filter based on pre-treatment / treatment naive status. 
# For patients with more than 1 pre-treatment sample, take the earliest one.
metaCristianoTxNaive <- metaCristiano %>% filter(grepl("Preoperative treatment naive|Pre-treatment", subject_id)) %>% droplevels()
metaCristianoTxNaive$patientWithoutDay <- gsub("_D.+","",metaCristianoTxNaive$alias)
metaCristianoTxNaive$dayPreTx <- as.numeric(gsub("CG.+ Day ","",metaCristianoTxNaive$subject_id))

# The following code finds patients with more than 1 sample, finds which sample was taken earliest,
# and then creates a list of the samples to remove
patientsWithMultipleSamples <- metaCristianoTxNaive[duplicated(metaCristianoTxNaive$patientWithoutDay)|duplicated(metaCristianoTxNaive$patientWithoutDay, fromLast = TRUE),
                     c("patientWithoutDay","dayPreTx")] 
patientsWithMultipleSamplesFixed <- patientsWithMultipleSamples %>% 
  arrange(patientWithoutDay, dayPreTx) %>%
  distinct(patientWithoutDay, .keep_all = TRUE)
samplesToRemove <- rownames(patientsWithMultipleSamples[!(rownames(patientsWithMultipleSamples) %in% rownames(patientsWithMultipleSamplesFixed)),])

# Create final filtered metadata table
metaCristianoTxNaiveFilt <- droplevels(metaCristianoTxNaive[!(rownames(metaCristianoTxNaive) %in% samplesToRemove),])
metaCristianoTxNaiveFilt$phenotype <- factor(metaCristianoTxNaiveFilt$phenotype)
metaCristianoTxNaiveFilt$cancer_status <- factor(ifelse(metaCristianoTxNaiveFilt$phenotype == "Healthy", 
                                                 yes="Healthy",no="Cancer"), levels = c("Cancer","Healthy"))
metaCristianoTxNaiveFilt$stage[metaCristianoTxNaiveFilt$phenotype == "Healthy"] <- "Healthy"
## Import read count data
cristiano_rep200Data_BIOM <- read_biom(biom_file = "Input_data/Qiita_results/rep200_OGU_plasma_validation_cohort_2.biom")
cristiano_rep200Data <- t(as(biom_data(cristiano_rep200Data_BIOM), "matrix"))
cristiano_rep200Data_Filt <- cristiano_rep200Data[rownames(metaCristianoTxNaiveFilt),]
dim(cristiano_rep200Data_Filt) # 491 7418

## Extract only fungal features
cristiano_rep200Data_Filt_Fungi <- cristiano_rep200Data_Filt[,colnames(cristiano_rep200Data_Filt) %in% fungiOGUs]
dim(cristiano_rep200Data_Filt_Fungi) # 491 296
## Extract only bacterial features
cristiano_rep200Data_Filt_Bacteria <- cristiano_rep200Data_Filt[,colnames(cristiano_rep200Data_Filt) %in% bacteriaOGUs]
dim(cristiano_rep200Data_Filt_Bacteria) # 491 6908

#----------------------------------------------------------------------------------------------#
# Create fungi phyloseq object to summarize reads at various taxa levels and match to Weizmann features
#----------------------------------------------------------------------------------------------#

rep200TaxSplit_Fungi_Paired_to_Weizmann <- read.csv("Supporting_data/rep200TaxSplit_Fungi_Paired_To_Weizmann_Final.csv", stringsAsFactors = FALSE, row.names = 1)

# Build phyloseq object
psCristianoFungi <- phyloseq(otu_table(cristiano_rep200Data_Filt_Fungi, taxa_are_rows = FALSE), 
                             tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), 
                             sample_data(metaCristianoTxNaiveFilt))

## Aggregate counts
psCristianoFungi_phylum = aggregate_taxa(psCristianoFungi, "phylum")
psCristianoFungi_class = aggregate_taxa(psCristianoFungi, "class")
psCristianoFungi_order = aggregate_taxa(psCristianoFungi, "order")
psCristianoFungi_family = aggregate_taxa(psCristianoFungi, "family")
psCristianoFungi_genus = aggregate_taxa(psCristianoFungi, "genus")
psCristianoFungi_species = cristiano_rep200Data_Filt_Fungi
colnames(psCristianoFungi_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(psCristianoFungi_species),"species"]

## Load shared features with Weizmann
load("Interim_data/shared_fungi_features_at_each_taxa_level_13Sep21.RData")

## Subset taxa to those shared with the Weizmann data
psCristianoFungi_phylum_shared <- subset_taxa(psCristianoFungi_phylum, phylum %in% sharedPhylum)
psCristianoFungi_class_shared <- subset_taxa(psCristianoFungi_class, class %in% sharedClass)
psCristianoFungi_order_shared <- subset_taxa(psCristianoFungi_order, order %in% sharedOrder)
psCristianoFungi_family_shared <- subset_taxa(psCristianoFungi_family, family %in% sharedFamily)
psCristianoFungi_genus_shared <- subset_taxa(psCristianoFungi_genus, genus %in% sharedGenus)
psCristianoFungi_species_shared <- psCristianoFungi_species[,colnames(psCristianoFungi_species) %in% sharedSpecies]

## Create data.frames of summarized data with and without feature intersection
# Without feature intersection
cristianoRep200FungiPhylum <- data.frame(t(otu_table(psCristianoFungi_phylum)))
cristianoRep200FungiClass <- data.frame(t(otu_table(psCristianoFungi_class)))
cristianoRep200FungiOrder <- data.frame(t(otu_table(psCristianoFungi_order)))
cristianoRep200FungiFamily <- data.frame(t(otu_table(psCristianoFungi_family)))
cristianoRep200FungiGenus <- data.frame(t(otu_table(psCristianoFungi_genus)))
cristianoRep200FungiSpecies <- data.frame(psCristianoFungi_species)
cristianoRep200FungiSpecies[1:3,1:3]
# With feature intersection
cristianoRep200FungiPhylumShared <- data.frame(t(otu_table(psCristianoFungi_phylum_shared)))
cristianoRep200FungiClassShared <- data.frame(t(otu_table(psCristianoFungi_class_shared)))
cristianoRep200FungiOrderShared <- data.frame(t(otu_table(psCristianoFungi_order_shared)))
cristianoRep200FungiFamilyShared <- data.frame(t(otu_table(psCristianoFungi_family_shared)))
cristianoRep200FungiGenusShared <- data.frame(t(otu_table(psCristianoFungi_genus_shared)))
cristianoRep200FungiSpeciesShared <- data.frame(psCristianoFungi_species_shared)
cristianoRep200FungiGenusShared[1:3,1:3]

#-----------------------------------------------#
# Decontaminate using TCGA and Nature Val data
#-----------------------------------------------#

### NEED TO IMPLEMENT THIS STEP AFTER WRITING THE VAL #1 SCRIPT (16 SEP 21) ###

load("Interim_data/contaminants_fungi_OGUs_UCSD_TCGA_25Sep21.RData") # for "contaminantsFungiUCSDAndPlateCenterTCGA" object
cristiano_rep200Data_Filt_Fungi_Decontam <- cristiano_rep200Data_Filt_Fungi[,!(colnames(cristiano_rep200Data_Filt_Fungi) %in% contaminantsFungiUCSDAndPlateCenterTCGA)]

#----------------------------------------------------------------------------------------------#
# Create bacteria phyloseq object to compare performance of fungi to bacteria
#----------------------------------------------------------------------------------------------#
# Build phyloseq object
psCristianoBacteria <- phyloseq(otu_table(cristiano_rep200Data_Filt_Bacteria, taxa_are_rows = FALSE), 
                             tax_table(as.matrix(rep200TaxSplit_Bacteria)), 
                             sample_data(metaCristianoTxNaiveFilt))

## Aggregate counts
psCristianoBacteria_phylum = aggregate_taxa(psCristianoBacteria, "Phylum")
psCristianoBacteria_class = aggregate_taxa(psCristianoBacteria, "Class")
psCristianoBacteria_order = aggregate_taxa(psCristianoBacteria, "Order")
psCristianoBacteria_family = aggregate_taxa(psCristianoBacteria, "Family")
psCristianoBacteria_genus = aggregate_taxa(psCristianoBacteria, "Genus")
psCristianoBacteria_species = cristiano_rep200Data_Filt_Bacteria
colnames(psCristianoBacteria_species) <- rep200TaxSplit_Bacteria[colnames(psCristianoBacteria_species),"Species"]

## Create data.frames of summarized data
cristianoRep200BacteriaPhylum <- data.frame(t(otu_table(psCristianoBacteria_phylum)))
cristianoRep200BacteriaClass <- data.frame(t(otu_table(psCristianoBacteria_class)))
cristianoRep200BacteriaOrder <- data.frame(t(otu_table(psCristianoBacteria_order)))
cristianoRep200BacteriaFamily <- data.frame(t(otu_table(psCristianoBacteria_family)))
cristianoRep200BacteriaGenus <- data.frame(t(otu_table(psCristianoBacteria_genus)))
cristianoRep200BacteriaSpecies <- data.frame(psCristianoBacteria_species)
cristianoRep200BacteriaSpecies[1:3,1:3]

#-----------------------------------------#
# Machine learning: Healthy vs. Cancer
# - Full rep200 vs. 
# - Bacteria only vs. 
# - Fungi only vs. 
# - Intersected fungi only
#-----------------------------------------#

source("00-Functions.R") # for mlCristiano() function

cristianoCancerVsHealthy_FullRep200 <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                   countData = cristiano_rep200Data_Filt,
                                                   dataString = "full_rep200",
                                                   col2Predict = "cancer_status",
                                                   varImpFlag = FALSE)
cristianoCancerVsHealthy_BacteriaOnly <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                     countData = cristianoRep200BacteriaSpecies,
                                                     dataString = "bacteria_only",
                                                     col2Predict = "cancer_status",
                                                     varImpFlag = FALSE)
cristianoCancerVsHealthy_FungiOnly <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                  countData = cristianoRep200FungiSpecies,
                                                  dataString = "fungi_only",
                                                  col2Predict = "cancer_status",
                                                  varImpFlag = FALSE)
cristianoCancerVsHealthy_FungiDecontam <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                      countData = cristiano_rep200Data_Filt_Fungi_Decontam,
                                                      dataString = "fungi_decontam",
                                                      col2Predict = "cancer_status",
                                                      varImpFlag = FALSE)
cristianoCancerVsHealthy_FungiIntersect <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                       countData = cristianoRep200FungiSpeciesShared,
                                                       dataString = "fungi_intersected_with_Weizmann",
                                                       col2Predict = "cancer_status",
                                                       varImpFlag = FALSE)
cristianoCancerVsHealthyResults <- rbind(cristianoCancerVsHealthy_FullRep200$rep_perf,
                                         cristianoCancerVsHealthy_BacteriaOnly$rep_perf,
                                         cristianoCancerVsHealthy_FungiDecontam$rep_perf,
                                         cristianoCancerVsHealthy_FungiOnly$rep_perf,
                                         cristianoCancerVsHealthy_FungiIntersect$rep_perf)

colnames(cristianoCancerVsHealthyResults)[1:2] <- c("AUROC","AUPR")
cristianoCancerVsHealthyResults$dataString[cristianoCancerVsHealthyResults$dataString=="full_rep200"] <- "Full dataset (rep200)"
cristianoCancerVsHealthyResults$dataString[cristianoCancerVsHealthyResults$dataString=="bacteria_only"] <- "Bacteria all (rep200)"
cristianoCancerVsHealthyResults$dataString[cristianoCancerVsHealthyResults$dataString=="fungi_only"] <- "Fungi all (rep200)"
cristianoCancerVsHealthyResults$dataString[cristianoCancerVsHealthyResults$dataString=="fungi_decontam"] <- "Fungi decontaminated (rep200)"
cristianoCancerVsHealthyResults$dataString[cristianoCancerVsHealthyResults$dataString=="fungi_intersected_with_Weizmann"] <- "Intersected fungi only (rep200+Weizmann)"
cristianoCancerVsHealthyResults$dataString <- factor(cristianoCancerVsHealthyResults$dataString, levels = c("Full dataset (rep200)",
                                                                                                            "Bacteria all (rep200)",
                                                                                                            "Fungi all (rep200)",
                                                                                                            "Fungi decontaminated (rep200)",
                                                                                                            "Intersected fungi only (rep200+Weizmann)"))
require(ggrepel)
source("Supporting_scripts/S05-SummarySE.R")
# Facet by AUROC vs AUPR; color by full dataset vs fungi (dataString)
cristianoCancerVsHealthyResults %>%
  reshape2::melt(id.vars = c("dataString","rep","diseaseType", "col2Predict")) %>%
  summarySE(measurevar = "value", groupvars = c("diseaseType","dataString","variable")) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=dataString)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=3) + 
  ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + ggtitle("Plasma validation cohort #2 cancer vs. healthy (231 cancer | 260 healthy)") +
  theme(plot.title = element_text(hjust = 0), legend.position = "right") + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "Cancer vs. Healthy") +
  rotate_x_text(0) + scale_color_nejm(name = "Dataset") + geom_hline(yintercept = 1, linetype="dashed") + 
  geom_label_repel(aes(label=round(value,2)), size=3, box.padding = 0.75, point.padding = 6.5, direction = "y", position = position_dodge(width = 1), show.legend = FALSE) +
  ggsave("Figures/Figure_5/figure_5_F_cristiano_all_cancer_vs_healthy_all_datasets_16Sep21.jpeg", 
         dpi = "retina", width = 8, height = 4.5, units = "in")
cristianoCancerVsHealthyResults %>% write.csv("Figures_data/Figure_5/figure_5_F_cristiano_all_cancer_vs_healthy_all_datasets_16Sep21.csv")

#---------------------------------------------------------------------------------------------------------------------------#
# ML perf each cancer vs. healthy -- 4 datasets
#---------------------------------------------------------------------------------------------------------------------------#
source("00-Functions.R") # for the ml1VsAllCristiano10kRep1_Iterate() function
cristianoPerCancerIterate <- ml1VsAllCristiano10kRep1_Iterate()

cristianoPerCancerIterate_results <- cristianoPerCancerIterate$rep_perf
colnames(cristianoPerCancerIterate_results)[1:2] <- c("AUROC","AUPR")
cristianoPerCancerIterate_results$dataString[cristianoPerCancerIterate_results$dataString=="full_rep200"] <- "Full dataset (rep200)"
cristianoPerCancerIterate_results$dataString[cristianoPerCancerIterate_results$dataString=="bacteria_only"] <- "Bacteria all (rep200)"
cristianoPerCancerIterate_results$dataString[cristianoPerCancerIterate_results$dataString=="fungi_only"] <- "Fungi all (rep200)"
cristianoPerCancerIterate_results$dataString[cristianoPerCancerIterate_results$dataString=="fungi_decontam"] <- "Fungi decontaminated (rep200)"
cristianoPerCancerIterate_results$dataString[cristianoPerCancerIterate_results$dataString=="fungi_intersected_with_Weizmann"] <- "Intersected fungi only (rep200+Weizmann)"
cristianoPerCancerIterate_results$dataString <- factor(cristianoPerCancerIterate_results$dataString, levels = c("Full dataset (rep200)",
                                                                                                            "Bacteria all (rep200)",
                                                                                                            "Fungi all (rep200)",
                                                                                                            "Fungi decontaminated (rep200)",
                                                                                                            "Intersected fungi only (rep200+Weizmann)"))
cristianoPerCancerIterate_results$diseaseType[cristianoPerCancerIterate_results$diseaseType=="Bile Duct Cancer"] <- "Bile Duct"
cristianoPerCancerIterate_results$diseaseType[cristianoPerCancerIterate_results$diseaseType=="Breast Cancer"] <- "Breast"
cristianoPerCancerIterate_results$diseaseType[cristianoPerCancerIterate_results$diseaseType=="Colorectal Cancer"] <- "Colorectal"
cristianoPerCancerIterate_results$diseaseType[cristianoPerCancerIterate_results$diseaseType=="Gastric cancer"] <- "Gastric"
cristianoPerCancerIterate_results$diseaseType[cristianoPerCancerIterate_results$diseaseType=="Lung Cancer"] <- "Lung"
cristianoPerCancerIterate_results$diseaseType[cristianoPerCancerIterate_results$diseaseType=="Ovarian Cancer"] <- "Ovarian"
cristianoPerCancerIterate_results$diseaseType[cristianoPerCancerIterate_results$diseaseType=="Pancreatic Cancer"] <- "Pancreatic"

save(cristianoPerCancerIterate_results,
     file = "Interim_data/cristianoPerCancerIterate_results_25Sep21.RData")

source("Supporting_scripts/S05-SummarySE.R")
# Facet by AUROC vs AUPR; color by full dataset vs fungi (dataString)
cristianoPerCancerIterate_results %>%
  reshape2::melt(id.vars = c("dataString","rep","diseaseType", "col2Predict")) %>%
  summarySE(measurevar = "value", groupvars = c("diseaseType","dataString","variable")) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=dataString)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=3) + xlab("Individual cancer type versus all healthy samples") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + ggtitle("Plasma validation cohort #2: Per cancer type vs. healthy performance") +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  rotate_x_text(0) + scale_color_nejm(name = "Dataset") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_5/figure_5_G_cristiano_per_cancer_type_vs_healthy_all_datasets_16Sep21.png", 
         dpi = "retina", width = 12, height = 4.5, units = "in")
cristianoPerCancerIterate_results %>% write.csv("Figures_data/Figure_5/figure_5_G_cristiano_per_cancer_type_vs_healthy_all_datasets_16Sep21.csv")

#-----------------------------Test ML perf at varying stages levels-----------------------------#

source("00-Functions.R") # for the ml1VsAllCristiano10kRep1_Iterate_Stage() function
cristianoStageIterate <- ml1VsAllCristiano10kRep1_Iterate_Stage()

cristianoStageIterate_results <- cristianoStageIterate$rep_perf
colnames(cristianoStageIterate_results)[1:2] <- c("AUROC","AUPR")
cristianoStageIterate_results$dataString[cristianoStageIterate_results$dataString=="full_rep200"] <- "Full dataset (rep200)"
cristianoStageIterate_results$dataString[cristianoStageIterate_results$dataString=="bacteria_only"] <- "Bacteria all (rep200)"
cristianoStageIterate_results$dataString[cristianoStageIterate_results$dataString=="fungi_only"] <- "Fungi all (rep200)"
cristianoStageIterate_results$dataString[cristianoStageIterate_results$dataString=="fungi_decontam"] <- "Fungi decontaminated (rep200)"
cristianoStageIterate_results$dataString[cristianoStageIterate_results$dataString=="fungi_intersected_with_Weizmann"] <- "Intersected fungi only (rep200+Weizmann)"
cristianoStageIterate_results$dataString <- factor(cristianoStageIterate_results$dataString, levels = c("Full dataset (rep200)",
                                                                                                                "Bacteria all (rep200)",
                                                                                                                "Fungi all (rep200)",
                                                                                                                "Fungi decontaminated (rep200)",
                                                                                                                "Intersected fungi only (rep200+Weizmann)"))
cristianoStageIterate_results$stageNum[cristianoStageIterate_results$stageNum=="I"] <- "Stage I"
cristianoStageIterate_results$stageNum[cristianoStageIterate_results$stageNum=="II"] <- "Stage II"
cristianoStageIterate_results$stageNum[cristianoStageIterate_results$stageNum=="III"] <- "Stage III"
cristianoStageIterate_results$stageNum[cristianoStageIterate_results$stageNum=="IV"] <- "Stage IV"
cristianoStageIterate_results$stageNum <- ordered(cristianoStageIterate_results$stageNum,
                                                  levels = c("Stage I","Stage II","Stage III","Stage IV"))

save(cristianoStageIterate_results,
     file = "Interim_data/cristianoStageIterate_results_25Sep21.RData")

source("Supporting_scripts/S05-SummarySE.R")
# Facet by AUROC vs AUPR; color by full dataset vs fungi (dataString)
cristianoStageIterate_results %>%
  reshape2::melt(id.vars = c("dataString","rep","stageNum", "diseaseType", "col2Predict")) %>%
  summarySE(measurevar = "value", groupvars = c("stageNum","dataString","diseaseType","variable")) %>%
  ggplot(aes(stageNum,value, color=dataString)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=3) + xlab("All cancer types grouped within a particular stage versus all healthy samples") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + ggtitle("Plasma validation cohort #2: per cancer stage vs. healthy performance") +
  theme(plot.title = element_text(hjust = 0.5)) + #ylim(c(0,1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  rotate_x_text(0) + scale_color_nejm(name = "Dataset") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_5/figure_5_H_cristiano_per_cancer_stage_vs_healthy_all_datasets_17Sep21.jpeg", 
         dpi = "retina", width = 10, height = 6, units = "in")
cristianoStageIterate_results %>% write.csv("Figures_data/Figure_5/figure_5_H_cristiano_per_cancer_stage_vs_healthy_all_datasets_17Sep21.csv")

#-----------------------------ML perf at varying taxa levels-----------------------------#

source("00-Functions.R") # for mlCristiano() function

## Bacteria
cristianoBacteriaPhylum <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                    countData = cristianoRep200BacteriaPhylum,
                                    col2Predict = "cancer_status",
                                    dataString = "bacteria_phylum",
                                    varImpFlag = FALSE)
cristianoBacteriaClass <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                   countData = cristianoRep200BacteriaClass,
                                   col2Predict = "cancer_status",
                                   dataString = "bacteria_class",
                                   varImpFlag = FALSE)
cristianoBacteriaOrder <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                   countData = cristianoRep200BacteriaOrder,
                                   col2Predict = "cancer_status",
                                   dataString = "bacteria_order",
                                   varImpFlag = FALSE)
cristianoBacteriaFamily <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                    countData = cristianoRep200BacteriaFamily,
                                    col2Predict = "cancer_status",
                                    dataString = "bacteria_family",
                                    varImpFlag = FALSE)
cristianoBacteriaGenus <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                   countData = cristianoRep200BacteriaGenus,
                                   col2Predict = "cancer_status",
                                   dataString = "bacteria_genus",
                                   varImpFlag = TRUE)
cristianoBacteriaSpecies <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                     countData = cristianoRep200BacteriaSpecies,
                                     col2Predict = "cancer_status",
                                     dataString = "bacteria_species",
                                     varImpFlag = TRUE)

## Fungi
cristianoFungiPhylum <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                     countData = cristianoRep200FungiPhylum,
                                                     col2Predict = "cancer_status",
                                                    dataString = "fungi_phylum",
                                                     varImpFlag = FALSE)
cristianoFungiClass <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                    countData = cristianoRep200FungiClass,
                                                    col2Predict = "cancer_status",
                                   dataString = "fungi_class",
                                                    varImpFlag = FALSE)
cristianoFungiOrder <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                    countData = cristianoRep200FungiOrder,
                                                    col2Predict = "cancer_status",
                                   dataString = "fungi_order",
                                                    varImpFlag = FALSE)
cristianoFungiFamily <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                     countData = cristianoRep200FungiFamily,
                                                     col2Predict = "cancer_status",
                                    dataString = "fungi_family",
                                                     varImpFlag = FALSE)
cristianoFungiGenus <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                    countData = cristianoRep200FungiGenus,
                                                    col2Predict = "cancer_status",
                                   dataString = "fungi_genus",
                                                    varImpFlag = TRUE)
cristianoFungiSpecies <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                      countData = cristianoRep200FungiSpecies,
                                                      col2Predict = "cancer_status",
                                     dataString = "fungi_species",
                                                      varImpFlag = TRUE)

cristianoTaxaLevel_results <- cbind(rbind(cristianoBacteriaPhylum$rep_perf,
                                          cristianoBacteriaClass$rep_perf,
                                          cristianoBacteriaOrder$rep_perf,
                                          cristianoBacteriaFamily$rep_perf,
                                          cristianoBacteriaGenus$rep_perf,
                                          cristianoBacteriaSpecies$rep_perf,
                                          cristianoFungiPhylum$rep_perf,
                                          cristianoFungiClass$rep_perf,
                                          cristianoFungiOrder$rep_perf,
                                          cristianoFungiFamily$rep_perf,
                                          cristianoFungiGenus$rep_perf,
                                          cristianoFungiSpecies$rep_perf),
                                    taxaLevel = c(rep(c("Phylum","Class","Order","Family","Genus","Species"),each=10),
                                                  rep(c("Phylum","Class","Order","Family","Genus","Species"),each=10)))
colnames(cristianoTaxaLevel_results)[1:2] <- c("AUROC","AUPR")
cristianoTaxaLevel_results$dataString <- factor(ifelse(grepl("bacteria",cristianoTaxaLevel_results$dataString), yes = "Bacteria", no = "Fungi"),
                                                levels = c("Bacteria","Fungi"))
cristianoTaxaLevel_results$taxaLevel <- factor(cristianoTaxaLevel_results$taxaLevel,
                                               levels = c("Phylum","Class","Order","Family","Genus","Species"))

source("Supporting_scripts/S05-SummarySE.R")
# Facet by AUROC vs AUPR; color by full dataset vs fungi (dataString)
cristianoTaxaLevel_results %>%
  reshape2::melt(id.vars = c("dataString","rep","taxaLevel","diseaseType", "col2Predict")) %>%
  summarySE(measurevar = "value", groupvars = c("taxaLevel","diseaseType","dataString","variable")) %>%
  ggplot(aes(taxaLevel,value, color=dataString)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=3) + xlab("") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + ggtitle("Plasma validation cohort #2: All cancer (n=231) vs. healthy (n=260) at varying taxa levels\n(Only fungi overlapping Weizmann data at each taxa level used)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "top") +
  rotate_x_text(0) + scale_color_nejm(name = "Dataset") + geom_hline(yintercept = 1, linetype="dashed") + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  # geom_label_repel(aes(label=round(value,2)), size=3, box.padding = 1, point.padding = 10, position = position_dodge(width = 1), show.legend = FALSE) +
  ggsave("Figures/Supplementary_Figures/cristiano_all_cancer_vs_healthy_bacteria_vs_fungi_varying_taxa_levels_17Sep21.jpeg", dpi = "retina",
         width = 10, height = 4.5, units = "in")
cristianoTaxaLevel_results %>% write.csv("Figures_data/Supplementary_Figures/cristiano_all_cancer_vs_healthy_bacteria_vs_fungi_varying_taxa_levels_17Sep21.csv")


