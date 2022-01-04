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
patientsWithMultipleSamples <- metaCristianoTxNaive[duplicated(metaCristianoTxNaive$patientWithoutDay)|
                                                      duplicated(metaCristianoTxNaive$patientWithoutDay, fromLast = TRUE),
                                                    c("patientWithoutDay","dayPreTx")] 
patientsWithMultipleSamplesFixed <- patientsWithMultipleSamples %>% 
  arrange(patientWithoutDay, dayPreTx) %>%
  distinct(patientWithoutDay, .keep_all = TRUE)
samplesToRemove <- rownames(patientsWithMultipleSamples[!(rownames(patientsWithMultipleSamples) %in% 
                                                            rownames(patientsWithMultipleSamplesFixed)),])

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
# Decontaminate using TCGA data
#-----------------------------------------------#

load("Interim_data/decontamResultsV2_13Oct21.RData")

cristiano_rep200Data_Filt_Fungi_Decontam <- cristiano_rep200Data_Filt_Fungi[,colnames(cristiano_rep200Data_Filt_Fungi) %in%
                                                                              rownames(decontamResultsV2[decontamResultsV2$decision == "KEEP",])]

# Build phyloseq object
psCristianoFungiDecontam <- phyloseq(otu_table(cristiano_rep200Data_Filt_Fungi_Decontam, taxa_are_rows = FALSE), 
                             tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), 
                             sample_data(metaCristianoTxNaiveFilt))

## Aggregate counts
psCristianoFungiDecontam_phylum = aggregate_taxa(psCristianoFungiDecontam, "phylum")
psCristianoFungiDecontam_class = aggregate_taxa(psCristianoFungiDecontam, "class")
psCristianoFungiDecontam_order = aggregate_taxa(psCristianoFungiDecontam, "order")
psCristianoFungiDecontam_family = aggregate_taxa(psCristianoFungiDecontam, "family")
psCristianoFungiDecontam_genus = aggregate_taxa(psCristianoFungiDecontam, "genus")
psCristianoFungiDecontam_species = cristiano_rep200Data_Filt_Fungi_Decontam
colnames(psCristianoFungiDecontam_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[
  colnames(psCristianoFungiDecontam_species),"species"]

## Create data.frames of summarized data
cristianoRep200FungiDecontamPhylum <- data.frame(t(otu_table(psCristianoFungiDecontam_phylum)))
cristianoRep200FungiDecontamClass <- data.frame(t(otu_table(psCristianoFungiDecontam_class)))
cristianoRep200FungiDecontamOrder <- data.frame(t(otu_table(psCristianoFungiDecontam_order)))
cristianoRep200FungiDecontamFamily <- data.frame(t(otu_table(psCristianoFungiDecontam_family)))
cristianoRep200FungiDecontamGenus <- data.frame(t(otu_table(psCristianoFungiDecontam_genus)))
cristianoRep200FungiDecontamSpecies <- data.frame(psCristianoFungiDecontam_species)

# ### NEED TO IMPLEMENT THIS STEP AFTER WRITING THE VAL #1 SCRIPT (16 SEP 21) ###
# load("Interim_data/contaminants_fungi_OGUs_UCSD_TCGA_25Sep21.RData") # for "contaminantsFungiUCSDAndPlateCenterTCGA" object
# cristiano_rep200Data_Filt_Fungi_Decontam <- cristiano_rep200Data_Filt_Fungi[,!(colnames(cristiano_rep200Data_Filt_Fungi) %in% contaminantsFungiUCSDAndPlateCenterTCGA)]

#----------------------------------------------------------------------------------------------#
# Create bacteria phyloseq object to compare performance of fungi to bacteria
#----------------------------------------------------------------------------------------------#
# Modify bacterial taxa table
rep200TaxSplit_Bacteria_Formatted <- apply(rep200TaxSplit_Bacteria, 2, function(x) gsub("^[k|p|c|o|f|g|s]__","",x))
rep200TaxSplit_Bacteria_Formatted[rep200TaxSplit_Bacteria_Formatted == ""] <- "other"

# Build phyloseq object
psCristianoBacteria <- phyloseq(otu_table(cristiano_rep200Data_Filt_Bacteria, taxa_are_rows = FALSE), 
                             tax_table(as.matrix(rep200TaxSplit_Bacteria_Formatted)), 
                             sample_data(metaCristianoTxNaiveFilt))

## Aggregate counts
psCristianoBacteria_phylum = aggregate_taxa(psCristianoBacteria, "Phylum")
psCristianoBacteria_class = aggregate_taxa(psCristianoBacteria, "Class")
psCristianoBacteria_order = aggregate_taxa(psCristianoBacteria, "Order")
psCristianoBacteria_family = aggregate_taxa(psCristianoBacteria, "Family")
psCristianoBacteria_genus = aggregate_taxa(psCristianoBacteria, "Genus")
psCristianoBacteria_species = aggregate_taxa(psCristianoBacteria, "Species") #cristiano_rep200Data_Filt_Bacteria
# colnames(psCristianoBacteria_species) <- rep200TaxSplit_Bacteria_Formatted[colnames(psCristianoBacteria_species),"Species"]

## Create data.frames of summarized data
cristianoRep200BacteriaPhylum <- data.frame(t(otu_table(psCristianoBacteria_phylum)))
cristianoRep200BacteriaClass <- data.frame(t(otu_table(psCristianoBacteria_class)))
cristianoRep200BacteriaOrder <- data.frame(t(otu_table(psCristianoBacteria_order)))
cristianoRep200BacteriaFamily <- data.frame(t(otu_table(psCristianoBacteria_family)))
cristianoRep200BacteriaGenus <- data.frame(t(otu_table(psCristianoBacteria_genus)))
cristianoRep200BacteriaSpecies <- data.frame(t(otu_table(psCristianoBacteria_species))) #data.frame(psCristianoBacteria_species)
cristianoRep200BacteriaSpecies[1:3,1:3]

## Subset taxa to those shared with the Weizmann data
load("Interim_data/shared_bacterial_features_at_each_taxa_level_15Oct21.RData")

psCristianoBacteria_phylum_shared <- subset_taxa(psCristianoBacteria_phylum, Phylum %in% sharedPhylumBacteria$intersectedTaxa)
psCristianoBacteria_class_shared <- subset_taxa(psCristianoBacteria_class, Class %in% sharedClassBacteria$intersectedTaxa)
psCristianoBacteria_order_shared <- subset_taxa(psCristianoBacteria_order, Order %in% sharedOrderBacteria$intersectedTaxa)
psCristianoBacteria_family_shared <- subset_taxa(psCristianoBacteria_family, Family %in% sharedFamilyBacteria$intersectedTaxa)
psCristianoBacteria_genus_shared <- subset_taxa(psCristianoBacteria_genus, Genus %in% sharedGenusBacteria$intersectedTaxa)
psCristianoBacteria_species_shared <- subset_taxa(psCristianoBacteria_species, Species %in% sharedSpeciesBacteria$intersectedTaxa)
# psCristianoBacteria_species_shared <- psCristianoBacteria_species[,colnames(psCristianoBacteria_species) %in% sharedSpeciesBacteria$intersectedTaxa]

## Create data.frames of summarized data with feature intersection
cristianoRep200BacteriaPhylumShared <- data.frame(t(otu_table(psCristianoBacteria_phylum_shared)))
cristianoRep200BacteriaClassShared <- data.frame(t(otu_table(psCristianoBacteria_class_shared)))
cristianoRep200BacteriaOrderShared <- data.frame(t(otu_table(psCristianoBacteria_order_shared)))
cristianoRep200BacteriaFamilyShared <- data.frame(t(otu_table(psCristianoBacteria_family_shared)))
cristianoRep200BacteriaGenusShared <- data.frame(t(otu_table(psCristianoBacteria_genus_shared)))
cristianoRep200BacteriaSpeciesShared <- data.frame(t(otu_table(psCristianoBacteria_species_shared))) # data.frame(psCristianoBacteria_species_shared)
cristianoRep200BacteriaSpeciesShared[1:3,1:3]

#-----------------------------------------#
# Calculating log ratios based on MMvec findings
#-----------------------------------------#

f1Genera <- c("Malassezia","Trichosporon","Ramularia")
f2Genera <- c("Aspergillus","Candida")
f3Genera <- c("Colletotrichum","Fusarium","Cutaneotrichosporon","Phialocephala","Trichoderma","Talaromyces",
              "Yarrowia","Stereum","Aureobasidium","Hyphopichia","Dissoconium","Agaricus","Exophiala",
              "Alternaria","Tilletiopsis","Cryptococcus","Penicillium","Puccinia")
# Subset with a pseudocount (to avoid dropping many samples)
cristianoRep200FungiGenusSubset <- cristianoRep200FungiGenus[,c(f1Genera,f2Genera,f3Genera)]+1 

lrMetaHopkins <- metaCristianoTxNaiveFilt
lrMetaHopkins$lr_labels <- gsub(" cancer| Cancer","", lrMetaHopkins$phenotype)
lrMetaHopkins$lr_f1_f2 <- log10(rowSums(cristianoRep200FungiGenusSubset[,f1Genera])/rowSums(cristianoRep200FungiGenusSubset[,f2Genera]))
lrMetaHopkins$lr_f1_f3 <- log10(rowSums(cristianoRep200FungiGenusSubset[,f1Genera])/rowSums(cristianoRep200FungiGenusSubset[,f3Genera]))
lrMetaHopkins$lr_f2_f3 <- log10(rowSums(cristianoRep200FungiGenusSubset[,f2Genera])/rowSums(cristianoRep200FungiGenusSubset[,f3Genera]))

#------------Between cancer types------------#
# NOTE: Duodenal Cancer only has 1 sample and is removed below
lrMetaHopkins %>%
  filter(phenotype != "Healthy") %>%
  filter(phenotype != "Duodenal Cancer") %>% 
  filter(is.finite(lr_f1_f2)) %>%
  filter(!is.na(lr_f1_f2)) %>%
  ggplot(aes(reorder(lr_labels, -lr_f1_f2, FUN=median),lr_f1_f2, fill=lr_labels)) +
  geom_boxplot(notch = TRUE) + scale_fill_nejm() + theme_pubr() + geom_point(alpha=0.4) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  labs(x = "Cancer Types", y = "log(F1/F2)", title = "Log ratio of F1/F2 in plasma data") +
  stat_compare_means(method = "anova", label.y = 2.2)
ggsave(filename = "Figures/Figure_5/log_ratio_cristiano_f1_f2.svg", dpi = "retina",
       units = "in", width = 6.5, height = 4)

lrMetaHopkins %>%
  filter(phenotype != "Healthy") %>%
  filter(phenotype != "Duodenal Cancer") %>% 
  filter(is.finite(lr_f1_f3)) %>%
  filter(!is.na(lr_f1_f3)) %>%
  ggplot(aes(reorder(lr_labels, -lr_f1_f3, FUN=median),lr_f1_f3, fill=lr_labels)) +
  geom_boxplot(notch = TRUE) + scale_fill_nejm() + theme_pubr() + geom_point(alpha=0.4) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  labs(x = "Cancer Types", y = "log(F1/F3)", title = "Log ratio of F1/F3 in plasma data") +
  stat_compare_means(method = "anova", label.y = 1.5)
ggsave(filename = "Figures/Figure_5/log_ratio_cristiano_f1_f3.svg", dpi = "retina",
       units = "in", width = 6.5, height = 4)

lrMetaHopkins %>%
  filter(phenotype != "Healthy") %>%
  filter(phenotype != "Duodenal Cancer") %>% 
  filter(is.finite(lr_f2_f3)) %>%
  filter(!is.na(lr_f2_f3)) %>%
  ggplot(aes(reorder(lr_labels, -lr_f2_f3, FUN=mean),lr_f2_f3, fill=lr_labels)) +
  geom_boxplot(notch = TRUE) + scale_fill_nejm() + theme_pubr() + geom_point(alpha=0.4) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  labs(x = "Cancer Types", y = "log(F2/F3)", title = "Log ratio of F2/F3 in plasma data") +
  stat_compare_means(method = "anova", label.y = 0.5)
ggsave(filename = "Figures/Figure_5/log_ratio_cristiano_f2_f3.svg", dpi = "retina",
       units = "in", width = 6.5, height = 4)

#------------Cancer vs Healthy------------#
lrMetaHopkins %>%
  filter(is.finite(lr_f1_f2)) %>%
  filter(!is.na(lr_f1_f2)) %>%
  ggplot(aes(reorder(cancer_status, -lr_f1_f2, FUN=median),lr_f1_f2, fill=cancer_status)) +
  geom_boxplot(notch = TRUE) + scale_fill_nejm() + theme_pubr() + geom_point(alpha=0.4) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  labs(x = "Disease Status", y = "log(F1/F2)", title = "Hopkins: Log ratio of F1/F2") +
  stat_compare_means(method = "wilcox.test", label.y = 2.2, comparisons = list(c("Cancer","Healthy")))
ggsave(filename = "Figures/Figure_5/log_ratio_cristiano_f1_f2_cancer_vs_healthy.svg", dpi = "retina",
       units = "in", width = 3, height = 4)

lrMetaHopkins %>%
  filter(is.finite(lr_f1_f3)) %>%
  filter(!is.na(lr_f1_f3)) %>%
  ggplot(aes(reorder(cancer_status, -lr_f1_f3, FUN=median),lr_f1_f3, fill=cancer_status)) +
  geom_boxplot(notch = TRUE) + scale_fill_nejm() + theme_pubr() + geom_point(alpha=0.4) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  labs(x = "Disease Status", y = "log(F1/F3)", title = "Hopkins: Log ratio of F1/F3") +
  stat_compare_means(method = "wilcox.test", label.y = 2.2, comparisons = list(c("Cancer","Healthy")))
ggsave(filename = "Figures/Figure_5/log_ratio_cristiano_f1_f3_cancer_vs_healthy.svg", dpi = "retina",
       units = "in", width = 3, height = 4)

lrMetaHopkins %>%
  filter(is.finite(lr_f2_f3)) %>%
  filter(!is.na(lr_f2_f3)) %>%
  ggplot(aes(reorder(cancer_status, -lr_f2_f3, FUN=median),lr_f2_f3, fill=cancer_status)) +
  geom_boxplot(notch = TRUE) + scale_fill_nejm() + theme_pubr() + geom_point(alpha=0.4) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  labs(x = "Disease Status", y = "log(F2/F3)", title = "Hopkins: Log ratio of F2/F3") +
  stat_compare_means(method = "wilcox.test", label.y = 1.2, comparisons = list(c("Cancer","Healthy")))
ggsave(filename = "Figures/Figure_5/log_ratio_cristiano_f2_f3_cancer_vs_healthy.svg", dpi = "retina",
       units = "in", width = 3, height = 4)

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
                                                  varImpFlag = TRUE)
cristianoCancerVsHealthy_FungiDecontam <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                      countData = cristiano_rep200Data_Filt_Fungi_Decontam,
                                                      dataString = "fungi_decontam",
                                                      col2Predict = "cancer_status",
                                                      varImpFlag = TRUE)
cristianoCancerVsHealthy_FungiIntersect <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                       countData = cristianoRep200FungiSpeciesShared,
                                                       dataString = "fungi_intersected_with_Weizmann",
                                                       col2Predict = "cancer_status",
                                                       varImpFlag = FALSE)
cristianoCancerVsHealthy_BacteriaIntersect <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                       countData = cristianoRep200BacteriaSpeciesShared,
                                                       dataString = "bacteria_intersected_with_Weizmann",
                                                       col2Predict = "cancer_status",
                                                       varImpFlag = TRUE)
cristianoCancerVsHealthy_FungiBacteriaIntersect <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                                          countData = cbind(cristianoRep200FungiSpeciesShared, 
                                                                            cristianoRep200BacteriaSpeciesShared),
                                                          dataString = "fungi_and_bacteria_intersected_with_Weizmann",
                                                          col2Predict = "cancer_status",
                                                          varImpFlag = FALSE)
# Compile results
cristianoCancerVsHealthyResults <- rbind(cristianoCancerVsHealthy_FullRep200$rep_perf,
                                         cristianoCancerVsHealthy_BacteriaOnly$rep_perf,
                                         cristianoCancerVsHealthy_FungiOnly$rep_perf,
                                         cristianoCancerVsHealthy_FungiDecontam$rep_perf,
                                         cristianoCancerVsHealthy_FungiIntersect$rep_perf,
                                         cristianoCancerVsHealthy_BacteriaIntersect$rep_perf,
                                         cristianoCancerVsHealthy_FungiBacteriaIntersect$rep_perf)

colnames(cristianoCancerVsHealthyResults)[1:2] <- c("AUROC","AUPR")
cristianoCancerVsHealthyResults$nullAUROC <- 0.5
# In the line below, "minorityClassSize" isn't technically the minority but rather
# the positive class (cancer), so it still serves to infer the null AUPR
cristianoCancerVsHealthyResults$nullAUPR <- cristianoCancerVsHealthyResults$minorityClassSize/
  (cristianoCancerVsHealthyResults$minorityClassSize+cristianoCancerVsHealthyResults$majorityClassSize)

cristianoCancerVsHealthyResults$dataString[cristianoCancerVsHealthyResults$dataString=="full_rep200"] <- "Full multikingdom database (rep200)"
cristianoCancerVsHealthyResults$dataString[cristianoCancerVsHealthyResults$dataString=="bacteria_only"] <- "Bacteria all (Species)"
cristianoCancerVsHealthyResults$dataString[cristianoCancerVsHealthyResults$dataString=="fungi_only"] <- "Fungi all (Species)"
cristianoCancerVsHealthyResults$dataString[cristianoCancerVsHealthyResults$dataString=="fungi_decontam"] <- "Fungi decontaminated (Species)"
cristianoCancerVsHealthyResults$dataString[cristianoCancerVsHealthyResults$dataString=="fungi_intersected_with_Weizmann"] <- "Fungi ∩ WIS (Species)"
cristianoCancerVsHealthyResults$dataString[cristianoCancerVsHealthyResults$dataString=="bacteria_intersected_with_Weizmann"] <- "Bacteria ∩ WIS (Species)"
cristianoCancerVsHealthyResults$dataString[cristianoCancerVsHealthyResults$dataString=="fungi_and_bacteria_intersected_with_Weizmann"] <- "Fungi+bacteria ∩ WIS (Species)"
cristianoCancerVsHealthyResults$dataString <- factor(cristianoCancerVsHealthyResults$dataString, levels = c("Full multikingdom database (rep200)",
                                                                                                            "Fungi+bacteria ∩ WIS (Species)",
                                                                                                            "Bacteria all (Species)",
                                                                                                            "Fungi all (Species)",
                                                                                                            "Fungi decontaminated (Species)",
                                                                                                            "Fungi ∩ WIS (Species)",
                                                                                                            "Bacteria ∩ WIS (Species)"))
require(ggrepel)
source("Supporting_scripts/S05-SummarySE.R")
# Facet by AUROC vs AUPR; color by full dataset vs fungi (dataString)
cristianoCancerVsHealthyResults %>%
  filter(grepl("Full|\\+|decontam",dataString)) %>%
  reshape2::melt(id.vars = c("dataString","rep","diseaseType", "col2Predict","minorityClassSize","majorityClassSize","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("diseaseType","dataString","variable","minorityClassSize","majorityClassSize","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=dataString)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + 
  ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + 
  ggtitle("Hopkins plasma validation cohort:\nCancer vs. healthy") +
  theme(plot.title = element_text(hjust = 0), legend.position = "right") + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "Cancer vs. Healthy") +
  rotate_x_text(0) + scale_color_aaas(name = "Dataset") + geom_hline(yintercept = 1, linetype="dashed") #+ 
  # geom_text_repel(aes(label=round(value,2)), size=3, box.padding = 0.75, force = 100, 
  #                  point.padding = 6.5, direction = "y", position = position_dodge(width = 1), show.legend = FALSE) #+
ggsave("Figures/Figure_5/cristiano_all_cancer_vs_healthy_full_WIS_fungi_15Nov21.svg", 
       dpi = "retina", width = 6, height = 3, units = "in")
cristianoCancerVsHealthyResults %>% write.csv("Figures_data/Figure_5/cristiano_all_cancer_vs_healthy_full_WIS_fungi_15Nov21.csv")

#---------------------------------------------------------------------------------------------------------------------------#
# ML repeated CV perf --> ROC plot with 99% CIs
#---------------------------------------------------------------------------------------------------------------------------#

mlFullDB_CV <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                           countData = cristiano_rep200Data_Filt,
                           dataString = "full_rep200",
                           col2Predict = "cancer_status",
                           numResampleIter = 10,
                           varImpFlag = FALSE)
mlFungiBacteriaWISIntersect_CV <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                             countData = cbind(cristianoRep200FungiSpeciesShared, 
                                                               cristianoRep200BacteriaSpeciesShared),
                                             dataString = "fungi_and_bacteria_intersected_with_Weizmann",
                                             col2Predict = "cancer_status",
                                             numResampleIter = 10,
                                             varImpFlag = FALSE)
save(mlFungiBacteriaWISIntersect_CV, file = "Interim_data/hopkins_mlFungiBacteriaWISIntersect_CV_15Nov21.RData")

mlBacteriaWISIntersect_CV <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                              countData = cristianoRep200BacteriaSpeciesShared,
                                              dataString = "bacteria_intersected_with_Weizmann",
                                              col2Predict = "cancer_status",
                                              numResampleIter = 10,
                                              varImpFlag = FALSE)

mlFungiWISIntersect_CV <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                         countData = cristianoRep200FungiSpeciesShared,
                                         dataString = "fungi_intersected_with_Weizmann",
                                         col2Predict = "cancer_status",
                                         numResampleIter = 10,
                                         varImpFlag = FALSE)

mlFungiDecontam_CV <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                               countData = cristiano_rep200Data_Filt_Fungi_Decontam,
                               dataString = "fungi_intersected_with_Weizmann",
                               col2Predict = "cancer_status",
                               numResampleIter = 10,
                               varImpFlag = FALSE)
#----------Identify top X signature----------#
load("Interim_data/decontamResultsV2_13Oct21.RData")
topXNum <- 20
topXHopkinsSig <- decontamResultsV2[head(mlFungiDecontam_CV$varImp$Taxa,topXNum),c("species","reason")]
save(topXHopkinsSig, file = "Interim_data/topXHopkinsSig_13Nov21.RData")
topXHopkinsSig %>% write.csv(file = "Figures/Figure_5/topXHopkinsSig.csv")

topXHopkinsData <- cristiano_rep200Data_Filt_Fungi_Decontam[,rownames(topXHopkinsSig)]

mlFungiDecontam_CV_topFeat <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                    countData = topXHopkinsData,
                    dataString = "fungi_intersected_with_Weizmann",
                    col2Predict = "cancer_status",
                    numResampleIter = 10,
                    varImpFlag = FALSE)
hopkinsFungi_topX <- mlFungiDecontam_CV_topFeat
save(hopkinsFungi_topX, file = "Interim_data/cristianoFungi_topX_ML_model_15Nov21.RData")
#--------------------------------------------#

plots_mlFullDB_CV <- plotMLWithCIs(mlFullDB_CV, 
                                   showRepCurves=TRUE,
                                   sizeAvgCurve=0.5, 
                                   sizeRepCurves = 0.25,
                                   ciAlpha = 0.4,
                                   colorAvgCurve = "#3B4992FF",
                                   ciFillColor = "blue",
                                   ciLevel = 0.99)
plots_mlFungiBacteriaWISIntersect_CV <- plotMLWithCIs(mlFungiBacteriaWISIntersect_CV, 
                                                      showRepCurves=TRUE, 
                                                      sizeAvgCurve=0.5, 
                                                      sizeRepCurves = 0.25,
                                                      ciAlpha = 0.4, 
                                                      colorAvgCurve = "#BB0021FF",
                                                      ciFillColor = "red",
                                                      ciLevel = 0.99)
plots_mlBacteriaWISIntersect_CV <- plotMLWithCIs(mlBacteriaWISIntersect_CV, 
                                                      showRepCurves=TRUE, 
                                                      sizeAvgCurve=0.5, 
                                                      sizeRepCurves = 0.25,
                                                      ciAlpha = 0.4, 
                                                      colorAvgCurve = "#631879",
                                                      ciFillColor = "purple",
                                                      ciLevel = 0.99)
plots_mlFungiWISIntersect_CV <- plotMLWithCIs(mlFungiWISIntersect_CV, 
                                                 showRepCurves=TRUE, 
                                                 sizeAvgCurve=0.5, 
                                                 sizeRepCurves = 0.25,
                                                 ciAlpha = 0.4, 
                                                 colorAvgCurve = "#631879",
                                                 ciFillColor = "purple",
                                                 ciLevel = 0.99)
plots_mlFungiDecontam_CV <- plotMLWithCIs(mlFungiDecontam_CV,
                                          sizeAvgCurve=0.5, 
                                          showRepCurves=TRUE, 
                                          ciAlpha = 0.4, 
                                          sizeRepCurves = 0.25,
                                          colorAvgCurve = "#008B45FF",
                                          ciFillColor = "lightgreen",
                                          ciLevel = 0.99)

plots_mlFungiDecontam_CV_topFeat <- plotMLWithCIs(mlFungiDecontam_CV_topFeat,
                                          sizeAvgCurve=0.5, 
                                          showRepCurves=TRUE, 
                                          ciAlpha = 0.4, 
                                          sizeRepCurves = 0.25,
                                          colorAvgCurve = "#631879",
                                          ciFillColor = "purple",
                                          ciLevel = 0.99)

# save(mlFullDB_CV,
#      mlFungiDecontam_CV,
#      mlFungiBacteriaWISIntersect_CV,
#      plots_mlFullDB_CV,
#      plots_mlFungiBacteriaWISIntersect_CV,
#      plots_mlFungiDecontam_CV,
#      file = "Interim_data/cristiano_CV10_overlay_data_and_plots_11Nov21.RData")

# Combine ROC curves
rocCombined <- plots_mlFungiBacteriaWISIntersect_CV$rocPlot
for(ii in 1:10){
  rocCombined <- rocCombined + 
    geom_path(data = plots_mlFullDB_CV$rocCurveData[[ii]],aes(x=fpr,y=tpr), color = "lightgray", size = 0.5) +
    geom_path(data = plots_mlFungiDecontam_CV$rocCurveData[[ii]],aes(x=fpr,y=tpr), color = "lightgray", size = 0.5)
}
rocCombinedAll <- rocCombined +
  geom_ribbon(data = plots_mlFullDB_CV$interpROCYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
              fill = "lightblue", alpha = 0.6, inherit.aes = F) +
  geom_path(data = plots_mlFullDB_CV$interpROCYDf_CI, aes(x = xval, y = Estimate), color = "#3B4992FF", size = 0.5) +
  geom_ribbon(data = plots_mlFungiDecontam_CV$interpROCYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
              fill = "lightgreen", alpha = 0.3, inherit.aes = F) + theme_bw() +
  geom_path(data = plots_mlFungiDecontam_CV$interpROCYDf_CI, aes(x = xval, y = Estimate), color = "#008B45FF", size = 0.5) +
  annotate("segment", x = 0.35, xend = 0.45, y = 0.2, yend = 0.2, color = "#3B4992FF") +
  annotate("text", x = 0.47, y = 0.2, color = "#3B4992FF", label = paste0("AUROC 99% CI: [", 
                                                     paste0(100*round(plots_mlFullDB_CV$aurocCI[2],4),
                                                            ", ",100*round(plots_mlFullDB_CV$aurocCI[3],4)),"]"), hjust = 0) +
  annotate("segment", x = 0.35, xend = 0.45, y = 0.15, yend = 0.15, color = "#BB0021FF") +
  annotate("text", x = 0.47, y = 0.15, color = "#BB0021FF", label = paste0("AUROC 99% CI: [", 
                                                                           paste0(100*round(plots_mlFungiBacteriaWISIntersect_CV$aurocCI[2],4),
                                                                                  ", ",100*round(plots_mlFungiBacteriaWISIntersect_CV$aurocCI[3],4)),"]"), hjust = 0) +
  annotate("segment", x = 0.35, xend = 0.45, y = 0.1, yend = 0.1, color = "#008B45FF") +
  annotate("text", x = 0.47, y = 0.1, color = "#008B45FF", label = paste0("AUROC 99% CI: [", 
                                                                          paste0(100*round(plots_mlFungiDecontam_CV$aurocCI[2],4),
                                                                                 ", ",100*round(plots_mlFungiDecontam_CV$aurocCI[3],4)),"]"), hjust = 0)

# Combine PR curves
prCombined <- plots_mlFungiBacteriaWISIntersect_CV$prPlot
for(ii in 1:10){
  prCombined <- prCombined + 
    geom_path(data = plots_mlFullDB_CV$prCurveData[[ii]],aes(x=recall,y=precision), color = "lightgray", size = 0.5) +
    geom_path(data = plots_mlFungiDecontam_CV$prCurveData[[ii]],aes(x=recall,y=precision), color = "lightgray", size = 0.5)
}
prCombinedAll <- prCombined +
  geom_ribbon(data = plots_mlFullDB_CV$interpPRYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
              fill = "lightblue", alpha = 0.6, inherit.aes = F) +
  geom_path(data = plots_mlFullDB_CV$interpPRYDf_CI, aes(x = xval, y = Estimate), color = "#3B4992FF", size = 0.5) +
  geom_ribbon(data = plots_mlFungiDecontam_CV$interpPRYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
              fill = "lightgreen", alpha = 0.3, inherit.aes = F) + theme_bw() +
  geom_path(data = plots_mlFungiDecontam_CV$interpPRYDf_CI, aes(x = xval, y = Estimate), color = "#008B45FF", size = 0.5) +
  annotate("segment", x = 0.15, xend = 0.25, y = 0.2, yend = 0.2, color = "#3B4992FF") +
  annotate("text", x = 0.27, y = 0.2, color = "#3B4992FF", label = paste0("AUPR 99% CI: [", 
                                                                          paste0(100*round(plots_mlFullDB_CV$auprCI[2],4),
                                                                                 ", ",100*round(plots_mlFullDB_CV$auprCI[3],4)),"]"), hjust = 0) +
  annotate("segment", x = 0.15, xend = 0.25, y = 0.15, yend = 0.15, color = "#BB0021FF") +
  annotate("text", x = 0.27, y = 0.15, color = "#BB0021FF", label = paste0("AUPR 99% CI: [", 
                                                                           paste0(100*round(plots_mlFungiBacteriaWISIntersect_CV$auprCI[2],4),
                                                                                  ", ",100*round(plots_mlFungiBacteriaWISIntersect_CV$auprCI[3],4)),"]"), hjust = 0) +
  annotate("segment", x = 0.15, xend = 0.25, y = 0.1, yend = 0.1, color = "#008B45FF") +
  annotate("text", x = 0.27, y = 0.1, color = "#008B45FF", label = paste0("AUPR 99% CI: [", 
                                                                          paste0(100*round(plots_mlFungiDecontam_CV$auprCI[2],4),
                                                                                 ", ",100*round(plots_mlFungiDecontam_CV$auprCI[3],4)),"]"), hjust = 0)
combinedPlotTitle <- paste0("Hopkins plasma cohort: Cancer vs Healthy (10-fold CV repeated 10 times)")
combinedOverlayPlots <- ggarrange(rocCombinedAll, prCombinedAll, ncol = 2)
combinedOverlayPlotsAnnotated <- annotate_figure(combinedOverlayPlots, 
                                                      top = text_grob(combinedPlotTitle, 
                                                                      color = "black", face = "bold", size = 14))
print(combinedOverlayPlotsAnnotated)
ggsave(filename = "Figures/Figure_5/roc_pr_CV_overlay_cristiano_full_WIS_fungionly.svg",
       units = "in", width = 10, height = 5)


#--------------------Overlay topX model with decontaminated model--------------------#
# Combine ROC curves
rocCombinedTopX <- plots_mlFungiDecontam_CV_topFeat$rocPlot
for(ii in 1:10){
  rocCombinedTopX <- rocCombinedTopX + 
    geom_path(data = plots_mlFungiDecontam_CV$rocCurveData[[ii]],aes(x=fpr,y=tpr), color = "lightgray", size = 0.25)
}
rocCombinedTopXAll <- rocCombinedTopX +
  geom_ribbon(data = plots_mlFungiDecontam_CV$interpROCYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
              fill = "lightgreen", alpha = 0.3, inherit.aes = F) + theme_bw() +
  geom_path(data = plots_mlFungiDecontam_CV$interpROCYDf_CI, aes(x = xval, y = Estimate), color = "#008B45FF", size = 0.5) +
  annotate("segment", x = 0.25, xend = 0.35, y = 0.15, yend = 0.15, color = "#631879") +
  annotate("text", x = 0.37, y = 0.15, color = "#631879", label = paste0("AUROC 99% CI: [", 
                                                                           paste0(100*round(plots_mlFungiDecontam_CV_topFeat$aurocCI[2],4),
                                                                                  ", ",100*round(plots_mlFungiDecontam_CV_topFeat$aurocCI[3],4)),"]"), hjust = 0) +
  annotate("segment", x = 0.25, xend = 0.35, y = 0.1, yend = 0.1, color = "#008B45FF") +
  annotate("text", x = 0.37, y = 0.1, color = "#008B45FF", label = paste0("AUROC 99% CI: [", 
                                                                          paste0(100*round(plots_mlFungiDecontam_CV$aurocCI[2],4),
                                                                                 ", ",100*round(plots_mlFungiDecontam_CV$aurocCI[3],4)),"]"), hjust = 0)

# Combine PR curves
prCombinedTopX <- plots_mlFungiDecontam_CV_topFeat$prPlot
for(ii in 1:10){
  prCombinedTopX <- prCombinedTopX + 
    geom_path(data = plots_mlFungiDecontam_CV$prCurveData[[ii]],aes(x=recall,y=precision), color = "lightgray", size = 0.25)
}
prCombinedTopXAll <- prCombinedTopX +
  geom_ribbon(data = plots_mlFungiDecontam_CV$interpPRYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
              fill = "lightgreen", alpha = 0.3, inherit.aes = F) + theme_bw() +
  geom_path(data = plots_mlFungiDecontam_CV$interpPRYDf_CI, aes(x = xval, y = Estimate), color = "#008B45FF", size = 0.5) +
  annotate("segment", x = 0.1, xend = 0.2, y = 0.15, yend = 0.15, color = "#631879") +
  annotate("text", x = 0.22, y = 0.15, color = "#631879", label = paste0("AUPR 99% CI: [", 
                                                                           paste0(100*round(plots_mlFungiDecontam_CV_topFeat$auprCI[2],4),
                                                                                  ", ",100*round(plots_mlFungiDecontam_CV_topFeat$auprCI[3],4)),"]"), hjust = 0) +
  annotate("segment", x = 0.1, xend = 0.2, y = 0.1, yend = 0.1, color = "#008B45FF") +
  annotate("text", x = 0.22, y = 0.1, color = "#008B45FF", label = paste0("AUPR 99% CI: [", 
                                                                          paste0(100*round(plots_mlFungiDecontam_CV$auprCI[2],4),
                                                                                 ", ",100*round(plots_mlFungiDecontam_CV$auprCI[3],4)),"]"), hjust = 0)
combinedPlotTitleTopX <- paste0("Hopkins plasma cohort: Cancer vs Healthy (10-fold CV repeated 10 times)\n",
                                "Top ",topXNum," fungi used (purple) vs. decontaminated dataset with 209 fungi (green)")
combinedOverlayPlotsTopX <- ggarrange(rocCombinedTopXAll, prCombinedTopXAll, ncol = 2)
combinedOverlayPlotsAnnotatedTopX <- annotate_figure(combinedOverlayPlotsTopX, 
                                                 top = text_grob(combinedPlotTitleTopX, 
                                                                 color = "black", face = "bold", size = 14))
print(combinedOverlayPlotsAnnotatedTopX)
ggsave(filename = "Figures/Figure_5/roc_pr_CV_overlay_cristiano_fungi_decontam_vs_topX.svg",
       units = "in", width = 10, height = 5)

#--------------------Overlay WIS fungi, bacteria, fungi+bacteria--------------------#

synergyWISDf <- data.frame(AUROC = c(plots_mlFungiWISIntersect_CV$auroc,
                                     plots_mlBacteriaWISIntersect_CV$auroc,
                                     plots_mlFungiBacteriaWISIntersect_CV$auroc),
                           AUPR = c(plots_mlFungiWISIntersect_CV$aupr,
                                     plots_mlBacteriaWISIntersect_CV$aupr,
                                     plots_mlFungiBacteriaWISIntersect_CV$aupr),
                           Dataset = c(rep("Fungi",10),
                                       rep("Bacteria",10),
                                       rep("Fungi+Bacteria",10)))

source("Supporting_scripts/S05-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals
synergyWISDf %>%
  reshape2::melt(id.vars = c("Dataset")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","Dataset")) %>%
  ggplot(aes(reorder(Dataset, value, FUN=median),value, color=Dataset)) +
  facet_wrap(~variable) + 
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),
                width=0.2,size=0.6,position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=0.5) + 
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() + scale_color_nejm()

require(rstatix)
synergyWISDf %>%
  wilcox_test(AUROC ~ Dataset) %>%
  adjust_pvalue() %>%
  add_significance() %>%
  add_xy_position(group = "Dataset") -> roc.stat

synergyWISDf %>%
  ggerrorplot(x = "Dataset",
              y = "AUROC",
              add = "jitter",
              color = "Dataset",
              legend = "none",
              xlab = "WIS ∩ feature set",
              size = 0.1,
              ci = 0.99,
              palette = "nejm",
              add.params = list(alpha=0.4)) +
  stat_pvalue_manual(data = roc.stat,
                     label = "Wilcoxon, q = {p.adj}") -> roc.plot
synergyWISDf %>%
  filter(!grepl("Fungi$",Dataset)) %>%
  ggerrorplot(x = "Dataset",
              y = "AUROC",
              add = "jitter",
              color = "Dataset",
              legend = "none",
              xlab = "",
              ylab = "",
              ylim = c(0.945, 0.968),
              size = 0.5,
              ci = 0.99,
              palette = c("#0072B5FF","#E18727FF"),
              add.params = list(alpha=0.4)) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("Bacteria","Fungi+Bacteria")),
                     label = "p.signif") -> roc.plot.inset
print(roc.plot.inset)
ggsave(filename = "Figures/Supplementary_Figures/hopkins_synergy_CvsH_plots_roc_inset.svg",
       dpi = "retina", units = "in", width = 3, height = 4)

synergyWISDf %>%
  wilcox_test(AUPR ~ Dataset) %>%
  adjust_pvalue() %>%
  add_significance() %>%
  add_xy_position(group = "Dataset") -> pr.stat

synergyWISDf %>%
  ggerrorplot(x = "Dataset",
              y = "AUPR",
              add = "jitter",
              color = "Dataset",
              legend = "none",
              xlab = "WIS ∩ feature set",
              size = 0.1,
              ci = 0.99,
              palette = "nejm",
              add.params = list(alpha=0.4)) +
  stat_pvalue_manual(data = pr.stat,
                     label = "Wilcoxon, q = {p.adj}") -> pr.plot
synergyWISDf %>%
  filter(!grepl("Fungi$",Dataset)) %>%
  ggerrorplot(x = "Dataset",
              y = "AUPR",
              add = "jitter",
              color = "Dataset",
              legend = "none",
              xlab = "",
              ylab = "",
              size = 0.5,
              ci = 0.99,
              ylim = c(0.955, 0.970),
              palette = c("#0072B5FF","#E18727FF"),
              add.params = list(alpha=0.4)) +
  stat_compare_means(method = "wilcox", 
                     comparisons = list(c("Bacteria","Fungi+Bacteria")),
                     label = "p.signif") -> pr.plot.inset
print(pr.plot.inset)
ggsave(filename = "Figures/Supplementary_Figures/hopkins_synergy_CvsH_plots_pr_inset.svg",
       dpi = "retina", units = "in", width = 3, height = 4)

combinedSynergyWISPlot <- ggarrange(roc.plot, pr.plot, ncol = 2)
print(combinedSynergyWISPlot)
ggsave(filename = "Figures/Supplementary_Figures/hopkins_synergy_CvsH_plots.svg",
       dpi = "retina", units = "in", width = 8, height = 6)

  
roc.stat.test <- cristianoStageIterate_results %>%
  filter(grepl("Full|\\+|decontam|Top",dataString)) %>%
  group_by(dataString) %>%
  anova_test(AUROC ~ stageNum) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  data.frame() %>%
  mutate(p.adj = round(p.adj,4)) %>%
  mutate(group1="stageNum", group2 = "stageNum") %>%
  mutate(xmin = 2)
# Combine ROC curves
rocCombinedSynergy <- plots_mlBacteriaWISIntersect_CV$rocPlot
for(ii in 1:10){
  rocCombinedSynergy <- rocCombinedTopX + 
    geom_path(data = plots_mlFungiBacteriaWISIntersect_CV$rocCurveData[[ii]],aes(x=fpr,y=tpr), color = "lightgray", size = 0.25)
}
rocCombinedSynergyAll <- rocCombinedSynergy +
  geom_ribbon(data = plots_mlFungiBacteriaWISIntersect_CV$interpROCYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
              fill = "lightgreen", alpha = 0.3, inherit.aes = F) + theme_bw() +
  geom_path(data = plots_mlFungiBacteriaWISIntersect_CV$interpROCYDf_CI, aes(x = xval, y = Estimate), color = "#008B45FF", size = 0.5) +
  annotate("segment", x = 0.25, xend = 0.35, y = 0.15, yend = 0.15, color = "#631879") +
  annotate("text", x = 0.37, y = 0.15, color = "#631879", label = paste0("AUROC 99% CI: [", 
                                                                         paste0(100*round(plots_mlBacteriaWISIntersect_CV$aurocCI[2],4),
                                                                                ", ",100*round(plots_mlBacteriaWISIntersect_CV$aurocCI[3],4)),"]"), hjust = 0) +
  annotate("segment", x = 0.25, xend = 0.35, y = 0.1, yend = 0.1, color = "#008B45FF") +
  annotate("text", x = 0.37, y = 0.1, color = "#008B45FF", label = paste0("AUROC 99% CI: [", 
                                                                          paste0(100*round(plots_mlFungiBacteriaWISIntersect_CV$aurocCI[2],4),
                                                                                 ", ",100*round(plots_mlFungiBacteriaWISIntersect_CV$aurocCI[3],4)),"]"), hjust = 0)

# Combine PR curves
prCombinedSynergy <- plots_mlBacteriaWISIntersect_CV$prPlot
for(ii in 1:10){
  prCombinedSynergy <- prCombinedSynergy + 
    geom_path(data = plots_mlFungiBacteriaWISIntersect_CV$prCurveData[[ii]],aes(x=recall,y=precision), color = "lightgray", size = 0.25)
}
prCombinedSynergyAll <- prCombinedSynergy +
  geom_ribbon(data = plots_mlFungiBacteriaWISIntersect_CV$interpPRYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
              fill = "lightgreen", alpha = 0.3, inherit.aes = F) + theme_bw() +
  geom_path(data = plots_mlFungiBacteriaWISIntersect_CV$interpPRYDf_CI, aes(x = xval, y = Estimate), color = "#008B45FF", size = 0.5) +
  annotate("segment", x = 0.1, xend = 0.2, y = 0.15, yend = 0.15, color = "#631879") +
  annotate("text", x = 0.22, y = 0.15, color = "#631879", label = paste0("AUPR 99% CI: [", 
                                                                         paste0(100*round(plots_mlBacteriaWISIntersect_CV$auprCI[2],4),
                                                                                ", ",100*round(plots_mlBacteriaWISIntersect_CV$auprCI[3],4)),"]"), hjust = 0) +
  annotate("segment", x = 0.1, xend = 0.2, y = 0.1, yend = 0.1, color = "#008B45FF") +
  annotate("text", x = 0.22, y = 0.1, color = "#008B45FF", label = paste0("AUPR 99% CI: [", 
                                                                          paste0(100*round(plots_mlFungiBacteriaWISIntersect_CV$auprCI[2],4),
                                                                                 ", ",100*round(plots_mlFungiBacteriaWISIntersect_CV$auprCI[3],4)),"]"), hjust = 0)
combinedPlotTitleSynergy <- paste0("Hopkins plasma cohort: Cancer vs Healthy (10-fold CV repeated 10 times)\n",
                                "WIS synergistic performance\n",
                                "WIS fungi (purple) vs. WIS bacteria vs. WIS fungi+bacteria (green)")
combinedOverlayPlotsSynergy <- ggarrange(rocCombinedSynergyAll, prCombinedSynergyAll, ncol = 2)
combinedOverlayPlotsSynergyAnnotated <- annotate_figure(combinedOverlayPlotsSynergy, 
                                                     top = text_grob(combinedPlotTitleSynergy, 
                                                                     color = "black", face = "bold", size = 14))
print(combinedOverlayPlotsSynergyAnnotated)
ggsave(filename = "Figures/Figure_5/roc_pr_CV_overlay_cristiano_WIS_synergy.svg",
       units = "in", width = 10, height = 5)

#---------------------------------------------------------------------------------------------------------------------------#
# ML perf each cancer vs. healthy -- all datasets
#---------------------------------------------------------------------------------------------------------------------------#
source("00-Functions.R") # for the ml1VsAllCristiano10kRep1_Iterate() function
cristianoPerCancerIterate <- ml1VsAllCristiano10kRep1_Iterate()

cristianoPerCancerIterate_results <- cristianoPerCancerIterate$rep_perf
cristianoPerCancerIterate_results$nullAUROC <- 0.5
cristianoPerCancerIterate_results$nullAUPR <- cristianoPerCancerIterate_results$minorityClassSize/
  (cristianoPerCancerIterate_results$minorityClassSize+cristianoPerCancerIterate_results$majorityClassSize)

colnames(cristianoPerCancerIterate_results)[1:2] <- c("AUROC","AUPR")
cristianoPerCancerIterate_results$dataString[cristianoPerCancerIterate_results$dataString=="full_rep200"] <- "Full multikingdom database (rep200)"
cristianoPerCancerIterate_results$dataString[cristianoPerCancerIterate_results$dataString=="bacteria_only"] <- "Bacteria all (Species)"
cristianoPerCancerIterate_results$dataString[cristianoPerCancerIterate_results$dataString=="fungi_only"] <- "Fungi all (Species)"
cristianoPerCancerIterate_results$dataString[cristianoPerCancerIterate_results$dataString=="fungi_decontam"] <- "Fungi decontaminated (Species)"
cristianoPerCancerIterate_results$dataString[cristianoPerCancerIterate_results$dataString=="fungi_intersected_with_Weizmann"] <- "Fungi ∩ WIS (Species)"
cristianoPerCancerIterate_results$dataString[cristianoPerCancerIterate_results$dataString=="bacteria_intersected_with_Weizmann"] <- "Bacteria ∩ WIS (Species)"
cristianoPerCancerIterate_results$dataString[cristianoPerCancerIterate_results$dataString=="fungi_and_bacteria_intersected_with_Weizmann"] <- "Fungi+bacteria ∩ WIS (Species)"
cristianoPerCancerIterate_results$dataString[cristianoPerCancerIterate_results$dataString=="topX_fungi"] <- "Top 20 fungi (Species)"
cristianoPerCancerIterate_results$dataString <- factor(cristianoPerCancerIterate_results$dataString, levels = c("Full multikingdom database (rep200)",
                                                                                                                "Fungi+bacteria ∩ WIS (Species)",
                                                                                                                "Bacteria all (Species)",
                                                                                                                "Fungi all (Species)",
                                                                                                                "Fungi decontaminated (Species)",
                                                                                                                "Fungi ∩ WIS (Species)",
                                                                                                                "Bacteria ∩ WIS (Species)",
                                                                                                                "Top 20 fungi (Species)"))
cristianoPerCancerIterate_results$diseaseType <- gsub(" cancer| Cancer","",cristianoPerCancerIterate_results$diseaseType)

save(cristianoPerCancerIterate_results,
     file = "Interim_data/cristianoPerCancerIterate_results_13Nov21.RData")

## Overlay multiple data types on the same graph
source("Supporting_scripts/S05-SummarySE.R")
cristianoPerCancerIterate_results %>%
  filter(grepl("Full|\\+|decontam|Top",dataString)) %>%
  # filter(grepl("\\+|Bacteria ∩",dataString)) %>%
  reshape2::melt(id.vars = c("rep","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=dataString)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid", position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Hopkins plasma validation cohort:\nPer cancer type vs. healthy performance") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Figure_5/figure_5_cristiano_per_cancer_type_vs_healthy_full_WIS_fungi_topX_datasets_13Nov21.svg", 
       dpi = "retina", width = 8, height = 4.5, units = "in")

cristianoPerCancerIterate_results %>% 
  write.csv("Figures_data/Figure_5/figure_5_G_cristiano_per_cancer_type_vs_healthy_all_datasets_13Nov21.csv")

# #----------Testing WIS fungi+bacteria vs. WIS bacteria
# cristianoPerCancerIterate_results %>%
#   filter(grepl("\\+|Bacteria ∩",dataString)) %>%
#   filter(diseaseType %in% c("Breast","Gastric")) %>%
#   ggboxplot(x = "dataString",
#             y = "AUPR",
#             fill = "dataString",
#             notch = TRUE) +
#   stat_compare_means(method = "wilcox", label.y = 1.1)

## Subset to decontaminated fungi
cristianoPerCancerIterate_results %>%
  filter(grepl("decontam",dataString)) %>%
  reshape2::melt(id.vars = c("rep","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  # guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Hopkins plasma validation cohort:\nPer cancer type vs. healthy performance using fungi") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Figure_5/figure_5_G_cristiano_per_cancer_type_vs_healthy_fungi_decontam_only_13Nov21.svg", 
       dpi = "retina", width = 7, height = 4.5, units = "in")

cristianoPerCancerIterate_results %>%
  filter(grepl("decontam",dataString)) %>% 
  pull(AUROC) %>% ci()
# Estimate   CI lower   CI upper Std. Error 
# 0.77668498 0.74385524 0.80951472 0.01645646 

cristianoPerCancerIterate_results %>%
  filter(grepl("decontam",dataString)) %>% 
  filter(grepl("Breast",diseaseType)) %>% 
  pull(AUROC) %>% ci()
# Estimate   CI lower   CI upper Std. Error 
# 0.8746154  0.8139667  0.9352640  0.0268101 

cristianoPerCancerIterate_results %>%
  filter(grepl("Top",dataString)) %>% 
  pull(AUROC) %>% ci()
# Estimate   CI lower   CI upper Std. Error 
# 0.78840659 0.75768774 0.81912545 0.01539834

## Subset to decontaminated and topX fungi
cristianoPerCancerIterate_results %>%
  filter(grepl("decontam|Top",dataString)) %>%
  reshape2::melt(id.vars = c("rep","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=dataString)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  # guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Hopkins plasma validation cohort:\nPer cancer type vs. healthy performance using fungi") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) +
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Figure_5/hopkins_per_cancer_type_vs_healthy_fungi_decontam_and_TopX_13Nov21.svg", 
       dpi = "retina", width = 7, height = 4.5, units = "in")

#---------------------------------------------------------------------------------------------------------------------------#
# ML perf between cancer types -- all datasets
#---------------------------------------------------------------------------------------------------------------------------#
source("00-Functions.R") # for the ml1VsAllCristiano10kRep1_Iterate() function
# NOTE: Remove healthy samples and the 1 duodenal cancer sample
metaCristianoTxNaiveFilt_CancerOnly <- metaCristianoTxNaiveFilt %>% 
  # filter(!(phenotype %in% c("Healthy","Duodenal Cancer"))) %>%
  filter(!(phenotype %in% c("Healthy"))) %>%
  droplevels()
cristianoBtwnCancerIterate <- ml1VsAllCristiano10kRep1_Iterate(metaData = metaCristianoTxNaiveFilt_CancerOnly,
                                                               col2Predict = "phenotype")

cristianoBtwnCancerIterate_results <- cristianoBtwnCancerIterate$rep_perf
cristianoBtwnCancerIterate_results$nullAUROC <- 0.5
cristianoBtwnCancerIterate_results$nullAUPR <- cristianoBtwnCancerIterate_results$minorityClassSize/
  (cristianoBtwnCancerIterate_results$minorityClassSize+cristianoBtwnCancerIterate_results$majorityClassSize)

colnames(cristianoBtwnCancerIterate_results)[1:2] <- c("AUROC","AUPR")
cristianoBtwnCancerIterate_results$dataString[cristianoBtwnCancerIterate_results$dataString=="full_rep200"] <- "Full multikingdom database (rep200)"
cristianoBtwnCancerIterate_results$dataString[cristianoBtwnCancerIterate_results$dataString=="bacteria_only"] <- "Bacteria all (Species)"
cristianoBtwnCancerIterate_results$dataString[cristianoBtwnCancerIterate_results$dataString=="fungi_only"] <- "Fungi all (Species)"
cristianoBtwnCancerIterate_results$dataString[cristianoBtwnCancerIterate_results$dataString=="fungi_decontam"] <- "Fungi decontaminated (Species)"
cristianoBtwnCancerIterate_results$dataString[cristianoBtwnCancerIterate_results$dataString=="fungi_intersected_with_Weizmann"] <- "Fungi ∩ WIS (Species)"
cristianoBtwnCancerIterate_results$dataString[cristianoBtwnCancerIterate_results$dataString=="bacteria_intersected_with_Weizmann"] <- "Bacteria ∩ WIS (Species)"
cristianoBtwnCancerIterate_results$dataString[cristianoBtwnCancerIterate_results$dataString=="fungi_and_bacteria_intersected_with_Weizmann"] <- "Fungi+bacteria ∩ WIS (Species)"
cristianoBtwnCancerIterate_results$dataString[cristianoBtwnCancerIterate_results$dataString=="topX_fungi"] <- "Top 20 fungi (Species)"
cristianoBtwnCancerIterate_results$dataString <- factor(cristianoBtwnCancerIterate_results$dataString, levels = c("Full multikingdom database (rep200)",
                                                                                                                "Fungi+bacteria ∩ WIS (Species)",
                                                                                                                "Bacteria all (Species)",
                                                                                                                "Fungi all (Species)",
                                                                                                                "Fungi decontaminated (Species)",
                                                                                                                "Fungi ∩ WIS (Species)",
                                                                                                                "Bacteria ∩ WIS (Species)",
                                                                                                                "Top 20 fungi (Species)"))
cristianoBtwnCancerIterate_results$diseaseType <- gsub(" cancer| Cancer","",cristianoBtwnCancerIterate_results$diseaseType)

save(cristianoBtwnCancerIterate_results,
     file = "Interim_data/cristianoBtwnCancerIterate_results_13Nov21.RData")

## Overlay multiple data types on the same graph
# NOTE: TopX fungi were selected based on cancer vs healthy comparisons
# so they are not plotted here
source("Supporting_scripts/S05-SummarySE.R")
cristianoBtwnCancerIterate_results %>%
  filter(grepl("Full|\\+|decontam",dataString)) %>%
  # filter(grepl("\\+|Bacteria ∩",dataString)) %>%
  reshape2::melt(id.vars = c("rep","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=dataString)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid", position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Hopkins plasma validation cohort:\nOne cancer type vs all others") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Figure_5/hopkins_between_cancer_types_full_WIS_fungi_datasets_13Nov21.svg", 
       dpi = "retina", width = 8, height = 5, units = "in")

cristianoBtwnCancerIterate_results %>% 
  write.csv("Figures_data/Figure_5/hopkins_between_cancer_types_full_WIS_fungi_datasets_13Nov21.csv")

## Subset to decontaminated fungi
cristianoBtwnCancerIterate_results %>%
  filter(grepl("decontam",dataString)) %>%
  reshape2::melt(id.vars = c("rep","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  # guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Hopkins plasma validation cohort:\nOne cancer type vs all others") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Figure_5/hopkins_between_cancer_types_fungi_only_13Nov21.svg", 
       dpi = "retina", width = 6.5, height = 4, units = "in")

cristianoBtwnCancerIterate_results %>%
  filter(grepl("Full",dataString)) %>% 
  pull(AUROC) %>% ci()
# Estimate   CI lower   CI upper Std. Error 
# 0.90471525 0.87310032 0.93633018 0.01584752

cristianoBtwnCancerIterate_results %>%
  filter(grepl("\\+",dataString)) %>% 
  pull(AUROC) %>% ci()
# Estimate   CI lower   CI upper Std. Error 
# 0.78742888 0.75300854 0.82184923 0.01725378 

cristianoBtwnCancerIterate_results %>%
  filter(grepl("decontam",dataString)) %>% 
  pull(AUROC) %>% ci()
# Estimate   CI lower   CI upper Std. Error 
# 0.68447299 0.64259177 0.72635421 0.02099367

#---------------------------------------------------------------------------------------------------------------------------#
# Test ML perf at varying stages levels
#---------------------------------------------------------------------------------------------------------------------------#
source("00-Functions.R") # for the ml1VsAllCristiano10kRep1_Iterate_Stage() function
cristianoStageIterate <- ml1VsAllCristiano10kRep1_Iterate_Stage()

cristianoStageIterate_results <- cristianoStageIterate$rep_perf
colnames(cristianoStageIterate_results)[1:2] <- c("AUROC","AUPR")
cristianoStageIterate_results$nullAUROC <- 0.5
cristianoStageIterate_results$nullAUPR <- cristianoStageIterate_results$minorityClassSize/
  (cristianoStageIterate_results$minorityClassSize+cristianoStageIterate_results$majorityClassSize)

cristianoStageIterate_results$dataString[cristianoStageIterate_results$dataString=="full_rep200"] <- "Full multikingdom database (rep200)"
cristianoStageIterate_results$dataString[cristianoStageIterate_results$dataString=="bacteria_only"] <- "Bacteria all (Species)"
cristianoStageIterate_results$dataString[cristianoStageIterate_results$dataString=="fungi_only"] <- "Fungi all (Species)"
cristianoStageIterate_results$dataString[cristianoStageIterate_results$dataString=="fungi_decontam"] <- "Fungi decontaminated (Species)"
cristianoStageIterate_results$dataString[cristianoStageIterate_results$dataString=="fungi_intersected_with_Weizmann"] <- "Fungi ∩ WIS (Species)"
cristianoStageIterate_results$dataString[cristianoStageIterate_results$dataString=="bacteria_intersected_with_Weizmann"] <- "Bacteria ∩ WIS (Species)"
cristianoStageIterate_results$dataString[cristianoStageIterate_results$dataString=="fungi_and_bacteria_intersected_with_Weizmann"] <- "Fungi+bacteria ∩ WIS (Species)"
cristianoStageIterate_results$dataString[cristianoStageIterate_results$dataString=="topX_fungi"] <- "Top 20 fungi (Species)"
cristianoStageIterate_results$dataString <- factor(cristianoStageIterate_results$dataString, levels = c("Full multikingdom database (rep200)",
                                                                                                        "Fungi+bacteria ∩ WIS (Species)",
                                                                                                        "Bacteria all (Species)",
                                                                                                        "Fungi all (Species)",
                                                                                                        "Fungi decontaminated (Species)",
                                                                                                        "Fungi ∩ WIS (Species)",
                                                                                                        "Bacteria ∩ WIS (Species)",
                                                                                                        "Top 20 fungi (Species)"))
cristianoStageIterate_results$diseaseType <- gsub(" cancer| Cancer","",cristianoStageIterate_results$diseaseType)
cristianoStageIterate_results$stageNum[cristianoStageIterate_results$stageNum=="I"] <- "Stage I"
cristianoStageIterate_results$stageNum[cristianoStageIterate_results$stageNum=="II"] <- "Stage II"
cristianoStageIterate_results$stageNum[cristianoStageIterate_results$stageNum=="III"] <- "Stage III"
cristianoStageIterate_results$stageNum[cristianoStageIterate_results$stageNum=="IV"] <- "Stage IV"
cristianoStageIterate_results$stageNum <- ordered(cristianoStageIterate_results$stageNum,
                                                  levels = c("Stage I","Stage II","Stage III","Stage IV"))

save(cristianoStageIterate_results,
     file = "Interim_data/cristianoStageIterate_results_13Nov21.RData")

## Full, WIS, fungi
source("Supporting_scripts/S05-SummarySE.R")
cristianoStageIterate_results %>%
  # filter(grepl("Full|\\+|decontam|Top",dataString)) %>%
  filter(grepl("\\+|Bacteria ∩",dataString)) %>%
  reshape2::melt(id.vars = c("rep","dataString","stageNum","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","stageNum","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(stageNum,value, color=dataString)) +
  # ggplot(aes(stageNum,value, color=dataString)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Hopkins plasma validation cohort:\nPer cancer stage vs. healthy performance") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Figure_5/figure_5_cristiano_per_stage_vs_healthy_full_WIS_fungi_TopX_13Nov21.svg", 
       dpi = "retina", width = 7, height = 4.5, units = "in")

cristianoStageIterate_results %>% write.csv("Figures_data/Figure_5/figure_5_cristiano_per_cancer_stage_vs_healthy_all_datasets_13Nov21.csv")

# cristianoStageIterate_results %>%
#   # filter(grepl("Full|\\+|decontam|Top",dataString)) %>%
#   filter(grepl("\\+|Bacteria ∩",dataString)) %>%
#   # filter(stageNum == "Stage III") %>%
#   ggboxplot(x = "dataString",
#             y = "AUROC") +
#   stat_compare_means(method = "wilcox")

## Subset to decontaminated fungi
cristianoStageIterate_results %>%
  filter(grepl("decontam",dataString)) %>%
  reshape2::melt(id.vars = c("rep","dataString","stageNum","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","stageNum","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  # ggplot(aes(reorder(stageNum, value, FUN=median),value, color=variable)) +
  ggplot(aes(stageNum,value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  # guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Hopkins plasma validation cohort:\nPer cancer stage vs. healthy") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Figure_5/figure_5_cristiano_per_stage_fungi_decontam_only_13Nov21.svg", 
       dpi = "retina", width = 4, height = 3.5, units = "in")

## Subset to decontaminated fungi
cristianoStageIterate_results %>%
  filter(grepl("decontam|Top",dataString)) %>%
  reshape2::melt(id.vars = c("rep","dataString","stageNum","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","stageNum","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  # ggplot(aes(reorder(stageNum, value, FUN=median),value, color=variable)) +
  ggplot(aes(stageNum,value, color=dataString)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer stage") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  # guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Hopkins plasma validation cohort:\nPer cancer stage vs. healthy") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Figure_5/figure_5_cristiano_per_stage_fungi_decontam_TopX_13Nov21.svg", 
       dpi = "retina", width = 7, height = 4, units = "in")

## Plot perf across stages
# Calculate one-way ANOVAs with adjusted p-vals
# NOTE: rstatix syntax does not play nicely with anova_test()
# since 2 groups are technically not compared, so the test
# results are converted to a df followed by manual plotting
require(rstatix)
roc.stat.test <- cristianoStageIterate_results %>%
  filter(grepl("Full|\\+|decontam|Top",dataString)) %>%
  group_by(dataString) %>%
  anova_test(AUROC ~ stageNum) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  data.frame() %>%
  mutate(p.adj = round(p.adj,4)) %>%
  mutate(group1="stageNum", group2 = "stageNum") %>%
  mutate(xmin = 2)

cristianoStageIterate_results %>%
  filter(grepl("Full|\\+|decontam|Top",dataString)) %>%
  ggboxplot(x = "stageNum",
            y = "AUROC",
            fill = "stageNum",
            xlab = "Cancer stage",
            title = "Hopkins: Cancer stage vs. healthy",
            palette = "nejm",
            facet.by = "dataString",
            legend = "none",
            add = "jitter",
            ylim = c(0.4, 1.2),
            add.params = list(alpha=0.4),
            notch = FALSE) +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_pvalue_manual(roc.stat.test, 
                     x = "xmin",
                     remove.bracket = TRUE,
                     label = "ANOVA padj = {p.adj}",
                     y.position = 1.1)
ggsave(filename = "Figures/Figure_5/hopkins_auroc_across_stages_16Nov21.svg",
       dpi = "retina", units = "in", width = 6.5, height = 5)

pr.stat.test <- cristianoStageIterate_results %>%
  filter(grepl("Full|\\+|decontam|Top",dataString)) %>%
  group_by(dataString) %>%
  anova_test(AUPR ~ stageNum) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  data.frame() %>%
  mutate(p.adj = round(p.adj,4)) %>%
  mutate(group1="stageNum", group2 = "stageNum") %>%
  mutate(xmin = 2)

cristianoStageIterate_results %>%
  filter(grepl("Full|\\+|decontam|Top",dataString)) %>%
  ggboxplot(x = "stageNum",
            y = "AUPR",
            fill = "stageNum",
            xlab = "Cancer stage",
            title = "Hopkins: Cancer stage vs. healthy",
            palette = "nejm",
            facet.by = "dataString",
            legend = "none",
            add = "jitter",
            ylim = c(0, 1.2),
            add.params = list(alpha=0.4),
            notch = FALSE) +
  theme(plot.title = element_text(hjust=0.5)) +
  stat_pvalue_manual(pr.stat.test, 
                     x = "xmin",
                     remove.bracket = TRUE,
                     label = "ANOVA padj = {p.adj}",
                     y.position = 1.1)
ggsave(filename = "Figures/Figure_5/hopkins_aupr_across_stages_16Nov21.svg",
       dpi = "retina", units = "in", width = 6.5, height = 5)
#---------------------------------------------------------------------------------------------------------------------------#
# Stage ML repeated CV perf --> ROC plot with 95% CIs
#---------------------------------------------------------------------------------------------------------------------------#
source("00-Functions.R")
stage_mlFullDB_CV <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                           countData = cristiano_rep200Data_Filt,
                           dataString = "full_rep200",
                           col2Predict = "stage",
                           stageNum = "I",
                           numResampleIter = 10,
                           varImpFlag = FALSE)
stage_mlFungiBacteriaWISIntersect_CV <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                              countData = cbind(cristianoRep200FungiSpeciesShared, 
                                                                cristianoRep200BacteriaSpeciesShared),
                                              dataString = "fungi_and_bacteria_intersected_with_Weizmann",
                                              col2Predict = "stage",
                                              stageNum = "I",
                                              numResampleIter = 10,
                                              varImpFlag = FALSE)
stage_mlFungiDecontam_CV <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                  countData = cristiano_rep200Data_Filt_Fungi_Decontam,
                                  dataString = "fungi_intersected_with_Weizmann",
                                  col2Predict = "stage",
                                  stageNum = "I",
                                  numResampleIter = 10,
                                  varImpFlag = FALSE)

source("00-Functions.R")
plots_stage_mlFullDB_CV <- plotMLWithCIs(stage_mlFullDB_CV, 
                                   showRepCurves=TRUE,
                                   sizeAvgCurve=0.5, 
                                   sizeRepCurves = 0.25,
                                   ciAlpha = 0.3,
                                   colorAvgCurve = "#3B4992FF",
                                   ciFillColor = "blue",
                                   colorRepCurves="lightblue",
                                   ciLevel = 0.99, 
                                   positiveClass = "I",
                                   negativeClass = "Healthy")
plots_stage_mlFungiBacteriaWISIntersect_CV <- plotMLWithCIs(stage_mlFungiBacteriaWISIntersect_CV, 
                                                      showRepCurves=TRUE, 
                                                      sizeAvgCurve=0.5, 
                                                      sizeRepCurves = 0.25,
                                                      ciAlpha = 0.3, 
                                                      colorAvgCurve = "#BB0021FF",
                                                      ciFillColor = "red",
                                                      colorRepCurves="darksalmon",
                                                      ciLevel = 0.99,
                                                      positiveClass = "I",
                                                      negativeClass = "Healthy")
plots_stage_mlFungiDecontam_CV <- plotMLWithCIs(stage_mlFungiDecontam_CV,
                                          sizeAvgCurve=0.5, 
                                          showRepCurves=TRUE, 
                                          ciAlpha = 0.3, 
                                          sizeRepCurves = 0.25,
                                          colorAvgCurve = "#008B45FF",
                                          ciFillColor = "lightgreen",
                                          colorRepCurves="lightgreen",
                                          ciLevel = 0.99,
                                          positiveClass = "I",
                                          negativeClass = "Healthy")

# save(stage_mlFullDB_CV,
#      stage_mlFungiDecontam_CV,
#      stage_mlFungiBacteriaWISIntersect_CV,
#      plots_mlFullDB_CV,
#      plots_stage_mlFungiBacteriaWISIntersect_CV,
#      plots_mlFungiDecontam_CV,
#      file = "Interim_data/cristiano_stageI_CV10_overlay_data_and_plots_11Nov21.RData")

# Combine ROC curves
rocCombinedStage <- plots_stage_mlFungiBacteriaWISIntersect_CV$rocPlot
for(ii in 1:10){
  rocCombinedStage <- rocCombinedStage + 
    geom_path(data = plots_stage_mlFullDB_CV$rocCurveData[[ii]],aes(x=fpr,y=tpr), color = "lightblue", size = 0.25) +
    geom_path(data = plots_stage_mlFungiDecontam_CV$rocCurveData[[ii]],aes(x=fpr,y=tpr), color = "lightgreen", size = 0.25)
}
rocCombinedStageAll <- rocCombinedStage +
  geom_ribbon(data = plots_stage_mlFullDB_CV$interpROCYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
              fill = "lightblue", alpha = 0.3, inherit.aes = F) +
  geom_path(data = plots_stage_mlFullDB_CV$interpROCYDf_CI, aes(x = xval, y = Estimate), color = "#3B4992FF", size = 0.5) +
  geom_ribbon(data = plots_stage_mlFungiDecontam_CV$interpROCYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
              fill = "lightgreen", alpha = 0.3, inherit.aes = F) + theme_bw() +
  geom_path(data = plots_stage_mlFungiDecontam_CV$interpROCYDf_CI, aes(x = xval, y = Estimate), color = "#008B45FF", size = 0.5) +
  annotate("segment", x = 0.35, xend = 0.45, y = 0.2, yend = 0.2, color = "#3B4992FF") +
  annotate("text", x = 0.47, y = 0.2, color = "#3B4992FF", label = paste0("AUROC 99% CI: [", 
                                                                          paste0(100*round(plots_stage_mlFullDB_CV$aurocCI[2],4),
                                                                                 ", ",100*round(plots_stage_mlFullDB_CV$aurocCI[3],4)),"]"), hjust = 0) +
  annotate("segment", x = 0.35, xend = 0.45, y = 0.15, yend = 0.15, color = "#BB0021FF") +
  annotate("text", x = 0.47, y = 0.15, color = "#BB0021FF", label = paste0("AUROC 99% CI: [", 
                                                                           paste0(100*round(plots_stage_mlFungiBacteriaWISIntersect_CV$aurocCI[2],4),
                                                                                  ", ",100*round(plots_stage_mlFungiBacteriaWISIntersect_CV$aurocCI[3],4)),"]"), hjust = 0) +
  annotate("segment", x = 0.35, xend = 0.45, y = 0.1, yend = 0.1, color = "#008B45FF") +
  annotate("text", x = 0.47, y = 0.1, color = "#008B45FF", label = paste0("AUROC 99% CI: [", 
                                                                          paste0(100*round(plots_stage_mlFungiDecontam_CV$aurocCI[2],4),
                                                                                 ", ",100*round(plots_stage_mlFungiDecontam_CV$aurocCI[3],4)),"]"), hjust = 0)
rocCombinedStageAll

# Combine PR curves
prCombinedStage <- plots_stage_mlFungiBacteriaWISIntersect_CV$prPlot
for(ii in 1:10){
  prCombinedStage <- prCombinedStage + 
    geom_path(data = plots_stage_mlFullDB_CV$prCurveData[[ii]],aes(x=recall,y=precision), color = "lightblue", size = 0.25) +
    geom_path(data = plots_stage_mlFungiDecontam_CV$prCurveData[[ii]],aes(x=recall,y=precision), color = "lightgreen", size = 0.25)
}
prCombinedStageAll <- prCombinedStage +
  geom_ribbon(data = plots_stage_mlFullDB_CV$interpPRYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
              fill = "lightblue", alpha = 0.3, inherit.aes = F) +
  geom_path(data = plots_stage_mlFullDB_CV$interpPRYDf_CI, aes(x = xval, y = Estimate), color = "#3B4992FF", size = 0.5) +
  geom_ribbon(data = plots_stage_mlFungiDecontam_CV$interpPRYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
              fill = "lightgreen", alpha = 0.3, inherit.aes = F) + theme_bw() +
  geom_path(data = plots_stage_mlFungiDecontam_CV$interpPRYDf_CI, aes(x = xval, y = Estimate), color = "#008B45FF", size = 0.5) +
  annotate("segment", x = 0.15, xend = 0.25, y = 0.2, yend = 0.2, color = "#3B4992FF") +
  annotate("text", x = 0.27, y = 0.2, color = "#3B4992FF", label = paste0("AUPR 99% CI: [", 
                                                                          paste0(100*round(plots_stage_mlFullDB_CV$auprCI[2],4),
                                                                                 ", ",100*round(plots_stage_mlFullDB_CV$auprCI[3],4)),"]"), hjust = 0) +
  annotate("segment", x = 0.15, xend = 0.25, y = 0.15, yend = 0.15, color = "#BB0021FF") +
  annotate("text", x = 0.27, y = 0.15, color = "#BB0021FF", label = paste0("AUPR 99% CI: [", 
                                                                           paste0(100*round(plots_stage_mlFungiBacteriaWISIntersect_CV$auprCI[2],4),
                                                                                  ", ",100*round(plots_stage_mlFungiBacteriaWISIntersect_CV$auprCI[3],4)),"]"), hjust = 0) +
  annotate("segment", x = 0.15, xend = 0.25, y = 0.1, yend = 0.1, color = "#008B45FF") +
  annotate("text", x = 0.27, y = 0.1, color = "#008B45FF", label = paste0("AUPR 99% CI: [", 
                                                                          paste0(100*round(plots_stage_mlFungiDecontam_CV$auprCI[2],4),
                                                                                 ", ",100*round(plots_stage_mlFungiDecontam_CV$auprCI[3],4)),"]"), hjust = 0)
combinedStageOverlayPlots <- ggarrange(rocCombinedStageAll, prCombinedStageAll, ncol = 2)
combinedStagePlotTitle <- paste0("Hopkins plasma cohort: Stage I vs Healthy (10-fold CV repeated 10 times)")
combinedStageOverlayPlotsAnnotated <- annotate_figure(combinedStageOverlayPlots, 
                                         top = text_grob(combinedStagePlotTitle, 
                                                         color = "black", face = "bold", size = 14))
print(combinedStageOverlayPlotsAnnotated)
ggsave(filename = "Figures/Figure_5/roc_pr_CV_overlay_cristiano_stageI_vs_healthy_full_WIS_fungionly.svg",
       units = "in", width = 10, height = 5)

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
## Fungi decontaminated
cristianoFungiDecontamPhylum <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                    countData = cristianoRep200FungiDecontamPhylum,
                                    col2Predict = "cancer_status",
                                    dataString = "fungi_decontam_phylum",
                                    varImpFlag = FALSE)
cristianoFungiDecontamClass <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                   countData = cristianoRep200FungiDecontamClass,
                                   col2Predict = "cancer_status",
                                   dataString = "fungi_decontam_class",
                                   varImpFlag = FALSE)
cristianoFungiDecontamOrder <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                   countData = cristianoRep200FungiDecontamOrder,
                                   col2Predict = "cancer_status",
                                   dataString = "fungi_decontam_order",
                                   varImpFlag = FALSE)
cristianoFungiDecontamFamily <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                    countData = cristianoRep200FungiDecontamFamily,
                                    col2Predict = "cancer_status",
                                    dataString = "fungi_decontam_family",
                                    varImpFlag = FALSE)
cristianoFungiDecontamGenus <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                   countData = cristianoRep200FungiDecontamGenus,
                                   col2Predict = "cancer_status",
                                   dataString = "fungi_decontam_genus",
                                   varImpFlag = TRUE)
cristianoFungiDecontamSpecies <- mlCristiano(metaData = metaCristianoTxNaiveFilt,
                                     countData = cristianoRep200FungiDecontamSpecies,
                                     col2Predict = "cancer_status",
                                     dataString = "fungi_decontam_species",
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
                                          cristianoFungiSpecies$rep_perf,
                                          
                                          cristianoFungiDecontamPhylum$rep_perf,
                                          cristianoFungiDecontamClass$rep_perf,
                                          cristianoFungiDecontamOrder$rep_perf,
                                          cristianoFungiDecontamFamily$rep_perf,
                                          cristianoFungiDecontamGenus$rep_perf,
                                          cristianoFungiDecontamSpecies$rep_perf),
                                    taxaLevel = c(rep(c("Phylum","Class","Order","Family","Genus","Species"),each=10),
                                                  rep(c("Phylum","Class","Order","Family","Genus","Species"),each=10),
                                                  rep(c("Phylum","Class","Order","Family","Genus","Species"),each=10)))
colnames(cristianoTaxaLevel_results)[1:2] <- c("AUROC","AUPR")
cristianoTaxaLevel_results$nullAUROC <- 0.5
cristianoTaxaLevel_results$nullAUPR <- cristianoTaxaLevel_results$minorityClassSize/
  (cristianoTaxaLevel_results$minorityClassSize+cristianoTaxaLevel_results$majorityClassSize)
cristianoTaxaLevel_results$dataString[grepl("bacteria",cristianoTaxaLevel_results$dataString)] <- "Bacteria (Raw)"
cristianoTaxaLevel_results$dataString[grepl("fungi_decontam",cristianoTaxaLevel_results$dataString)] <- "Fungi decontaminated"
cristianoTaxaLevel_results$dataString[grepl("fungi_",cristianoTaxaLevel_results$dataString)] <- "Fungi (Raw)"
cristianoTaxaLevel_results$dataString <- factor(cristianoTaxaLevel_results$dataString,
                                                levels = c("Bacteria (Raw)",
                                                           "Fungi (Raw)",
                                                           "Fungi decontaminated"))
cristianoTaxaLevel_results$taxaLevel <- factor(cristianoTaxaLevel_results$taxaLevel,
                                               levels = c("Phylum","Class","Order","Family","Genus","Species"))

## All data types
cristianoTaxaLevel_results %>%
  # filter(grepl("decontam",dataString)) %>%
  reshape2::melt(id.vars = c("rep","dataString","taxaLevel","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","taxaLevel","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  # ggplot(aes(reorder(stageNum, value, FUN=median),value, color=variable)) +
  ggplot(aes(taxaLevel,value, color=dataString)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  # guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Hopkins plasma validation cohort:\nPer taxa level cancer vs. healthy") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Figure_5/figure_5_cristiano_per_taxa_level_bacteria_and_fungi_13Nov21.svg", 
       dpi = "retina", width = 8, height = 4, units = "in")

## Subset to decontaminated fungi
cristianoTaxaLevel_results %>%
  filter(grepl("decontam",dataString)) %>%
  reshape2::melt(id.vars = c("rep","dataString","taxaLevel","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","taxaLevel","dataString","minorityClassSize","majorityClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  # ggplot(aes(reorder(stageNum, value, FUN=median),value, color=variable)) +
  ggplot(aes(taxaLevel,value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  # guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Hopkins plasma validation cohort:\nPer taxa level cancer vs. healthy") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Figure_5/figure_5_cristiano_per_taxa_level_bacteria_and_fungi_13Nov21.svg", 
       dpi = "retina", width = 6, height = 4, units = "in")

#----------------------------------------------------------------------------------#
# Scrambled and shuffled negative controls
# For healthy vs cancer
#----------------------------------------------------------------------------------#

#-----------------------ML with scrambled and shuffled data-----------------------#
source("00-Functions.R")
hopkinsScrambledIterate <- ml1VsAllCristiano10kRep1_Iterate(metaData = metaCristianoTxNaiveFilt,
                                                            listOfData = list(cristiano_rep200Data_Filt_Fungi_Decontam),
                                                            dataStringList = c("scrambled"),
                                                            scrambleFlag = TRUE)
hopkinsShuffledIterate <- ml1VsAllCristiano10kRep1_Iterate(metaData = metaCristianoTxNaiveFilt,
                                                           listOfData = list(cristiano_rep200Data_Filt_Fungi_Decontam),
                                                           dataStringList = c("shuffled"),
                                                           shuffleFlag = TRUE)
#-----------------------Format data-----------------------#

hopkinsScrambledIterate_results <- hopkinsScrambledIterate$rep_perf
hopkinsShuffledIterate_results <- hopkinsShuffledIterate$rep_perf

hopkinsScrambledIterate_results$nullAUROC <- 0.5
hopkinsScrambledIterate_results$nullAUPR <- hopkinsScrambledIterate_results$minorityClassSize/
  (hopkinsScrambledIterate_results$minorityClassSize+hopkinsScrambledIterate_results$majorityClassSize)
hopkinsShuffledIterate_results$nullAUROC <- 0.5
hopkinsShuffledIterate_results$nullAUPR <- hopkinsShuffledIterate_results$minorityClassSize/
  (hopkinsShuffledIterate_results$minorityClassSize+hopkinsShuffledIterate_results$majorityClassSize)

colnames(hopkinsScrambledIterate_results)[1:2] <- c("AUROC","AUPR")
colnames(hopkinsShuffledIterate_results)[1:2] <- c("AUROC","AUPR")
hopkinsScrambledIterate_results$dataString[hopkinsScrambledIterate_results$dataString=="scrambled"] <- "Scrambled"
hopkinsShuffledIterate_results$dataString[hopkinsShuffledIterate_results$dataString=="shuffled"] <- "Shuffled"

hopkinsScrambledIterate_results$diseaseType <- gsub(" cancer| Cancer","",hopkinsScrambledIterate_results$diseaseType)
hopkinsShuffledIterate_results$diseaseType <- gsub(" cancer| Cancer","",hopkinsShuffledIterate_results$diseaseType)

#-----------------------Combine data from actual and controls-----------------------#
hopkinsOverlay_HvsC_ActualvsControls <- rbind(cristianoPerCancerIterate_results,
                                              hopkinsScrambledIterate_results,
                                              hopkinsShuffledIterate_results) %>%
  filter(grepl("decontaminated|Top|Scrambled|Shuffled",dataString)) %>% droplevels()
hopkinsOverlay_HvsC_ActualvsControls$statGroups <- ifelse(grepl("shuffled|scrambled",
                                                                hopkinsOverlay_HvsC_ActualvsControls$dataString),
                                                         yes = "Control", no = "Actual")
hopkinsOverlay_HvsC_ActualvsControls$dataString <- as.character(hopkinsOverlay_HvsC_ActualvsControls$dataString)
hopkinsOverlay_HvsC_ActualvsControls$dataString[grepl("Top",hopkinsOverlay_HvsC_ActualvsControls$dataString)] <- "Top 20 fungi"
hopkinsOverlay_HvsC_ActualvsControls$dataString[grepl("Fungi",hopkinsOverlay_HvsC_ActualvsControls$dataString)] <- "Decontaminated\nfungi"
hopkinsOverlay_HvsC_ActualvsControls$dataString <- factor(hopkinsOverlay_HvsC_ActualvsControls$dataString,
                                                       levels = c("Decontaminated\nfungi",
                                                                  "Top 20 fungi",
                                                                  "Scrambled",
                                                                  "Shuffled"))
## Make grouped comparisons (all actual samples vs controls)
require(rstatix)
## AUROC
hopkinsOverlay_HvsC_ActualvsControls %>%
  distinct() %>% droplevels() %>%
  group_by(dataString) %>% data.frame() %>%
  wilcox_test(AUROC ~ dataString, exact = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "dataString") -> roc.tumor.stat.test

hopkinsOverlay_HvsC_ActualvsControls %>%
  distinct() %>% droplevels() %>%
  ggboxplot(x = "dataString",
            y = "AUROC",
            notch = TRUE,
            add = "jitter",
            add.params = list(alpha=0.4),
            legend = "none",
            fill = "dataString",
            xlab = "Data type",
            palette = "nejm") +
  rotate_x_text(30) + 
  stat_pvalue_manual(roc.tumor.stat.test, label = "q = {p.adj}", 
                     # y.position = c(1.03, 1.12, 1.06), 
                     y.position = c(1.03, 1.10,1.17, 1.24, 1.31, 1.45),
                     size = 3) +
  geom_hline(yintercept = 0.5, linetype="dotted") + 
  labs(fill = "Dataset") + theme(plot.title = element_text(hjust=0.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1.5)) -> aggregatedTumorROC
## AUPR
hopkinsOverlay_HvsC_ActualvsControls %>%
  distinct() %>% droplevels() %>%
  group_by(dataString) %>% data.frame() %>%
  wilcox_test(AUPR ~ dataString, exact = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "dataString") -> pr.tumor.stat.test

hopkinsOverlay_HvsC_ActualvsControls %>%
  distinct() %>% droplevels() %>%
  ggboxplot(x = "dataString",
            y = "AUPR",
            notch = TRUE,
            add = "jitter",
            add.params = list(alpha=0.4),
            legend = "none",
            fill = "dataString",
            xlab = "Data type",
            palette = "nejm") +
  rotate_x_text(30) + 
  stat_pvalue_manual(pr.tumor.stat.test, label = "q = {p.adj}", 
                     # y.position = c(1.03, 1.12, 1.06),
                     y.position = c(1.03, 1.10,1.17, 1.24, 1.31, 1.45),
                     size = 3) +
  labs(fill = "Dataset") + theme(plot.title = element_text(hjust=0.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1.5)) -> aggregatedTumorPR

combinedAggregatedTumorPlot <- ggarrange(aggregatedTumorROC, aggregatedTumorPR, ncol = 2) 
combinedAggregatedTumorPlotAnnotated <- annotate_figure(combinedAggregatedTumorPlot, 
                                                        top = text_grob("Hopkins: Aggregated cancer vs. healthy\nperformance with controls\n(Decontaminated fungi only)", 
                                                                        color = "black", face = "bold", size = 14))
print(combinedAggregatedTumorPlotAnnotated)
ggsave(filename = paste0("Figures/Figure_5/controls_hopkins_auroc_aupr.svg"),
       plot = combinedAggregatedTumorPlotAnnotated,
       dpi = "retina", units = "in", width = 8, height = 6)
