#-----------------------------------------------------------------------------
# 06-Perform-machine-learning-on-Weizmann-data.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Perform machine learning on Weizmann data
#-----------------------------------------------------------------------------

#----------------------------------------------------------#
# Load environments
#----------------------------------------------------------#
require(doMC)
require(plyr)
require(dplyr)
require(phyloseq)
require(microbiome)
require(vegan)
require(tibble)
require(biomformat)
require(rhdf5)
require(ggpubr)
require(ggsci)
require(scales)

numCores <- detectCores()
registerDoMC(cores=numCores)
#----------------------------------------------------------#
# Import data
#----------------------------------------------------------#

load("Interim_data/summarized_weizmann_data_various_taxa_levels_29Mar22.RData",
     verbose = TRUE)

#----------------------------------------------------------#
# Convert data to relative abundance and binary formats
# and add a taxa diversity column
#----------------------------------------------------------#

##----------------------All features----------------------##
# Relative abundances
ra_weizmannPhylum <- weizmannPhylum/rowSums(weizmannPhylum)
ra_weizmannClass <- weizmannClass/rowSums(weizmannClass)
ra_weizmannOrder <- weizmannOrder/rowSums(weizmannOrder)
ra_weizmannFamily <- weizmannFamily/rowSums(weizmannFamily)
ra_weizmannGenus <- weizmannGenus/rowSums(weizmannGenus)
ra_weizmannSpecies <- weizmannSpecies/rowSums(weizmannSpecies)
# Binary formats (based on presence/absence)
binary_weizmannPhylum <- data.frame((weizmannPhylum>0)*1)
binary_weizmannClass <- data.frame((weizmannClass>0)*1)
binary_weizmannOrder <- data.frame((weizmannOrder>0)*1)
binary_weizmannFamily <- data.frame((weizmannFamily>0)*1)
binary_weizmannGenus <- data.frame((weizmannGenus>0)*1)
binary_weizmannSpecies <- data.frame((weizmannSpecies>0)*1)
# Add diversity columns
ra_weizmannPhylum$diversity <- rowSums(ra_weizmannPhylum>0)
ra_weizmannClass$diversity <- rowSums(ra_weizmannClass>0)
ra_weizmannOrder$diversity <- rowSums(ra_weizmannOrder>0)
ra_weizmannFamily$diversity <- rowSums(ra_weizmannFamily>0)
ra_weizmannGenus$diversity <- rowSums(ra_weizmannGenus>0)
ra_weizmannSpecies$diversity <- rowSums(ra_weizmannSpecies>0)
binary_weizmannPhylum$diversity <- rowSums(binary_weizmannPhylum)
binary_weizmannClass$diversity <- rowSums(binary_weizmannClass)
binary_weizmannOrder$diversity <- rowSums(binary_weizmannOrder)
binary_weizmannFamily$diversity <- rowSums(binary_weizmannFamily)
binary_weizmannGenus$diversity <- rowSums(binary_weizmannGenus)
binary_weizmannSpecies$diversity <- rowSums(binary_weizmannSpecies)

##----------------------Shared features----------------------##
# Relative abundances
ra_weizmannPhylumShared_Nonzero <- weizmannPhylumShared_Nonzero/rowSums(weizmannPhylumShared_Nonzero)
ra_weizmannClassShared_Nonzero <- weizmannClassShared_Nonzero/rowSums(weizmannClassShared_Nonzero)
ra_weizmannOrderShared_Nonzero <- weizmannOrderShared_Nonzero/rowSums(weizmannOrderShared_Nonzero)
ra_weizmannFamilyShared_Nonzero <- weizmannFamilyShared_Nonzero/rowSums(weizmannFamilyShared_Nonzero)
ra_weizmannGenusShared_Nonzero <- weizmannGenusShared_Nonzero/rowSums(weizmannGenusShared_Nonzero)
ra_weizmannSpeciesShared_Nonzero <- weizmannSpeciesShared_Nonzero/rowSums(weizmannSpeciesShared_Nonzero)
# Binary formats (based on presence/absence)
binary_weizmannPhylumShared_Nonzero <- data.frame((weizmannPhylumShared_Nonzero>0)*1)
binary_weizmannClassShared_Nonzero <- data.frame((weizmannClassShared_Nonzero>0)*1)
binary_weizmannOrderShared_Nonzero <- data.frame((weizmannOrderShared_Nonzero>0)*1)
binary_weizmannFamilyShared_Nonzero <- data.frame((weizmannFamilyShared_Nonzero>0)*1)
binary_weizmannGenusShared_Nonzero <- data.frame((weizmannGenusShared_Nonzero>0)*1)
binary_weizmannSpeciesShared_Nonzero <- data.frame((weizmannSpeciesShared_Nonzero>0)*1)
# Add diversity columns
ra_weizmannPhylumShared_Nonzero$diversity <- rowSums(ra_weizmannPhylumShared_Nonzero>0)
ra_weizmannClassShared_Nonzero$diversity <- rowSums(ra_weizmannClassShared_Nonzero>0)
ra_weizmannOrderShared_Nonzero$diversity <- rowSums(ra_weizmannOrderShared_Nonzero>0)
ra_weizmannFamilyShared_Nonzero$diversity <- rowSums(ra_weizmannFamilyShared_Nonzero>0)
ra_weizmannGenusShared_Nonzero$diversity <- rowSums(ra_weizmannGenusShared_Nonzero>0)
ra_weizmannSpeciesShared_Nonzero$diversity <- rowSums(ra_weizmannSpeciesShared_Nonzero>0)
binary_weizmannPhylumShared_Nonzero$diversity <- rowSums(binary_weizmannPhylumShared_Nonzero)
binary_weizmannClassShared_Nonzero$diversity <- rowSums(binary_weizmannClassShared_Nonzero)
binary_weizmannOrderShared_Nonzero$diversity <- rowSums(binary_weizmannOrderShared_Nonzero)
binary_weizmannFamilyShared_Nonzero$diversity <- rowSums(binary_weizmannFamilyShared_Nonzero)
binary_weizmannGenusShared_Nonzero$diversity <- rowSums(binary_weizmannGenusShared_Nonzero)
binary_weizmannSpeciesShared_Nonzero$diversity <- rowSums(binary_weizmannSpeciesShared_Nonzero)

##----------------------Add diversity columns to count data----------------------##
## Add diversity columns -- CANNOT add before the above code (makes the the RA calculations incorrect)
weizmannPhylum$diversity <- rowSums(weizmannPhylum>0)
weizmannClass$diversity <- rowSums(weizmannClass>0)
weizmannOrder$diversity <- rowSums(weizmannOrder>0)
weizmannFamily$diversity <- rowSums(weizmannFamily>0)
weizmannGenus$diversity <- rowSums(weizmannGenus>0)
weizmannSpecies$diversity <- rowSums(weizmannSpecies>0)

weizmannPhylumShared_Nonzero$diversity <- rowSums(weizmannPhylumShared_Nonzero>0)
weizmannClassShared_Nonzero$diversity <- rowSums(weizmannClassShared_Nonzero>0)
weizmannOrderShared_Nonzero$diversity <- rowSums(weizmannOrderShared_Nonzero>0)
weizmannFamilyShared_Nonzero$diversity <- rowSums(weizmannFamilyShared_Nonzero>0)
weizmannGenusShared_Nonzero$diversity <- rowSums(weizmannGenusShared_Nonzero>0)
weizmannSpeciesShared_Nonzero$diversity <- rowSums(weizmannSpeciesShared_Nonzero>0)

save(weizmannPhylum,
     ra_weizmannPhylum,
     binary_weizmannPhylum,
     weizmannMetaPhylum,
     
     weizmannClass,
     ra_weizmannClass,
     binary_weizmannClass,
     weizmannMetaClass,
     
     weizmannOrder,
     ra_weizmannOrder,
     binary_weizmannOrder,
     weizmannMetaOrder,
     
     weizmannFamily,
     ra_weizmannFamily,
     binary_weizmannFamily,
     weizmannMetaFamily,
     
     weizmannGenus,
     ra_weizmannGenus,
     binary_weizmannGenus,
     weizmannMetaGenus,
     
     weizmannSpecies,
     ra_weizmannSpecies,
     binary_weizmannSpecies,
     weizmannMetaSpecies,
     
     weizmannPhylumShared_Nonzero,
     ra_weizmannPhylumShared_Nonzero,
     binary_weizmannPhylumShared_Nonzero,
     weizmannMetaPhylumShared_Nonzero,
     
     weizmannClassShared_Nonzero,
     ra_weizmannClassShared_Nonzero,
     binary_weizmannClassShared_Nonzero,
     weizmannMetaClassShared_Nonzero,
     
     weizmannOrderShared_Nonzero,
     ra_weizmannOrderShared_Nonzero,
     binary_weizmannOrderShared_Nonzero,
     weizmannMetaOrderShared_Nonzero,
     
     weizmannFamilyShared_Nonzero,
     ra_weizmannFamilyShared_Nonzero,
     binary_weizmannFamilyShared_Nonzero,
     weizmannMetaFamilyShared_Nonzero,
     
     weizmannGenusShared_Nonzero,
     ra_weizmannGenusShared_Nonzero,
     binary_weizmannGenusShared_Nonzero,
     weizmannMetaGenusShared_Nonzero,
     
     weizmannSpeciesShared_Nonzero,
     ra_weizmannSpeciesShared_Nonzero,
     binary_weizmannSpeciesShared_Nonzero,
     weizmannMetaSpeciesShared_Nonzero,
     
     file = "Interim_data/data_for_ml_weizmann_2Apr22.RData")

#----------------------------------------------------------#
# Run ML on Weizmann data
#----------------------------------------------------------#

# Use "Interim_data/data_for_ml_weizmann_2Apr22" to run
# "Supporting_scripts/S16R-ML-fungi-10k-rep1-weizmann.R"
# Which was run on a Slurm compute cluster using the submission script
# "Supporting_scripts/S16R-Submit-job.sh"
# Then download the final CSV file
# which has been placed under the "Interim_data" folder

# source('test-S11-ML-fungi-10k-rep1-weizmann.R')

#----------------------------------------------------------#
# Plot ML performance
#----------------------------------------------------------#

wzMLres <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_weizmann_with_and_without_shared_feat_ALL_2Apr22.csv")
wzMLres <- wzMLres[,!(colnames(wzMLres) == "X")]
colnames(wzMLres)[1:2] <- c("AUROC","AUPR")
# Phylum
wzMLres$datasetName[wzMLres$datasetName == "ra_weizmannPhylum"] <- "Phylum (RA)"
wzMLres$datasetName[wzMLres$datasetName == "binary_weizmannPhylum"] <- "Phylum (Binary)"
wzMLres$datasetName[wzMLres$datasetName == "weizmannPhylum"] <- "Phylum (Raw Counts)"
wzMLres$datasetName[wzMLres$datasetName == "ra_weizmannPhylumShared_Nonzero"] <- "Phylum Shared (RA)"
wzMLres$datasetName[wzMLres$datasetName == "binary_weizmannPhylumShared_Nonzero"] <- "Phylum Shared (Binary)"
wzMLres$datasetName[wzMLres$datasetName == "weizmannPhylumShared_Nonzero"] <- "Phylum Shared (Raw Counts)"
# Class
wzMLres$datasetName[wzMLres$datasetName == "ra_weizmannClass"] <- "Class (RA)"
wzMLres$datasetName[wzMLres$datasetName == "binary_weizmannClass"] <- "Class (Binary)"
wzMLres$datasetName[wzMLres$datasetName == "weizmannClass"] <- "Class (Raw Counts)"
wzMLres$datasetName[wzMLres$datasetName == "ra_weizmannClassShared_Nonzero"] <- "Class Shared (RA)"
wzMLres$datasetName[wzMLres$datasetName == "binary_weizmannClassShared_Nonzero"] <- "Class Shared (Binary)"
wzMLres$datasetName[wzMLres$datasetName == "weizmannClassShared_Nonzero"] <- "Class Shared (Raw Counts)"
# Order
wzMLres$datasetName[wzMLres$datasetName == "ra_weizmannOrder"] <- "Order (RA)"
wzMLres$datasetName[wzMLres$datasetName == "binary_weizmannOrder"] <- "Order (Binary)"
wzMLres$datasetName[wzMLres$datasetName == "weizmannOrder"] <- "Order (Raw Counts)"
wzMLres$datasetName[wzMLres$datasetName == "ra_weizmannOrderShared_Nonzero"] <- "Order Shared (RA)"
wzMLres$datasetName[wzMLres$datasetName == "binary_weizmannOrderShared_Nonzero"] <- "Order Shared (Binary)"
wzMLres$datasetName[wzMLres$datasetName == "weizmannOrderShared_Nonzero"] <- "Order Shared (Raw Counts)"
# Family
wzMLres$datasetName[wzMLres$datasetName == "ra_weizmannFamily"] <- "Family (RA)"
wzMLres$datasetName[wzMLres$datasetName == "binary_weizmannFamily"] <- "Family (Binary)"
wzMLres$datasetName[wzMLres$datasetName == "weizmannFamily"] <- "Family (Raw Counts)"
wzMLres$datasetName[wzMLres$datasetName == "ra_weizmannFamilyShared_Nonzero"] <- "Family Shared (RA)"
wzMLres$datasetName[wzMLres$datasetName == "binary_weizmannFamilyShared_Nonzero"] <- "Family Shared (Binary)"
wzMLres$datasetName[wzMLres$datasetName == "weizmannFamilyShared_Nonzero"] <- "Family Shared (Raw Counts)"
# Genus
wzMLres$datasetName[wzMLres$datasetName == "ra_weizmannGenus"] <- "Genus (RA)"
wzMLres$datasetName[wzMLres$datasetName == "binary_weizmannGenus"] <- "Genus (Binary)"
wzMLres$datasetName[wzMLres$datasetName == "weizmannGenus"] <- "Genus (Raw Counts)"
wzMLres$datasetName[wzMLres$datasetName == "ra_weizmannGenusShared_Nonzero"] <- "Genus Shared (RA)"
wzMLres$datasetName[wzMLres$datasetName == "binary_weizmannGenusShared_Nonzero"] <- "Genus Shared (Binary)"
wzMLres$datasetName[wzMLres$datasetName == "weizmannGenusShared_Nonzero"] <- "Genus Shared (Raw Counts)"
# Species
wzMLres$datasetName[wzMLres$datasetName == "ra_weizmannSpecies"] <- "Species (RA)"
wzMLres$datasetName[wzMLres$datasetName == "binary_weizmannSpecies"] <- "Species (Binary)"
wzMLres$datasetName[wzMLres$datasetName == "weizmannSpecies"] <- "Species (Raw Counts)"
wzMLres$datasetName[wzMLres$datasetName == "ra_weizmannSpeciesShared_Nonzero"] <- "Species Shared (RA)"
wzMLres$datasetName[wzMLres$datasetName == "binary_weizmannSpeciesShared_Nonzero"] <- "Species Shared (Binary)"
wzMLres$datasetName[wzMLres$datasetName == "weizmannSpeciesShared_Nonzero"] <- "Species Shared (Raw Counts)"

wzMLres$datasetName <- factor(wzMLres$datasetName,
                                           levels = c("Phylum (Raw Counts)", "Phylum (RA)", "Phylum (Binary)",
                                                      "Phylum Shared (Raw Counts)", "Phylum Shared (RA)", "Phylum Shared (Binary)",
                                                      
                                                      "Class (Raw Counts)", "Class (RA)", "Class (Binary)",
                                                      "Class Shared (Raw Counts)", "Class Shared (RA)", "Class Shared (Binary)",
                                                      
                                                      "Order (Raw Counts)", "Order (RA)", "Order (Binary)",
                                                      "Order Shared (Raw Counts)", "Order Shared (RA)", "Order Shared (Binary)",
                                                      
                                                      "Family (Raw Counts)", "Family (RA)", "Family (Binary)",
                                                      "Family Shared (Raw Counts)", "Family Shared (RA)", "Family Shared (Binary)",
                                                      
                                                      "Genus (Raw Counts)", "Genus (RA)", "Genus (Binary)",
                                                      "Genus Shared (Raw Counts)", "Genus Shared (RA)", "Genus Shared (Binary)",
                                                      
                                                      "Species (Raw Counts)", "Species (RA)", "Species (Binary)",
                                                      "Species Shared (Raw Counts)", "Species Shared (RA)", "Species Shared (Binary)"))
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
# There are lots of options to play with here. For purposes of the paper, we are focusing on the
# best performing option: full binary (presence/absence) data at the species level
source('Supporting_scripts/S00-SummarySE.R')
wzMLres %>%
  filter(sampleType == "tumor") %>%
  filter(!grepl("Shared",datasetName)) %>%
  filter(grepl("Binary",datasetName)) %>%
  filter(grepl("Species",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("diseaseType","datasetName","variable")) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Tumor | 1 Vs All | Binary (Presence/Absence) data | Species Level") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") +
  geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/ml_weizmann_tumor_1vsAll_species_binary.pdf", dpi = "retina",
         width = 10, height = 6, units = "in")

# There are lots of options to play with here. For purposes of the paper, we are focusing on the
# best performing option: full binary (presence/absence) data at the species level
source('Supporting_scripts/S00-SummarySE.R')
wzMLres %>%
  filter(sampleType == "tumor") %>%
  filter(!grepl("Shared",datasetName)) %>%
  # filter(grepl("Binary",datasetName)) %>%
  filter(grepl("Species",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("diseaseType","datasetName","variable")) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Tumor | 1 Vs All | Binary (Presence/Absence) data | Species Level") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") +
  geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/ml_weizmann_tumor_1vsAll_species_counts_ra_binary.pdf", dpi = "retina",
         width = 10, height = 6, units = "in")

#-------------------------Plot primary tumor vs. nat performance-------------------------#
# No options appeared to differentiate tumor vs. nat well here, so all are shown
source('Supporting_scripts/S00-SummarySE.R')
wzMLres %>%
  filter(sampleType == "tumor vs nat") %>%
  filter(!grepl("Shared",datasetName)) %>%
  # filter(grepl("Binary",datasetName)) %>%
  # filter(grepl("Species",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("diseaseType","datasetName","variable")) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Tumor vs. NAT | Phylum, Class, Order, Family, Genus, Species Levels") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_igv(name = "Features") +
  geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/ml_weizmann_tumor_vs_nat_all.pdf", dpi = "retina",
         width = 10, height = 6, units = "in")

