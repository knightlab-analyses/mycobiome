#-----------------------------------------------------------------------------
# 03-Intersect-TCGA-data-with-Weizmann-data.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Intersect TCGA fungi data with Weizmann fungi data at each taxa level
# - Calculate and plot shared overlap at genus and species levels
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

# NOTE: "fungal_norm_unique_hit.rds" contains a list with 7 phyloseq objects
# containing fungal data summarized at all_rank, phylum, class, order, family,
# genus, and species levels
wzRDS <- readRDS("Input_data/Weizmann_data/fungal_norm_unique_hit.rds")
pzWz_allRank <- wzRDS$all_rank
psWz_phylum <- wzRDS$phylum_phy
psWz_class <- wzRDS$class_phy
psWz_order <- wzRDS$order_phy
psWz_family <- wzRDS$family_phy
psWz_genus <- wzRDS$genus_phy
psWz_species <- wzRDS$species_phy

# Subset to non-control samples to compare overlap
pzWz_allRank_Bio <- subset_samples(pzWz_allRank, type.detail %in% c("normal", "nat", "tumor"))
psWz_phylum_Bio <- subset_samples(psWz_phylum, type.detail %in% c("normal", "nat", "tumor"))
psWz_class_Bio <- subset_samples(psWz_class, type.detail %in% c("normal", "nat", "tumor"))
psWz_order_Bio <- subset_samples(psWz_order, type.detail %in% c("normal", "nat", "tumor"))
psWz_family_Bio <- subset_samples(psWz_family, type.detail %in% c("normal", "nat", "tumor"))
psWz_genus_Bio <- subset_samples(psWz_genus, type.detail %in% c("normal", "nat", "tumor"))
psWz_species_Bio <- subset_samples(psWz_species, type.detail %in% c("normal", "nat", "tumor"))

## Load TCGA data
# metaQiitaCombined_Nonzero_WithBamcounts, 
# rep200Data_WGS_RNA_Matched_Bacteria,
# rep200Data_WGS_RNA_Matched_Fungi,
load("Interim_data/metaQiitaCombined_Nonzero_WithBamcounts_and_Data_13Sep21.RData") # To load the metadata and raw data objects

# Load the revised taxa mapping file. The following changes were made to compare to the Weizmann taxa mapping:
# - Names for taxa levels were edited (e.g. "Domain" --> "kingdom") and lowercased
# - Empty entries (e.g. "c__") were replaced with "other" (e.g. "other")
# - Entries with brackets (e.g. "[Candida] ...") had their brackets removed (e.g. "Candida ...")
# - The preceding taxa level symbol and double underscore (e.g. "c__XXX" for class) were removed (e.g. "XXX")
# - The entries under "kingdom" were changed from "Eukaryota" to "Fungi"
# - All spaces were converted to underscores
rep200TaxSplit_Fungi_Paired_to_Weizmann <- read.csv("Supporting_data/rep200TaxSplit_Fungi_Paired_To_Weizmann_Final.csv", stringsAsFactors = FALSE, row.names = 1)

#----------------------------------------------------------#
# Create TCGA phyloseq object and aggregate counts at various taxa levels
#----------------------------------------------------------#

psFungiHiSeqFungi_Paired2Wz <- phyloseq(otu_table(rep200Data_WGS_RNA_Matched_Fungi, taxa_are_rows = FALSE), 
                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_WithBamcounts))

## Aggregate counts for TCGA data
psFungiHiSeqFungi_Paired2Wz_phylum = aggregate_taxa(psFungiHiSeqFungi_Paired2Wz, "phylum")
psFungiHiSeqFungi_Paired2Wz_class = aggregate_taxa(psFungiHiSeqFungi_Paired2Wz, "class")
psFungiHiSeqFungi_Paired2Wz_order = aggregate_taxa(psFungiHiSeqFungi_Paired2Wz, "order")
psFungiHiSeqFungi_Paired2Wz_family = aggregate_taxa(psFungiHiSeqFungi_Paired2Wz, "family")
psFungiHiSeqFungi_Paired2Wz_genus = aggregate_taxa(psFungiHiSeqFungi_Paired2Wz, "genus")
psFungiHiSeqFungi_Paired2Wz_species = rep200Data_WGS_RNA_Matched_Fungi
colnames(psFungiHiSeqFungi_Paired2Wz_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(psFungiHiSeqFungi_Paired2Wz_species),"species"]
psFungiHiSeqFungi_Paired2Wz_species[1:3,1:3]

#----------------------------------------------------------#
# Intersect features at varying taxa levels
#----------------------------------------------------------#

source("00-Functions.R") # for findSharedTaxa() function

sharedPhylum <- findSharedTaxa(psFungiHiSeqFungi_Paired2Wz_phylum, psWz_phylum_Bio, "phylum")
sharedClass <- findSharedTaxa(psFungiHiSeqFungi_Paired2Wz_class, psWz_class_Bio, "class")
sharedOrder <- findSharedTaxa(psFungiHiSeqFungi_Paired2Wz_order, psWz_order_Bio, "order")
sharedFamily <- findSharedTaxa(psFungiHiSeqFungi_Paired2Wz_family, psWz_family_Bio, "family")
sharedGenus <- findSharedTaxa(psFungiHiSeqFungi_Paired2Wz_genus, psWz_genus_Bio, "genus")
sharedSpecies <- findSharedTaxa(psFungiHiSeqFungi_Paired2Wz_species, psWz_species_Bio, "species")

save(sharedPhylum,
     sharedClass,
     sharedOrder,
     sharedFamily,
     sharedGenus,
     sharedSpecies,
     file = "Interim_data/shared_fungi_features_at_each_taxa_level_13Sep21.RData")

#----------------------------------------------------------#
# Subset TCGA and Weizmann data by intersecting features at varying taxa levels
#----------------------------------------------------------#
## Subset TCGA data
psFungiHiSeqFungi_Paired2Wz_phylum_shared <- subset_taxa(psFungiHiSeqFungi_Paired2Wz_phylum, phylum %in% sharedPhylum)
psFungiHiSeqFungi_Paired2Wz_class_shared <- subset_taxa(psFungiHiSeqFungi_Paired2Wz_class, class %in% sharedClass)
psFungiHiSeqFungi_Paired2Wz_order_shared <- subset_taxa(psFungiHiSeqFungi_Paired2Wz_order, order %in% sharedOrder)
psFungiHiSeqFungi_Paired2Wz_family_shared <- subset_taxa(psFungiHiSeqFungi_Paired2Wz_family, family %in% sharedFamily)
psFungiHiSeqFungi_Paired2Wz_genus_shared <- subset_taxa(psFungiHiSeqFungi_Paired2Wz_genus, genus %in% sharedGenus)
psFungiHiSeqFungi_Paired2Wz_species_shared <- psFungiHiSeqFungi_Paired2Wz_species[,colnames(psFungiHiSeqFungi_Paired2Wz_species) %in% sharedSpecies]

## Extract TCGA count tables with shared features
metaQiitaCombined_Nonzero_shared <- metaQiitaCombined_Nonzero_WithBamcounts
rep200FungiPhylumShared <- data.frame(t(otu_table(psFungiHiSeqFungi_Paired2Wz_phylum_shared)))
rep200FungiClassShared <- data.frame(t(otu_table(psFungiHiSeqFungi_Paired2Wz_class_shared)))
rep200FungiOrderShared <- data.frame(t(otu_table(psFungiHiSeqFungi_Paired2Wz_order_shared)))
rep200FungiFamilyShared <- data.frame(t(otu_table(psFungiHiSeqFungi_Paired2Wz_family_shared)))
rep200FungiGenusShared <- data.frame(t(otu_table(psFungiHiSeqFungi_Paired2Wz_genus_shared)))
rep200FungiSpeciesShared <- psFungiHiSeqFungi_Paired2Wz_species_shared
rep200FungiSpeciesShared[1:3,1:3]

## After subsetting features, remove samples with 0 counts. Note that the metadata will have to be subsetted as well
# Phylum
rep200FungiPhylumShared_Nonzero <- rep200FungiPhylumShared[rowSums(rep200FungiPhylumShared) != 0,]
metaQiitaCombined_Nonzero_PhylumShared <- droplevels(metaQiitaCombined_Nonzero_WithBamcounts[rowSums(rep200FungiPhylumShared) != 0,])
# Class
rep200FungiClassShared_Nonzero <- rep200FungiClassShared[rowSums(rep200FungiClassShared) != 0,]
metaQiitaCombined_Nonzero_ClassShared <- droplevels(metaQiitaCombined_Nonzero_WithBamcounts[rowSums(rep200FungiClassShared) != 0,])
# Order
rep200FungiOrderShared_Nonzero <- rep200FungiOrderShared[rowSums(rep200FungiOrderShared) != 0,]
metaQiitaCombined_Nonzero_OrderShared <- droplevels(metaQiitaCombined_Nonzero_WithBamcounts[rowSums(rep200FungiOrderShared) != 0,])
# Family
rep200FungiFamilyShared_Nonzero <- rep200FungiFamilyShared[rowSums(rep200FungiFamilyShared) != 0,]
metaQiitaCombined_Nonzero_FamilyShared <- droplevels(metaQiitaCombined_Nonzero_WithBamcounts[rowSums(rep200FungiFamilyShared) != 0,])
# Genus
rep200FungiGenusShared_Nonzero <- rep200FungiGenusShared[rowSums(rep200FungiGenusShared) != 0,]
metaQiitaCombined_Nonzero_GenusShared <- droplevels(metaQiitaCombined_Nonzero_WithBamcounts[rowSums(rep200FungiGenusShared) != 0,])
# Species
rep200FungiSpeciesShared_Nonzero <- rep200FungiSpeciesShared[rowSums(rep200FungiSpeciesShared) != 0,]
metaQiitaCombined_Nonzero_SpeciesShared <- droplevels(metaQiitaCombined_Nonzero_WithBamcounts[rowSums(rep200FungiSpeciesShared) != 0,])


save(rep200FungiPhylumShared_Nonzero,
     metaQiitaCombined_Nonzero_PhylumShared,
     rep200FungiClassShared_Nonzero,
     metaQiitaCombined_Nonzero_ClassShared,
     rep200FungiOrderShared_Nonzero,
     metaQiitaCombined_Nonzero_OrderShared,
     rep200FungiFamilyShared_Nonzero,
     metaQiitaCombined_Nonzero_FamilyShared,
     rep200FungiGenusShared_Nonzero,
     metaQiitaCombined_Nonzero_GenusShared,
     rep200FungiSpeciesShared_Nonzero,
     metaQiitaCombined_Nonzero_SpeciesShared,
     file = "Interim_data/tcga_data_taxa_levels_features_shared_with_Weizmann_13Sep21.RData")

## Subset Weizmann data
psWz_phylum_Bio_shared <- subset_taxa(psWz_phylum_Bio, phylum %in% sharedPhylum)
psWz_class_Bio_shared <- subset_taxa(psWz_class_Bio, class %in% sharedClass)
psWz_order_Bio_shared <- subset_taxa(psWz_order_Bio, order %in% sharedOrder)
psWz_family_Bio_shared <- subset_taxa(psWz_family_Bio, family %in% sharedFamily)
psWz_genus_Bio_shared <- subset_taxa(psWz_genus_Bio, genus %in% sharedGenus)
psWz_species_Bio_shared <- subset_taxa(psWz_species_Bio, species %in% sharedSpecies)

## Convert Weizmann data with shared features to data frames and save corresponding metadata
## Need to remove samples with 0 counts after subseting features
# Phylum
weizmannPhylumShared <- data.frame(t(otu_table(psWz_phylum_Bio_shared)))
colnames(weizmannPhylumShared) <- data.frame(tax_table(psWz_phylum_Bio_shared))[colnames(weizmannPhylumShared),"phylum"]
weizmannMetaPhylumShared <- data.frame(sample_data(psWz_phylum_Bio_shared))
weizmannPhylumShared_Nonzero <- weizmannPhylumShared[rowSums(weizmannPhylumShared) != 0,]
weizmannMetaPhylumShared_Nonzero <- droplevels(weizmannMetaPhylumShared[rowSums(weizmannPhylumShared) != 0,])
# Class
weizmannClassShared <- data.frame(t(otu_table(psWz_class_Bio_shared)))
colnames(weizmannClassShared) <- data.frame(tax_table(psWz_class_Bio_shared))[colnames(weizmannClassShared),"class"]
weizmannMetaClassShared <- data.frame(sample_data(psWz_class_Bio_shared))
weizmannClassShared_Nonzero <- weizmannClassShared[rowSums(weizmannClassShared) != 0,]
weizmannMetaClassShared_Nonzero <- droplevels(weizmannMetaClassShared[rowSums(weizmannClassShared) != 0,])
# Order
weizmannOrderShared <- data.frame(t(otu_table(psWz_order_Bio_shared)))
colnames(weizmannOrderShared) <- data.frame(tax_table(psWz_order_Bio_shared))[colnames(weizmannOrderShared),"order"]
weizmannMetaOrderShared <- data.frame(sample_data(psWz_order_Bio_shared))
weizmannOrderShared_Nonzero <- weizmannOrderShared[rowSums(weizmannOrderShared) != 0,]
weizmannMetaOrderShared_Nonzero <- droplevels(weizmannMetaOrderShared[rowSums(weizmannOrderShared) != 0,])
# Family
weizmannFamilyShared <- data.frame(t(otu_table(psWz_family_Bio_shared)))
colnames(weizmannFamilyShared) <- data.frame(tax_table(psWz_family_Bio_shared))[colnames(weizmannFamilyShared),"family"]
weizmannMetaFamilyShared <- data.frame(sample_data(psWz_family_Bio_shared))
weizmannFamilyShared_Nonzero <- weizmannFamilyShared[rowSums(weizmannFamilyShared) != 0,]
weizmannMetaFamilyShared_Nonzero <- droplevels(weizmannMetaFamilyShared[rowSums(weizmannFamilyShared) != 0,])
# Genus
weizmannGenusShared <- data.frame(t(otu_table(psWz_genus_Bio_shared)))
colnames(weizmannGenusShared) <- data.frame(tax_table(psWz_genus_Bio_shared))[colnames(weizmannGenusShared),"genus"]
weizmannMetaGenusShared <- data.frame(sample_data(psWz_genus_Bio_shared))
weizmannGenusShared_Nonzero <- weizmannGenusShared[rowSums(weizmannGenusShared) != 0,]
weizmannMetaGenusShared_Nonzero <- droplevels(weizmannMetaGenusShared[rowSums(weizmannGenusShared) != 0,])
# Species
weizmannSpeciesShared <- data.frame(t(otu_table(psWz_species_Bio_shared)))
colnames(weizmannSpeciesShared) <- data.frame(tax_table(psWz_species_Bio_shared))[colnames(weizmannSpeciesShared),"species"]
weizmannMetaSpeciesShared <- data.frame(sample_data(psWz_species_Bio_shared))
weizmannSpeciesShared_Nonzero <- weizmannSpeciesShared[rowSums(weizmannSpeciesShared) != 0,]
weizmannMetaSpeciesShared_Nonzero <- droplevels(weizmannMetaSpeciesShared[rowSums(weizmannSpeciesShared) != 0,])

## Convert Weizmann data with ALL features to data frames and save corresponding metadata
# Phylum
weizmannPhylum <- data.frame(t(otu_table(psWz_phylum_Bio)))
colnames(weizmannPhylum) <- data.frame(tax_table(psWz_phylum_Bio))[colnames(weizmannPhylum),"phylum"]
weizmannMetaPhylum <- data.frame(sample_data(psWz_phylum_Bio))
# Class
weizmannClass <- data.frame(t(otu_table(psWz_class_Bio)))
colnames(weizmannClass) <- data.frame(tax_table(psWz_class_Bio))[colnames(weizmannClass),"class"]
weizmannMetaClass <- data.frame(sample_data(psWz_class_Bio))
# Order
weizmannOrder <- data.frame(t(otu_table(psWz_order_Bio)))
colnames(weizmannOrder) <- data.frame(tax_table(psWz_order_Bio))[colnames(weizmannOrder),"order"]
weizmannMetaOrder <- data.frame(sample_data(psWz_order_Bio))
# Family
weizmannFamily <- data.frame(t(otu_table(psWz_family_Bio)))
colnames(weizmannFamily) <- data.frame(tax_table(psWz_family_Bio))[colnames(weizmannFamily),"family"]
weizmannMetaFamily <- data.frame(sample_data(psWz_family_Bio))
# Genus
weizmannGenus <- data.frame(t(otu_table(psWz_genus_Bio)))
colnames(weizmannGenus) <- data.frame(tax_table(psWz_genus_Bio))[colnames(weizmannGenus),"genus"]
weizmannMetaGenus <- data.frame(sample_data(psWz_genus_Bio))
# Species
weizmannSpecies <- data.frame(t(otu_table(psWz_species_Bio)))
colnames(weizmannSpecies) <- data.frame(tax_table(psWz_species_Bio))[colnames(weizmannSpecies),"species"]
weizmannMetaSpecies <- data.frame(sample_data(psWz_species_Bio))


save(weizmannPhylumShared_Nonzero,
     weizmannClassShared_Nonzero,
     weizmannOrderShared_Nonzero,
     weizmannFamilyShared_Nonzero,
     weizmannGenusShared_Nonzero,
     weizmannSpeciesShared_Nonzero,
     weizmannMetaPhylumShared_Nonzero,
     weizmannMetaClassShared_Nonzero,
     weizmannMetaOrderShared_Nonzero,
     weizmannMetaFamilyShared_Nonzero,
     weizmannMetaGenusShared_Nonzero,
     weizmannMetaSpeciesShared_Nonzero,
     # Shared feature objects above and ALL feature objects below
     weizmannPhylum,
     weizmannClass,
     weizmannOrder,
     weizmannFamily,
     weizmannGenus,
     weizmannSpecies,
     weizmannMetaPhylum,
     weizmannMetaClass,
     weizmannMetaOrder,
     weizmannMetaFamily,
     weizmannMetaGenus,
     weizmannMetaSpecies,
     file = "Interim_data/summarized_weizmann_data_various_taxa_levels_15Sep21.RData")
# The "summarized_weizmann_data_various_taxa_levels_15Sep21.RData" file will be loaded
# in the "06-Perform-machine-learning-on-Weizmann-data.R" script and processed there for
# machine learning

#----------------------------------------------------------------------------------#
# Extract 8 cancer types to directly compare with Weizmann cancers and find shared features
#----------------------------------------------------------------------------------#
# Cancer types to include: breast, lung, melanoma, colon, GBM, pancreas, ovary, bone
# Combine: LUAD and LUSC, COAD and READ

metaQiitaCombined_Nonzero_8cancer_shared <- metaQiitaCombined_Nonzero_WithBamcounts %>%
  filter(investigation %in% c("TCGA-BRCA","TCGA-LUAD","TCGA-LUSC","TCGA-SKCM","TCGA-COAD","TCGA-READ","TCGA-GBM","TCGA-PAAD","TCGA-OV","TCGA-SARC")) %>%
  droplevels()
dim(metaQiitaCombined_Nonzero_8cancer_shared) # 5892   41
metaQiitaCombined_Nonzero_8cancer_shared$disease_type <- as.character(metaQiitaCombined_Nonzero_8cancer_shared$disease_type) 
metaQiitaCombined_Nonzero_8cancer_shared$disease_type[metaQiitaCombined_Nonzero_8cancer_shared$disease_type == "Lung Squamous Cell Carcinoma"] <- "Lung Cancer"
metaQiitaCombined_Nonzero_8cancer_shared$disease_type[metaQiitaCombined_Nonzero_8cancer_shared$disease_type == "Lung Adenocarcinoma"] <- "Lung Cancer"
metaQiitaCombined_Nonzero_8cancer_shared$disease_type[metaQiitaCombined_Nonzero_8cancer_shared$disease_type == "Colon Adenocarcinoma"] <- "Colorectal Cancer"
metaQiitaCombined_Nonzero_8cancer_shared$disease_type[metaQiitaCombined_Nonzero_8cancer_shared$disease_type == "Rectum Adenocarcinoma"] <- "Colorectal Cancer"
metaQiitaCombined_Nonzero_8cancer_shared$disease_type[metaQiitaCombined_Nonzero_8cancer_shared$disease_type == "Breast Invasive Carcinoma"] <- "Breast Cancer"
metaQiitaCombined_Nonzero_8cancer_shared$disease_type[metaQiitaCombined_Nonzero_8cancer_shared$disease_type == "Glioblastoma Multiforme"] <- "Glioblastoma"
metaQiitaCombined_Nonzero_8cancer_shared$disease_type[metaQiitaCombined_Nonzero_8cancer_shared$disease_type == "Pancreatic Adenocarcinoma"] <- "Pancreatic Cancer"
metaQiitaCombined_Nonzero_8cancer_shared$disease_type[metaQiitaCombined_Nonzero_8cancer_shared$disease_type == "Ovarian Serous Cystadenocarcinoma"] <- "Ovarian Cancer"
metaQiitaCombined_Nonzero_8cancer_shared$disease_type[metaQiitaCombined_Nonzero_8cancer_shared$disease_type == "Skin Cutaneous Melanoma"] <- "Melanoma"
metaQiitaCombined_Nonzero_8cancer_shared$disease_type[metaQiitaCombined_Nonzero_8cancer_shared$disease_type == "Sarcoma"] <- "Bone Cancer"
table(metaQiitaCombined_Nonzero_8cancer_shared$disease_type)

## TCGA data
rep200FungiPhylumShared_8cancer <- rep200FungiPhylumShared[rownames(metaQiitaCombined_Nonzero_8cancer_shared),]
rep200FungiClassShared_8cancer <- rep200FungiClassShared[rownames(metaQiitaCombined_Nonzero_8cancer_shared),]
rep200FungiOrderShared_8cancer <- rep200FungiOrderShared[rownames(metaQiitaCombined_Nonzero_8cancer_shared),]
rep200FungiFamilyShared_8cancer <- rep200FungiFamilyShared[rownames(metaQiitaCombined_Nonzero_8cancer_shared),]
rep200FungiGenusShared_8cancer <- rep200FungiGenusShared[rownames(metaQiitaCombined_Nonzero_8cancer_shared),]
rep200FungiSpeciesShared_8cancer <- rep200FungiSpeciesShared[rownames(metaQiitaCombined_Nonzero_8cancer_shared),]
rep200FungiGenusShared_8cancer[1:3,1:3]

## After subsetting features, remove samples with 0 counts. Note that the metadata will have to be subsetted as well
# Phylum
rep200FungiPhylumShared_8cancer_Nonzero <- rep200FungiPhylumShared_8cancer[rowSums(rep200FungiPhylumShared_8cancer) != 0,]
metaQiitaCombined_Nonzero_8cancer_PhylumShared <- droplevels(metaQiitaCombined_Nonzero_8cancer_shared[rowSums(rep200FungiPhylumShared_8cancer) != 0,])
# Class
rep200FungiClassShared_8cancer_Nonzero <- rep200FungiClassShared_8cancer[rowSums(rep200FungiClassShared_8cancer) != 0,]
metaQiitaCombined_Nonzero_8cancer_ClassShared <- droplevels(metaQiitaCombined_Nonzero_8cancer_shared[rowSums(rep200FungiClassShared_8cancer) != 0,])
# Order
rep200FungiOrderShared_8cancer_Nonzero <- rep200FungiOrderShared_8cancer[rowSums(rep200FungiOrderShared_8cancer) != 0,]
metaQiitaCombined_Nonzero_8cancer_OrderShared <- droplevels(metaQiitaCombined_Nonzero_8cancer_shared[rowSums(rep200FungiOrderShared_8cancer) != 0,])
# Family
rep200FungiFamilyShared_8cancer_Nonzero <- rep200FungiFamilyShared_8cancer[rowSums(rep200FungiFamilyShared_8cancer) != 0,]
metaQiitaCombined_Nonzero_8cancer_FamilyShared <- droplevels(metaQiitaCombined_Nonzero_8cancer_shared[rowSums(rep200FungiFamilyShared_8cancer) != 0,])
# Genus
rep200FungiGenusShared_8cancer_Nonzero <- rep200FungiGenusShared_8cancer[rowSums(rep200FungiGenusShared_8cancer) != 0,]
metaQiitaCombined_Nonzero_8cancer_GenusShared <- droplevels(metaQiitaCombined_Nonzero_8cancer_shared[rowSums(rep200FungiGenusShared_8cancer) != 0,])
# Species
rep200FungiSpeciesShared_8cancer_Nonzero <- rep200FungiSpeciesShared_8cancer[rowSums(rep200FungiSpeciesShared_8cancer) != 0,]
metaQiitaCombined_Nonzero_8cancer_SpeciesShared <- droplevels(metaQiitaCombined_Nonzero_8cancer_shared[rowSums(rep200FungiSpeciesShared_8cancer) != 0,])

save(rep200FungiPhylumShared_8cancer_Nonzero,
     metaQiitaCombined_Nonzero_8cancer_PhylumShared,
     rep200FungiClassShared_8cancer_Nonzero,
     metaQiitaCombined_Nonzero_8cancer_ClassShared,
     rep200FungiOrderShared_8cancer_Nonzero,
     metaQiitaCombined_Nonzero_8cancer_OrderShared,
     rep200FungiFamilyShared_8cancer_Nonzero,
     metaQiitaCombined_Nonzero_8cancer_FamilyShared,
     rep200FungiGenusShared_8cancer_Nonzero,
     metaQiitaCombined_Nonzero_8cancer_GenusShared,
     rep200FungiSpeciesShared_8cancer_Nonzero,
     metaQiitaCombined_Nonzero_8cancer_SpeciesShared,
     file = "Interim_data/data_tcga_8_cancers_features_matched_to_weizmann_15Sep21.RData")

#-----------------------------------------#
# Intersect presence/absence per CT at genus level
#-----------------------------------------#

source("00-Functions.R") # for the following functions:
# findPrevGenus() and findPrevSpecies()

seqCentersAll <- names(table(metaQiitaCombined_Nonzero_8cancer_shared$data_submitting_center_label))
metaQiitaCombined_Nonzero_8cancer_WGS <- metaQiitaCombined_Nonzero_8cancer_shared %>% filter(experimental_strategy == "WGS") %>% droplevels()
seqCentersWGS <- names(table(metaQiitaCombined_Nonzero_8cancer_WGS$data_submitting_center_label))
metaQiitaCombined_Nonzero_8cancer_RNA <- metaQiitaCombined_Nonzero_8cancer_shared %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()
seqCentersRNA <- names(table(metaQiitaCombined_Nonzero_8cancer_RNA$data_submitting_center_label))

#---------------------------WGS comparisons at genus level---------------------------#
source("00-Functions.R")
# bone: 11/14 (Weizmann) | 11/37 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Bone Cancer", taxaFlag = FALSE, 
              seqCenter = seqCentersWGS, weizmannIntersectFlag=TRUE, weizmannTissue="bone")
# breast: weizmann: 37/38 (Weizmann) | 37/51 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Breast Cancer", taxaFlag = FALSE, 
              seqCenter = seqCentersWGS, weizmannIntersectFlag=TRUE, weizmannTissue="breast")
# colon: weizmann:  10/10 (Weizmann) | 10/52 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Colorectal Cancer", taxaFlag = FALSE, 
              seqCenter = seqCentersWGS, weizmannIntersectFlag=TRUE, weizmannTissue="colon")
# gbm: weizmann: 10/10 (Weizmann) | 10/54 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Glioblastoma", taxaFlag = FALSE, 
              seqCenter = seqCentersWGS, weizmannIntersectFlag=TRUE, weizmannTissue="gbm")
# lung: weizmann: 28/28 (Weizmann) | 28/54 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Lung Cancer", taxaFlag = FALSE, 
              seqCenter = seqCentersWGS, weizmannIntersectFlag=TRUE, weizmannTissue="lung")
# melanoma: 10/15 (Weizmann) | 10/29 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Melanoma", taxaFlag = FALSE, 
              seqCenter = seqCentersWGS, weizmannIntersectFlag=TRUE, weizmannTissue="melanoma")
# ovary: weizmann: 21/24 (Weizmann) | 21/44 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Ovarian Cancer", taxaFlag = FALSE, 
              seqCenter = seqCentersWGS, weizmannIntersectFlag=TRUE, weizmannTissue="ovary")
# Note that pancreas is only available via RNA-Seq

#---------------------------RNA comparisons at genus level---------------------------#
# Note that bone cancer / sarcoma is only available via WGS
# breast: weizmann: 24/38 (Weizmann) | 24/30 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Breast Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersRNA, weizmannIntersectFlag=TRUE, weizmannTissue="breast")
# colon: weizmann:  8/10 (Weizmann) | 8/26 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Colorectal Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersRNA, weizmannIntersectFlag=TRUE, weizmannTissue="colon")
# gbm: weizmann: 10/10 (Weizmann) | 10/54 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Glioblastoma", taxaFlag = TRUE, 
              seqCenter = seqCentersRNA, weizmannIntersectFlag=TRUE, weizmannTissue="gbm")
# lung: weizmann: 27/28 (Weizmann) | 27/50 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Lung Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersRNA, weizmannIntersectFlag=TRUE, weizmannTissue="lung")
# melanoma: weizmann: 9/15 (Weizmann) | 9/19 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Melanoma", taxaFlag = TRUE, 
              seqCenter = seqCentersRNA, weizmannIntersectFlag=TRUE, weizmannTissue="melanoma")
# ovary: weizmann: 24/24 (Weizmann) | 24/54 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Ovarian Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersRNA, weizmannIntersectFlag=TRUE, weizmannTissue="ovary")
# pancreas: weizmann: 8/13 (Weizmann) | 8/22 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Pancreatic Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersRNA, weizmannIntersectFlag=TRUE, weizmannTissue="pancreas")

#---------------------------WGS+RNA comparisons at genus level---------------------------#
# bone: weizmann: 11/14 (Weizmann) | 11/37 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Bone Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="bone")
# breast: weizmann: 37/38 (Weizmann) | 37/51 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Breast Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="breast")
# colon: weizmann:  10/10 (Weizmann) | 10/52 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Colorectal Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="colon")
# gbm: weizmann: 10/10 (Weizmann) | 10/54 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Glioblastoma", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="gbm")
# lung: weizmann: 28/28 (Weizmann) | 28/54 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Lung Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="lung")
# melanoma: weizmann: 13/15 (Weizmann) | 13/34 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Melanoma", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="melanoma")
# ovary: weizmann: 24/24 (Weizmann) | 24/54 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Ovarian Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="ovary")
# pancreas: weizmann: 8/13 (Weizmann) | 8/22 (TCGA)
findPrevGenus(countData=rep200FungiGenusShared_8cancer, cancerType = "Pancreatic Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="pancreas")

#---------------------------Plot comparisons using genus---------------------------#
# Create data frame with data
# NOTE: Each pair of numerical entries form a *stacked* barplot, so the first number
# reflects the number of taxa *only* found in the Weizmann data while the second number
# reflect the number of taxa *shared* between the Weizmann data and TCGA.
# ALSO NOTE: The limit of these intersections is defined by the overlap between the rep200 database
# and the Weizmann data. This is how the upper bound of the Weizmann data taxa has been calculated,
# since the intersection can be no more than the maximum intersection with all of rep200.
wgsOverlapGenus <- data.frame(dataset=rep(c("Weizmann only", "Found in Weizmann and TCGA"), 7),
                              cancer_type=rep(c("Bone", "Breast", "Colorectal", "Glioblastoma","Lung","Melanoma","Ovarian"),each=2),
                              genus = c(3, 11, 1, 37, 0, 10, 0, 10, 0, 28, 5, 10, 3, 21))
ggbarplot(wgsOverlapGenus, "cancer_type", "genus", 
          fill = "dataset", color = "dataset", palette = "Paired", ylim=c(0,50),
          xlab = "Cancer/tissue type (tumor and NAT combined when available)", ylab = "Genus count",
          title = "Overlapping fungal genera between TCGA WGS and Weizmann datasets\n(note: intersection among features found in rep200 database)",
          label = TRUE, lab.col = "black", lab.pos = "out") +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank()) +
  ggsave(file = "Figures/Supplementary_Figures/wgs_overlap_rep200_reduced_feature_space_genus.png",
         dpi = "retina", width = 9, height = 4, units = "in")

rnaOverlapGenus <- data.frame(dataset=rep(c("Weizmann only", "Found in Weizmann and TCGA"), 7),
                              cancer_type=rep(c("Breast", "Colorectal", "Glioblastoma","Lung","Melanoma","Ovarian","Pancreatic"),each=2),
                              genus = c(14, 24, 2, 8, 0, 10, 1, 28, 6, 9, 0, 24, 5, 8))
ggbarplot(rnaOverlapGenus, "cancer_type", "genus", 
          fill = "dataset", color = "dataset", palette = "Paired", ylim=c(0,50),
          xlab = "Cancer/tissue type (tumor and NAT combined when available)", ylab = "Genus count",
          title = "Overlapping fungal genera between TCGA RNA-Seq and Weizmann datasets\n(note: intersection among features found in rep200 database)",
          label = TRUE, lab.col = "black", lab.pos = "out") +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank()) +
  ggsave(file = "Figures/Supplementary_Figures/rnaseq_overlap_rep200_reduced_feature_space_genus.png",
         dpi = "retina", width = 9, height = 4, units = "in")

allOverlapGenus <- data.frame(dataset=rep(c("Weizmann only", "Found in Weizmann and TCGA"), 8),
                              cancer_type=rep(c("Bone","Breast", "Colorectal", "Glioblastoma","Lung","Melanoma","Ovarian","Pancreatic"),each=2),
                              genus = c(3, 11, 1, 37, 0, 10, 0, 10, 0, 28, 2, 13, 0, 24, 5, 8))
ggbarplot(allOverlapGenus, "cancer_type", "genus", 
          fill = "dataset", color = "dataset", palette = "Paired", ylim=c(0,50),
          xlab = "Cancer/tissue type (tumor and NAT combined when available)", ylab = "Genus count",
          title = "Overlapping fungal genera between TCGA WGS and RNA-Seq and Weizmann datasets\n(note: intersection among features found in rep200 database)",
          label = TRUE, lab.col = "black", lab.pos = "out") +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank()) +
  ggsave(file = "Figures/Supplementary_Figures/wgs_and_rnaseq_overlap_rep200_reduced_feature_space_genus.png",
         dpi = "retina", width = 9, height = 4, units = "in")

#-----------------------------------------#
# Intersect presence/absence per CT at species level
#-----------------------------------------#

source("00-Functions.R") # for the following functions:
# findPrevGenus() and findPrevSpecies()

#---------------------------WGS comparisons at species level---------------------------#
source("00-Functions.R")
# bone: 9/10 (Weizmann) | 9/27 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Bone Cancer", taxaFlag = FALSE, 
              seqCenter = seqCentersWGS, weizmannIntersectFlag=TRUE, weizmannTissue="bone")
# breast: weizmann: 19/20 (Weizmann) | 19/32 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Breast Cancer", taxaFlag = FALSE, 
              seqCenter = seqCentersWGS, weizmannIntersectFlag=TRUE, weizmannTissue="breast")
# colon: weizmann:  7/8 (Weizmann) | 7/30 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Colorectal Cancer", taxaFlag = FALSE, 
              seqCenter = seqCentersWGS, weizmannIntersectFlag=TRUE, weizmannTissue="colon")
# gbm: weizmann: 6/6 (Weizmann) | 6/34 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Glioblastoma", taxaFlag = FALSE, 
              seqCenter = seqCentersWGS, weizmannIntersectFlag=TRUE, weizmannTissue="gbm")
# lung: weizmann: 21/21 (Weizmann) | 21/34 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Lung Cancer", taxaFlag = FALSE, 
              seqCenter = seqCentersWGS, weizmannIntersectFlag=TRUE, weizmannTissue="lung")
# melanoma: 6/11 (Weizmann) | 6/15 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Melanoma", taxaFlag = FALSE, 
              seqCenter = seqCentersWGS, weizmannIntersectFlag=TRUE, weizmannTissue="melanoma")
# ovary: weizmann: 10/11 (Weizmann) | 10/27 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Ovarian Cancer", taxaFlag = FALSE, 
              seqCenter = seqCentersWGS, weizmannIntersectFlag=TRUE, weizmannTissue="ovary")
# Note that pancreas is only available via RNA-Seq

#---------------------------RNA comparisons at species level---------------------------#
# Note that bone cancer / sarcoma is only available via WGS
# breast: weizmann: 7/20 (Weizmann) | 7/16 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Breast Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersRNA, weizmannIntersectFlag=TRUE, weizmannTissue="breast")
# colon: weizmann:  3/8 (Weizmann) | 3/11 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Colorectal Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersRNA, weizmannIntersectFlag=TRUE, weizmannTissue="colon")
# gbm: weizmann: 6/6 (Weizmann) | 6/34 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Glioblastoma", taxaFlag = TRUE, 
              seqCenter = seqCentersRNA, weizmannIntersectFlag=TRUE, weizmannTissue="gbm")
# lung: weizmann: 21/21 (Weizmann) | 21/34 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Lung Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersRNA, weizmannIntersectFlag=TRUE, weizmannTissue="lung")
# melanoma: weizmann: 5/11 (Weizmann) | 5/11 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Melanoma", taxaFlag = TRUE, 
              seqCenter = seqCentersRNA, weizmannIntersectFlag=TRUE, weizmannTissue="melanoma")
# ovary: weizmann: 11/11 (Weizmann) | 11/34 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Ovarian Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersRNA, weizmannIntersectFlag=TRUE, weizmannTissue="ovary")
# pancreas: weizmann: 2/7 (Weizmann) | 2/9 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Pancreatic Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersRNA, weizmannIntersectFlag=TRUE, weizmannTissue="pancreas")

#---------------------------WGS+RNA comparisons at species level---------------------------#
# bone: weizmann: 9/10 (Weizmann) | 9/27 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Bone Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="bone")
# breast: weizmann: 19/20 (Weizmann) | 19/32 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Breast Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="breast")
# colon: weizmann:  7/8 (Weizmann) | 7/30 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Colorectal Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="colon")
# gbm: weizmann: 6/6 (Weizmann) | 6/34 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Glioblastoma", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="gbm")
# lung: weizmann: 21/21 (Weizmann) | 21/34 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Lung Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="lung")
# melanoma: weizmann: 7/11 (Weizmann) | 7/19 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Melanoma", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="melanoma")
# ovary: weizmann: 11/11 (Weizmann) | 11/34 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Ovarian Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="ovary")
# pancreas: weizmann: 2/7 (Weizmann) | 2/9 (TCGA)
findPrevSpecies(countData=rep200FungiSpeciesShared_8cancer, cancerType = "Pancreatic Cancer", taxaFlag = TRUE, 
              seqCenter = seqCentersAll, weizmannIntersectFlag=TRUE, weizmannTissue="pancreas")

#---------------------------Plot comparisons using species---------------------------#
# Create data frame with data
# NOTE: Each pair of numerical entries form a *stacked* barplot, so the first number
# reflects the number of taxa *only* found in the Weizmann data while the second number
# reflect the number of taxa *shared* between the Weizmann data and TCGA.
# ALSO NOTE: The limit of these intersections is defined by the overlap between the rep200 database
# and the Weizmann data. This is how the upper bound of the Weizmann data taxa has been calculated,
# since the intersection can be no more than the maximum intersection with all of rep200.
wgsOverlapSpecies <- data.frame(dataset=rep(c("Weizmann only", "Found in Weizmann and TCGA"), 7),
                                cancer_type=rep(c("Bone", "Breast", "Colorectal", "Glioblastoma","Lung","Melanoma","Ovarian"),each=2),
                                species = c(1, 9, 1, 19, 1, 7, 0, 6, 0, 21, 5, 6, 1, 10))
ggbarplot(wgsOverlapSpecies, "cancer_type", "species", 
          fill = "dataset", color = "dataset", palette = "Paired", ylim=c(0,25),
          xlab = "Cancer/tissue type (tumor and NAT combined when available)", ylab = "Species count",
          title = "Overlapping fungal species between TCGA WGS and Weizmann datasets\n(note: intersection among features found in rep200 database)",
          label = TRUE, lab.col = "black", lab.pos = "out") +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank()) +
  ggsave(file = "Figures/Supplementary_Figures/wgs_overlap_rep200_reduced_feature_space_species.png",
         dpi = "retina", width = 9, height = 4, units = "in")

rnaOverlapSpecies <- data.frame(dataset=rep(c("Weizmann only", "Found in Weizmann and TCGA"), 7),
                                cancer_type=rep(c("Breast", "Colorectal", "Glioblastoma","Lung","Melanoma","Ovarian","Pancreatic"),each=2),
                                species = c(13, 7, 5, 3, 0, 6, 0, 21, 6, 5, 0, 11, 5, 2))
ggbarplot(rnaOverlapSpecies, "cancer_type", "species", 
          fill = "dataset", color = "dataset", palette = "Paired", ylim=c(0,25),
          xlab = "Cancer/tissue type (tumor and NAT combined when available)", ylab = "Species count",
          title = "Overlapping fungal species between TCGA RNA-Seq and Weizmann datasets\n(note: intersection among features found in rep200 database)",
          label = TRUE, lab.col = "black", lab.pos = "out") +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank()) +
  ggsave(file = "Figures/Supplementary_Figures/rnaseq_overlap_rep200_reduced_feature_space_species.png",
         dpi = "retina", width = 9, height = 4, units = "in")

allOverlapSpecies <- data.frame(dataset=rep(c("Weizmann only", "Found in Weizmann and TCGA"), 8),
                                cancer_type=rep(c("Bone","Breast", "Colorectal", "Glioblastoma","Lung","Melanoma","Ovarian","Pancreatic"),each=2),
                                species = c(1, 9, 1, 19, 1, 7, 0, 6, 0, 21, 4, 7, 0, 11, 5, 2))
ggbarplot(allOverlapSpecies, "cancer_type", "species", 
          fill = "dataset", color = "dataset", palette = "Paired", ylim=c(0,25),
          xlab = "Cancer/tissue type (tumor and NAT combined when available)", ylab = "Species count",
          title = "Overlapping fungal species between TCGA WGS and RNA-Seq and Weizmann datasets\n(note: intersection among features found in rep200 database)",
          label = TRUE, lab.col = "black", lab.pos = "out") +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank()) +
  ggsave(file = "Figures/Supplementary_Figures/wgs_and_rnaseq_overlap_rep200_reduced_feature_space_species.png",
         dpi = "retina", width = 9, height = 4, units = "in")

