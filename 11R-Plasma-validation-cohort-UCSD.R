#-----------------------------------------------------------------------------
# 11R-Plasma-validation-cohort-UCSD.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Import rep200 data for plasma validation cohort #1 (Poore & Kopylova et al. 2020. Nature)
# - Decontaminate using frequency and prevalence modes of decontam
# - Batch correct for known differences in age and gender between groups
# - Perform machine learning between cancer vs. healthy groups
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
# https://qiita.ucsd.edu/analysis/description/47212/

## Import Qiita metadata
metaPoore <- read.csv("Input_data/qiita_metadata_plasma_validation_cohort_UCSD.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(metaPoore)[1] <- "sample_name"
# Load consolidated metadata from Poore & Kopylova et al. 2020. Nature
metaUCSD <- read.csv("Input_data/poore_et_al_plasma_validation_metadata_consolidated.txt", sep = "\t", stringsAsFactors = FALSE, row.names = 1)

sum(metaPoore$tube_id %in% metaUCSD$tube_id) # 169
sum(metaUCSD$tube_id %in% metaPoore$tube_id) # 169

## Intersect tube IDs and subset
metaPooreFilt <- droplevels(metaPoore[metaPoore$tube_id %in% metaUCSD$tube_id,])
metaUCSD$sampleID <- rownames(metaUCSD)

metaUCSDJoined <- left_join(metaUCSD,
                            metaPooreFilt[,c("sample_name","tube_id")],
                            by = "tube_id") %>% column_to_rownames("sample_name")

#-------------------------Import count data-------------------------#

ucsd_rep200Data_BIOM <- read_biom(biom_file = "Input_data/Qiita_results/rep200_OGU_plasma_validation_cohort_UCSD.biom")
ucsd_rep200Data <- t(as(biom_data(ucsd_rep200Data_BIOM), "matrix"))

sum(rownames(ucsd_rep200Data) %in% rownames(metaUCSDJoined)) # 169
sum(rownames(metaUCSDJoined)%in% rownames(ucsd_rep200Data)) # 169

ucsd_rep200Data_Filt <- ucsd_rep200Data[rownames(metaUCSDJoined),]
dim(ucsd_rep200Data_Filt) # 169 7777

## Extract only fungal and bacterial features
ucsd_rep200Data_Filt_Fungi <- ucsd_rep200Data_Filt[,colnames(ucsd_rep200Data_Filt) %in% fungiOGUs]
ucsd_rep200Data_Filt_Bacteria <- ucsd_rep200Data_Filt[,colnames(ucsd_rep200Data_Filt) %in% bacteriaOGUs]

## Remove samples with 0 counts after filtering to Fungi: {"12667.X2051123", "12691.PC22", "12691.PC47"}
zeroSumFungiSamples <- names(which(rowSums(ucsd_rep200Data_Filt_Fungi) == 0))
ucsd_rep200Data_Filt_Fungi_Nonzero <- ucsd_rep200Data_Filt_Fungi[!(rownames(ucsd_rep200Data_Filt_Fungi) %in% zeroSumFungiSamples),]
metaUCSDJoined_Fungi_Nonzero <- droplevels(metaUCSDJoined[!(rownames(metaUCSDJoined) %in% zeroSumFungiSamples),])
all(rownames(metaUCSDJoined_Fungi_Nonzero) == rownames(ucsd_rep200Data_Filt_Fungi_Nonzero)) # TRUE

## NOTE: No bacteria samples have 0 sum. Nonetheless, we will create a matched table to directly compare fungi vs. bacteria
ucsd_rep200Data_Filt_Bacteria_Matched <- ucsd_rep200Data_Filt_Bacteria[!(rownames(ucsd_rep200Data_Filt_Bacteria) %in% zeroSumFungiSamples),]
metaUCSDJoined_Bacteria_Matched <- droplevels(metaUCSDJoined[!(rownames(metaUCSDJoined) %in% zeroSumFungiSamples),])
all(rownames(metaUCSDJoined_Bacteria_Matched) == rownames(metaUCSDJoined_Fungi_Nonzero)) # TRUE

## NOTE: No full rep200 samples have 0 sum. Nonetheless, we will create a matched table to directly compare to fungi and bacteria
ucsd_rep200Data_Filt_Matched <- ucsd_rep200Data_Filt[!(rownames(ucsd_rep200Data_Filt) %in% zeroSumFungiSamples),]
metaUCSDJoined_FullRep200_Matched <- droplevels(metaUCSDJoined[!(rownames(metaUCSDJoined) %in% zeroSumFungiSamples),])
all(rownames(metaUCSDJoined_FullRep200_Matched) == rownames(metaUCSDJoined_Fungi_Nonzero)) # TRUE

#----------------------------------------------------------------------------------------------#
# Summarize bacteria data at each taxonomic level
#----------------------------------------------------------------------------------------------#
# Modify bacterial taxa table
rep200TaxSplit_Bacteria_Formatted <- apply(rep200TaxSplit_Bacteria, 2, function(x) gsub("^[k|p|c|o|f|g|s]__","",x))
rep200TaxSplit_Bacteria_Formatted[rep200TaxSplit_Bacteria_Formatted == ""] <- "other"

# Build phyloseq object
psUCSDBacteria <- phyloseq(otu_table(ucsd_rep200Data_Filt_Bacteria_Matched, taxa_are_rows = FALSE), 
                        tax_table(as.matrix(rep200TaxSplit_Bacteria_Formatted)), 
                        sample_data(metaUCSDJoined_Bacteria_Matched))

## Aggregate counts
psUCSDBacteria_phylum = aggregate_taxa(psUCSDBacteria, "Phylum")
psUCSDBacteria_class = aggregate_taxa(psUCSDBacteria, "Class")
psUCSDBacteria_order = aggregate_taxa(psUCSDBacteria, "Order")
psUCSDBacteria_family = aggregate_taxa(psUCSDBacteria, "Family")
psUCSDBacteria_genus = aggregate_taxa(psUCSDBacteria, "Genus")
psUCSDBacteria_species = aggregate_taxa(psUCSDBacteria, "Species")
# colnames(psUCSDBacteria_species) <- rep200TaxSplit_Bacteria_Formatted[colnames(psUCSDBacteria_species),"Species"]

## Create data.frames of summarized data
ucsdRep200BacteriaPhylum <- data.frame(t(otu_table(psUCSDBacteria_phylum)))
ucsdRep200BacteriaClass <- data.frame(t(otu_table(psUCSDBacteria_class)))
ucsdRep200BacteriaOrder <- data.frame(t(otu_table(psUCSDBacteria_order)))
ucsdRep200BacteriaFamily <- data.frame(t(otu_table(psUCSDBacteria_family)))
ucsdRep200BacteriaGenus <- data.frame(t(otu_table(psUCSDBacteria_genus)))
ucsdRep200BacteriaSpecies <- data.frame(t(otu_table(psUCSDBacteria_species))) # data.frame(psUCSDBacteria_species)
ucsdRep200BacteriaSpecies[1:3,1:3]

## Subset taxa to those shared with the Weizmann data
load("Interim_data/shared_bacterial_features_at_each_taxa_level_29Mar22.RData")

psUCSDBacteria_phylum_shared <- subset_taxa(psUCSDBacteria_phylum, Phylum %in% sharedPhylumBacteria$intersectedTaxa)
psUCSDBacteria_class_shared <- subset_taxa(psUCSDBacteria_class, Class %in% sharedClassBacteria$intersectedTaxa)
psUCSDBacteria_order_shared <- subset_taxa(psUCSDBacteria_order, Order %in% sharedOrderBacteria$intersectedTaxa)
psUCSDBacteria_family_shared <- subset_taxa(psUCSDBacteria_family, Family %in% sharedFamilyBacteria$intersectedTaxa)
psUCSDBacteria_genus_shared <- subset_taxa(psUCSDBacteria_genus, Genus %in% sharedGenusBacteria$intersectedTaxa)
psUCSDBacteria_species_shared <- subset_taxa(psUCSDBacteria_species, Species %in% sharedSpeciesBacteria$intersectedTaxa)

# psUCSDBacteria_species_shared <- psUCSDBacteria_species[,colnames(psUCSDBacteria_species) %in% sharedSpeciesBacteria$intersectedTaxa]

## Create data.frames of summarized data with feature intersection
ucsdRep200BacteriaPhylumShared <- data.frame(t(otu_table(psUCSDBacteria_phylum_shared)))
ucsdRep200BacteriaClassShared <- data.frame(t(otu_table(psUCSDBacteria_class_shared)))
ucsdRep200BacteriaOrderShared <- data.frame(t(otu_table(psUCSDBacteria_order_shared)))
ucsdRep200BacteriaFamilyShared <- data.frame(t(otu_table(psUCSDBacteria_family_shared)))
ucsdRep200BacteriaGenusShared <- data.frame(t(otu_table(psUCSDBacteria_genus_shared)))
ucsdRep200BacteriaSpeciesShared <- data.frame(t(otu_table(psUCSDBacteria_species_shared)))
ucsdRep200BacteriaSpeciesShared[1:3,1:3]

#----------------------------------------------------------------------------------------------#
# Summarize fungi data at each taxonomic level and intersect with Weizmann features
#----------------------------------------------------------------------------------------------#

rep200TaxSplit_Fungi_Paired_to_Weizmann <- read.csv("Supporting_data/rep200TaxSplit_Fungi_Paired_To_Weizmann_Final.csv", stringsAsFactors = FALSE, row.names = 1)

# Build phyloseq object
psUCSDFungi <- phyloseq(otu_table(ucsd_rep200Data_Filt_Fungi_Nonzero, taxa_are_rows = FALSE), 
                        tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), 
                        sample_data(metaUCSDJoined_Fungi_Nonzero))

## Aggregate counts
psUCSDFungi_phylum = aggregate_taxa(psUCSDFungi, "phylum")
psUCSDFungi_class = aggregate_taxa(psUCSDFungi, "class")
psUCSDFungi_order = aggregate_taxa(psUCSDFungi, "order")
psUCSDFungi_family = aggregate_taxa(psUCSDFungi, "family")
psUCSDFungi_genus = aggregate_taxa(psUCSDFungi, "genus")
psUCSDFungi_species = ucsd_rep200Data_Filt_Fungi_Nonzero
colnames(psUCSDFungi_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(psUCSDFungi_species),"species"]

## Load shared features with Weizmann
load("Interim_data/shared_fungi_features_at_each_taxa_level_29Mar22.RData")

## Subset taxa to those shared with the Weizmann data
psUCSDFungi_phylum_shared <- subset_taxa(psUCSDFungi_phylum, phylum %in% sharedPhylum)
psUCSDFungi_class_shared <- subset_taxa(psUCSDFungi_class, class %in% sharedClass)
psUCSDFungi_order_shared <- subset_taxa(psUCSDFungi_order, order %in% sharedOrder)
psUCSDFungi_family_shared <- subset_taxa(psUCSDFungi_family, family %in% sharedFamily)
psUCSDFungi_genus_shared <- subset_taxa(psUCSDFungi_genus, genus %in% sharedGenus)
psUCSDFungi_species_shared <- psUCSDFungi_species[,colnames(psUCSDFungi_species) %in% sharedSpecies]

## Create data.frames of summarized data with and without feature intersection
# Without feature intersection
ucsdRep200FungiPhylum <- data.frame(t(otu_table(psUCSDFungi_phylum)))
ucsdRep200FungiClass <- data.frame(t(otu_table(psUCSDFungi_class)))
ucsdRep200FungiOrder <- data.frame(t(otu_table(psUCSDFungi_order)))
ucsdRep200FungiFamily <- data.frame(t(otu_table(psUCSDFungi_family)))
ucsdRep200FungiGenus <- data.frame(t(otu_table(psUCSDFungi_genus)))
ucsdRep200FungiSpecies <- data.frame(psUCSDFungi_species)
ucsdRep200FungiSpecies[1:3,1:3]
# With feature intersection
ucsdRep200FungiPhylumShared <- data.frame(t(otu_table(psUCSDFungi_phylum_shared)))
ucsdRep200FungiClassShared <- data.frame(t(otu_table(psUCSDFungi_class_shared)))
ucsdRep200FungiOrderShared <- data.frame(t(otu_table(psUCSDFungi_order_shared)))
ucsdRep200FungiFamilyShared <- data.frame(t(otu_table(psUCSDFungi_family_shared)))
ucsdRep200FungiGenusShared <- data.frame(t(otu_table(psUCSDFungi_genus_shared)))
ucsdRep200FungiSpeciesShared <- data.frame(psUCSDFungi_species_shared)
ucsdRep200FungiGenusShared[1:3,1:3]

# The following taxa levels had zero sum samples after intersecting features: Family, Genus, Species
# NOTE: All samples without feature intersection had nonzero sample sums
zeroFungiFamilyShared <- unname(which(rowSums(ucsdRep200FungiFamilyShared)==0))
zeroFungiGenusShared <- unname(which(rowSums(ucsdRep200FungiGenusShared)==0))
zeroFungiSpeciesShared <- unname(which(rowSums(ucsdRep200FungiSpeciesShared)==0))

ucsdRep200FungiPhylumShared_Nonzero <- ucsdRep200FungiPhylumShared
ucsdRep200FungiClassShared_Nonzero <- ucsdRep200FungiClassShared
ucsdRep200FungiOrderShared_Nonzero <- ucsdRep200FungiOrderShared
ucsdRep200FungiFamilyShared_Nonzero <- ucsdRep200FungiFamilyShared[-zeroFungiFamilyShared,]
ucsdRep200FungiGenusShared_Nonzero <- ucsdRep200FungiGenusShared[-zeroFungiGenusShared,]
ucsdRep200FungiSpeciesShared_Nonzero <- ucsdRep200FungiSpeciesShared[-zeroFungiSpeciesShared,]

metaUCSDJoined_Fungi_Nonzero_PhylumShared <- metaUCSDJoined_Fungi_Nonzero
metaUCSDJoined_Fungi_Nonzero_ClassShared <- metaUCSDJoined_Fungi_Nonzero
metaUCSDJoined_Fungi_Nonzero_OrderShared <- metaUCSDJoined_Fungi_Nonzero
metaUCSDJoined_Fungi_Nonzero_FamilyShared <- droplevels(metaUCSDJoined_Fungi_Nonzero[-zeroFungiFamilyShared,])
metaUCSDJoined_Fungi_Nonzero_GenusShared <- droplevels(metaUCSDJoined_Fungi_Nonzero[-zeroFungiGenusShared,])
metaUCSDJoined_Fungi_Nonzero_SpeciesShared <- droplevels(metaUCSDJoined_Fungi_Nonzero[-zeroFungiSpeciesShared,])

#-----------------------------------------------#
# Decontaminate using controls and TCGA data
#-----------------------------------------------#
require(decontam)
require(gtools)
require(readr)

# NOTE: Three Qiita projects ((12667 (HIV-free); 12691 (PC); 12692 (LC and SKCM))
# contain the Poore & Kopylova et al. 2020 Nature plasma data even though they were sequenced
# in a single sequencing run. To make it possible for anyone to analyze data from any one Qiita project, 
# the experimental contamination controls were added to each of the three Qiita studies. When aggregating
# the three datasets, as we have done here, only one set of contamination controls from a single project
# needs to be used. For simplicity, we will use controls from Qiita study 12667 (HIV-free) for
# decontamination, although the others would be equally suitable.

## Subset metadata and data using sample_name identifiers containing Qiita project IDs
metaPooreDecontamDedup <- droplevels(metaPoore[-grep("[12692|12691]\\.Control*", metaPoore$sample_name),])
ucsd_rep200Data_decontam_dedup <- ucsd_rep200Data[metaPooreDecontamDedup$sample_name,]
ucsd_rep200Data_decontam_dedup_fungi <- ucsd_rep200Data_decontam_dedup[,colnames(ucsd_rep200Data_decontam_dedup) %in% fungiOGUs]

## Remove samples with 0 fungi counts prior to decontam
# NOTE: This step is different than the one above because all samples are being fed into
# decontam whereas the above steps separated out only the biological plasma samples of interest,
# which we will decontaminate and normalize prior to machine learning
zeroCountSamplesDecontamFungi <- unname(which(rowSums(ucsd_rep200Data_decontam_dedup_fungi)==0))
metaPooreDecontamDedup_Nonzero <- droplevels(metaPooreDecontamDedup[-zeroCountSamplesDecontamFungi,])
ucsd_rep200Data_decontam_dedup_fungi_nonzero <- ucsd_rep200Data_decontam_dedup_fungi[-zeroCountSamplesDecontamFungi,]

#--------------------------Decontaminate using prevalence--------------------------#
require(decontam)
# Prevalence based decontamination using all control samples (extraction blanks and well blanks)
contamdf.prev.fungi.ucsd <- isContaminant(seqtab = as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero), 
                                     method = "prevalence",
                                     neg = grepl("control",metaPooreDecontamDedup_Nonzero$sample_type),
                                     threshold = 0.5)
# **Quoted from the decontam tutorial:**
# "In the prevalence test there is a special value worth knowing, 
# threshold=0.5, that will identify as contaminants all sequences 
# thare are more prevalent in negative controls than in positive samples."
table(contamdf.prev.fungi.ucsd$contaminant) # TRUE=30 | FALSE=197
contamdf.prev.fungi.ucsd[which(contamdf.prev.fungi.ucsd$contaminant),]

notContamSumPrev <- colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero)[,!contamdf.prev.fungi.ucsd$contaminant])
contamSumPrev <- colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero)[,contamdf.prev.fungi.ucsd$contaminant])
sum(contamSumPrev)/sum(colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero))) #--> 0.08865087

#--------------------------Decontaminate using frequency--------------------------#
# First remove samples with DNA concentrations of 0 or decontam will fail
metaPooreDecontamDedup_Nonzero_Freq <- droplevels(metaPooreDecontamDedup_Nonzero[metaPooreDecontamDedup_Nonzero$well_conc != 0,])
ucsd_rep200Data_decontam_dedup_fungi_nonzero_freq <- ucsd_rep200Data_decontam_dedup_fungi_nonzero[metaPooreDecontamDedup_Nonzero_Freq$sample_name,]
dim(ucsd_rep200Data_decontam_dedup_fungi_nonzero_freq) # 346 227

contamdf.freq.fungi.ucsd <- isContaminant(seqtab = as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero_freq), 
                                     method = "frequency",
                                     conc = metaPooreDecontamDedup_Nonzero_Freq$well_conc,
                                     threshold = 0.1)
table(contamdf.freq.fungi.ucsd$contaminant) # TRUE=4 | FALSE=223

notContamSumFreq <- colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero_freq)[,!contamdf.freq.fungi.ucsd$contaminant])
contamSumFreq <- colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero_freq)[,contamdf.freq.fungi.ucsd$contaminant])
sum(contamSumFreq)/sum(colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero_freq))) #--> 0.002012228 (Apr 4 2022)

#--------------------------Combine prevalence and frequency contaminants and calculate % read removal--------------------------#
putativeContaminantsFungiUCSD <- unique(c(rownames(contamdf.prev.fungi.ucsd[which(contamdf.prev.fungi.ucsd$contaminant),]),
                                       rownames(contamdf.freq.fungi.ucsd[which(contamdf.freq.fungi.ucsd$contaminant),])))
length(putativeContaminantsFungiUCSD) # 32

# Calculate total % read count removal using prevalence and frequency decontamination
notContamSumTotal <- colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero)[,!(colnames(ucsd_rep200Data_decontam_dedup_fungi_nonzero) %in% putativeContaminantsFungiUCSD)])
contamSumTotal <- colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero)[,putativeContaminantsFungiUCSD])
sum(contamSumTotal)/sum(colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero))) #--> 0.08948648 (Apr 4 2022)

#--------------------------Combine decontam and TCGA contaminants--------------------------#
load("Interim_data/decontamResultsV2_25Mar22.RData")
putativeContaminantsFungiUCSD_decontamV2Results <- decontamResultsV2[putativeContaminantsFungiUCSD,]
# Remove fungi with unknown human associations predicted to be contaminants
indx2Keep <- which(grepl("Unknown human association|No known human association",putativeContaminantsFungiUCSD_decontamV2Results$reason))
putativeContaminantsFungiUCSD_decontamV2ResultsFilt <- putativeContaminantsFungiUCSD_decontamV2Results[indx2Keep,]

notContamSumTotalWithTCGA <- colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero)[,!(colnames(ucsd_rep200Data_decontam_dedup_fungi_nonzero) %in% rownames(putativeContaminantsFungiUCSD_decontamV2ResultsFilt))])
contamSumTotalWithTCGA <- colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero)[,colnames(ucsd_rep200Data_decontam_dedup_fungi_nonzero) %in% rownames(putativeContaminantsFungiUCSD_decontamV2ResultsFilt)])
sum(contamSumTotalWithTCGA)/sum(colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero))) #--> 0.03783045 (Apr 4 2022)

#--------------------------Remove contaminants from UCSD plasma cohort--------------------------#

ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam <- ucsd_rep200Data_Filt_Fungi_Nonzero[,!(colnames(ucsd_rep200Data_Filt_Fungi_Nonzero) %in% 
                                                                                       rownames(putativeContaminantsFungiUCSD_decontamV2ResultsFilt))]
dim(ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam) # 166 215
# Sanity check
all(rownames(metaUCSDJoined_Fungi_Nonzero) == rownames(ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam)) # TRUE

#--------------------------Subset to topX signature from Hopkins cohort--------------------------#
load("Interim_data/topXHopkinsSig_13Nov21.RData")
ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_TopX <- ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam[,colnames(ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam) %in%
                                                                                                  rownames(topXHopkinsSig)]
dim(ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_TopX) # 166  18 
# Remove 6 zero sum samples after subsetting to top features from Hopkins cohort
ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_TopX_Nonzero <- ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_TopX[-which(rowSums(ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_TopX)==0),]
# Sanity check before subsetting metadata
all(rownames(metaUCSDJoined_Fungi_Nonzero) == rownames(ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_TopX)) # TRUE
metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero <- metaUCSDJoined_Fungi_Nonzero[-which(rowSums(ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_TopX)==0),]

#----------------------------------------------------------------------------------------------#
# Summarize decontaminated fungi data at each taxonomic level
#----------------------------------------------------------------------------------------------#

rep200TaxSplit_Fungi_Paired_to_Weizmann <- read.csv("Supporting_data/rep200TaxSplit_Fungi_Paired_To_Weizmann_Final.csv", stringsAsFactors = FALSE, row.names = 1)

# Build phyloseq object
psUCSDFungiDecontam <- phyloseq(otu_table(ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam, taxa_are_rows = FALSE), 
                        tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), 
                        sample_data(metaUCSDJoined_Fungi_Nonzero))

## Aggregate counts
psUCSDFungiDecontam_phylum = aggregate_taxa(psUCSDFungiDecontam, "phylum")
psUCSDFungiDecontam_class = aggregate_taxa(psUCSDFungiDecontam, "class")
psUCSDFungiDecontam_order = aggregate_taxa(psUCSDFungiDecontam, "order")
psUCSDFungiDecontam_family = aggregate_taxa(psUCSDFungiDecontam, "family")
psUCSDFungiDecontam_genus = aggregate_taxa(psUCSDFungiDecontam, "genus")
psUCSDFungiDecontam_species = ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam
colnames(psUCSDFungiDecontam_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(psUCSDFungiDecontam_species),"species"]

## Create data.frames of summarized data
ucsdRep200FungiDecontamPhylum <- data.frame(t(otu_table(psUCSDFungiDecontam_phylum)))
ucsdRep200FungiDecontamClass <- data.frame(t(otu_table(psUCSDFungiDecontam_class)))
ucsdRep200FungiDecontamOrder <- data.frame(t(otu_table(psUCSDFungiDecontam_order)))
ucsdRep200FungiDecontamFamily <- data.frame(t(otu_table(psUCSDFungiDecontam_family)))
ucsdRep200FungiDecontamGenus <- data.frame(t(otu_table(psUCSDFungiDecontam_genus)))
ucsdRep200FungiDecontamSpecies <- data.frame(psUCSDFungiDecontam_species)
ucsdRep200FungiDecontamGenus[1:3,1:3]

#-----------------------------------------------
#           VSNM Normalization                 #
#-----------------------------------------------

source("00-Functions.R") # for vsnmFunctionUCSD()

# Full rep200 data
snmDataUCSD_FullRep200_OGU <- vsnmFunctionUCSD(qcData = ucsd_rep200Data_Filt_Matched, qcMetadata = metaUCSDJoined_FullRep200_Matched)

# All bacteria data by taxa levels
snmData_ucsdRep200BacteriaPhylum <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaPhylum, qcMetadata = metaUCSDJoined_Bacteria_Matched)
snmData_ucsdRep200BacteriaClass <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaClass, qcMetadata = metaUCSDJoined_Bacteria_Matched)
snmData_ucsdRep200BacteriaOrder <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaOrder, qcMetadata = metaUCSDJoined_Bacteria_Matched)
snmData_ucsdRep200BacteriaFamily <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaFamily, qcMetadata = metaUCSDJoined_Bacteria_Matched)
snmData_ucsdRep200BacteriaGenus <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaGenus, qcMetadata = metaUCSDJoined_Bacteria_Matched)
snmData_ucsdRep200BacteriaSpecies <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaSpecies, qcMetadata = metaUCSDJoined_Bacteria_Matched)
snmDataUCSD_Bacteria_OGU <- vsnmFunctionUCSD(qcData = ucsd_rep200Data_Filt_Bacteria_Matched, qcMetadata = metaUCSDJoined_Bacteria_Matched)

# # Shared bacteria data by taxa levels
# snmData_ucsdRep200BacteriaPhylumShared <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaPhylumShared, qcMetadata = metaUCSDJoined_Bacteria_Matched)
# # snmData_ucsdRep200BacteriaClassShared <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaClassShared, qcMetadata = metaUCSDJoined_Bacteria_Matched)
# # snmData_ucsdRep200BacteriaOrderShared <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaOrderShared, qcMetadata = metaUCSDJoined_Bacteria_Matched)
# snmData_ucsdRep200BacteriaFamilyShared <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaFamilyShared, qcMetadata = metaUCSDJoined_Bacteria_Matched)
# snmData_ucsdRep200BacteriaGenusShared <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaGenusShared, qcMetadata = metaUCSDJoined_Bacteria_Matched)
snmData_ucsdRep200BacteriaSpeciesShared <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaSpeciesShared, qcMetadata = metaUCSDJoined_Bacteria_Matched)
# 
# # Fungi data by taxa levels
snmData_ucsdRep200FungiPhylum <- vsnmFunctionUCSD(qcData = ucsdRep200FungiPhylum, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiClass <- vsnmFunctionUCSD(qcData = ucsdRep200FungiClass, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiOrder <- vsnmFunctionUCSD(qcData = ucsdRep200FungiOrder, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiFamily <- vsnmFunctionUCSD(qcData = ucsdRep200FungiFamily, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiGenus <- vsnmFunctionUCSD(qcData = ucsdRep200FungiGenus, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiSpecies <- vsnmFunctionUCSD(qcData = ucsdRep200FungiSpecies, qcMetadata = metaUCSDJoined_Fungi_Nonzero)

# Decontaminated fungi data
snmDataUCSD_Fungi_Decontam_OGU <- vsnmFunctionUCSD(qcData = ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam, qcMetadata = metaUCSDJoined_Fungi_Nonzero)

# TopX fungi data
snmDataUCSD_Fungi_TopX <- vsnmFunctionUCSD(qcData = ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_TopX_Nonzero, 
                                           qcMetadata = metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero)

# # Decontaminated fungi data by taxa levels
snmData_ucsdRep200FungiDecontamPhylum <- vsnmFunctionUCSD(qcData = ucsdRep200FungiDecontamPhylum, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiDecontamClass <- vsnmFunctionUCSD(qcData = ucsdRep200FungiDecontamClass, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiDecontamOrder <- vsnmFunctionUCSD(qcData = ucsdRep200FungiDecontamOrder, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiDecontamFamily <- vsnmFunctionUCSD(qcData = ucsdRep200FungiDecontamFamily, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiDecontamGenus <- vsnmFunctionUCSD(qcData = ucsdRep200FungiDecontamGenus, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiDecontamSpecies <- vsnmFunctionUCSD(qcData = ucsdRep200FungiDecontamSpecies, qcMetadata = metaUCSDJoined_Fungi_Nonzero)

# Fungi data intersected with Weizmann features
# # NOTE: The following taxa levels did not converge and are not shown: Phylum
# # snmData_ucsdRep200FungiPhylumShared <- vsnmFunctionUCSD(qcData = ucsdRep200FungiPhylumShared_Nonzero, qcMetadata = metaUCSDJoined_Fungi_Nonzero_PhylumShared)
# snmData_ucsdRep200FungiClassShared <- vsnmFunctionUCSD(qcData = ucsdRep200FungiClassShared_Nonzero, qcMetadata = metaUCSDJoined_Fungi_Nonzero_ClassShared)
# snmData_ucsdRep200FungiOrderShared <- vsnmFunctionUCSD(qcData = ucsdRep200FungiOrderShared_Nonzero, qcMetadata = metaUCSDJoined_Fungi_Nonzero_OrderShared)
# snmData_ucsdRep200FungiFamilyShared <- vsnmFunctionUCSD(qcData = ucsdRep200FungiFamilyShared_Nonzero, qcMetadata = metaUCSDJoined_Fungi_Nonzero_FamilyShared)
# snmData_ucsdRep200FungiGenusShared <- vsnmFunctionUCSD(qcData = ucsdRep200FungiGenusShared_Nonzero, qcMetadata = metaUCSDJoined_Fungi_Nonzero_GenusShared)
snmData_ucsdRep200FungiSpeciesShared <- vsnmFunctionUCSD(qcData = ucsdRep200FungiSpeciesShared_Nonzero, qcMetadata = metaUCSDJoined_Fungi_Nonzero_SpeciesShared)

#----------------------------------------------------------------------------------#
# Machine learning: Aggregated cancer vs. Healthy @ OGU/species level
# NOTE: For fungal data, there is only 1 OGU per species
#----------------------------------------------------------------------------------#

source("00-Functions.R") # for ml1VsAllUCSD() function
metaUCSDJoined_Fungi_Nonzero %>% count(HvsC)
# HvsC  n
# 1  Cancer 98
# 2 Control 68

#---------------------------------Full data---------------------------------#
ml_HvsC_Fullrep200_OGU <- ml1VsAllUCSD(metaData = metaUCSDJoined_FullRep200_Matched,
                                      snmData = snmDataUCSD_FullRep200_OGU,
                                      col2Predict = "HvsC",
                                      dataString = "full_rep200",
                                      varImpFlag = FALSE)

#---------------------------------Bacteria---------------------------------#
ml_HvsC_Bacteria_OGU <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                   snmData = snmDataUCSD_Bacteria_OGU,
                                   col2Predict = "HvsC",
                                   dataString = "bacteria_only",
                                   varImpFlag = FALSE)

#---------------------------------Fungi---------------------------------#
ml_HvsC_Fungi_Species <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                           snmData = snmData_ucsdRep200FungiSpecies,
                                           col2Predict = "HvsC",
                                           dataString = "fungi_only",
                                           varImpFlag = FALSE)

#---------------------------------Decontam Fungi---------------------------------#
ml_HvsC_Decontam_Fungi_OGU <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                   snmData = snmDataUCSD_Fungi_Decontam_OGU,
                                   col2Predict = "HvsC",
                                   dataString = "fungi_decontam",
                                   varImpFlag = FALSE)

#---------------------------------Intersected fungi---------------------------------#
ml_HvsC_Shared_Fungi_Species <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                   snmData = snmData_ucsdRep200FungiSpeciesShared,
                                   col2Predict = "HvsC",
                                   dataString = "fungi_intersected_with_Weizmann",
                                   varImpFlag = FALSE)

#---------------------------------Intersected bacteria---------------------------------#
ml_HvsC_Shared_Bacteria_Species <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                             snmData = snmData_ucsdRep200BacteriaSpeciesShared,
                                             col2Predict = "HvsC",
                                             dataString = "bacteria_intersected_with_Weizmann",
                                             varImpFlag = FALSE)

#---------------------------------Intersected fungi+bacteria---------------------------------#

ml_HvsC_Shared_FungiBacteria_Species <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                             snmData = cbind(snmData_ucsdRep200FungiSpeciesShared,
                                                             snmData_ucsdRep200BacteriaSpeciesShared[rownames(snmData_ucsdRep200FungiSpeciesShared),]),
                                             col2Predict = "HvsC",
                                             dataString = "fungi_and_bacteria_intersected_with_Weizmann",
                                             varImpFlag = FALSE)

#---------------------------------TopX fungi---------------------------------#

ml_HvsC_TopX <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero,
                             snmData = snmDataUCSD_Fungi_TopX,
                             col2Predict = "HvsC",
                             dataString = "topX_fungi",
                             varImpFlag = FALSE)

#---------------------------------Combine results---------------------------------#
ucsdCancerVsHealthyResults <- rbind(ml_HvsC_Fullrep200_OGU$rep_perf,
                                    ml_HvsC_Bacteria_OGU$rep_perf,
                                    ml_HvsC_Fungi_Species$rep_perf,
                                    ml_HvsC_Decontam_Fungi_OGU$rep_perf,
                                    ml_HvsC_Shared_Fungi_Species$rep_perf,
                                    ml_HvsC_Shared_Bacteria_Species$rep_perf,
                                    ml_HvsC_Shared_FungiBacteria_Species$rep_perf,
                                    ml_HvsC_TopX$rep_perf)

colnames(ucsdCancerVsHealthyResults)[1:2] <- c("AUROC","AUPR")
ucsdCancerVsHealthyResults$nullAUROC <- 0.5
ucsdCancerVsHealthyResults$nullAUPR <- ucsdCancerVsHealthyResults$positiveClassSize/
  (ucsdCancerVsHealthyResults$positiveClassSize+ucsdCancerVsHealthyResults$negativeClassSize)

ucsdCancerVsHealthyResults$dataString[ucsdCancerVsHealthyResults$dataString=="full_rep200"] <- "Full multikingdom database (rep200)"
ucsdCancerVsHealthyResults$dataString[ucsdCancerVsHealthyResults$dataString=="bacteria_only"] <- "Bacteria all (Species)"
ucsdCancerVsHealthyResults$dataString[ucsdCancerVsHealthyResults$dataString=="fungi_only"] <- "Fungi all (Species)"
ucsdCancerVsHealthyResults$dataString[ucsdCancerVsHealthyResults$dataString=="fungi_decontam"] <- "Fungi decontaminated (Species)"
ucsdCancerVsHealthyResults$dataString[ucsdCancerVsHealthyResults$dataString=="fungi_intersected_with_Weizmann"] <- "Fungi ∩ WIS (Species)"
ucsdCancerVsHealthyResults$dataString[ucsdCancerVsHealthyResults$dataString=="bacteria_intersected_with_Weizmann"] <- "Bacteria ∩ WIS (Species)"
ucsdCancerVsHealthyResults$dataString[ucsdCancerVsHealthyResults$dataString=="fungi_and_bacteria_intersected_with_Weizmann"] <- "Fungi+bacteria ∩ WIS (Species)"
ucsdCancerVsHealthyResults$dataString[ucsdCancerVsHealthyResults$dataString=="topX_fungi"] <- "Top 20 fungi in Hopkins cohort (Species)"

ucsdCancerVsHealthyResults$dataString <- factor(ucsdCancerVsHealthyResults$dataString, levels = c("Full multikingdom database (rep200)",
                                                                                                  "Fungi+bacteria ∩ WIS (Species)",
                                                                                                  "Bacteria all (Species)",
                                                                                                  "Fungi all (Species)",
                                                                                                  "Fungi decontaminated (Species)",
                                                                                                  "Fungi ∩ WIS (Species)",
                                                                                                  "Bacteria ∩ WIS (Species)",
                                                                                                  "Top 20 fungi in Hopkins cohort (Species)"))

source("Supporting_scripts/S00-SummarySE.R")
# Facet by AUROC vs AUPR; color by full dataset vs fungi (dataString)
ucsdCancerVsHealthyResults %>%
  filter(grepl("Full|\\+|decontam|Top",dataString)) %>%
  # filter(grepl("\\+|Bacteria ∩",dataString)) %>%
  reshape2::melt(id.vars = c("rep","dataString","positiveClassSize","negativeClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","dataString","positiveClassSize","negativeClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=dataString)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid", position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + 
  ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + ggtitle("UCSD plasma validation cohort:\nCancer vs. healthy") +
  theme(plot.title = element_text(hjust = 0), legend.position = "right") + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "Cancer vs. Healthy") +
  rotate_x_text(0) + scale_color_nejm(name = "Dataset") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Main_Figures/ucsd_all_cancer_vs_healthy_full_WIS_fungi_TopX_4Apr22.svg", 
         dpi = "retina", width = 6, height = 3, units = "in")

#-----------------------------------------#
# Calculating log ratios based on MMvec findings
#-----------------------------------------#

f1Genera <- c("Malassezia","Trichosporon","Ramularia")
f2Genera <- c("Aspergillus","Candida")
f3Genera <- c("Colletotrichum","Fusarium","Cutaneotrichosporon","Phialocephala","Trichoderma","Talaromyces",
              "Yarrowia","Stereum","Aureobasidium","Hyphopichia","Dissoconium","Agaricus","Exophiala",
              "Alternaria","Tilletiopsis","Cryptococcus","Penicillium","Puccinia")

# Subset with the pseudocount (to avoid dropping many samples)
ucsdRep200FungiGenusSubset <- ucsdRep200FungiGenus[,colnames(ucsdRep200FungiGenus) %in% c(f1Genera,f2Genera,f3Genera)]+1 

lrMetaUCSD <- metaUCSDJoined_Fungi_Nonzero
lrMetaUCSD$lr_labels <- lrMetaUCSD$disease_type_consol
lrMetaUCSD$lr_f1_f2 <- log10(rowSums(ucsdRep200FungiGenusSubset[,colnames(ucsdRep200FungiGenusSubset) %in% f1Genera])/
                           rowSums(ucsdRep200FungiGenusSubset[,colnames(ucsdRep200FungiGenusSubset) %in% f2Genera]))
lrMetaUCSD$lr_f1_f3 <- log10(rowSums(ucsdRep200FungiGenusSubset[,colnames(ucsdRep200FungiGenusSubset) %in% f1Genera])/
                           rowSums(ucsdRep200FungiGenusSubset[,colnames(ucsdRep200FungiGenusSubset) %in% f3Genera]))
lrMetaUCSD$lr_f2_f3 <- log10(rowSums(ucsdRep200FungiGenusSubset[,colnames(ucsdRep200FungiGenusSubset) %in% f2Genera])/
                           rowSums(ucsdRep200FungiGenusSubset[,colnames(ucsdRep200FungiGenusSubset) %in%f3Genera]))

#------------Between cancer types------------#
lrMetaUCSD %>%
  filter(lr_labels != "Control") %>%
  filter(is.finite(lr_f1_f2)) %>%
  filter(!is.na(lr_f1_f2)) %>%
  ggplot(aes(reorder(lr_labels, -lr_f1_f2, FUN=median),lr_f1_f2, fill=lr_labels)) +
  geom_boxplot(notch = TRUE) + scale_fill_nejm() + theme_pubr() + geom_point(alpha=0.4) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  labs(x = "Cancer Types", y = "log(F1/F2)", title = "UCSD: Log ratio of F1/F2 in plasma data") +
  stat_compare_means(method = "anova", label.y = 2.2)
ggsave(filename = "Figures/Other_Figures/log_ratio_ucsd_f1_f2.pdf", dpi = "retina",
       units = "in", width = 4, height = 4)

lrMetaUCSD %>%
  filter(lr_labels != "Control") %>%
  filter(is.finite(lr_f1_f3)) %>%
  filter(!is.na(lr_f1_f3)) %>%
  ggplot(aes(reorder(lr_labels, -lr_f1_f3, FUN=median),lr_f1_f3, fill=lr_labels)) +
  geom_boxplot(notch = TRUE) + scale_fill_nejm() + theme_pubr() + geom_point(alpha=0.4) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  labs(x = "Cancer Types", y = "log(F1/F3)", title = "UCSD: Log ratio of F1/F3 in plasma data") +
  stat_compare_means(method = "anova", label.y = 1.5)
ggsave(filename = "Figures/Other_Figures/log_ratio_ucsd_f1_f3.pdf", dpi = "retina",
       units = "in", width = 4, height = 4)

lrMetaUCSD %>%
  filter(lr_labels != "Control") %>%
  filter(is.finite(lr_f2_f3)) %>%
  filter(!is.na(lr_f2_f3)) %>%
  ggplot(aes(reorder(lr_labels, -lr_f2_f3, FUN=mean),lr_f2_f3, fill=lr_labels)) +
  geom_boxplot(notch = TRUE) + scale_fill_nejm() + theme_pubr() + geom_point(alpha=0.4) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  labs(x = "Cancer Types", y = "log(F2/F3)", title = "UCSD: Log ratio of F2/F3 in plasma data") +
  stat_compare_means(method = "anova", label.y = 0.5)
ggsave(filename = "Figures/Other_Figures/log_ratio_ucsd_f2_f3.pdf", dpi = "retina",
       units = "in", width = 4, height = 4)

#------------Cancer vs Healthy------------#
lrMetaUCSD %>%
  filter(is.finite(lr_f1_f2)) %>%
  filter(!is.na(lr_f1_f2)) %>%
  ggplot(aes(reorder(HvsC, -lr_f1_f2, FUN=median),lr_f1_f2, fill=HvsC)) +
  geom_boxplot(notch = TRUE) + scale_fill_nejm() + theme_pubr() + geom_point(alpha=0.4) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  labs(x = "Disease Status", y = "log(F1/F2)", title = "UCSD: Log ratio of F1/F2 in plasma data") +
  stat_compare_means(method = "wilcox.test", label.y = 2.5, comparisons = list(c("Cancer","Control")))
ggsave(filename = "Figures/Other_Figures/log_ratio_ucsd_f1_f2_cancer_vs_healthy.pdf", dpi = "retina",
       units = "in", width = 4, height = 4)

lrMetaUCSD %>%
  filter(is.finite(lr_f1_f3)) %>%
  filter(!is.na(lr_f1_f3)) %>%
  ggplot(aes(reorder(HvsC, -lr_f1_f3, FUN=median),lr_f1_f3, fill=HvsC)) +
  geom_boxplot(notch = TRUE) + scale_fill_nejm() + theme_pubr() + geom_point(alpha=0.4) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  labs(x = "Disease Status", y = "log(F1/F3)", title = "UCSD: Log ratio of F1/F3 in plasma data") +
  stat_compare_means(method = "wilcox.test", label.y = 1.5, comparisons = list(c("Cancer","Control")))
ggsave(filename = "Figures/Supplementary_Figures/log_ratio_ucsd_f1_f3_cancer_vs_healthy.pdf", dpi = "retina",
       units = "in", width = 4, height = 4)

lrMetaUCSD %>%
  filter(is.finite(lr_f2_f3)) %>%
  filter(!is.na(lr_f2_f3)) %>%
  ggplot(aes(reorder(HvsC, -lr_f2_f3, FUN=mean),lr_f2_f3, fill=HvsC)) +
  geom_boxplot(notch = TRUE) + scale_fill_nejm() + theme_pubr() + geom_point(alpha=0.4) +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  labs(x = "Disease Status", y = "log(F2/F3)", title = "UCSD: Log ratio of F2/F3 in plasma data") +
  stat_compare_means(method = "wilcox.test", label.y = 0.85, comparisons = list(c("Cancer","Control")))
ggsave(filename = "Figures/Other_Figures/log_ratio_ucsd_f2_f3_cancer_vs_healthy.pdf", dpi = "retina",
       units = "in", width = 4, height = 4)

#----------------------------------------------------------------------------------#
# ML repeated CV perf --> ROC plot with 99% CIs
# (1) TopX between Hopkins and UCSD
# (2) WIS-overlapping fungi+bacteria between Hopkins and UCSD
#----------------------------------------------------------------------------------#

#-------------------TopX between Hopkins and UCSD-------------------#
ucsdFungi_TopX <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero,
                              snmData = snmDataUCSD_Fungi_TopX,
                              col2Predict = "HvsC",
                              numResampleIter = 10,
                              dataString = "fungi_topX",
                              varImpFlag = TRUE)
# Load hopkinsFungi_topX object
load("Interim_data/cristianoFungi_topX_ML_model_15Nov21.RData")

plots_hopkins_topX_CV <- plotMLWithCIs(hopkinsFungi_topX, 
                                       sizeAvgCurve=0.5, 
                                       showRepCurves=TRUE, 
                                       ciAlpha = 0.4, 
                                       sizeRepCurves = 0.25,
                                       colorAvgCurve = "#631879FF",
                                       # colorRepCurves = "#631879FF",
                                       ciFillColor = "purple",
                                       ciLevel = 0.99,
                                       positiveClass="Cancer", 
                                       negativeClass="Healthy")
plots_ucsd_topX_CV <- plotMLWithCIs(ucsdFungi_TopX,
                                    sizeAvgCurve=0.5, 
                                    showRepCurves=TRUE, 
                                    ciAlpha = 0.4, 
                                    sizeRepCurves = 0.25,
                                    colorAvgCurve = "#79AF97FF", # "#A20056FF",
                                    ciFillColor = "#79AF97FF",
                                    ciLevel = 0.99,
                                    positiveClass="Cancer", 
                                    negativeClass="Control")
# Combine ROC curves
rocCombined <- plots_hopkins_topX_CV$rocPlot
for(ii in 1:10){
  rocCombined <- rocCombined + 
    geom_path(data = plots_ucsd_topX_CV$rocCurveData[[ii]],aes(x=fpr,y=tpr), color = "lightgray", size = 0.25)
}
rocCombinedAll <- rocCombined +
  geom_ribbon(data = plots_ucsd_topX_CV$interpROCYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
              fill = "#79AF97FF", alpha = 0.3, inherit.aes = F) +
  geom_path(data = plots_ucsd_topX_CV$interpROCYDf_CI, aes(x = xval, y = Estimate), color = "#79AF97FF", size = 0.5) +
  annotate("segment", x = 0.15, xend = 0.25, y = 0.15, yend = 0.15, color = "#79AF97FF") +
  annotate("text", x = 0.27, y = 0.15, color = "#79AF97FF", label = paste0("UCSD: AUROC 99% CI: [", 
                                                                          paste0(100*round(plots_ucsd_topX_CV$aurocCI[2],4),
                                                                                 ", ",100*round(plots_ucsd_topX_CV$aurocCI[3],4)),"]"), hjust = 0) +
  annotate("segment", x = 0.15, xend = 0.25, y = 0.10, yend = 0.10, color = "#631879FF") +
  annotate("text", x = 0.27, y = 0.10, color = "#631879FF", label = paste0("Hopkins: AUROC 99% CI: [", 
                                                                           paste0(100*round(plots_hopkins_topX_CV$aurocCI[2],4),
                                                                                  ", ",100*round(plots_hopkins_topX_CV$aurocCI[3],4)),"]"), hjust = 0)

# Combine PR curves
prCombined <- plots_hopkins_topX_CV$prPlot
for(ii in 1:10){
  prCombined <- prCombined + 
    geom_path(data = plots_ucsd_topX_CV$prCurveData[[ii]],aes(x=recall,y=precision), color = "lightgray", size = 0.25)
}
prCombinedAll <- prCombined +
  geom_ribbon(data = plots_ucsd_topX_CV$interpPRYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
              fill = "#79AF97FF", alpha = 0.3, inherit.aes = F) +
  geom_path(data = plots_ucsd_topX_CV$interpPRYDf_CI, aes(x = xval, y = Estimate), color = "#79AF97FF", size = 0.5) +
  annotate("segment", x = 0.05, xend = 0.15, y = 0.15, yend = 0.15, color = "#79AF97FF") +
  annotate("text", x = 0.17, y = 0.15, color = "#79AF97FF", label = paste0("UCSD: AUPR 99% CI: [", 
                                                                          paste0(100*round(plots_ucsd_topX_CV$auprCI[2],4),
                                                                                 ", ",100*round(plots_ucsd_topX_CV$auprCI[3],4)),"]"), hjust = 0) +
  annotate("segment", x = 0.05, xend = 0.15, y = 0.1, yend = 0.1, color = "#631879FF") +
  annotate("text", x = 0.17, y = 0.10, color = "#631879FF", label = paste0("Hopkins: AUPR 99% CI: [", 
                                                                           paste0(100*round(plots_hopkins_topX_CV$auprCI[2],4),
                                                                                  ", ",100*round(plots_hopkins_topX_CV$auprCI[3],4)),"]"), hjust = 0)
combinedPlotTitle <- paste0("Hopkins and UCSD plasma cohorts using top 20 fungi from Hopkins cohort\nCancer vs Healthy (10-fold CV repeated 10 times)")
combinedOverlayPlots <- ggarrange(rocCombinedAll, prCombinedAll, ncol = 2)
combinedOverlayPlotsAnnotated <- annotate_figure(combinedOverlayPlots, 
                                                 top = text_grob(combinedPlotTitle, 
                                                                 color = "black", face = "bold", size = 14))
print(combinedOverlayPlotsAnnotated)
ggsave(filename = "Figures/Main_Figures/roc_pr_CV_overlay_hopkins_ucsd_topX_fungi.svg",
       units = "in", width = 10, height = 5)

#----------------------------------------------------------------------------------#
# Machine learning: Individual cancer vs. Healthy
#----------------------------------------------------------------------------------#

metaUCSDJoined_Fungi_Nonzero %>% count(disease_type_consol)
# disease_type_consol  n
# 1             Control 68
# 2               NSCLC 25
# 3                PRAD 57
# 4                SKCM 16

#---------------------------------Full data---------------------------------#
ml_HvsC_Fullrep200_OGU_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_FullRep200_Matched,
                                       snmData = snmDataUCSD_FullRep200_OGU,
                                       col2Predict = "one_cancer_vs_controls",
                                       dzOfInterest = "NSCLC",
                                       dataString = "full_rep200",
                                       varImpFlag = FALSE)
ml_HvsC_Fullrep200_OGU_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_FullRep200_Matched,
                                       snmData = snmDataUCSD_FullRep200_OGU,
                                       col2Predict = "one_cancer_vs_controls",
                                       dzOfInterest = "PRAD",
                                       dataString = "full_rep200",
                                       varImpFlag = FALSE)
ml_HvsC_Fullrep200_OGU_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_FullRep200_Matched,
                                            snmData = snmDataUCSD_FullRep200_OGU,
                                            col2Predict = "one_cancer_vs_controls",
                                            dzOfInterest = "SKCM",
                                            dataString = "full_rep200",
                                            varImpFlag = FALSE)

#---------------------------------Bacteria---------------------------------#
ml_HvsC_Bacteria_OGU_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                     snmData = snmDataUCSD_Bacteria_OGU,
                                     col2Predict = "one_cancer_vs_controls",
                                     dzOfInterest = "NSCLC",
                                     dataString = "bacteria_only",
                                     varImpFlag = FALSE)
ml_HvsC_Bacteria_OGU_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                     snmData = snmDataUCSD_Bacteria_OGU,
                                     col2Predict = "one_cancer_vs_controls",
                                     dzOfInterest = "PRAD",
                                     dataString = "bacteria_only",
                                     varImpFlag = FALSE)
ml_HvsC_Bacteria_OGU_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                     snmData = snmDataUCSD_Bacteria_OGU,
                                     col2Predict = "one_cancer_vs_controls",
                                     dzOfInterest = "SKCM",
                                     dataString = "bacteria_only",
                                     varImpFlag = FALSE)

#---------------------------------Fungi---------------------------------#
ml_HvsC_Fungi_Species_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                                 snmData = snmData_ucsdRep200FungiSpecies,
                                                 col2Predict = "one_cancer_vs_controls",
                                                 dzOfInterest = "NSCLC",
                                                 dataString = "fungi_only",
                                                 varImpFlag = FALSE)
ml_HvsC_Fungi_Species_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                                snmData = snmData_ucsdRep200FungiSpecies,
                                                col2Predict = "one_cancer_vs_controls",
                                                dzOfInterest = "PRAD",
                                                dataString = "fungi_only",
                                                varImpFlag = FALSE)
ml_HvsC_Fungi_Species_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                                snmData = snmData_ucsdRep200FungiSpecies,
                                                col2Predict = "one_cancer_vs_controls",
                                                dzOfInterest = "SKCM",
                                                dataString = "fungi_only",
                                                varImpFlag = FALSE)

#---------------------------------Decontam Fungi---------------------------------#
ml_HvsC_Decontam_Fungi_OGU_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                           snmData = snmDataUCSD_Fungi_Decontam_OGU,
                                           col2Predict = "one_cancer_vs_controls",
                                           dzOfInterest = "NSCLC",
                                           dataString = "fungi_decontam",
                                           varImpFlag = FALSE)
ml_HvsC_Decontam_Fungi_OGU_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                           snmData = snmDataUCSD_Fungi_Decontam_OGU,
                                           col2Predict = "one_cancer_vs_controls",
                                           dzOfInterest = "PRAD",
                                           dataString = "fungi_decontam",
                                           varImpFlag = FALSE)
ml_HvsC_Decontam_Fungi_OGU_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                           snmData = snmDataUCSD_Fungi_Decontam_OGU,
                                           col2Predict = "one_cancer_vs_controls",
                                           dzOfInterest = "SKCM",
                                           dataString = "fungi_decontam",
                                           varImpFlag = FALSE)

#---------------------------------Intersected fungi---------------------------------#
ml_HvsC_Shared_Fungi_Species_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                             snmData = snmData_ucsdRep200FungiSpeciesShared,
                                             col2Predict = "one_cancer_vs_controls",
                                             dzOfInterest = "NSCLC",
                                             dataString = "fungi_intersected_with_Weizmann",
                                             varImpFlag = FALSE)
ml_HvsC_Shared_Fungi_Species_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                             snmData = snmData_ucsdRep200FungiSpeciesShared,
                                             col2Predict = "one_cancer_vs_controls",
                                             dzOfInterest = "PRAD",
                                             dataString = "fungi_intersected_with_Weizmann",
                                             varImpFlag = FALSE)
ml_HvsC_Shared_Fungi_Species_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                             snmData = snmData_ucsdRep200FungiSpeciesShared,
                                             col2Predict = "one_cancer_vs_controls",
                                             dzOfInterest = "SKCM",
                                             dataString = "fungi_intersected_with_Weizmann",
                                             varImpFlag = FALSE)

#---------------------------------Intersected bacteria---------------------------------#
ml_HvsC_Shared_Bacteria_Species_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                                   snmData = snmData_ucsdRep200BacteriaSpeciesShared,
                                                   col2Predict = "one_cancer_vs_controls",
                                                   dzOfInterest = "NSCLC",
                                                   dataString = "bacteria_intersected_with_Weizmann",
                                                   varImpFlag = FALSE)
ml_HvsC_Shared_Bacteria_Species_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                                  snmData = snmData_ucsdRep200BacteriaSpeciesShared,
                                                  col2Predict = "one_cancer_vs_controls",
                                                  dzOfInterest = "PRAD",
                                                  dataString = "bacteria_intersected_with_Weizmann",
                                                  varImpFlag = FALSE)
ml_HvsC_Shared_Bacteria_Species_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                                  snmData = snmData_ucsdRep200BacteriaSpeciesShared,
                                                  col2Predict = "one_cancer_vs_controls",
                                                  dzOfInterest = "SKCM",
                                                  dataString = "bacteria_intersected_with_Weizmann",
                                                  varImpFlag = FALSE)

#---------------------------------Intersected fungi+bacteria---------------------------------#
ml_HvsC_Shared_FungiBacteria_Species_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                                      snmData = cbind(snmData_ucsdRep200FungiSpeciesShared,
                                                                      snmData_ucsdRep200BacteriaSpeciesShared[rownames(snmData_ucsdRep200FungiSpeciesShared),]),
                                                      col2Predict = "one_cancer_vs_controls",
                                                      dzOfInterest = "NSCLC",
                                                      dataString = "fungi_and_bacteria_intersected_with_Weizmann",
                                                      varImpFlag = FALSE)
ml_HvsC_Shared_FungiBacteria_Species_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                                     snmData = cbind(snmData_ucsdRep200FungiSpeciesShared,
                                                                     snmData_ucsdRep200BacteriaSpeciesShared[rownames(snmData_ucsdRep200FungiSpeciesShared),]),
                                                     col2Predict = "one_cancer_vs_controls",
                                                     dzOfInterest = "PRAD",
                                                     dataString = "fungi_and_bacteria_intersected_with_Weizmann",
                                                     varImpFlag = FALSE)
ml_HvsC_Shared_FungiBacteria_Species_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                                     snmData = cbind(snmData_ucsdRep200FungiSpeciesShared,
                                                                     snmData_ucsdRep200BacteriaSpeciesShared[rownames(snmData_ucsdRep200FungiSpeciesShared),]),
                                                     col2Predict = "one_cancer_vs_controls",
                                                     dzOfInterest = "SKCM",
                                                     dataString = "fungi_and_bacteria_intersected_with_Weizmann",
                                                     varImpFlag = FALSE)

#---------------------------------TopX fungi---------------------------------#
ml_HvsC_TopX_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero,
                                   snmData = snmDataUCSD_Fungi_TopX,
                                   col2Predict = "one_cancer_vs_controls",
                                   dzOfInterest = "NSCLC",
                                   dataString = "topX_fungi",
                                   varImpFlag = FALSE)
ml_HvsC_TopX_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero,
                                  snmData = snmDataUCSD_Fungi_TopX,
                                  col2Predict = "one_cancer_vs_controls",
                                  dzOfInterest = "PRAD",
                                  dataString = "topX_fungi",
                                  varImpFlag = FALSE)
ml_HvsC_TopX_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero,
                                  snmData = snmDataUCSD_Fungi_TopX,
                                  col2Predict = "one_cancer_vs_controls",
                                  dzOfInterest = "SKCM",
                                  dataString = "topX_fungi",
                                  varImpFlag = FALSE)

#---------------------------------Combine results---------------------------------#
ucsdPerCancerVsHealthyResults <- rbind(ml_HvsC_Fullrep200_OGU_NSCLC$rep_perf,
                                    ml_HvsC_Fullrep200_OGU_PRAD$rep_perf,
                                    ml_HvsC_Fullrep200_OGU_SKCM$rep_perf,
                                    
                                    ml_HvsC_Bacteria_OGU_NSCLC$rep_perf,
                                    ml_HvsC_Bacteria_OGU_PRAD$rep_perf,
                                    ml_HvsC_Bacteria_OGU_SKCM$rep_perf,
                                    
                                    ml_HvsC_Fungi_Species_NSCLC$rep_perf,
                                    ml_HvsC_Fungi_Species_PRAD$rep_perf,
                                    ml_HvsC_Fungi_Species_SKCM$rep_perf,
                                    
                                    ml_HvsC_Decontam_Fungi_OGU_NSCLC$rep_perf,
                                    ml_HvsC_Decontam_Fungi_OGU_PRAD$rep_perf,
                                    ml_HvsC_Decontam_Fungi_OGU_SKCM$rep_perf,
                                    
                                    ml_HvsC_Shared_Fungi_Species_NSCLC$rep_perf,
                                    ml_HvsC_Shared_Fungi_Species_PRAD$rep_perf,
                                    ml_HvsC_Shared_Fungi_Species_SKCM$rep_perf,
                                    
                                    ml_HvsC_Shared_Bacteria_Species_NSCLC$rep_perf,
                                    ml_HvsC_Shared_Bacteria_Species_PRAD$rep_perf,
                                    ml_HvsC_Shared_Bacteria_Species_SKCM$rep_perf,
                                    
                                    ml_HvsC_Shared_FungiBacteria_Species_NSCLC$rep_perf,
                                    ml_HvsC_Shared_FungiBacteria_Species_PRAD$rep_perf,
                                    ml_HvsC_Shared_FungiBacteria_Species_SKCM$rep_perf,
                                    
                                    ml_HvsC_TopX_NSCLC$rep_perf,
                                    ml_HvsC_TopX_PRAD$rep_perf,
                                    ml_HvsC_TopX_SKCM$rep_perf)

colnames(ucsdPerCancerVsHealthyResults)[1:2] <- c("AUROC","AUPR")
ucsdPerCancerVsHealthyResults$nullAUROC <- 0.5
ucsdPerCancerVsHealthyResults$nullAUPR <- ucsdPerCancerVsHealthyResults$positiveClassSize/
  (ucsdPerCancerVsHealthyResults$positiveClassSize+ucsdPerCancerVsHealthyResults$negativeClassSize)

ucsdPerCancerVsHealthyResults$dataString[ucsdPerCancerVsHealthyResults$dataString=="full_rep200"] <- "Full multikingdom database (rep200)"
ucsdPerCancerVsHealthyResults$dataString[ucsdPerCancerVsHealthyResults$dataString=="bacteria_only"] <- "Bacteria all (Species)"
ucsdPerCancerVsHealthyResults$dataString[ucsdPerCancerVsHealthyResults$dataString=="fungi_only"] <- "Fungi all (Species)"
ucsdPerCancerVsHealthyResults$dataString[ucsdPerCancerVsHealthyResults$dataString=="fungi_decontam"] <- "Fungi decontaminated (Species)"
ucsdPerCancerVsHealthyResults$dataString[ucsdPerCancerVsHealthyResults$dataString=="fungi_intersected_with_Weizmann"] <- "Fungi ∩ WIS (Species)"
ucsdPerCancerVsHealthyResults$dataString[ucsdPerCancerVsHealthyResults$dataString=="bacteria_intersected_with_Weizmann"] <- "Bacteria ∩ WIS (Species)"
ucsdPerCancerVsHealthyResults$dataString[ucsdPerCancerVsHealthyResults$dataString=="fungi_and_bacteria_intersected_with_Weizmann"] <- "Fungi+bacteria ∩ WIS (Species)"
ucsdPerCancerVsHealthyResults$dataString[ucsdPerCancerVsHealthyResults$dataString=="topX_fungi"] <- "Top 20 fungi in Hopkins cohort (Species)"

ucsdPerCancerVsHealthyResults$dataString <- factor(ucsdPerCancerVsHealthyResults$dataString, levels = c("Full multikingdom database (rep200)",
                                                                                                  "Fungi+bacteria ∩ WIS (Species)",
                                                                                                  "Bacteria all (Species)",
                                                                                                  "Fungi all (Species)",
                                                                                                  "Fungi decontaminated (Species)",
                                                                                                  "Fungi ∩ WIS (Species)",
                                                                                                  "Bacteria ∩ WIS (Species)",
                                                                                                  "Top 20 fungi in Hopkins cohort (Species)"))

save(ucsdPerCancerVsHealthyResults,
     file = "Interim_data/ucsdPerCancerVsHealthyResults_4Apr22.RData")

source("Supporting_scripts/S00-SummarySE.R")
# Facet by AUROC vs AUPR; color by full dataset vs fungi (dataString)
ucsdPerCancerVsHealthyResults %>%
  filter(grepl("Full|\\+|decontam|Top",dataString)) %>%
  # filter(grepl("\\+|Bacteria ∩",dataString)) %>%
  reshape2::melt(id.vars = c("dataString","rep","diseaseType", "col2Predict","nullAUPR","nullAUROC","positiveClassSize","negativeClassSize")) %>%
  summarySE(measurevar = "value", groupvars = c("diseaseType","dataString","variable","nullAUPR","nullAUROC","positiveClassSize","negativeClassSize")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=dataString)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Individual cancer type versus all healthy samples") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + ggtitle("UCSD plasma validation cohort:\nPer cancer type vs. healthy") +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  rotate_x_text(0) + scale_color_nejm(name = "Dataset") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Supplementary_Figures/ucsd_per_cancer_type_vs_healthy_full_WIS_fungi_topX_4Apr22.svg", 
       dpi = "retina", width = 7, height = 4, units = "in")

#----------------------------------------------------------------------------------#
# Machine learning: Between cancer types
#----------------------------------------------------------------------------------#

metaUCSDJoined_Fungi_Nonzero %>% count(disease_type_consol)
# disease_type_consol  n
# 1             Control 68
# 2               NSCLC 25
# 3                PRAD 57
# 4                SKCM 16

source("00-Functions.R")
#---------------------------------Full data---------------------------------#
ml_1vsAll_Fullrep200_OGU_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_FullRep200_Matched,
                                             snmData = snmDataUCSD_FullRep200_OGU,
                                             col2Predict = "one_cancer_vs_others",
                                             dzOfInterest = "NSCLC",
                                             dataString = "full_rep200",
                                             varImpFlag = FALSE)
ml_1vsAll_Fullrep200_OGU_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_FullRep200_Matched,
                                            snmData = snmDataUCSD_FullRep200_OGU,
                                            col2Predict = "one_cancer_vs_others",
                                            dzOfInterest = "PRAD",
                                            dataString = "full_rep200",
                                            varImpFlag = FALSE)
ml_1vsAll_Fullrep200_OGU_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_FullRep200_Matched,
                                            snmData = snmDataUCSD_FullRep200_OGU,
                                            col2Predict = "one_cancer_vs_others",
                                            dzOfInterest = "SKCM",
                                            dataString = "full_rep200",
                                            varImpFlag = FALSE)

#---------------------------------Bacteria---------------------------------#
ml_1vsAll_Bacteria_OGU_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                           snmData = snmDataUCSD_Bacteria_OGU,
                                           col2Predict = "one_cancer_vs_others",
                                           dzOfInterest = "NSCLC",
                                           dataString = "bacteria_only",
                                           varImpFlag = FALSE)
ml_1vsAll_Bacteria_OGU_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                          snmData = snmDataUCSD_Bacteria_OGU,
                                          col2Predict = "one_cancer_vs_others",
                                          dzOfInterest = "PRAD",
                                          dataString = "bacteria_only",
                                          varImpFlag = FALSE)
ml_1vsAll_Bacteria_OGU_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                          snmData = snmDataUCSD_Bacteria_OGU,
                                          col2Predict = "one_cancer_vs_others",
                                          dzOfInterest = "SKCM",
                                          dataString = "bacteria_only",
                                          varImpFlag = FALSE)

#---------------------------------Fungi---------------------------------#
ml_1vsAll_Fungi_Species_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                            snmData = snmData_ucsdRep200FungiSpecies,
                                            col2Predict = "one_cancer_vs_others",
                                            dzOfInterest = "NSCLC",
                                            dataString = "fungi_only",
                                            varImpFlag = FALSE)
ml_1vsAll_Fungi_Species_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                           snmData = snmData_ucsdRep200FungiSpecies,
                                           col2Predict = "one_cancer_vs_others",
                                           dzOfInterest = "PRAD",
                                           dataString = "fungi_only",
                                           varImpFlag = FALSE)
ml_1vsAll_Fungi_Species_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                           snmData = snmData_ucsdRep200FungiSpecies,
                                           col2Predict = "one_cancer_vs_others",
                                           dzOfInterest = "SKCM",
                                           dataString = "fungi_only",
                                           varImpFlag = FALSE)

#---------------------------------Decontam Fungi---------------------------------#
ml_1vsAll_Decontam_Fungi_OGU_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                                 snmData = snmDataUCSD_Fungi_Decontam_OGU,
                                                 col2Predict = "one_cancer_vs_others",
                                                 dzOfInterest = "NSCLC",
                                                 dataString = "fungi_decontam",
                                                 varImpFlag = FALSE)
ml_1vsAll_Decontam_Fungi_OGU_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                                snmData = snmDataUCSD_Fungi_Decontam_OGU,
                                                col2Predict = "one_cancer_vs_others",
                                                dzOfInterest = "PRAD",
                                                dataString = "fungi_decontam",
                                                varImpFlag = FALSE)
ml_1vsAll_Decontam_Fungi_OGU_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                                snmData = snmDataUCSD_Fungi_Decontam_OGU,
                                                col2Predict = "one_cancer_vs_others",
                                                dzOfInterest = "SKCM",
                                                dataString = "fungi_decontam",
                                                varImpFlag = FALSE)

#---------------------------------Intersected fungi---------------------------------#
ml_1vsAll_Shared_Fungi_Species_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                                   snmData = snmData_ucsdRep200FungiSpeciesShared,
                                                   col2Predict = "one_cancer_vs_others",
                                                   dzOfInterest = "NSCLC",
                                                   dataString = "fungi_intersected_with_Weizmann",
                                                   varImpFlag = FALSE)
ml_1vsAll_Shared_Fungi_Species_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                                  snmData = snmData_ucsdRep200FungiSpeciesShared,
                                                  col2Predict = "one_cancer_vs_others",
                                                  dzOfInterest = "PRAD",
                                                  dataString = "fungi_intersected_with_Weizmann",
                                                  varImpFlag = FALSE)
ml_1vsAll_Shared_Fungi_Species_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                                  snmData = snmData_ucsdRep200FungiSpeciesShared,
                                                  col2Predict = "one_cancer_vs_others",
                                                  dzOfInterest = "SKCM",
                                                  dataString = "fungi_intersected_with_Weizmann",
                                                  varImpFlag = FALSE)

#---------------------------------Intersected bacteria---------------------------------#
ml_1vsAll_Shared_Bacteria_Species_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                                      snmData = snmData_ucsdRep200BacteriaSpeciesShared,
                                                      col2Predict = "one_cancer_vs_others",
                                                      dzOfInterest = "NSCLC",
                                                      dataString = "bacteria_intersected_with_Weizmann",
                                                      varImpFlag = FALSE)
ml_1vsAll_Shared_Bacteria_Species_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                                     snmData = snmData_ucsdRep200BacteriaSpeciesShared,
                                                     col2Predict = "one_cancer_vs_others",
                                                     dzOfInterest = "PRAD",
                                                     dataString = "bacteria_intersected_with_Weizmann",
                                                     varImpFlag = FALSE)
ml_1vsAll_Shared_Bacteria_Species_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                                     snmData = snmData_ucsdRep200BacteriaSpeciesShared,
                                                     col2Predict = "one_cancer_vs_others",
                                                     dzOfInterest = "SKCM",
                                                     dataString = "bacteria_intersected_with_Weizmann",
                                                     varImpFlag = FALSE)

#---------------------------------Intersected fungi+bacteria---------------------------------#
ml_1vsAll_Shared_FungiBacteria_Species_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                                           snmData = cbind(snmData_ucsdRep200FungiSpeciesShared,
                                                                           snmData_ucsdRep200BacteriaSpeciesShared[rownames(snmData_ucsdRep200FungiSpeciesShared),]),
                                                           col2Predict = "one_cancer_vs_others",
                                                           dzOfInterest = "NSCLC",
                                                           dataString = "fungi_and_bacteria_intersected_with_Weizmann",
                                                           varImpFlag = FALSE)
ml_1vsAll_Shared_FungiBacteria_Species_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                                          snmData = cbind(snmData_ucsdRep200FungiSpeciesShared,
                                                                          snmData_ucsdRep200BacteriaSpeciesShared[rownames(snmData_ucsdRep200FungiSpeciesShared),]),
                                                          col2Predict = "one_cancer_vs_others",
                                                          dzOfInterest = "PRAD",
                                                          dataString = "fungi_and_bacteria_intersected_with_Weizmann",
                                                          varImpFlag = FALSE)
ml_1vsAll_Shared_FungiBacteria_Species_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                                          snmData = cbind(snmData_ucsdRep200FungiSpeciesShared,
                                                                          snmData_ucsdRep200BacteriaSpeciesShared[rownames(snmData_ucsdRep200FungiSpeciesShared),]),
                                                          col2Predict = "one_cancer_vs_others",
                                                          dzOfInterest = "SKCM",
                                                          dataString = "fungi_and_bacteria_intersected_with_Weizmann",
                                                          varImpFlag = FALSE)

#---------------------------------TopX fungi---------------------------------#
ml_1vsAll_TopX_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero,
                                   snmData = snmDataUCSD_Fungi_TopX,
                                   col2Predict = "one_cancer_vs_others",
                                   dzOfInterest = "NSCLC",
                                   dataString = "topX_fungi",
                                   varImpFlag = FALSE)
ml_1vsAll_TopX_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero,
                                  snmData = snmDataUCSD_Fungi_TopX,
                                  col2Predict = "one_cancer_vs_others",
                                  dzOfInterest = "PRAD",
                                  dataString = "topX_fungi",
                                  varImpFlag = FALSE)
ml_1vsAll_TopX_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero,
                                  snmData = snmDataUCSD_Fungi_TopX,
                                  col2Predict = "one_cancer_vs_others",
                                  dzOfInterest = "SKCM",
                                  dataString = "topX_fungi",
                                  varImpFlag = FALSE)

#---------------------------------Combine results---------------------------------#
ucsd1VsAllResults <- rbind(ml_1vsAll_Fullrep200_OGU_NSCLC$rep_perf,
                           ml_1vsAll_Fullrep200_OGU_PRAD$rep_perf,
                           ml_1vsAll_Fullrep200_OGU_SKCM$rep_perf,
                           
                           ml_1vsAll_Bacteria_OGU_NSCLC$rep_perf,
                           ml_1vsAll_Bacteria_OGU_PRAD$rep_perf,
                           ml_1vsAll_Bacteria_OGU_SKCM$rep_perf,
                           
                           ml_1vsAll_Fungi_Species_NSCLC$rep_perf,
                           ml_1vsAll_Fungi_Species_PRAD$rep_perf,
                           ml_1vsAll_Fungi_Species_SKCM$rep_perf,
                           
                           ml_1vsAll_Decontam_Fungi_OGU_NSCLC$rep_perf,
                           ml_1vsAll_Decontam_Fungi_OGU_PRAD$rep_perf,
                           ml_1vsAll_Decontam_Fungi_OGU_SKCM$rep_perf,
                           
                           ml_1vsAll_Shared_Fungi_Species_NSCLC$rep_perf,
                           ml_1vsAll_Shared_Fungi_Species_PRAD$rep_perf,
                           ml_1vsAll_Shared_Fungi_Species_SKCM$rep_perf,
                           
                           ml_1vsAll_Shared_Bacteria_Species_NSCLC$rep_perf,
                           ml_1vsAll_Shared_Bacteria_Species_PRAD$rep_perf,
                           ml_1vsAll_Shared_Bacteria_Species_SKCM$rep_perf,
                           
                           ml_1vsAll_Shared_FungiBacteria_Species_NSCLC$rep_perf,
                           ml_1vsAll_Shared_FungiBacteria_Species_PRAD$rep_perf,
                           ml_1vsAll_Shared_FungiBacteria_Species_SKCM$rep_perf,
                           
                           ml_1vsAll_TopX_NSCLC$rep_perf,
                           ml_1vsAll_TopX_PRAD$rep_perf,
                           ml_1vsAll_TopX_SKCM$rep_perf)

colnames(ucsd1VsAllResults)[1:2] <- c("AUROC","AUPR")
ucsd1VsAllResults$nullAUROC <- 0.5
ucsd1VsAllResults$nullAUPR <- ucsd1VsAllResults$positiveClassSize/
  (ucsd1VsAllResults$positiveClassSize+ucsd1VsAllResults$negativeClassSize)

ucsd1VsAllResults$dataString[ucsd1VsAllResults$dataString=="full_rep200"] <- "Full multikingdom database (rep200)"
ucsd1VsAllResults$dataString[ucsd1VsAllResults$dataString=="bacteria_only"] <- "Bacteria all (Species)"
ucsd1VsAllResults$dataString[ucsd1VsAllResults$dataString=="fungi_only"] <- "Fungi all (Species)"
ucsd1VsAllResults$dataString[ucsd1VsAllResults$dataString=="fungi_decontam"] <- "Fungi decontaminated (Species)"
ucsd1VsAllResults$dataString[ucsd1VsAllResults$dataString=="fungi_intersected_with_Weizmann"] <- "Fungi ∩ WIS (Species)"
ucsd1VsAllResults$dataString[ucsd1VsAllResults$dataString=="bacteria_intersected_with_Weizmann"] <- "Bacteria ∩ WIS (Species)"
ucsd1VsAllResults$dataString[ucsd1VsAllResults$dataString=="fungi_and_bacteria_intersected_with_Weizmann"] <- "Fungi+bacteria ∩ WIS (Species)"
ucsd1VsAllResults$dataString[ucsd1VsAllResults$dataString=="topX_fungi"] <- "Top 20 fungi in Hopkins cohort (Species)"

ucsd1VsAllResults$dataString <- factor(ucsd1VsAllResults$dataString, levels = c("Full multikingdom database (rep200)",
                                                                                                        "Fungi+bacteria ∩ WIS (Species)",
                                                                                                        "Bacteria all (Species)",
                                                                                                        "Fungi all (Species)",
                                                                                                        "Fungi decontaminated (Species)",
                                                                                                        "Fungi ∩ WIS (Species)",
                                                                                                        "Bacteria ∩ WIS (Species)",
                                                                                                        "Top 20 fungi in Hopkins cohort (Species)"))
save(ucsd1VsAllResults,
     file = "Interim_data/ucsd1VsAllResults_4Apr22.RData")

source("Supporting_scripts/S00-SummarySE.R")
# NOTE: TopX fungi were selected based on cancer vs healthy comparisons
# so they are not plotted here
ucsd1VsAllResults %>%
  filter(grepl("Full|\\+|decontam",dataString)) %>%
  reshape2::melt(id.vars = c("dataString","rep","diseaseType", "col2Predict","nullAUPR","nullAUROC","positiveClassSize","negativeClassSize")) %>%
  summarySE(measurevar = "value", groupvars = c("diseaseType","dataString","variable","nullAUPR","nullAUROC","positiveClassSize","negativeClassSize")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=dataString)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer types") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + ggtitle("UCSD plasma validation cohort:\nOne cancer type vs all others") +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  rotate_x_text(0) + scale_color_nejm(name = "Dataset") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Supplementary_Figures/ucsd_between_cancer_types_full_WIS_fungi_4Apr22.svg", 
       dpi = "retina", width = 7, height = 4, units = "in")

# ucsd1VsAllResults %>%
  # # filter(grepl("Full|\\+|decontam",dataString)) %>%
  # filter(grepl("\\+|Bacteria ∩",dataString)) %>%
  # filter(diseaseType == "SKCM") %>%
  # ggboxplot(x = "dataString",
  #           y = "AUPR") +
  # stat_compare_means(method = "wilcox")

ucsd1VsAllResults %>%
  filter(grepl("Full",dataString)) %>% 
  pull(AUROC) %>% ci()
# Estimate   CI lower   CI upper Std. Error 
# 0.85538624 0.78539109 0.92538140 0.03422362 

ucsd1VsAllResults %>%
  filter(grepl("\\+",dataString)) %>% 
  pull(AUROC) %>% ci()
# Estimate   CI lower   CI upper Std. Error 
# 0.8146190  0.7414845  0.8877536  0.0357586 

ucsd1VsAllResults %>%
  filter(grepl("decontam",dataString)) %>% 
  pull(AUROC) %>% ci()
# Estimate   CI lower   CI upper Std. Error 
# 0.7762884  0.7028096  0.8497671  0.0359269 

#----------------------------------------------------------------------------------#
# Machine learning: Aggregated cancer vs. healthy for varying taxa levels for bacteria and fungi
#----------------------------------------------------------------------------------#

#---------------------------------Bacteria---------------------------------#
ml_HvsC_Bacteria_Phylum <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                     snmData = snmData_ucsdRep200BacteriaPhylum,
                                     col2Predict = "HvsC",
                                     dataString = "Bacteria",
                                     varImpFlag = FALSE)
ml_HvsC_Bacteria_Class <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                     snmData = snmData_ucsdRep200BacteriaClass,
                                     col2Predict = "HvsC",
                                     dataString = "Bacteria",
                                     varImpFlag = FALSE)
ml_HvsC_Bacteria_Order <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                     snmData = snmData_ucsdRep200BacteriaOrder,
                                     col2Predict = "HvsC",
                                     dataString = "Bacteria",
                                     varImpFlag = FALSE)
ml_HvsC_Bacteria_Family <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                     snmData = snmData_ucsdRep200BacteriaFamily,
                                     col2Predict = "HvsC",
                                     dataString = "Bacteria",
                                     varImpFlag = FALSE)
ml_HvsC_Bacteria_Genus <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                     snmData = snmData_ucsdRep200BacteriaGenus,
                                     col2Predict = "HvsC",
                                     dataString = "Bacteria",
                                     varImpFlag = FALSE)
ml_HvsC_Bacteria_Species <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                     snmData = snmData_ucsdRep200BacteriaSpecies,
                                     col2Predict = "HvsC",
                                     dataString = "Bacteria",
                                     varImpFlag = FALSE)

#---------------------------------Fungi---------------------------------#
ml_HvsC_Fungi_Phylum <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                        snmData = snmData_ucsdRep200FungiPhylum,
                                        col2Predict = "HvsC",
                                        dataString = "Fungi_full",
                                        varImpFlag = FALSE)
ml_HvsC_Fungi_Class <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                       snmData = snmData_ucsdRep200FungiClass,
                                       col2Predict = "HvsC",
                                       dataString = "Fungi_full",
                                       varImpFlag = FALSE)
ml_HvsC_Fungi_Order <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                       snmData = snmData_ucsdRep200FungiOrder,
                                       col2Predict = "HvsC",
                                       dataString = "Fungi_full",
                                       varImpFlag = FALSE)
ml_HvsC_Fungi_Family <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                        snmData = snmData_ucsdRep200FungiFamily,
                                        col2Predict = "HvsC",
                                        dataString = "Fungi_full",
                                        varImpFlag = FALSE)
ml_HvsC_Fungi_Genus <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                       snmData = snmData_ucsdRep200FungiGenus,
                                       col2Predict = "HvsC",
                                       dataString = "Fungi_full",
                                       varImpFlag = FALSE)
ml_HvsC_Fungi_Species <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                         snmData = snmData_ucsdRep200FungiSpecies,
                                         col2Predict = "HvsC",
                                         dataString = "Fungi_full",
                                         varImpFlag = FALSE)

#---------------------------------Decontaminated Fungi---------------------------------#
ml_HvsC_FungiDecontam_Phylum <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                     snmData = snmData_ucsdRep200FungiDecontamPhylum,
                                     col2Predict = "HvsC",
                                     dataString = "Fungi_decontam",
                                     varImpFlag = FALSE)
ml_HvsC_FungiDecontam_Class <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                    snmData = snmData_ucsdRep200FungiDecontamClass,
                                    col2Predict = "HvsC",
                                    dataString = "Fungi_decontam",
                                    varImpFlag = FALSE)
ml_HvsC_FungiDecontam_Order <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                    snmData = snmData_ucsdRep200FungiDecontamOrder,
                                    col2Predict = "HvsC",
                                    dataString = "Fungi_decontam",
                                    varImpFlag = FALSE)
ml_HvsC_FungiDecontam_Family <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                     snmData = snmData_ucsdRep200FungiDecontamFamily,
                                     col2Predict = "HvsC",
                                     dataString = "Fungi_decontam",
                                     varImpFlag = FALSE)
ml_HvsC_FungiDecontam_Genus <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                    snmData = snmData_ucsdRep200FungiDecontamGenus,
                                    col2Predict = "HvsC",
                                    dataString = "Fungi_decontam",
                                    varImpFlag = FALSE)
ml_HvsC_FungiDecontam_Species <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                      snmData = snmData_ucsdRep200FungiDecontamSpecies,
                                      col2Predict = "HvsC",
                                      dataString = "Fungi_decontam",
                                      varImpFlag = FALSE)

#---------------------------------Combine results---------------------------------#
ucsdTaxaLevel_results <- cbind(rbind(ml_HvsC_Bacteria_Phylum$rep_perf,
                                     ml_HvsC_Bacteria_Class$rep_perf,
                                     ml_HvsC_Bacteria_Order$rep_perf,
                                     ml_HvsC_Bacteria_Family$rep_perf,
                                     ml_HvsC_Bacteria_Genus$rep_perf,
                                     ml_HvsC_Bacteria_Species$rep_perf,
                                     
                                     ml_HvsC_Fungi_Phylum$rep_perf,
                                     ml_HvsC_Fungi_Class$rep_perf,
                                     ml_HvsC_Fungi_Order$rep_perf,
                                     ml_HvsC_Fungi_Family$rep_perf,
                                     ml_HvsC_Fungi_Genus$rep_perf,
                                     ml_HvsC_Fungi_Species$rep_perf,
                                     
                                     ml_HvsC_FungiDecontam_Phylum$rep_perf,
                                     ml_HvsC_FungiDecontam_Class$rep_perf,
                                     ml_HvsC_FungiDecontam_Order$rep_perf,
                                     ml_HvsC_FungiDecontam_Family$rep_perf,
                                     ml_HvsC_FungiDecontam_Genus$rep_perf,
                                     ml_HvsC_FungiDecontam_Species$rep_perf),
                                    taxaLevel = c(rep(c("Phylum","Class","Order","Family","Genus","Species"),each=10),
                                                  rep(c("Phylum","Class","Order","Family","Genus","Species"),each=10),
                                                  rep(c("Phylum","Class","Order","Family","Genus","Species"),each=10)))

colnames(ucsdTaxaLevel_results)[1:2] <- c("AUROC","AUPR")
ucsdTaxaLevel_results$nullAUROC <- 0.5
ucsdTaxaLevel_results$nullAUPR <- ucsdTaxaLevel_results$positiveClassSize/
  (ucsdTaxaLevel_results$positiveClassSize+ucsdTaxaLevel_results$negativeClassSize)
ucsdTaxaLevel_results$dataString[grepl("Bacteria",ucsdTaxaLevel_results$dataString)] <- "Bacteria (Raw)"
ucsdTaxaLevel_results$dataString[grepl("Fungi_decontam",ucsdTaxaLevel_results$dataString)] <- "Fungi decontaminated"
ucsdTaxaLevel_results$dataString[grepl("Fungi_",ucsdTaxaLevel_results$dataString)] <- "Fungi (Raw)"
ucsdTaxaLevel_results$dataString <- factor(ucsdTaxaLevel_results$dataString,
                                                levels = c("Bacteria (Raw)",
                                                           "Fungi (Raw)",
                                                           "Fungi decontaminated"))
ucsdTaxaLevel_results$taxaLevel <- factor(ucsdTaxaLevel_results$taxaLevel,
                                               levels = c("Phylum","Class","Order","Family","Genus","Species"))

source("Supporting_scripts/S00-SummarySE.R")
## All data types
ucsdTaxaLevel_results %>%
  # filter(grepl("decontam",dataString)) %>%
  reshape2::melt(id.vars = c("rep","dataString","taxaLevel","positiveClassSize","negativeClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","taxaLevel","dataString","positiveClassSize","negativeClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(taxaLevel,value, color=dataString)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  # guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("UCSD plasma validation cohort:\nPer taxa level cancer vs. healthy") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/ucsd_per_taxa_level_bacteria_and_fungi_4Apr22.svg", 
       dpi = "retina", width = 8, height = 4, units = "in")

## Subset to decontaminated fungi
ucsdTaxaLevel_results %>%
  filter(grepl("decontam",dataString)) %>%
  reshape2::melt(id.vars = c("rep","dataString","taxaLevel","positiveClassSize","negativeClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("variable","taxaLevel","dataString","positiveClassSize","negativeClassSize","diseaseType","col2Predict","nullAUPR","nullAUROC")) %>%
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
  ggtitle("UCSD plasma validation cohort:\nPer taxa level cancer vs. healthy") + theme(plot.title = element_text(hjust = 0.5)) +
  # rotate_x_text(90) + 
  scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/ucsd_per_taxa_level_fungi_decontam_4Apr22.svg", 
       dpi = "retina", width = 6, height = 4, units = "in")

#----------------------------------------------------------------------------------#
# Scrambled and shuffled negative controls
# For healthy vs cancer
#----------------------------------------------------------------------------------#

#-----------------------ML with scrambled data-----------------------#
source("00-Functions.R")
scrambled_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                               snmData = snmDataUCSD_Fungi_Decontam_OGU,
                               col2Predict = "one_cancer_vs_controls",
                               dzOfInterest = "NSCLC",
                               dataString = "scrambled",
                               scrambleFlag=TRUE,
                               varImpFlag = FALSE)
scrambled_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                snmData = snmDataUCSD_Fungi_Decontam_OGU,
                                col2Predict = "one_cancer_vs_controls",
                                dzOfInterest = "PRAD",
                                dataString = "scrambled",
                                scrambleFlag=TRUE,
                                varImpFlag = FALSE)
scrambled_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                snmData = snmDataUCSD_Fungi_Decontam_OGU,
                                col2Predict = "one_cancer_vs_controls",
                                dzOfInterest = "SKCM",
                                dataString = "scrambled",
                                scrambleFlag=TRUE,
                                varImpFlag = FALSE)
#-----------------------ML with shuffled data-----------------------#
shuffled_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                snmData = snmDataUCSD_Fungi_Decontam_OGU,
                                col2Predict = "one_cancer_vs_controls",
                                dzOfInterest = "NSCLC",
                                dataString = "shuffled",
                                shuffleFlag=TRUE,
                                varImpFlag = FALSE)
shuffled_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                               snmData = snmDataUCSD_Fungi_Decontam_OGU,
                               col2Predict = "one_cancer_vs_controls",
                               dzOfInterest = "PRAD",
                               dataString = "shuffled",
                               shuffleFlag=TRUE,
                               varImpFlag = FALSE)
shuffled_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                               snmData = snmDataUCSD_Fungi_Decontam_OGU,
                               col2Predict = "one_cancer_vs_controls",
                               dzOfInterest = "SKCM",
                               dataString = "shuffled",
                               shuffleFlag=TRUE,
                               varImpFlag = FALSE)

#-----------------------Format data-----------------------#

ucsdControls <- rbind(scrambled_NSCLC$rep_perf,
                      scrambled_PRAD$rep_perf,
                      scrambled_SKCM$rep_perf,
                     
                      shuffled_NSCLC$rep_perf,
                      shuffled_PRAD$rep_perf,
                      shuffled_SKCM$rep_perf)

colnames(ucsdControls)[1:2] <- c("AUROC","AUPR")
ucsdControls$nullAUROC <- 0.5
ucsdControls$nullAUPR <- ucsdControls$positiveClassSize/
  (ucsdControls$positiveClassSize+ucsdControls$negativeClassSize)

colnames(ucsdControls)[1:2] <- c("AUROC","AUPR")
ucsdControls$dataString[ucsdControls$dataString=="scrambled"] <- "Scrambled"
ucsdControls$dataString[ucsdControls$dataString=="shuffled"] <- "Shuffled"

ucsdControls$diseaseType <- gsub(" cancer| Cancer","",ucsdControls$diseaseType)

#-----------------------Combine data from actual and controls-----------------------#

ucsdOverlay_HvsC_ActualvsControls <- rbind(ucsdPerCancerVsHealthyResults,
                                           ucsdControls) %>%
  filter(grepl("decontaminated|Top|Scrambled|Shuffled",dataString)) %>% droplevels()
ucsdOverlay_HvsC_ActualvsControls$statGroups <- ifelse(grepl("shuffled|scrambled",
                                                                ucsdOverlay_HvsC_ActualvsControls$dataString),
                                                          yes = "Control", no = "Actual")
ucsdOverlay_HvsC_ActualvsControls$dataString <- as.character(ucsdOverlay_HvsC_ActualvsControls$dataString)
ucsdOverlay_HvsC_ActualvsControls$dataString[grepl("Top",ucsdOverlay_HvsC_ActualvsControls$dataString)] <- "Top 20 fungi"
ucsdOverlay_HvsC_ActualvsControls$dataString[grepl("decontaminated",ucsdOverlay_HvsC_ActualvsControls$dataString)] <- "Decontaminated\nfungi"
ucsdOverlay_HvsC_ActualvsControls$dataString <- factor(ucsdOverlay_HvsC_ActualvsControls$dataString,
                                                       levels = c("Decontaminated\nfungi",
                                                                  "Top 20 fungi",
                                                                  "Scrambled",
                                                                  "Shuffled"))

## Make grouped comparisons (all actual samples vs controls)
require(rstatix)
## AUROC
ucsdOverlay_HvsC_ActualvsControls %>%
  distinct() %>% droplevels() %>%
  group_by(dataString) %>% data.frame() %>%
  wilcox_test(AUROC ~ dataString, exact = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "dataString") -> roc.tumor.stat.test.ucsd
data.frame(roc.tumor.stat.test.ucsd)

ucsdOverlay_HvsC_ActualvsControls %>%
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
  stat_pvalue_manual(roc.tumor.stat.test.ucsd, 
                     # y.position = c(1.03, 1.12, 1.06),
                     y.position = c(1.03, 1.10,1.17, 1.24, 1.31, 1.45),
                     label = "q = {p.adj}", 
                     size = 3) +
  geom_hline(yintercept = 0.5, linetype="dotted") + 
  labs(fill = "Dataset") + theme(plot.title = element_text(hjust=0.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1.5)) -> aggregatedTumorROC
## AUPR
ucsdOverlay_HvsC_ActualvsControls %>%
  distinct() %>% droplevels() %>%
  group_by(dataString) %>% data.frame() %>%
  wilcox_test(AUPR ~ dataString, exact = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "dataString") -> pr.tumor.stat.test.ucsd

ucsdOverlay_HvsC_ActualvsControls %>%
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
  stat_pvalue_manual(pr.tumor.stat.test.ucsd, label = "q = {p.adj}", 
                     # y.position = c(1.03, 1.12, 1.06), 
                     y.position = c(1.03, 1.10,1.17, 1.24, 1.31, 1.45),
                     size = 3) +
  labs(fill = "Dataset") + theme(plot.title = element_text(hjust=0.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1.5)) -> aggregatedTumorPR

combinedAggregatedTumorPlot <- ggarrange(aggregatedTumorROC, aggregatedTumorPR, ncol = 2) 
combinedAggregatedTumorPlotAnnotated <- annotate_figure(combinedAggregatedTumorPlot, 
                                                        top = text_grob("UCSD: Aggregated cancer vs. healthy\nperformance with controls\n(Decontaminated fungi only)", 
                                                                        color = "black", face = "bold", size = 14))
print(combinedAggregatedTumorPlotAnnotated)
ggsave(filename = paste0("Figures/Supplementary_Figures/controls_ucsd_auroc_aupr.svg"),
       plot = combinedAggregatedTumorPlotAnnotated,
       dpi = "retina", units = "in", width = 8, height = 6)

#----------------------------------------------------------------------------------#
# Exploratory immunotherapy response predictions using the mycobiome
#----------------------------------------------------------------------------------#

ucsdMetaIOResponse <- metadataSandip <- read.csv("Input_data/ucsd_plasma_cohort_immunotherapy_response.csv",
                                                 header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
# Reformat data to be compatible with metadataPSMatchedDPQCFiltered
ucsdMetaIOResponse$De.identified.ID <- gsub("T","",ucsdMetaIOResponse$De.identified.ID)
ucsdMetaIOResponse$De.identified.ID <- gsub("-","-0",ucsdMetaIOResponse$De.identified.ID)

# Sanity check
sum(ucsdMetaIOResponse$De.identified.ID %in% metaUCSDJoined_Fungi_Nonzero$tube_id) # 41

# Filter and merge metadata
ucsdMetaIOResponseFilt <- droplevels(ucsdMetaIOResponse[,c("Age.at.Treatment", "Gender",
                                                           "Diagnosis","Responder","De.identified.ID")])
metaUCSDJoined_Fungi_Nonzero_IOR <- left_join(metaUCSDJoined_Fungi_Nonzero,
                                              ucsdMetaIOResponseFilt,
                                              by = c("tube_id" = "De.identified.ID"))
rownames(metaUCSDJoined_Fungi_Nonzero_IOR) <- rownames(metaUCSDJoined_Fungi_Nonzero)
metaUCSDJoined_Fungi_Nonzero_IOR$Responder <- factor(metaUCSDJoined_Fungi_Nonzero_IOR$Responder,
                                                     levels = c("Yes","No", "Indeterminate"))

# Sanity check — look at merged dataframe
droplevels(metaUCSDJoined_Fungi_Nonzero_IOR[metaUCSDJoined_Fungi_Nonzero_IOR$disease_type_consol %in% c("SKCM","NSCLC"),])

# Check row order
all(rownames(metaUCSDJoined_Fungi_Nonzero_IOR) == rownames(ucsd_rep200Data_Filt_Matched)) # TRUE
all(rownames(metaUCSDJoined_Fungi_Nonzero_IOR) == rownames(ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam)) # TRUE

# Subset metadata for WIS intersection and again check row order
metaUCSDJoined_Fungi_Nonzero_SpeciesShared_IOR <- droplevels(metaUCSDJoined_Fungi_Nonzero_IOR[rownames(metaUCSDJoined_Fungi_Nonzero_SpeciesShared),])
ucsdWIS_VSNM_FungiBacteria <- cbind(snmData_ucsdRep200FungiSpeciesShared,
      snmData_ucsdRep200BacteriaSpeciesShared[rownames(snmData_ucsdRep200FungiSpeciesShared),])
all(rownames(metaUCSDJoined_Fungi_Nonzero_SpeciesShared_IOR) == rownames(ucsdWIS_VSNM_FungiBacteria)) # TRUE

# Subset metadata for topX intersection and again check row order
metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero_IOR <- droplevels(metaUCSDJoined_Fungi_Nonzero_IOR[rownames(metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero),])
all(rownames(metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero_IOR) == rownames(metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero)) # TRUE

# Full data: ucsd_rep200Data_Filt_Matched | snmDataUCSD_FullRep200_OGU
# WIS fungi+bacteria: cbind(snmData_ucsdRep200FungiSpeciesShared, snmData_ucsdRep200BacteriaSpeciesShared)
# Fungi decontaminated: ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam | snmDataUCSD_Fungi_Decontam_OGU

#-----------------------Batch correct using responder information-----------------------#
# Subset metadata to samples with responder information
metaUCSDJoined_Fungi_Nonzero_IOR_Filt <- metaUCSDJoined_Fungi_Nonzero_IOR %>% filter(Responder %in% c("Yes","No")) %>% droplevels()
metaUCSDJoined_Fungi_Nonzero_SpeciesShared_IOR_Filt <- metaUCSDJoined_Fungi_Nonzero_SpeciesShared_IOR %>%
  filter(Responder %in% c("Yes","No")) %>% droplevels()
metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero_IOR_Filt <- metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero_IOR %>% 
  filter(Responder %in% c("Yes","No")) %>% droplevels()
# Subset count data
ucsd_rep200Data_Filt_Matched_IOR_Filt <- ucsd_rep200Data_Filt_Matched[rownames(metaUCSDJoined_Fungi_Nonzero_IOR_Filt),]
ucsdRep200BacteriaSpeciesShared_IOR_Filt <- ucsdRep200BacteriaSpeciesShared[rownames(metaUCSDJoined_Fungi_Nonzero_IOR_Filt),]
ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_IOR_Filt <- ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam[rownames(metaUCSDJoined_Fungi_Nonzero_IOR_Filt),]
ucsdRep200FungiSpeciesShared_Nonzero_IOR_Filt <- ucsdRep200FungiSpeciesShared_Nonzero[rownames(metaUCSDJoined_Fungi_Nonzero_SpeciesShared_IOR_Filt),]
ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_TopX_Nonzero_IOR_Filt <- ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_TopX_Nonzero[rownames(metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero_IOR_Filt),]
ucsdRep200BacteriaGenus_IOR_Filt <- ucsdRep200BacteriaGenus[rownames(metaUCSDJoined_Fungi_Nonzero_IOR_Filt),]

source("00-Functions.R") # for vsnmFunctionUCSD_IOR()

# Full rep200 data
snmDataUCSD_FullRep200_OGU_IOR <- vsnmFunctionUCSD_IOR(qcData = ucsd_rep200Data_Filt_Matched_IOR_Filt, 
                                                       qcMetadata = metaUCSDJoined_Fungi_Nonzero_IOR_Filt)
# Bacteria at genus level
snmData_ucsdRep200BacteriaGenus_IOR <- vsnmFunctionUCSD_IOR(qcData = ucsdRep200BacteriaGenus_IOR_Filt, 
                                                      qcMetadata = metaUCSDJoined_Fungi_Nonzero_IOR_Filt)

# Shared bacteria data by taxa levels
snmData_ucsdRep200BacteriaSpeciesShared_IOR <- vsnmFunctionUCSD_IOR(qcData = ucsdRep200BacteriaSpeciesShared_IOR_Filt, 
                                                                    qcMetadata = metaUCSDJoined_Fungi_Nonzero_IOR_Filt)

# Decontaminated fungi data
snmDataUCSD_Fungi_Decontam_OGU_IOR <- vsnmFunctionUCSD_IOR(qcData = ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_IOR_Filt, 
                                                           qcMetadata = metaUCSDJoined_Fungi_Nonzero_IOR_Filt)
# TopX fungi data
snmDataUCSD_Fungi_TopX_IOR <- vsnmFunctionUCSD(qcData = ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_TopX_Nonzero_IOR_Filt, 
                                           qcMetadata = metaUCSDJoined_Fungi_Nonzero_TopX_Nonzero_IOR_Filt)

# Fungi data intersected with Weizmann features
snmData_ucsdRep200FungiSpeciesShared_IOR <- vsnmFunctionUCSD_IOR(qcData = ucsdRep200FungiSpeciesShared_Nonzero_IOR_Filt, 
                                                                 qcMetadata = metaUCSDJoined_Fungi_Nonzero_SpeciesShared_IOR_Filt)

# Join WIS-overlapping fungi+bacterial data
ucsdWIS_VSNM_FungiBacteria_IOR_Filt <- cbind(snmData_ucsdRep200FungiSpeciesShared_IOR,
                                             snmData_ucsdRep200BacteriaSpeciesShared_IOR[rownames(snmData_ucsdRep200FungiSpeciesShared_IOR),])
#-----------------------Machine learning-----------------------#

source("00-Functions.R")
loocvIOR_FungiShared_SKCM <- loocvIOR(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared_IOR_Filt,
                                      snmData = snmData_ucsdRep200FungiSpeciesShared_IOR, 
                                      dzType = "SKCM",
                                      savePlotFlag = TRUE,
                                      dataStringForPlotFilename = "WIS_Shared_Fungi")
# AUROC: 0.71 | AUPR: 0.71
