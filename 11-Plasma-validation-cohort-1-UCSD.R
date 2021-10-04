#-----------------------------------------------------------------------------
# 11-Plasma-validation-cohort-1-UCSD.R
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
metaPoore <- read.csv("Input_data/qiita_metadata_plasma_validation_cohort_1.txt", sep = "\t", stringsAsFactors = FALSE)
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

ucsd_rep200Data_BIOM <- read_biom(biom_file = "Input_data/Qiita_results/rep200_OGU_plasma_validation_cohort_1.biom")
ucsd_rep200Data <- t(as(biom_data(ucsd_rep200Data_BIOM), "matrix"))

sum(rownames(ucsd_rep200Data) %in% rownames(metaUCSDJoined)) # 169
sum(rownames(metaUCSDJoined)%in% rownames(ucsd_rep200Data)) # 169

ucsd_rep200Data_Filt <- ucsd_rep200Data[rownames(metaUCSDJoined),]
dim(ucsd_rep200Data_Filt) # 169 7777

## Extract only fungal and bacterial features
ucsd_rep200Data_Filt_Fungi <- ucsd_rep200Data_Filt[,colnames(ucsd_rep200Data_Filt) %in% fungiOGUs]
ucsd_rep200Data_Filt_Bacteria <- ucsd_rep200Data_Filt[,colnames(ucsd_rep200Data_Filt) %in% bacteriaOGUs]

## Remove samples with 0 counts after filtering to Fungi: "12667.X2051123" "12691.PC22"     "12691.PC47"
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

# Build phyloseq object
psUCSDBacteria <- phyloseq(otu_table(ucsd_rep200Data_Filt_Bacteria_Matched, taxa_are_rows = FALSE), 
                        tax_table(as.matrix(rep200TaxSplit_Bacteria)), 
                        sample_data(metaUCSDJoined_Bacteria_Matched))

## Aggregate counts
psUCSDBacteria_phylum = aggregate_taxa(psUCSDBacteria, "Phylum")
psUCSDBacteria_class = aggregate_taxa(psUCSDBacteria, "Class")
psUCSDBacteria_order = aggregate_taxa(psUCSDBacteria, "Order")
psUCSDBacteria_family = aggregate_taxa(psUCSDBacteria, "Family")
psUCSDBacteria_genus = aggregate_taxa(psUCSDBacteria, "Genus")
psUCSDBacteria_species = ucsd_rep200Data_Filt_Bacteria_Matched
colnames(psUCSDBacteria_species) <- rep200TaxSplit_Bacteria[colnames(psUCSDBacteria_species),"Species"]

## Create data.frames of summarized data
ucsdRep200BacteriaPhylum <- data.frame(t(otu_table(psUCSDBacteria_phylum)))
ucsdRep200BacteriaClass <- data.frame(t(otu_table(psUCSDBacteria_class)))
ucsdRep200BacteriaOrder <- data.frame(t(otu_table(psUCSDBacteria_order)))
ucsdRep200BacteriaFamily <- data.frame(t(otu_table(psUCSDBacteria_family)))
ucsdRep200BacteriaGenus <- data.frame(t(otu_table(psUCSDBacteria_genus)))
ucsdRep200BacteriaSpecies <- data.frame(psUCSDBacteria_species)
ucsdRep200BacteriaSpecies[1:3,1:3]

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
load("Interim_data/shared_fungi_features_at_each_taxa_level_13Sep21.RData")

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
# Decontaminate using TCGA and Nature Val data
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
# decontamination, although the others would be equally suitable too.

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
# **From the decontam tutorial:**
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
sum(contamSumFreq)/sum(colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero_freq))) #--> 0.002012228

#--------------------------Combine prevalence and frequency contaminants and calculate % read removal--------------------------#
contaminantsFungiUCSD <- unique(c(rownames(contamdf.prev.fungi.ucsd[which(contamdf.prev.fungi.ucsd$contaminant),]),
                                       rownames(contamdf.freq.fungi.ucsd[which(contamdf.freq.fungi.ucsd$contaminant),])))
length(contaminantsFungiUCSD) # 32

# Calculate total % read count removal using prevalence and frequency decontamination
notContamSumTotal <- colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero)[,!(colnames(ucsd_rep200Data_decontam_dedup_fungi_nonzero) %in% contaminantsFungiUCSD)])
contamSumTotal <- colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero)[,contaminantsFungiUCSD])
sum(contamSumTotal)/sum(colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero))) #--> 0.08948648

#--------------------------Combine decontam and TCGA contaminants--------------------------#
load("Interim_data/contamdf.freq.fungi.plateCenter_13Sep21.RData") # loads "contamdf.freq.fungi.plateCenter" object
contaminantsFungiUCSDAndPlateCenterTCGA <- unique(c(contaminantsFungiUCSD,
                                  rownames(contamdf.freq.fungi.plateCenter[which(contamdf.freq.fungi.plateCenter$contaminant),])))
save(contaminantsFungiUCSDAndPlateCenterTCGA, 
     contaminantsFungiUCSD,
     file = "Interim_data/contaminants_fungi_OGUs_UCSD_TCGA_25Sep21.RData")

# Calculate total % read count removal using decontam and TCGA decontamination
notContamSumTotalWithTCGA <- colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero)[,!(colnames(ucsd_rep200Data_decontam_dedup_fungi_nonzero) %in% contaminantsFungiUCSDAndPlateCenterTCGA)])
contamSumTotalWithTCGA <- colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero)[,colnames(ucsd_rep200Data_decontam_dedup_fungi_nonzero) %in% contaminantsFungiUCSDAndPlateCenterTCGA])
sum(contamSumTotalWithTCGA)/sum(colSums(as.matrix(ucsd_rep200Data_decontam_dedup_fungi_nonzero))) #--> 0.1892535

#--------------------------Remove contaminants from UCSD plasma cohort--------------------------#

ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam <- ucsd_rep200Data_Filt_Fungi_Nonzero[,!(colnames(ucsd_rep200Data_Filt_Fungi_Nonzero) %in% contaminantsFungiUCSDAndPlateCenterTCGA)]
dim(ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam) # 166 145

# There is 1 zero sum sample after decontamination (12691.PC48), so it should be removed
zeroSumFungiSamplesPostDecontam <- names(which(rowSums(ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam) == 0))
ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_Final <- ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam[!(rownames(ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam) %in% zeroSumFungiSamplesPostDecontam),]
metaUCSDJoined_Fungi_Nonzero_Final <- droplevels(metaUCSDJoined_Fungi_Nonzero[!(rownames(metaUCSDJoined_Fungi_Nonzero) %in% zeroSumFungiSamplesPostDecontam),])
all(rownames(metaUCSDJoined_Fungi_Nonzero_Final) == rownames(ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_Final)) # TRUE
dim(ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_Final) # 165  145

#-----------------------------------------------
#           VSNM Normalization                 #
#-----------------------------------------------

source("00-Functions.R") # for vsnmFunctionUCSD()

# Full rep200 data
snmDataUCSD_FullRep200_OGU <- vsnmFunctionUCSD(qcData = ucsd_rep200Data_Filt_Matched, qcMetadata = metaUCSDJoined_FullRep200_Matched)

# Bacteria data by taxa levels
snmData_ucsdRep200BacteriaPhylum <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaPhylum, qcMetadata = metaUCSDJoined_Bacteria_Matched)
snmData_ucsdRep200BacteriaClass <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaClass, qcMetadata = metaUCSDJoined_Bacteria_Matched)
snmData_ucsdRep200BacteriaOrder <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaOrder, qcMetadata = metaUCSDJoined_Bacteria_Matched)
snmData_ucsdRep200BacteriaFamily <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaFamily, qcMetadata = metaUCSDJoined_Bacteria_Matched)
snmData_ucsdRep200BacteriaGenus <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaGenus, qcMetadata = metaUCSDJoined_Bacteria_Matched)
snmData_ucsdRep200BacteriaSpecies <- vsnmFunctionUCSD(qcData = ucsdRep200BacteriaSpecies, qcMetadata = metaUCSDJoined_Bacteria_Matched)
snmDataUCSD_Bacteria_OGU <- vsnmFunctionUCSD(qcData = ucsd_rep200Data_Filt_Bacteria_Matched, qcMetadata = metaUCSDJoined_Bacteria_Matched)

# Fungi data by taxa levels
snmData_ucsdRep200FungiPhylum <- vsnmFunctionUCSD(qcData = ucsdRep200FungiPhylum, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiClass <- vsnmFunctionUCSD(qcData = ucsdRep200FungiClass, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiOrder <- vsnmFunctionUCSD(qcData = ucsdRep200FungiOrder, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiFamily <- vsnmFunctionUCSD(qcData = ucsdRep200FungiFamily, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiGenus <- vsnmFunctionUCSD(qcData = ucsdRep200FungiGenus, qcMetadata = metaUCSDJoined_Fungi_Nonzero)
snmData_ucsdRep200FungiSpecies <- vsnmFunctionUCSD(qcData = ucsdRep200FungiSpecies, qcMetadata = metaUCSDJoined_Fungi_Nonzero)

# Decontaminated fungi data
snmDataUCSD_Fungi_Decontam_OGU <- vsnmFunctionUCSD(qcData = ucsd_rep200Data_Filt_Fungi_Nonzero_Decontam_Final, qcMetadata = metaUCSDJoined_Fungi_Nonzero_Final)

# Fungi data intersected with Weizmann features
# NOTE: The following taxa levels did not converge and are not shown: Phylum
# snmData_ucsdRep200FungiPhylumShared <- vsnmFunctionUCSD(qcData = ucsdRep200FungiPhylumShared_Nonzero, qcMetadata = metaUCSDJoined_Fungi_Nonzero_PhylumShared)
snmData_ucsdRep200FungiClassShared <- vsnmFunctionUCSD(qcData = ucsdRep200FungiClassShared_Nonzero, qcMetadata = metaUCSDJoined_Fungi_Nonzero_ClassShared)
snmData_ucsdRep200FungiOrderShared <- vsnmFunctionUCSD(qcData = ucsdRep200FungiOrderShared_Nonzero, qcMetadata = metaUCSDJoined_Fungi_Nonzero_OrderShared)
snmData_ucsdRep200FungiFamilyShared <- vsnmFunctionUCSD(qcData = ucsdRep200FungiFamilyShared_Nonzero, qcMetadata = metaUCSDJoined_Fungi_Nonzero_FamilyShared)
snmData_ucsdRep200FungiGenusShared <- vsnmFunctionUCSD(qcData = ucsdRep200FungiGenusShared_Nonzero, qcMetadata = metaUCSDJoined_Fungi_Nonzero_GenusShared)
snmData_ucsdRep200FungiSpeciesShared <- vsnmFunctionUCSD(qcData = ucsdRep200FungiSpeciesShared_Nonzero, qcMetadata = metaUCSDJoined_Fungi_Nonzero_SpeciesShared)

#----------------------------------------------------------------------------------#
# Machine learning: Aggregated cancer vs. Healthy @ OGU/species level
# NOTE: For fungal data, there is only 1 OGU per species
#----------------------------------------------------------------------------------#

source("00-Functions.R") # for ml1VsAllUCSD() function
metaUCSDJoined_Fungi_Nonzero_Final %>% count(HvsC)
# HvsC  n
# 1  Cancer 97
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
                                   dataString = "bacteria",
                                   varImpFlag = FALSE)

#---------------------------------Fungi---------------------------------#
ml_HvsC_Fungi_Species <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                           snmData = snmData_ucsdRep200FungiSpecies,
                                           col2Predict = "HvsC",
                                           dataString = "fungi",
                                           varImpFlag = FALSE)

#---------------------------------Decontam Fungi---------------------------------#
ml_HvsC_Decontam_Fungi_OGU <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_Final,
                                   snmData = snmDataUCSD_Fungi_Decontam_OGU,
                                   col2Predict = "HvsC",
                                   dataString = "decontam_fungi",
                                   varImpFlag = FALSE)

#---------------------------------Intersected fungi---------------------------------#
ml_HvsC_Shared_Fungi_Species <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                   snmData = snmData_ucsdRep200FungiSpeciesShared,
                                   col2Predict = "HvsC",
                                   dataString = "intersected_fungi",
                                   varImpFlag = FALSE)

#---------------------------------Combine results---------------------------------#
ucsdCancerVsHealthyResults <- rbind(ml_HvsC_Fullrep200_OGU$rep_perf,
                                    ml_HvsC_Bacteria_OGU$rep_perf,
                                    ml_HvsC_Fungi_Species$rep_perf,
                                    ml_HvsC_Decontam_Fungi_OGU$rep_perf,
                                    ml_HvsC_Shared_Fungi_Species$rep_perf)

colnames(ucsdCancerVsHealthyResults)[1:2] <- c("AUROC","AUPR")
ucsdCancerVsHealthyResults$dataString[ucsdCancerVsHealthyResults$dataString=="full_rep200"] <- "Full dataset (rep200)"
ucsdCancerVsHealthyResults$dataString[ucsdCancerVsHealthyResults$dataString=="bacteria"] <- "Bacteria all (rep200)"
ucsdCancerVsHealthyResults$dataString[ucsdCancerVsHealthyResults$dataString=="fungi"] <- "Fungi all (rep200)"
ucsdCancerVsHealthyResults$dataString[ucsdCancerVsHealthyResults$dataString=="decontam_fungi"] <- "Fungi decontaminated (rep200)"
ucsdCancerVsHealthyResults$dataString[ucsdCancerVsHealthyResults$dataString=="intersected_fungi"] <- "Intersected fungi only (rep200+Weizmann)"
ucsdCancerVsHealthyResults$dataString <- factor(ucsdCancerVsHealthyResults$dataString, levels = c("Full dataset (rep200)",
                                                                                                  "Bacteria all (rep200)",
                                                                                                  "Fungi all (rep200)",
                                                                                                  "Fungi decontaminated (rep200)",
                                                                                                  "Intersected fungi only (rep200+Weizmann)"))
require(ggrepel)
source("Supporting_scripts/S05-SummarySE.R")
# Facet by AUROC vs AUPR; color by full dataset vs fungi (dataString)
ucsdCancerVsHealthyResults %>%
  reshape2::melt(id.vars = c("dataString","rep","diseaseType", "col2Predict")) %>%
  summarySE(measurevar = "value", groupvars = c("diseaseType","dataString","variable")) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=dataString)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=3) + 
  ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + ggtitle("Plasma validation cohort #1 cancer vs. healthy (97 cancer | 68 healthy)") +
  theme(plot.title = element_text(hjust = 0), legend.position = "right") + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "Cancer vs. Healthy") +
  rotate_x_text(0) + scale_color_nejm(name = "Dataset") + geom_hline(yintercept = 1, linetype="dashed") + 
  geom_label_repel(aes(label=round(value,2)), size=3, box.padding = 0.75, point.padding = 6.5, direction = "y", position = position_dodge(width = 1), show.legend = FALSE) +
  ggsave("Figures/Figure_5/figure_5D_ucsd_all_cancer_vs_healthy_all_datasets_25Sep21.jpeg", 
         dpi = "retina", width = 8, height = 4.5, units = "in")
ucsdCancerVsHealthyResults %>% write.csv("Figures_data/Figure_5/figure_5_D_ucsd_all_cancer_vs_healthy_all_datasets_25Sep21.csv")


#----------------------------------------------------------------------------------#
# Machine learning: Individual cancer vs. Healthy
#----------------------------------------------------------------------------------#

metaUCSDJoined_Fungi_Nonzero_Final %>% count(disease_type_consol)
# disease_type_consol  n
# 1             Control 68
# 2               NSCLC 25
# 3                PRAD 56
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
                                     dataString = "bacteria",
                                     varImpFlag = FALSE)
ml_HvsC_Bacteria_OGU_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                     snmData = snmDataUCSD_Bacteria_OGU,
                                     col2Predict = "one_cancer_vs_controls",
                                     dzOfInterest = "PRAD",
                                     dataString = "bacteria",
                                     varImpFlag = FALSE)
ml_HvsC_Bacteria_OGU_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Bacteria_Matched,
                                     snmData = snmDataUCSD_Bacteria_OGU,
                                     col2Predict = "one_cancer_vs_controls",
                                     dzOfInterest = "SKCM",
                                     dataString = "bacteria",
                                     varImpFlag = FALSE)

#---------------------------------Fungi---------------------------------#
ml_HvsC_Fungi_Species_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                                 snmData = snmData_ucsdRep200FungiSpecies,
                                                 col2Predict = "one_cancer_vs_controls",
                                                 dzOfInterest = "NSCLC",
                                                 dataString = "fungi",
                                                 varImpFlag = FALSE)
ml_HvsC_Fungi_Species_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                                snmData = snmData_ucsdRep200FungiSpecies,
                                                col2Predict = "one_cancer_vs_controls",
                                                dzOfInterest = "PRAD",
                                                dataString = "fungi",
                                                varImpFlag = FALSE)
ml_HvsC_Fungi_Species_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                                snmData = snmData_ucsdRep200FungiSpecies,
                                                col2Predict = "one_cancer_vs_controls",
                                                dzOfInterest = "SKCM",
                                                dataString = "fungi",
                                                varImpFlag = FALSE)

#---------------------------------Decontam Fungi---------------------------------#
ml_HvsC_Decontam_Fungi_OGU_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_Final,
                                           snmData = snmDataUCSD_Fungi_Decontam_OGU,
                                           col2Predict = "one_cancer_vs_controls",
                                           dzOfInterest = "NSCLC",
                                           dataString = "decontam_fungi",
                                           varImpFlag = FALSE)
ml_HvsC_Decontam_Fungi_OGU_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_Final,
                                           snmData = snmDataUCSD_Fungi_Decontam_OGU,
                                           col2Predict = "one_cancer_vs_controls",
                                           dzOfInterest = "PRAD",
                                           dataString = "decontam_fungi",
                                           varImpFlag = FALSE)
ml_HvsC_Decontam_Fungi_OGU_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_Final,
                                           snmData = snmDataUCSD_Fungi_Decontam_OGU,
                                           col2Predict = "one_cancer_vs_controls",
                                           dzOfInterest = "SKCM",
                                           dataString = "decontam_fungi",
                                           varImpFlag = FALSE)

#---------------------------------Intersected fungi---------------------------------#
ml_HvsC_Shared_Fungi_Species_NSCLC <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                             snmData = snmData_ucsdRep200FungiSpeciesShared,
                                             col2Predict = "one_cancer_vs_controls",
                                             dzOfInterest = "NSCLC",
                                             dataString = "intersected_fungi",
                                             varImpFlag = FALSE)
ml_HvsC_Shared_Fungi_Species_PRAD <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                             snmData = snmData_ucsdRep200FungiSpeciesShared,
                                             col2Predict = "one_cancer_vs_controls",
                                             dzOfInterest = "PRAD",
                                             dataString = "intersected_fungi",
                                             varImpFlag = FALSE)
ml_HvsC_Shared_Fungi_Species_SKCM <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero_SpeciesShared,
                                             snmData = snmData_ucsdRep200FungiSpeciesShared,
                                             col2Predict = "one_cancer_vs_controls",
                                             dzOfInterest = "SKCM",
                                             dataString = "intersected_fungi",
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
                                    ml_HvsC_Shared_Fungi_Species_SKCM$rep_perf)

colnames(ucsdPerCancerVsHealthyResults)[1:2] <- c("AUROC","AUPR")
ucsdPerCancerVsHealthyResults$dataString[ucsdPerCancerVsHealthyResults$dataString=="full_rep200"] <- "Full dataset (rep200)"
ucsdPerCancerVsHealthyResults$dataString[ucsdPerCancerVsHealthyResults$dataString=="bacteria"] <- "Bacteria all (rep200)"
ucsdPerCancerVsHealthyResults$dataString[ucsdPerCancerVsHealthyResults$dataString=="fungi"] <- "Fungi all (rep200)"
ucsdPerCancerVsHealthyResults$dataString[ucsdPerCancerVsHealthyResults$dataString=="decontam_fungi"] <- "Fungi decontaminated (rep200)"
ucsdPerCancerVsHealthyResults$dataString[ucsdPerCancerVsHealthyResults$dataString=="intersected_fungi"] <- "Intersected fungi only (rep200+Weizmann)"
ucsdPerCancerVsHealthyResults$dataString <- factor(ucsdPerCancerVsHealthyResults$dataString, levels = c("Full dataset (rep200)",
                                                                                                  "Bacteria all (rep200)",
                                                                                                  "Fungi all (rep200)",
                                                                                                  "Fungi decontaminated (rep200)",
                                                                                                  "Intersected fungi only (rep200+Weizmann)"))

source("Supporting_scripts/S05-SummarySE.R")
# Facet by AUROC vs AUPR; color by full dataset vs fungi (dataString)
ucsdPerCancerVsHealthyResults %>%
  reshape2::melt(id.vars = c("dataString","rep","diseaseType", "col2Predict")) %>%
  summarySE(measurevar = "value", groupvars = c("diseaseType","dataString","variable")) %>%
  ggplot(aes(reorder(diseaseType, value, FUN=median),value, color=dataString)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=3) + xlab("Individual cancer type versus all healthy samples") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + ggtitle("Plasma validation cohort #1: Per cancer type vs. healthy performance") +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  rotate_x_text(0) + scale_color_nejm(name = "Dataset") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_5/figure_5_E_ucsd_per_cancer_type_vs_healthy_all_datasets_25Sep21.png", 
         dpi = "retina", width = 12, height = 4.5, units = "in")
ucsdPerCancerVsHealthyResults %>% write.csv("Figures_data/Figure_5/figure_5_E_ucsd_per_cancer_type_vs_healthy_all_datasets_25Sep21.csv")

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
                                        dataString = "Fungi",
                                        varImpFlag = FALSE)
ml_HvsC_Fungi_Class <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                       snmData = snmData_ucsdRep200FungiClass,
                                       col2Predict = "HvsC",
                                       dataString = "Fungi",
                                       varImpFlag = FALSE)
ml_HvsC_Fungi_Order <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                       snmData = snmData_ucsdRep200FungiOrder,
                                       col2Predict = "HvsC",
                                       dataString = "Fungi",
                                       varImpFlag = FALSE)
ml_HvsC_Fungi_Family <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                        snmData = snmData_ucsdRep200FungiFamily,
                                        col2Predict = "HvsC",
                                        dataString = "Fungi",
                                        varImpFlag = FALSE)
ml_HvsC_Fungi_Genus <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                       snmData = snmData_ucsdRep200FungiGenus,
                                       col2Predict = "HvsC",
                                       dataString = "Fungi",
                                       varImpFlag = FALSE)
ml_HvsC_Fungi_Species <- ml1VsAllUCSD(metaData = metaUCSDJoined_Fungi_Nonzero,
                                         snmData = snmData_ucsdRep200FungiSpecies,
                                         col2Predict = "HvsC",
                                         dataString = "Fungi",
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
                                     ml_HvsC_Fungi_Species$rep_perf),
                                    taxaLevel = c(rep(c("Phylum","Class","Order","Family","Genus","Species"),each=10),
                                                  rep(c("Phylum","Class","Order","Family","Genus","Species"),each=10)))
colnames(ucsdTaxaLevel_results)[1:2] <- c("AUROC","AUPR")
ucsdTaxaLevel_results$dataString <- factor(ucsdTaxaLevel_results$dataString, levels = c("Bacteria","Fungi"))
ucsdTaxaLevel_results$taxaLevel <- factor(ucsdTaxaLevel_results$taxaLevel,
                                               levels = c("Phylum","Class","Order","Family","Genus","Species"))

source("Supporting_scripts/S05-SummarySE.R")
# Facet by AUROC vs AUPR; color by full dataset vs fungi (dataString)
ucsdTaxaLevel_results %>%
  reshape2::melt(id.vars = c("dataString","rep","taxaLevel","diseaseType", "col2Predict")) %>%
  summarySE(measurevar = "value", groupvars = c("taxaLevel","diseaseType","dataString","variable")) %>%
  ggplot(aes(taxaLevel,value, color=dataString)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=3) + xlab("") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + ggtitle("Plasma validation cohort #1: All cancer (n=97) vs. healthy (n=68) at varying taxa levels") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "top") +
  rotate_x_text(0) + scale_color_nejm(name = "Dataset") + geom_hline(yintercept = 1, linetype="dashed") + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  # geom_label_repel(aes(label=round(value,2)), size=3, box.padding = 1, point.padding = 10, position = position_dodge(width = 1), show.legend = FALSE) +
  ggsave("Figures/Supplementary_Figures/ucsd_all_cancer_vs_healthy_bacteria_vs_fungi_varying_taxa_levels_25Sep21.jpeg", dpi = "retina",
         width = 10, height = 4.5, units = "in")
ucsdTaxaLevel_results %>% write.csv("Figures_data/Supplementary_Figures/ucsd_all_cancer_vs_healthy_bacteria_vs_fungi_varying_taxa_levels_25Sep21.csv")



