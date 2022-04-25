#-----------------------------------------------------------------------------
# 05B-Prepare-TCGA-data-for-machine-learning.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Aggregate counts to higher taxa levels
# - Prepare TCGA data for machine learning analyses
# - Plot machine learning results
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
# Load TCGA fungi data and convert to phyloseq object
#----------------------------------------------------------#

load("Interim_data/snmDataFungi_DecontamV2_25Mar22.RData", verbose = T)
rep200TaxSplit <- read.csv("Supporting_data/rep200_lineage_map_split.csv", stringsAsFactors = FALSE, row.names = 1)

#--------------------------------------------------------------------------------------------------------------------#
# Separate raw data into seq center-experimental strategy groups (to preclude needing batch correction)
#--------------------------------------------------------------------------------------------------------------------#
metaQiitaCombined_Nonzero_DecontamV2 %>% count(data_submitting_center_label, experimental_strategy)

#--------------------Subset metadata--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaQiitaCombined_Nonzero_DecontamV2_HMS <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(data_submitting_center_label == "Harvard Medical School") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_BCM <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(data_submitting_center_label == "Baylor College of Medicine") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_MDA <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_WashU <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(data_submitting_center_label == "Washington University School of Medicine") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS <- metaQiitaCombined_Nonzero_DecontamV2 %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
metaQiitaCombined_Nonzero_DecontamV2_UNC <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(data_submitting_center_label == "University of North Carolina") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_CMS <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA <- metaQiitaCombined_Nonzero_DecontamV2 %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "RNA-Seq") %>% droplevels()

#--------------------Subset count data--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
rep200_HiSeq_Fungi_DecontamV2_HMS <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_HMS),]
rep200_HiSeq_Fungi_DecontamV2_BCM <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_BCM),]
rep200_HiSeq_Fungi_DecontamV2_MDA <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_MDA),]
rep200_HiSeq_Fungi_DecontamV2_WashU <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_WashU),]
rep200_HiSeq_Fungi_DecontamV2_Broad_WGS <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
rep200_HiSeq_Fungi_DecontamV2_UNC <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_UNC),]
rep200_HiSeq_Fungi_DecontamV2_CMS <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_CMS),]
rep200_HiSeq_Fungi_DecontamV2_Broad_RNA <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA),]
#--------------------Intersect with fungi with >1% aggregate genome coverage--------------------#
# metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts, 
# rep200Data_WGS_RNA_Matched,
# rep200Data_WGS_RNA_Matched_Bacteria,
# rep200Data_WGS_RNA_Matched_Fungi,
load("Interim_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_25Mar22.RData")
coverageFungiAllSamples <- read.csv("Input_data/fungi_filt_updated_29Sep21_coverage_all_wgs_and_rna_output.csv", stringsAsFactors = FALSE)
coverageFungiAllSamples_1Percent <- coverageFungiAllSamples %>% filter(coverage_ratio >= 0.01)
coverageFungiAllSamples_1Percent_OGUs <- coverageFungiAllSamples_1Percent$gotu
# WGS
rep200_HiSeq_Fungi_HMS_Cov <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_HMS),colnames(rep200Data_WGS_RNA_Matched_Fungi) %in% coverageFungiAllSamples_1Percent_OGUs]
rep200_HiSeq_Fungi_BCM_Cov <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_BCM),colnames(rep200Data_WGS_RNA_Matched_Fungi) %in% coverageFungiAllSamples_1Percent_OGUs]
rep200_HiSeq_Fungi_MDA_Cov <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_MDA),colnames(rep200Data_WGS_RNA_Matched_Fungi) %in% coverageFungiAllSamples_1Percent_OGUs]
rep200_HiSeq_Fungi_WashU_Cov <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_WashU),colnames(rep200Data_WGS_RNA_Matched_Fungi) %in% coverageFungiAllSamples_1Percent_OGUs]
rep200_HiSeq_Fungi_Broad_WGS_Cov <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS),colnames(rep200Data_WGS_RNA_Matched_Fungi) %in% coverageFungiAllSamples_1Percent_OGUs]
# RNA
rep200_HiSeq_Fungi_UNC_Cov <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_UNC),colnames(rep200Data_WGS_RNA_Matched_Fungi) %in% coverageFungiAllSamples_1Percent_OGUs]
rep200_HiSeq_Fungi_CMS_Cov <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_CMS),colnames(rep200Data_WGS_RNA_Matched_Fungi) %in% coverageFungiAllSamples_1Percent_OGUs]
rep200_HiSeq_Fungi_Broad_RNA_Cov <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA),colnames(rep200Data_WGS_RNA_Matched_Fungi) %in% coverageFungiAllSamples_1Percent_OGUs]

## Remove zero sum samples -- Broad_WGS and UNC have zero sum samples; for naming convention, all will be processed
# e.g. ummary(rowSums(rep200_HiSeq_Fungi_Broad_WGS_Cov)==0)
rep200_HiSeq_Fungi_HMS_Cov_Nonzero <- rep200_HiSeq_Fungi_HMS_Cov[!rowSums(rep200_HiSeq_Fungi_HMS_Cov)==0,]
rep200_HiSeq_Fungi_BCM_Cov_Nonzero <- rep200_HiSeq_Fungi_BCM_Cov[!rowSums(rep200_HiSeq_Fungi_BCM_Cov)==0,]
rep200_HiSeq_Fungi_MDA_Cov_Nonzero <- rep200_HiSeq_Fungi_MDA_Cov[!rowSums(rep200_HiSeq_Fungi_MDA_Cov)==0,]
rep200_HiSeq_Fungi_WashU_Cov_Nonzero <- rep200_HiSeq_Fungi_WashU_Cov[!rowSums(rep200_HiSeq_Fungi_WashU_Cov)==0,]
rep200_HiSeq_Fungi_Broad_WGS_Cov_Nonzero <- rep200_HiSeq_Fungi_Broad_WGS_Cov[!rowSums(rep200_HiSeq_Fungi_Broad_WGS_Cov)==0,]
rep200_HiSeq_Fungi_UNC_Cov_Nonzero <- rep200_HiSeq_Fungi_UNC_Cov[!rowSums(rep200_HiSeq_Fungi_UNC_Cov)==0,]
rep200_HiSeq_Fungi_CMS_Cov_Nonzero <- rep200_HiSeq_Fungi_CMS_Cov[!rowSums(rep200_HiSeq_Fungi_CMS_Cov)==0,]
rep200_HiSeq_Fungi_Broad_RNA_Cov_Nonzero <- rep200_HiSeq_Fungi_Broad_RNA_Cov[!rowSums(rep200_HiSeq_Fungi_Broad_RNA_Cov)==0,]
# Match metadata
metaQiitaCombined_Nonzero_DecontamV2_HMS_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_HMS[rownames(rep200_HiSeq_Fungi_HMS_Cov_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_BCM_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_BCM[rownames(rep200_HiSeq_Fungi_BCM_Cov_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_MDA_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_MDA[rownames(rep200_HiSeq_Fungi_MDA_Cov_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_WashU_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_WashU[rownames(rep200_HiSeq_Fungi_WashU_Cov_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS[rownames(rep200_HiSeq_Fungi_Broad_WGS_Cov_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_UNC_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_UNC[rownames(rep200_HiSeq_Fungi_UNC_Cov_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_CMS_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_CMS[rownames(rep200_HiSeq_Fungi_CMS_Cov_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA[rownames(rep200_HiSeq_Fungi_Broad_RNA_Cov_Nonzero),])

#--------------------Save data for ML--------------------#

save(rep200_HiSeq_Fungi_DecontamV2_HMS,
     rep200_HiSeq_Fungi_DecontamV2_BCM,
     rep200_HiSeq_Fungi_DecontamV2_MDA,
     rep200_HiSeq_Fungi_DecontamV2_WashU,
     rep200_HiSeq_Fungi_DecontamV2_Broad_WGS,
     rep200_HiSeq_Fungi_DecontamV2_UNC,
     rep200_HiSeq_Fungi_DecontamV2_CMS,
     rep200_HiSeq_Fungi_DecontamV2_Broad_RNA,
     
     metaQiitaCombined_Nonzero_DecontamV2_HMS,
     metaQiitaCombined_Nonzero_DecontamV2_BCM,
     metaQiitaCombined_Nonzero_DecontamV2_MDA,
     metaQiitaCombined_Nonzero_DecontamV2_WashU,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS,
     metaQiitaCombined_Nonzero_DecontamV2_UNC,
     metaQiitaCombined_Nonzero_DecontamV2_CMS,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA,
     
     rep200_HiSeq_Fungi_HMS_Cov_Nonzero,
     rep200_HiSeq_Fungi_BCM_Cov_Nonzero,
     rep200_HiSeq_Fungi_MDA_Cov_Nonzero,
     rep200_HiSeq_Fungi_WashU_Cov_Nonzero,
     rep200_HiSeq_Fungi_Broad_WGS_Cov_Nonzero,
     rep200_HiSeq_Fungi_UNC_Cov_Nonzero,
     rep200_HiSeq_Fungi_CMS_Cov_Nonzero,
     rep200_HiSeq_Fungi_Broad_RNA_Cov_Nonzero,
     
     metaQiitaCombined_Nonzero_DecontamV2_HMS_Cov_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_BCM_Cov_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_MDA_Cov_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_WashU_Cov_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_Cov_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_UNC_Cov_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_CMS_Cov_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_Cov_Nonzero,
     file = "Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy_with_coverage_filter_decontamV2_2Apr22.RData")

# Scripts: S03R

#--------------------Aggregate count data at taxa levels--------------------#

rep200TaxSplit_Fungi_Paired_to_Weizmann <- read.csv("Supporting_data/rep200TaxSplit_Fungi_Paired_To_Weizmann_Final.csv", stringsAsFactors = FALSE, row.names = 1)

# Build phyloseq object
psRep200_HiSeq_Fungi_DecontamV2_HMS <- phyloseq(otu_table(rep200_HiSeq_Fungi_DecontamV2_HMS, taxa_are_rows = FALSE),
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_HMS))
psRep200_HiSeq_Fungi_DecontamV2_BCM <- phyloseq(otu_table(rep200_HiSeq_Fungi_DecontamV2_BCM, taxa_are_rows = FALSE),
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_BCM))
psRep200_HiSeq_Fungi_DecontamV2_MDA <- phyloseq(otu_table(rep200_HiSeq_Fungi_DecontamV2_MDA, taxa_are_rows = FALSE),
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_MDA))
psRep200_HiSeq_Fungi_DecontamV2_WashU <- phyloseq(otu_table(rep200_HiSeq_Fungi_DecontamV2_WashU, taxa_are_rows = FALSE),
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_WashU))
psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS <- phyloseq(otu_table(rep200_HiSeq_Fungi_DecontamV2_Broad_WGS, taxa_are_rows = FALSE),
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS))
psRep200_HiSeq_Fungi_DecontamV2_UNC <- phyloseq(otu_table(rep200_HiSeq_Fungi_DecontamV2_UNC, taxa_are_rows = FALSE),
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_UNC))
psRep200_HiSeq_Fungi_DecontamV2_CMS <- phyloseq(otu_table(rep200_HiSeq_Fungi_DecontamV2_CMS, taxa_are_rows = FALSE),
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_CMS))
psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA <- phyloseq(otu_table(rep200_HiSeq_Fungi_DecontamV2_Broad_RNA, taxa_are_rows = FALSE),
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA))

## Aggregate counts - HMS
psRep200_HiSeq_Fungi_DecontamV2_HMS_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_HMS, "phylum")
psRep200_HiSeq_Fungi_DecontamV2_HMS_class = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_HMS, "class")
psRep200_HiSeq_Fungi_DecontamV2_HMS_order = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_HMS, "order")
psRep200_HiSeq_Fungi_DecontamV2_HMS_family = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_HMS, "family")
psRep200_HiSeq_Fungi_DecontamV2_HMS_genus = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_HMS, "genus")
psRep200_HiSeq_Fungi_DecontamV2_HMS_species = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_HMS, "species")
## Aggregate counts - BCM
psRep200_HiSeq_Fungi_DecontamV2_BCM_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_BCM, "phylum")
psRep200_HiSeq_Fungi_DecontamV2_BCM_class = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_BCM, "class")
psRep200_HiSeq_Fungi_DecontamV2_BCM_order = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_BCM, "order")
psRep200_HiSeq_Fungi_DecontamV2_BCM_family = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_BCM, "family")
psRep200_HiSeq_Fungi_DecontamV2_BCM_genus = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_BCM, "genus")
psRep200_HiSeq_Fungi_DecontamV2_BCM_species = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_BCM, "species")
## Aggregate counts - MDA
psRep200_HiSeq_Fungi_DecontamV2_MDA_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_MDA, "phylum")
psRep200_HiSeq_Fungi_DecontamV2_MDA_class = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_MDA, "class")
psRep200_HiSeq_Fungi_DecontamV2_MDA_order = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_MDA, "order")
psRep200_HiSeq_Fungi_DecontamV2_MDA_family = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_MDA, "family")
psRep200_HiSeq_Fungi_DecontamV2_MDA_genus = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_MDA, "genus")
psRep200_HiSeq_Fungi_DecontamV2_MDA_species = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_MDA, "species")
## Aggregate counts - WashU
psRep200_HiSeq_Fungi_DecontamV2_WashU_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_WashU, "phylum")
psRep200_HiSeq_Fungi_DecontamV2_WashU_class = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_WashU, "class")
psRep200_HiSeq_Fungi_DecontamV2_WashU_order = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_WashU, "order")
psRep200_HiSeq_Fungi_DecontamV2_WashU_family = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_WashU, "family")
psRep200_HiSeq_Fungi_DecontamV2_WashU_genus = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_WashU, "genus")
psRep200_HiSeq_Fungi_DecontamV2_WashU_species = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_WashU, "species")
## Aggregate counts - Broad_WGS
psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS, "phylum")
psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_class = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS, "class")
psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_order = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS, "order")
psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_family = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS, "family")
psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_genus = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS, "genus")
psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_species = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS, "species")
## Aggregate counts - UNC
psRep200_HiSeq_Fungi_DecontamV2_UNC_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_UNC, "phylum")
psRep200_HiSeq_Fungi_DecontamV2_UNC_class = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_UNC, "class")
psRep200_HiSeq_Fungi_DecontamV2_UNC_order = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_UNC, "order")
psRep200_HiSeq_Fungi_DecontamV2_UNC_family = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_UNC, "family")
psRep200_HiSeq_Fungi_DecontamV2_UNC_genus = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_UNC, "genus")
psRep200_HiSeq_Fungi_DecontamV2_UNC_species = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_UNC, "species")
## Aggregate counts - HMS
psRep200_HiSeq_Fungi_DecontamV2_CMS_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_CMS, "phylum")
psRep200_HiSeq_Fungi_DecontamV2_CMS_class = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_CMS, "class")
psRep200_HiSeq_Fungi_DecontamV2_CMS_order = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_CMS, "order")
psRep200_HiSeq_Fungi_DecontamV2_CMS_family = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_CMS, "family")
psRep200_HiSeq_Fungi_DecontamV2_CMS_genus = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_CMS, "genus")
psRep200_HiSeq_Fungi_DecontamV2_CMS_species = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_CMS, "species")
## Aggregate counts - Broad_RNA
psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA, "phylum")
psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_class = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA, "class")
psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_order = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA, "order")
psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_family = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA, "family")
psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_genus = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA, "genus")
psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_species = aggregate_taxa(psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA, "species")

#--------------------Extract data frames of aggregated counts--------------------#
## HMS
df_psRep200_HiSeq_Fungi_DecontamV2_HMS_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_HMS_phylum)))
df_psRep200_HiSeq_Fungi_DecontamV2_HMS_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_HMS_class)))
df_psRep200_HiSeq_Fungi_DecontamV2_HMS_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_HMS_order)))
df_psRep200_HiSeq_Fungi_DecontamV2_HMS_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_HMS_family)))
df_psRep200_HiSeq_Fungi_DecontamV2_HMS_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_HMS_genus)))
df_psRep200_HiSeq_Fungi_DecontamV2_HMS_species <- rep200_HiSeq_Fungi_DecontamV2_HMS
colnames(df_psRep200_HiSeq_Fungi_DecontamV2_HMS_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_DecontamV2_HMS_species), "species"]
## BCM
df_psRep200_HiSeq_Fungi_DecontamV2_BCM_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_BCM_phylum)))
df_psRep200_HiSeq_Fungi_DecontamV2_BCM_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_BCM_class)))
df_psRep200_HiSeq_Fungi_DecontamV2_BCM_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_BCM_order)))
df_psRep200_HiSeq_Fungi_DecontamV2_BCM_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_BCM_family)))
df_psRep200_HiSeq_Fungi_DecontamV2_BCM_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_BCM_genus)))
df_psRep200_HiSeq_Fungi_DecontamV2_BCM_species <- rep200_HiSeq_Fungi_DecontamV2_BCM
colnames(df_psRep200_HiSeq_Fungi_DecontamV2_BCM_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_DecontamV2_BCM_species), "species"]
## MDA
df_psRep200_HiSeq_Fungi_DecontamV2_MDA_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_MDA_phylum)))
df_psRep200_HiSeq_Fungi_DecontamV2_MDA_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_MDA_class)))
df_psRep200_HiSeq_Fungi_DecontamV2_MDA_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_MDA_order)))
df_psRep200_HiSeq_Fungi_DecontamV2_MDA_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_MDA_family)))
df_psRep200_HiSeq_Fungi_DecontamV2_MDA_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_MDA_genus)))
df_psRep200_HiSeq_Fungi_DecontamV2_MDA_species <- rep200_HiSeq_Fungi_DecontamV2_MDA
colnames(df_psRep200_HiSeq_Fungi_DecontamV2_MDA_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_DecontamV2_MDA_species), "species"]
## WashU
df_psRep200_HiSeq_Fungi_DecontamV2_WashU_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_WashU_phylum)))
df_psRep200_HiSeq_Fungi_DecontamV2_WashU_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_WashU_class)))
df_psRep200_HiSeq_Fungi_DecontamV2_WashU_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_WashU_order)))
df_psRep200_HiSeq_Fungi_DecontamV2_WashU_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_WashU_family)))
df_psRep200_HiSeq_Fungi_DecontamV2_WashU_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_WashU_genus)))
df_psRep200_HiSeq_Fungi_DecontamV2_WashU_species <- rep200_HiSeq_Fungi_DecontamV2_WashU
colnames(df_psRep200_HiSeq_Fungi_DecontamV2_WashU_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_DecontamV2_WashU_species), "species"]
## Broad_WGS
df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_phylum)))
df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_class)))
df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_order)))
df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_family)))
df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_genus)))
df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_species <- rep200_HiSeq_Fungi_DecontamV2_Broad_WGS
colnames(df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_species), "species"]
## UNC
df_psRep200_HiSeq_Fungi_DecontamV2_UNC_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_UNC_phylum)))
df_psRep200_HiSeq_Fungi_DecontamV2_UNC_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_UNC_class)))
df_psRep200_HiSeq_Fungi_DecontamV2_UNC_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_UNC_order)))
df_psRep200_HiSeq_Fungi_DecontamV2_UNC_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_UNC_family)))
df_psRep200_HiSeq_Fungi_DecontamV2_UNC_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_UNC_genus)))
df_psRep200_HiSeq_Fungi_DecontamV2_UNC_species <- rep200_HiSeq_Fungi_DecontamV2_UNC
colnames(df_psRep200_HiSeq_Fungi_DecontamV2_UNC_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_DecontamV2_UNC_species), "species"]
## CMS
df_psRep200_HiSeq_Fungi_DecontamV2_CMS_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_CMS_phylum)))
df_psRep200_HiSeq_Fungi_DecontamV2_CMS_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_CMS_class)))
df_psRep200_HiSeq_Fungi_DecontamV2_CMS_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_CMS_order)))
df_psRep200_HiSeq_Fungi_DecontamV2_CMS_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_CMS_family)))
df_psRep200_HiSeq_Fungi_DecontamV2_CMS_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_CMS_genus)))
df_psRep200_HiSeq_Fungi_DecontamV2_CMS_species <- rep200_HiSeq_Fungi_DecontamV2_CMS
colnames(df_psRep200_HiSeq_Fungi_DecontamV2_CMS_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_DecontamV2_CMS_species), "species"]
## Broad_RNA
df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_phylum)))
df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_class)))
df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_order)))
df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_family)))
df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_genus)))
df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_species <- rep200_HiSeq_Fungi_DecontamV2_Broad_RNA
colnames(df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_species), "species"]

#----------------------------------------------------------------------#
#--------------------Intersect with Weizmann cohort--------------------#
# DIFFERENCE FROM ABOVE IS THAT INTERSECTION IS DONE ON FULL DATA NOT DECONTAMINATED DATA
#----------------------------------------------------------------------#

## Load shared features with Weizmann
load("Interim_data/shared_fungi_features_at_each_taxa_level_29Mar22.RData")

rep200TaxSplit_Fungi_Paired_to_Weizmann <- read.csv("Supporting_data/rep200TaxSplit_Fungi_Paired_To_Weizmann_Final.csv", stringsAsFactors = FALSE, row.names = 1)

# Build phyloseq object
psRep200_HiSeq_Fungi_HMS <- phyloseq(otu_table(rep200Data_WGS_RNA_Matched_Fungi, taxa_are_rows = FALSE), 
                                                tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_HMS))
psRep200_HiSeq_Fungi_BCM <- phyloseq(otu_table(rep200Data_WGS_RNA_Matched_Fungi, taxa_are_rows = FALSE), 
                                                tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_BCM))
psRep200_HiSeq_Fungi_MDA <- phyloseq(otu_table(rep200Data_WGS_RNA_Matched_Fungi, taxa_are_rows = FALSE), 
                                                tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_MDA))
psRep200_HiSeq_Fungi_WashU <- phyloseq(otu_table(rep200Data_WGS_RNA_Matched_Fungi, taxa_are_rows = FALSE), 
                                                  tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_WashU))
psRep200_HiSeq_Fungi_Broad_WGS <- phyloseq(otu_table(rep200Data_WGS_RNA_Matched_Fungi, taxa_are_rows = FALSE), 
                                                      tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS))
psRep200_HiSeq_Fungi_UNC <- phyloseq(otu_table(rep200Data_WGS_RNA_Matched_Fungi, taxa_are_rows = FALSE), 
                                                tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_UNC))
psRep200_HiSeq_Fungi_CMS <- phyloseq(otu_table(rep200Data_WGS_RNA_Matched_Fungi, taxa_are_rows = FALSE), 
                                                tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_CMS))
psRep200_HiSeq_Fungi_Broad_RNA <- phyloseq(otu_table(rep200Data_WGS_RNA_Matched_Fungi, taxa_are_rows = FALSE), 
                                                      tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA))
## Aggregate counts - HMS
psRep200_HiSeq_Fungi_HMS_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_HMS, "phylum")
psRep200_HiSeq_Fungi_HMS_class = aggregate_taxa(psRep200_HiSeq_Fungi_HMS, "class")
psRep200_HiSeq_Fungi_HMS_order = aggregate_taxa(psRep200_HiSeq_Fungi_HMS, "order")
psRep200_HiSeq_Fungi_HMS_family = aggregate_taxa(psRep200_HiSeq_Fungi_HMS, "family")
psRep200_HiSeq_Fungi_HMS_genus = aggregate_taxa(psRep200_HiSeq_Fungi_HMS, "genus")
psRep200_HiSeq_Fungi_HMS_species = aggregate_taxa(psRep200_HiSeq_Fungi_HMS, "species")
## Aggregate counts - BCM
psRep200_HiSeq_Fungi_BCM_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_BCM, "phylum")
psRep200_HiSeq_Fungi_BCM_class = aggregate_taxa(psRep200_HiSeq_Fungi_BCM, "class")
psRep200_HiSeq_Fungi_BCM_order = aggregate_taxa(psRep200_HiSeq_Fungi_BCM, "order")
psRep200_HiSeq_Fungi_BCM_family = aggregate_taxa(psRep200_HiSeq_Fungi_BCM, "family")
psRep200_HiSeq_Fungi_BCM_genus = aggregate_taxa(psRep200_HiSeq_Fungi_BCM, "genus")
psRep200_HiSeq_Fungi_BCM_species = aggregate_taxa(psRep200_HiSeq_Fungi_BCM, "species")
## Aggregate counts - MDA
psRep200_HiSeq_Fungi_MDA_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_MDA, "phylum")
psRep200_HiSeq_Fungi_MDA_class = aggregate_taxa(psRep200_HiSeq_Fungi_MDA, "class")
psRep200_HiSeq_Fungi_MDA_order = aggregate_taxa(psRep200_HiSeq_Fungi_MDA, "order")
psRep200_HiSeq_Fungi_MDA_family = aggregate_taxa(psRep200_HiSeq_Fungi_MDA, "family")
psRep200_HiSeq_Fungi_MDA_genus = aggregate_taxa(psRep200_HiSeq_Fungi_MDA, "genus")
psRep200_HiSeq_Fungi_MDA_species = aggregate_taxa(psRep200_HiSeq_Fungi_MDA, "species")
## Aggregate counts - WashU
psRep200_HiSeq_Fungi_WashU_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_WashU, "phylum")
psRep200_HiSeq_Fungi_WashU_class = aggregate_taxa(psRep200_HiSeq_Fungi_WashU, "class")
psRep200_HiSeq_Fungi_WashU_order = aggregate_taxa(psRep200_HiSeq_Fungi_WashU, "order")
psRep200_HiSeq_Fungi_WashU_family = aggregate_taxa(psRep200_HiSeq_Fungi_WashU, "family")
psRep200_HiSeq_Fungi_WashU_genus = aggregate_taxa(psRep200_HiSeq_Fungi_WashU, "genus")
psRep200_HiSeq_Fungi_WashU_species = aggregate_taxa(psRep200_HiSeq_Fungi_WashU, "species")
## Aggregate counts - Broad_WGS
psRep200_HiSeq_Fungi_Broad_WGS_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_Broad_WGS, "phylum")
psRep200_HiSeq_Fungi_Broad_WGS_class = aggregate_taxa(psRep200_HiSeq_Fungi_Broad_WGS, "class")
psRep200_HiSeq_Fungi_Broad_WGS_order = aggregate_taxa(psRep200_HiSeq_Fungi_Broad_WGS, "order")
psRep200_HiSeq_Fungi_Broad_WGS_family = aggregate_taxa(psRep200_HiSeq_Fungi_Broad_WGS, "family")
psRep200_HiSeq_Fungi_Broad_WGS_genus = aggregate_taxa(psRep200_HiSeq_Fungi_Broad_WGS, "genus")
psRep200_HiSeq_Fungi_Broad_WGS_species = aggregate_taxa(psRep200_HiSeq_Fungi_Broad_WGS, "species")
## Aggregate counts - UNC
psRep200_HiSeq_Fungi_UNC_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_UNC, "phylum")
psRep200_HiSeq_Fungi_UNC_class = aggregate_taxa(psRep200_HiSeq_Fungi_UNC, "class")
psRep200_HiSeq_Fungi_UNC_order = aggregate_taxa(psRep200_HiSeq_Fungi_UNC, "order")
psRep200_HiSeq_Fungi_UNC_family = aggregate_taxa(psRep200_HiSeq_Fungi_UNC, "family")
psRep200_HiSeq_Fungi_UNC_genus = aggregate_taxa(psRep200_HiSeq_Fungi_UNC, "genus")
psRep200_HiSeq_Fungi_UNC_species = aggregate_taxa(psRep200_HiSeq_Fungi_UNC, "species")
## Aggregate counts - HMS
psRep200_HiSeq_Fungi_CMS_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_CMS, "phylum")
psRep200_HiSeq_Fungi_CMS_class = aggregate_taxa(psRep200_HiSeq_Fungi_CMS, "class")
psRep200_HiSeq_Fungi_CMS_order = aggregate_taxa(psRep200_HiSeq_Fungi_CMS, "order")
psRep200_HiSeq_Fungi_CMS_family = aggregate_taxa(psRep200_HiSeq_Fungi_CMS, "family")
psRep200_HiSeq_Fungi_CMS_genus = aggregate_taxa(psRep200_HiSeq_Fungi_CMS, "genus")
psRep200_HiSeq_Fungi_CMS_species = aggregate_taxa(psRep200_HiSeq_Fungi_CMS, "species")
## Aggregate counts - Broad_RNA
psRep200_HiSeq_Fungi_Broad_RNA_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_Broad_RNA, "phylum")
psRep200_HiSeq_Fungi_Broad_RNA_class = aggregate_taxa(psRep200_HiSeq_Fungi_Broad_RNA, "class")
psRep200_HiSeq_Fungi_Broad_RNA_order = aggregate_taxa(psRep200_HiSeq_Fungi_Broad_RNA, "order")
psRep200_HiSeq_Fungi_Broad_RNA_family = aggregate_taxa(psRep200_HiSeq_Fungi_Broad_RNA, "family")
psRep200_HiSeq_Fungi_Broad_RNA_genus = aggregate_taxa(psRep200_HiSeq_Fungi_Broad_RNA, "genus")
psRep200_HiSeq_Fungi_Broad_RNA_species = aggregate_taxa(psRep200_HiSeq_Fungi_Broad_RNA, "species")

#--------------------Extract data frames of aggregated counts--------------------#
## HMS
df_psRep200_HiSeq_Fungi_HMS_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_HMS_phylum)))
df_psRep200_HiSeq_Fungi_HMS_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_HMS_class)))
df_psRep200_HiSeq_Fungi_HMS_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_HMS_order)))
df_psRep200_HiSeq_Fungi_HMS_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_HMS_family)))
df_psRep200_HiSeq_Fungi_HMS_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_HMS_genus)))
df_psRep200_HiSeq_Fungi_HMS_species <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_HMS),]
colnames(df_psRep200_HiSeq_Fungi_HMS_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_HMS_species), "species"]
## BCM
df_psRep200_HiSeq_Fungi_BCM_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_BCM_phylum)))
df_psRep200_HiSeq_Fungi_BCM_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_BCM_class)))
df_psRep200_HiSeq_Fungi_BCM_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_BCM_order)))
df_psRep200_HiSeq_Fungi_BCM_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_BCM_family)))
df_psRep200_HiSeq_Fungi_BCM_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_BCM_genus)))
df_psRep200_HiSeq_Fungi_BCM_species <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_BCM),]
colnames(df_psRep200_HiSeq_Fungi_BCM_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_BCM_species), "species"]
## MDA
df_psRep200_HiSeq_Fungi_MDA_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_MDA_phylum)))
df_psRep200_HiSeq_Fungi_MDA_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_MDA_class)))
df_psRep200_HiSeq_Fungi_MDA_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_MDA_order)))
df_psRep200_HiSeq_Fungi_MDA_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_MDA_family)))
df_psRep200_HiSeq_Fungi_MDA_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_MDA_genus)))
df_psRep200_HiSeq_Fungi_MDA_species <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_MDA),]
colnames(df_psRep200_HiSeq_Fungi_MDA_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_MDA_species), "species"]
## WashU
df_psRep200_HiSeq_Fungi_WashU_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_WashU_phylum)))
df_psRep200_HiSeq_Fungi_WashU_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_WashU_class)))
df_psRep200_HiSeq_Fungi_WashU_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_WashU_order)))
df_psRep200_HiSeq_Fungi_WashU_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_WashU_family)))
df_psRep200_HiSeq_Fungi_WashU_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_WashU_genus)))
df_psRep200_HiSeq_Fungi_WashU_species <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_WashU),]
colnames(df_psRep200_HiSeq_Fungi_WashU_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_WashU_species), "species"]
## Broad_WGS
df_psRep200_HiSeq_Fungi_Broad_WGS_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Broad_WGS_phylum)))
df_psRep200_HiSeq_Fungi_Broad_WGS_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Broad_WGS_class)))
df_psRep200_HiSeq_Fungi_Broad_WGS_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Broad_WGS_order)))
df_psRep200_HiSeq_Fungi_Broad_WGS_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Broad_WGS_family)))
df_psRep200_HiSeq_Fungi_Broad_WGS_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Broad_WGS_genus)))
df_psRep200_HiSeq_Fungi_Broad_WGS_species <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS),]
colnames(df_psRep200_HiSeq_Fungi_Broad_WGS_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_Broad_WGS_species), "species"]
## UNC
df_psRep200_HiSeq_Fungi_UNC_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_UNC_phylum)))
df_psRep200_HiSeq_Fungi_UNC_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_UNC_class)))
df_psRep200_HiSeq_Fungi_UNC_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_UNC_order)))
df_psRep200_HiSeq_Fungi_UNC_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_UNC_family)))
df_psRep200_HiSeq_Fungi_UNC_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_UNC_genus)))
df_psRep200_HiSeq_Fungi_UNC_species <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_UNC),]
colnames(df_psRep200_HiSeq_Fungi_UNC_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_UNC_species), "species"]
## CMS
df_psRep200_HiSeq_Fungi_CMS_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_CMS_phylum)))
df_psRep200_HiSeq_Fungi_CMS_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_CMS_class)))
df_psRep200_HiSeq_Fungi_CMS_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_CMS_order)))
df_psRep200_HiSeq_Fungi_CMS_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_CMS_family)))
df_psRep200_HiSeq_Fungi_CMS_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_CMS_genus)))
df_psRep200_HiSeq_Fungi_CMS_species <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_CMS),]
colnames(df_psRep200_HiSeq_Fungi_CMS_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_CMS_species), "species"]
## Broad_RNA
df_psRep200_HiSeq_Fungi_Broad_RNA_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Broad_RNA_phylum)))
df_psRep200_HiSeq_Fungi_Broad_RNA_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Broad_RNA_class)))
df_psRep200_HiSeq_Fungi_Broad_RNA_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Broad_RNA_order)))
df_psRep200_HiSeq_Fungi_Broad_RNA_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Broad_RNA_family)))
df_psRep200_HiSeq_Fungi_Broad_RNA_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Broad_RNA_genus)))
df_psRep200_HiSeq_Fungi_Broad_RNA_species <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA),]
colnames(df_psRep200_HiSeq_Fungi_Broad_RNA_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_Broad_RNA_species), "species"]

#---------------------------------------#
## HMS -- subset features
df_psRep200_HiSeq_Fungi_HMS_phylum_Shared <- df_psRep200_HiSeq_Fungi_HMS_phylum[,colnames(df_psRep200_HiSeq_Fungi_HMS_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_HMS_class_Shared <- df_psRep200_HiSeq_Fungi_HMS_class[,colnames(df_psRep200_HiSeq_Fungi_HMS_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_HMS_order_Shared <- df_psRep200_HiSeq_Fungi_HMS_order[,colnames(df_psRep200_HiSeq_Fungi_HMS_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_HMS_family_Shared <- df_psRep200_HiSeq_Fungi_HMS_family[,colnames(df_psRep200_HiSeq_Fungi_HMS_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_HMS_genus_Shared <- df_psRep200_HiSeq_Fungi_HMS_genus[,colnames(df_psRep200_HiSeq_Fungi_HMS_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_HMS_species_Shared <- df_psRep200_HiSeq_Fungi_HMS_species[,colnames(df_psRep200_HiSeq_Fungi_HMS_species) %in% sharedSpecies]
## HMS -- remove zero sum samples
df_psRep200_HiSeq_Fungi_HMS_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_HMS_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_HMS_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_HMS_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_HMS_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_HMS_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_HMS_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_HMS_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_HMS_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_HMS_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_HMS_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_HMS_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_HMS_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_HMS_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_HMS_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_HMS_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_HMS_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_HMS_species_Shared)==0,]
## HMS - metadata
metaQiitaCombined_Nonzero_DecontamV2_HMS_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_HMS[rownames(df_psRep200_HiSeq_Fungi_HMS_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_HMS_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_HMS[rownames(df_psRep200_HiSeq_Fungi_HMS_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_HMS_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_HMS[rownames(df_psRep200_HiSeq_Fungi_HMS_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_HMS_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_HMS[rownames(df_psRep200_HiSeq_Fungi_HMS_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_HMS_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_HMS[rownames(df_psRep200_HiSeq_Fungi_HMS_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_HMS_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_HMS[rownames(df_psRep200_HiSeq_Fungi_HMS_species_Shared_Nonzero),])

## BCM -- subset features
df_psRep200_HiSeq_Fungi_BCM_phylum_Shared <- df_psRep200_HiSeq_Fungi_BCM_phylum[,colnames(df_psRep200_HiSeq_Fungi_BCM_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_BCM_class_Shared <- df_psRep200_HiSeq_Fungi_BCM_class[,colnames(df_psRep200_HiSeq_Fungi_BCM_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_BCM_order_Shared <- df_psRep200_HiSeq_Fungi_BCM_order[,colnames(df_psRep200_HiSeq_Fungi_BCM_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_BCM_family_Shared <- df_psRep200_HiSeq_Fungi_BCM_family[,colnames(df_psRep200_HiSeq_Fungi_BCM_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_BCM_genus_Shared <- df_psRep200_HiSeq_Fungi_BCM_genus[,colnames(df_psRep200_HiSeq_Fungi_BCM_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_BCM_species_Shared <- df_psRep200_HiSeq_Fungi_BCM_species[,colnames(df_psRep200_HiSeq_Fungi_BCM_species) %in% sharedSpecies]
## BCM -- remove zero sum samples
df_psRep200_HiSeq_Fungi_BCM_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_BCM_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_BCM_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_BCM_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_BCM_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_BCM_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_BCM_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_BCM_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_BCM_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_BCM_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_BCM_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_BCM_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_BCM_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_BCM_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_BCM_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_BCM_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_BCM_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_BCM_species_Shared)==0,]
## BCM - metadata
metaQiitaCombined_Nonzero_DecontamV2_BCM_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_BCM[rownames(df_psRep200_HiSeq_Fungi_BCM_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_BCM_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_BCM[rownames(df_psRep200_HiSeq_Fungi_BCM_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_BCM_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_BCM[rownames(df_psRep200_HiSeq_Fungi_BCM_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_BCM_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_BCM[rownames(df_psRep200_HiSeq_Fungi_BCM_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_BCM_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_BCM[rownames(df_psRep200_HiSeq_Fungi_BCM_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_BCM_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_BCM[rownames(df_psRep200_HiSeq_Fungi_BCM_species_Shared_Nonzero),])

## MDA -- subset features
df_psRep200_HiSeq_Fungi_MDA_phylum_Shared <- df_psRep200_HiSeq_Fungi_MDA_phylum[,colnames(df_psRep200_HiSeq_Fungi_MDA_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_MDA_class_Shared <- df_psRep200_HiSeq_Fungi_MDA_class[,colnames(df_psRep200_HiSeq_Fungi_MDA_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_MDA_order_Shared <- df_psRep200_HiSeq_Fungi_MDA_order[,colnames(df_psRep200_HiSeq_Fungi_MDA_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_MDA_family_Shared <- df_psRep200_HiSeq_Fungi_MDA_family[,colnames(df_psRep200_HiSeq_Fungi_MDA_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_MDA_genus_Shared <- df_psRep200_HiSeq_Fungi_MDA_genus[,colnames(df_psRep200_HiSeq_Fungi_MDA_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_MDA_species_Shared <- df_psRep200_HiSeq_Fungi_MDA_species[,colnames(df_psRep200_HiSeq_Fungi_MDA_species) %in% sharedSpecies]
## MDA -- remove zero sum samples
df_psRep200_HiSeq_Fungi_MDA_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_MDA_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_MDA_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_MDA_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_MDA_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_MDA_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_MDA_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_MDA_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_MDA_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_MDA_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_MDA_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_MDA_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_MDA_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_MDA_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_MDA_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_MDA_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_MDA_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_MDA_species_Shared)==0,]
## MDA - metadata
metaQiitaCombined_Nonzero_DecontamV2_MDA_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_MDA[rownames(df_psRep200_HiSeq_Fungi_MDA_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_MDA_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_MDA[rownames(df_psRep200_HiSeq_Fungi_MDA_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_MDA_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_MDA[rownames(df_psRep200_HiSeq_Fungi_MDA_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_MDA_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_MDA[rownames(df_psRep200_HiSeq_Fungi_MDA_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_MDA_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_MDA[rownames(df_psRep200_HiSeq_Fungi_MDA_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_MDA_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_MDA[rownames(df_psRep200_HiSeq_Fungi_MDA_species_Shared_Nonzero),])

## WashU -- subset features
df_psRep200_HiSeq_Fungi_WashU_phylum_Shared <- df_psRep200_HiSeq_Fungi_WashU_phylum[,colnames(df_psRep200_HiSeq_Fungi_WashU_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_WashU_class_Shared <- df_psRep200_HiSeq_Fungi_WashU_class[,colnames(df_psRep200_HiSeq_Fungi_WashU_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_WashU_order_Shared <- df_psRep200_HiSeq_Fungi_WashU_order[,colnames(df_psRep200_HiSeq_Fungi_WashU_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_WashU_family_Shared <- df_psRep200_HiSeq_Fungi_WashU_family[,colnames(df_psRep200_HiSeq_Fungi_WashU_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_WashU_genus_Shared <- df_psRep200_HiSeq_Fungi_WashU_genus[,colnames(df_psRep200_HiSeq_Fungi_WashU_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_WashU_species_Shared <- df_psRep200_HiSeq_Fungi_WashU_species[,colnames(df_psRep200_HiSeq_Fungi_WashU_species) %in% sharedSpecies]
## WashU -- remove zero sum samples
df_psRep200_HiSeq_Fungi_WashU_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_WashU_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_WashU_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_WashU_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_WashU_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_WashU_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_WashU_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_WashU_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_WashU_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_WashU_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_WashU_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_WashU_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_WashU_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_WashU_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_WashU_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_WashU_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_WashU_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_WashU_species_Shared)==0,]
## WashU - metadata
metaQiitaCombined_Nonzero_DecontamV2_WashU_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_WashU[rownames(df_psRep200_HiSeq_Fungi_WashU_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_WashU_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_WashU[rownames(df_psRep200_HiSeq_Fungi_WashU_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_WashU_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_WashU[rownames(df_psRep200_HiSeq_Fungi_WashU_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_WashU_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_WashU[rownames(df_psRep200_HiSeq_Fungi_WashU_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_WashU_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_WashU[rownames(df_psRep200_HiSeq_Fungi_WashU_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_WashU_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_WashU[rownames(df_psRep200_HiSeq_Fungi_WashU_species_Shared_Nonzero),])

## Broad_WGS -- subset features
df_psRep200_HiSeq_Fungi_Broad_WGS_phylum_Shared <- df_psRep200_HiSeq_Fungi_Broad_WGS_phylum[,colnames(df_psRep200_HiSeq_Fungi_Broad_WGS_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_Broad_WGS_class_Shared <- df_psRep200_HiSeq_Fungi_Broad_WGS_class[,colnames(df_psRep200_HiSeq_Fungi_Broad_WGS_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_Broad_WGS_order_Shared <- df_psRep200_HiSeq_Fungi_Broad_WGS_order[,colnames(df_psRep200_HiSeq_Fungi_Broad_WGS_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_Broad_WGS_family_Shared <- df_psRep200_HiSeq_Fungi_Broad_WGS_family[,colnames(df_psRep200_HiSeq_Fungi_Broad_WGS_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_Broad_WGS_genus_Shared <- df_psRep200_HiSeq_Fungi_Broad_WGS_genus[,colnames(df_psRep200_HiSeq_Fungi_Broad_WGS_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_Broad_WGS_species_Shared <- df_psRep200_HiSeq_Fungi_Broad_WGS_species[,colnames(df_psRep200_HiSeq_Fungi_Broad_WGS_species) %in% sharedSpecies]
## Broad_WGS -- remove zero sum samples
df_psRep200_HiSeq_Fungi_Broad_WGS_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_WGS_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Broad_WGS_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_Broad_WGS_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_WGS_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Broad_WGS_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_Broad_WGS_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_WGS_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Broad_WGS_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_Broad_WGS_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_WGS_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Broad_WGS_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_Broad_WGS_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_WGS_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Broad_WGS_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_Broad_WGS_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_WGS_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Broad_WGS_species_Shared)==0,]
## Broad_WGS - metadata
metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Broad_WGS_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Broad_WGS_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Broad_WGS_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Broad_WGS_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Broad_WGS_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Broad_WGS_species_Shared_Nonzero),])

## UNC -- subset features
df_psRep200_HiSeq_Fungi_UNC_phylum_Shared <- df_psRep200_HiSeq_Fungi_UNC_phylum[,colnames(df_psRep200_HiSeq_Fungi_UNC_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_UNC_class_Shared <- df_psRep200_HiSeq_Fungi_UNC_class[,colnames(df_psRep200_HiSeq_Fungi_UNC_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_UNC_order_Shared <- df_psRep200_HiSeq_Fungi_UNC_order[,colnames(df_psRep200_HiSeq_Fungi_UNC_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_UNC_family_Shared <- df_psRep200_HiSeq_Fungi_UNC_family[,colnames(df_psRep200_HiSeq_Fungi_UNC_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_UNC_genus_Shared <- df_psRep200_HiSeq_Fungi_UNC_genus[,colnames(df_psRep200_HiSeq_Fungi_UNC_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_UNC_species_Shared <- df_psRep200_HiSeq_Fungi_UNC_species[,colnames(df_psRep200_HiSeq_Fungi_UNC_species) %in% sharedSpecies]
## UNC -- remove zero sum samples
df_psRep200_HiSeq_Fungi_UNC_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_UNC_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_UNC_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_UNC_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_UNC_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_UNC_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_UNC_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_UNC_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_UNC_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_UNC_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_UNC_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_UNC_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_UNC_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_UNC_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_UNC_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_UNC_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_UNC_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_UNC_species_Shared)==0,]
## UNC - metadata
metaQiitaCombined_Nonzero_DecontamV2_UNC_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_UNC[rownames(df_psRep200_HiSeq_Fungi_UNC_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_UNC_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_UNC[rownames(df_psRep200_HiSeq_Fungi_UNC_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_UNC_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_UNC[rownames(df_psRep200_HiSeq_Fungi_UNC_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_UNC_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_UNC[rownames(df_psRep200_HiSeq_Fungi_UNC_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_UNC_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_UNC[rownames(df_psRep200_HiSeq_Fungi_UNC_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_UNC_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_UNC[rownames(df_psRep200_HiSeq_Fungi_UNC_species_Shared_Nonzero),])

## CMS -- subset features
df_psRep200_HiSeq_Fungi_CMS_phylum_Shared <- df_psRep200_HiSeq_Fungi_CMS_phylum[,colnames(df_psRep200_HiSeq_Fungi_CMS_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_CMS_class_Shared <- df_psRep200_HiSeq_Fungi_CMS_class[,colnames(df_psRep200_HiSeq_Fungi_CMS_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_CMS_order_Shared <- df_psRep200_HiSeq_Fungi_CMS_order[,colnames(df_psRep200_HiSeq_Fungi_CMS_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_CMS_family_Shared <- df_psRep200_HiSeq_Fungi_CMS_family[,colnames(df_psRep200_HiSeq_Fungi_CMS_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_CMS_genus_Shared <- df_psRep200_HiSeq_Fungi_CMS_genus[,colnames(df_psRep200_HiSeq_Fungi_CMS_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_CMS_species_Shared <- df_psRep200_HiSeq_Fungi_CMS_species[,colnames(df_psRep200_HiSeq_Fungi_CMS_species) %in% sharedSpecies]
## CMS -- remove zero sum samples
df_psRep200_HiSeq_Fungi_CMS_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_CMS_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_CMS_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_CMS_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_CMS_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_CMS_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_CMS_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_CMS_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_CMS_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_CMS_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_CMS_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_CMS_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_CMS_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_CMS_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_CMS_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_CMS_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_CMS_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_CMS_species_Shared)==0,]
## CMS - metadata
metaQiitaCombined_Nonzero_DecontamV2_CMS_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_CMS[rownames(df_psRep200_HiSeq_Fungi_CMS_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_CMS_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_CMS[rownames(df_psRep200_HiSeq_Fungi_CMS_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_CMS_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_CMS[rownames(df_psRep200_HiSeq_Fungi_CMS_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_CMS_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_CMS[rownames(df_psRep200_HiSeq_Fungi_CMS_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_CMS_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_CMS[rownames(df_psRep200_HiSeq_Fungi_CMS_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_CMS_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_CMS[rownames(df_psRep200_HiSeq_Fungi_CMS_species_Shared_Nonzero),])

## Broad_RNA -- subset features
df_psRep200_HiSeq_Fungi_Broad_RNA_phylum_Shared <- df_psRep200_HiSeq_Fungi_Broad_RNA_phylum[,colnames(df_psRep200_HiSeq_Fungi_Broad_RNA_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_Broad_RNA_class_Shared <- df_psRep200_HiSeq_Fungi_Broad_RNA_class[,colnames(df_psRep200_HiSeq_Fungi_Broad_RNA_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_Broad_RNA_order_Shared <- df_psRep200_HiSeq_Fungi_Broad_RNA_order[,colnames(df_psRep200_HiSeq_Fungi_Broad_RNA_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_Broad_RNA_family_Shared <- df_psRep200_HiSeq_Fungi_Broad_RNA_family[,colnames(df_psRep200_HiSeq_Fungi_Broad_RNA_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_Broad_RNA_genus_Shared <- df_psRep200_HiSeq_Fungi_Broad_RNA_genus[,colnames(df_psRep200_HiSeq_Fungi_Broad_RNA_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_Broad_RNA_species_Shared <- df_psRep200_HiSeq_Fungi_Broad_RNA_species[,colnames(df_psRep200_HiSeq_Fungi_Broad_RNA_species) %in% sharedSpecies]
## Broad_RNA -- remove zero sum samples
df_psRep200_HiSeq_Fungi_Broad_RNA_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_RNA_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Broad_RNA_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_Broad_RNA_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_RNA_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Broad_RNA_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_Broad_RNA_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_RNA_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Broad_RNA_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_Broad_RNA_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_RNA_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Broad_RNA_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_Broad_RNA_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_RNA_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Broad_RNA_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_Broad_RNA_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_RNA_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Broad_RNA_species_Shared)==0,]
## Broad_RNA - metadata
metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Broad_RNA_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Broad_RNA_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Broad_RNA_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Broad_RNA_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Broad_RNA_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Broad_RNA_species_Shared_Nonzero),])

#--------------------Save data for ML--------------------#

save(df_psRep200_HiSeq_Fungi_DecontamV2_HMS_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_HMS_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_HMS_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_HMS_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_HMS_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_HMS_species,
     df_psRep200_HiSeq_Fungi_HMS_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_HMS_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_HMS_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_HMS_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_HMS_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_HMS_species_Shared_Nonzero,
     # BCM
     df_psRep200_HiSeq_Fungi_DecontamV2_BCM_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_BCM_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_BCM_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_BCM_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_BCM_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_BCM_species,
     df_psRep200_HiSeq_Fungi_BCM_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_BCM_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_BCM_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_BCM_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_BCM_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_BCM_species_Shared_Nonzero,
     # MDA
     df_psRep200_HiSeq_Fungi_DecontamV2_MDA_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_MDA_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_MDA_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_MDA_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_MDA_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_MDA_species,
     df_psRep200_HiSeq_Fungi_MDA_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_MDA_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_MDA_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_MDA_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_MDA_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_MDA_species_Shared_Nonzero,
     # WashU
     df_psRep200_HiSeq_Fungi_DecontamV2_WashU_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_WashU_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_WashU_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_WashU_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_WashU_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_WashU_species,
     df_psRep200_HiSeq_Fungi_WashU_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_WashU_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_WashU_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_WashU_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_WashU_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_WashU_species_Shared_Nonzero,
     # Broad_WGS
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_species,
     df_psRep200_HiSeq_Fungi_Broad_WGS_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Broad_WGS_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Broad_WGS_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Broad_WGS_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Broad_WGS_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Broad_WGS_species_Shared_Nonzero,
     # UNC
     df_psRep200_HiSeq_Fungi_DecontamV2_UNC_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_UNC_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_UNC_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_UNC_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_UNC_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_UNC_species,
     df_psRep200_HiSeq_Fungi_UNC_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_UNC_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_UNC_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_UNC_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_UNC_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_UNC_species_Shared_Nonzero,
     # CMS
     df_psRep200_HiSeq_Fungi_DecontamV2_CMS_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_CMS_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_CMS_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_CMS_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_CMS_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_CMS_species,
     df_psRep200_HiSeq_Fungi_CMS_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_CMS_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_CMS_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_CMS_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_CMS_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_CMS_species_Shared_Nonzero,
     # Broad_RNA
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_phylum,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_class,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_order,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_family,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_genus,
     df_psRep200_HiSeq_Fungi_DecontamV2_Broad_RNA_species,
     df_psRep200_HiSeq_Fungi_Broad_RNA_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Broad_RNA_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Broad_RNA_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Broad_RNA_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Broad_RNA_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Broad_RNA_species_Shared_Nonzero,
     # Metadata - HMS
     metaQiitaCombined_Nonzero_DecontamV2_HMS,
     metaQiitaCombined_Nonzero_DecontamV2_HMS_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_HMS_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_HMS_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_HMS_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_HMS_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_HMS_species_Shared_Nonzero,
     # Metadata - BCM
     metaQiitaCombined_Nonzero_DecontamV2_BCM,
     metaQiitaCombined_Nonzero_DecontamV2_BCM_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_BCM_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_BCM_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_BCM_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_BCM_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_BCM_species_Shared_Nonzero,
     # Metadata - MDA
     metaQiitaCombined_Nonzero_DecontamV2_MDA,
     metaQiitaCombined_Nonzero_DecontamV2_MDA_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_MDA_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_MDA_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_MDA_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_MDA_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_MDA_species_Shared_Nonzero,
     # Metadata - WashU
     metaQiitaCombined_Nonzero_DecontamV2_WashU,
     metaQiitaCombined_Nonzero_DecontamV2_WashU_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_WashU_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_WashU_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_WashU_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_WashU_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_WashU_species_Shared_Nonzero,
     # Metadata - Broad_WGS
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_species_Shared_Nonzero,
     # Metadata - UNC
     metaQiitaCombined_Nonzero_DecontamV2_UNC,
     metaQiitaCombined_Nonzero_DecontamV2_UNC_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_UNC_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_UNC_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_UNC_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_UNC_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_UNC_species_Shared_Nonzero,
     # Metadata - CMS
     metaQiitaCombined_Nonzero_DecontamV2_CMS,
     metaQiitaCombined_Nonzero_DecontamV2_CMS_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_CMS_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_CMS_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_CMS_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_CMS_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_CMS_species_Shared_Nonzero,
     # Metadata - Broad_RNA
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_species_Shared_Nonzero,
     file = "Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy_taxa_levels_and_wz_intersect_decontamV2_2Apr22.RData")

# Scripts: S04R

#--------------------------------------------------------------------------------------------------------------------#
# Separate WGS and RNA-Seq data for ML (to compare performance between them after batch correcting)
# Using CT Voom-SNM data
#--------------------------------------------------------------------------------------------------------------------#

metaQiitaCombined_Nonzero_DecontamV2_WGS <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(experimental_strategy == "WGS") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_RNA <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()

snmDataOGUFungiDecontamV2_WGS <- snmDataOGUFungiDecontamV2[rownames(metaQiitaCombined_Nonzero_DecontamV2_WGS),]
snmDataOGUFungiDecontamV2_RNA <- snmDataOGUFungiDecontamV2[rownames(metaQiitaCombined_Nonzero_DecontamV2_RNA),]

save(metaQiitaCombined_Nonzero_DecontamV2,
     metaQiitaCombined_Nonzero_DecontamV2_WGS,
     metaQiitaCombined_Nonzero_DecontamV2_RNA,
     snmDataOGUFungiDecontamV2,
     snmDataOGUFungiDecontamV2_WGS,
     snmDataOGUFungiDecontamV2_RNA,
     file = "Interim_data/data_for_ml_tcga_wgs_vs_rna_decontamV2_2Apr22.RData")

# Scripts: S05R

#----------------------------------------------------------#
# Create phyloseq objects and aggregate counts at various taxa levels
#----------------------------------------------------------#

rep200TaxSplit_Fungi_Paired_to_Weizmann <- read.csv("Supporting_data/rep200TaxSplit_Fungi_Paired_To_Weizmann_Final.csv", stringsAsFactors = FALSE, row.names = 1)

psFungiHiSeqFungi_DecontamV2 <- phyloseq(otu_table(rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero, taxa_are_rows = FALSE), 
                                tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_DecontamV2))

## Aggregate counts
psFungiHiSeqFungi_DecontamV2_phylum = aggregate_taxa(psFungiHiSeqFungi_DecontamV2, "phylum")
psFungiHiSeqFungi_DecontamV2_class = aggregate_taxa(psFungiHiSeqFungi_DecontamV2, "class")
psFungiHiSeqFungi_DecontamV2_order = aggregate_taxa(psFungiHiSeqFungi_DecontamV2, "order")
psFungiHiSeqFungi_DecontamV2_family = aggregate_taxa(psFungiHiSeqFungi_DecontamV2, "family")
psFungiHiSeqFungi_DecontamV2_genus = aggregate_taxa(psFungiHiSeqFungi_DecontamV2, "genus")
psFungiHiSeqFungi_DecontamV2_species = aggregate_taxa(psFungiHiSeqFungi_DecontamV2, "species")

#----------------------------------------------------------#
# Extract aggregated and subset feature matrices
#----------------------------------------------------------#
rep200FungiDecontamV2Phylum <- data.frame(t(otu_table(psFungiHiSeqFungi_DecontamV2_phylum)))
rep200FungiDecontamV2Class <- data.frame(t(otu_table(psFungiHiSeqFungi_DecontamV2_class)))
rep200FungiDecontamV2Order <- data.frame(t(otu_table(psFungiHiSeqFungi_DecontamV2_order)))
rep200FungiDecontamV2Family <- data.frame(t(otu_table(psFungiHiSeqFungi_DecontamV2_family)))
rep200FungiDecontamV2Genus <- data.frame(t(otu_table(psFungiHiSeqFungi_DecontamV2_genus)))
rep200FungiDecontamV2Species <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero
colnames(rep200FungiDecontamV2Species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(rep200FungiDecontamV2Species), "species"]
rep200FungiDecontamV2Species[1:3,1:3]

save(metaQiitaCombined_Nonzero_DecontamV2,
     rep200FungiDecontamV2Phylum,
     rep200FungiDecontamV2Class,
     rep200FungiDecontamV2Order,
     rep200FungiDecontamV2Family,
     rep200FungiDecontamV2Genus,
     rep200FungiDecontamV2Species,
     file = "Interim_data/raw_data_tcga_all_cancers_decontamV2_2Apr22.RData")

#-----------------------------------------#
# Extract 8 cancer types to directly compare with Weizmann cancers
#-----------------------------------------#
# Cancer types to include: breast, lung, melanoma, colon, GBM, pancreas, ovary, bone
# Combine: LUAD and LUSC into LC, COAD and READ into CRC

metaQiitaCombined_Nonzero_DecontamV2_8cancer <- metaQiitaCombined_Nonzero_DecontamV2 %>%
  filter(investigation %in% c("TCGA-BRCA","TCGA-LUAD","TCGA-LUSC","TCGA-SKCM","TCGA-COAD",
                              "TCGA-READ","TCGA-GBM","TCGA-PAAD","TCGA-OV","TCGA-SARC")) %>% droplevels()
dim(metaQiitaCombined_Nonzero_DecontamV2_8cancer) # 5875   41
metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type <- as.character(metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type) 
metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type[metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type == "Lung Squamous Cell Carcinoma"] <- "Lung Cancer"
metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type[metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type == "Lung Adenocarcinoma"] <- "Lung Cancer"
metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type[metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type == "Colon Adenocarcinoma"] <- "Colorectal Cancer"
metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type[metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type == "Rectum Adenocarcinoma"] <- "Colorectal Cancer"
metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type[metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type == "Breast Invasive Carcinoma"] <- "Breast Cancer"
metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type[metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type == "Glioblastoma Multiforme"] <- "Glioblastoma"
metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type[metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type == "Pancreatic Adenocarcinoma"] <- "Pancreatic Cancer"
metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type[metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type == "Ovarian Serous Cystadenocarcinoma"] <- "Ovarian Cancer"
metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type[metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type == "Skin Cutaneous Melanoma"] <- "Melanoma"
metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type[metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type == "Sarcoma"] <- "Bone Cancer"
table(metaQiitaCombined_Nonzero_DecontamV2_8cancer$disease_type)

rep200FungiDecontamV2Phylum_8cancer <- rep200FungiDecontamV2Phylum[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer),]
rep200FungiDecontamV2Class_8cancer <- rep200FungiDecontamV2Class[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer),]
rep200FungiDecontamV2Order_8cancer <- rep200FungiDecontamV2Order[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer),]
rep200FungiDecontamV2Family_8cancer <- rep200FungiDecontamV2Family[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer),]
rep200FungiDecontamV2Genus_8cancer <- rep200FungiDecontamV2Genus[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer),]
rep200FungiDecontamV2Species_8cancer <- rep200FungiDecontamV2Species[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer),]

save(metaQiitaCombined_Nonzero_DecontamV2_8cancer,
     rep200FungiDecontamV2Phylum_8cancer,
     rep200FungiDecontamV2Class_8cancer,
     rep200FungiDecontamV2Order_8cancer,
     rep200FungiDecontamV2Family_8cancer,
     rep200FungiDecontamV2Genus_8cancer,
     rep200FungiDecontamV2Species_8cancer,
     file = "Interim_data/raw_data_tcga_8_cancers_decontamV2_2Apr22.RData")

#----------------------------------------------------------#
# Run VSNM for each taxa level
#----------------------------------------------------------#

source("00-Functions.R") # for vsnmFunctionTCGA() function
# NOTE: The Voom-SNM batch effect correction can utilize multiple
# "biological" and "technical" variables. It operates to iteratively
# remove variance due to "technical" variables while maintaining variance
# attributed to "biological" factors. It does NOT artificially increase variance
# attributed to "biological" factors. Rather, by removing the noise from "technical"
# factors in a supervised manner, the goal is to allow the underlying signal to be observed.
# Also note that the adding more variables in the Voom-SNM normalization can sometimes preclude
# convergence of the algorithm, particularly when fewer features are available. This lack of convergence
# was observed for some higher taxa levels when adding "disease_type" as a "biological" variable. Thus,
# all possible options are shown below. Including "disease_type" as a "biological" variable in addition to
# "sample_type" improves the overall modelling and downstream machine learning, but it is not necessary to
# show that the mycobiome varies between cancer types and sample types.

## TCGA full data - biological variable: sample_type | technical variables: data_submitting_center_label + experimental_strategy
# NOTE: Phylum and Class levels were attempted but did not converge, so they are not shown here
rep200FungiDecontamV2OrderVSNM_Obj <- vsnmFunctionTCGA(rep200FungiDecontamV2Order) # converges
rep200FungiDecontamV2OrderVSNM <- rep200FungiDecontamV2OrderVSNM_Obj$snmData
rep200FungiDecontamV2FamilyVSNM_Obj <- vsnmFunctionTCGA(rep200FungiDecontamV2Family) # converges
rep200FungiDecontamV2FamilyVSNM <- rep200FungiDecontamV2FamilyVSNM_Obj$snmData
rep200FungiDecontamV2GenusVSNM_Obj <- vsnmFunctionTCGA(rep200FungiDecontamV2Genus) # converges
rep200FungiDecontamV2GenusVSNM <- rep200FungiDecontamV2GenusVSNM_Obj$snmData
rep200FungiDecontamV2SpeciesVSNM_Obj <- vsnmFunctionTCGA(rep200FungiDecontamV2Species) # converges
rep200FungiDecontamV2SpeciesVSNM <- rep200FungiDecontamV2SpeciesVSNM_Obj$snmData

## Summarized 8 cancer types - biological variable: sample_type | technical variables: data_submitting_center_label + experimental_strategy
# NOTE: Phylum and Class levels were attempted but did not converge, so they are not shown here
rep200FungiDecontamV2Order_8cancer_VSNM_Obj <- vsnmFunctionTCGA(rep200FungiDecontamV2Order_8cancer, qcMetadata = metaQiitaCombined_Nonzero_DecontamV2_8cancer) # converges
rep200FungiDecontamV2Order_8cancer_VSNM <- rep200FungiDecontamV2Order_8cancer_VSNM_Obj$snmData
rep200FungiDecontamV2Family_8cancer_VSNM_Obj <- vsnmFunctionTCGA(rep200FungiDecontamV2Family_8cancer, qcMetadata = metaQiitaCombined_Nonzero_DecontamV2_8cancer) # converges
rep200FungiDecontamV2Family_8cancer_VSNM <- rep200FungiDecontamV2Family_8cancer_VSNM_Obj$snmData
rep200FungiDecontamV2Genus_8cancer_VSNM_Obj <- vsnmFunctionTCGA(rep200FungiDecontamV2Genus_8cancer, qcMetadata = metaQiitaCombined_Nonzero_DecontamV2_8cancer) # converges
rep200FungiDecontamV2Genus_8cancer_VSNM <- rep200FungiDecontamV2Genus_8cancer_VSNM_Obj$snmData
rep200FungiDecontamV2Species_8cancer_VSNM_Obj <- vsnmFunctionTCGA(rep200FungiDecontamV2Species_8cancer, qcMetadata = metaQiitaCombined_Nonzero_DecontamV2_8cancer) # converges
rep200FungiDecontamV2Species_8cancer_VSNM <- rep200FungiDecontamV2Species_8cancer_VSNM_Obj$snmData

# NOT CURRENTLY RUN:
## Summarized 8 cancer types - biological variables: disease_type + sample_type | technical variables: data_submitting_center_label + experimental_strategy
# NOTE: Phylum, Class, Order, and Family levels were attempted but did not converge, so they are not shown here
# rep200FungiDecontamV2Genus_8cancer_VSNM_Obj_CT <- vsnmFunctionTCGA(rep200FungiDecontamV2Genus_8cancer, qcMetadata = metaQiitaCombined_Nonzero_DecontamV2_8cancer, cancerTypeFlag = TRUE) # converges
# rep200FungiDecontamV2Genus_8cancer_VSNM_CT <- rep200FungiDecontamV2Genus_8cancer_VSNM_Obj_CT$snmData
# rep200FungiDecontamV2Species_8cancer_VSNM_Obj_CT <- vsnmFunctionTCGA(rep200FungiDecontamV2Species_8cancer, qcMetadata = metaQiitaCombined_Nonzero_DecontamV2_8cancer, cancerTypeFlag = TRUE) # converges
# rep200FungiDecontamV2Species_8cancer_VSNM_CT <- rep200FungiDecontamV2Species_8cancer_VSNM_Obj_CT$snmData

#----------------------------------------------------------#
# Save data for running PVCA on each subsetted data to examine batch correction quality
#----------------------------------------------------------#

save(metaQiitaCombined_Nonzero_DecontamV2,
     rep200FungiDecontamV2OrderVSNM_Obj,
     rep200FungiDecontamV2FamilyVSNM_Obj,
     rep200FungiDecontamV2GenusVSNM_Obj,
     rep200FungiDecontamV2SpeciesVSNM_Obj,
     # Full TCGA above, Weizmann matched cancers below
     metaQiitaCombined_Nonzero_DecontamV2_8cancer,
     rep200FungiDecontamV2Order_8cancer_VSNM_Obj,
     rep200FungiDecontamV2Family_8cancer_VSNM_Obj,
     rep200FungiDecontamV2Genus_8cancer_VSNM_Obj,
     rep200FungiDecontamV2Species_8cancer_VSNM_Obj,
     file = "Interim_data/data_for_pvca_tcga_taxa_levels_decontamV2_2Apr22.RData")

#--------------------------------------------------------------------------------------------------------------------#
# Separate WGS and RNA-Seq data for ML (to compare performance between them after batch correcting)
#--------------------------------------------------------------------------------------------------------------------#
load("Interim_data/data_for_pvca_tcga_taxa_levels_decontamV2_2Apr22.RData")

metaQiitaCombined_Nonzero_DecontamV2_WGS <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(experimental_strategy == "WGS") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_RNA <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()

rep200FungiDecontamV2SpeciesVSNM <- rep200FungiDecontamV2SpeciesVSNM_Obj$snmData

rep200FungiDecontamV2SpeciesVSNM_WGS <- rep200FungiDecontamV2SpeciesVSNM[rownames(metaQiitaCombined_Nonzero_DecontamV2_WGS),]
rep200FungiDecontamV2SpeciesVSNM_RNA <- rep200FungiDecontamV2SpeciesVSNM[rownames(metaQiitaCombined_Nonzero_DecontamV2_RNA),]

save(metaQiitaCombined_Nonzero_DecontamV2,
     metaQiitaCombined_Nonzero_DecontamV2_WGS,
     metaQiitaCombined_Nonzero_DecontamV2_RNA,
     rep200FungiDecontamV2SpeciesVSNM,
     rep200FungiDecontamV2SpeciesVSNM_WGS,
     rep200FungiDecontamV2SpeciesVSNM_RNA,
     file = "Interim_data/data_for_ml_tcga_wgs_vs_rna_sample_type_VSNM_decontamV2_2Apr22.RData")

# Scripts: S06R

#--------------------Format data and regression plot--------------------#
source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_full_WGS_RNA_v2 <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_full_wgs_rna_vsnm_STonly_ALL_DecontamV2_2Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_full_WGS_RNA_v2$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_full_WGS_RNA_v2$diseaseType,"abbrev"]
mlPerfAll10k_full_WGS_RNA_v2 <- mlPerfAll10k_full_WGS_RNA_v2[,!(colnames(mlPerfAll10k_full_WGS_RNA_v2) == "X")]
colnames(mlPerfAll10k_full_WGS_RNA_v2)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_full_WGS_RNA_v2$nullAUPR <- ifelse(mlPerfAll10k_full_WGS_RNA_v2$minorityClassName == "SolidTissueNormal",
                                             yes=mlPerfAll10k_full_WGS_RNA_v2$majorityClassSize/(mlPerfAll10k_full_WGS_RNA_v2$minorityClassSize+mlPerfAll10k_full_WGS_RNA_v2$majorityClassSize),
                                             no=mlPerfAll10k_full_WGS_RNA_v2$minorityClassSize/(mlPerfAll10k_full_WGS_RNA_v2$minorityClassSize+mlPerfAll10k_full_WGS_RNA_v2$majorityClassSize))
mlPerfAll10k_full_WGS_RNA_v2$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
mlPerfAll10k_full_WGS_RNA_v2$datasetName[mlPerfAll10k_full_WGS_RNA_v2$datasetName == "rep200FungiDecontamV2SpeciesVSNM"] <- "WGS+RNA-Seq Species (Decontaminated)"
mlPerfAll10k_full_WGS_RNA_v2$datasetName[mlPerfAll10k_full_WGS_RNA_v2$datasetName == "rep200FungiDecontamV2SpeciesVSNM_WGS"] <- "WGS Species (Decontaminated)"
mlPerfAll10k_full_WGS_RNA_v2$datasetName[mlPerfAll10k_full_WGS_RNA_v2$datasetName == "rep200FungiDecontamV2SpeciesVSNM_RNA"] <- "RNA-Seq Species (Decontaminated)"

mlPerfAll10k_full_WGS_RNA_v2$datasetName <- factor(mlPerfAll10k_full_WGS_RNA_v2$datasetName,
                                                levels = c("WGS+RNA-Seq Species (Decontaminated)",
                                                           "WGS Species (Decontaminated)",
                                                           "RNA-Seq Species (Decontaminated)"))
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
mlPerfAll10k_full_WGS_RNA_v2 %>%
  filter(sampleType == "Primary Tumor") %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor | 1 Vs All | Species | WGS vs. RNA") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_wgs_vs_rna_VSNM_STonly_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

#-------------------------Plot primary tumor vs. adjacent normal performance-------------------------#
mlPerfAll10k_full_WGS_RNA_v2 %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species | WGS vs. RNA") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_STN_wgs_vs_rna_VSNM_STonly_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

#-------------------------Linear regression of WGS-only OR RNA-Seq-only vs. VSNM data-------------------------#

perfSummarizedFullvsWGSvsRNA <- mlPerfAll10k_full_WGS_RNA_v2 %>% distinct() %>% 
  filter(sampleType == "Primary Tumor") %>% droplevels() %>%
  group_by(sampleType, diseaseType, datasetName) %>% 
  mutate(avgROC = mean(AUROC), stdROC=sd(AUROC), avgPR = mean(AUPR), stdPR=sd(AUPR),) %>%
  select(sampleType, diseaseType, datasetName, abbrev, avgROC, stdROC, avgPR, stdPR, minorityClassSize, majorityClassSize) %>% 
  unique() %>% data.frame()

perfSummarizedFullvsWGSvsRNA_Full <- perfSummarizedFullvsWGSvsRNA %>%
  filter(datasetName=="WGS+RNA-Seq Species (Decontaminated)") %>% droplevels()
perfSummarizedFullvsWGSvsRNA_WGSonly <- perfSummarizedFullvsWGSvsRNA %>%
  filter(datasetName=="WGS Species (Decontaminated)") %>% droplevels() %>%
  rename(avgROC_WGS = avgROC, stdROC_WGS = stdROC, avgPR_WGS = avgPR, stdPR_WGS = stdPR)
perfSummarizedFullvsWGSvsRNA_RNAonly <- perfSummarizedFullvsWGSvsRNA %>%
  filter(datasetName=="RNA-Seq Species (Decontaminated)") %>% droplevels() %>%
  rename(avgROC_RNA = avgROC, stdROC_RNA = stdROC, avgPR_RNA = avgPR, stdPR_RNA = stdPR)

perfSummarizedFullvsWGSvsRNA_Full_Joined <- perfSummarizedFullvsWGSvsRNA_Full %>%
  left_join(perfSummarizedFullvsWGSvsRNA_WGSonly[,grep("abbrev|ROC|PR",colnames(perfSummarizedFullvsWGSvsRNA_WGSonly))], by = "abbrev") %>%
  left_join(perfSummarizedFullvsWGSvsRNA_RNAonly[,grep("abbrev|ROC|PR",colnames(perfSummarizedFullvsWGSvsRNA_RNAonly))], by = "abbrev")

require(ggpmisc)
require(ggrepel)
# WGS vs. Full | AUROC
ggplot(perfSummarizedFullvsWGSvsRNA_Full_Joined, 
       aes(x = avgROC, y =  avgROC_WGS, label = abbrev)) +
  geom_point(alpha = 0.8) +
  coord_equal() + xlim(c(0, 1)) + ylim(c(0,1)) +
  theme_bw() +
  geom_abline(linetype = 2) +
  labs(x = "Average classifier performances on full data\n(primary tumor 1 cancer type vs all others)",
       y = "Average classifier performances on WGS-only data\n(primary tumor 1 cancer type vs all others)",
       title = "WGS-only vs. full data | AUROC") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  geom_ribbon(stat='smooth', method = "lm", fill = "orange", se=TRUE, alpha=0.2, fullrange = TRUE) +
  geom_line(stat='smooth', method = "lm", color = "darkorange", alpha=1, fullrange = TRUE) +
  geom_text_repel(force = 120, size=3) +
  stat_poly_eq(formula = y ~ x, label.x = 0.9, label.y = 0.05,
               aes(label = paste(..eq.label.., sep = "~~~")), 
               parse = TRUE)
ggsave(filename = "Figures/Other_Figures/regr_auroc_tcga_wgs_vs_full_PT_VSNM_STonly.pdf",
         dpi = "retina", units = "in", height = 5, width = 5)
# WGS vs. Full | AUPR
ggplot(perfSummarizedFullvsWGSvsRNA_Full_Joined, 
       aes(x = avgPR, y =  avgPR_WGS, label = abbrev)) +
  geom_point(alpha = 0.8) +
  coord_equal() + xlim(c(0, 1)) + ylim(c(0,1)) +
  theme_bw() +
  geom_abline(linetype = 2) +
  labs(x = "Average classifier performances on full data\n(primary tumor 1 cancer type vs all others)",
       y = "Average classifier performances on WGS-only data\n(primary tumor 1 cancer type vs all others)",
       title = "WGS-only vs. full data | AUPR") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  geom_ribbon(stat='smooth', method = "lm", fill = "orange", se=TRUE, alpha=0.2, fullrange = TRUE) +
  geom_line(stat='smooth', method = "lm", color = "darkorange", alpha=1, fullrange = TRUE) +
  geom_text_repel(size=3) +
  stat_poly_eq(formula = y ~ x, label.x = 0.9, label.y = 0.05,
               aes(label = paste(..eq.label.., sep = "~~~")), 
               parse = TRUE)
ggsave(filename = "Figures/Other_Figures/regr_aupr_tcga_wgs_vs_full_PT_VSNM_STonly.pdf",
         dpi = "retina", units = "in", height = 5, width = 5)
# RNA vs. Full | AUROC
ggplot(perfSummarizedFullvsWGSvsRNA_Full_Joined, 
       aes(x = avgROC, y =  avgROC_RNA, label = abbrev)) +
  geom_point(alpha = 0.8) +
  coord_equal() + xlim(c(0, 1)) + ylim(c(0,1)) +
  theme_bw() +
  geom_abline(linetype = 2) +
  labs(x = "Average classifier performances on full data\n(primary tumor 1 cancer type vs all others)",
       y = "Average classifier performances on RNA-only data\n(primary tumor 1 cancer type vs all others)",
       title = "RNA-only vs. full data | AUROC") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  geom_ribbon(stat='smooth', method = "lm", fill = "orange", se=TRUE, alpha=0.2, fullrange = TRUE) +
  geom_line(stat='smooth', method = "lm", color = "darkorange", alpha=1, fullrange = TRUE) +
  geom_text_repel(force = 120, size=3) +
  stat_poly_eq(formula = y ~ x, label.x = 0.7, label.y = 0.05,
               aes(label = paste(..eq.label.., sep = "~~~")), 
               parse = TRUE)
ggsave(filename = "Figures/Other_Figures/regr_auroc_tcga_rna_vs_full_PT_VSNM_STonly.pdf",
         dpi = "retina", units = "in", height = 5, width = 5)
# RNA vs. Full | AUPR
ggplot(perfSummarizedFullvsWGSvsRNA_Full_Joined, 
       aes(x = avgPR, y =  avgPR_RNA, label = abbrev)) +
  geom_point(alpha = 0.8) +
  coord_equal() + xlim(c(0, 1)) + ylim(c(0,1)) +
  theme_bw() +
  geom_abline(linetype = 2) +
  labs(x = "Average classifier performances on full data\n(primary tumor 1 cancer type vs all others)",
       y = "Average classifier performances on RNA-only data\n(primary tumor 1 cancer type vs all others)",
       title = "RNA-only vs. full data | AUPR") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  geom_ribbon(stat='smooth', method = "lm", fill = "orange", se=TRUE, alpha=0.2, fullrange = TRUE) +
  geom_line(stat='smooth', method = "lm", color = "darkorange", alpha=1, fullrange = TRUE) +
  geom_text_repel(force = 150, size=3) +
  stat_poly_eq(formula = y ~ x, label.x = 0.05, label.y = 0.9,
               aes(label = paste(..eq.label.., sep = "~~~")), 
               parse = TRUE)
ggsave(filename = "Figures/Other_Figures/regr_aupr_tcga_rna_vs_full_PT_VSNM_STonly.pdf",
         dpi = "retina", units = "in", height = 5, width = 5)

#--------------------------------------------------------------------------------------------------------------------#
# Test Hopkins top 20 fungi signature on TCGA blood samples
#--------------------------------------------------------------------------------------------------------------------#
load("Interim_data/topXHopkinsSig_13Nov21.RData")

#---------------------Subset VSNM data---------------------#
load("Interim_data/data_for_pvca_tcga_taxa_levels_decontamV2_2Apr22.RData")

rep200FungiDecontamV2SpeciesVSNM <- rep200FungiDecontamV2SpeciesVSNM_Obj$snmData
rep200FungiDecontamV2SpeciesVSNM_TopX <- rep200FungiDecontamV2SpeciesVSNM[,colnames(rep200FungiDecontamV2SpeciesVSNM) %in%
                                                                            gsub("\\s","_",topXHopkinsSig$species)]
dim(rep200FungiDecontamV2SpeciesVSNM_TopX) # 14495    20

save(metaQiitaCombined_Nonzero_DecontamV2,
     rep200FungiDecontamV2SpeciesVSNM_TopX,
     file = "Interim_data/data_for_ml_tcga_vsnm_topX_2Apr22.RData")

# Scripts: S07R

#---------------------Subset raw data---------------------#
load("Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy_with_coverage_filter_decontamV2_2Apr22.RData")

rep200_HiSeq_Fungi_DecontamV2_TopX_HMS <- rep200_HiSeq_Fungi_DecontamV2_HMS[,rownames(topXHopkinsSig)]
rep200_HiSeq_Fungi_DecontamV2_TopX_BCM <- rep200_HiSeq_Fungi_DecontamV2_BCM[,rownames(topXHopkinsSig)]
rep200_HiSeq_Fungi_DecontamV2_TopX_MDA <- rep200_HiSeq_Fungi_DecontamV2_MDA[,rownames(topXHopkinsSig)]
rep200_HiSeq_Fungi_DecontamV2_TopX_WashU <- rep200_HiSeq_Fungi_DecontamV2_WashU[,rownames(topXHopkinsSig)]
rep200_HiSeq_Fungi_DecontamV2_TopX_Broad_WGS <- rep200_HiSeq_Fungi_DecontamV2_Broad_WGS[,rownames(topXHopkinsSig)]
rep200_HiSeq_Fungi_DecontamV2_TopX_UNC <- rep200_HiSeq_Fungi_DecontamV2_UNC[,rownames(topXHopkinsSig)]
rep200_HiSeq_Fungi_DecontamV2_TopX_CMS <- rep200_HiSeq_Fungi_DecontamV2_CMS[,rownames(topXHopkinsSig)]
rep200_HiSeq_Fungi_DecontamV2_TopX_Broad_RNA <- rep200_HiSeq_Fungi_DecontamV2_Broad_RNA[,rownames(topXHopkinsSig)]

# Out of the above, the following have 1 or more zero sum samples after subsetting to topX features:
# rep200_HiSeq_Fungi_DecontamV2_TopX_Broad_WGS (13 zero sum samples)
# rep200_HiSeq_Fungi_DecontamV2_TopX_UNC (901 zero sum samples)
# rep200_HiSeq_Fungi_DecontamV2_TopX_CMS (18 zero sum samples)
# These will need to have their metadata subsetted
# Broad WGS
rep200_HiSeq_Fungi_DecontamV2_TopX_Broad_WGS_Nonzero <- rep200_HiSeq_Fungi_DecontamV2_TopX_Broad_WGS[-which(rowSums(rep200_HiSeq_Fungi_DecontamV2_TopX_Broad_WGS)==0),]
metaQiitaCombined_Nonzero_DecontamV2_TopX_Broad_WGS_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS[rownames(rep200_HiSeq_Fungi_DecontamV2_TopX_Broad_WGS_Nonzero),])
# UNC
rep200_HiSeq_Fungi_DecontamV2_TopX_UNC_Nonzero <- rep200_HiSeq_Fungi_DecontamV2_TopX_UNC[-which(rowSums(rep200_HiSeq_Fungi_DecontamV2_TopX_UNC)==0),]
metaQiitaCombined_Nonzero_DecontamV2_TopX_UNC_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_UNC[rownames(rep200_HiSeq_Fungi_DecontamV2_TopX_UNC_Nonzero),])
# CMS
rep200_HiSeq_Fungi_DecontamV2_TopX_CMS_Nonzero <- rep200_HiSeq_Fungi_DecontamV2_TopX_CMS[-which(rowSums(rep200_HiSeq_Fungi_DecontamV2_TopX_CMS)==0),]
metaQiitaCombined_Nonzero_DecontamV2_TopX_CMS_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_CMS[rownames(rep200_HiSeq_Fungi_DecontamV2_TopX_CMS_Nonzero),])

save(rep200_HiSeq_Fungi_DecontamV2_TopX_HMS,
     rep200_HiSeq_Fungi_DecontamV2_TopX_BCM,
     rep200_HiSeq_Fungi_DecontamV2_TopX_MDA,
     rep200_HiSeq_Fungi_DecontamV2_TopX_WashU,
     rep200_HiSeq_Fungi_DecontamV2_TopX_Broad_WGS_Nonzero,
     rep200_HiSeq_Fungi_DecontamV2_TopX_UNC_Nonzero,
     rep200_HiSeq_Fungi_DecontamV2_TopX_CMS_Nonzero,
     rep200_HiSeq_Fungi_DecontamV2_TopX_Broad_RNA,
     
     metaQiitaCombined_Nonzero_DecontamV2_HMS,
     metaQiitaCombined_Nonzero_DecontamV2_BCM,
     metaQiitaCombined_Nonzero_DecontamV2_MDA,
     metaQiitaCombined_Nonzero_DecontamV2_WashU,
     metaQiitaCombined_Nonzero_DecontamV2_TopX_Broad_WGS_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_TopX_UNC_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_TopX_CMS_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA,
     file = "Interim_data/data_for_ml_tcga_raw_topX_per_seq_center_2Apr22.RData")

# Scripts: S08R

#----------------------------------------------------------#
# Separate 8 cancer data into WGS and RNA to compare ML performance
#----------------------------------------------------------#

metaQiitaCombined_Nonzero_DecontamV2_8cancer_WGS <- metaQiitaCombined_Nonzero_DecontamV2_8cancer %>% filter(experimental_strategy=="WGS") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_8cancer_RNA <- metaQiitaCombined_Nonzero_DecontamV2_8cancer %>% filter(experimental_strategy=="RNA-Seq") %>% droplevels()

# Subset WGS data
rep200FungiDecontamV2Order_8cancer_VSNM_WGS <- rep200FungiDecontamV2Order_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer_WGS),]
rep200FungiDecontamV2Family_8cancer_VSNM_WGS <- rep200FungiDecontamV2Family_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer_WGS),]
rep200FungiDecontamV2Genus_8cancer_VSNM_WGS <- rep200FungiDecontamV2Genus_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer_WGS),]
rep200FungiDecontamV2Species_8cancer_VSNM_WGS <- rep200FungiDecontamV2Species_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer_WGS),]

# Subset RNA data
rep200FungiDecontamV2Order_8cancer_VSNM_RNA <- rep200FungiDecontamV2Order_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer_RNA),]
rep200FungiDecontamV2Family_8cancer_VSNM_RNA <- rep200FungiDecontamV2Family_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer_RNA),]
rep200FungiDecontamV2Genus_8cancer_VSNM_RNA <- rep200FungiDecontamV2Genus_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer_RNA),]
rep200FungiDecontamV2Species_8cancer_VSNM_RNA <- rep200FungiDecontamV2Species_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer_RNA),]

save(metaQiitaCombined_Nonzero_DecontamV2_8cancer_WGS,
     rep200FungiDecontamV2Order_8cancer_VSNM_WGS,
     rep200FungiDecontamV2Family_8cancer_VSNM_WGS,
     rep200FungiDecontamV2Genus_8cancer_VSNM_WGS,
     rep200FungiDecontamV2Species_8cancer_VSNM_WGS,
     
     metaQiitaCombined_Nonzero_DecontamV2_8cancer_RNA,
     rep200FungiDecontamV2Order_8cancer_VSNM_RNA,
     rep200FungiDecontamV2Family_8cancer_VSNM_RNA,
     rep200FungiDecontamV2Genus_8cancer_VSNM_RNA,
     rep200FungiDecontamV2Species_8cancer_VSNM_RNA,
     file = "Interim_data/data_for_ml_8_cancers_wgs_vs_rna_decontamV2_2Apr22.RData")

# Scripts: S09R

#----------------------------------------------------------#
# Load TCGA data matched to Weizmann by cancer types and shared features to run VSNM
#----------------------------------------------------------#

load("Interim_data/data_tcga_8_cancers_features_matched_to_weizmann_29Mar22.RData")
load("Interim_data/tcga_data_taxa_levels_features_shared_with_Weizmann_29Mar22.RData")

## Full TCGA intersected features
source("00-Functions.R")
# NOTE: Phylum, Class, and Order levels were attempted but did not converge, so they are not shown here
rep200FungiFamilyShared_NonzeroVSNM_Obj <- vsnmFunctionTCGA(rep200FungiFamilyShared_Nonzero, qcMetadata = metaQiitaCombined_Nonzero_FamilyShared) # worked
rep200FungiFamilyShared_NonzeroVSNM <- rep200FungiFamilyShared_NonzeroVSNM_Obj$snmData
rep200FungiGenusShared_NonzeroVSNM_Obj <- vsnmFunctionTCGA(rep200FungiGenusShared_Nonzero, qcMetadata = metaQiitaCombined_Nonzero_GenusShared) # worked
rep200FungiGenusShared_NonzeroVSNM <- rep200FungiGenusShared_NonzeroVSNM_Obj$snmData
rep200FungiSpeciesShared_NonzeroVSNM_Obj <- vsnmFunctionTCGA(rep200FungiSpeciesShared_Nonzero, qcMetadata = metaQiitaCombined_Nonzero_SpeciesShared) # worked
rep200FungiSpeciesShared_NonzeroVSNM <- rep200FungiSpeciesShared_NonzeroVSNM_Obj$snmData

## Partial TCGA intersected by cancer types and by features
# NOTE: Class and Order levels were attempted but did not converge, so they are not shown here
rep200FungiPhylumShared_8cancer_Nonzero_VSNM_Obj <- vsnmFunctionTCGA(rep200FungiPhylumShared_8cancer_Nonzero, qcMetadata = metaQiitaCombined_Nonzero_8cancer_PhylumShared) # doesn't work
rep200FungiPhylumShared_8cancer_Nonzero_VSNM <- rep200FungiPhylumShared_8cancer_Nonzero_VSNM_Obj$snmData
rep200FungiFamilyShared_8cancer_Nonzero_VSNM_Obj <- vsnmFunctionTCGA(rep200FungiFamilyShared_8cancer_Nonzero, qcMetadata = metaQiitaCombined_Nonzero_8cancer_FamilyShared) # worked
rep200FungiFamilyShared_8cancer_Nonzero_VSNM <- rep200FungiFamilyShared_8cancer_Nonzero_VSNM_Obj$snmData
rep200FungiGenusShared_8cancer_Nonzero_VSNM_Obj <- vsnmFunctionTCGA(rep200FungiGenusShared_8cancer_Nonzero, qcMetadata = metaQiitaCombined_Nonzero_8cancer_GenusShared) # worked
rep200FungiGenusShared_8cancer_Nonzero_VSNM <- rep200FungiGenusShared_8cancer_Nonzero_VSNM_Obj$snmData
rep200FungiSpeciesShared_8cancer_Nonzero_VSNM_Obj <- vsnmFunctionTCGA(rep200FungiSpeciesShared_8cancer_Nonzero, qcMetadata = metaQiitaCombined_Nonzero_8cancer_SpeciesShared) # worked
rep200FungiSpeciesShared_8cancer_Nonzero_VSNM <- rep200FungiSpeciesShared_8cancer_Nonzero_VSNM_Obj$snmData

#----------------------------------------------------------#
# Subset full TCGA data and 8 WIS cancers to features with ≥1% aggregate coverage to run VSNM
#----------------------------------------------------------#

# metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts, 
# rep200Data_WGS_RNA_Matched,
# rep200Data_WGS_RNA_Matched_Bacteria,
# rep200Data_WGS_RNA_Matched_Fungi,
load("Interim_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_25Mar22.RData")
coverageFungiAllSamples <- read.csv("Input_data/fungi_filt_updated_29Sep21_coverage_all_wgs_and_rna_output.csv", stringsAsFactors = FALSE)
coverageFungiAllSamples_1Percent <- coverageFungiAllSamples %>% filter(coverage_ratio >= 0.01)
coverageFungiAllSamples_1Percent_OGUs <- coverageFungiAllSamples_1Percent$gotu
#----------------------------Full TCGA data----------------------------#
# Subset count data
rep200Data_WGS_RNA_Matched_Fungi_Cov <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2), coverageFungiAllSamples_1Percent_OGUs] 
# Remove zero sum samples
rep200Data_WGS_RNA_Matched_Fungi_Cov_Nonzero <- rep200Data_WGS_RNA_Matched_Fungi_Cov[rowSums(rep200Data_WGS_RNA_Matched_Fungi_Cov) != 0,]
metaQiitaCombined_Nonzero_DecontamV2_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2[rowSums(rep200Data_WGS_RNA_Matched_Fungi_Cov) != 0,])
# Filter out samples to enable Voom-SNM to converge
metaQiitaCombined_Nonzero_DecontamV2_Cov_Nonzero_Filt <- metaQiitaCombined_Nonzero_DecontamV2_Cov_Nonzero %>% 
  filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal", "Blood Derived Normal")) %>% droplevels()
rep200Data_WGS_RNA_Matched_Fungi_Cov_Nonzero_Filt <- rep200Data_WGS_RNA_Matched_Fungi_Cov_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_Cov_Nonzero_Filt),]
# Voom-SNM
rep200FungiCovSpeciesVSNM_Obj <- vsnmFunctionTCGA(rep200Data_WGS_RNA_Matched_Fungi_Cov_Nonzero_Filt, qcMetadata = metaQiitaCombined_Nonzero_DecontamV2_Cov_Nonzero_Filt) # worked
rep200FungiCovSpeciesVSNM <- rep200FungiCovSpeciesVSNM_Obj$snmData
#----------------------------TCGA matching WIS 8 cancers data----------------------------#
# Subset count data
rep200Fungi_8cancer_Cov <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer), coverageFungiAllSamples_1Percent_OGUs] 
# Remove zero sum samples
rep200Fungi_8cancer_Cov_Nonzero <- rep200Fungi_8cancer_Cov[rowSums(rep200Fungi_8cancer_Cov) != 0,]
metaQiitaCombined_Nonzero_DecontamV2_8cancer_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_8cancer[rowSums(rep200Fungi_8cancer_Cov) != 0,])
# Filter out samples to enable Voom-SNM to converge
metaQiitaCombined_Nonzero_DecontamV2_8cancer_Cov_Nonzero_Filt <- metaQiitaCombined_Nonzero_DecontamV2_8cancer_Cov_Nonzero %>% 
  filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal", "Blood Derived Normal")) %>% droplevels()
rep200Fungi_8cancer_Cov_Nonzero_Filt <- rep200Fungi_8cancer_Cov_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_8cancer_Cov_Nonzero_Filt),]
# Voom-SNM
rep200FungiSpecies_8cancer_Cov_Nonzero_VSNM_Obj <- vsnmFunctionTCGA(rep200Fungi_8cancer_Cov_Nonzero_Filt, qcMetadata = metaQiitaCombined_Nonzero_DecontamV2_8cancer_Cov_Nonzero_Filt) # worked
rep200FungiSpecies_8cancer_Cov_Nonzero_VSNM <- rep200FungiSpecies_8cancer_Cov_Nonzero_VSNM_Obj$snmData

#----------------------------------------------------------#
# Save data for running machine learning
#----------------------------------------------------------#

save(metaQiitaCombined_Nonzero_DecontamV2_8cancer,
     rep200FungiDecontamV2Order_8cancer_VSNM,
     rep200FungiDecontamV2Family_8cancer_VSNM,
     rep200FungiDecontamV2Genus_8cancer_VSNM,
     rep200FungiDecontamV2Species_8cancer_VSNM,
     file = "Interim_data/data_for_ml_8_cancers_decontamV2_2Apr22.RData")

# Scripts: S10R

save(metaQiitaCombined_Nonzero_DecontamV2,
     rep200FungiDecontamV2OrderVSNM,
     rep200FungiDecontamV2FamilyVSNM,
     rep200FungiDecontamV2GenusVSNM,
     rep200FungiDecontamV2SpeciesVSNM,
     snmDataOGUFungiDecontamV2,
     file = "Interim_data/data_for_ml_tcga_decontamV2_2Apr22.RData")

# Scripts: S11R

save(rep200FungiFamilyShared_NonzeroVSNM,
     metaQiitaCombined_Nonzero_FamilyShared,
     rep200FungiGenusShared_NonzeroVSNM,
     metaQiitaCombined_Nonzero_GenusShared,
     rep200FungiSpeciesShared_NonzeroVSNM,
     metaQiitaCombined_Nonzero_SpeciesShared,
     file = "Interim_data/data_for_ml_tcga_shared_features_decontamV2_2Apr22.RData")

# Scripts: S12R

save(rep200FungiPhylumShared_8cancer_Nonzero_VSNM,
     metaQiitaCombined_Nonzero_8cancer_PhylumShared,
     rep200FungiFamilyShared_8cancer_Nonzero_VSNM,
     metaQiitaCombined_Nonzero_8cancer_FamilyShared,
     rep200FungiGenusShared_8cancer_Nonzero_VSNM,
     metaQiitaCombined_Nonzero_8cancer_GenusShared,
     rep200FungiSpeciesShared_8cancer_Nonzero_VSNM,
     metaQiitaCombined_Nonzero_8cancer_SpeciesShared,
     file = "Interim_data/data_for_ml_8_cancers_shared_features_decontamV2_2Apr22.RData")

# Scripts: S13R

save(metaQiitaCombined_Nonzero_DecontamV2_Cov_Nonzero_Filt,
     rep200FungiCovSpeciesVSNM,
     file = "Interim_data/data_for_ml_tcga_high_coverage_2Apr22.RData")

# Scripts: S14R

save(metaQiitaCombined_Nonzero_DecontamV2_8cancer_Cov_Nonzero_Filt,
     rep200FungiSpecies_8cancer_Cov_Nonzero_VSNM,
     file = "Interim_data/data_for_ml_tcga_8cancers_high_coverage_2Apr22.RData")

# Scripts: S15R

#----------------------------------------------------------#
# Subset to EukDetect and HMP taxa at species level
# using VSNM and raw data
#----------------------------------------------------------#

load("Interim_data/shared_fungi_eukdetect_filtAndUnfilt_1Apr22.RData", verbose = TRUE)
load("Interim_data/decontamResultsV2_25Mar22.RData", verbose = TRUE)
decontamResultsV2_HMP <- decontamResultsV2 %>% filter(in_hmp_gut_mycobiome_metagenomic_data == "YES") %>% droplevels()
decontamResultsV2_HMP$species <- gsub(" ","_",decontamResultsV2_HMP$species)

#-------------------VSNM data subset--------------------#

rep200FungiDecontamV2SpeciesVSNM_EukDetectFiltered <- rep200FungiDecontamV2SpeciesVSNM[,colnames(rep200FungiDecontamV2SpeciesVSNM) %in%
                                                                                         sharedSpeciesEukDetectFiltered$intersectedTaxa]
rep200FungiDecontamV2SpeciesVSNM_EukDetectUnfiltered <- rep200FungiDecontamV2SpeciesVSNM[,colnames(rep200FungiDecontamV2SpeciesVSNM) %in%
                                                                                         sharedSpeciesEukDetectUnfiltered$intersectedTaxa]
rep200FungiDecontamV2SpeciesVSNM_HMP <- rep200FungiDecontamV2SpeciesVSNM[,colnames(rep200FungiDecontamV2SpeciesVSNM) %in%
                                                                           decontamResultsV2_HMP$species]
dim(rep200FungiDecontamV2SpeciesVSNM_EukDetectFiltered) # 14495    20
dim(rep200FungiDecontamV2SpeciesVSNM_EukDetectUnfiltered) # 14495    38
dim(rep200FungiDecontamV2SpeciesVSNM_HMP) # 14495    135

#-------------------Per sequencing center subset--------------------#
## HMS -- subset features
df_psRep200_HiSeq_Fungi_HMS_species_EukDetectFiltered <- df_psRep200_HiSeq_Fungi_HMS_species[,colnames(df_psRep200_HiSeq_Fungi_HMS_species) %in% sharedSpeciesEukDetectFiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_HMS_species_EukDetectUnfiltered <- df_psRep200_HiSeq_Fungi_HMS_species[,colnames(df_psRep200_HiSeq_Fungi_HMS_species) %in% sharedSpeciesEukDetectUnfiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_HMS_species_HMP <- df_psRep200_HiSeq_Fungi_HMS_species[,colnames(df_psRep200_HiSeq_Fungi_HMS_species) %in% decontamResultsV2_HMP$species]
## HMS -- remove zero sum samples
df_psRep200_HiSeq_Fungi_HMS_species_EukDetectFiltered_Nonzero <- df_psRep200_HiSeq_Fungi_HMS_species_EukDetectFiltered[!rowSums(df_psRep200_HiSeq_Fungi_HMS_species_EukDetectFiltered)==0,]
df_psRep200_HiSeq_Fungi_HMS_species_EukDetectUnfiltered_Nonzero <- df_psRep200_HiSeq_Fungi_HMS_species_EukDetectUnfiltered[!rowSums(df_psRep200_HiSeq_Fungi_HMS_species_EukDetectUnfiltered)==0,]
df_psRep200_HiSeq_Fungi_HMS_species_HMP_Nonzero <- df_psRep200_HiSeq_Fungi_HMS_species_HMP[!rowSums(df_psRep200_HiSeq_Fungi_HMS_species_HMP)==0,]
## HMS - metadata
metaQiitaCombined_Nonzero_DecontamV2_HMS_species_EukDetectFiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_HMS[rownames(df_psRep200_HiSeq_Fungi_HMS_species_EukDetectFiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_HMS_species_EukDetectUnfiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_HMS[rownames(df_psRep200_HiSeq_Fungi_HMS_species_EukDetectUnfiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_HMS_species_HMP_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_HMS[rownames(df_psRep200_HiSeq_Fungi_HMS_species_HMP_Nonzero),])

## BCM -- subset features
df_psRep200_HiSeq_Fungi_BCM_species_EukDetectFiltered <- df_psRep200_HiSeq_Fungi_BCM_species[,colnames(df_psRep200_HiSeq_Fungi_BCM_species) %in% sharedSpeciesEukDetectFiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_BCM_species_EukDetectUnfiltered <- df_psRep200_HiSeq_Fungi_BCM_species[,colnames(df_psRep200_HiSeq_Fungi_BCM_species) %in% sharedSpeciesEukDetectUnfiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_BCM_species_HMP <- df_psRep200_HiSeq_Fungi_BCM_species[,colnames(df_psRep200_HiSeq_Fungi_BCM_species) %in% decontamResultsV2_HMP$species]
## BCM -- remove zero sum samples
df_psRep200_HiSeq_Fungi_BCM_species_EukDetectFiltered_Nonzero <- df_psRep200_HiSeq_Fungi_BCM_species_EukDetectFiltered[!rowSums(df_psRep200_HiSeq_Fungi_BCM_species_EukDetectFiltered)==0,]
df_psRep200_HiSeq_Fungi_BCM_species_EukDetectUnfiltered_Nonzero <- df_psRep200_HiSeq_Fungi_BCM_species_EukDetectUnfiltered[!rowSums(df_psRep200_HiSeq_Fungi_BCM_species_EukDetectUnfiltered)==0,]
df_psRep200_HiSeq_Fungi_BCM_species_HMP_Nonzero <- df_psRep200_HiSeq_Fungi_BCM_species_HMP[!rowSums(df_psRep200_HiSeq_Fungi_BCM_species_HMP)==0,]
## BCM - metadata
metaQiitaCombined_Nonzero_DecontamV2_BCM_species_EukDetectFiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_BCM[rownames(df_psRep200_HiSeq_Fungi_BCM_species_EukDetectFiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_BCM_species_EukDetectUnfiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_BCM[rownames(df_psRep200_HiSeq_Fungi_BCM_species_EukDetectUnfiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_BCM_species_HMP_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_BCM[rownames(df_psRep200_HiSeq_Fungi_BCM_species_HMP_Nonzero),])

## MDA -- subset features
df_psRep200_HiSeq_Fungi_MDA_species_EukDetectFiltered <- df_psRep200_HiSeq_Fungi_MDA_species[,colnames(df_psRep200_HiSeq_Fungi_MDA_species) %in% sharedSpeciesEukDetectFiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_MDA_species_EukDetectUnfiltered <- df_psRep200_HiSeq_Fungi_MDA_species[,colnames(df_psRep200_HiSeq_Fungi_MDA_species) %in% sharedSpeciesEukDetectUnfiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_MDA_species_HMP <- df_psRep200_HiSeq_Fungi_MDA_species[,colnames(df_psRep200_HiSeq_Fungi_MDA_species) %in% decontamResultsV2_HMP$species]
## MDA -- remove zero sum samples
df_psRep200_HiSeq_Fungi_MDA_species_EukDetectFiltered_Nonzero <- df_psRep200_HiSeq_Fungi_MDA_species_EukDetectFiltered[!rowSums(df_psRep200_HiSeq_Fungi_MDA_species_EukDetectFiltered)==0,]
df_psRep200_HiSeq_Fungi_MDA_species_EukDetectUnfiltered_Nonzero <- df_psRep200_HiSeq_Fungi_MDA_species_EukDetectUnfiltered[!rowSums(df_psRep200_HiSeq_Fungi_MDA_species_EukDetectUnfiltered)==0,]
df_psRep200_HiSeq_Fungi_MDA_species_HMP_Nonzero <- df_psRep200_HiSeq_Fungi_MDA_species_HMP[!rowSums(df_psRep200_HiSeq_Fungi_MDA_species_HMP)==0,]
## MDA - metadata
metaQiitaCombined_Nonzero_DecontamV2_MDA_species_EukDetectFiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_MDA[rownames(df_psRep200_HiSeq_Fungi_MDA_species_EukDetectFiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_MDA_species_EukDetectUnfiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_MDA[rownames(df_psRep200_HiSeq_Fungi_MDA_species_EukDetectUnfiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_MDA_species_HMP_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_MDA[rownames(df_psRep200_HiSeq_Fungi_MDA_species_HMP_Nonzero),])

## WashU -- subset features
df_psRep200_HiSeq_Fungi_WashU_species_EukDetectFiltered <- df_psRep200_HiSeq_Fungi_WashU_species[,colnames(df_psRep200_HiSeq_Fungi_WashU_species) %in% sharedSpeciesEukDetectFiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_WashU_species_EukDetectUnfiltered <- df_psRep200_HiSeq_Fungi_WashU_species[,colnames(df_psRep200_HiSeq_Fungi_WashU_species) %in% sharedSpeciesEukDetectUnfiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_WashU_species_HMP <- df_psRep200_HiSeq_Fungi_WashU_species[,colnames(df_psRep200_HiSeq_Fungi_WashU_species) %in% decontamResultsV2_HMP$species]
## WashU -- remove zero sum samples
df_psRep200_HiSeq_Fungi_WashU_species_EukDetectFiltered_Nonzero <- df_psRep200_HiSeq_Fungi_WashU_species_EukDetectFiltered[!rowSums(df_psRep200_HiSeq_Fungi_WashU_species_EukDetectFiltered)==0,]
df_psRep200_HiSeq_Fungi_WashU_species_EukDetectUnfiltered_Nonzero <- df_psRep200_HiSeq_Fungi_WashU_species_EukDetectUnfiltered[!rowSums(df_psRep200_HiSeq_Fungi_WashU_species_EukDetectUnfiltered)==0,]
df_psRep200_HiSeq_Fungi_WashU_species_HMP_Nonzero <- df_psRep200_HiSeq_Fungi_WashU_species_HMP[!rowSums(df_psRep200_HiSeq_Fungi_WashU_species_HMP)==0,]
## WashU - metadata
metaQiitaCombined_Nonzero_DecontamV2_WashU_species_EukDetectFiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_WashU[rownames(df_psRep200_HiSeq_Fungi_WashU_species_EukDetectFiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_WashU_species_EukDetectUnfiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_WashU[rownames(df_psRep200_HiSeq_Fungi_WashU_species_EukDetectUnfiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_WashU_species_HMP_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_WashU[rownames(df_psRep200_HiSeq_Fungi_WashU_species_HMP_Nonzero),])

## Broad_WGS -- subset features
df_psRep200_HiSeq_Fungi_Broad_WGS_species_EukDetectFiltered <- df_psRep200_HiSeq_Fungi_Broad_WGS_species[,colnames(df_psRep200_HiSeq_Fungi_Broad_WGS_species) %in% sharedSpeciesEukDetectFiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_Broad_WGS_species_EukDetectUnfiltered <- df_psRep200_HiSeq_Fungi_Broad_WGS_species[,colnames(df_psRep200_HiSeq_Fungi_Broad_WGS_species) %in% sharedSpeciesEukDetectUnfiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_Broad_WGS_species_HMP <- df_psRep200_HiSeq_Fungi_Broad_WGS_species[,colnames(df_psRep200_HiSeq_Fungi_Broad_WGS_species) %in% decontamResultsV2_HMP$species]
## Broad_WGS -- remove zero sum samples
df_psRep200_HiSeq_Fungi_Broad_WGS_species_EukDetectFiltered_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_WGS_species_EukDetectFiltered[!rowSums(df_psRep200_HiSeq_Fungi_Broad_WGS_species_EukDetectFiltered)==0,]
df_psRep200_HiSeq_Fungi_Broad_WGS_species_EukDetectUnfiltered_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_WGS_species_EukDetectUnfiltered[!rowSums(df_psRep200_HiSeq_Fungi_Broad_WGS_species_EukDetectUnfiltered)==0,]
df_psRep200_HiSeq_Fungi_Broad_WGS_species_HMP_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_WGS_species_HMP[!rowSums(df_psRep200_HiSeq_Fungi_Broad_WGS_species_HMP)==0,]
## Broad_WGS - metadata
metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_species_EukDetectFiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Broad_WGS_species_EukDetectFiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_species_EukDetectUnfiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Broad_WGS_species_EukDetectUnfiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_species_HMP_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Broad_WGS_species_HMP_Nonzero),])

## UNC -- subset features
df_psRep200_HiSeq_Fungi_UNC_species_EukDetectFiltered <- df_psRep200_HiSeq_Fungi_UNC_species[,colnames(df_psRep200_HiSeq_Fungi_UNC_species) %in% sharedSpeciesEukDetectFiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_UNC_species_EukDetectUnfiltered <- df_psRep200_HiSeq_Fungi_UNC_species[,colnames(df_psRep200_HiSeq_Fungi_UNC_species) %in% sharedSpeciesEukDetectUnfiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_UNC_species_HMP <- df_psRep200_HiSeq_Fungi_UNC_species[,colnames(df_psRep200_HiSeq_Fungi_UNC_species) %in% decontamResultsV2_HMP$species]
## UNC -- remove zero sum samples
df_psRep200_HiSeq_Fungi_UNC_species_EukDetectFiltered_Nonzero <- df_psRep200_HiSeq_Fungi_UNC_species_EukDetectFiltered[!rowSums(df_psRep200_HiSeq_Fungi_UNC_species_EukDetectFiltered)==0,]
df_psRep200_HiSeq_Fungi_UNC_species_EukDetectUnfiltered_Nonzero <- df_psRep200_HiSeq_Fungi_UNC_species_EukDetectUnfiltered[!rowSums(df_psRep200_HiSeq_Fungi_UNC_species_EukDetectUnfiltered)==0,]
df_psRep200_HiSeq_Fungi_UNC_species_HMP_Nonzero <- df_psRep200_HiSeq_Fungi_UNC_species_HMP[!rowSums(df_psRep200_HiSeq_Fungi_UNC_species_HMP)==0,]
## UNC - metadata
metaQiitaCombined_Nonzero_DecontamV2_UNC_species_EukDetectFiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_UNC[rownames(df_psRep200_HiSeq_Fungi_UNC_species_EukDetectFiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_UNC_species_EukDetectUnfiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_UNC[rownames(df_psRep200_HiSeq_Fungi_UNC_species_EukDetectUnfiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_UNC_species_HMP_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_UNC[rownames(df_psRep200_HiSeq_Fungi_UNC_species_HMP_Nonzero),])

## CMS -- subset features
df_psRep200_HiSeq_Fungi_CMS_species_EukDetectFiltered <- df_psRep200_HiSeq_Fungi_CMS_species[,colnames(df_psRep200_HiSeq_Fungi_CMS_species) %in% sharedSpeciesEukDetectFiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_CMS_species_EukDetectUnfiltered <- df_psRep200_HiSeq_Fungi_CMS_species[,colnames(df_psRep200_HiSeq_Fungi_CMS_species) %in% sharedSpeciesEukDetectUnfiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_CMS_species_HMP <- df_psRep200_HiSeq_Fungi_CMS_species[,colnames(df_psRep200_HiSeq_Fungi_CMS_species) %in% decontamResultsV2_HMP$species]
## CMS -- remove zero sum samples
df_psRep200_HiSeq_Fungi_CMS_species_EukDetectFiltered_Nonzero <- df_psRep200_HiSeq_Fungi_CMS_species_EukDetectFiltered[!rowSums(df_psRep200_HiSeq_Fungi_CMS_species_EukDetectFiltered)==0,]
df_psRep200_HiSeq_Fungi_CMS_species_EukDetectUnfiltered_Nonzero <- df_psRep200_HiSeq_Fungi_CMS_species_EukDetectUnfiltered[!rowSums(df_psRep200_HiSeq_Fungi_CMS_species_EukDetectUnfiltered)==0,]
df_psRep200_HiSeq_Fungi_CMS_species_HMP_Nonzero <- df_psRep200_HiSeq_Fungi_CMS_species_HMP[!rowSums(df_psRep200_HiSeq_Fungi_CMS_species_HMP)==0,]
## CMS - metadata
metaQiitaCombined_Nonzero_DecontamV2_CMS_species_EukDetectFiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_CMS[rownames(df_psRep200_HiSeq_Fungi_CMS_species_EukDetectFiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_CMS_species_EukDetectUnfiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_CMS[rownames(df_psRep200_HiSeq_Fungi_CMS_species_EukDetectUnfiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_CMS_species_HMP_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_CMS[rownames(df_psRep200_HiSeq_Fungi_CMS_species_HMP_Nonzero),])

## Broad_RNA -- subset features
df_psRep200_HiSeq_Fungi_Broad_RNA_species_EukDetectFiltered <- df_psRep200_HiSeq_Fungi_Broad_RNA_species[,colnames(df_psRep200_HiSeq_Fungi_Broad_RNA_species) %in% sharedSpeciesEukDetectFiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_Broad_RNA_species_EukDetectUnfiltered <- df_psRep200_HiSeq_Fungi_Broad_RNA_species[,colnames(df_psRep200_HiSeq_Fungi_Broad_RNA_species) %in% sharedSpeciesEukDetectUnfiltered$intersectedTaxa]
df_psRep200_HiSeq_Fungi_Broad_RNA_species_HMP <- df_psRep200_HiSeq_Fungi_Broad_RNA_species[,colnames(df_psRep200_HiSeq_Fungi_Broad_RNA_species) %in% decontamResultsV2_HMP$species]
## Broad_RNA -- remove zero sum samples
df_psRep200_HiSeq_Fungi_Broad_RNA_species_EukDetectFiltered_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_RNA_species_EukDetectFiltered[!rowSums(df_psRep200_HiSeq_Fungi_Broad_RNA_species_EukDetectFiltered)==0,]
df_psRep200_HiSeq_Fungi_Broad_RNA_species_EukDetectUnfiltered_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_RNA_species_EukDetectUnfiltered[!rowSums(df_psRep200_HiSeq_Fungi_Broad_RNA_species_EukDetectUnfiltered)==0,]
df_psRep200_HiSeq_Fungi_Broad_RNA_species_HMP_Nonzero <- df_psRep200_HiSeq_Fungi_Broad_RNA_species_HMP[!rowSums(df_psRep200_HiSeq_Fungi_Broad_RNA_species_HMP)==0,]
## Broad_RNA - metadata
metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_species_EukDetectFiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Broad_RNA_species_EukDetectFiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_species_EukDetectUnfiltered_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Broad_RNA_species_EukDetectUnfiltered_Nonzero),])
metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_species_HMP_Nonzero <- droplevels(metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Broad_RNA_species_HMP_Nonzero),])

#--------------------Save data for ML--------------------#

save(# VSNM
     rep200FungiDecontamV2SpeciesVSNM_EukDetectFiltered,
     rep200FungiDecontamV2SpeciesVSNM_EukDetectUnfiltered,
     rep200FungiDecontamV2SpeciesVSNM_HMP,
     # HMS
     df_psRep200_HiSeq_Fungi_HMS_species_EukDetectFiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_HMS_species_EukDetectUnfiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_HMS_species_HMP_Nonzero,
     # BCM
     df_psRep200_HiSeq_Fungi_BCM_species_EukDetectFiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_BCM_species_EukDetectUnfiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_BCM_species_HMP_Nonzero,
     # MDA
     df_psRep200_HiSeq_Fungi_MDA_species_EukDetectFiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_MDA_species_EukDetectUnfiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_MDA_species_HMP_Nonzero,
     # WashU
     df_psRep200_HiSeq_Fungi_WashU_species_EukDetectFiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_WashU_species_EukDetectUnfiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_WashU_species_HMP_Nonzero,
     # Broad_WGS
     df_psRep200_HiSeq_Fungi_Broad_WGS_species_EukDetectFiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_Broad_WGS_species_EukDetectUnfiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_Broad_WGS_species_HMP_Nonzero,
     # UNC
     df_psRep200_HiSeq_Fungi_UNC_species_EukDetectFiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_UNC_species_EukDetectUnfiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_UNC_species_HMP_Nonzero,
     # CMS
     df_psRep200_HiSeq_Fungi_CMS_species_EukDetectFiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_CMS_species_EukDetectUnfiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_CMS_species_HMP_Nonzero,
     # Broad_RNA
     df_psRep200_HiSeq_Fungi_Broad_RNA_species_EukDetectFiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_Broad_RNA_species_EukDetectUnfiltered_Nonzero,
     df_psRep200_HiSeq_Fungi_Broad_RNA_species_HMP_Nonzero,
     # Metadata - VSNM
     metaQiitaCombined_Nonzero_DecontamV2,
     # Metadata - HMS
     metaQiitaCombined_Nonzero_DecontamV2_HMS_species_EukDetectFiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_HMS_species_EukDetectUnfiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_HMS_species_HMP_Nonzero,
     # Metadata - BCM
     metaQiitaCombined_Nonzero_DecontamV2_BCM_species_EukDetectFiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_BCM_species_EukDetectUnfiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_BCM_species_HMP_Nonzero,
     # Metadata - MDA
     metaQiitaCombined_Nonzero_DecontamV2_MDA_species_EukDetectFiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_MDA_species_EukDetectUnfiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_MDA_species_HMP_Nonzero,
     # Metadata - WashU
     metaQiitaCombined_Nonzero_DecontamV2_WashU_species_EukDetectFiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_WashU_species_EukDetectUnfiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_WashU_species_HMP_Nonzero,
     # Metadata - Broad_WGS
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_species_EukDetectFiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_species_EukDetectUnfiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS_species_HMP_Nonzero,
     # Metadata - UNC
     metaQiitaCombined_Nonzero_DecontamV2_UNC_species_EukDetectFiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_UNC_species_EukDetectUnfiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_UNC_species_HMP_Nonzero,
     # Metadata - CMS
     metaQiitaCombined_Nonzero_DecontamV2_CMS_species_EukDetectFiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_CMS_species_EukDetectUnfiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_CMS_species_HMP_Nonzero,
     # Metadata - Broad_RNA
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_species_EukDetectFiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_species_EukDetectUnfiltered_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2_Broad_RNA_species_HMP_Nonzero,
     file = "Interim_data/data_for_ml_tcga_eukdetect_and_hmp_by_seq_center_and_taxa_levels_2Apr22.RData")

# Scripts: S17R

#----------------------------------------------------------#
# Plot EukDetect and HMP machine learning performances for all TCGA cancers
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_ED_HMP_VSNM <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_eukdetect_hmp_seq_center_taxa_level_ALL_2Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_ED_HMP_VSNM$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_ED_HMP_VSNM$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_ED_HMP_VSNM <- mlPerfAll10k_Allcancer_ED_HMP_VSNM[,!(colnames(mlPerfAll10k_Allcancer_ED_HMP_VSNM) == "X")]
colnames(mlPerfAll10k_Allcancer_ED_HMP_VSNM)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_Allcancer_ED_HMP_VSNM$nullAUPR <- ifelse(mlPerfAll10k_Allcancer_ED_HMP_VSNM$minorityClassName == "SolidTissueNormal",
                                                    yes=mlPerfAll10k_Allcancer_ED_HMP_VSNM$majorityClassSize/(mlPerfAll10k_Allcancer_ED_HMP_VSNM$minorityClassSize+mlPerfAll10k_Allcancer_ED_HMP_VSNM$majorityClassSize),
                                                    no=mlPerfAll10k_Allcancer_ED_HMP_VSNM$minorityClassSize/(mlPerfAll10k_Allcancer_ED_HMP_VSNM$minorityClassSize+mlPerfAll10k_Allcancer_ED_HMP_VSNM$majorityClassSize))
mlPerfAll10k_Allcancer_ED_HMP_VSNM$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
# VSNM
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "rep200FungiDecontamV2SpeciesVSNM_EukDetectFiltered"] <- "Species ∩ EukDetect"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "rep200FungiDecontamV2SpeciesVSNM_EukDetectUnfiltered"] <- "Species ∩ EukDetect Unfiltered"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "rep200FungiDecontamV2SpeciesVSNM_HMP"] <- "Species ∩ HMP"
# HMS
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_HMS_species_EukDetectFiltered_Nonzero"] <- "HMS species ∩ EukDetect (WGS)"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_HMS_species_EukDetectUnfiltered_Nonzero"] <- "HMS species ∩ EukDetect Unfiltered (WGS)"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_HMS_species_HMP_Nonzero"] <- "HMS species ∩ HMP (WGS)"
# BCM
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_BCM_species_EukDetectFiltered_Nonzero"] <- "BCM species ∩ EukDetect (WGS)"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_BCM_species_EukDetectUnfiltered_Nonzero"] <- "BCM species ∩ EukDetect Unfiltered (WGS)"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_BCM_species_HMP_Nonzero"] <- "BCM species ∩ HMP (WGS)"
# MDA
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_MDA_species_EukDetectFiltered_Nonzero"] <- "MDA species ∩ EukDetect (WGS)"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_MDA_species_EukDetectUnfiltered_Nonzero"] <- "MDA species ∩ EukDetect Unfiltered (WGS)"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_MDA_species_HMP_Nonzero"] <- "MDA species ∩ HMP (WGS)"

# WashU
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_WashU_species_EukDetectFiltered_Nonzero"] <- "WashU species ∩ EukDetect (WGS)"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_WashU_species_EukDetectUnfiltered_Nonzero"] <- "WashU species ∩ EukDetect Unfiltered (WGS)"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_WashU_species_HMP_Nonzero"] <- "WashU species ∩ HMP (WGS)"

# Broad_WGS
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_Broad_WGS_species_EukDetectFiltered_Nonzero"] <- "Broad species ∩ EukDetect (WGS)"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_Broad_WGS_species_EukDetectUnfiltered_Nonzero"] <- "Broad species ∩ EukDetect Unfiltered (WGS)"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_Broad_WGS_species_HMP_Nonzero"] <- "Broad species ∩ HMP (WGS)"

# UNC
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_UNC_species_EukDetectFiltered_Nonzero"] <- "UNC species ∩ EukDetect (RNA-Seq)"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_UNC_species_EukDetectUnfiltered_Nonzero"] <- "UNC species ∩ EukDetect Unfiltered (RNA-Seq)"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_UNC_species_HMP_Nonzero"] <- "UNC species ∩ HMP (RNA-Seq)"

# CMS
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_CMS_species_EukDetectFiltered_Nonzero"] <- "CMS species ∩ EukDetect (RNA-Seq)"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_CMS_species_EukDetectUnfiltered_Nonzero"] <- "CMS species ∩ EukDetect Unfiltered (RNA-Seq)"
mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName[mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName == "df_psRep200_HiSeq_Fungi_CMS_species_HMP_Nonzero"] <- "CMS species ∩ HMP (RNA-Seq)"

## NOTE: Broad_RNA only included GBM tumors, so no ML comparisons were made

mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName <- factor(mlPerfAll10k_Allcancer_ED_HMP_VSNM$datasetName,
                                                                          levels = c(#VSNM
                                                                                     "Species ∩ EukDetect",
                                                                                     "Species ∩ EukDetect Unfiltered",
                                                                                     "Species ∩ HMP",
                                                                                     # HMS
                                                                                     "HMS species ∩ EukDetect (WGS)",
                                                                                     "HMS species ∩ EukDetect Unfiltered (WGS)",
                                                                                     "HMS species ∩ HMP (WGS)",
                                                                                     # BCM
                                                                                     "BCM species ∩ EukDetect (WGS)",
                                                                                     "BCM species ∩ EukDetect Unfiltered (WGS)",
                                                                                     "BCM species ∩ HMP (WGS)",
                                                                                     # MDA
                                                                                     "MDA species ∩ EukDetect (WGS)",
                                                                                     "MDA species ∩ EukDetect Unfiltered (WGS)",
                                                                                     "MDA species ∩ HMP (WGS)",
                                                                                     # WashU
                                                                                     "WashU species ∩ EukDetect (WGS)",
                                                                                     "WashU species ∩ EukDetect Unfiltered (WGS)",
                                                                                     "WashU species ∩ HMP (WGS)",
                                                                                     # Broad_WGS
                                                                                     "Broad species ∩ EukDetect (WGS)",
                                                                                     "Broad species ∩ EukDetect Unfiltered (WGS)",
                                                                                     "Broad species ∩ HMP (WGS)",
                                                                                     # UNC
                                                                                     "UNC species ∩ EukDetect (RNA-Seq)",
                                                                                     "UNC species ∩ EukDetect Unfiltered (RNA-Seq)",
                                                                                     "UNC species ∩ HMP (RNA-Seq)",
                                                                                     # CMS
                                                                                     "CMS species ∩ EukDetect (RNA-Seq)",
                                                                                     "CMS species ∩ EukDetect Unfiltered (RNA-Seq)",
                                                                                     "CMS species ∩ HMP (RNA-Seq)"))

save(mlPerfAll10k_Allcancer_ED_HMP_VSNM,
     file = "Interim_data/mlPerfAll10k_Allcancer_ED_HMP_VSNM_2Apr22.RData")

# #-------------------------Plot ML performances-------------------------#
#----------------------PT----------------------#
# EukDetect & HMP
mlPerfAll10k_Allcancer_ED_HMP_VSNM %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | EukDetect & HMP Fungi | Primary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_EukDetect_HMP.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# EukDetect only
mlPerfAll10k_Allcancer_ED_HMP_VSNM %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species ∩ EukDetect",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | EukDetect | Primary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_EukDetect_only.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

#----------------------Primary Tumor vs NAT----------------------#
# EukDetect & HMP
mlPerfAll10k_Allcancer_ED_HMP_VSNM %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("EukDetect & HMP Fungi | Primary Tumor vs NAT") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_NAT_EukDetect_HMP.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# EukDetect only
mlPerfAll10k_Allcancer_ED_HMP_VSNM %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species ∩ EukDetect",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("EukDetect | Primary Tumor vs NAT") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_NAT_EukDetect_only.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

#----------------------BDN----------------------#
# EukDetect & HMP
mlPerfAll10k_Allcancer_ED_HMP_VSNM %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | EukDetect & HMP Fungi | Blood Derived Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_BDN_EukDetect_HMP.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# EukDetect only
mlPerfAll10k_Allcancer_ED_HMP_VSNM %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species ∩ EukDetect",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("EukDetect | Blood Derived Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_BDN_EukDetect_only.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

require(gmodels)
mlPerfAll10k_Allcancer_ED_HMP_VSNM %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species ∩ EukDetect",datasetName)) %>%
  distinct() %>% droplevels() %>%
  pull(AUROC) %>% ci()
# Estimate    CI lower    CI upper  Std. Error 
# 0.836275118 0.828741944 0.843808292 0.003836247

mlPerfAll10k_Allcancer_ED_HMP_VSNM %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species ∩ EukDetect",datasetName)) %>%
  distinct() %>% droplevels() %>%
  pull(AUROC) %>% ci()
# Estimate    CI lower    CI upper  Std. Error 
# 0.770848506 0.752618351 0.789078662 0.009263622

mlPerfAll10k_Allcancer_ED_HMP_VSNM %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species ∩ EukDetect",datasetName)) %>%
  distinct() %>% droplevels() %>%
  pull(AUROC) %>% ci()
# Estimate    CI lower    CI upper  Std. Error 
# 0.906508619 0.899480660 0.913536577 0.003574882

#----------------------------------------------------------#
# Plot machine learning performances for 8 cancers
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_8cancer <- read.csv("Interim_data/rep_perfFungi_10k_rep1_8cancer_higher_taxa_and_intersected_ALL_decontamV2_2Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_8cancer <- read.csv("Supporting_data/tcga_abbreviations_8cancer.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_8cancer$abbrev <- abbreviationsTCGA_8cancer[mlPerfAll10k_8cancer$diseaseType,"abbrev"]
mlPerfAll10k_8cancer <- mlPerfAll10k_8cancer[,!(colnames(mlPerfAll10k_8cancer) == "X")]
colnames(mlPerfAll10k_8cancer)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_8cancer$nullAUPR <- ifelse(mlPerfAll10k_8cancer$minorityClassName == "SolidTissueNormal",
                                          yes=mlPerfAll10k_8cancer$majorityClassSize/(mlPerfAll10k_8cancer$minorityClassSize+mlPerfAll10k_8cancer$majorityClassSize),
                                          no=mlPerfAll10k_8cancer$minorityClassSize/(mlPerfAll10k_8cancer$minorityClassSize+mlPerfAll10k_8cancer$majorityClassSize))
mlPerfAll10k_8cancer$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
mlPerfAll10k_8cancer$datasetName[mlPerfAll10k_8cancer$datasetName == "rep200FungiDecontamV2Order_8cancer_VSNM"] <- "Order (Decontaminated)"
mlPerfAll10k_8cancer$datasetName[mlPerfAll10k_8cancer$datasetName == "rep200FungiDecontamV2Family_8cancer_VSNM"] <- "Family (Decontaminated)"
mlPerfAll10k_8cancer$datasetName[mlPerfAll10k_8cancer$datasetName == "rep200FungiDecontamV2Genus_8cancer_VSNM"] <- "Genus (Decontaminated)"
mlPerfAll10k_8cancer$datasetName[mlPerfAll10k_8cancer$datasetName == "rep200FungiDecontamV2Species_8cancer_VSNM"] <- "Species (Decontaminated)"

mlPerfAll10k_8cancer$datasetName <- factor(mlPerfAll10k_8cancer$datasetName,
                                           levels = c("Order (Decontaminated)",
                                                      "Family (Decontaminated)",
                                                      "Genus (Decontaminated)",
                                                      "Species (Decontaminated)"))
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#

# All taxa levels and their intersections
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor") %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>% 
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid", position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor | 1 Vs All | Order, Family, Genus, Species Levels") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_8cancers_PT_order_family_genus_species_decontamV2.pdf", dpi = "retina",
       width = 10, height = 6, units = "in")

#-------------------------Plot primary tumor vs. adjacent tissue normal performance-------------------------#
# All taxa levels and their intersections
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>% mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid", position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs. Adjacent Normal | Order, Family, Genus, Species Levels") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_8cancers_PT_vs_STN_order_family_genus_species_decontamV2.pdf", dpi = "retina",
       width = 10, height = 5, units = "in")

#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
# All taxa levels and their intersections
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>% mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position=position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + 
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Blood Derived Normal | 1 Vs All | Order, Family, Genus, Species Levels") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_8cancers_BDN_order_family_genus_species_decontamV2.pdf", dpi = "retina",
       width = 10, height = 6, units = "in")

#----------------------------------------------------------#
# Plot machine learning performances for 8 cancers - WGS vs RNA
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_8cancer_WGS_RNA <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_8cancers_wgs_vs_rna_taxa_levels_ALL_decontamV2_2Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_8cancer <- read.csv("Supporting_data/tcga_abbreviations_8cancer.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_8cancer_WGS_RNA$abbrev <- abbreviationsTCGA_8cancer[mlPerfAll10k_8cancer_WGS_RNA$diseaseType,"abbrev"]
mlPerfAll10k_8cancer_WGS_RNA <- mlPerfAll10k_8cancer_WGS_RNA[,!(colnames(mlPerfAll10k_8cancer_WGS_RNA) == "X")]
colnames(mlPerfAll10k_8cancer_WGS_RNA)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_8cancer_WGS_RNA$nullAUPR <- ifelse(mlPerfAll10k_8cancer_WGS_RNA$minorityClassName == "SolidTissueNormal",
                                        yes=mlPerfAll10k_8cancer_WGS_RNA$majorityClassSize/(mlPerfAll10k_8cancer_WGS_RNA$minorityClassSize+mlPerfAll10k_8cancer_WGS_RNA$majorityClassSize),
                                        no=mlPerfAll10k_8cancer_WGS_RNA$minorityClassSize/(mlPerfAll10k_8cancer_WGS_RNA$minorityClassSize+mlPerfAll10k_8cancer_WGS_RNA$majorityClassSize))
mlPerfAll10k_8cancer_WGS_RNA$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
table(mlPerfAll10k_8cancer_WGS_RNA$datasetName)
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiDecontamV2Family_8cancer_VSNM_WGS"] <- "Family WGS (Decontaminated)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiDecontamV2Family_8cancer_VSNM_RNA"] <- "Family RNA (Decontaminated)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiDecontamV2Genus_8cancer_VSNM_WGS"] <- "Genus WGS (Decontaminated)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiDecontamV2Genus_8cancer_VSNM_RNA"] <- "Genus RNA (Decontaminated)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiDecontamV2Species_8cancer_VSNM_WGS"] <- "Species WGS (Decontaminated)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiDecontamV2Species_8cancer_VSNM_RNA"] <- "Species RNA (Decontaminated)"

mlPerfAll10k_8cancer_WGS_RNA$datasetName <- factor(mlPerfAll10k_8cancer_WGS_RNA$datasetName,
                                           levels = c("Family WGS (Decontaminated)",
                                                      "Family RNA (Decontaminated)",
                                                      "Genus WGS (Decontaminated)",
                                                      "Genus RNA (Decontaminated)",
                                                      "Species WGS (Decontaminated)",
                                                      "Species RNA (Decontaminated)"))
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#

mlPerfAll10k_8cancer_WGS_RNA %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","metadataName")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor | 1 Vs All | Species | WGS vs. RNA-Seq") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_8cancers_PT_species_wgs_vs_rna_decontamV2.pdf", dpi = "retina",
         width = 10, height = 6, units = "in")

#-------------------------Plot primary tumor vs. normal performance-------------------------#

mlPerfAll10k_8cancer_WGS_RNA %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","metadataName")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs. NAT | Species | WGS vs. RNA-Seq") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_8cancers_PT_vs_NAT_species_wgs_vs_rna_decontamV2.pdf", dpi = "retina",
       width = 10, height = 6, units = "in")

#----------------------------------------------------------#
# Plot machine learning performances for 8 cancers - shared features
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_8cancer_Shared <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_8cancers_shared_features_ALL_decontamV2_2Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_8cancer <- read.csv("Supporting_data/tcga_abbreviations_8cancer.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_8cancer_Shared$abbrev <- abbreviationsTCGA_8cancer[mlPerfAll10k_8cancer_Shared$diseaseType,"abbrev"]
mlPerfAll10k_8cancer_Shared <- mlPerfAll10k_8cancer_Shared[,!(colnames(mlPerfAll10k_8cancer_Shared) == "X")]
colnames(mlPerfAll10k_8cancer_Shared)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_8cancer_Shared$nullAUPR <- ifelse(mlPerfAll10k_8cancer_Shared$minorityClassName == "SolidTissueNormal",
                                                yes=mlPerfAll10k_8cancer_Shared$majorityClassSize/(mlPerfAll10k_8cancer_Shared$minorityClassSize+mlPerfAll10k_8cancer_Shared$majorityClassSize),
                                                no=mlPerfAll10k_8cancer_Shared$minorityClassSize/(mlPerfAll10k_8cancer_Shared$minorityClassSize+mlPerfAll10k_8cancer_Shared$majorityClassSize))
mlPerfAll10k_8cancer_Shared$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
mlPerfAll10k_8cancer_Shared$datasetName[mlPerfAll10k_8cancer_Shared$datasetName == "rep200FungiPhylumShared_8cancer_Nonzero_VSNM"] <- "Phylum ∩ Weizmann"
mlPerfAll10k_8cancer_Shared$datasetName[mlPerfAll10k_8cancer_Shared$datasetName == "rep200FungiFamilyShared_8cancer_Nonzero_VSNM"] <- "Family ∩ Weizmann"
mlPerfAll10k_8cancer_Shared$datasetName[mlPerfAll10k_8cancer_Shared$datasetName == "rep200FungiGenusShared_8cancer_Nonzero_VSNM"] <- "Genus ∩ Weizmann"
mlPerfAll10k_8cancer_Shared$datasetName[mlPerfAll10k_8cancer_Shared$datasetName == "rep200FungiSpeciesShared_8cancer_Nonzero_VSNM"] <- "Species ∩ Weizmann"

mlPerfAll10k_8cancer_Shared$datasetName <- factor(mlPerfAll10k_8cancer_Shared$datasetName,
                                                   levels = c("Phylum ∩ Weizmann",
                                                              "Family ∩ Weizmann",
                                                              "Genus ∩ Weizmann",
                                                              "Species ∩ Weizmann"))
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#

mlPerfAll10k_8cancer_Shared %>%
  filter(sampleType == "Primary Tumor") %>%
  # filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","metadataName")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor | 1 Vs All | Species ∩ Weizmann") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_8cancers_PT_species_shared_decontamV2.svg", dpi = "retina",
       width = 10, height = 6, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#

mlPerfAll10k_8cancer_Shared %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  # filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","metadataName")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs NAT | Phylum, Family, Genus, & Species ∩ Weizmann") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_8cancers_PT_vs_NAT_phylum_family_genus_species_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 6, units = "in")

#-------------------------Plot blood derived normal 1 vs all others performance-------------------------#

mlPerfAll10k_8cancer_Shared %>%
  filter(sampleType == "Blood Derived Normal") %>%
  # filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","metadataName")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Blood Derived Normal | 1 vs All Others | Phylum, Family, Genus, & Species ∩ Weizmann") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_8cancers_BDN_phylum_family_genus_species_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 6, units = "in")

#----------------------------------------------------------#
# Overlay plot machine learning performances for for 8 cancers:
# - Full features
# - WGS vs RNA
# - Shared features
#----------------------------------------------------------#

#--------------------------First load high coverage data--------------------------#
mlPerfAll10k_8cancer_Cov <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_8cancers_high_coverage_ALL_decontamV2_2Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_8cancer <- read.csv("Supporting_data/tcga_abbreviations_8cancer.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_8cancer_Cov$abbrev <- abbreviationsTCGA_8cancer[mlPerfAll10k_8cancer_Cov$diseaseType,"abbrev"]
mlPerfAll10k_8cancer_Cov <- mlPerfAll10k_8cancer_Cov[,!(colnames(mlPerfAll10k_8cancer_Cov) == "X")]
colnames(mlPerfAll10k_8cancer_Cov)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_8cancer_Cov$nullAUPR <- ifelse(mlPerfAll10k_8cancer_Cov$minorityClassName == "SolidTissueNormal",
                                              yes=mlPerfAll10k_8cancer_Cov$majorityClassSize/(mlPerfAll10k_8cancer_Cov$minorityClassSize+mlPerfAll10k_8cancer_Cov$majorityClassSize),
                                              no=mlPerfAll10k_8cancer_Cov$minorityClassSize/(mlPerfAll10k_8cancer_Cov$minorityClassSize+mlPerfAll10k_8cancer_Cov$majorityClassSize))
mlPerfAll10k_8cancer_Cov$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
mlPerfAll10k_8cancer_Cov$datasetName[mlPerfAll10k_8cancer_Cov$datasetName == "rep200FungiSpecies_8cancer_Cov_Nonzero_VSNM"] <- "Species high coverage"
mlPerfAll10k_8cancer_Cov$datasetName <- factor(mlPerfAll10k_8cancer_Cov$datasetName,
                                                 levels = c("Species high coverage"))

# Make sure to run the above 3 sections above to correctly format the following objects:
# mlPerfAll10k_8cancer, mlPerfAll10k_8cancer_WGS_RNA, mlPerfAll10k_8cancer_Shared
# Note that the mlPerfAll10k_8cancer is missing the "metadataName" column, so one is added prior to rbind
mlPerfAll10k_8cancer_Overlay <- rbind(cbind(mlPerfAll10k_8cancer, metadataName="metaQiitaCombined_Nonzero_DecontamV2_8cancer"),
                                      mlPerfAll10k_8cancer_WGS_RNA,
                                      mlPerfAll10k_8cancer_Cov,
                                      mlPerfAll10k_8cancer_Shared)

#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#

mlPerfAll10k_8cancer_Overlay %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Primary Tumor | 1 Vs All | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/overlay_mlPerfAll10k_rep1_8cancers_PT_species_decontamV2.svg", dpi = "retina",
       width = 10, height = 6, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#

mlPerfAll10k_8cancer_Overlay %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","metadataName")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Primary Tumor vs NAT | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/overlay_mlPerfAll10k_rep1_8cancers_PT_vs_NAT_species_decontamV2.svg", dpi = "retina",
         width = 10, height = 6, units = "in")

#-------------------------Plot BDN 1 vs. all others performance-------------------------#

mlPerfAll10k_8cancer_Overlay %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species [WGS|∩|high]",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","metadataName")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  # guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Blood Derived Normal | 1 Vs All | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/overlay_mlPerfAll10k_rep1_8cancers_BDN_species_decontamV2.svg", dpi = "retina",
         width = 10, height = 6, units = "in")

#----------------------------------------------------------#
# Plot machine learning performances for all TCGA cancers
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_higher_taxa_and_intersected_ALL_decontamV2_2Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer <- mlPerfAll10k_Allcancer[,!(colnames(mlPerfAll10k_Allcancer) == "X")]
colnames(mlPerfAll10k_Allcancer)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_Allcancer$nullAUPR <- ifelse(mlPerfAll10k_Allcancer$minorityClassName == "SolidTissueNormal",
                                          yes=mlPerfAll10k_Allcancer$majorityClassSize/(mlPerfAll10k_Allcancer$minorityClassSize+mlPerfAll10k_Allcancer$majorityClassSize),
                                          no=mlPerfAll10k_Allcancer$minorityClassSize/(mlPerfAll10k_Allcancer$minorityClassSize+mlPerfAll10k_Allcancer$majorityClassSize))
mlPerfAll10k_Allcancer$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
mlPerfAll10k_Allcancer$datasetName[mlPerfAll10k_Allcancer$datasetName == "rep200FungiDecontamV2OrderVSNM"] <- "Order (Decontaminated)"
mlPerfAll10k_Allcancer$datasetName[mlPerfAll10k_Allcancer$datasetName == "rep200FungiDecontamV2FamilyVSNM"] <- "Family (Decontaminated)"
mlPerfAll10k_Allcancer$datasetName[mlPerfAll10k_Allcancer$datasetName == "rep200FungiDecontamV2GenusVSNM"] <- "Genus (Decontaminated)"
mlPerfAll10k_Allcancer$datasetName[mlPerfAll10k_Allcancer$datasetName == "rep200FungiDecontamV2SpeciesVSNM"] <- "Species (Decontaminated)"

mlPerfAll10k_Allcancer$datasetName <- factor(mlPerfAll10k_Allcancer$datasetName,
                                           levels = c("Order (Decontaminated)",
                                                      "Family (Decontaminated)",
                                                      "Genus (Decontaminated)",
                                                      "Species (Decontaminated)"))

save(mlPerfAll10k_Allcancer, file = "Interim_data/mlPerfAll10k_Allcancer_2Apr22.RData")

#-------------------------Regress avg performance vs minority class size-------------------------#
ptPerfSummarizedMinClassSize <- mlPerfAll10k_Allcancer %>% distinct() %>% filter(grepl("Species",datasetName)) %>% 
                            filter(grepl("Decontaminated",datasetName)) %>%
                            filter(sampleType == "Primary Tumor") %>%
                            group_by(sampleType, diseaseType) %>% 
                            mutate(avgROC = mean(AUROC), avgPR = mean(AUPR)) %>%
                            select(sampleType, diseaseType, abbrev, avgROC, avgPR, minorityClassSize, majorityClassSize) %>% 
                            mutate(logAvgMinClass = log10(minorityClassSize), ratioMinMajClass=minorityClassSize/majorityClassSize) %>%
                            unique() %>% data.frame()

## Regress using minority class size
summary(lm(avgPR ~ minorityClassSize, ptPerfSummarizedMinClassSize))
summary(lm(avgROC ~ minorityClassSize, ptPerfSummarizedMinClassSize))

require(ggrepel)
ptPerfSummarizedMinClassSize %>%
  ggplot(aes(x=minorityClassSize, y=avgPR, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.x = 0, label.y = 0.6) + 
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  geom_text_repel() +
  labs(x = "Minority class size (per cancer type)", y = "Average AUPR", 
       title = "Correlating AUPR and minority class size\nfor primary tumor TCGA models")
ggsave(filename = "Figures/Supplementary_Figures/aupr_vs_minority_class_size_all_TCGA_PT.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

ptPerfSummarizedMinClassSize %>%
  ggplot(aes(x=minorityClassSize, y=avgROC, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.x = 600, label.y = 0.7) + 
  geom_text_repel() +
  theme_bw()+ theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  labs(x = "Minority class size (per cancer type)", y = "Average AUROC", 
       title = "Correlating AUROC and minority class size\nfor primary tumor TCGA models")
ggsave(filename = "Figures/Supplementary_Figures/auroc_vs_minority_class_size_all_TCGA_PT.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

## Regress using minority class *ratio*
summary(lm(avgPR ~ ratioMinMajClass, ptPerfSummarizedMinClassSize))
summary(lm(avgROC ~ ratioMinMajClass, ptPerfSummarizedMinClassSize))

ptPerfSummarizedMinClassSize %>%
  ggplot(aes(x=ratioMinMajClass, y=avgPR, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.y = 1.1) + 
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  geom_text_repel() +
  labs(x = "Ratio of minority to majority class size (per cancer type)", y = "Average AUPR", 
       title = "Correlating AUPR and ratio of minority to majority class size\nfor primary tumor TCGA models")
ggsave(filename = "Figures/Supplementary_Figures/aupr_vs_ratio_minority_to_majority_class_size_all_TCGA_PT.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

ptPerfSummarizedMinClassSize %>%
  ggplot(aes(x=ratioMinMajClass, y=avgROC, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.y = 1.05) + 
  geom_text_repel() +
  theme_bw()+ theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  labs(x = "Ratio of minority to majority class size (per cancer type)", y = "Average AUROC", 
       title = "Correlating AUROC and ratio of minority to majority class size\nfor primary tumor TCGA models")
ggsave(filename = "Figures/Supplementary_Figures/auroc_vs_ratio_minority_to_majority_class_size_all_TCGA_PT.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

#-------------------------Regress avg performance vs average read depth-------------------------#
load("Interim_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_25Mar22.RData",
     verbose = TRUE)

metaAvgReadDepthVSNM <- droplevels(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts[rownames(metaQiitaCombined_Nonzero_DecontamV2),])
metaAvgReadDepthVSNM %>% 
  filter(sample_type == "Primary Tumor") %>%
  mutate(abbrev = gsub("^TCGA-","",investigation)) %>% group_by(abbrev) %>%
  summarize(Mean_bam_total_reads = mean(bam_total_reads), 
            Mean_bam_mapped_reads = mean(bam_mapped_reads),
            Mean_bam_unmapped_reads = mean(bam_unmapped_reads),
            Mean_bam_ratio_unmapped = mean(bam_ratio_unmapped)) %>% data.frame() -> metaAvgReadDepthVSNMFormatted

ptPerfSummarizedReadDepth <- mlPerfAll10k_Allcancer %>% distinct() %>% filter(grepl("Species",datasetName)) %>% 
  filter(grepl("Decontaminated",datasetName)) %>%
  filter(sampleType == "Primary Tumor") %>%
  group_by(sampleType, diseaseType) %>% 
  mutate(avgROC = mean(AUROC), avgPR = mean(AUPR)) %>%
  select(sampleType, diseaseType, abbrev, avgROC, avgPR, minorityClassSize, majorityClassSize) %>% 
  mutate(logAvgMinClass = log10(minorityClassSize), ratioMinMajClass=minorityClassSize/majorityClassSize) %>%
  unique() %>% 
  left_join(metaAvgReadDepthVSNMFormatted, by = "abbrev") %>% 
  data.frame()

## Mean total bam reads
summary(lm(avgPR ~ Mean_bam_total_reads, ptPerfSummarizedReadDepth))
summary(lm(avgROC ~ Mean_bam_total_reads, ptPerfSummarizedReadDepth))

require(ggrepel)
ptPerfSummarizedReadDepth %>%
  ggplot(aes(x=Mean_bam_total_reads, y=avgPR, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.x = 1e9, label.y = 1) + 
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  geom_text_repel() +
  labs(x = "Average Total Read Depth Per\nCancer Type (Primary Tumor)", y = "Average AUPR", 
       title = "Correlating AUPR and read depth\nfor primary tumor TCGA models")
ggsave(filename = "Figures/Supplementary_Figures/aupr_vs_read_depth_all_TCGA_PT.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

ptPerfSummarizedReadDepth %>%
  ggplot(aes(x=Mean_bam_total_reads, y=avgROC, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.x = 1e9, label.y = 1.05) + 
  geom_text_repel() +
  theme_bw()+ theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  labs(x = "Average Total Read Depth Per\nCancer Type (Primary Tumor)", y = "Average AUROC", 
       title = "Correlating AUROC and read depth\nfor primary tumor TCGA models")
ggsave(filename = "Figures/Supplementary_Figures/auroc_vs_read_depth_all_TCGA_PT.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

## Mean unmapped bam reads
summary(lm(avgPR ~ Mean_bam_unmapped_reads, ptPerfSummarizedReadDepth))
summary(lm(avgROC ~ Mean_bam_unmapped_reads, ptPerfSummarizedReadDepth))

ptPerfSummarizedReadDepth %>%
  ggplot(aes(x=Mean_bam_unmapped_reads, y=avgPR, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.y = 0.82) + 
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  geom_text_repel() +
  labs(x = "Average Unmapped Reads Per\nCancer Type (Primary Tumor)", y = "Average AUPR", 
       title = "Correlating AUPR and unmapped reads\nfor primary tumor TCGA models")
ggsave(filename = "Figures/Supplementary_Figures/aupr_vs_unmapped_reads_all_TCGA_PT.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

ptPerfSummarizedReadDepth %>%
  ggplot(aes(x=Mean_bam_unmapped_reads, y=avgROC, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red",label.y = 1.05) + 
  geom_text_repel() +
  theme_bw()+ theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  labs(x = "Average Unmapped Reads Per\nCancer Type (Primary Tumor)", y = "Average AUROC", 
       title = "Correlating AUROC and unmapped reads\nfor primary tumor TCGA models")
ggsave(filename = "Figures/Supplementary_Figures/auroc_vs_unmapped_reads_all_TCGA_PT.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
# Order level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Order",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","minorityClassSize","majorityClassSize","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor | 1 Vs All | Order") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_order_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

# Family level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Family",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","minorityClassSize","majorityClassSize","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor | 1 Vs All | Family") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_family_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

# Genus level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Genus",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","minorityClassSize","majorityClassSize","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor | 1 Vs All | Genus") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_genus_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

# Species level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","minorityClassSize","majorityClassSize","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_species_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

require(gmodels)
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  pull(AUROC) %>% ci()
  
#-------------------------Plot primary tumor vs. adjacent normal tissue performance-------------------------#
# Order level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Order",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","minorityClassSize","majorityClassSize","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor vs. Adjacent Normal | Order") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_STN_order_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

# Family level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Family",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","minorityClassSize","majorityClassSize","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor vs. Adjacent Normal | Family") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_STN_family_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

# Genus level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Genus",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","minorityClassSize","majorityClassSize","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor vs. Adjacent Normal | Genus") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_STN_genus_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

# Species level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","minorityClassSize","majorityClassSize","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor vs. Adjacent Normal | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_STN_species_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
# Order level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Order",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","minorityClassSize","majorityClassSize","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Blood Derived Normal | 1 Vs All | Order") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_BDN_order_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

# Family level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Family",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","minorityClassSize","majorityClassSize","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Blood Derived Normal | 1 Vs All | Family") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_BDN_family_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

# Genus level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Genus",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","minorityClassSize","majorityClassSize","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Blood Derived Normal | 1 Vs All | Genus") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_BDN_genus_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

# Species level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","minorityClassSize","majorityClassSize","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_BDN_species_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

require(gmodels)
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  pull(AUROC) %>% ci()

#----------------------------------------------------------#
# Plot machine learning performances for all TCGA cancers - shared features
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Shared <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_shared_features_ALL_decontamV2_2Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_Shared$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_Shared$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_Shared <- mlPerfAll10k_Allcancer_Shared[,!(colnames(mlPerfAll10k_Allcancer_Shared) == "X")]
colnames(mlPerfAll10k_Allcancer_Shared)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_Allcancer_Shared$nullAUPR <- ifelse(mlPerfAll10k_Allcancer_Shared$minorityClassName == "SolidTissueNormal",
                                               yes=mlPerfAll10k_Allcancer_Shared$majorityClassSize/(mlPerfAll10k_Allcancer_Shared$minorityClassSize+mlPerfAll10k_Allcancer_Shared$majorityClassSize),
                                               no=mlPerfAll10k_Allcancer_Shared$minorityClassSize/(mlPerfAll10k_Allcancer_Shared$minorityClassSize+mlPerfAll10k_Allcancer_Shared$majorityClassSize))
mlPerfAll10k_Allcancer_Shared$nullAUROC <- 0.5

# Rename entries in the "datasetName" column

mlPerfAll10k_Allcancer_Shared$datasetName[mlPerfAll10k_Allcancer_Shared$datasetName == "rep200FungiFamilyShared_NonzeroVSNM"] <- "Family ∩ Weizmann"
mlPerfAll10k_Allcancer_Shared$datasetName[mlPerfAll10k_Allcancer_Shared$datasetName == "rep200FungiGenusShared_NonzeroVSNM"] <- "Genus ∩ Weizmann"
mlPerfAll10k_Allcancer_Shared$datasetName[mlPerfAll10k_Allcancer_Shared$datasetName == "rep200FungiSpeciesShared_NonzeroVSNM"] <- "Species ∩ Weizmann"
mlPerfAll10k_Allcancer_Shared$datasetName <- factor(mlPerfAll10k_Allcancer_Shared$datasetName,
                                             levels = c("Family ∩ Weizmann",
                                                        "Genus ∩ Weizmann",
                                                        "Species ∩ Weizmann"))

save(mlPerfAll10k_Allcancer_Shared, file = "Interim_data/mlPerfAll10k_Allcancer_Shared_2Apr22.RData")
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
mlPerfAll10k_Allcancer_Shared %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor | 1 Vs All | Species ∩ Weizmann") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_species_shared_decontamV2.svg", dpi = "retina",
         width = 12, height = 6, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#
mlPerfAll10k_Allcancer_Shared %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor vs NAT | Species ∩ Weizmann") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_NAT_species_shared_decontamV2.svg", dpi = "retina",
         width = 12, height = 6, units = "in")

#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
mlPerfAll10k_Allcancer_Shared %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Blood Derived Normal | 1 Vs All | Species ∩ Weizmann") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_BDN_species_shared_decontamV2.svg", dpi = "retina",
         width = 12, height = 6, units = "in")

#----------------------------------------------------------#
# Plot machine learning performances for WGS vs RNA-Seq
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_full_WGS_RNA <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_full_wgs_rna_vsnm_ALL_DecontamV2_2Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_full_WGS_RNA$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_full_WGS_RNA$diseaseType,"abbrev"]
mlPerfAll10k_full_WGS_RNA <- mlPerfAll10k_full_WGS_RNA[,!(colnames(mlPerfAll10k_full_WGS_RNA) == "X")]
colnames(mlPerfAll10k_full_WGS_RNA)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_full_WGS_RNA$nullAUPR <- ifelse(mlPerfAll10k_full_WGS_RNA$minorityClassName == "SolidTissueNormal",
                                                 yes=mlPerfAll10k_full_WGS_RNA$majorityClassSize/(mlPerfAll10k_full_WGS_RNA$minorityClassSize+mlPerfAll10k_full_WGS_RNA$majorityClassSize),
                                                 no=mlPerfAll10k_full_WGS_RNA$minorityClassSize/(mlPerfAll10k_full_WGS_RNA$minorityClassSize+mlPerfAll10k_full_WGS_RNA$majorityClassSize))
mlPerfAll10k_full_WGS_RNA$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
mlPerfAll10k_full_WGS_RNA$datasetName[mlPerfAll10k_full_WGS_RNA$datasetName == "snmDataOGUFungiDecontamV2"] <- "WGS+RNA-Seq Species (Decontaminated)"
mlPerfAll10k_full_WGS_RNA$datasetName[mlPerfAll10k_full_WGS_RNA$datasetName == "snmDataOGUFungiDecontamV2_WGS"] <- "WGS Species (Decontaminated)"
mlPerfAll10k_full_WGS_RNA$datasetName[mlPerfAll10k_full_WGS_RNA$datasetName == "snmDataOGUFungiDecontamV2_RNA"] <- "RNA-Seq Species (Decontaminated)"

mlPerfAll10k_full_WGS_RNA$datasetName <- factor(mlPerfAll10k_full_WGS_RNA$datasetName,
                                             levels = c("WGS+RNA-Seq Species (Decontaminated)",
                                                        "WGS Species (Decontaminated)",
                                                        "RNA-Seq Species (Decontaminated)"))
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
mlPerfAll10k_full_WGS_RNA %>%
  filter(sampleType == "Primary Tumor") %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>% 
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor | 1 Vs All | Species | WGS vs. RNA") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_wgs_vs_rna_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

#-------------------------Plot primary tumor vs. adjacent normal performance-------------------------#
mlPerfAll10k_full_WGS_RNA %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>% 
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species | WGS vs. RNA") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_STN_wgs_vs_rna_decontamV2.pdf", dpi = "retina",
         width = 12, height = 6, units = "in")

## NOTE: Only WGS was run on Blood Derived Normal TCGA samples, so there is no comparison between them

#-------------------------Linear regression of WGS-only OR RNA-Seq-only vs. VSNM data-------------------------#

perfSummarizedFullvsWGSvsRNA <- mlPerfAll10k_full_WGS_RNA %>% distinct() %>%
  filter(sampleType == "Primary Tumor") %>%
  group_by(sampleType, diseaseType, datasetName) %>% 
  mutate(avgROC = mean(AUROC), stdROC=sd(AUROC), avgPR = mean(AUPR), stdPR=sd(AUPR),) %>%
  select(sampleType, diseaseType, datasetName, abbrev, avgROC, stdROC, avgPR, stdPR, minorityClassSize, majorityClassSize) %>% 
  unique() %>% data.frame()

summary(lm(avgPR ~ minorityClassSize, perfSummarizedFullvsWGSvsRNA))
summary(lm(avgROC ~ minorityClassSize, perfSummarizedFullvsWGSvsRNA))

perfSummarizedFullvsWGSvsRNA_Full <- perfSummarizedFullvsWGSvsRNA %>%
  filter(datasetName=="WGS+RNA-Seq Species (Decontaminated)") %>% droplevels()
perfSummarizedFullvsWGSvsRNA_WGSonly <- perfSummarizedFullvsWGSvsRNA %>%
  filter(datasetName=="WGS Species (Decontaminated)") %>% droplevels() %>%
  rename(avgROC_WGS = avgROC, stdROC_WGS = stdROC, avgPR_WGS = avgPR, stdPR_WGS = stdPR)
perfSummarizedFullvsWGSvsRNA_RNAonly <- perfSummarizedFullvsWGSvsRNA %>%
  filter(datasetName=="RNA-Seq Species (Decontaminated)") %>% droplevels() %>%
  rename(avgROC_RNA = avgROC, stdROC_RNA = stdROC, avgPR_RNA = avgPR, stdPR_RNA = stdPR)

perfSummarizedFullvsWGSvsRNA_Full_Joined <- perfSummarizedFullvsWGSvsRNA_Full %>%
  left_join(perfSummarizedFullvsWGSvsRNA_WGSonly[,grep("abbrev|ROC|PR",colnames(perfSummarizedFullvsWGSvsRNA_WGSonly))], by = "abbrev") %>%
  left_join(perfSummarizedFullvsWGSvsRNA_RNAonly[,grep("abbrev|ROC|PR",colnames(perfSummarizedFullvsWGSvsRNA_RNAonly))], by = "abbrev")

require(ggpmisc)
# WGS vs. Full | AUROC
ggplot(perfSummarizedFullvsWGSvsRNA_Full_Joined, 
       aes(x = avgROC, y =  avgROC_WGS, label = abbrev)) +
  geom_point(alpha = 0.8) +
  coord_equal() + xlim(c(0, 1)) + ylim(c(0,1)) +
  theme_bw() +
  geom_abline(linetype = 2) +
  labs(x = "Average classifier performances on full data\n(primary tumor 1 cancer type vs all others)",
       y = "Average classifier performances on WGS-only data\n(primary tumor 1 cancer type vs all others)",
       title = "WGS-only vs. full data | AUROC") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  geom_ribbon(stat='smooth', method = "lm", fill = "orange", se=TRUE, alpha=0.2, fullrange = TRUE) +
  geom_line(stat='smooth', method = "lm", color = "darkorange", alpha=1, fullrange = TRUE) +
  geom_text_repel(force = 120, size=3) +
  stat_poly_eq(formula = y ~ x, label.x = 0.9, label.y = 0.05,
               aes(label = paste(..eq.label.., sep = "~~~")), 
               parse = TRUE)
ggsave(filename = "Figures/Other_Figures/regr_auroc_tcga_wgs_vs_full_PT_VSNM.pdf",
         dpi = "retina", units = "in", height = 5, width = 5)
# WGS vs. Full | AUPR
ggplot(perfSummarizedFullvsWGSvsRNA_Full_Joined, 
       aes(x = avgPR, y =  avgPR_WGS, label = abbrev)) +
  geom_point(alpha = 0.8) +
  coord_equal() + xlim(c(0, 1)) + ylim(c(0,1)) +
  theme_bw() +
  geom_abline(linetype = 2) +
  labs(x = "Average classifier performances on full data\n(primary tumor 1 cancer type vs all others)",
       y = "Average classifier performances on WGS-only data\n(primary tumor 1 cancer type vs all others)",
       title = "WGS-only vs. full data | AUPR") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  geom_ribbon(stat='smooth', method = "lm", fill = "orange", se=TRUE, alpha=0.2, fullrange = TRUE) +
  geom_line(stat='smooth', method = "lm", color = "darkorange", alpha=1, fullrange = TRUE) +
  geom_text_repel(size=3) +
  stat_poly_eq(formula = y ~ x, label.x = 0.9, label.y = 0.05,
               aes(label = paste(..eq.label.., sep = "~~~")), 
               parse = TRUE)
ggsave(filename = "Figures/Other_Figures/regr_aupr_tcga_wgs_vs_full_PT_VSNM.pdf",
         dpi = "retina", units = "in", height = 5, width = 5)
# RNA vs. Full | AUROC
ggplot(perfSummarizedFullvsWGSvsRNA_Full_Joined, 
       aes(x = avgROC, y =  avgROC_RNA, label = abbrev)) +
  geom_point(alpha = 0.8) +
  coord_equal() + xlim(c(0, 1)) + ylim(c(0,1)) +
  theme_bw() +
  geom_abline(linetype = 2) +
  labs(x = "Average classifier performances on full data\n(primary tumor 1 cancer type vs all others)",
       y = "Average classifier performances on RNA-only data\n(primary tumor 1 cancer type vs all others)",
       title = "RNA-only vs. full data | AUROC") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  geom_ribbon(stat='smooth', method = "lm", fill = "orange", se=TRUE, alpha=0.2, fullrange = TRUE) +
  geom_line(stat='smooth', method = "lm", color = "darkorange", alpha=1, fullrange = TRUE) +
  geom_text_repel(force = 120, size=3) +
  stat_poly_eq(formula = y ~ x, label.x = 0.7, label.y = 0.05,
               aes(label = paste(..eq.label.., sep = "~~~")), 
               parse = TRUE)
ggsave(filename = "Figures/Other_Figures/regr_auroc_tcga_rna_vs_full_PT_VSNM.pdf",
         dpi = "retina", units = "in", height = 5, width = 5)
# RNA vs. Full | AUPR
ggplot(perfSummarizedFullvsWGSvsRNA_Full_Joined, 
       aes(x = avgPR, y =  avgPR_RNA, label = abbrev)) +
  geom_point(alpha = 0.8) +
  coord_equal() + xlim(c(0, 1)) + ylim(c(0,1)) +
  theme_bw() +
  geom_abline(linetype = 2) +
  labs(x = "Average classifier performances on full data\n(primary tumor 1 cancer type vs all others)",
       y = "Average classifier performances on RNA-only data\n(primary tumor 1 cancer type vs all others)",
       title = "RNA-only vs. full data | AUPR") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  geom_ribbon(stat='smooth', method = "lm", fill = "orange", se=TRUE, alpha=0.2, fullrange = TRUE) +
  geom_line(stat='smooth', method = "lm", color = "darkorange", alpha=1, fullrange = TRUE) +
  geom_text_repel(force = 150, size=3) +
  stat_poly_eq(formula = y ~ x, label.x = 0.05, label.y = 0.9,
               aes(label = paste(..eq.label.., sep = "~~~")), 
               parse = TRUE)
ggsave(filename = "Figures/Other_Figures/regr_aupr_tcga_rna_vs_full_PT_VSNM.pdf",
         dpi = "retina", units = "in", height = 5, width = 5)

#----------------------------------------------------------#
# Overlay plots machine learning performances for all TCGA cancers:
# - VSNM data
# - WGS vs. RNA
# - Shared features
# - High coverage features
#----------------------------------------------------------#

#--------------------------First load high coverage data--------------------------#
mlPerfAll10k_Allcancer_Cov <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_high_coverage_ALL_decontamV2_2Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_Cov$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_Cov$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_Cov <- mlPerfAll10k_Allcancer_Cov[,!(colnames(mlPerfAll10k_Allcancer_Cov) == "X")]
colnames(mlPerfAll10k_Allcancer_Cov)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_Allcancer_Cov$nullAUPR <- ifelse(mlPerfAll10k_Allcancer_Cov$minorityClassName == "SolidTissueNormal",
                                              yes=mlPerfAll10k_Allcancer_Cov$majorityClassSize/(mlPerfAll10k_Allcancer_Cov$minorityClassSize+mlPerfAll10k_Allcancer_Cov$majorityClassSize),
                                              no=mlPerfAll10k_Allcancer_Cov$minorityClassSize/(mlPerfAll10k_Allcancer_Cov$minorityClassSize+mlPerfAll10k_Allcancer_Cov$majorityClassSize))
mlPerfAll10k_Allcancer_Cov$nullAUROC <- 0.5

# Rename entries in the "datasetName" column

mlPerfAll10k_Allcancer_Cov$datasetName[mlPerfAll10k_Allcancer_Cov$datasetName == "rep200FungiCovSpeciesVSNM"] <- "Species high coverage"

mlPerfAll10k_Allcancer_Cov$datasetName <- factor(mlPerfAll10k_Allcancer_Cov$datasetName,
                                                 levels = c("Species high coverage"))

save(mlPerfAll10k_Allcancer_Cov, file = "Interim_data/mlPerfAll10k_Allcancer_Cov_2Apr22.RData")

#------------------------------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

# NOTE: There are too many cancer types to include more than two errorbars per cancer type,
# so the WGS vs. RNA comparison is not included in this overlay. If you would like to include
# it, you can remove the commented line below.
mlPerfAll10k_Allcancer_Overlay <- rbind(cbind(mlPerfAll10k_Allcancer, metadataName=NA),
                                        mlPerfAll10k_Allcancer_Cov,
                                        # cbind(mlPerfAll10k_WGS_RNA, metadataName=NA),
                                        mlPerfAll10k_Allcancer_Shared,
                                        mlPerfAll10k_Allcancer_ED_HMP_VSNM)
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
mlPerfAll10k_Allcancer_Overlay %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  mutate(datasetName = factor(datasetName, levels = c("Species high coverage", 
                                                      "Species ∩ Weizmann", 
                                                      "Species (Decontaminated)"))) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","metadataName")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  # ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  ggplot(aes(abbrev,value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.6,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All TCGA Cancers | Primary Tumor | 1 Vs All | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Main_Figures/overlay_mlPerfAll10k_rep1_Allcancers_PT_species_decontamV2_alphabetical.svg", dpi = "retina",
         width = 12, height = 4, units = "in")

mlPerfAll10k_Allcancer_Overlay %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  mutate(datasetName = factor(datasetName, levels = c("Species ∩ EukDetect",
                                                      "Species high coverage", 
                                                      "Species ∩ Weizmann", 
                                                      "Species (Decontaminated)"))) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","metadataName")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  # ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  ggplot(aes(abbrev,value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.6,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All TCGA Cancers | Primary Tumor | 1 Vs All | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Main_Figures/overlay_mlPerfAll10k_rep1_Allcancers_PT_species_decontamV2_alphabetical_EukDetect.svg", dpi = "retina",
       width = 12, height = 4, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#
mlPerfAll10k_Allcancer_Overlay %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","metadataName")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs NAT | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Supplementary_Figures/overlay_mlPerfAll10k_rep1_Allcancers_PT_vs_NAT_species_decontamV2.svg", dpi = "retina",
         width = 10, height = 6, units = "in")

mlPerfAll10k_Allcancer_Overlay %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  mutate(datasetName = factor(datasetName, levels = c("Species ∩ EukDetect",
                                                      "Species high coverage", 
                                                      "Species ∩ Weizmann", 
                                                      "Species (Decontaminated)"))) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","metadataName")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  # ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  ggplot(aes(abbrev,value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.6,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs NAT | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/overlay_mlPerfAll10k_rep1_Allcancers_PT_vs_NAT_species_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 6, units = "in")
#-------------------------Plot BDN 1 vs. all others performance-------------------------#
mlPerfAll10k_Allcancer_Overlay %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  mutate(datasetName = factor(datasetName, levels = c("Species high coverage", "Species ∩ Weizmann", "Species (Decontaminated)"))) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","metadataName")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  # ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  ggplot(aes(abbrev,value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Blood Derived Normal | 1 Vs All | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Main_Figures/overlay_mlPerfAll10k_rep1_Allcancers_BDN_species_decontamV2_alphabetical.svg", dpi = "retina",
         width = 12, height = 4, units = "in")

mlPerfAll10k_Allcancer_Overlay %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  mutate(datasetName = factor(datasetName, levels = c("Species ∩ EukDetect",
                                                      "Species high coverage", 
                                                      "Species ∩ Weizmann", 
                                                      "Species (Decontaminated)"))) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC","metadataName")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  # ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  ggplot(aes(abbrev,value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.6,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=1) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Blood Derived Normal | 1 Vs All | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Main_Figures/overlay_mlPerfAll10k_rep1_Allcancers_BDN_species_decontamV2_alphabetical_EukDetect.svg", dpi = "retina",
       width = 12, height = 4, units = "in")

#-------------------------Linear regression of EukDetect perf-------------------------#

#--------------------Primary Tumor--------------------#
perfSummarizedEukDetectPT <- mlPerfAll10k_Allcancer_Overlay %>% distinct() %>% 
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  droplevels() %>%
  group_by(sampleType, diseaseType, datasetName) %>% 
  mutate(avgROC = mean(AUROC), stdROC=sd(AUROC), avgPR = mean(AUPR), stdPR=sd(AUPR),) %>%
  select(sampleType, diseaseType, datasetName, abbrev, avgROC, stdROC, avgPR, stdPR, minorityClassSize, majorityClassSize) %>% 
  distinct() %>% data.frame()

perfSummarizedEukDetectPT_Cov <- perfSummarizedEukDetectPT %>%
  filter(datasetName=="Species high coverage") %>% droplevels()
perfSummarizedEukDetectPT_EukDetect <- perfSummarizedEukDetectPT %>%
  filter(datasetName=="Species ∩ EukDetect") %>% droplevels() %>%
  rename(avgROC_ED = avgROC, stdROC_ED = stdROC, avgPR_ED = avgPR, stdPR_ED = stdPR)

perfSummarizedEukDetectPT_Joined <- perfSummarizedEukDetectPT %>% select(abbrev) %>%
  left_join(perfSummarizedEukDetectPT_Cov, by = "abbrev") %>%
  left_join(perfSummarizedEukDetectPT_EukDetect[,grep("abbrev|ROC|PR",colnames(perfSummarizedEukDetectPT_EukDetect))], by = "abbrev")

summary(lm(avgROC ~ avgROC_ED, data=perfSummarizedEukDetectPT_Joined))
summary(lm(avgPR ~ avgPR_ED, data=perfSummarizedEukDetectPT_Joined))

require(ggpmisc)
require(ggrepel)
# AUROC
ggplot(perfSummarizedEukDetectPT_Joined, 
       aes(x = avgROC, y =  avgROC_ED)) +
  geom_point(alpha = 0.8) +
  coord_equal() + xlim(c(0, 1)) + ylim(c(0,1)) +
  theme_bw() +
  geom_abline(linetype = 2) +
  labs(x = "Average classifier performances on EukDetect fungi\n(primary tumor 1 cancer type vs all others)",
       y = "Average classifier performances on high coverage fungi\n(primary tumor 1 cancer type vs all others)",
       title = "High coverage vs. EukDetect Fungi | AUROC\nPrimary Tumor") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  geom_ribbon(stat='smooth', method = "lm", fill = "orange", se=TRUE, alpha=0.2, fullrange = TRUE) +
  geom_line(stat='smooth', method = "lm", color = "darkorange", alpha=1, fullrange = TRUE) +
  stat_poly_eq(formula = y ~ x, label.x = 0.9, label.y = 0.05,
               aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~")), 
               parse = TRUE)
ggsave(filename = "Figures/Supplementary_Figures/regr_auroc_tcga_VSNM_PT_high_cov_vs_EukDetect.pdf",
       dpi = "retina", units = "in", height = 5, width = 5)

# AUPR
ggplot(perfSummarizedEukDetectPT_Joined, 
       aes(x = avgPR, y =  avgPR_ED)) +
  geom_point(alpha = 0.8) +
  coord_equal() + xlim(c(0, 1)) + ylim(c(0,1)) +
  theme_bw() +
  geom_abline(linetype = 2) +
  labs(x = "Average classifier performances on EukDetect fungi\n(primary tumor 1 cancer type vs all others)",
       y = "Average classifier performances on high coverage fungi\n(primary tumor 1 cancer type vs all others)",
       title = "High coverage vs. EukDetect Fungi | AUPR\nPrimary Tumor") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  geom_ribbon(stat='smooth', method = "lm", fill = "orange", se=TRUE, alpha=0.2, fullrange = TRUE) +
  geom_line(stat='smooth', method = "lm", color = "darkorange", alpha=1, fullrange = TRUE) +
  stat_poly_eq(formula = y ~ x, label.x = 0.9, label.y = 0.05,
               aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~")), 
               parse = TRUE)
ggsave(filename = "Figures/Supplementary_Figures/regr_aupr_tcga_PT_high_cov_vs_EukDetect.pdf",
       dpi = "retina", units = "in", height = 5, width = 5)

#--------------------Blood Derived Normal--------------------#
perfSummarizedEukDetectBDN <- mlPerfAll10k_Allcancer_Overlay %>% distinct() %>% 
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  droplevels() %>%
  group_by(sampleType, diseaseType, datasetName) %>% 
  mutate(avgROC = mean(AUROC), stdROC=sd(AUROC), avgPR = mean(AUPR), stdPR=sd(AUPR),) %>%
  select(sampleType, diseaseType, datasetName, abbrev, avgROC, stdROC, avgPR, stdPR, minorityClassSize, majorityClassSize) %>% 
  distinct() %>% data.frame()

perfSummarizedEukDetectBDN_Cov <- perfSummarizedEukDetectBDN %>%
  filter(datasetName=="Species high coverage") %>% droplevels()
perfSummarizedEukDetectBDN_EukDetect <- perfSummarizedEukDetectBDN %>%
  filter(datasetName=="Species ∩ EukDetect") %>% droplevels() %>%
  rename(avgROC_ED = avgROC, stdROC_ED = stdROC, avgPR_ED = avgPR, stdPR_ED = stdPR)

perfSummarizedEukDetectBDN_Joined <- perfSummarizedEukDetectBDN %>% select(abbrev) %>%
  left_join(perfSummarizedEukDetectBDN_Cov, by = "abbrev") %>%
  left_join(perfSummarizedEukDetectBDN_EukDetect[,grep("abbrev|ROC|PR",colnames(perfSummarizedEukDetectBDN_EukDetect))], by = "abbrev")

summary(lm(avgROC ~ avgROC_ED, data=perfSummarizedEukDetectBDN_Joined))
summary(lm(avgPR ~ avgPR_ED, data=perfSummarizedEukDetectBDN_Joined))

require(ggpmisc)
require(ggrepel)
# AUROC
ggplot(perfSummarizedEukDetectBDN_Joined, 
       aes(x = avgROC, y =  avgROC_ED)) +
  geom_point(alpha = 0.8) +
  coord_equal() + xlim(c(0, 1)) + ylim(c(0,1)) +
  theme_bw() +
  geom_abline(linetype = 2) +
  labs(x = "Average classifier performances on EukDetect fungi\n(primary tumor 1 cancer type vs all others)",
       y = "Average classifier performances on high coverage fungi\n(primary tumor 1 cancer type vs all others)",
       title = "High coverage vs. EukDetect Fungi | AUROC\nBlood Derived Normal") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  geom_ribbon(stat='smooth', method = "lm", fill = "orange", se=TRUE, alpha=0.2, fullrange = TRUE) +
  geom_line(stat='smooth', method = "lm", color = "darkorange", alpha=1, fullrange = TRUE) +
  stat_poly_eq(formula = y ~ x, label.x = 0.9, label.y = 0.05,
               aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~")), 
               parse = TRUE)
ggsave(filename = "Figures/Supplementary_Figures/regr_auroc_tcga_VSNM_BDN_high_cov_vs_EukDetect.pdf",
       dpi = "retina", units = "in", height = 5, width = 5)

# AUPR
ggplot(perfSummarizedEukDetectBDN_Joined, 
       aes(x = avgPR, y =  avgPR_ED)) +
  geom_point(alpha = 0.8) +
  coord_equal() + xlim(c(0, 1)) + ylim(c(0,1)) +
  theme_bw() +
  geom_abline(linetype = 2) +
  labs(x = "Average classifier performances on EukDetect fungi\n(primary tumor 1 cancer type vs all others)",
       y = "Average classifier performances on high coverage fungi\n(primary tumor 1 cancer type vs all others)",
       title = "High coverage vs. EukDetect Fungi | AUPR\nBlood Derived Normal") +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) +
  geom_ribbon(stat='smooth', method = "lm", fill = "orange", se=TRUE, alpha=0.2, fullrange = TRUE) +
  geom_line(stat='smooth', method = "lm", color = "darkorange", alpha=1, fullrange = TRUE) +
  stat_poly_eq(formula = y ~ x, label.x = 0.9, label.y = 0.05,
               aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~")), 
               parse = TRUE)
ggsave(filename = "Figures/Supplementary_Figures/regr_aupr_tcga_BDN_high_cov_vs_EukDetect.pdf",
       dpi = "retina", units = "in", height = 5, width = 5)

#----------------------------------------------------------#
# Plot machine learning performances for all TCGA cancers on raw data split by sequencing center
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Raw <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_by_seq_center_ALL_DecontamV2_2Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_Raw$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_Raw$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_Raw <- mlPerfAll10k_Allcancer_Raw[,!(colnames(mlPerfAll10k_Allcancer_Raw) == "X")]
colnames(mlPerfAll10k_Allcancer_Raw)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_Allcancer_Raw$nullAUPR <- ifelse(mlPerfAll10k_Allcancer_Raw$minorityClassName == "SolidTissueNormal",
                                        yes=mlPerfAll10k_Allcancer_Raw$majorityClassSize/(mlPerfAll10k_Allcancer_Raw$minorityClassSize+mlPerfAll10k_Allcancer_Raw$majorityClassSize),
                                        no=mlPerfAll10k_Allcancer_Raw$minorityClassSize/(mlPerfAll10k_Allcancer_Raw$minorityClassSize+mlPerfAll10k_Allcancer_Raw$majorityClassSize))
mlPerfAll10k_Allcancer_Raw$nullAUROC <- 0.5

# Rename entries in the "datasetName" column

mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_HMS"] <- "HMS species decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_HMS_Cov_Nonzero"] <- "HMS species high coverage (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_MDA"] <- "MDA species decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_MDA_Cov_Nonzero"] <- "MDA species high coverage (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_BCM"] <- "BCM species decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_BCM_Cov_Nonzero"] <- "BCM species high coverage (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_WashU"] <- "WashU species decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_WashU_Cov_Nonzero"] <- "WashU species high coverage (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_Broad_WGS"] <- "Broad species decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Broad_WGS_Cov_Nonzero"] <- "Broad species high coverage (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_UNC"] <- "UNC species decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_UNC_Cov_Nonzero"] <- "UNC species high coverage (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_CMS"] <- "CMS species decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_CMS_Cov_Nonzero"] <- "CMS species high coverage (RNA-Seq)"
## NOTE: Broad_RNA only included GBM tumors, so no ML comparisons were made

mlPerfAll10k_Allcancer_Raw$datasetName <- factor(mlPerfAll10k_Allcancer_Raw$datasetName,
                                                    levels = c("HMS species decontaminated (WGS)",
                                                               "HMS species high coverage (WGS)",
                                                               "MDA species decontaminated (WGS)",
                                                               "MDA species high coverage (WGS)",
                                                               "BCM species decontaminated (WGS)",
                                                               "BCM species high coverage (WGS)",
                                                               "WashU species decontaminated (WGS)",
                                                               "WashU species high coverage (WGS)",
                                                               "Broad species decontaminated (WGS)",
                                                               "Broad species high coverage (WGS)",
                                                               "UNC species decontaminated (RNA-Seq)",
                                                               "UNC species high coverage (RNA-Seq)",
                                                               "CMS species decontaminated (RNA-Seq)",
                                                               "CMS species high coverage (RNA-Seq)"))

save(mlPerfAll10k_Allcancer_Raw,
     file = "Interim_data/mlPerfAll10k_Allcancer_Raw_2Apr22.RData")
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
# HMS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_HMS_PT_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# MDA
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("MDA",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_MDA_PT_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# BCM
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("BCM",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_BCM_PT_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# WashU
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("WashU",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Washington University (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_WashU_PT_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# Broad WGS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Broad",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Broad_WGS_PT_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# UNC
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("UNC",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University of North Carolina (RNA-Seq) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_UNC_PT_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# CMS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("CMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_CMS_PT_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#
# HMS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/figure_4_XX_mlPerfAll10k_rep1_HMS_PT_vs_NAT_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# MDA does not have sufficient PT and NAT samples to compare
# BCM
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_BCM_PT_vs_NAT_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# WashU does not have sufficient PT and NAT samples to compare
# Broad WGS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# UNC
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("UNC",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University of North Carolina (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_UNC_PT_vs_NAT_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# CMS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("CMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_CMS_PT_vs_NAT_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")

#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
# HMS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_HMS_BDN_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# MDA
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("MDA",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_MDA_BDN_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# BCM
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_BCM_BDN_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# WashU
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("WashU",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Washington University (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_WashU_BDN_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# Broad WGS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Broad_WGS_BDN_species_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# UNC does not have any BDN samples
# CMS does not have any BDN samples

#----------------------------------------------------------#
# Plot machine learning performances for all TCGA cancers on raw data split by sequencing center
# AND aggregated to each taxa level AND ∩ Weizmann data
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_by_seq_center_taxa_level_and_WIS_intersect_ALL_decontamV2_2Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann <- mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann[,!(colnames(mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann) == "X")]
colnames(mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$nullAUPR <- ifelse(mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$minorityClassName == "SolidTissueNormal",
                                              yes=mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$majorityClassSize/(mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$minorityClassSize+mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$majorityClassSize),
                                              no=mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$minorityClassSize/(mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$minorityClassSize+mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$majorityClassSize))
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$nullAUROC <- 0.5

# Rename entries in the "datasetName" column

# HMS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_phylum"] <- "HMS phylum decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_class"] <- "HMS class decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_order"] <- "HMS order decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_family"] <- "HMS family decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_genus"] <- "HMS genus decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_HMS_species"] <- "HMS species decontaminated (WGS)"
# HMS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_HMS_phylum_Shared_Nonzero"] <- "HMS phylum ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_HMS_class_Shared_Nonzero"] <- "HMS class ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_HMS_order_Shared_Nonzero"] <- "HMS order ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_HMS_family_Shared_Nonzero"] <- "HMS family ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_HMS_genus_Shared_Nonzero"] <- "HMS genus ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_HMS_species_Shared_Nonzero"] <- "HMS species ∩ WIS (WGS)"
# BCM - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_phylum"] <- "BCM phylum decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_class"] <- "BCM class decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_order"] <- "BCM order decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_family"] <- "BCM family decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_genus"] <- "BCM genus decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_BCM_species"] <- "BCM species decontaminated (WGS)"
# BCM - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_BCM_phylum_Shared_Nonzero"] <- "BCM phylum ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_BCM_class_Shared_Nonzero"] <- "BCM class ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_BCM_order_Shared_Nonzero"] <- "BCM order ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_BCM_family_Shared_Nonzero"] <- "BCM family ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_BCM_genus_Shared_Nonzero"] <- "BCM genus ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_BCM_species_Shared_Nonzero"] <- "BCM species ∩ WIS (WGS)"
# MDA - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_phylum"] <- "MDA phylum decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_class"] <- "MDA class decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_order"] <- "MDA order decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_family"] <- "MDA family decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_genus"] <- "MDA genus decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_MDA_species"] <- "MDA species decontaminated (WGS)"
# MDA - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_MDA_phylum_Shared_Nonzero"] <- "MDA phylum ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_MDA_class_Shared_Nonzero"] <- "MDA class ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_MDA_order_Shared_Nonzero"] <- "MDA order ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_MDA_family_Shared_Nonzero"] <- "MDA family ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_MDA_genus_Shared_Nonzero"] <- "MDA genus ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_MDA_species_Shared_Nonzero"] <- "MDA species ∩ WIS (WGS)"
# WashU - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_phylum"] <- "WashU phylum decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_class"] <- "WashU class decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_order"] <- "WashU order decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_family"] <- "WashU family decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_genus"] <- "WashU genus decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_WashU_species"] <- "WashU species decontaminated (WGS)"
# WashU - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_WashU_phylum_Shared_Nonzero"] <- "WashU phylum ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_WashU_class_Shared_Nonzero"] <- "WashU class ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_WashU_order_Shared_Nonzero"] <- "WashU order ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_WashU_family_Shared_Nonzero"] <- "WashU family ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_WashU_genus_Shared_Nonzero"] <- "WashU genus ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_WashU_species_Shared_Nonzero"] <- "WashU species ∩ WIS (WGS)"
# Broad_WGS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_phylum"] <- "Broad phylum decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_class"] <- "Broad class decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_order"] <- "Broad order decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_family"] <- "Broad family decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_genus"] <- "Broad genus decontaminated (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_Broad_WGS_species"] <- "Broad species decontaminated (WGS)"
# Broad_WGS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Broad_WGS_phylum_Shared_Nonzero"] <- "Broad phylum ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Broad_WGS_class_Shared_Nonzero"] <- "Broad class ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Broad_WGS_order_Shared_Nonzero"] <- "Broad order ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Broad_WGS_family_Shared_Nonzero"] <- "Broad family ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Broad_WGS_genus_Shared_Nonzero"] <- "Broad genus ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Broad_WGS_species_Shared_Nonzero"] <- "Broad species ∩ WIS (WGS)"
# UNC - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_phylum"] <- "UNC phylum decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_class"] <- "UNC class decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_order"] <- "UNC order decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_family"] <- "UNC family decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_genus"] <- "UNC genus decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_UNC_species"] <- "UNC species decontaminated (RNA-Seq)"
# UNC - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_UNC_phylum_Shared_Nonzero"] <- "UNC phylum ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_UNC_class_Shared_Nonzero"] <- "UNC class ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_UNC_order_Shared_Nonzero"] <- "UNC order ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_UNC_family_Shared_Nonzero"] <- "UNC family ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_UNC_genus_Shared_Nonzero"] <- "UNC genus ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_UNC_species_Shared_Nonzero"] <- "UNC species ∩ WIS (RNA-Seq)"
# CMS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_phylum"] <- "CMS phylum decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_class"] <- "CMS class decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_order"] <- "CMS order decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_family"] <- "CMS family decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_genus"] <- "CMS genus decontaminated (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_DecontamV2_CMS_species"] <- "CMS species decontaminated (RNA-Seq)"
# CMS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_CMS_phylum_Shared_Nonzero"] <- "CMS phylum ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_CMS_class_Shared_Nonzero"] <- "CMS class ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_CMS_order_Shared_Nonzero"] <- "CMS order ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_CMS_family_Shared_Nonzero"] <- "CMS family ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_CMS_genus_Shared_Nonzero"] <- "CMS genus ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_CMS_species_Shared_Nonzero"] <- "CMS species ∩ WIS (RNA-Seq)"
## NOTE: Broad_RNA only included GBM tumors, so no ML comparisons were made

mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName <- factor(mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName,
                                                 levels = c("HMS phylum decontaminated (WGS)",
                                                            "HMS class decontaminated (WGS)",
                                                            "HMS order decontaminated (WGS)",
                                                            "HMS family decontaminated (WGS)",
                                                            "HMS genus decontaminated (WGS)",
                                                            "HMS species decontaminated (WGS)",
                                                            "HMS phylum ∩ WIS (WGS)",
                                                            "HMS class ∩ WIS (WGS)",
                                                            "HMS order ∩ WIS (WGS)",
                                                            "HMS family ∩ WIS (WGS)",
                                                            "HMS genus ∩ WIS (WGS)",
                                                            "HMS species ∩ WIS (WGS)",
                                                            # BCM
                                                            "BCM phylum decontaminated (WGS)",
                                                            "BCM class decontaminated (WGS)",
                                                            "BCM order decontaminated (WGS)",
                                                            "BCM family decontaminated (WGS)",
                                                            "BCM genus decontaminated (WGS)",
                                                            "BCM species decontaminated (WGS)",
                                                            "BCM phylum ∩ WIS (WGS)",
                                                            "BCM class ∩ WIS (WGS)",
                                                            "BCM order ∩ WIS (WGS)",
                                                            "BCM family ∩ WIS (WGS)",
                                                            "BCM genus ∩ WIS (WGS)",
                                                            "BCM species ∩ WIS (WGS)",
                                                            # MDA
                                                            "MDA phylum decontaminated (WGS)",
                                                            "MDA class decontaminated (WGS)",
                                                            "MDA order decontaminated (WGS)",
                                                            "MDA family decontaminated (WGS)",
                                                            "MDA genus decontaminated (WGS)",
                                                            "MDA species decontaminated (WGS)",
                                                            "MDA phylum ∩ WIS (WGS)",
                                                            "MDA class ∩ WIS (WGS)",
                                                            "MDA order ∩ WIS (WGS)",
                                                            "MDA family ∩ WIS (WGS)",
                                                            "MDA genus ∩ WIS (WGS)",
                                                            "MDA species ∩ WIS (WGS)",
                                                            # WashU
                                                            "WashU phylum decontaminated (WGS)",
                                                            "WashU class decontaminated (WGS)",
                                                            "WashU order decontaminated (WGS)",
                                                            "WashU family decontaminated (WGS)",
                                                            "WashU genus decontaminated (WGS)",
                                                            "WashU species decontaminated (WGS)",
                                                            "WashU phylum ∩ WIS (WGS)",
                                                            "WashU class ∩ WIS (WGS)",
                                                            "WashU order ∩ WIS (WGS)",
                                                            "WashU family ∩ WIS (WGS)",
                                                            "WashU genus ∩ WIS (WGS)",
                                                            "WashU species ∩ WIS (WGS)",
                                                            # Broad
                                                            "Broad phylum decontaminated (WGS)",
                                                            "Broad class decontaminated (WGS)",
                                                            "Broad order decontaminated (WGS)",
                                                            "Broad family decontaminated (WGS)",
                                                            "Broad genus decontaminated (WGS)",
                                                            "Broad species decontaminated (WGS)",
                                                            "Broad phylum ∩ WIS (WGS)",
                                                            "Broad class ∩ WIS (WGS)",
                                                            "Broad order ∩ WIS (WGS)",
                                                            "Broad family ∩ WIS (WGS)",
                                                            "Broad genus ∩ WIS (WGS)",
                                                            "Broad species ∩ WIS (WGS)",
                                                            # UNC
                                                            "UNC phylum decontaminated (RNA-Seq)",
                                                            "UNC class decontaminated (RNA-Seq)",
                                                            "UNC order decontaminated (RNA-Seq)",
                                                            "UNC family decontaminated (RNA-Seq)",
                                                            "UNC genus decontaminated (RNA-Seq)",
                                                            "UNC species decontaminated (RNA-Seq)",
                                                            "UNC phylum ∩ WIS (RNA-Seq)",
                                                            "UNC class ∩ WIS (RNA-Seq)",
                                                            "UNC order ∩ WIS (RNA-Seq)",
                                                            "UNC family ∩ WIS (RNA-Seq)",
                                                            "UNC genus ∩ WIS (RNA-Seq)",
                                                            "UNC species ∩ WIS (RNA-Seq)",
                                                            # CMS
                                                            "CMS phylum decontaminated (RNA-Seq)",
                                                            "CMS class decontaminated (RNA-Seq)",
                                                            "CMS order decontaminated (RNA-Seq)",
                                                            "CMS family decontaminated (RNA-Seq)",
                                                            "CMS genus decontaminated (RNA-Seq)",
                                                            "CMS species decontaminated (RNA-Seq)",
                                                            "CMS phylum ∩ WIS (RNA-Seq)",
                                                            "CMS class ∩ WIS (RNA-Seq)",
                                                            "CMS order ∩ WIS (RNA-Seq)",
                                                            "CMS family ∩ WIS (RNA-Seq)",
                                                            "CMS genus ∩ WIS (RNA-Seq)",
                                                            "CMS species ∩ WIS (RNA-Seq)",
                                                            # Broad
                                                            "Broad phylum decontaminated (RNA-Seq)",
                                                            "Broad class decontaminated (RNA-Seq)",
                                                            "Broad order decontaminated (RNA-Seq)",
                                                            "Broad family decontaminated (RNA-Seq)",
                                                            "Broad genus decontaminated (RNA-Seq)",
                                                            "Broad species decontaminated (RNA-Seq)",
                                                            "Broad phylum ∩ WIS (RNA-Seq)",
                                                            "Broad class ∩ WIS (RNA-Seq)",
                                                            "Broad order ∩ WIS (RNA-Seq)",
                                                            "Broad family ∩ WIS (RNA-Seq)",
                                                            "Broad genus ∩ WIS (RNA-Seq)",
                                                            "Broad species ∩ WIS (RNA-Seq)"))

save(mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann, 
     file = "Interim_data/mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann_2Apr22.RData")
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
# HMS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data\nPrimary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_HMS_PT_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 5, units = "in")
# HMS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data ∩ WIS Features\nPrimary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_HMS_PT_taxa_level_shared_decontamV2.svg", dpi = "retina",
       width = 10, height = 5, units = "in")

# BCM - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("BCM",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data\nPrimary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_BCM_PT_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 5, units = "in")
# BCM - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("BCM",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data ∩ WIS Features  | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_BCM_PT_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

# MDA - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("MDA",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data\nPrimary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_MDA_PT_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 5, units = "in")
# MDA - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("MDA",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data ∩ WIS Features  | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_MDA_PT_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

# WashU - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("WashU",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data\nPrimary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_WashU_PT_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 5, units = "in")
# WashU - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("WashU",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data ∩ WIS Features  | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_WashU_PT_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

# Broad_WGS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Broad",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data\nPrimary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Broad_WGS_PT_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 5, units = "in")
# Broad_WGS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Broad",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data ∩ WIS Features  | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Broad_WGS_PT_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

# UNC - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("UNC",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data\nPrimary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_UNC_PT_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 18, height = 6, units = "in")
# UNC - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("UNC",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data ∩ WIS Features  | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_UNC_PT_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

# CMS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("CMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data\nPrimary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_CMS_PT_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 5, units = "in")
# CMS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("CMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data ∩ WIS Features  | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_CMS_PT_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#
# HMS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_HMS_PT_vs_NAT_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# HMS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data ∩ WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") 
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_HMS_PT_vs_NAT_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

# BCM - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_BCM_PT_vs_NAT_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# BCM - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data ∩ WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_BCM_PT_vs_NAT_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

## NOTE: Neither MDA nor WashU has enough tumor vs. normal samples to plot

# Broad_WGS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# Broad_WGS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data ∩ WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

# UNC - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("UNC",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_UNC_PT_vs_NAT_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# UNC - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("UNC",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) + 
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),color="#0072B5FF",lty="dotted") + 
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="#BC3C29FF",lty="dotted") +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data ∩ WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_UNC_PT_vs_NAT_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

# CMS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("CMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_CMS_PT_vs_NAT_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 4, units = "in")
# CMS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("CMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data ∩ WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_CMS_PT_vs_NAT_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
# HMS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Blood Derived Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_HMS_BDN_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 5, units = "in")
# HMS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data ∩ WIS Features  | Blood Derived Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_HMS_BDN_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

# BCM - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Blood Derived Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_BCM_BDN_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 5, units = "in")
# BCM - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data ∩ WIS Features  | Blood Derived Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_BCM_BDN_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

# MDA - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("MDA",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Blood Derived Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_MDA_BDN_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 5, units = "in")
# MDA - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("MDA",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data ∩ WIS Features  | Blood Derived Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_MDA_BDN_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

# WashU - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("WashU",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data | Blood Derived Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_WashU_BDN_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 5, units = "in")
# WashU - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("WashU",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data ∩ WIS Features  | Blood Derived Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_WashU_BDN_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

# Broad_WGS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Blood Derived Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Broad_WGS_BDN_taxa_level_decontam_decontamV2.pdf", dpi = "retina",
         width = 10, height = 5, units = "in")
# Broad_WGS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data ∩ WIS Features  | Blood Derived Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Broad_WGS_BDN_taxa_level_shared_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")

#----------------------------------------------------------#
# Overlay plots machine learning performances for all TCGA cancers
# using RAW data per seq center:
# - Decontam species data
# - High coverage species
# - Intersected species with WIS
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

# NOTE: There are too many cancer types to include more than two errorbars per cancer type,
# so the WGS vs. RNA comparison is not included in this overlay. If you would like to include
# it, you can remove the commented line below. ## mlPerfAll10k_Allcancer_Raw
mlPerfAll10k_Allcancer_Raw_Overlay <- rbind(mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann,
                                            mlPerfAll10k_Allcancer_Raw,
                                            mlPerfAll10k_Allcancer_ED_HMP_VSNM)
mlPerfAll10k_Allcancer_Raw_Overlay_Filt <- mlPerfAll10k_Allcancer_Raw_Overlay %>%
  filter(grepl("species",datasetName)) %>% droplevels()
mlPerfAll10k_Allcancer_Raw_Overlay_Filt$datasetName <- factor(mlPerfAll10k_Allcancer_Raw_Overlay_Filt$datasetName,
                                                              levels = rev(levels(mlPerfAll10k_Allcancer_Raw_Overlay_Filt$datasetName)))

#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
# HMS - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Main_Figures/mlPerfAll10k_rep1_HMS_PT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# HMS - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Main_Figures/mlPerfAll10k_rep1_HMS_PT_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# BCM - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("BCM",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_BCM_PT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# BCM - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("BCM",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_BCM_PT_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# MDA - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("MDA",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_MDA_PT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# MDA - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("MDA",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_MDA_PT_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# WashU - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("WashU",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_WashU_PT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# WashU - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("WashU",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_WashU_PT_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# Broad_WGS - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Broad",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Broad_WGS_PT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# Broad_WGS - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Broad",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Broad_WGS_PT_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# UNC - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("UNC",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=0.8) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_UNC_PT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# UNC - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("UNC",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University of North Carolina (RNA-Seq) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_UNC_PT_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# CMS - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("CMS",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data\n Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_CMS_PT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# CMS - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("CMS",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data\n Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_CMS_PT_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#
# HMS - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_HMS_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# HMS - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_HMS_PT_vs_NAT_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# BCM - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_BCM_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# BCM - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_BCM_PT_vs_NAT_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")
## NOTE: Neither MDA nor WashU had sufficient samples to plot primary tumor vs. NAT performance

# Broad_WGS - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# Broad_WGS - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# UNC - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("UNC",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data\nPrimary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_UNC_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# UNC - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("UNC",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University of North Carolina (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_UNC_PT_vs_NAT_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# CMS - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("CMS",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data\nPrimary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_CMS_PT_vs_NAT_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# CMS - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("CMS",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data\nPrimary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_CMS_PT_vs_NAT_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
# HMS - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(30) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Main_Figures/mlPerfAll10k_rep1_HMS_BDN_species_overlay_decontamV2.svg", dpi = "retina",
         width = 8.5, height = 4, units = "in")
# HMS - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Main_Figures/mlPerfAll10k_rep1_HMS_BDN_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# BCM - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_BCM_BDN_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# BCM - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_BCM_BDN_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# MDA - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("MDA",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_MDA_BDN_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# MDA - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("MDA",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_MDA_BDN_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# WashU - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("WashU",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_WashU_BDN_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# WashU - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("WashU",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_WashU_BDN_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

# Broad_WGS - overlay without EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  filter(!grepl("EukDetect|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Broad_WGS_BDN_species_overlay_decontamV2.svg", dpi = "retina",
         width = 10, height = 4, units = "in")
# Broad_WGS - overlay with EukDetect
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  filter(!grepl("Unfiltered|HMP",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),color="darkgray",lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Broad_WGS_BDN_species_overlay_decontamV2_EukDetect.svg", dpi = "retina",
       width = 10, height = 4, units = "in")

#-------------------------Regress avg performance vs minority class size using raw data-------------------------#

ptRawDataPerfSummarizedMinClassSize <- mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>% distinct() %>%
  filter(grepl("decontaminated",datasetName)) %>%
  filter(sampleType == "Primary Tumor") %>%
  group_by(sampleType, metadataName, diseaseType) %>% 
  mutate(avgROC = mean(AUROC), avgPR = mean(AUPR)) %>%
  select(sampleType, diseaseType, abbrev, avgROC, avgPR, minorityClassSize, majorityClassSize, metadataName) %>% 
  mutate(logAvgMinClass = log10(minorityClassSize), ratioMinMajClass=minorityClassSize/majorityClassSize) %>%
  unique() %>% data.frame()

ptRawDataPerfSummarizedMinClassSize %>% filter(abbrev == "OV")
# Add seq center names
ptRawDataPerfSummarizedMinClassSize$data_submitting_center_label <- NA
ptRawDataPerfSummarizedMinClassSize$data_submitting_center_label[grepl("HMS",ptRawDataPerfSummarizedMinClassSize$metadataName)] <- "Harvard Medical School"
ptRawDataPerfSummarizedMinClassSize$data_submitting_center_label[grepl("BCM",ptRawDataPerfSummarizedMinClassSize$metadataName)] <- "Baylor College of Medicine"
ptRawDataPerfSummarizedMinClassSize$data_submitting_center_label[grepl("MDA",ptRawDataPerfSummarizedMinClassSize$metadataName)] <- "MD Anderson - Institute for Applied Cancer Science"
ptRawDataPerfSummarizedMinClassSize$data_submitting_center_label[grepl("WashU",ptRawDataPerfSummarizedMinClassSize$metadataName)] <- "Washington University School of Medicine"
ptRawDataPerfSummarizedMinClassSize$data_submitting_center_label[grepl("Broad_WGS",ptRawDataPerfSummarizedMinClassSize$metadataName)] <- "Broad Institute of MIT and Harvard"
ptRawDataPerfSummarizedMinClassSize$data_submitting_center_label[grepl("UNC",ptRawDataPerfSummarizedMinClassSize$metadataName)] <- "University of North Carolina"
ptRawDataPerfSummarizedMinClassSize$data_submitting_center_label[grepl("CMS",ptRawDataPerfSummarizedMinClassSize$metadataName)] <- "Canada's Michael Smith Genome Sciences Centre"
# Add exp strategy names
ptRawDataPerfSummarizedMinClassSize$experimental_strategy <- NA
ptRawDataPerfSummarizedMinClassSize$experimental_strategy[grepl("HMS|BCM|MDA|WashU|Broad_WGS",ptRawDataPerfSummarizedMinClassSize$metadataName)] <- "WGS"
ptRawDataPerfSummarizedMinClassSize$experimental_strategy[grepl("UNC|CMS",ptRawDataPerfSummarizedMinClassSize$metadataName)] <- "RNA-Seq"

## Regress using minority class size
summary(lm(avgPR ~ minorityClassSize, ptRawDataPerfSummarizedMinClassSize))
summary(lm(avgROC ~ minorityClassSize, ptRawDataPerfSummarizedMinClassSize))

require(ggrepel)
ptRawDataPerfSummarizedMinClassSize %>%
  # filter(grepl("HMS|MDA|BCM|WashU|WGS",metadataName)) %>%
  ggplot(aes(x=minorityClassSize, y=avgPR, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.x = 400, label.y = 0.02) +
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  geom_text_repel() +
  labs(x = "Minority class size (per cancer type)", y = "Average AUPR", 
       title = "Correlating AUPR and minority class size\nfor primary tumor TCGA models (Raw data)")
ggsave(filename = "Figures/Other_Figures/aupr_vs_minority_class_size_all_TCGA_PT_raw_data_AllSeqCenters.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

ptRawDataPerfSummarizedMinClassSize %>%
  ggplot(aes(x=minorityClassSize, y=avgROC, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.x = 400, label.y = 0.02) +
  geom_text_repel() +
  theme_bw()+ theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  labs(x = "Minority class size (per cancer type)", y = "Average AUROC", 
       title = "Correlating AUROC and minority class size\nfor primary tumor TCGA models (Raw data)")
ggsave(filename = "Figures/Other_Figures/auroc_vs_minority_class_size_all_TCGA_PT_raw_data_AllSeqCenters.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

## Regress using minority class *ratio*
summary(lm(avgPR ~ ratioMinMajClass, ptRawDataPerfSummarizedMinClassSize))
summary(lm(avgROC ~ ratioMinMajClass, ptRawDataPerfSummarizedMinClassSize))

ptRawDataPerfSummarizedMinClassSize %>%
  # filter(grepl("HMS|MDA|BCM|WashU|WGS",metadataName)) %>%
  ggplot(aes(x=ratioMinMajClass, y=avgPR, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.y = 1.1) + 
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  geom_text_repel() +
  labs(x = "Ratio of minority to majority class size (per cancer type)", y = "Average AUPR", 
       title = "Correlating AUPR and ratio of minority to majority class\nsize for primary tumor TCGA models (Raw data)")
ggsave(filename = "Figures/Supplementary_Figures/aupr_vs_ratio_minority_to_majority_class_size_all_TCGA_PT_raw_data_AllSeqCenters.pdf",
         dpi = "retina", units = "in", width = 5, height = 5)

ptRawDataPerfSummarizedMinClassSize %>%
  ggplot(aes(x=ratioMinMajClass, y=avgROC, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.y = 1.05) + 
  geom_text_repel() +
  theme_bw()+ theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  labs(x = "Ratio of minority to majority class size (per cancer type)", y = "Average AUROC", 
       title = "Correlating AUROC and ratio of minority to majority class size\nfor primary tumor TCGA models (Raw data)")
ggsave(filename = "Figures/Supplementary_Figures/auroc_vs_ratio_minority_to_majority_class_size_all_TCGA_PT_raw_data_AllSeqCenters.pdf",
         dpi = "retina", units = "in", width = 5, height = 5)

#-------------------------Regress avg performance vs average read depth-------------------------#
load("Interim_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_25Mar22.RData",
     verbose = TRUE)

metaAvgReadDepthRaw <- droplevels(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts[rownames(metaQiitaCombined_Nonzero_DecontamV2),])
metaAvgReadDepthRaw %>% 
  filter(sample_type == "Primary Tumor") %>%
  mutate(abbrev = gsub("^TCGA-","",investigation)) %>% group_by(data_submitting_center_label, experimental_strategy, abbrev) %>%
  summarize(Mean_bam_total_reads = mean(bam_total_reads), 
            Mean_bam_mapped_reads = mean(bam_mapped_reads),
            Mean_bam_unmapped_reads = mean(bam_unmapped_reads),
            Mean_bam_ratio_unmapped = mean(bam_ratio_unmapped)) %>% data.frame() -> metaAvgReadDepthRawFormatted

ptRawDataPerfSummarizedReadDepth <- ptRawDataPerfSummarizedMinClassSize %>%
  left_join(metaAvgReadDepthRawFormatted, by = c("data_submitting_center_label","experimental_strategy","abbrev")) %>% 
  data.frame()

## Mean total bam reads
summary(lm(avgPR ~ Mean_bam_total_reads, ptRawDataPerfSummarizedReadDepth))
summary(lm(avgROC ~ Mean_bam_total_reads, ptRawDataPerfSummarizedReadDepth))

require(ggrepel)
ptRawDataPerfSummarizedReadDepth %>%
  ggplot(aes(x=Mean_bam_total_reads, y=avgPR, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(aes(color = NULL), method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.x = 5e8, label.y = 1) + 
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  geom_text_repel() +
  labs(x = "Average Total Read Depth Per\nSeq Center & Cancer Type (Primary Tumor)", 
       y = "Average AUPR", 
       title = "Correlating AUPR and read depth\nfor primary tumor TCGA models (Raw data)")
ggsave(filename = "Figures/Other_Figures/aupr_vs_read_depth_all_TCGA_PT_raw_data_AllSeqCenters.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

ptRawDataPerfSummarizedReadDepth %>%
  ggplot(aes(x=Mean_bam_total_reads, y=avgROC, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.x = 5e8, label.y = 1.05) + 
  geom_text_repel() +
  theme_bw()+ theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  labs(x = "Average Total Read Depth Per\nSeq Center & Cancer Type (Primary Tumor)", 
       y = "Average AUROC", 
       title = "Correlating AUROC and read depth\nfor primary tumor TCGA models (Raw data)")
ggsave(filename = "Figures/Other_Figures/auroc_vs_read_depth_all_TCGA_PT_raw_data_AllSeqCenters.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

## Mean unmapped bam reads
summary(lm(avgPR ~ Mean_bam_unmapped_reads, ptRawDataPerfSummarizedReadDepth))
summary(lm(avgROC ~ Mean_bam_unmapped_reads, ptRawDataPerfSummarizedReadDepth))

ptRawDataPerfSummarizedReadDepth %>%
  ggplot(aes(x=Mean_bam_unmapped_reads, y=avgPR, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red", label.x = 5e7, label.y = 0.15) + 
  theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  geom_text_repel() +
  labs(x = "Average Unmapped Reads Per\nSeq Center & Cancer Type (Primary Tumor)", 
       y = "Average AUPR", 
       title = "Correlating AUPR and unmapped reads\nfor primary tumor TCGA models (Raw data)")
ggsave(filename = "Figures/Other_Figures/aupr_vs_unmapped_reads_all_TCGA_PT_raw_data_AllSeqCenters.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

ptRawDataPerfSummarizedReadDepth %>%
  ggplot(aes(x=Mean_bam_unmapped_reads, y=avgROC, label=abbrev)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + 
  stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE, color="red",label.y = 1.05) + 
  geom_text_repel() +
  theme_bw()+ theme(aspect.ratio=1, plot.title = element_text(hjust=0.5)) + coord_fixed() +
  labs(x = "Average Unmapped Reads Per\nSeq Center & Cancer Type (Primary Tumor)", 
       y = "Average AUROC", 
       title = "Correlating AUROC and unmapped reads\nfor primary tumor TCGA models (Raw data)")
ggsave(filename = "Figures/Other_Figures/auroc_vs_unmapped_reads_all_TCGA_PT_raw_data_AllSeqCenters.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

#-------------------------Plot WGS vs RNA-Seq performance-------------------------#
ptRawDataPerfSummarizedMinClassSize %>%
  ggboxplot(x = "experimental_strategy", 
            y = "avgROC",
            notch = TRUE,
            fill = "experimental_strategy",
            xlab = "Experimental strategy",
            ylab = "Average ROC\n(one per cancer type per center)",
            add = "jitter",
            legend = "none",
            add.params = list(alpha=0.4),
            palette = "nejm") +
  stat_compare_means(comparisons = list(c("WGS","RNA-Seq"))) -> avgROCPlot

ptRawDataPerfSummarizedMinClassSize %>%
  ggboxplot(x = "experimental_strategy",
            y = "avgPR",
            notch = TRUE,
            fill = "experimental_strategy",
            xlab = "Experimental strategy",
            ylab = "Average PR\n(one per cancer type per center)",
            legend = "none",
            add = "jitter",
            add.params = list(alpha=0.4),
            palette = "nejm") +
  stat_compare_means(comparisons = list(c("WGS","RNA-Seq"))) -> avgPRPlot

ggarrange(avgROCPlot, avgPRPlot, ncol = 2)
ggsave(filename = "Figures/Supplementary_Figures/auroc_aupr_raw_data_wgs_vs_rna.pdf",
         dpi = "retina", units = "in", width = 5, height = 4)

#----------------------------------------------------------#
# Plot topX machine learning performances for all TCGA cancers - VSNM data
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_TopX_VSNM <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_vsnm_topX_ALL_2Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_TopX_VSNM$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_TopX_VSNM$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_TopX_VSNM <- mlPerfAll10k_Allcancer_TopX_VSNM[,!(colnames(mlPerfAll10k_Allcancer_TopX_VSNM) == "X")]
colnames(mlPerfAll10k_Allcancer_TopX_VSNM)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_Allcancer_TopX_VSNM$nullAUPR <- ifelse(mlPerfAll10k_Allcancer_TopX_VSNM$minorityClassName == "SolidTissueNormal",
                                          yes=mlPerfAll10k_Allcancer_TopX_VSNM$majorityClassSize/(mlPerfAll10k_Allcancer_TopX_VSNM$minorityClassSize+mlPerfAll10k_Allcancer_TopX_VSNM$majorityClassSize),
                                          no=mlPerfAll10k_Allcancer_TopX_VSNM$minorityClassSize/(mlPerfAll10k_Allcancer_TopX_VSNM$minorityClassSize+mlPerfAll10k_Allcancer_TopX_VSNM$majorityClassSize))
mlPerfAll10k_Allcancer_TopX_VSNM$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
table(mlPerfAll10k_Allcancer_TopX_VSNM$datasetName)
mlPerfAll10k_Allcancer_TopX_VSNM$datasetName[mlPerfAll10k_Allcancer_TopX_VSNM$datasetName == "rep200FungiDecontamV2SpeciesVSNM_TopX"] <- "Top 20 fungi from Hopkins cohort"

mlPerfAll10k_Allcancer_TopX_VSNM$datasetName <- factor(mlPerfAll10k_Allcancer_TopX_VSNM$datasetName,
                                             levels = c("Top 20 fungi from Hopkins cohort"))

#-------------------------Plot ML performances-------------------------#
## PT
mlPerfAll10k_Allcancer_TopX_VSNM %>%
  filter(sampleType == "Primary Tumor") %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Top 20 fungi from Hopkins cohort | Primary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_TopX.pdf", dpi = "retina",
         width = 9, height = 4, units = "in")

## Primary Tumor vs NAT
mlPerfAll10k_Allcancer_TopX_VSNM %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Top 20 fungi from Hopkins cohort | Primary Tumor vs NAT") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_NAT_TopX.pdf", dpi = "retina",
       width = 8, height = 4, units = "in")

## BDN
mlPerfAll10k_Allcancer_TopX_VSNM %>%
  filter(sampleType == "Blood Derived Normal") %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Top 20 fungi from Hopkins cohort | Blood Derived Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_BDN_TopX.pdf", dpi = "retina",
       width = 8, height = 4, units = "in")

require(gmodels)
mlPerfAll10k_Allcancer_TopX_VSNM %>%
  filter(sampleType == "Blood Derived Normal") %>%
  distinct() %>% droplevels() %>%
  pull(AUROC) %>% ci()
# Estimate   CI lower   CI upper Std. Error 
# 0.88775374 0.87757196 0.89793553 0.00516329

#----------------------------------------------------------#
# Plot topX machine learning performances for all TCGA cancers - Raw data
#----------------------------------------------------------#

source("Supporting_scripts/S00-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_TopX_Raw <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_raw_topX_by_seq_center_ALL_DecontamV2_2Apr22.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_TopX_Raw$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_TopX_Raw$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_TopX_Raw <- mlPerfAll10k_Allcancer_TopX_Raw[,!(colnames(mlPerfAll10k_Allcancer_TopX_Raw) == "X")]
colnames(mlPerfAll10k_Allcancer_TopX_Raw)[1:2] <- c("AUROC","AUPR")
# Add null perf values. Note: AUPR null is prevalence of **positive class**
# For 1-vs-all-others, prevalence is (minority class)/(total samples)
# For PT vs. NAT, prevalence is (majority class [PT])/(total samples)
mlPerfAll10k_Allcancer_TopX_Raw$nullAUPR <- ifelse(mlPerfAll10k_Allcancer_TopX_Raw$minorityClassName == "SolidTissueNormal",
                                                    yes=mlPerfAll10k_Allcancer_TopX_Raw$majorityClassSize/(mlPerfAll10k_Allcancer_TopX_Raw$minorityClassSize+mlPerfAll10k_Allcancer_TopX_Raw$majorityClassSize),
                                                    no=mlPerfAll10k_Allcancer_TopX_Raw$minorityClassSize/(mlPerfAll10k_Allcancer_TopX_Raw$minorityClassSize+mlPerfAll10k_Allcancer_TopX_Raw$majorityClassSize))
mlPerfAll10k_Allcancer_TopX_Raw$nullAUROC <- 0.5

# Rename entries in the "datasetName" column
table(mlPerfAll10k_Allcancer_TopX_Raw$datasetName)
mlPerfAll10k_Allcancer_TopX_Raw$datasetName[mlPerfAll10k_Allcancer_TopX_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_TopX_HMS"] <- "HMS top 20 fungi from Hopkins cohort"
mlPerfAll10k_Allcancer_TopX_Raw$datasetName[mlPerfAll10k_Allcancer_TopX_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_TopX_BCM"] <- "BCM top 20 fungi from Hopkins cohort"
mlPerfAll10k_Allcancer_TopX_Raw$datasetName[mlPerfAll10k_Allcancer_TopX_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_TopX_MDA"] <- "MDA top 20 fungi from Hopkins cohort"
mlPerfAll10k_Allcancer_TopX_Raw$datasetName[mlPerfAll10k_Allcancer_TopX_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_TopX_WashU"] <- "WashU top 20 fungi from Hopkins cohort"
mlPerfAll10k_Allcancer_TopX_Raw$datasetName[mlPerfAll10k_Allcancer_TopX_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_TopX_Broad_WGS_Nonzero"] <- "Broad top 20 fungi from Hopkins cohort"
mlPerfAll10k_Allcancer_TopX_Raw$datasetName[mlPerfAll10k_Allcancer_TopX_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_TopX_UNC_Nonzero"] <- "UNC top 20 fungi from Hopkins cohort"
mlPerfAll10k_Allcancer_TopX_Raw$datasetName[mlPerfAll10k_Allcancer_TopX_Raw$datasetName == "rep200_HiSeq_Fungi_DecontamV2_TopX_CMS_Nonzero"] <- "CMS top 20 fungi from Hopkins cohort"

mlPerfAll10k_Allcancer_TopX_Raw$datasetName <- factor(mlPerfAll10k_Allcancer_TopX_Raw$datasetName,
                                                       levels = c("HMS top 20 fungi from Hopkins cohort",
                                                                  "BCM top 20 fungi from Hopkins cohort",
                                                                  "MDA top 20 fungi from Hopkins cohort",
                                                                  "WashU top 20 fungi from Hopkins cohort",
                                                                  "Broad top 20 fungi from Hopkins cohort",
                                                                  "UNC top 20 fungi from Hopkins cohort",
                                                                  "CMS top 20 fungi from Hopkins cohort"))

#-------------------------Plot PT ML performances-------------------------#
# HMS
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("HMS | Raw Data | Top 20 fungi from Hopkins cohort\nPrimary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_TopX_Raw_HMS.pdf", dpi = "retina",
       width = 7, height = 4, units = "in")
# BCM
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("BCM",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("BCM | Raw Data | Top 20 fungi from Hopkins cohort\nPrimary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_TopX_Raw_BCM.pdf", dpi = "retina",
       width = 7, height = 4, units = "in")

# MDA
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("MDA",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MDA | Raw Data | Top 20 fungi from Hopkins cohort\nPrimary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_TopX_Raw_MDA.pdf", dpi = "retina",
       width = 7, height = 4, units = "in")

# WashU
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("WashU",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU | Raw Data | Top 20 fungi from Hopkins cohort\nPrimary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_TopX_Raw_WashU.pdf", dpi = "retina",
       width = 7, height = 4, units = "in")

# UNC
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("UNC",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("UNC | Raw Data | Top 20 fungi from Hopkins cohort\nPrimary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_TopX_Raw_UNC.pdf", dpi = "retina",
       width = 9, height = 4, units = "in")

# CMS
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("CMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("CMS | Raw Data | Top 20 fungi from Hopkins cohort\nPrimary Tumor | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PT_TopX_Raw_CMS.pdf", dpi = "retina",
       width = 7, height = 4, units = "in")

#-------------------------Plot Primary Tumor vs NAT ML performances-------------------------#
# HMS
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("HMS | Raw Data | Top 20 fungi from Hopkins cohort\nPrimary Tumor vs Solid Tissue Normal") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PTvsNAT_TopX_Raw_HMS.pdf", dpi = "retina",
       width = 6, height = 4, units = "in")
# BCM
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("BCM | Raw Data | Top 20 fungi from Hopkins cohort\nPrimary Tumor vs Solid Tissue Normal") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PTvsNAT_TopX_Raw_BCM.pdf", dpi = "retina",
       width = 6, height = 4, units = "in")

# MDA and WashU did not have sufficient samples

# UNC
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("UNC",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("UNC | Raw Data | Top 20 fungi from Hopkins cohort\nPrimary Tumor vs Solid Tissue Normal") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PTvsNAT_TopX_Raw_UNC.pdf", dpi = "retina",
       width = 7, height = 4, units = "in")

# CMS
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("CMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("CMS | Raw Data | Top 20 fungi from Hopkins cohort\nPrimary Tumor vs Solid Tissue Normal") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(0) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_PTvsNAT_TopX_Raw_CMS.pdf", dpi = "retina",
       width = 5, height = 4, units = "in")

#-------------------------Plot BDN ML performances-------------------------#
# HMS
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("HMS | Raw Data | Top 20 fungi from Hopkins cohort\nBlood Derived Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_BDN_TopX_Raw_HMS.pdf", dpi = "retina",
       width = 6, height = 4, units = "in")
# BCM
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("BCM | Raw Data | Top 20 fungi from Hopkins cohort\nBlood Derived Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_BDN_TopX_Raw_BCM.pdf", dpi = "retina",
       width = 6, height = 4, units = "in")

# MDA
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("MDA",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MDA | Raw Data | Top 20 fungi from Hopkins cohort\nBlood Derived Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_BDN_TopX_Raw_MDA.pdf", dpi = "retina",
       width = 6, height = 4, units = "in")

# WashU
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("WashU",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="dotted",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU | Raw Data | Top 20 fungi from Hopkins cohort\nBlood Derived Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_BDN_TopX_Raw_WashU.pdf", dpi = "retina",
       width = 6, height = 4, units = "in")

# All
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  # filter(grepl("WashU",datasetName)) %>%
  distinct() %>% droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","metadataName","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName","nullAUPR","nullAUROC")) %>%
  mutate(nullAUPR = ifelse(variable=="AUROC",NA,nullAUPR), nullAUROC = ifelse(variable=="AUPR",NA,nullAUROC)) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=ifelse(value-ci<0,0,value-ci), ymax=ifelse(value+ci>1,1,value+ci)),width=0.4,size=0.6,position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUPR,ymin=nullAUPR,ymax=nullAUPR),lty="solid",position = position_dodge(0.9)) +
  geom_errorbar(aes(y=nullAUROC,ymin=nullAUROC,ymax=nullAUROC),lty="dotted",position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=1) +
  xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All Seq Centers | Raw Data | Top 20 fungi from Hopkins cohort\nBlood Derived Normal | 1 Vs All") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed")
ggsave("Figures/Other_Figures/mlPerfAll10k_rep1_Allcancers_BDN_TopX_Raw_AllSeqCenters.pdf", dpi = "retina",
       width = 10, height = 5, units = "in")

require(gmodels)
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  distinct() %>% droplevels() %>%
  pull(AUROC) %>% ci()
# Estimate    CI lower    CI upper  Std. Error 
# 0.737488652 0.718456862 0.756520442 0.009674554

require(gmodels)
mlPerfAll10k_Allcancer_TopX_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  distinct() %>% droplevels() %>%
  pull(AUROC) %>% ci()
# Estimate    CI lower    CI upper  Std. Error 
# 0.724949184 0.714574859 0.735323510 0.005283249
