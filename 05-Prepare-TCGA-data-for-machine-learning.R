#-----------------------------------------------------------------------------
# 05-Prepare-TCGA-data-for-machine-learning.R
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
# snmDataOGUFungi,
# vdge_dataE,
# rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero,
# metaQiitaCombined_Nonzero,
load("Interim_data/snmDataFungi_13Sep21.RData") # To load the metadata and raw data objects
rep200TaxSplit <- read.csv("Supporting_data/rep200_lineage_map_split.csv", stringsAsFactors = FALSE, row.names = 1)

#--------------------------------------------------------------------------------------------------------------------#
# Separate raw data into seq center-experimental strategy groups (to preclude needing batch correction)
#--------------------------------------------------------------------------------------------------------------------#
metaQiitaCombined_Nonzero %>% count(data_submitting_center_label, experimental_strategy)

#--------------------Subset metadata--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
metaQiitaCombined_Nonzero_HMS <- metaQiitaCombined_Nonzero %>% filter(data_submitting_center_label == "Harvard Medical School") %>% droplevels()
metaQiitaCombined_Nonzero_BCM <- metaQiitaCombined_Nonzero %>% filter(data_submitting_center_label == "Baylor College of Medicine") %>% droplevels()
metaQiitaCombined_Nonzero_MDA <- metaQiitaCombined_Nonzero %>% filter(data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science") %>% droplevels()
metaQiitaCombined_Nonzero_WashU <- metaQiitaCombined_Nonzero %>% filter(data_submitting_center_label == "Washington University School of Medicine") %>% droplevels()
metaQiitaCombined_Nonzero_Broad_WGS <- metaQiitaCombined_Nonzero %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "WGS") %>% droplevels()

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
metaQiitaCombined_Nonzero_UNC <- metaQiitaCombined_Nonzero %>% filter(data_submitting_center_label == "University of North Carolina") %>% droplevels()
metaQiitaCombined_Nonzero_CMS <- metaQiitaCombined_Nonzero %>% filter(data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre") %>% droplevels()
metaQiitaCombined_Nonzero_Broad_RNA <- metaQiitaCombined_Nonzero %>% 
  filter(data_submitting_center_label == "Broad Institute of MIT and Harvard") %>%
  filter(experimental_strategy == "RNA-Seq") %>% droplevels()

#--------------------Subset count data--------------------#
# WGS (note that Broad has both WGS and RNA, so separate objects are made for both)
rep200_HiSeq_Fungi_Decontam_HMS <- rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero[rownames(metaQiitaCombined_Nonzero_HMS),]
rep200_HiSeq_Fungi_Decontam_BCM <- rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero[rownames(metaQiitaCombined_Nonzero_BCM),]
rep200_HiSeq_Fungi_Decontam_MDA <- rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero[rownames(metaQiitaCombined_Nonzero_MDA),]
rep200_HiSeq_Fungi_Decontam_WashU <- rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero[rownames(metaQiitaCombined_Nonzero_WashU),]
rep200_HiSeq_Fungi_Decontam_Broad_WGS <- rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero[rownames(metaQiitaCombined_Nonzero_Broad_WGS),]

# RNA-Seq (note that Broad has both WGS and RNA-Seq, so separate objects are made for both)
rep200_HiSeq_Fungi_Decontam_UNC <- rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero[rownames(metaQiitaCombined_Nonzero_UNC),]
rep200_HiSeq_Fungi_Decontam_CMS <- rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero[rownames(metaQiitaCombined_Nonzero_CMS),]
rep200_HiSeq_Fungi_Decontam_Broad_RNA <- rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero[rownames(metaQiitaCombined_Nonzero_Broad_RNA),]
#--------------------Intersect with fungi with >1% aggregate genome coverage--------------------#

coverageFungiAllSamples <- read.csv("Input_data/fungi_filt_updated_29Sep21_coverage_all_wgs_and_rna_output.csv", stringsAsFactors = FALSE)
coverageFungiAllSamples_1Percent <- coverageFungiAllSamples %>% filter(coverage_ratio >= 0.01)
coverageFungiAllSamples_1Percent_OGUs <- coverageFungiAllSamples_1Percent$gotu
# WGS
rep200_HiSeq_Fungi_Decontam_HMS_Cov <- rep200_HiSeq_Fungi_Decontam_HMS[,colnames(rep200_HiSeq_Fungi_Decontam_HMS) %in% coverageFungiAllSamples_1Percent_OGUs]
rep200_HiSeq_Fungi_Decontam_BCM_Cov <- rep200_HiSeq_Fungi_Decontam_BCM[,colnames(rep200_HiSeq_Fungi_Decontam_BCM) %in% coverageFungiAllSamples_1Percent_OGUs]
rep200_HiSeq_Fungi_Decontam_MDA_Cov <- rep200_HiSeq_Fungi_Decontam_MDA[,colnames(rep200_HiSeq_Fungi_Decontam_MDA) %in% coverageFungiAllSamples_1Percent_OGUs]
rep200_HiSeq_Fungi_Decontam_WashU_Cov <- rep200_HiSeq_Fungi_Decontam_WashU[,colnames(rep200_HiSeq_Fungi_Decontam_WashU) %in% coverageFungiAllSamples_1Percent_OGUs]
rep200_HiSeq_Fungi_Decontam_Broad_WGS_Cov <- rep200_HiSeq_Fungi_Decontam_Broad_WGS[,colnames(rep200_HiSeq_Fungi_Decontam_Broad_WGS) %in% coverageFungiAllSamples_1Percent_OGUs]
# RNA
rep200_HiSeq_Fungi_Decontam_UNC_Cov <- rep200_HiSeq_Fungi_Decontam_UNC[,colnames(rep200_HiSeq_Fungi_Decontam_UNC) %in% coverageFungiAllSamples_1Percent_OGUs]
rep200_HiSeq_Fungi_Decontam_CMS_Cov <- rep200_HiSeq_Fungi_Decontam_CMS[,colnames(rep200_HiSeq_Fungi_Decontam_CMS) %in% coverageFungiAllSamples_1Percent_OGUs]
rep200_HiSeq_Fungi_Decontam_Broad_RNA_Cov <- rep200_HiSeq_Fungi_Decontam_Broad_RNA[,colnames(rep200_HiSeq_Fungi_Decontam_Broad_RNA) %in% coverageFungiAllSamples_1Percent_OGUs]

## Remove zero sum samples -- Broad_WGS and UNC have zero sum samples; for naming convention, all will be processed
# e.g. ummary(rowSums(rep200_HiSeq_Fungi_Decontam_Broad_WGS_Cov)==0)
rep200_HiSeq_Fungi_Decontam_HMS_Cov_Nonzero <- rep200_HiSeq_Fungi_Decontam_HMS_Cov[!rowSums(rep200_HiSeq_Fungi_Decontam_HMS_Cov)==0,]
rep200_HiSeq_Fungi_Decontam_BCM_Cov_Nonzero <- rep200_HiSeq_Fungi_Decontam_BCM_Cov[!rowSums(rep200_HiSeq_Fungi_Decontam_BCM_Cov)==0,]
rep200_HiSeq_Fungi_Decontam_MDA_Cov_Nonzero <- rep200_HiSeq_Fungi_Decontam_MDA_Cov[!rowSums(rep200_HiSeq_Fungi_Decontam_MDA_Cov)==0,]
rep200_HiSeq_Fungi_Decontam_WashU_Cov_Nonzero <- rep200_HiSeq_Fungi_Decontam_WashU_Cov[!rowSums(rep200_HiSeq_Fungi_Decontam_WashU_Cov)==0,]
rep200_HiSeq_Fungi_Decontam_Broad_WGS_Cov_Nonzero <- rep200_HiSeq_Fungi_Decontam_Broad_WGS_Cov[!rowSums(rep200_HiSeq_Fungi_Decontam_Broad_WGS_Cov)==0,]
rep200_HiSeq_Fungi_Decontam_UNC_Cov_Nonzero <- rep200_HiSeq_Fungi_Decontam_UNC_Cov[!rowSums(rep200_HiSeq_Fungi_Decontam_UNC_Cov)==0,]
rep200_HiSeq_Fungi_Decontam_CMS_Cov_Nonzero <- rep200_HiSeq_Fungi_Decontam_CMS_Cov[!rowSums(rep200_HiSeq_Fungi_Decontam_CMS_Cov)==0,]
rep200_HiSeq_Fungi_Decontam_Broad_RNA_Cov_Nonzero <- rep200_HiSeq_Fungi_Decontam_Broad_RNA_Cov[!rowSums(rep200_HiSeq_Fungi_Decontam_Broad_RNA_Cov)==0,]
# Match metadata
metaQiitaCombined_Nonzero_HMS_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_HMS[rownames(rep200_HiSeq_Fungi_Decontam_HMS_Cov_Nonzero),])
metaQiitaCombined_Nonzero_BCM_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_BCM[rownames(rep200_HiSeq_Fungi_Decontam_BCM_Cov_Nonzero),])
metaQiitaCombined_Nonzero_MDA_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_MDA[rownames(rep200_HiSeq_Fungi_Decontam_MDA_Cov_Nonzero),])
metaQiitaCombined_Nonzero_WashU_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_WashU[rownames(rep200_HiSeq_Fungi_Decontam_WashU_Cov_Nonzero),])
metaQiitaCombined_Nonzero_Broad_WGS_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_Broad_WGS[rownames(rep200_HiSeq_Fungi_Decontam_Broad_WGS_Cov_Nonzero),])
metaQiitaCombined_Nonzero_UNC_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_UNC[rownames(rep200_HiSeq_Fungi_Decontam_UNC_Cov_Nonzero),])
metaQiitaCombined_Nonzero_CMS_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_CMS[rownames(rep200_HiSeq_Fungi_Decontam_CMS_Cov_Nonzero),])
metaQiitaCombined_Nonzero_Broad_RNA_Cov_Nonzero <- droplevels(metaQiitaCombined_Nonzero_Broad_RNA[rownames(rep200_HiSeq_Fungi_Decontam_Broad_RNA_Cov_Nonzero),])

#--------------------Save data for ML--------------------#

save(rep200_HiSeq_Fungi_Decontam_HMS,
     rep200_HiSeq_Fungi_Decontam_BCM,
     rep200_HiSeq_Fungi_Decontam_MDA,
     rep200_HiSeq_Fungi_Decontam_WashU,
     rep200_HiSeq_Fungi_Decontam_Broad_WGS,
     rep200_HiSeq_Fungi_Decontam_UNC,
     rep200_HiSeq_Fungi_Decontam_CMS,
     rep200_HiSeq_Fungi_Decontam_Broad_RNA,
     
     metaQiitaCombined_Nonzero_HMS,
     metaQiitaCombined_Nonzero_BCM,
     metaQiitaCombined_Nonzero_MDA,
     metaQiitaCombined_Nonzero_WashU,
     metaQiitaCombined_Nonzero_Broad_WGS,
     metaQiitaCombined_Nonzero_UNC,
     metaQiitaCombined_Nonzero_CMS,
     metaQiitaCombined_Nonzero_Broad_RNA,
     
     rep200_HiSeq_Fungi_Decontam_HMS_Cov_Nonzero,
     rep200_HiSeq_Fungi_Decontam_BCM_Cov_Nonzero,
     rep200_HiSeq_Fungi_Decontam_MDA_Cov_Nonzero,
     rep200_HiSeq_Fungi_Decontam_WashU_Cov_Nonzero,
     rep200_HiSeq_Fungi_Decontam_Broad_WGS_Cov_Nonzero,
     rep200_HiSeq_Fungi_Decontam_UNC_Cov_Nonzero,
     rep200_HiSeq_Fungi_Decontam_CMS_Cov_Nonzero,
     rep200_HiSeq_Fungi_Decontam_Broad_RNA_Cov_Nonzero,
     
     metaQiitaCombined_Nonzero_HMS_Cov_Nonzero,
     metaQiitaCombined_Nonzero_BCM_Cov_Nonzero,
     metaQiitaCombined_Nonzero_MDA_Cov_Nonzero,
     metaQiitaCombined_Nonzero_WashU_Cov_Nonzero,
     metaQiitaCombined_Nonzero_Broad_WGS_Cov_Nonzero,
     metaQiitaCombined_Nonzero_UNC_Cov_Nonzero,
     metaQiitaCombined_Nonzero_CMS_Cov_Nonzero,
     metaQiitaCombined_Nonzero_Broad_RNA_Cov_Nonzero,
     file = "Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy_with_coverage_filter_02Oct21.RData")

#--------------------Aggregate count data at taxa levels--------------------#

rep200TaxSplit_Fungi_Paired_to_Weizmann <- read.csv("Supporting_data/rep200TaxSplit_Fungi_Paired_To_Weizmann_Final.csv", stringsAsFactors = FALSE, row.names = 1)

# Build phyloseq object
psRep200_HiSeq_Fungi_Decontam_HMS <- phyloseq(otu_table(rep200_HiSeq_Fungi_Decontam_HMS, taxa_are_rows = FALSE), 
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_HMS))
psRep200_HiSeq_Fungi_Decontam_BCM <- phyloseq(otu_table(rep200_HiSeq_Fungi_Decontam_BCM, taxa_are_rows = FALSE), 
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_BCM))
psRep200_HiSeq_Fungi_Decontam_MDA <- phyloseq(otu_table(rep200_HiSeq_Fungi_Decontam_MDA, taxa_are_rows = FALSE), 
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_MDA))
psRep200_HiSeq_Fungi_Decontam_WashU <- phyloseq(otu_table(rep200_HiSeq_Fungi_Decontam_WashU, taxa_are_rows = FALSE), 
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_WashU))
psRep200_HiSeq_Fungi_Decontam_Broad_WGS <- phyloseq(otu_table(rep200_HiSeq_Fungi_Decontam_Broad_WGS, taxa_are_rows = FALSE), 
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_Broad_WGS))
psRep200_HiSeq_Fungi_Decontam_UNC <- phyloseq(otu_table(rep200_HiSeq_Fungi_Decontam_UNC, taxa_are_rows = FALSE), 
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_UNC))
psRep200_HiSeq_Fungi_Decontam_CMS <- phyloseq(otu_table(rep200_HiSeq_Fungi_Decontam_CMS, taxa_are_rows = FALSE), 
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_CMS))
psRep200_HiSeq_Fungi_Decontam_Broad_RNA <- phyloseq(otu_table(rep200_HiSeq_Fungi_Decontam_Broad_RNA, taxa_are_rows = FALSE), 
                                              tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), sample_data(metaQiitaCombined_Nonzero_Broad_RNA))

## Aggregate counts - HMS
psRep200_HiSeq_Fungi_Decontam_HMS_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_HMS, "phylum")
psRep200_HiSeq_Fungi_Decontam_HMS_class = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_HMS, "class")
psRep200_HiSeq_Fungi_Decontam_HMS_order = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_HMS, "order")
psRep200_HiSeq_Fungi_Decontam_HMS_family = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_HMS, "family")
psRep200_HiSeq_Fungi_Decontam_HMS_genus = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_HMS, "genus")
psRep200_HiSeq_Fungi_Decontam_HMS_species = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_HMS, "species")
## Aggregate counts - BCM
psRep200_HiSeq_Fungi_Decontam_BCM_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_BCM, "phylum")
psRep200_HiSeq_Fungi_Decontam_BCM_class = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_BCM, "class")
psRep200_HiSeq_Fungi_Decontam_BCM_order = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_BCM, "order")
psRep200_HiSeq_Fungi_Decontam_BCM_family = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_BCM, "family")
psRep200_HiSeq_Fungi_Decontam_BCM_genus = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_BCM, "genus")
psRep200_HiSeq_Fungi_Decontam_BCM_species = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_BCM, "species")
## Aggregate counts - MDA
psRep200_HiSeq_Fungi_Decontam_MDA_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_MDA, "phylum")
psRep200_HiSeq_Fungi_Decontam_MDA_class = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_MDA, "class")
psRep200_HiSeq_Fungi_Decontam_MDA_order = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_MDA, "order")
psRep200_HiSeq_Fungi_Decontam_MDA_family = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_MDA, "family")
psRep200_HiSeq_Fungi_Decontam_MDA_genus = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_MDA, "genus")
psRep200_HiSeq_Fungi_Decontam_MDA_species = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_MDA, "species")
## Aggregate counts - WashU
psRep200_HiSeq_Fungi_Decontam_WashU_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_WashU, "phylum")
psRep200_HiSeq_Fungi_Decontam_WashU_class = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_WashU, "class")
psRep200_HiSeq_Fungi_Decontam_WashU_order = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_WashU, "order")
psRep200_HiSeq_Fungi_Decontam_WashU_family = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_WashU, "family")
psRep200_HiSeq_Fungi_Decontam_WashU_genus = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_WashU, "genus")
psRep200_HiSeq_Fungi_Decontam_WashU_species = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_WashU, "species")
## Aggregate counts - Broad_WGS
psRep200_HiSeq_Fungi_Decontam_Broad_WGS_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_Broad_WGS, "phylum")
psRep200_HiSeq_Fungi_Decontam_Broad_WGS_class = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_Broad_WGS, "class")
psRep200_HiSeq_Fungi_Decontam_Broad_WGS_order = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_Broad_WGS, "order")
psRep200_HiSeq_Fungi_Decontam_Broad_WGS_family = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_Broad_WGS, "family")
psRep200_HiSeq_Fungi_Decontam_Broad_WGS_genus = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_Broad_WGS, "genus")
psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_Broad_WGS, "species")
## Aggregate counts - UNC
psRep200_HiSeq_Fungi_Decontam_UNC_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_UNC, "phylum")
psRep200_HiSeq_Fungi_Decontam_UNC_class = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_UNC, "class")
psRep200_HiSeq_Fungi_Decontam_UNC_order = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_UNC, "order")
psRep200_HiSeq_Fungi_Decontam_UNC_family = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_UNC, "family")
psRep200_HiSeq_Fungi_Decontam_UNC_genus = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_UNC, "genus")
psRep200_HiSeq_Fungi_Decontam_UNC_species = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_UNC, "species")
## Aggregate counts - HMS
psRep200_HiSeq_Fungi_Decontam_CMS_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_CMS, "phylum")
psRep200_HiSeq_Fungi_Decontam_CMS_class = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_CMS, "class")
psRep200_HiSeq_Fungi_Decontam_CMS_order = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_CMS, "order")
psRep200_HiSeq_Fungi_Decontam_CMS_family = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_CMS, "family")
psRep200_HiSeq_Fungi_Decontam_CMS_genus = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_CMS, "genus")
psRep200_HiSeq_Fungi_Decontam_CMS_species = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_CMS, "species")
## Aggregate counts - Broad_RNA
psRep200_HiSeq_Fungi_Decontam_Broad_RNA_phylum = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_Broad_RNA, "phylum")
psRep200_HiSeq_Fungi_Decontam_Broad_RNA_class = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_Broad_RNA, "class")
psRep200_HiSeq_Fungi_Decontam_Broad_RNA_order = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_Broad_RNA, "order")
psRep200_HiSeq_Fungi_Decontam_Broad_RNA_family = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_Broad_RNA, "family")
psRep200_HiSeq_Fungi_Decontam_Broad_RNA_genus = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_Broad_RNA, "genus")
psRep200_HiSeq_Fungi_Decontam_Broad_RNA_species = aggregate_taxa(psRep200_HiSeq_Fungi_Decontam_Broad_RNA, "species")

#--------------------Extract data frames of aggregated counts--------------------#
## HMS
df_psRep200_HiSeq_Fungi_Decontam_HMS_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_HMS_phylum)))
df_psRep200_HiSeq_Fungi_Decontam_HMS_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_HMS_class)))
df_psRep200_HiSeq_Fungi_Decontam_HMS_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_HMS_order)))
df_psRep200_HiSeq_Fungi_Decontam_HMS_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_HMS_family)))
df_psRep200_HiSeq_Fungi_Decontam_HMS_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_HMS_genus)))
df_psRep200_HiSeq_Fungi_Decontam_HMS_species <- rep200_HiSeq_Fungi_Decontam_HMS
colnames(df_psRep200_HiSeq_Fungi_Decontam_HMS_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_Decontam_HMS_species), "species"]
## BCM
df_psRep200_HiSeq_Fungi_Decontam_BCM_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_BCM_phylum)))
df_psRep200_HiSeq_Fungi_Decontam_BCM_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_BCM_class)))
df_psRep200_HiSeq_Fungi_Decontam_BCM_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_BCM_order)))
df_psRep200_HiSeq_Fungi_Decontam_BCM_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_BCM_family)))
df_psRep200_HiSeq_Fungi_Decontam_BCM_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_BCM_genus)))
df_psRep200_HiSeq_Fungi_Decontam_BCM_species <- rep200_HiSeq_Fungi_Decontam_BCM
colnames(df_psRep200_HiSeq_Fungi_Decontam_BCM_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_Decontam_BCM_species), "species"]
## MDA
df_psRep200_HiSeq_Fungi_Decontam_MDA_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_MDA_phylum)))
df_psRep200_HiSeq_Fungi_Decontam_MDA_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_MDA_class)))
df_psRep200_HiSeq_Fungi_Decontam_MDA_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_MDA_order)))
df_psRep200_HiSeq_Fungi_Decontam_MDA_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_MDA_family)))
df_psRep200_HiSeq_Fungi_Decontam_MDA_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_MDA_genus)))
df_psRep200_HiSeq_Fungi_Decontam_MDA_species <- rep200_HiSeq_Fungi_Decontam_MDA
colnames(df_psRep200_HiSeq_Fungi_Decontam_MDA_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_Decontam_MDA_species), "species"]
## WashU
df_psRep200_HiSeq_Fungi_Decontam_WashU_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_WashU_phylum)))
df_psRep200_HiSeq_Fungi_Decontam_WashU_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_WashU_class)))
df_psRep200_HiSeq_Fungi_Decontam_WashU_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_WashU_order)))
df_psRep200_HiSeq_Fungi_Decontam_WashU_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_WashU_family)))
df_psRep200_HiSeq_Fungi_Decontam_WashU_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_WashU_genus)))
df_psRep200_HiSeq_Fungi_Decontam_WashU_species <- rep200_HiSeq_Fungi_Decontam_WashU
colnames(df_psRep200_HiSeq_Fungi_Decontam_WashU_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_Decontam_WashU_species), "species"]
## Broad_WGS
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_Broad_WGS_phylum)))
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_Broad_WGS_class)))
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_Broad_WGS_order)))
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_Broad_WGS_family)))
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_Broad_WGS_genus)))
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species <- rep200_HiSeq_Fungi_Decontam_Broad_WGS
colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species), "species"]
## UNC
df_psRep200_HiSeq_Fungi_Decontam_UNC_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_UNC_phylum)))
df_psRep200_HiSeq_Fungi_Decontam_UNC_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_UNC_class)))
df_psRep200_HiSeq_Fungi_Decontam_UNC_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_UNC_order)))
df_psRep200_HiSeq_Fungi_Decontam_UNC_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_UNC_family)))
df_psRep200_HiSeq_Fungi_Decontam_UNC_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_UNC_genus)))
df_psRep200_HiSeq_Fungi_Decontam_UNC_species <- rep200_HiSeq_Fungi_Decontam_UNC
colnames(df_psRep200_HiSeq_Fungi_Decontam_UNC_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_Decontam_UNC_species), "species"]
## CMS
df_psRep200_HiSeq_Fungi_Decontam_CMS_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_CMS_phylum)))
df_psRep200_HiSeq_Fungi_Decontam_CMS_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_CMS_class)))
df_psRep200_HiSeq_Fungi_Decontam_CMS_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_CMS_order)))
df_psRep200_HiSeq_Fungi_Decontam_CMS_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_CMS_family)))
df_psRep200_HiSeq_Fungi_Decontam_CMS_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_CMS_genus)))
df_psRep200_HiSeq_Fungi_Decontam_CMS_species <- rep200_HiSeq_Fungi_Decontam_CMS
colnames(df_psRep200_HiSeq_Fungi_Decontam_CMS_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_Decontam_CMS_species), "species"]
## Broad_RNA
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_phylum <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_Broad_RNA_phylum)))
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_class <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_Broad_RNA_class)))
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_order <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_Broad_RNA_order)))
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_family <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_Broad_RNA_family)))
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_genus <- data.frame(t(otu_table(psRep200_HiSeq_Fungi_Decontam_Broad_RNA_genus)))
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_species <- rep200_HiSeq_Fungi_Decontam_Broad_RNA
colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_species) <- rep200TaxSplit_Fungi_Paired_to_Weizmann[colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_species), "species"]

#--------------------Intersect with Weizmann cohort--------------------#

## Load shared features with Weizmann
load("Interim_data/shared_fungi_features_at_each_taxa_level_13Sep21.RData")

## HMS -- subset features
df_psRep200_HiSeq_Fungi_Decontam_HMS_phylum_Shared <- df_psRep200_HiSeq_Fungi_Decontam_HMS_phylum[,colnames(df_psRep200_HiSeq_Fungi_Decontam_HMS_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_Decontam_HMS_class_Shared <- df_psRep200_HiSeq_Fungi_Decontam_HMS_class[,colnames(df_psRep200_HiSeq_Fungi_Decontam_HMS_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_Decontam_HMS_order_Shared <- df_psRep200_HiSeq_Fungi_Decontam_HMS_order[,colnames(df_psRep200_HiSeq_Fungi_Decontam_HMS_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_Decontam_HMS_family_Shared <- df_psRep200_HiSeq_Fungi_Decontam_HMS_family[,colnames(df_psRep200_HiSeq_Fungi_Decontam_HMS_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_Decontam_HMS_genus_Shared <- df_psRep200_HiSeq_Fungi_Decontam_HMS_genus[,colnames(df_psRep200_HiSeq_Fungi_Decontam_HMS_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_Decontam_HMS_species_Shared <- df_psRep200_HiSeq_Fungi_Decontam_HMS_species[,colnames(df_psRep200_HiSeq_Fungi_Decontam_HMS_species) %in% sharedSpecies]
## HMS -- remove zero sum samples
df_psRep200_HiSeq_Fungi_Decontam_HMS_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_HMS_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_HMS_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_HMS_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_HMS_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_HMS_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_HMS_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_HMS_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_HMS_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_HMS_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_HMS_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_HMS_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_HMS_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_HMS_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_HMS_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_HMS_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_HMS_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_HMS_species_Shared)==0,]
## HMS - metadata
metaQiitaCombined_Nonzero_HMS_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_HMS[rownames(df_psRep200_HiSeq_Fungi_Decontam_HMS_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_HMS_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_HMS[rownames(df_psRep200_HiSeq_Fungi_Decontam_HMS_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_HMS_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_HMS[rownames(df_psRep200_HiSeq_Fungi_Decontam_HMS_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_HMS_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_HMS[rownames(df_psRep200_HiSeq_Fungi_Decontam_HMS_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_HMS_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_HMS[rownames(df_psRep200_HiSeq_Fungi_Decontam_HMS_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_HMS_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_HMS[rownames(df_psRep200_HiSeq_Fungi_Decontam_HMS_species_Shared_Nonzero),])

## BCM -- subset features
df_psRep200_HiSeq_Fungi_Decontam_BCM_phylum_Shared <- df_psRep200_HiSeq_Fungi_Decontam_BCM_phylum[,colnames(df_psRep200_HiSeq_Fungi_Decontam_BCM_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_Decontam_BCM_class_Shared <- df_psRep200_HiSeq_Fungi_Decontam_BCM_class[,colnames(df_psRep200_HiSeq_Fungi_Decontam_BCM_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_Decontam_BCM_order_Shared <- df_psRep200_HiSeq_Fungi_Decontam_BCM_order[,colnames(df_psRep200_HiSeq_Fungi_Decontam_BCM_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_Decontam_BCM_family_Shared <- df_psRep200_HiSeq_Fungi_Decontam_BCM_family[,colnames(df_psRep200_HiSeq_Fungi_Decontam_BCM_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_Decontam_BCM_genus_Shared <- df_psRep200_HiSeq_Fungi_Decontam_BCM_genus[,colnames(df_psRep200_HiSeq_Fungi_Decontam_BCM_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_Decontam_BCM_species_Shared <- df_psRep200_HiSeq_Fungi_Decontam_BCM_species[,colnames(df_psRep200_HiSeq_Fungi_Decontam_BCM_species) %in% sharedSpecies]
## BCM -- remove zero sum samples
df_psRep200_HiSeq_Fungi_Decontam_BCM_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_BCM_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_BCM_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_BCM_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_BCM_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_BCM_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_BCM_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_BCM_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_BCM_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_BCM_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_BCM_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_BCM_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_BCM_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_BCM_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_BCM_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_BCM_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_BCM_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_BCM_species_Shared)==0,]
## BCM - metadata
metaQiitaCombined_Nonzero_BCM_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_BCM[rownames(df_psRep200_HiSeq_Fungi_Decontam_BCM_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_BCM_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_BCM[rownames(df_psRep200_HiSeq_Fungi_Decontam_BCM_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_BCM_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_BCM[rownames(df_psRep200_HiSeq_Fungi_Decontam_BCM_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_BCM_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_BCM[rownames(df_psRep200_HiSeq_Fungi_Decontam_BCM_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_BCM_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_BCM[rownames(df_psRep200_HiSeq_Fungi_Decontam_BCM_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_BCM_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_BCM[rownames(df_psRep200_HiSeq_Fungi_Decontam_BCM_species_Shared_Nonzero),])

## MDA -- subset features
df_psRep200_HiSeq_Fungi_Decontam_MDA_phylum_Shared <- df_psRep200_HiSeq_Fungi_Decontam_MDA_phylum[,colnames(df_psRep200_HiSeq_Fungi_Decontam_MDA_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_Decontam_MDA_class_Shared <- df_psRep200_HiSeq_Fungi_Decontam_MDA_class[,colnames(df_psRep200_HiSeq_Fungi_Decontam_MDA_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_Decontam_MDA_order_Shared <- df_psRep200_HiSeq_Fungi_Decontam_MDA_order[,colnames(df_psRep200_HiSeq_Fungi_Decontam_MDA_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_Decontam_MDA_family_Shared <- df_psRep200_HiSeq_Fungi_Decontam_MDA_family[,colnames(df_psRep200_HiSeq_Fungi_Decontam_MDA_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_Decontam_MDA_genus_Shared <- df_psRep200_HiSeq_Fungi_Decontam_MDA_genus[,colnames(df_psRep200_HiSeq_Fungi_Decontam_MDA_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_Decontam_MDA_species_Shared <- df_psRep200_HiSeq_Fungi_Decontam_MDA_species[,colnames(df_psRep200_HiSeq_Fungi_Decontam_MDA_species) %in% sharedSpecies]
## MDA -- remove zero sum samples
df_psRep200_HiSeq_Fungi_Decontam_MDA_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_MDA_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_MDA_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_MDA_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_MDA_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_MDA_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_MDA_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_MDA_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_MDA_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_MDA_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_MDA_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_MDA_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_MDA_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_MDA_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_MDA_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_MDA_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_MDA_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_MDA_species_Shared)==0,]
## MDA - metadata
metaQiitaCombined_Nonzero_MDA_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_MDA[rownames(df_psRep200_HiSeq_Fungi_Decontam_MDA_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_MDA_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_MDA[rownames(df_psRep200_HiSeq_Fungi_Decontam_MDA_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_MDA_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_MDA[rownames(df_psRep200_HiSeq_Fungi_Decontam_MDA_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_MDA_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_MDA[rownames(df_psRep200_HiSeq_Fungi_Decontam_MDA_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_MDA_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_MDA[rownames(df_psRep200_HiSeq_Fungi_Decontam_MDA_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_MDA_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_MDA[rownames(df_psRep200_HiSeq_Fungi_Decontam_MDA_species_Shared_Nonzero),])

## WashU -- subset features
df_psRep200_HiSeq_Fungi_Decontam_WashU_phylum_Shared <- df_psRep200_HiSeq_Fungi_Decontam_WashU_phylum[,colnames(df_psRep200_HiSeq_Fungi_Decontam_WashU_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_Decontam_WashU_class_Shared <- df_psRep200_HiSeq_Fungi_Decontam_WashU_class[,colnames(df_psRep200_HiSeq_Fungi_Decontam_WashU_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_Decontam_WashU_order_Shared <- df_psRep200_HiSeq_Fungi_Decontam_WashU_order[,colnames(df_psRep200_HiSeq_Fungi_Decontam_WashU_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_Decontam_WashU_family_Shared <- df_psRep200_HiSeq_Fungi_Decontam_WashU_family[,colnames(df_psRep200_HiSeq_Fungi_Decontam_WashU_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_Decontam_WashU_genus_Shared <- df_psRep200_HiSeq_Fungi_Decontam_WashU_genus[,colnames(df_psRep200_HiSeq_Fungi_Decontam_WashU_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_Decontam_WashU_species_Shared <- df_psRep200_HiSeq_Fungi_Decontam_WashU_species[,colnames(df_psRep200_HiSeq_Fungi_Decontam_WashU_species) %in% sharedSpecies]
## WashU -- remove zero sum samples
df_psRep200_HiSeq_Fungi_Decontam_WashU_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_WashU_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_WashU_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_WashU_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_WashU_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_WashU_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_WashU_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_WashU_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_WashU_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_WashU_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_WashU_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_WashU_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_WashU_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_WashU_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_WashU_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_WashU_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_WashU_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_WashU_species_Shared)==0,]
## WashU - metadata
metaQiitaCombined_Nonzero_WashU_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_WashU[rownames(df_psRep200_HiSeq_Fungi_Decontam_WashU_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_WashU_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_WashU[rownames(df_psRep200_HiSeq_Fungi_Decontam_WashU_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_WashU_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_WashU[rownames(df_psRep200_HiSeq_Fungi_Decontam_WashU_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_WashU_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_WashU[rownames(df_psRep200_HiSeq_Fungi_Decontam_WashU_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_WashU_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_WashU[rownames(df_psRep200_HiSeq_Fungi_Decontam_WashU_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_WashU_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_WashU[rownames(df_psRep200_HiSeq_Fungi_Decontam_WashU_species_Shared_Nonzero),])

## Broad_WGS -- subset features
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_phylum_Shared <- df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_phylum[,colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_class_Shared <- df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_class[,colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_order_Shared <- df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_order[,colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_family_Shared <- df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_family[,colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_genus_Shared <- df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_genus[,colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species_Shared <- df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species[,colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species) %in% sharedSpecies]
## Broad_WGS -- remove zero sum samples
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species_Shared)==0,]
## Broad_WGS - metadata
metaQiitaCombined_Nonzero_Broad_WGS_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_Broad_WGS_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_Broad_WGS_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_Broad_WGS_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_Broad_WGS_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_Broad_WGS_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_Broad_WGS[rownames(df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species_Shared_Nonzero),])

## UNC -- subset features
df_psRep200_HiSeq_Fungi_Decontam_UNC_phylum_Shared <- df_psRep200_HiSeq_Fungi_Decontam_UNC_phylum[,colnames(df_psRep200_HiSeq_Fungi_Decontam_UNC_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_Decontam_UNC_class_Shared <- df_psRep200_HiSeq_Fungi_Decontam_UNC_class[,colnames(df_psRep200_HiSeq_Fungi_Decontam_UNC_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_Decontam_UNC_order_Shared <- df_psRep200_HiSeq_Fungi_Decontam_UNC_order[,colnames(df_psRep200_HiSeq_Fungi_Decontam_UNC_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_Decontam_UNC_family_Shared <- df_psRep200_HiSeq_Fungi_Decontam_UNC_family[,colnames(df_psRep200_HiSeq_Fungi_Decontam_UNC_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_Decontam_UNC_genus_Shared <- df_psRep200_HiSeq_Fungi_Decontam_UNC_genus[,colnames(df_psRep200_HiSeq_Fungi_Decontam_UNC_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_Decontam_UNC_species_Shared <- df_psRep200_HiSeq_Fungi_Decontam_UNC_species[,colnames(df_psRep200_HiSeq_Fungi_Decontam_UNC_species) %in% sharedSpecies]
## UNC -- remove zero sum samples
df_psRep200_HiSeq_Fungi_Decontam_UNC_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_UNC_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_UNC_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_UNC_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_UNC_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_UNC_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_UNC_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_UNC_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_UNC_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_UNC_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_UNC_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_UNC_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_UNC_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_UNC_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_UNC_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_UNC_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_UNC_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_UNC_species_Shared)==0,]
## UNC - metadata
metaQiitaCombined_Nonzero_UNC_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_UNC[rownames(df_psRep200_HiSeq_Fungi_Decontam_UNC_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_UNC_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_UNC[rownames(df_psRep200_HiSeq_Fungi_Decontam_UNC_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_UNC_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_UNC[rownames(df_psRep200_HiSeq_Fungi_Decontam_UNC_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_UNC_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_UNC[rownames(df_psRep200_HiSeq_Fungi_Decontam_UNC_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_UNC_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_UNC[rownames(df_psRep200_HiSeq_Fungi_Decontam_UNC_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_UNC_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_UNC[rownames(df_psRep200_HiSeq_Fungi_Decontam_UNC_species_Shared_Nonzero),])

## CMS -- subset features
df_psRep200_HiSeq_Fungi_Decontam_CMS_phylum_Shared <- df_psRep200_HiSeq_Fungi_Decontam_CMS_phylum[,colnames(df_psRep200_HiSeq_Fungi_Decontam_CMS_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_Decontam_CMS_class_Shared <- df_psRep200_HiSeq_Fungi_Decontam_CMS_class[,colnames(df_psRep200_HiSeq_Fungi_Decontam_CMS_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_Decontam_CMS_order_Shared <- df_psRep200_HiSeq_Fungi_Decontam_CMS_order[,colnames(df_psRep200_HiSeq_Fungi_Decontam_CMS_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_Decontam_CMS_family_Shared <- df_psRep200_HiSeq_Fungi_Decontam_CMS_family[,colnames(df_psRep200_HiSeq_Fungi_Decontam_CMS_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_Decontam_CMS_genus_Shared <- df_psRep200_HiSeq_Fungi_Decontam_CMS_genus[,colnames(df_psRep200_HiSeq_Fungi_Decontam_CMS_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_Decontam_CMS_species_Shared <- df_psRep200_HiSeq_Fungi_Decontam_CMS_species[,colnames(df_psRep200_HiSeq_Fungi_Decontam_CMS_species) %in% sharedSpecies]
## CMS -- remove zero sum samples
df_psRep200_HiSeq_Fungi_Decontam_CMS_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_CMS_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_CMS_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_CMS_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_CMS_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_CMS_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_CMS_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_CMS_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_CMS_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_CMS_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_CMS_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_CMS_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_CMS_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_CMS_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_CMS_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_CMS_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_CMS_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_CMS_species_Shared)==0,]
## CMS - metadata
metaQiitaCombined_Nonzero_CMS_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_CMS[rownames(df_psRep200_HiSeq_Fungi_Decontam_CMS_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_CMS_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_CMS[rownames(df_psRep200_HiSeq_Fungi_Decontam_CMS_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_CMS_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_CMS[rownames(df_psRep200_HiSeq_Fungi_Decontam_CMS_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_CMS_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_CMS[rownames(df_psRep200_HiSeq_Fungi_Decontam_CMS_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_CMS_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_CMS[rownames(df_psRep200_HiSeq_Fungi_Decontam_CMS_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_CMS_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_CMS[rownames(df_psRep200_HiSeq_Fungi_Decontam_CMS_species_Shared_Nonzero),])

## Broad_RNA -- subset features
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_phylum_Shared <- df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_phylum[,colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_phylum) %in% sharedPhylum]
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_class_Shared <- df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_class[,colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_class) %in% sharedClass]
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_order_Shared <- df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_order[,colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_order) %in% sharedOrder]
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_family_Shared <- df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_family[,colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_family) %in% sharedFamily]
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_genus_Shared <- df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_genus[,colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_genus) %in% sharedGenus]
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_species_Shared <- df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_species[,colnames(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_species) %in% sharedSpecies]
## Broad_RNA -- remove zero sum samples
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_phylum_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_phylum_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_phylum_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_class_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_class_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_class_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_order_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_order_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_order_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_family_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_family_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_family_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_genus_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_genus_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_genus_Shared)==0,]
df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_species_Shared_Nonzero <- df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_species_Shared[!rowSums(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_species_Shared)==0,]
## Broad_RNA - metadata
metaQiitaCombined_Nonzero_Broad_RNA_phylum_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_phylum_Shared_Nonzero),])
metaQiitaCombined_Nonzero_Broad_RNA_class_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_class_Shared_Nonzero),])
metaQiitaCombined_Nonzero_Broad_RNA_order_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_order_Shared_Nonzero),])
metaQiitaCombined_Nonzero_Broad_RNA_family_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_family_Shared_Nonzero),])
metaQiitaCombined_Nonzero_Broad_RNA_genus_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_genus_Shared_Nonzero),])
metaQiitaCombined_Nonzero_Broad_RNA_species_Shared_Nonzero <- droplevels(metaQiitaCombined_Nonzero_Broad_RNA[rownames(df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_species_Shared_Nonzero),])

#--------------------Save data for ML--------------------#

save(df_psRep200_HiSeq_Fungi_Decontam_HMS_phylum,
     df_psRep200_HiSeq_Fungi_Decontam_HMS_class,
     df_psRep200_HiSeq_Fungi_Decontam_HMS_order,
     df_psRep200_HiSeq_Fungi_Decontam_HMS_family,
     df_psRep200_HiSeq_Fungi_Decontam_HMS_genus,
     df_psRep200_HiSeq_Fungi_Decontam_HMS_species,
     df_psRep200_HiSeq_Fungi_Decontam_HMS_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_HMS_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_HMS_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_HMS_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_HMS_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_HMS_species_Shared_Nonzero,
     # BCM
     df_psRep200_HiSeq_Fungi_Decontam_BCM_phylum,
     df_psRep200_HiSeq_Fungi_Decontam_BCM_class,
     df_psRep200_HiSeq_Fungi_Decontam_BCM_order,
     df_psRep200_HiSeq_Fungi_Decontam_BCM_family,
     df_psRep200_HiSeq_Fungi_Decontam_BCM_genus,
     df_psRep200_HiSeq_Fungi_Decontam_BCM_species,
     df_psRep200_HiSeq_Fungi_Decontam_BCM_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_BCM_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_BCM_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_BCM_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_BCM_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_BCM_species_Shared_Nonzero,
     # MDA
     df_psRep200_HiSeq_Fungi_Decontam_MDA_phylum,
     df_psRep200_HiSeq_Fungi_Decontam_MDA_class,
     df_psRep200_HiSeq_Fungi_Decontam_MDA_order,
     df_psRep200_HiSeq_Fungi_Decontam_MDA_family,
     df_psRep200_HiSeq_Fungi_Decontam_MDA_genus,
     df_psRep200_HiSeq_Fungi_Decontam_MDA_species,
     df_psRep200_HiSeq_Fungi_Decontam_MDA_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_MDA_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_MDA_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_MDA_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_MDA_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_MDA_species_Shared_Nonzero,
     # WashU
     df_psRep200_HiSeq_Fungi_Decontam_WashU_phylum,
     df_psRep200_HiSeq_Fungi_Decontam_WashU_class,
     df_psRep200_HiSeq_Fungi_Decontam_WashU_order,
     df_psRep200_HiSeq_Fungi_Decontam_WashU_family,
     df_psRep200_HiSeq_Fungi_Decontam_WashU_genus,
     df_psRep200_HiSeq_Fungi_Decontam_WashU_species,
     df_psRep200_HiSeq_Fungi_Decontam_WashU_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_WashU_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_WashU_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_WashU_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_WashU_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_WashU_species_Shared_Nonzero,
     # Broad_WGS
     df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_phylum,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_class,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_order,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_family,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_genus,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species_Shared_Nonzero,
     # UNC
     df_psRep200_HiSeq_Fungi_Decontam_UNC_phylum,
     df_psRep200_HiSeq_Fungi_Decontam_UNC_class,
     df_psRep200_HiSeq_Fungi_Decontam_UNC_order,
     df_psRep200_HiSeq_Fungi_Decontam_UNC_family,
     df_psRep200_HiSeq_Fungi_Decontam_UNC_genus,
     df_psRep200_HiSeq_Fungi_Decontam_UNC_species,
     df_psRep200_HiSeq_Fungi_Decontam_UNC_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_UNC_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_UNC_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_UNC_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_UNC_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_UNC_species_Shared_Nonzero,
     # CMS
     df_psRep200_HiSeq_Fungi_Decontam_CMS_phylum,
     df_psRep200_HiSeq_Fungi_Decontam_CMS_class,
     df_psRep200_HiSeq_Fungi_Decontam_CMS_order,
     df_psRep200_HiSeq_Fungi_Decontam_CMS_family,
     df_psRep200_HiSeq_Fungi_Decontam_CMS_genus,
     df_psRep200_HiSeq_Fungi_Decontam_CMS_species,
     df_psRep200_HiSeq_Fungi_Decontam_CMS_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_CMS_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_CMS_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_CMS_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_CMS_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_CMS_species_Shared_Nonzero,
     # Broad_RNA
     df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_phylum,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_class,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_order,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_family,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_genus,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_species,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_phylum_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_class_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_order_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_family_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_genus_Shared_Nonzero,
     df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_species_Shared_Nonzero,
     # Metadata - HMS
     metaQiitaCombined_Nonzero_HMS,
     metaQiitaCombined_Nonzero_HMS_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_HMS_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_HMS_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_HMS_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_HMS_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_HMS_species_Shared_Nonzero,
     # Metadata - BCM
     metaQiitaCombined_Nonzero_BCM,
     metaQiitaCombined_Nonzero_BCM_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_BCM_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_BCM_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_BCM_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_BCM_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_BCM_species_Shared_Nonzero,
     # Metadata - MDA
     metaQiitaCombined_Nonzero_MDA,
     metaQiitaCombined_Nonzero_MDA_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_MDA_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_MDA_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_MDA_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_MDA_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_MDA_species_Shared_Nonzero,
     # Metadata - WashU
     metaQiitaCombined_Nonzero_WashU,
     metaQiitaCombined_Nonzero_WashU_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_WashU_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_WashU_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_WashU_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_WashU_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_WashU_species_Shared_Nonzero,
     # Metadata - Broad_WGS
     metaQiitaCombined_Nonzero_Broad_WGS,
     metaQiitaCombined_Nonzero_Broad_WGS_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_Broad_WGS_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_Broad_WGS_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_Broad_WGS_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_Broad_WGS_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_Broad_WGS_species_Shared_Nonzero,
     # Metadata - UNC
     metaQiitaCombined_Nonzero_UNC,
     metaQiitaCombined_Nonzero_UNC_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_UNC_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_UNC_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_UNC_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_UNC_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_UNC_species_Shared_Nonzero,
     # Metadata - CMS
     metaQiitaCombined_Nonzero_CMS,
     metaQiitaCombined_Nonzero_CMS_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_CMS_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_CMS_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_CMS_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_CMS_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_CMS_species_Shared_Nonzero,
     # Metadata - Broad_RNA
     metaQiitaCombined_Nonzero_Broad_RNA,
     metaQiitaCombined_Nonzero_Broad_RNA_phylum_Shared_Nonzero,
     metaQiitaCombined_Nonzero_Broad_RNA_class_Shared_Nonzero,
     metaQiitaCombined_Nonzero_Broad_RNA_order_Shared_Nonzero,
     metaQiitaCombined_Nonzero_Broad_RNA_family_Shared_Nonzero,
     metaQiitaCombined_Nonzero_Broad_RNA_genus_Shared_Nonzero,
     metaQiitaCombined_Nonzero_Broad_RNA_species_Shared_Nonzero,
     file = "Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy_taxa_levels_and_wz_intersect_2Oct21.RData")

#--------------------------------------------------------------------------------------------------------------------#
# Separate WGS and RNA-Seq data for ML (to compare performance between them after batch correcting)
#--------------------------------------------------------------------------------------------------------------------#

metaQiitaCombined_Nonzero_WGS <- metaQiitaCombined_Nonzero %>% filter(experimental_strategy == "WGS") %>% droplevels()
metaQiitaCombined_Nonzero_RNA <- metaQiitaCombined_Nonzero %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()

snmDataOGUFungi_WGS <- snmDataOGUFungi[rownames(metaQiitaCombined_Nonzero_WGS),]
snmDataOGUFungi_RNA <- snmDataOGUFungi[rownames(metaQiitaCombined_Nonzero_RNA),]

save(metaQiitaCombined_Nonzero_WGS,
     metaQiitaCombined_Nonzero_RNA,
     snmDataOGUFungi_WGS,
     snmDataOGUFungi_RNA,
     file = "Interim_data/data_for_ml_tcga_wgs_vs_rna_13Sep21.RData")

# This .RData file is used to run the following scripts:
# - S06-ML-fungi-10k-rep1-tcga-wgs.R
# - S07-ML-fungi-10k-rep1-tcga-rna.R

#----------------------------------------------------------#
# Create phyloseq objects and aggregate counts at various taxa levels
#----------------------------------------------------------#

psFungiHiSeqFungi <- phyloseq(otu_table(rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero, taxa_are_rows = FALSE), 
                                tax_table(as.matrix(rep200TaxSplit)), sample_data(metaQiitaCombined_Nonzero))

## Aggregate counts
psFungiHiSeqFungi_phylum = aggregate_taxa(psFungiHiSeqFungi, "Phylum")
psFungiHiSeqFungi_class = aggregate_taxa(psFungiHiSeqFungi, "Class")
psFungiHiSeqFungi_order = aggregate_taxa(psFungiHiSeqFungi, "Order")
psFungiHiSeqFungi_family = aggregate_taxa(psFungiHiSeqFungi, "Family")
psFungiHiSeqFungi_genus = aggregate_taxa(psFungiHiSeqFungi, "Genus")
psFungiHiSeqFungi_species = aggregate_taxa(psFungiHiSeqFungi, "Species")

#----------------------------------------------------------#
# Extract aggregated and subsetted feature matrices
#----------------------------------------------------------#
rep200FungiPhylum <- data.frame(t(otu_table(psFungiHiSeqFungi_phylum)))
rep200FungiClass <- data.frame(t(otu_table(psFungiHiSeqFungi_class)))
rep200FungiOrder <- data.frame(t(otu_table(psFungiHiSeqFungi_order)))
rep200FungiFamily <- data.frame(t(otu_table(psFungiHiSeqFungi_family)))
rep200FungiGenus <- data.frame(t(otu_table(psFungiHiSeqFungi_genus)))
rep200FungiSpecies <- rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero
colnames(rep200FungiSpecies) <- rep200TaxSplit[colnames(rep200FungiSpecies), "Species"]
rep200FungiSpecies[1:3,1:3]

save(metaQiitaCombined_Nonzero,
     rep200FungiPhylum,
     rep200FungiClass,
     rep200FungiOrder,
     rep200FungiFamily,
     rep200FungiGenus,
     rep200FungiSpecies,
     file = "Interim_data/raw_data_tcga_all_cancers_16Sep21.RData")

#-----------------------------------------#
# Extract 8 cancer types to directly compare with Weizmann cancers
#-----------------------------------------#
# Cancer types to include: breast, lung, melanoma, colon, GBM, pancreas, ovary, bone
# Combine: LUAD and LUSC into LC, COAD and READ into CRC

metaQiitaCombined_Nonzero_8cancer <- metaQiitaCombined_Nonzero %>%
  filter(investigation %in% c("TCGA-BRCA","TCGA-LUAD","TCGA-LUSC","TCGA-SKCM","TCGA-COAD","TCGA-READ","TCGA-GBM","TCGA-PAAD","TCGA-OV","TCGA-SARC")) %>%
  droplevels()
dim(metaQiitaCombined_Nonzero_8cancer) # 5898   41
metaQiitaCombined_Nonzero_8cancer$disease_type <- as.character(metaQiitaCombined_Nonzero_8cancer$disease_type) 
metaQiitaCombined_Nonzero_8cancer$disease_type[metaQiitaCombined_Nonzero_8cancer$disease_type == "Lung Squamous Cell Carcinoma"] <- "Lung Cancer"
metaQiitaCombined_Nonzero_8cancer$disease_type[metaQiitaCombined_Nonzero_8cancer$disease_type == "Lung Adenocarcinoma"] <- "Lung Cancer"
metaQiitaCombined_Nonzero_8cancer$disease_type[metaQiitaCombined_Nonzero_8cancer$disease_type == "Colon Adenocarcinoma"] <- "Colorectal Cancer"
metaQiitaCombined_Nonzero_8cancer$disease_type[metaQiitaCombined_Nonzero_8cancer$disease_type == "Rectum Adenocarcinoma"] <- "Colorectal Cancer"
metaQiitaCombined_Nonzero_8cancer$disease_type[metaQiitaCombined_Nonzero_8cancer$disease_type == "Breast Invasive Carcinoma"] <- "Breast Cancer"
metaQiitaCombined_Nonzero_8cancer$disease_type[metaQiitaCombined_Nonzero_8cancer$disease_type == "Glioblastoma Multiforme"] <- "Glioblastoma"
metaQiitaCombined_Nonzero_8cancer$disease_type[metaQiitaCombined_Nonzero_8cancer$disease_type == "Pancreatic Adenocarcinoma"] <- "Pancreatic Cancer"
metaQiitaCombined_Nonzero_8cancer$disease_type[metaQiitaCombined_Nonzero_8cancer$disease_type == "Ovarian Serous Cystadenocarcinoma"] <- "Ovarian Cancer"
metaQiitaCombined_Nonzero_8cancer$disease_type[metaQiitaCombined_Nonzero_8cancer$disease_type == "Skin Cutaneous Melanoma"] <- "Melanoma"
metaQiitaCombined_Nonzero_8cancer$disease_type[metaQiitaCombined_Nonzero_8cancer$disease_type == "Sarcoma"] <- "Bone Cancer"
table(metaQiitaCombined_Nonzero_8cancer$disease_type)

rep200FungiPhylum_8cancer <- rep200FungiPhylum[rownames(metaQiitaCombined_Nonzero_8cancer),]
rep200FungiClass_8cancer <- rep200FungiClass[rownames(metaQiitaCombined_Nonzero_8cancer),]
rep200FungiOrder_8cancer <- rep200FungiOrder[rownames(metaQiitaCombined_Nonzero_8cancer),]
rep200FungiFamily_8cancer <- rep200FungiFamily[rownames(metaQiitaCombined_Nonzero_8cancer),]
rep200FungiGenus_8cancer <- rep200FungiGenus[rownames(metaQiitaCombined_Nonzero_8cancer),]
rep200FungiSpecies_8cancer <- rep200FungiSpecies[rownames(metaQiitaCombined_Nonzero_8cancer),]

save(metaQiitaCombined_Nonzero_8cancer,
     rep200FungiPhylum_8cancer,
     rep200FungiClass_8cancer,
     rep200FungiOrder_8cancer,
     rep200FungiFamily_8cancer,
     rep200FungiGenus_8cancer,
     rep200FungiSpecies_8cancer,
     file = "Interim_data/raw_data_tcga_8_cancers_16Sep21.RData")

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
rep200FungiOrderVSNM_Obj <- vsnmFunctionTCGA(rep200FungiOrder) # converges
rep200FungiOrderVSNM <- rep200FungiOrderVSNM_Obj$snmData
rep200FungiFamilyVSNM_Obj <- vsnmFunctionTCGA(rep200FungiFamily) # converges
rep200FungiFamilyVSNM <- rep200FungiFamilyVSNM_Obj$snmData
rep200FungiGenusVSNM_Obj <- vsnmFunctionTCGA(rep200FungiGenus) # converges
rep200FungiGenusVSNM <- rep200FungiGenusVSNM_Obj$snmData
rep200FungiSpeciesVSNM_Obj <- vsnmFunctionTCGA(rep200FungiSpecies) # converges
rep200FungiSpeciesVSNM <- rep200FungiSpeciesVSNM_Obj$snmData

## TCGA full data - biological variables: disease_type + sample_type | technical variables: data_submitting_center_label + experimental_strategy
# NOTE: Phylum, Class, Order, and Family levels were attempted but did not converge, so they are not shown here
rep200FungiGenusVSNM_Obj_CT <- vsnmFunctionTCGA(rep200FungiGenus, cancerTypeFlag = TRUE) # converges
rep200FungiGenusVSNM_CT <- rep200FungiGenusVSNM_Obj_CT$snmData
rep200FungiSpeciesVSNM_Obj_CT <- vsnmFunctionTCGA(rep200FungiSpecies, cancerTypeFlag = TRUE) # converges
rep200FungiSpeciesVSNM_CT <- rep200FungiSpeciesVSNM_Obj_CT$snmData

## Summarized 8 cancer types - biological variable: sample_type | technical variables: data_submitting_center_label + experimental_strategy
# NOTE: Phylum and Class levels were attempted but did not converge, so they are not shown here
rep200FungiOrder_8cancer_VSNM_Obj <- vsnmFunctionTCGA(rep200FungiOrder_8cancer, qcMetadata = metaQiitaCombined_Nonzero_8cancer) # converges
rep200FungiOrder_8cancer_VSNM <- rep200FungiOrder_8cancer_VSNM_Obj$snmData
rep200FungiFamily_8cancer_VSNM_Obj <- vsnmFunctionTCGA(rep200FungiFamily_8cancer, qcMetadata = metaQiitaCombined_Nonzero_8cancer) # converges
rep200FungiFamily_8cancer_VSNM <- rep200FungiFamily_8cancer_VSNM_Obj$snmData
rep200FungiGenus_8cancer_VSNM_Obj <- vsnmFunctionTCGA(rep200FungiGenus_8cancer, qcMetadata = metaQiitaCombined_Nonzero_8cancer) # converges
rep200FungiGenus_8cancer_VSNM <- rep200FungiGenus_8cancer_VSNM_Obj$snmData
rep200FungiSpecies_8cancer_VSNM_Obj <- vsnmFunctionTCGA(rep200FungiSpecies_8cancer, qcMetadata = metaQiitaCombined_Nonzero_8cancer) # converges
rep200FungiSpecies_8cancer_VSNM <- rep200FungiSpecies_8cancer_VSNM_Obj$snmData

## Summarized 8 cancer types - biological variables: disease_type + sample_type | technical variables: data_submitting_center_label + experimental_strategy
# NOTE: Phylum, Class, and Order levels were attempted but did not converge, so they are not shown here
rep200FungiFamily_8cancer_VSNM_Obj_CT <- vsnmFunctionTCGA(rep200FungiFamily_8cancer, qcMetadata = metaQiitaCombined_Nonzero_8cancer, cancerTypeFlag = TRUE) # converges
rep200FungiFamily_8cancer_VSNM_CT <- rep200FungiFamily_8cancer_VSNM_Obj_CT$snmData
rep200FungiGenus_8cancer_VSNM_Obj_CT <- vsnmFunctionTCGA(rep200FungiGenus_8cancer, qcMetadata = metaQiitaCombined_Nonzero_8cancer, cancerTypeFlag = TRUE) # converges
rep200FungiGenus_8cancer_VSNM_CT <- rep200FungiGenus_8cancer_VSNM_Obj_CT$snmData
rep200FungiSpecies_8cancer_VSNM_Obj_CT <- vsnmFunctionTCGA(rep200FungiSpecies_8cancer, qcMetadata = metaQiitaCombined_Nonzero_8cancer, cancerTypeFlag = TRUE) # converges
rep200FungiSpecies_8cancer_VSNM_CT <- rep200FungiSpecies_8cancer_VSNM_Obj_CT$snmData

#----------------------------------------------------------#
# Separate 8 cancer data into WGS and RNA to compare ML performance
#----------------------------------------------------------#

metaQiitaCombined_Nonzero_8cancer_WGS <- metaQiitaCombined_Nonzero_8cancer %>% filter(experimental_strategy=="WGS") %>% droplevels()
metaQiitaCombined_Nonzero_8cancer_RNA <- metaQiitaCombined_Nonzero_8cancer %>% filter(experimental_strategy=="RNA-Seq") %>% droplevels()

# Subset WGS data
rep200FungiOrder_8cancer_VSNM_WGS <- rep200FungiOrder_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_8cancer_WGS),]
rep200FungiFamily_8cancer_VSNM_WGS <- rep200FungiFamily_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_8cancer_WGS),]
rep200FungiGenus_8cancer_VSNM_WGS <- rep200FungiGenus_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_8cancer_WGS),]
rep200FungiSpecies_8cancer_VSNM_WGS <- rep200FungiSpecies_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_8cancer_WGS),]
rep200FungiFamily_8cancer_VSNM_CT_WGS <- rep200FungiFamily_8cancer_VSNM_CT[rownames(metaQiitaCombined_Nonzero_8cancer_WGS),]
rep200FungiGenus_8cancer_VSNM_CT_WGS <- rep200FungiGenus_8cancer_VSNM_CT[rownames(metaQiitaCombined_Nonzero_8cancer_WGS),]
rep200FungiSpecies_8cancer_VSNM_CT_WGS <- rep200FungiSpecies_8cancer_VSNM_CT[rownames(metaQiitaCombined_Nonzero_8cancer_WGS),]

# Subset RNA data
rep200FungiOrder_8cancer_VSNM_RNA <- rep200FungiOrder_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_8cancer_RNA),]
rep200FungiFamily_8cancer_VSNM_RNA <- rep200FungiFamily_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_8cancer_RNA),]
rep200FungiGenus_8cancer_VSNM_RNA <- rep200FungiGenus_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_8cancer_RNA),]
rep200FungiSpecies_8cancer_VSNM_RNA <- rep200FungiSpecies_8cancer_VSNM[rownames(metaQiitaCombined_Nonzero_8cancer_RNA),]
rep200FungiFamily_8cancer_VSNM_CT_RNA <- rep200FungiFamily_8cancer_VSNM_CT[rownames(metaQiitaCombined_Nonzero_8cancer_RNA),]
rep200FungiGenus_8cancer_VSNM_CT_RNA <- rep200FungiGenus_8cancer_VSNM_CT[rownames(metaQiitaCombined_Nonzero_8cancer_RNA),]
rep200FungiSpecies_8cancer_VSNM_CT_RNA <- rep200FungiSpecies_8cancer_VSNM_CT[rownames(metaQiitaCombined_Nonzero_8cancer_RNA),]

save(metaQiitaCombined_Nonzero_8cancer_WGS,
     rep200FungiOrder_8cancer_VSNM_WGS,
     rep200FungiFamily_8cancer_VSNM_WGS,
     rep200FungiGenus_8cancer_VSNM_WGS,
     rep200FungiSpecies_8cancer_VSNM_WGS,
     rep200FungiFamily_8cancer_VSNM_CT_WGS,
     rep200FungiGenus_8cancer_VSNM_CT_WGS,
     rep200FungiSpecies_8cancer_VSNM_CT_WGS,
     
     metaQiitaCombined_Nonzero_8cancer_RNA,
     rep200FungiOrder_8cancer_VSNM_RNA,
     rep200FungiFamily_8cancer_VSNM_RNA,
     rep200FungiGenus_8cancer_VSNM_RNA,
     rep200FungiSpecies_8cancer_VSNM_RNA,
     rep200FungiFamily_8cancer_VSNM_CT_RNA,
     rep200FungiGenus_8cancer_VSNM_CT_RNA,
     rep200FungiSpecies_8cancer_VSNM_CT_RNA,
     file = "Interim_data/data_for_ml_8_cancers_wgs_vs_rna_18Sep21.RData")
# Use "Interim_data/data_for_ml_8_cancers_wgs_vs_rna_18Sep21.RData" to run
# "Supporting_scripts/S13-ML-fungi-10k-rep1-tcga-8cancers-wgs-vs-rna-all-taxa-levels.R"
# Which was run on a Slurm compute cluster using the submission script
# "Supporting_scripts/S03-Submit-job.sh"
# Then download the final CSV file:
# rep_perfFungi_10k_rep1_8cancer_higher_taxa_and_intersected_ALL_13Sep21.csv
# which has been placed under the "Interim_data" folder

#----------------------------------------------------------#
# Load TCGA data matched to Weizmann by cancer types and shared features to run VSNM
#----------------------------------------------------------#
load("Interim_data/tcga_data_taxa_levels_features_shared_with_Weizmann_13Sep21.RData")
load("Interim_data/data_tcga_8_cancers_features_matched_to_weizmann_15Sep21.RData")

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
# Save data for running PVCA on each subsetted data to examine batch correction quality
#----------------------------------------------------------#

save(metaQiitaCombined_Nonzero,
     rep200FungiOrderVSNM_Obj,
     rep200FungiFamilyVSNM_Obj,
     rep200FungiGenusVSNM_Obj,
     rep200FungiSpeciesVSNM_Obj,
     rep200FungiGenusVSNM_Obj_CT,
     rep200FungiSpeciesVSNM_Obj_CT,
     # Full TCGA above, Weizmann matched cancers below
     metaQiitaCombined_Nonzero_8cancer,
     rep200FungiOrder_8cancer_VSNM_Obj,
     rep200FungiFamily_8cancer_VSNM_Obj,
     rep200FungiGenus_8cancer_VSNM_Obj,
     rep200FungiSpecies_8cancer_VSNM_Obj,
     rep200FungiFamily_8cancer_VSNM_Obj_CT,
     rep200FungiGenus_8cancer_VSNM_Obj_CT,
     rep200FungiSpecies_8cancer_VSNM_Obj_CT,
     file = "Interim_data/data_for_pvca_tcga_taxa_levels_15Sep21.RData")
# Use "Interim_data/data_for_pvca_tcga_taxa_levels_15Sep21.RData" to run
# "Supporting_scripts/S08-Run-pvca-all-taxa-levels.R"
# Which was run on a Slurm compute cluster using the submission script
# "Supporting_scripts/S08-Submit-job.sh"
# Figures are generated based on this script for each dataset tested,
# and the PVCA results for each dataset are saved as .RData files

#----------------------------------------------------------#
# Save data for running machine learning
#----------------------------------------------------------#

save(metaQiitaCombined_Nonzero_8cancer,
     rep200FungiOrder_8cancer_VSNM,
     rep200FungiFamily_8cancer_VSNM,
     rep200FungiGenus_8cancer_VSNM,
     rep200FungiSpecies_8cancer_VSNM,
     rep200FungiFamily_8cancer_VSNM_CT,
     rep200FungiGenus_8cancer_VSNM_CT,
     rep200FungiSpecies_8cancer_VSNM_CT,
     file = "Interim_data/data_for_ml_8_cancers_13Sep21.RData")
# Use "Interim_data/data_for_ml_8_cancers_13Sep21.RData" to run
# "Supporting_scripts/S03-ML-fungi-10k-rep1-higher-taxa-8cancer.R"
# Which was run on a Slurm compute cluster using the submission script
# "Supporting_scripts/S03-Submit-job.sh"
# Then download the final CSV file:
# rep_perfFungi_10k_rep1_8cancer_higher_taxa_and_intersected_ALL_13Sep21.csv
# which has been placed under the "Interim_data" folder

save(metaQiitaCombined_Nonzero,
     rep200FungiOrderVSNM,
     rep200FungiFamilyVSNM,
     rep200FungiGenusVSNM,
     rep200FungiSpeciesVSNM,
     rep200FungiGenusVSNM_CT,
     rep200FungiSpeciesVSNM_CT,
     snmDataOGUFungi,
     file = "Interim_data/data_for_ml_tcga_13Sep21.RData")
# Use "Interim_data/data_for_ml_8_cancers_13Sep21.RData" to run
# "Supporting_scripts/S04-ML-fungi-10k-rep1-higher-taxa-tcga.R"
# Which was run on a Slurm compute cluster using the submission script
# "Supporting_scripts/S04-Submit-job.sh"
# Then download the final CSV file:
# rep_perfFungi_10k_rep1_tcga_higher_taxa_and_intersected_ALL_13Sep21.csv
# which has been placed under the "Interim_data" folder

save(rep200FungiFamilyShared_NonzeroVSNM,
     metaQiitaCombined_Nonzero_FamilyShared,
     rep200FungiGenusShared_NonzeroVSNM,
     metaQiitaCombined_Nonzero_GenusShared,
     rep200FungiSpeciesShared_NonzeroVSNM,
     metaQiitaCombined_Nonzero_SpeciesShared,
     file = "Interim_data/data_for_ml_tcga_shared_features_15Sep21.RData")
# Use "Interim_data/data_for_ml_tcga_shared_features_15Sep21.RData" to run
# "Supporting_scripts/S09-ML-fungi-10k-rep1-tcga-shared-features.R"
# Which was run on a Slurm compute cluster using the submission script
# "Supporting_scripts/S09-Submit-job.sh"
# Then download the final CSV file
# which has been placed under the "Interim_data" folder


save(rep200FungiPhylumShared_8cancer_Nonzero_VSNM,
     metaQiitaCombined_Nonzero_8cancer_PhylumShared,
     rep200FungiFamilyShared_8cancer_Nonzero_VSNM,
     metaQiitaCombined_Nonzero_8cancer_FamilyShared,
     rep200FungiGenusShared_8cancer_Nonzero_VSNM,
     metaQiitaCombined_Nonzero_8cancer_GenusShared,
     rep200FungiSpeciesShared_8cancer_Nonzero_VSNM,
     metaQiitaCombined_Nonzero_8cancer_SpeciesShared,
     file = "Interim_data/data_for_ml_8_cancers_shared_features_15Sep21.RData")
# Use "Interim_data/data_for_ml_8_cancers_shared_features_15Sep21.RData" to run
# "Supporting_scripts/S10-ML-fungi-10k-rep1-tcga-8cancers-shared-features.R"
# Which was run on a Slurm compute cluster using the submission script
# "Supporting_scripts/S10-Submit-job.sh"
# Then download the final CSV file
# which has been placed under the "Interim_data" folder

#----------------------------------------------------------#
# Plot machine learning performances for 8 cancers
#----------------------------------------------------------#

source("Supporting_scripts/S05-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_8cancer <- read.csv("Interim_data/rep_perfFungi_10k_rep1_8cancer_higher_taxa_and_intersected_ALL_13Sep21.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_8cancer <- read.csv("Supporting_data/tcga_abbreviations_8cancer.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_8cancer$abbrev <- abbreviationsTCGA_8cancer[mlPerfAll10k_8cancer$diseaseType,"abbrev"]
mlPerfAll10k_8cancer <- mlPerfAll10k_8cancer[,!(colnames(mlPerfAll10k_8cancer) == "X")]
colnames(mlPerfAll10k_8cancer)[1:2] <- c("AUROC","AUPR")

# Rename entries in the "datasetName" column
# NOTE: These entries were originally listed in the "datasetListNames" object
# in the "Supporting_scripts/S03-ML-fungi-10k-rep1-higher-taxa-8cancer.R" script

mlPerfAll10k_8cancer$datasetName[mlPerfAll10k_8cancer$datasetName == "rep200FungiOrder_8cancer_VSNM"] <- "Order (Decontam)"
mlPerfAll10k_8cancer$datasetName[mlPerfAll10k_8cancer$datasetName == "rep200FungiFamily_8cancer_VSNM"] <- "Family (Decontam)"
mlPerfAll10k_8cancer$datasetName[mlPerfAll10k_8cancer$datasetName == "rep200FungiGenus_8cancer_VSNM"] <- "Genus (Decontam)"
mlPerfAll10k_8cancer$datasetName[mlPerfAll10k_8cancer$datasetName == "rep200FungiSpecies_8cancer_VSNM"] <- "Species (Decontam)"

mlPerfAll10k_8cancer$datasetName[mlPerfAll10k_8cancer$datasetName == "rep200FungiFamily_8cancer_VSNM_CT"] <- "Family (CT)"
mlPerfAll10k_8cancer$datasetName[mlPerfAll10k_8cancer$datasetName == "rep200FungiGenus_8cancer_VSNM_CT"] <- "Genus (CT)"
mlPerfAll10k_8cancer$datasetName[mlPerfAll10k_8cancer$datasetName == "rep200FungiSpecies_8cancer_VSNM_CT"] <- "Species (CT)"

mlPerfAll10k_8cancer$datasetName <- factor(mlPerfAll10k_8cancer$datasetName,
                                           levels = c("Class (Decontam)", 
                                                      "Order (Decontam)",
                                                      "Family (Decontam)",
                                                      "Genus (Decontam)",
                                                      "Species (Decontam)",
                                                      "Family (CT)",
                                                      "Genus (CT)",
                                                      "Species (CT)"))
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#

# All taxa levels and their intersections
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor | 1 Vs All | Order, Family, Genus, Species Levels") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
ggsave("Figures/Figure_4/figure_4_D__mlPerfAll10k_rep1_8cancers_PT_order_family_genus_species.jpeg", dpi = "retina",
       width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>% write.csv("Figures_data/Figure_4/figure_4_D__mlPerfAll10k_rep1_8cancers_PT_order_family_genus_species.csv")

# All taxa levels and their intersections - including CT Voom-SNM
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor") %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor | 1 Vs All | Order, Family, Genus, Species Levels") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_D__optional_mlPerfAll10k_rep1_8cancers_PT_order_family_genus_species.jpeg", dpi = "retina",
         width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>% write.csv("Figures_data/Figure_4/figure_4_D__optional_mlPerfAll10k_rep1_8cancers_PT_order_family_genus_species.csv")

# Order level
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Order",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor | 1 Vs All | Order") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") +
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_order.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")

# Family level
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Family",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor | 1 Vs All | Family") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") +
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_family.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")

# Genus level
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Genus",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor | 1 Vs All | Genus") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") +
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_genus.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")

# Species level
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") +
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")

#-------------------------Plot primary tumor vs. adjacent tissue normal performance-------------------------#
# All taxa levels and their intersections
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs. Adjacent Normal | Order, Family, Genus, Species Levels") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
ggsave("Figures/Figure_4/figure_4_E__mlPerfAll10k_rep1_8cancers_PT_vs_STN_order_family_genus_species.jpeg", dpi = "retina",
       width = 8, height = 5, units = "in")
# Save underlying data
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  droplevels() %>% write.csv("Figures_data/Figure_4/figure_4_E__mlPerfAll10k_rep1_8cancers_PT_vs_STN_order_family_genus_species.csv")

# All taxa levels and their intersections - including CT Voom-SNM
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs. Adjacent Normal | Order, Family, Genus, Species Levels") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_E__optional_mlPerfAll10k_rep1_8cancers_PT_vs_STN_order_family_genus_species.jpeg", dpi = "retina",
         width = 8, height = 5, units = "in")
# Save underlying data
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  droplevels() %>% write.csv("Figures_data/Figure_4/figure_4_E__optional_mlPerfAll10k_rep1_8cancers_PT_vs_STN_order_family_genus_species.csv")

# Order level
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Order",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs. Adjacent Normal | Order") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") +
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_vs_STN_order.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")

# Family level
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Family",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs. Adjacent Normal | Family") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") +
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_vs_STN_family.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")

# Genus level
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Genus",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs. Adjacent Normal | Genus") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") +
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_vs_STN_genus.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")

# Species level
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs. Adjacent Normal | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") +
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_vs_STN_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")


#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
# All taxa levels and their intersections
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Blood Derived Normal | 1 Vs All | Order, Family, Genus, Species Levels") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
ggsave("Figures/Figure_4/figure_4_F__mlPerfAll10k_rep1_8cancers_BDN_order_family_genus_species.jpeg", dpi = "retina",
       width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  droplevels() %>% write.csv("Figures_data/Figure_4/mlPerfAll10k_rep1_8cancers_BDN_order_family_genus_species.csv")

# All taxa levels and their intersections - including CT Voom-SNM
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Blood Derived Normal | 1 Vs All | Order, Family, Genus, Species Levels") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_F__optional_mlPerfAll10k_rep1_8cancers_BDN_order_family_genus_species.jpeg", dpi = "retina",
         width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  droplevels() %>% write.csv("Figures_data/Figure_4/optional_mlPerfAll10k_rep1_8cancers_BDN_order_family_genus_species.csv")

# Order level
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Order",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Blood Derived Normal | 1 Vs All | Order") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") +
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_BDN_order.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")

# Family level
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Family",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Blood Derived Normal | 1 Vs All | Family") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") +
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_BDN_family.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")

# Genus level
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Genus",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Blood Derived Normal | 1 Vs All | Genus") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") +
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_BDN_genus.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")

# Species level
mlPerfAll10k_8cancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") +
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_BDN_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")

#----------------------------------------------------------#
# Plot machine learning performances for 8 cancers - WGS vs RNA
#----------------------------------------------------------#

source("Supporting_scripts/S05-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_8cancer_WGS_RNA <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_8cancers_wgs_vs_rna_taxa_levels_ALL_18Sep21.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_8cancer <- read.csv("Supporting_data/tcga_abbreviations_8cancer.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_8cancer_WGS_RNA$abbrev <- abbreviationsTCGA_8cancer[mlPerfAll10k_8cancer_WGS_RNA$diseaseType,"abbrev"]
mlPerfAll10k_8cancer_WGS_RNA <- mlPerfAll10k_8cancer_WGS_RNA[,!(colnames(mlPerfAll10k_8cancer_WGS_RNA) == "X")]
colnames(mlPerfAll10k_8cancer_WGS_RNA)[1:2] <- c("AUROC","AUPR")

# Rename entries in the "datasetName" column

mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiOrder_8cancer_VSNM_WGS"] <- "Order WGS (Decontam)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiOrder_8cancer_VSNM_RNA"] <- "Order RNA (Decontam)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiFamily_8cancer_VSNM_WGS"] <- "Family WGS (Decontam)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiFamily_8cancer_VSNM_RNA"] <- "Family RNA (Decontam)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiGenus_8cancer_VSNM_WGS"] <- "Genus WGS (Decontam)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiGenus_8cancer_VSNM_RNA"] <- "Genus RNA (Decontam)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiSpecies_8cancer_VSNM_WGS"] <- "Species WGS (Decontam)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiSpecies_8cancer_VSNM_RNA"] <- "Species RNA (Decontam)"

mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiFamily_8cancer_VSNM_CT_WGS"] <- "Family WGS (CT)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiFamily_8cancer_VSNM_CT_RNA"] <- "Family RNA (CT)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiGenus_8cancer_VSNM_CT_WGS"] <- "Genus WGS (CT)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiGenus_8cancer_VSNM_CT_RNA"] <- "Genus RNA (CT)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiSpecies_8cancer_VSNM_CT_WGS"] <- "Species WGS (CT)"
mlPerfAll10k_8cancer_WGS_RNA$datasetName[mlPerfAll10k_8cancer_WGS_RNA$datasetName == "rep200FungiSpecies_8cancer_VSNM_CT_RNA"] <- "Species RNA (CT)"

mlPerfAll10k_8cancer_WGS_RNA$datasetName <- factor(mlPerfAll10k_8cancer_WGS_RNA$datasetName,
                                           levels = c("Order WGS (Decontam)",
                                                      "Order RNA (Decontam)",
                                                      "Family WGS (Decontam)",
                                                      "Family RNA (Decontam)",
                                                      "Genus WGS (Decontam)",
                                                      "Genus RNA (Decontam)",
                                                      "Species WGS (Decontam)",
                                                      "Species RNA (Decontam)",
                                                      "Family WGS (CT)",
                                                      "Family RNA (CT)",
                                                      "Genus WGS (CT)",
                                                      "Genus RNA (CT)",
                                                      "Species WGS (CT)",
                                                      "Species RNA (CT)"))
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#

mlPerfAll10k_8cancer_WGS_RNA %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","metadataName")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor | 1 Vs All | Species | WGS vs. RNA-Seq") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") #+ 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_species_wgs_vs_rna.jpeg", dpi = "retina",
         width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_8cancer_WGS_RNA %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>% write.csv("Figures_data/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_species_wgs_vs_rna.csv")

#-------------------------Plot primary tumor vs. normal performance-------------------------#

mlPerfAll10k_8cancer_WGS_RNA %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","metadataName")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs. NAT | Species | WGS vs. RNA-Seq") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_vs_NAT_species_wgs_vs_rna.jpeg", dpi = "retina",
       width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_8cancer_WGS_RNA %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>% write.csv("Figures_data/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_vs_NAT_species_wgs_vs_rna.jpeg")

#----------------------------------------------------------#
# Plot machine learning performances for 8 cancers - shared features
#----------------------------------------------------------#

source("Supporting_scripts/S05-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_8cancer_Shared <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_8cancers_shared_features_ALL_15Sep21.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_8cancer <- read.csv("Supporting_data/tcga_abbreviations_8cancer.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_8cancer_Shared$abbrev <- abbreviationsTCGA_8cancer[mlPerfAll10k_8cancer_Shared$diseaseType,"abbrev"]
mlPerfAll10k_8cancer_Shared <- mlPerfAll10k_8cancer_Shared[,!(colnames(mlPerfAll10k_8cancer_Shared) == "X")]
colnames(mlPerfAll10k_8cancer_Shared)[1:2] <- c("AUROC","AUPR")

# Rename entries in the "datasetName" column
mlPerfAll10k_8cancer_Shared$datasetName[mlPerfAll10k_8cancer_Shared$datasetName == "rep200FungiPhylumShared_8cancer_Nonzero_VSNM"] <- "Phylum intersected with Weizmann"
mlPerfAll10k_8cancer_Shared$datasetName[mlPerfAll10k_8cancer_Shared$datasetName == "rep200FungiFamilyShared_8cancer_Nonzero_VSNM"] <- "Family intersected with Weizmann"
mlPerfAll10k_8cancer_Shared$datasetName[mlPerfAll10k_8cancer_Shared$datasetName == "rep200FungiGenusShared_8cancer_Nonzero_VSNM"] <- "Genus intersected with Weizmann"
mlPerfAll10k_8cancer_Shared$datasetName[mlPerfAll10k_8cancer_Shared$datasetName == "rep200FungiSpeciesShared_8cancer_Nonzero_VSNM"] <- "Species intersected with Weizmann"

mlPerfAll10k_8cancer_Shared$datasetName <- factor(mlPerfAll10k_8cancer_Shared$datasetName,
                                                   levels = c("Phylum intersected with Weizmann",
                                                              "Family intersected with Weizmann",
                                                              "Genus intersected with Weizmann",
                                                              "Species intersected with Weizmann"))
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#

mlPerfAll10k_8cancer_Shared %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","metadataName")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor | 1 Vs All | Species intersected with Weizmann") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_species_shared.jpeg", dpi = "retina",
       width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_8cancer_WGS_RNA %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  droplevels() %>% write.csv("Figures_data/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_species_shared.csv")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#

mlPerfAll10k_8cancer_Shared %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  # filter(grepl("Species",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","metadataName")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs NAT | Phylum, Family, Genus, & Species intersected with Weizmann") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_vs_NAT_phylum_family_genus_species_shared.jpeg", dpi = "retina",
         width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_8cancer_WGS_RNA %>%
  filter(sampleType == "Primary Tumor") %>%
  # filter(grepl("Species",datasetName)) %>%
  droplevels() %>% write.csv("Figures_data/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_PT_vs_NAT_phylum_family_genus_species_shared.csv")

#-------------------------Plot blood derived normal 1 vs all others performance-------------------------#

mlPerfAll10k_8cancer_Shared %>%
  filter(sampleType == "Blood Derived Normal") %>%
  # filter(grepl("Species",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","metadataName")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Blood Derived Normal | 1 vs All Others | Phylum, Family, Genus, & Species intersected with Weizmann") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_BDN_phylum_family_genus_species_shared.jpeg", dpi = "retina",
         width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_8cancer_WGS_RNA %>%
  filter(sampleType == "Blood Derived Normal") %>%
  # filter(grepl("Species",datasetName)) %>%
  droplevels() %>% write.csv("Figures_data/Supplementary_Figures/mlPerfAll10k_rep1_8cancers_BDN_phylum_family_genus_species_shared.csv")

#----------------------------------------------------------#
# Overlay plot machine learning performances for for 8 cancers:
# - Full features
# - WGS vs RNA
# - Shared features
#----------------------------------------------------------#

# Make sure to run the above 3 sections above to correctly format the following objects:
# mlPerfAll10k_8cancer, mlPerfAll10k_8cancer_WGS_RNA, mlPerfAll10k_8cancer_Shared
# Note that the mlPerfAll10k_8cancer is missing the "metadataName" column, so one is added prior to rbind
mlPerfAll10k_8cancer_Overlay <- rbind(cbind(mlPerfAll10k_8cancer, metadataName="metaQiitaCombined_Nonzero_8cancer"),
                                      mlPerfAll10k_8cancer_WGS_RNA,
                                      mlPerfAll10k_8cancer_Shared)

#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#

mlPerfAll10k_8cancer_Overlay %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","metadataName")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor | 1 Vs All | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
ggsave("Figures/Figure_4/figure_4_D__overlay_mlPerfAll10k_rep1_8cancers_PT_species.jpeg", dpi = "retina",
       width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_8cancer_Overlay %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>% write.csv("Figures_data/Figure_4/figure_4_D__overlay_mlPerfAll10k_rep1_8cancers_PT_species.csv")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#

mlPerfAll10k_8cancer_Overlay %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","metadataName")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs NAT | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_E__overlay_mlPerfAll10k_rep1_8cancers_PT_vs_NAT_species.jpeg", dpi = "retina",
         width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_8cancer_Overlay %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>% write.csv("Figures_data/Figure_4/figure_4_E__overlay_mlPerfAll10k_rep1_8cancers_PT_vs_NAT_species.csv")

#-------------------------Plot BDN 1 vs. all others performance-------------------------#

mlPerfAll10k_8cancer_Overlay %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","metadataName")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Blood Derived Normal | 1 Vs All | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_F__overlay_mlPerfAll10k_rep1_8cancers_BDN_species.jpeg", dpi = "retina",
         width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_8cancer_Overlay %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>% write.csv("Figures_data/Figure_4/figure_4_F__overlay_mlPerfAll10k_rep1_8cancers_BDN_species.csv")

#----------------------------------------------------------#
# Plot machine learning performances for all TCGA cancers
#----------------------------------------------------------#

source("Supporting_scripts/S05-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_higher_taxa_and_intersected_ALL_13Sep21.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer <- mlPerfAll10k_Allcancer[,!(colnames(mlPerfAll10k_Allcancer) == "X")]
colnames(mlPerfAll10k_Allcancer)[1:2] <- c("AUROC","AUPR")

# Rename entries in the "datasetName" column
# NOTE: These entries were originally listed in the "datasetListNames" object
# in the "Supporting_scripts/S04-ML-fungi-10k-rep1-higher-taxa-tcga.R" script

mlPerfAll10k_Allcancer$datasetName[mlPerfAll10k_Allcancer$datasetName == "rep200FungiOrderVSNM"] <- "Order (Decontam)"
mlPerfAll10k_Allcancer$datasetName[mlPerfAll10k_Allcancer$datasetName == "rep200FungiFamilyVSNM"] <- "Family (Decontam)"
mlPerfAll10k_Allcancer$datasetName[mlPerfAll10k_Allcancer$datasetName == "rep200FungiGenusVSNM"] <- "Genus (Decontam)"
mlPerfAll10k_Allcancer$datasetName[mlPerfAll10k_Allcancer$datasetName == "rep200FungiSpeciesVSNM"] <- "Species (Decontam)"

mlPerfAll10k_Allcancer$datasetName[mlPerfAll10k_Allcancer$datasetName == "rep200FungiGenusVSNM_CT"] <- "Genus (CT)"
mlPerfAll10k_Allcancer$datasetName[mlPerfAll10k_Allcancer$datasetName == "rep200FungiSpeciesVSNM_CT"] <- "Species (CT)"
mlPerfAll10k_Allcancer$datasetName[mlPerfAll10k_Allcancer$datasetName == "snmDataOGUFungi"] <- "OGU (CT)"

mlPerfAll10k_Allcancer$datasetName <- factor(mlPerfAll10k_Allcancer$datasetName,
                                           levels = c("Order (Decontam)",
                                                      "Family (Decontam)",
                                                      "Genus (Decontam)",
                                                      "Species (Decontam)",
                                                      "Genus (CT)",
                                                      "Species (CT)",
                                                      "OGU (CT)"))
#-------------------------Regress avg performance vs minority class size-------------------------#
ptPerfSummarizedMinClassSize <- mlPerfAll10k_Allcancer %>% filter(grepl("Species",datasetName)) %>% 
                              filter(grepl("Decontam",datasetName)) %>%
                              filter(sampleType == "Primary Tumor") %>%
                              group_by(sampleType, diseaseType) %>% 
                              mutate(avgROC = mean(AUROC), avgPR = mean(AUPR), avgMinClass = mean(minorityClassSize)) %>%
                              select(sampleType, diseaseType, abbrev, avgROC, avgPR, avgMinClass) %>% 
                              mutate(logAvgMinClass = log10(avgMinClass)) %>%
                              unique() %>% data.frame()

summary(lm(avgPR ~ avgMinClass, ptPerfSummarizedMinClassSize))
summary(lm(avgROC ~ avgMinClass, ptPerfSummarizedMinClassSize))

ptPerfSummarizedMinClassSize %>%
  ggplot(aes(x=avgMinClass, y=avgPR)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE) + 
  theme_bw() + theme(aspect.ratio=1) + coord_fixed() #+
  ggsave(filename = "Figures/Figure_1/corr_fungal_and_bacterial_read_percentages_sample_type.jpeg",
         dpi = "retina", units = "in", width = 7, height = 7)

ptPerfSummarizedMinClassSize %>%
  ggplot(aes(x=logAvgMinClass, y=avgROC)) +
  geom_point(alpha = 0.4) + geom_smooth(method='lm') + stat_cor(method = "pearson", cor.coef.name = "R", show.legend = FALSE) + 
  theme_bw()


#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
# Order level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Order",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor | 1 Vs All | Order") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_order.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

# Family level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Family",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor | 1 Vs All | Family") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_family.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

# Genus level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Genus",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor | 1 Vs All | Genus") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_genus.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

# Species level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_species.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

#-------------------------Plot primary tumor vs. adjacent normal tissue performance-------------------------#
# Order level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Order",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor vs. Adjacent Normal | Order") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_STN_order.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

# Family level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Family",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor vs. Adjacent Normal | Family") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_STN_family.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

# Genus level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Genus",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor vs. Adjacent Normal | Genus") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") +
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_STN_genus.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

# Species level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor vs. Adjacent Normal | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_STN_species.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
# Order level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Order",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Blood Derived Normal | 1 Vs All | Order") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_BDN_order.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

# Family level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Family",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Blood Derived Normal | 1 Vs All | Family") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_BDN_family.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

# Genus level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Genus",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Blood Derived Normal | 1 Vs All | Genus") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_BDN_genus.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

# Species level
mlPerfAll10k_Allcancer %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_BDN_species.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

#----------------------------------------------------------#
# Plot machine learning performances for all TCGA cancers - shared features
#----------------------------------------------------------#

source("Supporting_scripts/S05-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Shared <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_shared_features_ALL_15Sep21.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_Shared$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_Shared$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_Shared <- mlPerfAll10k_Allcancer_Shared[,!(colnames(mlPerfAll10k_Allcancer_Shared) == "X")]
colnames(mlPerfAll10k_Allcancer_Shared)[1:2] <- c("AUROC","AUPR")

# Rename entries in the "datasetName" column

mlPerfAll10k_Allcancer_Shared$datasetName[mlPerfAll10k_Allcancer_Shared$datasetName == "rep200FungiFamilyShared_NonzeroVSNM"] <- "Family intersected with Weizmann"
mlPerfAll10k_Allcancer_Shared$datasetName[mlPerfAll10k_Allcancer_Shared$datasetName == "rep200FungiGenusShared_NonzeroVSNM"] <- "Genus intersected with Weizmann"
mlPerfAll10k_Allcancer_Shared$datasetName[mlPerfAll10k_Allcancer_Shared$datasetName == "rep200FungiSpeciesShared_NonzeroVSNM"] <- "Species intersected with Weizmann"
mlPerfAll10k_Allcancer_Shared$datasetName <- factor(mlPerfAll10k_Allcancer_Shared$datasetName,
                                             levels = c("Family intersected with Weizmann",
                                                        "Genus intersected with Weizmann",
                                                        "Species intersected with Weizmann"))
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
mlPerfAll10k_Allcancer_Shared %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor | 1 Vs All | Species intersected with Weizmann") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_species_shared.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#
mlPerfAll10k_Allcancer_Shared %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor vs NAT | Species intersected with Weizmann") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_NAT_species_shared.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
mlPerfAll10k_Allcancer_Shared %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=variable)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  # facet_wrap(~variable) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Blood Derived Normal | 1 Vs All | Species intersected with Weizmann") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_BDN_species_shared.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

#----------------------------------------------------------#
# Plot machine learning performances for WGS vs RNA-Seq
#----------------------------------------------------------#

source("Supporting_scripts/S05-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_WGS <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_wgs_ALL_13Sep21.csv", stringsAsFactors = FALSE)
mlPerfAll10k_RNA <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_rna_ALL_13Sep21.csv", stringsAsFactors = FALSE)
mlPerfAll10k_WGS_RNA <- rbind(mlPerfAll10k_WGS, mlPerfAll10k_RNA)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_WGS_RNA$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_WGS_RNA$diseaseType,"abbrev"]
mlPerfAll10k_WGS_RNA <- mlPerfAll10k_WGS_RNA[,!(colnames(mlPerfAll10k_WGS_RNA) == "X")]
colnames(mlPerfAll10k_WGS_RNA)[1:2] <- c("AUROC","AUPR")

# Rename entries in the "datasetName" column
# NOTE: These entries were originally listed in the "datasetListNames" object
# in the "Supporting_scripts/S04-ML-fungi-10k-rep1-higher-taxa-tcga.R" script

mlPerfAll10k_WGS_RNA$datasetName[mlPerfAll10k_WGS_RNA$datasetName == "snmDataOGUFungi_WGS"] <- "WGS Species (Decontam)"
mlPerfAll10k_WGS_RNA$datasetName[mlPerfAll10k_WGS_RNA$datasetName == "snmDataOGUFungi_RNA"] <- "RNA-Seq Species (Decontam)"

mlPerfAll10k_WGS_RNA$datasetName <- factor(mlPerfAll10k_WGS_RNA$datasetName,
                                             levels = c( "WGS Species (Decontam)",
                                                         "RNA-Seq Species (Decontam)"))
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
mlPerfAll10k_WGS_RNA %>%
  filter(sampleType == "Primary Tumor") %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor | 1 Vs All | OGUs | WGS vs. RNA") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_wgs_vs_rna.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

#-------------------------Plot primary tumor vs. adjacent normal performance-------------------------#
mlPerfAll10k_WGS_RNA %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All of TCGA | Primary Tumor vs Solid Tissue Normal | 1 Vs All | OGUs | WGS vs. RNA") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Supplementary_Figures/mlPerfAll10k_rep1_Allcancers_PT_vs_STN_wgs_vs_rna.jpeg", dpi = "retina",
         width = 12, height = 6, units = "in")

#----------------------------------------------------------------------------------------------------#
## NOTE: Only WGS was run on Blood Derived Normal TCGA samples, so there is no comparison between them
#----------------------------------------------------------------------------------------------------#

#----------------------------------------------------------#
# Overlay plots machine learning performances for all TCGA cancers:
# - VSNM data
# - WGS vs. RNA
# - Shared features
#----------------------------------------------------------#

source("Supporting_scripts/S05-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

# NOTE: There are too many cancer types to include more than two errorbars per cancer type,
# so the WGS vs. RNA comparison is not included in this overlay. If you would like to include
# it, you can remove the commented line below.
mlPerfAll10k_Allcancer_Overlay <- rbind(cbind(mlPerfAll10k_Allcancer, metadataName=NA),
                                        # cbind(mlPerfAll10k_WGS_RNA, metadataName=NA),
                                        mlPerfAll10k_Allcancer_Shared)

#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#

mlPerfAll10k_Allcancer_Overlay %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","metadataName")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("All TCGA Cancers | Primary Tumor | 1 Vs All | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_D__alternative_overlay_mlPerfAll10k_rep1_Allcancers_PT_species.jpeg", dpi = "retina",
         width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_Allcancer_Overlay %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>% write.csv("Figures_data/Figure_4/figure_4_D__alternative_overlay_mlPerfAll10k_rep1_Allcancers_PT_species.csv")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#

mlPerfAll10k_Allcancer_Overlay %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","metadataName")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Primary Tumor vs NAT | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_E__alternative_overlay_mlPerfAll10k_rep1_Allcancers_PT_vs_NAT_species.jpeg", dpi = "retina",
         width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_Allcancer_Overlay %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>% write.csv("Figures_data/Figure_4/figure_4_E__alternative_overlay_mlPerfAll10k_rep1_Allcancers_PT_vs_NAT_species.csv")

#-------------------------Plot BDN 1 vs. all others performance-------------------------#

mlPerfAll10k_Allcancer_Overlay %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","variable","abbrev","metadataName")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) + scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Blood Derived Normal | 1 Vs All | Species | Overlay") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_F__alternative_overlay_mlPerfAll10k_rep1_Allcancers_BDN_species.jpeg", dpi = "retina",
         width = 10, height = 6, units = "in")
# Save underlying data
mlPerfAll10k_8cancer_Overlay %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Species",datasetName)) %>%
  filter(!grepl("CT",datasetName)) %>%
  droplevels() %>% write.csv("Figures_data/Figure_4/figure_4_F__alternative_overlay_mlPerfAll10k_rep1_Allcancers_BDN_species.jpeg")

#----------------------------------------------------------#
# Plot machine learning performances for all TCGA cancers on raw data split by sequencing center
#----------------------------------------------------------#

source("Supporting_scripts/S05-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Raw <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_by_seq_center_ALL_02Oct21.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_Raw$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_Raw$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_Raw <- mlPerfAll10k_Allcancer_Raw[,!(colnames(mlPerfAll10k_Allcancer_Raw) == "X")]
colnames(mlPerfAll10k_Allcancer_Raw)[1:2] <- c("AUROC","AUPR")

# Rename entries in the "datasetName" column

mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Decontam_HMS"] <- "HMS species decontam (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Decontam_HMS_Cov_Nonzero"] <- "HMS species high coverage (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Decontam_MDA"] <- "MDA species decontam (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Decontam_MDA_Cov_Nonzero"] <- "MDA species high coverage (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Decontam_BCM"] <- "BCM species decontam (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Decontam_BCM_Cov_Nonzero"] <- "BCM species high coverage (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Decontam_WashU"] <- "WashU species decontam (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Decontam_WashU_Cov_Nonzero"] <- "WashU species high coverage (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Decontam_Broad_WGS"] <- "Broad species decontam (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Decontam_Broad_WGS_Cov_Nonzero"] <- "Broad species high coverage (WGS)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Decontam_UNC"] <- "UNC species decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Decontam_UNC_Cov_Nonzero"] <- "UNC species high coverage (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Decontam_CMS"] <- "CMS species decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw$datasetName[mlPerfAll10k_Allcancer_Raw$datasetName == "rep200_HiSeq_Fungi_Decontam_CMS_Cov_Nonzero"] <- "CMS species high coverage (RNA-Seq)"
## NOTE: Broad_RNA only included GBM tumors, so no ML comparisons were made

mlPerfAll10k_Allcancer_Raw$datasetName <- factor(mlPerfAll10k_Allcancer_Raw$datasetName,
                                                    levels = c("HMS species decontam (WGS)",
                                                               "HMS species high coverage (WGS)",
                                                               "MDA species decontam (WGS)",
                                                               "MDA species high coverage (WGS)",
                                                               "BCM species decontam (WGS)",
                                                               "BCM species high coverage (WGS)",
                                                               "WashU species decontam (WGS)",
                                                               "WashU species high coverage (WGS)",
                                                               "Broad species decontam (WGS)",
                                                               "Broad species high coverage (WGS)",
                                                               "UNC species decontam (RNA-Seq)",
                                                               "UNC species high coverage (RNA-Seq)",
                                                               "CMS species decontam (RNA-Seq)",
                                                               "CMS species high coverage (RNA-Seq)"))
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
# HMS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_HMS_PT_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")
# MDA
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("MDA",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_MDA_PT_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")
# BCM
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("BCM",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_BCM_PT_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")
# WashU
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("WashU",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Washington University (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_WashU_PT_species.jpeg", dpi = "retina",
         width = 5, height = 4, units = "in")
# Broad WGS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Broad",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_Broad_WGS_PT_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")
# UNC
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("UNC",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University of North Carolina (RNA-Seq) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_UNC_PT_species.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# CMS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("CMS",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_CMS_PT_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#
# HMS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_HMS_PT_vs_NAT_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")
# MDA does not have sufficient PT and NAT samples to compare
# BCM
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_BCM_PT_vs_NAT_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")
# WashU does not have sufficient PT and NAT samples to compare
# Broad WGS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")
# UNC
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("UNC",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University of North Carolina (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_UNC_PT_vs_NAT_species.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# CMS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("CMS",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_CMS_PT_vs_NAT_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")

#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
# HMS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_HMS_BDN_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")
# MDA
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("MDA",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_MDA_BDN_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")
# BCM
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_BCM_BDN_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")
# WashU
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("WashU",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Washington University (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_WashU_BDN_species.jpeg", dpi = "retina",
         width = 5, height = 4, units = "in")
# Broad WGS
mlPerfAll10k_Allcancer_Raw %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2.5) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_Broad_WGS_BDN_species.jpeg", dpi = "retina",
         width = 6, height = 4, units = "in")
# UNC does not have any BDN samples
# CMS does not have any BDN samples

#----------------------------------------------------------#
# Plot machine learning performances for all TCGA cancers on raw data split by sequencing center
# AND aggregated to each taxa level AND intersected with Weizmann data
#----------------------------------------------------------#

source("Supporting_scripts/S05-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann <- read.csv("Interim_data/rep_perfFungi_10k_rep1_tcga_by_seq_center_taxa_level_and_WIS_intersect_ALL_02Oct21.csv", stringsAsFactors = FALSE)
abbreviationsTCGA_Allcancer <- read.csv("Supporting_data/tcga_abbreviations.csv", stringsAsFactors = FALSE, row.names = 1)
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$abbrev <- abbreviationsTCGA_Allcancer[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$diseaseType,"abbrev"]
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann <- mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann[,!(colnames(mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann) == "X")]
colnames(mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann)[1:2] <- c("AUROC","AUPR")

# Rename entries in the "datasetName" column

# HMS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_HMS_phylum"] <- "HMS phylum decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_HMS_class"] <- "HMS class decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_HMS_order"] <- "HMS order decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_HMS_family"] <- "HMS family decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_HMS_genus"] <- "HMS genus decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_HMS_species"] <- "HMS species decontam (WGS)"
# HMS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_HMS_phylum_Shared_Nonzero"] <- "HMS phylum ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_HMS_class_Shared_Nonzero"] <- "HMS class ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_HMS_order_Shared_Nonzero"] <- "HMS order ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_HMS_family_Shared_Nonzero"] <- "HMS family ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_HMS_genus_Shared_Nonzero"] <- "HMS genus ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_HMS_species_Shared_Nonzero"] <- "HMS species ∩ WIS (WGS)"
# BCM - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_BCM_phylum"] <- "BCM phylum decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_BCM_class"] <- "BCM class decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_BCM_order"] <- "BCM order decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_BCM_family"] <- "BCM family decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_BCM_genus"] <- "BCM genus decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_BCM_species"] <- "BCM species decontam (WGS)"
# BCM - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_BCM_phylum_Shared_Nonzero"] <- "BCM phylum ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_BCM_class_Shared_Nonzero"] <- "BCM class ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_BCM_order_Shared_Nonzero"] <- "BCM order ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_BCM_family_Shared_Nonzero"] <- "BCM family ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_BCM_genus_Shared_Nonzero"] <- "BCM genus ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_BCM_species_Shared_Nonzero"] <- "BCM species ∩ WIS (WGS)"
# MDA - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_MDA_phylum"] <- "MDA phylum decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_MDA_class"] <- "MDA class decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_MDA_order"] <- "MDA order decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_MDA_family"] <- "MDA family decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_MDA_genus"] <- "MDA genus decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_MDA_species"] <- "MDA species decontam (WGS)"
# MDA - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_MDA_phylum_Shared_Nonzero"] <- "MDA phylum ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_MDA_class_Shared_Nonzero"] <- "MDA class ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_MDA_order_Shared_Nonzero"] <- "MDA order ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_MDA_family_Shared_Nonzero"] <- "MDA family ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_MDA_genus_Shared_Nonzero"] <- "MDA genus ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_MDA_species_Shared_Nonzero"] <- "MDA species ∩ WIS (WGS)"
# WashU - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_WashU_phylum"] <- "WashU phylum decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_WashU_class"] <- "WashU class decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_WashU_order"] <- "WashU order decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_WashU_family"] <- "WashU family decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_WashU_genus"] <- "WashU genus decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_WashU_species"] <- "WashU species decontam (WGS)"
# WashU - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_WashU_phylum_Shared_Nonzero"] <- "WashU phylum ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_WashU_class_Shared_Nonzero"] <- "WashU class ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_WashU_order_Shared_Nonzero"] <- "WashU order ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_WashU_family_Shared_Nonzero"] <- "WashU family ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_WashU_genus_Shared_Nonzero"] <- "WashU genus ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_WashU_species_Shared_Nonzero"] <- "WashU species ∩ WIS (WGS)"
# Broad_WGS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_phylum"] <- "Broad phylum decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_class"] <- "Broad class decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_order"] <- "Broad order decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_family"] <- "Broad family decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_genus"] <- "Broad genus decontam (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species"] <- "Broad species decontam (WGS)"
# Broad_WGS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_phylum_Shared_Nonzero"] <- "Broad phylum ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_class_Shared_Nonzero"] <- "Broad class ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_order_Shared_Nonzero"] <- "Broad order ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_family_Shared_Nonzero"] <- "Broad family ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_genus_Shared_Nonzero"] <- "Broad genus ∩ WIS (WGS)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_WGS_species_Shared_Nonzero"] <- "Broad species ∩ WIS (WGS)"
# UNC - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_UNC_phylum"] <- "UNC phylum decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_UNC_class"] <- "UNC class decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_UNC_order"] <- "UNC order decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_UNC_family"] <- "UNC family decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_UNC_genus"] <- "UNC genus decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_UNC_species"] <- "UNC species decontam (RNA-Seq)"
# UNC - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_UNC_phylum_Shared_Nonzero"] <- "UNC phylum ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_UNC_class_Shared_Nonzero"] <- "UNC class ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_UNC_order_Shared_Nonzero"] <- "UNC order ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_UNC_family_Shared_Nonzero"] <- "UNC family ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_UNC_genus_Shared_Nonzero"] <- "UNC genus ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_UNC_species_Shared_Nonzero"] <- "UNC species ∩ WIS (RNA-Seq)"
# CMS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_CMS_phylum"] <- "CMS phylum decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_CMS_class"] <- "CMS class decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_CMS_order"] <- "CMS order decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_CMS_family"] <- "CMS family decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_CMS_genus"] <- "CMS genus decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_CMS_species"] <- "CMS species decontam (RNA-Seq)"
# CMS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_CMS_phylum_Shared_Nonzero"] <- "CMS phylum ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_CMS_class_Shared_Nonzero"] <- "CMS class ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_CMS_order_Shared_Nonzero"] <- "CMS order ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_CMS_family_Shared_Nonzero"] <- "CMS family ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_CMS_genus_Shared_Nonzero"] <- "CMS genus ∩ WIS (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_CMS_species_Shared_Nonzero"] <- "CMS species ∩ WIS (RNA-Seq)"
# Broad_RNA - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_phylum"] <- "Broad phylum decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_class"] <- "Broad class decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_order"] <- "Broad order decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_family"] <- "Broad family decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_genus"] <- "Broad genus decontam (RNA-Seq)"
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName[mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName == "df_psRep200_HiSeq_Fungi_Decontam_Broad_RNA_species"] <- "Broad species decontam (RNA-Seq)"

## NOTE: Broad_RNA only included GBM tumors, so no ML comparisons were made

mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName <- factor(mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann$datasetName,
                                                 levels = c("HMS phylum decontam (WGS)",
                                                            "HMS class decontam (WGS)",
                                                            "HMS order decontam (WGS)",
                                                            "HMS family decontam (WGS)",
                                                            "HMS genus decontam (WGS)",
                                                            "HMS species decontam (WGS)",
                                                            "HMS phylum ∩ WIS (WGS)",
                                                            "HMS class ∩ WIS (WGS)",
                                                            "HMS order ∩ WIS (WGS)",
                                                            "HMS family ∩ WIS (WGS)",
                                                            "HMS genus ∩ WIS (WGS)",
                                                            "HMS species ∩ WIS (WGS)",
                                                            # BCM
                                                            "BCM phylum decontam (WGS)",
                                                            "BCM class decontam (WGS)",
                                                            "BCM order decontam (WGS)",
                                                            "BCM family decontam (WGS)",
                                                            "BCM genus decontam (WGS)",
                                                            "BCM species decontam (WGS)",
                                                            "BCM phylum ∩ WIS (WGS)",
                                                            "BCM class ∩ WIS (WGS)",
                                                            "BCM order ∩ WIS (WGS)",
                                                            "BCM family ∩ WIS (WGS)",
                                                            "BCM genus ∩ WIS (WGS)",
                                                            "BCM species ∩ WIS (WGS)",
                                                            # MDA
                                                            "MDA phylum decontam (WGS)",
                                                            "MDA class decontam (WGS)",
                                                            "MDA order decontam (WGS)",
                                                            "MDA family decontam (WGS)",
                                                            "MDA genus decontam (WGS)",
                                                            "MDA species decontam (WGS)",
                                                            "MDA phylum ∩ WIS (WGS)",
                                                            "MDA class ∩ WIS (WGS)",
                                                            "MDA order ∩ WIS (WGS)",
                                                            "MDA family ∩ WIS (WGS)",
                                                            "MDA genus ∩ WIS (WGS)",
                                                            "MDA species ∩ WIS (WGS)",
                                                            # WashU
                                                            "WashU phylum decontam (WGS)",
                                                            "WashU class decontam (WGS)",
                                                            "WashU order decontam (WGS)",
                                                            "WashU family decontam (WGS)",
                                                            "WashU genus decontam (WGS)",
                                                            "WashU species decontam (WGS)",
                                                            "WashU phylum ∩ WIS (WGS)",
                                                            "WashU class ∩ WIS (WGS)",
                                                            "WashU order ∩ WIS (WGS)",
                                                            "WashU family ∩ WIS (WGS)",
                                                            "WashU genus ∩ WIS (WGS)",
                                                            "WashU species ∩ WIS (WGS)",
                                                            # Broad
                                                            "Broad phylum decontam (WGS)",
                                                            "Broad class decontam (WGS)",
                                                            "Broad order decontam (WGS)",
                                                            "Broad family decontam (WGS)",
                                                            "Broad genus decontam (WGS)",
                                                            "Broad species decontam (WGS)",
                                                            "Broad phylum ∩ WIS (WGS)",
                                                            "Broad class ∩ WIS (WGS)",
                                                            "Broad order ∩ WIS (WGS)",
                                                            "Broad family ∩ WIS (WGS)",
                                                            "Broad genus ∩ WIS (WGS)",
                                                            "Broad species ∩ WIS (WGS)",
                                                            # UNC
                                                            "UNC phylum decontam (RNA-Seq)",
                                                            "UNC class decontam (RNA-Seq)",
                                                            "UNC order decontam (RNA-Seq)",
                                                            "UNC family decontam (RNA-Seq)",
                                                            "UNC genus decontam (RNA-Seq)",
                                                            "UNC species decontam (RNA-Seq)",
                                                            "UNC phylum ∩ WIS (RNA-Seq)",
                                                            "UNC class ∩ WIS (RNA-Seq)",
                                                            "UNC order ∩ WIS (RNA-Seq)",
                                                            "UNC family ∩ WIS (RNA-Seq)",
                                                            "UNC genus ∩ WIS (RNA-Seq)",
                                                            "UNC species ∩ WIS (RNA-Seq)",
                                                            # CMS
                                                            "CMS phylum decontam (RNA-Seq)",
                                                            "CMS class decontam (RNA-Seq)",
                                                            "CMS order decontam (RNA-Seq)",
                                                            "CMS family decontam (RNA-Seq)",
                                                            "CMS genus decontam (RNA-Seq)",
                                                            "CMS species decontam (RNA-Seq)",
                                                            "CMS phylum ∩ WIS (RNA-Seq)",
                                                            "CMS class ∩ WIS (RNA-Seq)",
                                                            "CMS order ∩ WIS (RNA-Seq)",
                                                            "CMS family ∩ WIS (RNA-Seq)",
                                                            "CMS genus ∩ WIS (RNA-Seq)",
                                                            "CMS species ∩ WIS (RNA-Seq)",
                                                            # Broad
                                                            "Broad phylum decontam (RNA-Seq)",
                                                            "Broad class decontam (RNA-Seq)",
                                                            "Broad order decontam (RNA-Seq)",
                                                            "Broad family decontam (RNA-Seq)",
                                                            "Broad genus decontam (RNA-Seq)",
                                                            "Broad species decontam (RNA-Seq)",
                                                            "Broad phylum ∩ WIS (RNA-Seq)",
                                                            "Broad class ∩ WIS (RNA-Seq)",
                                                            "Broad order ∩ WIS (RNA-Seq)",
                                                            "Broad family ∩ WIS (RNA-Seq)",
                                                            "Broad genus ∩ WIS (RNA-Seq)",
                                                            "Broad species ∩ WIS (RNA-Seq)"))
#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
# HMS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_HMS_PT_taxa_level_decontam.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# HMS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data Intersected with WIS Features  | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_HMS_PT_taxa_level_shared.jpeg", dpi = "retina",
       width = 8, height = 4, units = "in")

# BCM - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("BCM",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_BCM_PT_taxa_level_decontam.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# BCM - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("BCM",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data Intersected with WIS Features  | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_BCM_PT_taxa_level_shared.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")

# MDA - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("MDA",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_MDA_PT_taxa_level_decontam.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# MDA - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("MDA",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data Intersected with WIS Features  | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_MDA_PT_taxa_level_shared.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")

# WashU - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("WashU",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_WashU_PT_taxa_level_decontam.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# WashU - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("WashU",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data Intersected with WIS Features  | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_WashU_PT_taxa_level_shared.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")

# Broad_WGS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Broad",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_Broad_WGS_PT_taxa_level_decontam.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# Broad_WGS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Broad",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data Intersected with WIS Features  | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_Broad_WGS_PT_taxa_level_shared.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")

# UNC - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("UNC",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_UNC_PT_taxa_level_decontam.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# UNC - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("UNC",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data Intersected with WIS Features  | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_UNC_PT_taxa_level_shared.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")

# CMS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("CMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_CMS_PT_taxa_level_decontam.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# CMS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("CMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data Intersected with WIS Features  | Primary Tumor | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_CMS_PT_taxa_level_shared.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#
# HMS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_HMS_PT_vs_NAT_taxa_level_decontam.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# HMS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data Intersected with WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_HMS_PT_vs_NAT_taxa_level_shared.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")

# BCM - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_BCM_PT_vs_NAT_taxa_level_decontam.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# BCM - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data Intersected with WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_BCM_PT_vs_NAT_taxa_level_shared.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")

## NOTE: Neither MDA nor WashU has enough tumor vs. normal samples to plot

# Broad_WGS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT_taxa_level_decontam.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# Broad_WGS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data Intersected with WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT_taxa_level_shared.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")

# UNC - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("UNC",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_UNC_PT_vs_NAT_taxa_level_decontam.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# UNC - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("UNC",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data Intersected with WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_UNC_PT_vs_NAT_taxa_level_shared.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")

# CMS - decontam
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("CMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(!grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_CMS_PT_vs_NAT_taxa_level_decontam.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# CMS - shared
mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("CMS",datasetName)) %>%
  # filter(grepl("species",datasetName)) %>%
  filter(grepl("∩",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data Intersected with WIS Features  | Primary Tumor vs Solid Tissue Normal | 1 Vs All\n(Phylum, Class, Order, Family, Genus, Species)") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_CMS_PT_vs_NAT_taxa_level_shared.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")

#----------------------------------------------------------#
# Overlay plots machine learning performances for all TCGA cancers
# using RAW data per seq center:
# - Decontam species data
# - High coverage species
# - Intersected species with WIS
#----------------------------------------------------------#

source("Supporting_scripts/S05-SummarySE.R") # Contains a function that calculates std error and 95% confidence intervals

# NOTE: There are too many cancer types to include more than two errorbars per cancer type,
# so the WGS vs. RNA comparison is not included in this overlay. If you would like to include
# it, you can remove the commented line below.
mlPerfAll10k_Allcancer_Raw_Overlay <- rbind(mlPerfAll10k_Allcancer_Raw_Taxa_Levels_and_Weizmann,
                                            mlPerfAll10k_Allcancer_Raw)
mlPerfAll10k_Allcancer_Raw_Overlay_Filt <- mlPerfAll10k_Allcancer_Raw_Overlay %>%
  filter(grepl("species",datasetName)) %>% droplevels()
mlPerfAll10k_Allcancer_Raw_Overlay_Filt$datasetName <- factor(mlPerfAll10k_Allcancer_Raw_Overlay_Filt$datasetName,
                                                              levels = rev(levels(mlPerfAll10k_Allcancer_Raw_Overlay_Filt$datasetName)))

#-------------------------Plot primary tumor 1 vs. all others performance-------------------------#
# HMS - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("HMS",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_HMS_PT_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# BCM - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("BCM",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_BCM_PT_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# MDA - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("MDA",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_MDA_PT_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# WashU - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("WashU",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_WashU_PT_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# Broad_WGS - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("Broad",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_Broad_WGS_PT_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# UNC - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("UNC",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_UNC_PT_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# CMS - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor") %>%
  filter(grepl("CMS",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data | Primary Tumor | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_CMS_PT_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")

#-------------------------Plot primary tumor vs. NAT performance-------------------------#
# HMS - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_HMS_PT_vs_NAT_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# BCM - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_BCM_PT_vs_NAT_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
## NOTE: Neither MDA nor WashU had sufficient samples to plot primary tumor vs. NAT performance

# Broad_WGS - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_Broad_WGS_PT_vs_NAT_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# UNC - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("UNC",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("University North Carolina (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_UNC_PT_vs_NAT_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# CMS - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Primary Tumor vs Solid Tissue Normal") %>%
  filter(grepl("CMS",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Canada's Michael Smith Genome Sciences Centre (RNA-Seq) | Raw Data | Primary Tumor vs Solid Tissue Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_CMS_PT_vs_NAT_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")

#-------------------------Plot blood derived normal 1 vs. all others performance-------------------------#
# HMS - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("HMS",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Harvard Medical School (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_HMS_BDN_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# BCM - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("BCM",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Baylor College of Medicine (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_BCM_BDN_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# MDA - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("MDA",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("MD Anderson (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_MDA_BDN_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# WashU - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("WashU",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("WashU (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_WashU_BDN_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")
# Broad_WGS - overlay
mlPerfAll10k_Allcancer_Raw_Overlay_Filt %>%
  filter(sampleType == "Blood Derived Normal") %>%
  filter(grepl("Broad",datasetName)) %>%
  droplevels() %>%
  reshape2::melt(id.vars = c("rep","abbrev","diseaseType","sampleType","datasetName","metadataName","minorityClassSize","majorityClassSize","minorityClassName","majorityClassName")) %>%
  summarySE(measurevar = "value", groupvars = c("datasetName","metadataName","variable","abbrev")) %>%
  ggplot(aes(reorder(abbrev, value, FUN=median),value, color=datasetName)) +
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=0, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size=2) + xlab("Cancer type") + ylab("Area Under Curve") + theme_pubr() +
  facet_wrap(~variable) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  ggtitle("Broad Institute (WGS) | Raw Data | Blood Derived Normal | 1 Vs All | Species") + theme(plot.title = element_text(hjust = 0.5)) +
  rotate_x_text(90) + scale_color_nejm(name = "Features") + geom_hline(yintercept = 1, linetype="dashed") + 
  ggsave("Figures/Figure_4/figure_4_XX_mlPerfAll10k_rep1_Broad_WGS_BDN_species_overlay.jpeg", dpi = "retina",
         width = 8, height = 4, units = "in")

