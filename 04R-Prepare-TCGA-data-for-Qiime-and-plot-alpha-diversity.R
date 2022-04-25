#-----------------------------------------------------------------------------
# 04-Prepare-TCGA-data-for-Qiime-and-plot-alpha-diversity.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Prepare TCGA data for Qiime analyses to calculate alpha and beta diversities
# - Plot Qiime alpha diversity results (beta diversities to be shown using Emperor)
# - Prepare TCGA data for MMvec analyses to examine co-occurring fungi and bacteria
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
require(EnvStats) # to show sample sizes
require(biomformat)
require(Rhdf5lib)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Import data
#----------------------------------------------------------#

load("Interim_data/snmDataFungi_DecontamV2_25Mar22.RData", verbose = TRUE) # To load the metadata and raw data objects
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
load("Interim_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_25Mar22.RData", verbose = TRUE)
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
# e.g. summary(rowSums(rep200_HiSeq_Fungi_Broad_WGS_Cov)==0)
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
     file = "Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy_with_coverage_filter_decontamV2_29Mar22.RData")

# Scripts: S12B

# load("Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy_16Sep21.RData")

#----------------------------------------------------------#
# Prepare TCGA data for Qiime
# Formerly without decontamination but now with (Mar 29 2022)
#----------------------------------------------------------#
# Subset to primary tumor data
metaQiitaCombined_Nonzero_HMS_PT <- metaQiitaCombined_Nonzero_DecontamV2_HMS %>%
  filter(sample_type == "Primary Tumor") %>%
  rename("submitted_name" = "sample_name") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>% droplevels()
metaQiitaCombined_Nonzero_BCM_PT <- metaQiitaCombined_Nonzero_DecontamV2_BCM %>%
  filter(sample_type == "Primary Tumor") %>%
  rename("submitted_name" = "sample_name") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>% droplevels()
metaQiitaCombined_Nonzero_MDA_PT <- metaQiitaCombined_Nonzero_DecontamV2_MDA %>%
  filter(sample_type == "Primary Tumor") %>%
  rename("submitted_name" = "sample_name") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>% droplevels()
metaQiitaCombined_Nonzero_WashU_PT <- metaQiitaCombined_Nonzero_DecontamV2_WashU %>%
  filter(sample_type == "Primary Tumor") %>%
  rename("submitted_name" = "sample_name") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>% droplevels()
metaQiitaCombined_Nonzero_Broad_WGS_PT <- metaQiitaCombined_Nonzero_DecontamV2_Broad_WGS %>%
  filter(sample_type == "Primary Tumor") %>%
  rename("submitted_name" = "sample_name") %>%
  select(investigation, disease_type, experimental_strategy,data_submitting_center_label) %>%
  rownames_to_column("sampleid") %>% droplevels()

# Write metadata to text files
write.table(metaQiitaCombined_Nonzero_HMS_PT, file = "Qiime_data_and_scripts/Qiime_input_data/metaQiitaCombined_Nonzero_HMS_PT.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
write.table(metaQiitaCombined_Nonzero_BCM_PT, file = "Qiime_data_and_scripts/Qiime_input_data/metaQiitaCombined_Nonzero_BCM_PT.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
write.table(metaQiitaCombined_Nonzero_MDA_PT, file = "Qiime_data_and_scripts/Qiime_input_data/metaQiitaCombined_Nonzero_MDA_PT.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
write.table(metaQiitaCombined_Nonzero_WashU_PT, file = "Qiime_data_and_scripts/Qiime_input_data/metaQiitaCombined_Nonzero_WashU_PT.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)
write.table(metaQiitaCombined_Nonzero_Broad_WGS_PT, file = "Qiime_data_and_scripts/Qiime_input_data/metaQiitaCombined_Nonzero_Broad_WGS_PT.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)

## Subset count data
rep200_HiSeq_Fungi_Decontam_HMS_PT <- rep200_HiSeq_Fungi_DecontamV2_HMS[metaQiitaCombined_Nonzero_HMS_PT$sampleid,]
rep200_HiSeq_Fungi_Decontam_BCM_PT <- rep200_HiSeq_Fungi_DecontamV2_BCM[metaQiitaCombined_Nonzero_BCM_PT$sampleid,]
rep200_HiSeq_Fungi_Decontam_MDA_PT <- rep200_HiSeq_Fungi_DecontamV2_MDA[metaQiitaCombined_Nonzero_MDA_PT$sampleid,]
rep200_HiSeq_Fungi_Decontam_WashU_PT <- rep200_HiSeq_Fungi_DecontamV2_WashU[metaQiitaCombined_Nonzero_WashU_PT$sampleid,]
rep200_HiSeq_Fungi_Decontam_Broad_WGS_PT <- rep200_HiSeq_Fungi_DecontamV2_Broad_WGS[metaQiitaCombined_Nonzero_Broad_WGS_PT$sampleid,]

## Save count data as biom tables
rep200_HiSeq_Fungi_Decontam_HMS_PT_BIOM <- make_biom(t(rep200_HiSeq_Fungi_Decontam_HMS_PT))
write_biom(rep200_HiSeq_Fungi_Decontam_HMS_PT_BIOM, biom_file = "Qiime_data_and_scripts/Qiime_input_data/rep200_HiSeq_Fungi_Decontam_HMS_PT.biom")

rep200_HiSeq_Fungi_Decontam_BCM_PT_BIOM <- make_biom(t(rep200_HiSeq_Fungi_Decontam_BCM_PT))
write_biom(rep200_HiSeq_Fungi_Decontam_BCM_PT_BIOM, biom_file = "Qiime_data_and_scripts/Qiime_input_data/rep200_HiSeq_Fungi_Decontam_BCM_PT.biom")

rep200_HiSeq_Fungi_Decontam_MDA_PT_BIOM <- make_biom(t(rep200_HiSeq_Fungi_Decontam_MDA_PT))
write_biom(rep200_HiSeq_Fungi_Decontam_MDA_PT_BIOM, biom_file = "Qiime_data_and_scripts/Qiime_input_data/rep200_HiSeq_Fungi_Decontam_MDA_PT.biom")

rep200_HiSeq_Fungi_Decontam_WashU_PT_BIOM <- make_biom(t(rep200_HiSeq_Fungi_Decontam_WashU_PT))
write_biom(rep200_HiSeq_Fungi_Decontam_WashU_PT_BIOM, biom_file = "Qiime_data_and_scripts/Qiime_input_data/rep200_HiSeq_Fungi_Decontam_WashU_PT.biom")

rep200_HiSeq_Fungi_Decontam_Broad_WGS_PT_BIOM <- make_biom(t(rep200_HiSeq_Fungi_Decontam_Broad_WGS_PT))
write_biom(rep200_HiSeq_Fungi_Decontam_Broad_WGS_PT_BIOM, biom_file = "Qiime_data_and_scripts/Qiime_input_data/rep200_HiSeq_Fungi_Decontam_Broad_WGS_PT.biom")

#----------------------------------------------------------#
# Use Qiime to derive alpha and beta diversity results
#----------------------------------------------------------#

# Scripts for doing these steps are outlined in the 
# "Qiime_data_and_scripts" folder and some of their outputs
# will be used to create plots below for alpha diversity.

#----------------------------------------------------------#
# Import alpha diversity results from Qiime
#----------------------------------------------------------#
# NOTE: The loaded csv files were acquired by uploading the corresponding .qzv file
# to https://view.qiime2.org/ and then clicking "Download raw data as TSV" followed by
# (i) deleting the 2nd row (contains unnecessary column type info) and (ii) saving as a CSV file.

# HMS
alphaDivOtus_HMS <- read.csv("Qiime_data_and_scripts_resubmission_version/core_metrics_hms_5k/alpha_div_observed_features_hms_5k.csv", stringsAsFactors = FALSE)
alphaDivShannon_HMS <- read.csv("Qiime_data_and_scripts_resubmission_version/core_metrics_hms_5k/alpha_div_shannon_hms_5k.csv", stringsAsFactors = FALSE)
alphaDivOtus_HMS$investigation <- gsub("^TCGA\\-","",alphaDivOtus_HMS$investigation)
alphaDivShannon_HMS$investigation <- gsub("^TCGA\\-","",alphaDivShannon_HMS$investigation)
# BCM
alphaDivOtus_BCM <- read.csv("Qiime_data_and_scripts_resubmission_version/core_metrics_bcm_5k/alpha_div_observed_features_bcm_5k.csv", stringsAsFactors = FALSE)
alphaDivShannon_BCM <- read.csv("Qiime_data_and_scripts_resubmission_version/core_metrics_bcm_5k/alpha_div_shannon_bcm_5k.csv", stringsAsFactors = FALSE)
alphaDivOtus_BCM$investigation <- gsub("^TCGA\\-","",alphaDivOtus_BCM$investigation)
alphaDivShannon_BCM$investigation <- gsub("^TCGA\\-","",alphaDivShannon_BCM$investigation)
# MDA
alphaDivOtus_MDA <- read.csv("Qiime_data_and_scripts_resubmission_version/core_metrics_mda_5k/alpha_div_observed_features_mda_5k.csv", stringsAsFactors = FALSE)
alphaDivShannon_MDA <- read.csv("Qiime_data_and_scripts_resubmission_version/core_metrics_mda_5k/alpha_div_shannon_mda_5k.csv", stringsAsFactors = FALSE)
alphaDivOtus_MDA$investigation <- gsub("^TCGA\\-","",alphaDivOtus_MDA$investigation)
alphaDivShannon_MDA$investigation <- gsub("^TCGA\\-","",alphaDivShannon_MDA$investigation)
# WashU
alphaDivOtus_WashU <- read.csv("Qiime_data_and_scripts_resubmission_version/core_metrics_washU_5k/alpha_div_observed_features_washU_5k.csv", stringsAsFactors = FALSE)
alphaDivShannon_WashU <- read.csv("Qiime_data_and_scripts_resubmission_version/core_metrics_washU_5k/alpha_div_shannon_washU_5k.csv", stringsAsFactors = FALSE)
alphaDivOtus_WashU$investigation <- gsub("^TCGA\\-","",alphaDivOtus_WashU$investigation)
alphaDivShannon_WashU$investigation <- gsub("^TCGA\\-","",alphaDivShannon_WashU$investigation)
# Broad WGS
alphaDivOtus_Broad_WGS <- read.csv("Qiime_data_and_scripts_resubmission_version/core_metrics_broad_WGS_2k/alpha_div_observed_features_broad_WGS_2k.csv", stringsAsFactors = FALSE)
alphaDivShannon_Broad_WGS <- read.csv("Qiime_data_and_scripts_resubmission_version/core_metrics_broad_WGS_2k/alpha_div_shannon_broad_WGS_2k.csv", stringsAsFactors = FALSE)
alphaDivOtus_Broad_WGS$investigation <- gsub("^TCGA\\-","",alphaDivOtus_Broad_WGS$investigation)
alphaDivShannon_Broad_WGS$investigation <- gsub("^TCGA\\-","",alphaDivShannon_Broad_WGS$investigation)

#----------------------------------------------------------#
# Plot Qiime alpha diversity data for HMS and MDA with matching cancer type colors
#----------------------------------------------------------#
colors <- pal_d3("category20")(15)
table(alphaDivOtus_HMS$investigation)
# HMS observed features
alphaDivOtus_HMS %>%
  mutate(palColors = case_when(
      investigation == "BLCA" ~ colors[1],
      investigation == "BRCA" ~ colors[2],
      investigation == "CESC" ~ colors[3],
      investigation == "COAD" ~ colors[4],
      investigation == "HNSC" ~ colors[5],
      investigation == "LGG" ~ colors[6],
      investigation == "LUAD" ~ colors[7],
      investigation == "PRAD" ~ colors[8],
      investigation == "READ" ~ colors[9],
      investigation == "SKCM" ~ colors[10],
      investigation == "STAD" ~ colors[11],
      investigation == "THCA" ~ colors[12],
      investigation == "UCEC" ~ colors[13]
    )) -> alphaDivOtus_HMS_colored
# HMS shannon
alphaDivShannon_HMS %>%
  mutate(palColors = case_when(
    investigation == "BLCA" ~ colors[1],
    investigation == "BRCA" ~ colors[2],
    investigation == "CESC" ~ colors[3],
    investigation == "COAD" ~ colors[4],
    investigation == "HNSC" ~ colors[5],
    investigation == "LGG" ~ colors[6],
    investigation == "LUAD" ~ colors[7],
    investigation == "PRAD" ~ colors[8],
    investigation == "READ" ~ colors[9],
    investigation == "SKCM" ~ colors[10],
    investigation == "STAD" ~ colors[11],
    investigation == "THCA" ~ colors[12],
    investigation == "UCEC" ~ colors[13]
  )) -> alphaDivShannon_HMS_colored

table(alphaDivOtus_MDA$investigation)
# MDA observed features
alphaDivOtus_MDA %>%
  mutate(
    palColors = case_when(
      investigation == "BLCA" ~ colors[1],
      investigation == "BRCA" ~ colors[2],
      investigation == "CESC" ~ colors[3],
      # investigation == "COAD" ~ colors[4],
      investigation == "HNSC" ~ colors[5],
      # investigation == "LGG" ~ colors[6],
      # investigation == "LUAD" ~ colors[7],
      # investigation == "PRAD" ~ colors[8],
      # investigation == "READ" ~ colors[9],
      # investigation == "SKCM" ~ colors[10],
      investigation == "STAD" ~ colors[11],
      investigation == "THCA" ~ colors[12],
      # investigation == "UCEC" ~ colors[13]
      investigation == "ESCA" ~ colors[14],
      investigation == "UVM" ~ colors[15]
    )) -> alphaDivOtus_MDA_colored
# MDA shannon
alphaDivShannon_MDA %>%
  mutate(
    palColors = case_when(
      investigation == "BLCA" ~ colors[1],
      investigation == "BRCA" ~ colors[2],
      investigation == "CESC" ~ colors[3],
      # investigation == "COAD" ~ colors[4],
      investigation == "HNSC" ~ colors[5],
      # investigation == "LGG" ~ colors[6],
      # investigation == "LUAD" ~ colors[7],
      # investigation == "PRAD" ~ colors[8],
      # investigation == "READ" ~ colors[9],
      # investigation == "SKCM" ~ colors[10],
      investigation == "STAD" ~ colors[11],
      investigation == "THCA" ~ colors[12],
      # investigation == "UCEC" ~ colors[13]
      investigation == "ESCA" ~ colors[14],
      investigation == "UVM" ~ colors[15]
    )) -> alphaDivShannon_MDA_colored

# HMS observed features
alphaDivOtus_HMS_colored %>%
  ggplot(aes(reorder(investigation, observed_otus, FUN=median), observed_otus, fill=palColors)) +
  geom_boxplot() + xlab("TCGA cancer type (Harvard Med School)") + ylab("Observed features") + theme_pubr(legend = "none") +
  rotate_x_text(30) +
  scale_fill_manual(values = levels(factor(alphaDivOtus_HMS_colored$palColors))) +
  stat_compare_means(label.x.npc = 0.1, label.y.npc = 0.5) + 
  stat_n_text()
ggsave(filename = "Figures/Supplementary_Figures/alphaDivOtusPlot_HMS.pdf",
         dpi = "retina", units = "in", height = 4, width = 7)
alphaDivOtus_HMS_colored %>% write.csv("Figures_data/Supplementary_Figures/alphaDivOtusPlot_HMS.csv")

# HMS shannon
alphaDivShannon_HMS_colored %>%
  ggplot(aes(reorder(investigation, shannon, FUN=median), shannon, fill=palColors)) +
  geom_boxplot() + xlab("TCGA cancer type (Harvard Med School)") + ylab("Shannon diversity") + theme_pubr(legend = "none") +
  rotate_x_text(30) +
  scale_fill_manual(values = levels(factor(alphaDivShannon_HMS_colored$palColors))) +
  stat_compare_means(label.x.npc = 0.1, label.y.npc = 0.5) + 
  stat_n_text()
ggsave(filename = "Figures/Supplementary_Figures/alphaDivShannonPlot_HMS.pdf",
         dpi = "retina", units = "in", height = 4, width = 7)
alphaDivShannon_HMS_colored %>% write.csv("Figures_data/Supplementary_Figures/alphaDivOtusPlot_HMS.csv")

# MDA observed features
alphaDivOtus_MDA_colored %>%
  ggplot(aes(reorder(investigation, observed_otus, FUN=median), observed_otus, fill=palColors)) +
  geom_boxplot() + xlab("TCGA cancer type (MD Anderson)") + ylab("Observed features") + theme_pubr(legend = "none") +
  rotate_x_text(30) +
  scale_fill_manual(values = levels(factor(alphaDivOtus_MDA_colored$palColors))) +
  stat_compare_means(label.x.npc = 0.1, label.y.npc = 0.5) + 
  stat_n_text()
ggsave(filename = "Figures/Supplementary_Figures/alphaDivOtusPlot_MDA.pdf",
         dpi = "retina", units = "in", height = 4, width = 7)
alphaDivOtus_MDA_colored %>% write.csv("Figures_data/Supplementary_Figures/alphaDivOtusPlot_MDA.csv")

# MDA shannon
alphaDivShannon_MDA_colored %>%
  ggplot(aes(reorder(investigation, shannon, FUN=median), shannon, fill=palColors)) +
  geom_boxplot() + xlab("TCGA cancer type (MD Anderson)") + ylab("Shannon diversity") + theme_pubr(legend = "none") +
  rotate_x_text(30) +
  scale_fill_manual(values = levels(factor(alphaDivShannon_MDA_colored$palColors))) +
  stat_compare_means(label.x.npc = 0.1, label.y.npc = 0.5) + 
  stat_n_text()
ggsave(filename = "Figures/Supplementary_Figures/alphaDivShannonPlot_MDA.pdf",
         dpi = "retina", units = "in", height = 4, width = 7)
alphaDivShannon_MDA_colored %>% write.csv("Figures_data/Supplementary_Figures/alphaDivShannonPlot_MDA.csv")

# BCM observed features
alphaDivOtus_BCM %>%
  ggplot(aes(reorder(investigation, observed_otus, FUN=median), observed_otus, fill=investigation)) +
  geom_boxplot() + xlab("TCGA cancer type (BCM)") + ylab("Observed features") + theme_pubr(legend = "none") +
  rotate_x_text(30) +
  scale_fill_igv() +
  # scale_fill_manual(values = levels(factor(alphaDivOtus_BCM$palColors))) +
  stat_compare_means(label.x.npc = 0.1, label.y.npc = 0.5) + 
  stat_n_text()
ggsave(filename = "Figures/Supplementary_Figures/alphaDivOtusPlot_BCM.pdf",
         dpi = "retina", units = "in", height = 4, width = 7)
alphaDivOtus_BCM %>% write.csv("Figures_data/Supplementary_Figures/alphaDivOtusPlot_BCM.csv")

# BCM shannon
alphaDivShannon_BCM %>%
  ggplot(aes(reorder(investigation, shannon, FUN=median), shannon, fill=investigation)) +
  geom_boxplot() + xlab("TCGA cancer type (BCM)") + ylab("Shannon diversity") + theme_pubr(legend = "none") +
  rotate_x_text(30) +
  scale_fill_igv() +
  # scale_fill_manual(values = levels(factor(alphaDivShannon_BCM_colored$palColors))) +
  stat_compare_means(label.x.npc = 0.1, label.y.npc = 0.5) + 
  stat_n_text()
ggsave(filename = "Figures/Supplementary_Figures/alphaDivShannonPlot_BCM.pdf",
         dpi = "retina", units = "in", height = 4, width = 7)
alphaDivShannon_BCM %>% write.csv("Figures_data/Supplementary_Figures/alphaDivShannonPlot_BCM.csv")

# WashU observed features
alphaDivOtus_WashU %>%
  ggplot(aes(reorder(investigation, observed_otus, FUN=median), observed_otus, fill=investigation)) +
  geom_boxplot() + xlab("TCGA cancer type (WashU)") + ylab("Observed features") + theme_pubr(legend = "none") +
  rotate_x_text(30) +
  scale_fill_igv() +
  # scale_fill_manual(values = levels(factor(alphaDivOtus_WashU$palColors))) +
  stat_compare_means(label.x.npc = 0.1, label.y.npc = 0.5) + 
  stat_n_text()
ggsave(filename = "Figures/Supplementary_Figures/alphaDivOtusPlot_WashU.pdf",
         dpi = "retina", units = "in", height = 4, width = 7)
alphaDivOtus_WashU %>% write.csv("Figures_data/Supplementary_Figures/alphaDivOtusPlot_WashU.csv")

# WashU shannon
alphaDivShannon_WashU %>%
  ggplot(aes(reorder(investigation, shannon, FUN=median), shannon, fill=investigation)) +
  geom_boxplot() + xlab("TCGA cancer type (WashU)") + ylab("Shannon diversity") + theme_pubr(legend = "none") +
  rotate_x_text(30) +
  scale_fill_igv() +
  # scale_fill_manual(values = levels(factor(alphaDivShannon_WashU_colored$palColors))) +
  stat_compare_means(label.x.npc = 0.1, label.y.npc = 0.5) + 
  stat_n_text()
ggsave(filename = "Figures/Supplementary_Figures/alphaDivShannonPlot_WashU.pdf",
         dpi = "retina", units = "in", height = 4, width = 7)
alphaDivShannon_WashU %>% write.csv("Figures_data/Supplementary_Figures/alphaDivShannonPlot_WashU.csv")

# Broad WGS observed features
alphaDivOtus_Broad_WGS %>%
  ggplot(aes(reorder(investigation, observed_otus, FUN=median), observed_otus, fill=investigation)) +
  geom_boxplot() + xlab("TCGA cancer type (Broad WGS)") + ylab("Observed features") + theme_pubr(legend = "none") +
  rotate_x_text(30) +
  scale_fill_igv() +
  # scale_fill_manual(values = levels(factor(alphaDivOtus_Broad_WGS$palColors))) +
  stat_compare_means(label.x.npc = 0.1, label.y.npc = 0.8) + 
  stat_n_text()
ggsave(filename = "Figures/Supplementary_Figures/alphaDivOtusPlot_Broad_WGS.pdf",
         dpi = "retina", units = "in", height = 4, width = 7)
alphaDivOtus_Broad_WGS %>% write.csv("Figures_data/Supplementary_Figures/alphaDivOtusPlot_Broad_WGS.csv")

# Broad WGS shannon
alphaDivShannon_Broad_WGS %>%
  ggplot(aes(reorder(investigation, shannon, FUN=median), shannon, fill=investigation)) +
  geom_boxplot() + xlab("TCGA cancer type (Broad WGS)") + ylab("Shannon diversity") + theme_pubr(legend = "none") +
  rotate_x_text(30) +
  scale_fill_igv() +
  # scale_fill_manual(values = levels(factor(alphaDivShannon_Broad_WGS_colored$palColors))) +
  stat_compare_means(label.x.npc = 0.1, label.y.npc = 0.8) + 
  stat_n_text()
ggsave(filename = "Figures/Supplementary_Figures/alphaDivShannonPlot_Broad_WGS.pdf",
         dpi = "retina", units = "in", height = 4, width = 7)
alphaDivShannon_Broad_WGS %>% write.csv("Figures_data/Supplementary_Figures/alphaDivShannonPlot_Broad_WGS.csv")

