#-----------------------------------------------------------------------------
# 09-Prepare-TCGA-data-for-MMvec-and-ANCOM-BC.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Load data with Weizmann overlapping features
# - Format data for MMvec to examine co-occurences of (bacteria and fungi) and (immune cells and fungi)
# - Perform ANCOM-BC on subsetted data (without batch correction) and on batch corrected data to guide MMvec analyses
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
require(dplyr)
require(tidyr)
require(tibble)
require(reshape2)
require(ggpubr)
require(ggsci)
require(ANCOMBC)
require(biomformat)
require(Rhdf5lib)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Import data
#----------------------------------------------------------#

# snmDataOGUFungi,
# vdge_dataE,
# rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero,
# metaQiitaCombined_Nonzero_SpeciesShared,
load("Interim_data/snmDataFungi_13Sep21.RData")

# save(rep200FungiPhylumShared_Nonzero,
#      metaQiitaCombined_Nonzero_SpeciesShared_PhylumShared,
#      rep200FungiClassShared_Nonzero,
#      metaQiitaCombined_Nonzero_SpeciesShared_ClassShared,
#      rep200FungiOrderShared_Nonzero,
#      metaQiitaCombined_Nonzero_SpeciesShared_OrderShared,
#      rep200FungiFamilyShared_Nonzero,
#      metaQiitaCombined_Nonzero_SpeciesShared_FamilyShared,
#      rep200FungiGenusShared_Nonzero,
#      metaQiitaCombined_Nonzero_SpeciesShared_GenusShared,
#      rep200FungiSpeciesShared_Nonzero,
#      metaQiitaCombined_Nonzero_SpeciesShared_SpeciesShared,
#      file = "Interim_data/tcga_data_taxa_levels_features_shared_with_Weizmann_13Sep21.RData")
load("Interim_data/tcga_data_taxa_levels_features_shared_with_Weizmann_13Sep21.RData")

# metaQiitaCombined_Nonzero_WithBamcounts, 
# rep200Data_WGS_RNA_Matched_Bacteria,
# rep200Data_WGS_RNA_Matched_Fungi,
load("Interim_data/metaQiitaCombined_Nonzero_WithBamcounts_and_Data_13Sep21.RData")

#----------------------------------------------------------#
# Load cibersort immune cell data and merge to mycobiome metadata
#----------------------------------------------------------#

# save(gID2TaxAll,
#      rep200Data_WGS_RNA_Matched,
#      metaQiitaWGS_RNA_HiSeq_Filt_Nonzero,
#      file = "data_for_mmvec_tables_7Aug21.RData")
# load("data_for_mmvec_tables_7Aug21.RData")

# Load relative abundance cibersort table from Thorsson et al. 2018 Immunity
# Only TCGA case IDs were provided for the mapping, so those are used here to
# map to TCGA mycobiome data
cibersortTCGA <- read.csv("Supporting_data/thorsson_2018_immunity_tcga_cibersort_edited_14Aug21.csv",
                          row.names = 1, stringsAsFactors = FALSE)
cibersortTCGAFormatted <- cibersortTCGA %>% rownames_to_column("tcga_case_id")
## Merge by tcga_case_id in metaQiitaCombined_Nonzero_SpeciesShared
head(metaQiitaCombined_Nonzero_SpeciesShared)

# Check for non-unique overlap
sum(metaQiitaCombined_Nonzero_SpeciesShared$tcga_case_id %in% rownames(cibersortTCGA)) # 11720 (NB: is less than 14546 --> missing samples)
sum(rownames(cibersortTCGA) %in% metaQiitaCombined_Nonzero_SpeciesShared$tcga_case_id) # 7143 --> repeats exist

# Create mycobiome-matched cibersort immune cell table via merging
mycobiomeImmuneCells <- metaQiitaCombined_Nonzero_SpeciesShared %>% rownames_to_column("sampleid") %>%
  select(sampleid, tcga_case_id) %>% droplevels() %>%
  left_join(cibersortTCGAFormatted, by = "tcga_case_id") %>%
  drop_na() %>%
  select(-tcga_case_id) %>%
  column_to_rownames("sampleid") -> mycobiomeImmuneCellsMerged
colnames(mycobiomeImmuneCellsMerged) <- gsub('([[:punct:]])|\\s+','',colnames(mycobiomeImmuneCellsMerged))
dim(mycobiomeImmuneCellsMerged) # 11720    22
sum(is.na(mycobiomeImmuneCellsMerged)) # 0

# Subset metadata
metaImmune <- droplevels(metaQiitaCombined_Nonzero_SpeciesShared[rownames(mycobiomeImmuneCellsMerged),])
metaImmunePT <- metaImmune %>% filter(sample_type == "Primary Tumor") %>% droplevels()
metaImmunePT_WGS <- metaImmunePT %>% filter(experimental_strategy == "WGS") %>% droplevels()
metaImmunePT_RNA <- metaImmunePT %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()
metaImmunePT_WGS_HMS <- metaImmunePT_WGS %>% filter(data_submitting_center_label == "Harvard Medical School") %>% droplevels()
metaImmunePT_RNA_UNC <- metaImmunePT_RNA %>% filter(data_submitting_center_label == "University of North Carolina") %>% droplevels()

# Subset immune cell data
mycobiomeImmuneCellsMergedPT <- mycobiomeImmuneCellsMerged[rownames(metaImmunePT),]
mycobiomeImmuneCellsMergedPT_WGS <- mycobiomeImmuneCellsMerged[rownames(metaImmunePT_WGS),]
mycobiomeImmuneCellsMergedPT_RNA <- mycobiomeImmuneCellsMerged[rownames(metaImmunePT_RNA),]
mycobiomeImmuneCellsMergedPT_WGS_HMS <- mycobiomeImmuneCellsMerged[rownames(metaImmunePT_WGS_HMS),]
mycobiomeImmuneCellsMergedPT_RNA_UNC <- mycobiomeImmuneCellsMerged[rownames(metaImmunePT_RNA_UNC),]

#-----------------------Create matching mycobiome data-----------------------#

rep200TaxSplit_Fungi_Paired_to_Weizmann <- read.csv("Supporting_data/rep200TaxSplit_Fungi_Paired_To_Weizmann_Final.csv", stringsAsFactors = FALSE, row.names = 1)
rep200TaxSplit_Fungi_Paired_to_Weizmann_formatted <- rep200TaxSplit_Fungi_Paired_to_Weizmann %>% 
  rownames_to_column("OGU") %>% column_to_rownames("species")
colnames(rep200FungiSpeciesShared_Nonzero) <- rep200TaxSplit_Fungi_Paired_to_Weizmann_formatted[colnames(rep200FungiSpeciesShared_Nonzero), "OGU"]

# rep200TaxSplit <- read.csv("Supporting_data/rep200_lineage_map_split.csv", stringsAsFactors = FALSE, row.names = 1)
# rep200Kingdoms <- read.csv("Supporting_data/rep200_gOTU_kingdom_mapping.csv", stringsAsFactors = FALSE)
# rep200Kingdoms_Fungi <- rep200Kingdoms[which(rep200Kingdoms$kingdom == "fungi"),]
# rep200Kingdoms_Bacteria <- rep200Kingdoms[which(rep200Kingdoms$kingdom == "bacteria"),]
# 
# rep200TaxSplit_Fungi <- rep200TaxSplit[rownames(rep200TaxSplit) %in% rep200Kingdoms_Fungi$genomeID,,drop=FALSE]
# dim(rep200TaxSplit_Fungi) # 320   7
# 
# rep200TaxSplit_Bacteria <- rep200TaxSplit[rownames(rep200TaxSplit) %in% rep200Kingdoms_Bacteria$genomeID,,drop=FALSE]
# dim(rep200TaxSplit_Bacteria) # 11080     7
# 
# fungiOGUs <- rownames(rep200TaxSplit_Fungi)
# bacteriaOGUs <- rownames(rep200TaxSplit_Bacteria)

# # Adapted from code written on Aug 7 2021 (tcga_mycobiome_kingdom_analyses.R)
# gID2TaxBacteria <- gID2TaxAll[grepl("^k__Bacteria",gID2TaxAll$taxa),]
# dim(gID2TaxBacteria) # 11080     2
# gID2TaxFungi <- gID2TaxAll[grepl("mycota",gID2TaxAll$taxa),]
# dim(gID2TaxFungi) # 318   2

rep200Data_Matched2ImmunePT_Bacteria <- rep200Data_WGS_RNA_Matched_Bacteria[rownames(metaImmunePT),]
dim(rep200Data_Matched2ImmunePT_Bacteria) # 9157 11071

rep200Data_Matched2ImmunePT_Fungi <- rep200FungiSpeciesShared_Nonzero[rownames(metaImmunePT),]
dim(rep200Data_Matched2ImmunePT_Fungi) # 9157   34

# Sanity check
all(rownames(metaImmunePT) == rownames(rep200Data_Matched2ImmunePT_Fungi)) # TRUE
all(rownames(rep200Data_Matched2ImmunePT_Bacteria) == rownames(rep200Data_Matched2ImmunePT_Fungi)) # TRUE

# Write metadata to text file
metaImmunePT %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/metadata_immune_all_primary_tumors.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of full data
rep200Data_Matched2ImmunePT_Bacteria_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_BIOM, biom_file = "MMvec_analysis/immune_rep200_counts_bacteria_all_primary_tumors.biom")
rep200Data_Matched2ImmunePT_Fungi_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi))
write_biom(rep200Data_Matched2ImmunePT_Fungi_BIOM, biom_file = "MMvec_analysis/immune_rep200_counts_fungi_all_primary_tumors.biom")
mycobiomeImmuneCellsMergedPT_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT))
write_biom(mycobiomeImmuneCellsMergedPT_BIOM, biom_file = "MMvec_analysis/immune_cibersort_rel_abund_all_TCGA_primary_tumors.biom")

## Subset to HMS
rep200Data_Matched2ImmunePT_Bacteria_HMS_PT <- rep200Data_Matched2ImmunePT_Bacteria[rownames(metaImmunePT_WGS_HMS),]
rep200Data_Matched2ImmunePT_Fungi_HMS_PT <- rep200Data_Matched2ImmunePT_Fungi[rownames(metaImmunePT_WGS_HMS),]

metaImmunePT_WGS_HMS %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/metadata_immune_WGS_Harvard_Primary_Tumor.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of full data
rep200Data_Matched2ImmunePT_Bacteria_HMS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria_HMS_PT))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_HMS_PT_BIOM, biom_file = "MMvec_analysis/immune_rep200_counts_bacteria_TCGA_Harvard_WGS_Primary_Tumor.biom")
rep200Data_Matched2ImmunePT_Fungi_HMS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_HMS_PT))
write_biom(rep200Data_Matched2ImmunePT_Fungi_HMS_PT_BIOM, biom_file = "MMvec_analysis/immune_rep200_counts_fungi_TCGA_Harvard_WGS_Primary_Tumor.biom")
mycobiomeImmuneCellsMergedPT_WGS_HMS_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT_WGS_HMS))
write_biom(mycobiomeImmuneCellsMergedPT_WGS_HMS_BIOM, biom_file = "MMvec_analysis/immune_cibersort_rel_abund_TCGA_Harvard_WGS_Primary_Tumor.biom")

## Subset to UNC
rep200Data_Matched2ImmunePT_Bacteria_UNC_PT <- rep200Data_Matched2ImmunePT_Bacteria[rownames(metaImmunePT_RNA_UNC),]
rep200Data_Matched2ImmunePT_Fungi_UNC_PT <- rep200Data_Matched2ImmunePT_Fungi[rownames(metaImmunePT_RNA_UNC),]

metaImmunePT_RNA_UNC %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/metadata_immune_RNA_UNC_Primary_Tumor.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of full data
rep200Data_Matched2ImmunePT_Bacteria_UNC_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria_UNC_PT))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_UNC_PT_BIOM, biom_file = "MMvec_analysis/immune_rep200_counts_bacteria_TCGA_UNC_RNA_Primary_Tumor.biom")
rep200Data_Matched2ImmunePT_Fungi_UNC_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_UNC_PT))
write_biom(rep200Data_Matched2ImmunePT_Fungi_UNC_PT_BIOM, biom_file = "MMvec_analysis/immune_rep200_counts_fungi_TCGA_UNC_RNA_Primary_Tumor.biom")
mycobiomeImmuneCellsMergedPT_RNA_UNC_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT_RNA_UNC))
write_biom(mycobiomeImmuneCellsMergedPT_RNA_UNC_BIOM, biom_file = "MMvec_analysis/immune_cibersort_rel_abund_TCGA_UNC_RNA_Primary_Tumor.biom")

## Subset to all WGS
rep200Data_Matched2ImmunePT_Bacteria_WGS_PT <- rep200Data_Matched2ImmunePT_Bacteria[rownames(metaImmunePT_WGS),]
rep200Data_Matched2ImmunePT_Fungi_WGS_PT <- rep200Data_Matched2ImmunePT_Fungi[rownames(metaImmunePT_WGS),]

metaImmunePT_WGS %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/metadata_immune_WGS_AllSeqCenters_Primary_Tumor.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of full data
rep200Data_Matched2ImmunePT_Bacteria_WGS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria_WGS_PT))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_WGS_PT_BIOM, biom_file = "MMvec_analysis/immune_rep200_counts_bacteria_TCGA_AllSeqCenters_WGS_Primary_Tumor.biom")
rep200Data_Matched2ImmunePT_Fungi_WGS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_WGS_PT))
write_biom(rep200Data_Matched2ImmunePT_Fungi_WGS_PT_BIOM, biom_file = "MMvec_analysis/immune_rep200_counts_fungi_TCGA_AllSeqCenters_WGS_Primary_Tumor.biom")
mycobiomeImmuneCellsMergedPT_WGS_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT_WGS))
write_biom(mycobiomeImmuneCellsMergedPT_WGS_BIOM, biom_file = "MMvec_analysis/immune_cibersort_rel_abund_TCGA_AllSeqCenters_WGS_Primary_Tumor.biom")

