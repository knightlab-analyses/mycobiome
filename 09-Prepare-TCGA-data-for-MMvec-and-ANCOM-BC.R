#-----------------------------------------------------------------------------
# 09-Prepare-TCGA-data-for-MMvec-and-ANCOM-BC.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Load data with Weizmann overlapping features
# - Format data for MMvec to examine co-occurences of (bacteria and fungi) and (immune cells and fungi)
# - Perform ANCOM-BC on subsetted data (without batch correction) for tumor vs. NAT analyses
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
require(ggrepel)

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

rep200Data_Matched2ImmunePT_Fungi_species <- rep200FungiSpeciesShared_Nonzero[rownames(metaImmunePT),]
dim(rep200Data_Matched2ImmunePT_Fungi_species) # 9157   34

# Sanity check
all(rownames(metaImmunePT) == rownames(rep200Data_Matched2ImmunePT_Fungi_species)) # TRUE
all(rownames(rep200Data_Matched2ImmunePT_Bacteria) == rownames(rep200Data_Matched2ImmunePT_Fungi_species)) # TRUE

#-----------------------Write biom tables and metadata files-----------------------#

# Write metadata to text file
metaImmunePT %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/species_intersected_with_WIS/metadata_immune_all_primary_tumors.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of full data
rep200Data_Matched2ImmunePT_Bacteria_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_BIOM, biom_file = "MMvec_analysis/species_intersected_with_WIS/immune_rep200_counts_bacteria_all_primary_tumors.biom")
rep200Data_Matched2ImmunePT_Fungi_species_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_species))
write_biom(rep200Data_Matched2ImmunePT_Fungi_species_BIOM, biom_file = "MMvec_analysis/species_intersected_with_WIS/immune_rep200_counts_fungi_all_primary_tumors.biom")
mycobiomeImmuneCellsMergedPT_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT))
write_biom(mycobiomeImmuneCellsMergedPT_BIOM, biom_file = "MMvec_analysis/species_intersected_with_WIS/immune_cibersort_rel_abund_all_TCGA_primary_tumors.biom")

## Subset to HMS
rep200Data_Matched2ImmunePT_Bacteria_HMS_PT <- rep200Data_Matched2ImmunePT_Bacteria[rownames(metaImmunePT_WGS_HMS),]
rep200Data_Matched2ImmunePT_Fungi_species_HMS_PT <- rep200Data_Matched2ImmunePT_Fungi_species[rownames(metaImmunePT_WGS_HMS),]

metaImmunePT_WGS_HMS %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/species_intersected_with_WIS/metadata_immune_WGS_Harvard_Primary_Tumor.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of full data
rep200Data_Matched2ImmunePT_Bacteria_HMS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria_HMS_PT))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_HMS_PT_BIOM, biom_file = "MMvec_analysis/species_intersected_with_WIS/immune_rep200_counts_bacteria_TCGA_Harvard_WGS_Primary_Tumor.biom")
rep200Data_Matched2ImmunePT_Fungi_species_HMS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_species_HMS_PT))
write_biom(rep200Data_Matched2ImmunePT_Fungi_species_HMS_PT_BIOM, biom_file = "MMvec_analysis/species_intersected_with_WIS/immune_rep200_counts_fungi_TCGA_Harvard_WGS_Primary_Tumor.biom")
mycobiomeImmuneCellsMergedPT_WGS_HMS_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT_WGS_HMS))
write_biom(mycobiomeImmuneCellsMergedPT_WGS_HMS_BIOM, biom_file = "MMvec_analysis/species_intersected_with_WIS/immune_cibersort_rel_abund_TCGA_Harvard_WGS_Primary_Tumor.biom")

## Subset to UNC
rep200Data_Matched2ImmunePT_Bacteria_UNC_PT <- rep200Data_Matched2ImmunePT_Bacteria[rownames(metaImmunePT_RNA_UNC),]
rep200Data_Matched2ImmunePT_Fungi_species_UNC_PT <- rep200Data_Matched2ImmunePT_Fungi_species[rownames(metaImmunePT_RNA_UNC),]

metaImmunePT_RNA_UNC %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/species_intersected_with_WIS/metadata_immune_RNA_UNC_Primary_Tumor.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of full data
rep200Data_Matched2ImmunePT_Bacteria_UNC_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria_UNC_PT))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_UNC_PT_BIOM, biom_file = "MMvec_analysis/species_intersected_with_WIS/immune_rep200_counts_bacteria_TCGA_UNC_RNA_Primary_Tumor.biom")
rep200Data_Matched2ImmunePT_Fungi_species_UNC_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_species_UNC_PT))
write_biom(rep200Data_Matched2ImmunePT_Fungi_species_UNC_PT_BIOM, biom_file = "MMvec_analysis/species_intersected_with_WIS/immune_rep200_counts_fungi_TCGA_UNC_RNA_Primary_Tumor.biom")
mycobiomeImmuneCellsMergedPT_RNA_UNC_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT_RNA_UNC))
write_biom(mycobiomeImmuneCellsMergedPT_RNA_UNC_BIOM, biom_file = "MMvec_analysis/species_intersected_with_WIS/immune_cibersort_rel_abund_TCGA_UNC_RNA_Primary_Tumor.biom")

## Subset to all WGS
rep200Data_Matched2ImmunePT_Bacteria_WGS_PT <- rep200Data_Matched2ImmunePT_Bacteria[rownames(metaImmunePT_WGS),]
rep200Data_Matched2ImmunePT_Fungi_species_WGS_PT <- rep200Data_Matched2ImmunePT_Fungi_species[rownames(metaImmunePT_WGS),]

metaImmunePT_WGS %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/species_intersected_with_WIS/metadata_immune_WGS_AllSeqCenters_Primary_Tumor.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of full data
rep200Data_Matched2ImmunePT_Bacteria_WGS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria_WGS_PT))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_WGS_PT_BIOM, biom_file = "MMvec_analysis/species_intersected_with_WIS/immune_rep200_counts_bacteria_TCGA_AllSeqCenters_WGS_Primary_Tumor.biom")
rep200Data_Matched2ImmunePT_Fungi_species_WGS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_species_WGS_PT))
write_biom(rep200Data_Matched2ImmunePT_Fungi_species_WGS_PT_BIOM, biom_file = "MMvec_analysis/species_intersected_with_WIS/immune_rep200_counts_fungi_TCGA_AllSeqCenters_WGS_Primary_Tumor.biom")
mycobiomeImmuneCellsMergedPT_WGS_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT_WGS))
write_biom(mycobiomeImmuneCellsMergedPT_WGS_BIOM, biom_file = "MMvec_analysis/species_intersected_with_WIS/immune_cibersort_rel_abund_TCGA_AllSeqCenters_WGS_Primary_Tumor.biom")

#--------------------------------Write genus level files--------------------------------#
# There's one sample that is not shared between the metaImmunePT and rep200FungiGenusShared_Nonzero (13767.58cfa833e4b0c9d6adf6dc00),
# so they need to be subsetted

rep200TaxSplit <- read.csv("Supporting_data/rep200_lineage_map_split.csv", stringsAsFactors = FALSE, row.names = 1)
rep200Kingdoms <- read.csv("Supporting_data/rep200_gOTU_kingdom_mapping.csv", stringsAsFactors = FALSE)
rep200Kingdoms_Fungi <- rep200Kingdoms[which(rep200Kingdoms$kingdom == "fungi"),]
rep200Kingdoms_Bacteria <- rep200Kingdoms[which(rep200Kingdoms$kingdom == "bacteria"),]

# rep200TaxSplit_Fungi <- rep200TaxSplit[rownames(rep200TaxSplit) %in% rep200Kingdoms_Fungi$genomeID,,drop=FALSE]
# dim(rep200TaxSplit_Fungi) # 320   7

rep200TaxSplit_Bacteria <- rep200TaxSplit[rownames(rep200TaxSplit) %in% rep200Kingdoms_Bacteria$genomeID,,drop=FALSE]
dim(rep200TaxSplit_Bacteria) # 11080     7
## Subset
genusSamplesPT <- intersect(rownames(metaImmunePT), rownames(rep200FungiGenusShared_Nonzero))
metaImmunePT_genus <- droplevels(metaImmunePT[genusSamplesPT,])
rep200Data_Matched2ImmunePT_Fungi_genus <- rep200FungiGenusShared_Nonzero[genusSamplesPT,]
rep200Data_Matched2ImmunePT_Bacteria_Filt <- rep200Data_WGS_RNA_Matched_Bacteria[genusSamplesPT,]
dim(rep200Data_Matched2ImmunePT_Fungi_genus) # 9156   54
dim(rep200Data_Matched2ImmunePT_Bacteria_Filt) # 9156 11071

# Build phyloseq object
ps_rep200Data_Matched2ImmunePT_Bacteria_Filt <- phyloseq(otu_table(rep200Data_Matched2ImmunePT_Bacteria_Filt, taxa_are_rows = FALSE), 
                             tax_table(as.matrix(rep200TaxSplit_Bacteria)), 
                             sample_data(metaImmunePT_genus))

## Aggregate counts
# ps_rep200Data_Matched2ImmunePT_Bacteria_Filt_phylum = aggregate_taxa(ps_rep200Data_Matched2ImmunePT_Bacteria_Filt, "Phylum")
# ps_rep200Data_Matched2ImmunePT_Bacteria_Filt_class = aggregate_taxa(ps_rep200Data_Matched2ImmunePT_Bacteria_Filt, "Class")
# ps_rep200Data_Matched2ImmunePT_Bacteria_Filt_order = aggregate_taxa(ps_rep200Data_Matched2ImmunePT_Bacteria_Filt, "Order")
# ps_rep200Data_Matched2ImmunePT_Bacteria_Filt_family = aggregate_taxa(ps_rep200Data_Matched2ImmunePT_Bacteria_Filt, "Family")
ps_rep200Data_Matched2ImmunePT_Bacteria_Filt_genus = aggregate_taxa(ps_rep200Data_Matched2ImmunePT_Bacteria_Filt, "Genus")

rep200Data_Matched2ImmunePT_Bacteria_Filt_genus <- data.frame(t(otu_table(ps_rep200Data_Matched2ImmunePT_Bacteria_Filt_genus)))
dim(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus) # 9156 2675

#-----------------------Write biom tables and metadata files-----------------------#
metaImmunePT_genus_WGS <- metaImmunePT_genus %>% filter(experimental_strategy == "WGS") %>% droplevels()
metaImmunePT_genus_RNA <- metaImmunePT_genus %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()
metaImmunePT_genus_WGS_HMS <- metaImmunePT_genus_WGS %>% filter(data_submitting_center_label == "Harvard Medical School") %>% droplevels()
metaImmunePT_genus_RNA_UNC <- metaImmunePT_genus_RNA %>% filter(data_submitting_center_label == "University of North Carolina") %>% droplevels()

# Write metadata to text file
metaImmunePT_genus %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/genus_intersected_with_WIS/metadata_immune_all_primary_tumors.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of full data
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_bacteria_all_primary_tumors.biom")
rep200Data_Matched2ImmunePT_Fungi_genus_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_genus))
write_biom(rep200Data_Matched2ImmunePT_Fungi_genus_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_fungi_all_primary_tumors.biom")
mycobiomeImmuneCellsMergedPT_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT))
write_biom(mycobiomeImmuneCellsMergedPT_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_cibersort_rel_abund_all_TCGA_primary_tumors.biom")

## Subset to HMS
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_HMS_PT <- rep200Data_Matched2ImmunePT_Bacteria_Filt_genus[rownames(metaImmunePT_genus_WGS_HMS),]
rep200Data_Matched2ImmunePT_Fungi_genus_HMS_PT <- rep200Data_Matched2ImmunePT_Fungi_genus[rownames(metaImmunePT_genus_WGS_HMS),]

metaImmunePT_genus_WGS_HMS %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/genus_intersected_with_WIS/metadata_immune_WGS_Harvard_Primary_Tumor.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of HMS data
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_HMS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_HMS_PT))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_HMS_PT_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_bacteria_TCGA_Harvard_WGS_Primary_Tumor.biom")
rep200Data_Matched2ImmunePT_Fungi_genus_HMS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_genus_HMS_PT))
write_biom(rep200Data_Matched2ImmunePT_Fungi_genus_HMS_PT_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_fungi_TCGA_Harvard_WGS_Primary_Tumor.biom")
mycobiomeImmuneCellsMergedPT_WGS_HMS_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT_WGS_HMS))
write_biom(mycobiomeImmuneCellsMergedPT_WGS_HMS_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_cibersort_rel_abund_TCGA_Harvard_WGS_Primary_Tumor.biom")

## Subset to UNC
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_UNC_PT <- rep200Data_Matched2ImmunePT_Bacteria_Filt_genus[rownames(metaImmunePT_genus_RNA_UNC),]
rep200Data_Matched2ImmunePT_Fungi_genus_UNC_PT <- rep200Data_Matched2ImmunePT_Fungi_genus[rownames(metaImmunePT_genus_RNA_UNC),]

metaImmunePT_genus_RNA_UNC %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/genus_intersected_with_WIS/metadata_immune_RNA_UNC_Primary_Tumor.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of UNC data
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_UNC_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_UNC_PT))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_UNC_PT_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_bacteria_TCGA_UNC_RNA_Primary_Tumor.biom")
rep200Data_Matched2ImmunePT_Fungi_genus_UNC_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_genus_UNC_PT))
write_biom(rep200Data_Matched2ImmunePT_Fungi_genus_UNC_PT_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_fungi_TCGA_UNC_RNA_Primary_Tumor.biom")
mycobiomeImmuneCellsMergedPT_RNA_UNC_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT_RNA_UNC))
write_biom(mycobiomeImmuneCellsMergedPT_RNA_UNC_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_cibersort_rel_abund_TCGA_UNC_RNA_Primary_Tumor.biom")

## Subset to all WGS
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_WGS_PT <- rep200Data_Matched2ImmunePT_Bacteria_Filt_genus[rownames(metaImmunePT_genus_WGS),]
rep200Data_Matched2ImmunePT_Fungi_genus_WGS_PT <- rep200Data_Matched2ImmunePT_Fungi_genus[rownames(metaImmunePT_genus_WGS),]

metaImmunePT_genus_WGS %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/genus_intersected_with_WIS/metadata_immune_WGS_AllSeqCenters_Primary_Tumor.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of WGS data
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_WGS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_WGS_PT))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_WGS_PT_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_bacteria_TCGA_AllSeqCenters_WGS_Primary_Tumor.biom")
rep200Data_Matched2ImmunePT_Fungi_genus_WGS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_genus_WGS_PT))
write_biom(rep200Data_Matched2ImmunePT_Fungi_genus_WGS_PT_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_fungi_TCGA_AllSeqCenters_WGS_Primary_Tumor.biom")
mycobiomeImmuneCellsMergedPT_WGS_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT_WGS))
write_biom(mycobiomeImmuneCellsMergedPT_WGS_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_cibersort_rel_abund_TCGA_AllSeqCenters_WGS_Primary_Tumor.biom")


#----------------------------------------------------------#
# ANCOM-BC tumor vs NAT
# - Subset data to individual seq center and data type
#----------------------------------------------------------#

metaQiitaCombined_Nonzero_TvsNAT <- metaQiitaCombined_Nonzero %>% filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal")) %>% droplevels()
metaQiitaCombined_Nonzero_TvsNAT$sample_type <- relevel(metaQiitaCombined_Nonzero_TvsNAT$sample_type, ref = "Solid Tissue Normal")
rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero_TvsNAT <- rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero[rownames(metaQiitaCombined_Nonzero_TvsNAT),]

psFungiTCGADecontam_TvsNAT <- phyloseq(otu_table(rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero_TvsNAT, taxa_are_rows = FALSE), 
                             tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), 
                             sample_data(metaQiitaCombined_Nonzero_TvsNAT))
# psFungiTCGADecontam_TvsNAT <- subset_samples(psFungiTCGADecontam, sample_type %in% c("Primary Tumor", "Solid Tissue Normal"))
# Subset to seq center
psFungiTCGADecontam_TvsNAT_HMS <- subset_samples(psFungiTCGADecontam_TvsNAT, data_submitting_center_label == "Harvard Medical School")
psFungiTCGADecontam_TvsNAT_BCM <- subset_samples(psFungiTCGADecontam_TvsNAT, data_submitting_center_label == "Baylor College of Medicine")
psFungiTCGADecontam_TvsNAT_MDA <- subset_samples(psFungiTCGADecontam_TvsNAT, data_submitting_center_label == "MD Anderson - Institute for Applied Cancer Science")
psFungiTCGADecontam_TvsNAT_WashU <- subset_samples(psFungiTCGADecontam_TvsNAT, data_submitting_center_label == "Washington University School of Medicine")
psFungiTCGADecontam_TvsNAT_UNC <- subset_samples(psFungiTCGADecontam_TvsNAT, data_submitting_center_label == "University of North Carolina")
psFungiTCGADecontam_TvsNAT_CMS <- subset_samples(psFungiTCGADecontam_TvsNAT, data_submitting_center_label == "Canada's Michael Smith Genome Sciences Centre")
psFungiTCGADecontam_TvsNAT_Broad <- subset_samples(psFungiTCGADecontam_TvsNAT, data_submitting_center_label == "Broad Institute of MIT and Harvard")
psFungiTCGADecontam_TvsNAT_Broad_WGS <- subset_samples(psFungiTCGADecontam_TvsNAT_Broad, experimental_strategy == "WGS")
psFungiTCGADecontam_TvsNAT_Broad_RNA <- subset_samples(psFungiTCGADecontam_TvsNAT_Broad, experimental_strategy == "RNA-Seq")

## Write function

runAncomBC_TvsNAT <- function(psSeqCenter, ancombcLibCut = 1000){
  metaData <- data.frame(sample_data(psSeqCenter))
  SeqCenter <- as.character(metaData$data_submitting_center_label[1])
  SeqCenterFormatted <- gsub('([[:punct:]])|\\s+','',SeqCenter)
  cancerTypes <- as.character(unique(metaData$disease_type))
  
  
  for(ii in seq_along(cancerTypes)){
    Dz <- cancerTypes[ii]
    DzFormatted <- gsub('([[:punct:]])|\\s+','',Dz)
    
    metaDataFilt <- metaData %>% filter(disease_type == Dz)
    
    # SKIP CANCER TYPES THAT ONLY HAVE ONE CLASS OF A COMPARISON
    if(length(table(metaDataFilt$sample_type)) < 2){next}
    
    # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF 10 SAMPLES IN EITHER CLASS
    if(any(table(metaDataFilt$sample_type) < 10)){next}
    
    # If sufficient samples, then:
    print(Dz)
    print(sprintf("Number of samples (PT|NAT): %d | %d", 
                  unname(table(metaDataFilt$sample_type)[2]),
                  unname(table(metaDataFilt$sample_type)[1])))
    filtBoolean <- (sample_data(psSeqCenter)[["disease_type"]] == Dz)
    psFungiTCGADecontam_TvsNAT_X <- prune_samples(filtBoolean, psSeqCenter)
    # print(psFungiTCGADecontam_TvsNAT_X)
    print(sprintf("Read count distribution:"))
    print(summary(rowSums(otu_table(psFungiTCGADecontam_TvsNAT_X)))) # Sample read distribution
    print(sprintf("Now running ANCOM-BC..."))
    ancombc_Fungi_TvsNAT_X <- ancombc(phyloseq = psFungiTCGADecontam_TvsNAT_X, 
                                      formula = "sample_type", p_adj_method = "BH", zero_cut = 0.999, 
                                      lib_cut = ancombcLibCut, 
                                      # group = "sample_type", struc_zero = FALSE, neg_lb = FALSE,
                                      tol = 1e-5, max_iter = 100, conserve = FALSE, alpha = 0.05, global = FALSE)
    print(sprintf("Plotting and saving data..."))
    ancom_res_df_Fungi_X <- data.frame(
      beta = unlist(ancombc_Fungi_TvsNAT_X$res$beta),
      se = unlist(ancombc_Fungi_TvsNAT_X$res$se),
      W = unlist(ancombc_Fungi_TvsNAT_X$res$W),
      p_val = unlist(ancombc_Fungi_TvsNAT_X$res$p_val),
      q_val = unlist(ancombc_Fungi_TvsNAT_X$res$q_val),
      diff_abn = ifelse(unlist(ancombc_Fungi_TvsNAT_X$res$q_val)<=0.1, yes = TRUE, no = FALSE),
      genus = rep200TaxSplit_Fungi_Paired_to_Weizmann[row.names(ancombc_Fungi_TvsNAT_X$res$beta),"genus"],
      species = rep200TaxSplit_Fungi_Paired_to_Weizmann[row.names(ancombc_Fungi_TvsNAT_X$res$beta),"species"],
      OGUs = row.names(ancombc_Fungi_TvsNAT_X$res$beta),
      row.names = row.names(ancombc_Fungi_TvsNAT_X$res$beta))
    ancom_res_df_Fungi_X$diff_name_flag <- ancom_res_df_Fungi_X$diff_abn + 0
    ancom_res_df_Fungi_X_sorted <- ancom_res_df_Fungi_X[order(ancom_res_df_Fungi_X$q_val),]
    
    ancom_res_df_Fungi_X_sorted$diff_label <- ""
    ancom_res_df_Fungi_X_sorted$diff_label[1:5] <- ifelse(ancom_res_df_Fungi_X_sorted$diff_abn, 
                                                          yes = paste0(ancom_res_df_Fungi_X_sorted$OGUs,"\n(",ancom_res_df_Fungi_X_sorted$genus,")"),
                                                          no = "")[1:5]
    
    print(sprintf("Number of differentially abundant fungi: %d", sum(ancom_res_df_Fungi_X_sorted$diff_name_flag)))
    cat("\n")
    
    plotFilePath <- "Figures/Supplementary_Figures/"
    ancom_res_df_Fungi_X_sorted %>%
      ggplot(aes(x = beta, y = -log10(p_val), color = diff_abn, label = diff_label)) + geom_point() +
      theme_bw() + geom_hline(yintercept=-log10(0.05), col="black", linetype='dashed') +
      geom_vline(xintercept=c(-0.2, 0.2), col="black", linetype='dashed') + 
      theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + coord_fixed() +
      scale_color_aaas() +
      labs(x = "log-fold change\n(ANCOM-BC beta)", y = "-Log10(p-value)", color = "Differentially\nabundant (q<=0.1)", title = paste(SeqCenter, Dz, sep = "\n")) +
      geom_label_repel(size = 1.5, box.padding = 3, point.padding = 1e-04, show.legend = FALSE, color = "black", max.overlaps = 10) +
      # geom_label_repel(aes(label = ifelse(diff_name_flag == 1, as.character(species), "")), size = 2.5, show.legend = FALSE) +
      ggsave(filename = paste0(plotFilePath, "ancombc_TvsNAT_", SeqCenterFormatted,"_",DzFormatted,".pdf"), dpi = "retina",
             width = 6, height = 4, units = "in")
    # Write data to file
    dataFilePath <- "Figures_data/Supplementary_Figures/"
    ancom_res_df_Fungi_X_sorted %>% write.csv(file = paste0(dataFilePath, "ancombc_TvsNAT_", SeqCenterFormatted,"_",DzFormatted,".csv"))
  }
}
# runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_BCM)

runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_HMS)
runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_BCM)
runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_MDA)
runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_WashU)
runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_Broad_WGS)
runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_UNC, ancombcLibCut = 0)
runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_CMS, ancombcLibCut = 0)
runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_Broad_RNA, ancombcLibCut = 0)
