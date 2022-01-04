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

# save(snmDataOGUFungiDecontamV2,
#      vdge_dataE_DecontamV2,
#      rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero,
#      metaQiitaCombined_Nonzero_DecontamV2,
#      file = "Interim_data/snmDataFungi_DecontamV2_13Oct21.RData")
load("Interim_data/snmDataFungi_DecontamV2_13Oct21.RData")

# save(rep200FungiPhylumShared_Nonzero,
#      metaQiitaCombined_Nonzero_PhylumShared,
#      rep200FungiClassShared_Nonzero,
#      metaQiitaCombined_Nonzero_ClassShared,
#      rep200FungiOrderShared_Nonzero,
#      metaQiitaCombined_Nonzero_OrderShared,
#      rep200FungiFamilyShared_Nonzero,
#      metaQiitaCombined_Nonzero_FamilyShared,
#      rep200FungiGenusShared_Nonzero,
#      metaQiitaCombined_Nonzero_GenusShared,
#      rep200FungiSpeciesShared_Nonzero,
#      metaQiitaCombined_Nonzero_SpeciesShared,
#      file = "Interim_data/tcga_data_taxa_levels_features_shared_with_Weizmann_14Oct21.RData")
load("Interim_data/tcga_data_taxa_levels_features_shared_with_Weizmann_14Oct21.RData")

# # metaQiitaCombined_Nonzero_WithBamcounts, 
# # rep200Data_WGS_RNA_Matched_Bacteria,
# # rep200Data_WGS_RNA_Matched_Fungi,
# load("Interim_data/metaQiitaCombined_Nonzero_WithBamcounts_and_Data_13Sep21.RData")

# save(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts, 
#      rep200Data_WGS_RNA_Matched,
#      rep200Data_WGS_RNA_Matched_Bacteria,
#      rep200Data_WGS_RNA_Matched_Fungi,
#      file = "Interim_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_29Sep21.RData")
load("Interim_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_29Sep21.RData")

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
dim(metaQiitaCombined_Nonzero_SpeciesShared) # 12750    41

# Check for non-unique overlap
sum(metaQiitaCombined_Nonzero_SpeciesShared$tcga_case_id %in% rownames(cibersortTCGA)) # 11723 (NB: is less than 14546 --> missing samples)
sum(rownames(cibersortTCGA) %in% metaQiitaCombined_Nonzero_SpeciesShared$tcga_case_id) # 7145 --> repeats exist

# Create mycobiome-matched cibersort immune cell table via merging
mycobiomeImmuneCells <- metaQiitaCombined_Nonzero_SpeciesShared %>% rownames_to_column("sampleid") %>%
  select(sampleid, tcga_case_id) %>% droplevels() %>%
  left_join(cibersortTCGAFormatted, by = "tcga_case_id") %>%
  drop_na() %>%
  select(-tcga_case_id) %>%
  column_to_rownames("sampleid") -> mycobiomeImmuneCellsMerged
colnames(mycobiomeImmuneCellsMerged) <- gsub('([[:punct:]])|\\s+','',colnames(mycobiomeImmuneCellsMerged))
dim(mycobiomeImmuneCellsMerged) # 11723    22
sum(is.na(mycobiomeImmuneCellsMerged)) # 0

# Subset metadata
metaImmune <- droplevels(metaQiitaCombined_Nonzero_SpeciesShared[rownames(mycobiomeImmuneCellsMerged),])
metaImmunePT <- metaImmune %>% filter(sample_type == "Primary Tumor") %>% droplevels()
metaImmuneBDN <- metaImmune %>% filter(sample_type == "Blood Derived Normal") %>% droplevels()
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
# Load WIS matching fungi
rep200TaxSplit_Fungi_Paired_to_Weizmann <- read.csv("Supporting_data/rep200TaxSplit_Fungi_Paired_To_Weizmann_Final.csv", stringsAsFactors = FALSE, row.names = 1)
rep200TaxSplit_Fungi_Paired_to_Weizmann_formatted <- rep200TaxSplit_Fungi_Paired_to_Weizmann %>% 
  rownames_to_column("OGU") %>% column_to_rownames("species")
colnames(rep200FungiSpeciesShared_Nonzero) <- rep200TaxSplit_Fungi_Paired_to_Weizmann_formatted[colnames(rep200FungiSpeciesShared_Nonzero), "OGU"]
# Load WIS matching bacteria
load("Interim_data/shared_bacterial_features_at_each_taxa_level_15Oct21.RData")

rep200Data_Matched2ImmunePT_Bacteria <- rep200Data_WGS_RNA_Matched_Bacteria[rownames(metaImmunePT),sharedSpeciesBacteria$intersectedOGUs]
dim(rep200Data_Matched2ImmunePT_Bacteria) # 9159 267

rep200Data_Matched2ImmunePT_Fungi_species <- rep200FungiSpeciesShared_Nonzero[rownames(metaImmunePT),]
dim(rep200Data_Matched2ImmunePT_Fungi_species) # 9159   34

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

# Save biom files of HMS data
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

# Save biom files of UNC data
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

# Save biom files of WGS data
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
dim(rep200Data_Matched2ImmunePT_Fungi_genus) # 9158   54
dim(rep200Data_Matched2ImmunePT_Bacteria_Filt) # 9158 11071

all(rownames(rep200Data_Matched2ImmunePT_Fungi_genus) == rownames(rep200Data_Matched2ImmunePT_Bacteria_Filt)) # TRUE

# Build phyloseq object
ps_rep200Data_Matched2ImmunePT_Bacteria_Filt <- phyloseq(otu_table(rep200Data_Matched2ImmunePT_Bacteria_Filt, taxa_are_rows = FALSE), 
                             tax_table(as.matrix(rep200TaxSplit_Bacteria)), 
                             sample_data(metaImmunePT_genus))

## Aggregate counts
ps_rep200Data_Matched2ImmunePT_Bacteria_Filt_genus = aggregate_taxa(ps_rep200Data_Matched2ImmunePT_Bacteria_Filt, "Genus")

rep200Data_Matched2ImmunePT_Bacteria_Filt_genus <- data.frame(t(otu_table(ps_rep200Data_Matched2ImmunePT_Bacteria_Filt_genus)))
colnames(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus) <- gsub("^g__","",colnames(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus))
dim(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus) # 9158 2675

# Load WIS shared bacteria genera and subset
load("Interim_data/shared_bacterial_features_at_each_taxa_level_15Oct21.RData")
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared <- rep200Data_Matched2ImmunePT_Bacteria_Filt_genus[,colnames(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus) %in% unique(sharedGenusBacteria$intersectedTaxa)]
dim(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared) # 9158  196

## Save data
dfSharedSave <- data.frame(domain = c(rep("fungi",length(colnames(rep200Data_Matched2ImmunePT_Fungi_genus))),
                                      rep("bacteria",length(unique(sharedGenusBacteria$intersectedTaxa)))),
                           shared_genera = c(colnames(rep200Data_Matched2ImmunePT_Fungi_genus),
                                             unique(sharedGenusBacteria$intersectedTaxa)))
dfSharedSave %>% write.csv("Figures/Supplementary_Figures/mmvec_shared_genera_fungi_bacteria.csv", row.names = FALSE)

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
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_bacteria_all_primary_tumors.biom")
rep200Data_Matched2ImmunePT_Fungi_genus_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_genus))
write_biom(rep200Data_Matched2ImmunePT_Fungi_genus_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_fungi_all_primary_tumors.biom")
mycobiomeImmuneCellsMergedPT_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT))
write_biom(mycobiomeImmuneCellsMergedPT_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_cibersort_rel_abund_all_TCGA_primary_tumors.biom")

## Subset to HMS
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared_HMS_PT <- rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared[rownames(metaImmunePT_genus_WGS_HMS),]
rep200Data_Matched2ImmunePT_Fungi_genus_HMS_PT <- rep200Data_Matched2ImmunePT_Fungi_genus[rownames(metaImmunePT_genus_WGS_HMS),]

metaImmunePT_genus_WGS_HMS %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/genus_intersected_with_WIS/metadata_immune_WGS_Harvard_Primary_Tumor.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of HMS data
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared_HMS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared_HMS_PT))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared_HMS_PT_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_bacteria_TCGA_Harvard_WGS_Primary_Tumor.biom")
rep200Data_Matched2ImmunePT_Fungi_genus_HMS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_genus_HMS_PT))
write_biom(rep200Data_Matched2ImmunePT_Fungi_genus_HMS_PT_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_fungi_TCGA_Harvard_WGS_Primary_Tumor.biom")
mycobiomeImmuneCellsMergedPT_WGS_HMS_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT_WGS_HMS))
write_biom(mycobiomeImmuneCellsMergedPT_WGS_HMS_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_cibersort_rel_abund_TCGA_Harvard_WGS_Primary_Tumor.biom")

## Subset to UNC
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared_UNC_PT <- rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared[rownames(metaImmunePT_genus_RNA_UNC),]
rep200Data_Matched2ImmunePT_Fungi_genus_UNC_PT <- rep200Data_Matched2ImmunePT_Fungi_genus[rownames(metaImmunePT_genus_RNA_UNC),]

metaImmunePT_genus_RNA_UNC %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/genus_intersected_with_WIS/metadata_immune_RNA_UNC_Primary_Tumor.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of UNC data
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared_UNC_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared_UNC_PT))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared_UNC_PT_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_bacteria_TCGA_UNC_RNA_Primary_Tumor.biom")
rep200Data_Matched2ImmunePT_Fungi_genus_UNC_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_genus_UNC_PT))
write_biom(rep200Data_Matched2ImmunePT_Fungi_genus_UNC_PT_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_fungi_TCGA_UNC_RNA_Primary_Tumor.biom")
mycobiomeImmuneCellsMergedPT_RNA_UNC_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT_RNA_UNC))
write_biom(mycobiomeImmuneCellsMergedPT_RNA_UNC_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_cibersort_rel_abund_TCGA_UNC_RNA_Primary_Tumor.biom")

## Subset to all WGS
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared_WGS_PT <- rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared[rownames(metaImmunePT_genus_WGS),]
rep200Data_Matched2ImmunePT_Fungi_genus_WGS_PT <- rep200Data_Matched2ImmunePT_Fungi_genus[rownames(metaImmunePT_genus_WGS),]

metaImmunePT_genus_WGS %>% rownames_to_column("sampleid") %>%
  write.table(file = "MMvec_analysis/genus_intersected_with_WIS/metadata_immune_WGS_AllSeqCenters_Primary_Tumor.txt", 
              quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

# Save biom files of WGS data
rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared_WGS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared_WGS_PT))
write_biom(rep200Data_Matched2ImmunePT_Bacteria_Filt_genus_shared_WGS_PT_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_bacteria_TCGA_AllSeqCenters_WGS_Primary_Tumor.biom")
rep200Data_Matched2ImmunePT_Fungi_genus_WGS_PT_BIOM <- make_biom(t(rep200Data_Matched2ImmunePT_Fungi_genus_WGS_PT))
write_biom(rep200Data_Matched2ImmunePT_Fungi_genus_WGS_PT_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_rep200_counts_fungi_TCGA_AllSeqCenters_WGS_Primary_Tumor.biom")
mycobiomeImmuneCellsMergedPT_WGS_BIOM <- make_biom(t(mycobiomeImmuneCellsMergedPT_WGS))
write_biom(mycobiomeImmuneCellsMergedPT_WGS_BIOM, biom_file = "MMvec_analysis/genus_intersected_with_WIS/immune_cibersort_rel_abund_TCGA_AllSeqCenters_WGS_Primary_Tumor.biom")

#----------------------------------------------------------#
# Load all of Thorsson et al. 2018 Supplementary Table 1
#----------------------------------------------------------#

# Load Table S1 from Thorsson et al. 2018 Immunity *excluding* the CIBERSORT data
# Only TCGA case IDs were provided for the mapping, so those are used here to
# map to TCGA mycobiome data
thorssonTableS1 <- read.csv("Supporting_data/thorsson_2018_immunity_full_supplementary_table_S1.csv",
                          row.names = 1, stringsAsFactors = FALSE)
thorssonTableS1Formatted <- thorssonTableS1 %>% rownames_to_column("tcga_case_id")
## Merge by tcga_case_id in metaQiitaCombined_Nonzero_SpeciesShared
dim(metaQiitaCombined_Nonzero_SpeciesShared) # 12750    41

# Check for non-unique overlap
sum(metaQiitaCombined_Nonzero_SpeciesShared$tcga_case_id %in% rownames(thorssonTableS1)) # 12609 (NB: is less than 14546 --> missing samples)
sum(rownames(thorssonTableS1) %in% metaQiitaCombined_Nonzero_SpeciesShared$tcga_case_id) # 7576 --> repeats exist

# Create mycobiome-matched cibersort immune cell table via merging
metaQiitaCombined_Nonzero_SpeciesShared %>% rownames_to_column("sampleid") %>%
  select(sampleid, tcga_case_id) %>% droplevels() %>%
  left_join(thorssonTableS1Formatted, by = "tcga_case_id") %>%
  drop_na(Immune.Subtype) %>%
  select(-tcga_case_id) %>%
  column_to_rownames("sampleid") -> mycobiomeThorssonTableS1Merged
colnames(mycobiomeThorssonTableS1Merged) <- gsub('([[:punct:]])|\\s+','',colnames(mycobiomeThorssonTableS1Merged))
dim(mycobiomeThorssonTableS1Merged) # 11101    35

# Subset metadata
metaThorsson <- droplevels(metaQiitaCombined_Nonzero_SpeciesShared[rownames(mycobiomeThorssonTableS1Merged),])
metaThorssonPT <- metaThorsson %>% filter(sample_type == "Primary Tumor") %>% droplevels()
metaThorssonPT_WGS <- metaThorssonPT %>% filter(experimental_strategy == "WGS") %>% droplevels()
metaThorssonPT_RNA <- metaThorssonPT %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()
metaThorssonPT_WGS_HMS <- metaThorssonPT_WGS %>% filter(data_submitting_center_label == "Harvard Medical School") %>% droplevels()
metaThorssonPT_RNA_UNC <- metaThorssonPT_RNA %>% filter(data_submitting_center_label == "University of North Carolina") %>% droplevels()

# Subset immune cell data
mycobiomeThorssonTableS1MergedPT <- mycobiomeThorssonTableS1Merged[rownames(metaThorssonPT),]
mycobiomeThorssonTableS1MergedPT_WGS <- mycobiomeThorssonTableS1Merged[rownames(metaThorssonPT_WGS),]
mycobiomeThorssonTableS1MergedPT_RNA <- mycobiomeThorssonTableS1Merged[rownames(metaThorssonPT_RNA),]
mycobiomeThorssonTableS1MergedPT_WGS_HMS <- mycobiomeThorssonTableS1Merged[rownames(metaThorssonPT_WGS_HMS),]
mycobiomeThorssonTableS1MergedPT_RNA_UNC <- mycobiomeThorssonTableS1Merged[rownames(metaThorssonPT_RNA_UNC),]

mycobiomeThorssonTableS1MergedPT_WGS %>%
  rownames_to_column("sample_id") %>% 
  write.csv("MMvec_analysis/thorsson_et_al_2018_tableS1_other_metadata_02Nov21.csv", row.names = FALSE)

#----------------------------------------------------------#
# ANCOM-BC tumor vs NAT
# - Subset data to individual seq center and data type
#----------------------------------------------------------#

# load("Interim_data/snmDataFungi_DecontamV2_13Oct21.RData")
load("Interim_data/decontamResultsV2_13Oct21.RData")

metaQiitaCombined_Nonzero_DecontamV2_TvsNAT <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal")) %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_TvsNAT$sample_type <- relevel(metaQiitaCombined_Nonzero_DecontamV2_TvsNAT$sample_type, ref = "Solid Tissue Normal")
rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_TvsNAT <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_TvsNAT),]

psFungiTCGADecontam_TvsNAT <- phyloseq(otu_table(rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_TvsNAT, taxa_are_rows = FALSE), 
                             tax_table(as.matrix(rep200TaxSplit_Fungi_Paired_to_Weizmann)), 
                             sample_data(metaQiitaCombined_Nonzero_DecontamV2_TvsNAT))
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

## Source function
source("00-Functions.R") # for runAncomBC_TvsNAT() function
## Call function
runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_HMS)
runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_BCM)
runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_MDA)
runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_WashU)
runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_Broad_WGS)
runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_UNC, ancombcLibCut = 0)
runAncomBC_TvsNAT(psFungiTCGADecontam_TvsNAT_CMS, ancombcLibCut = 0)
# Broad_RNA does not have sufficient samples for comparing T vs NAT

#----------------------------------------------------------#
# ANCOM-BC tumor 1 cancer type vs all others: fungi
# - Subset data to individual seq center and data type
#----------------------------------------------------------#

load("Interim_data/decontamResultsV2_13Oct21.RData")

metaQiitaCombined_Nonzero_DecontamV2_PT <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(sample_type == "Primary Tumor") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_BDN <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(sample_type == "Blood Derived Normal") %>% droplevels()

rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_PT <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_PT),]
rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_BDN <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_BDN),]

# rep200FungiSpeciesShared_Nonzero,
# metaQiitaCombined_Nonzero_SpeciesShared,

# Run using WIS-intersected features. Need to reconvert species names back to OGUs. Ok since only 1 OGU/species
taxTableWISconversion <- rep200TaxSplit_Fungi_Paired_to_Weizmann %>% rownames_to_column("OGUs") %>% column_to_rownames("species")
metaQiitaCombined_Nonzero_SpeciesShared_PT <- metaQiitaCombined_Nonzero_SpeciesShared %>% filter(sample_type == "Primary Tumor") %>% droplevels()
rep200FungiSpeciesShared_Nonzero_PT <- rep200FungiSpeciesShared_Nonzero[rownames(metaQiitaCombined_Nonzero_SpeciesShared_PT),]
colnames(rep200FungiSpeciesShared_Nonzero_PT) <- taxTableWISconversion[colnames(rep200FungiSpeciesShared_Nonzero_PT), "OGUs"]

source("00-Functions.R") # for runAncomBC_1VsAll_Fungi() function
# Run for WIS intersected features
runAncomBC_1VsAll_Fungi(countData=rep200Data_Matched2ImmunePT_Fungi_species,
                        metaData = metaImmunePT,
                        taxTable = rep200TaxSplit_Fungi_Paired_to_Weizmann,
                        qvalCutoff = 0.05, showTopX = 3, decontamResV2 = decontamResultsV2,
                        ancombcLibCut = 100,
                        fileString = "WIS_fungi_ancombc_PT_1vsAll_",
                        taxaPlotLabel = "genus")
# Run for decontaminated data v2 PT
runAncomBC_1VsAll_Fungi(countData=rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_PT,
                        metaData = metaQiitaCombined_Nonzero_DecontamV2_PT,
                        taxTable = rep200TaxSplit_Fungi_Paired_to_Weizmann,
                        qvalCutoff = 0.05, showTopX = 3, decontamResV2 = decontamResultsV2,
                        ancombcLibCut = 1000,
                        fileString = "decontamV2_fungi_ancombc_PT_1vsAll_",
                        taxaPlotLabel = "genus")

# Run for decontaminated data v2 BDN
runAncomBC_1VsAll_Fungi(countData=rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_BDN,
                        metaData = metaQiitaCombined_Nonzero_DecontamV2_BDN,
                        taxTable = rep200TaxSplit_Fungi_Paired_to_Weizmann,
                        qvalCutoff = 0.05, showTopX = 3, decontamResV2 = decontamResultsV2,
                        ancombcLibCut = 1000,
                        fileString = "bdn_decontamV2_fungi_ancombc_PT_1vsAll_",
                        taxaPlotLabel = "genus")

#----------------------------------------------------------#
# ANCOM-BC tumor 1 cancer type vs all others: bacteria
# - Subset data to individual seq center and data type
#----------------------------------------------------------#
# save(rep200TaxSplit_Bacteria_Formatted,
#      sharedPhylumBacteria,
#      sharedClassBacteria,
#      sharedOrderBacteria,
#      sharedFamilyBacteria,
#      sharedGenusBacteria,
#      sharedSpeciesBacteria,
#      file = "Interim_data/shared_bacterial_features_at_each_taxa_level_15Oct21.RData")
load("Interim_data/shared_bacterial_features_at_each_taxa_level_15Oct21.RData")

# Sanity check
rep200Data_Matched2ImmunePT_Bacteria_WIS <- rep200Data_Matched2ImmunePT_Bacteria[,colnames(rep200Data_Matched2ImmunePT_Bacteria) %in%
                                                                                   sharedSpeciesBacteria$intersectedOGUs]
all(rownames(rep200Data_Matched2ImmunePT_Bacteria_WIS)==rownames(metaImmunePT)) # TRUE
rep200Data_Matched2ImmunePT_Bacteria_WIS_Nonzero <- rep200Data_Matched2ImmunePT_Bacteria_WIS[!rowSums(rep200Data_Matched2ImmunePT_Bacteria_WIS)==0,]
metaImmunePT_Bacteria_WIS_Nonzero <- droplevels(metaImmunePT[!rowSums(rep200Data_Matched2ImmunePT_Bacteria_WIS)==0,])

# For WIS intersected features
source("00-Functions.R") # for runAncomBC_1VsAll_Bacteria() function
runAncomBC_1VsAll_Bacteria(countData=rep200Data_Matched2ImmunePT_Bacteria_WIS_Nonzero,
                        metaData = metaImmunePT_Bacteria_WIS_Nonzero,
                        taxTable = rep200TaxSplit_Bacteria_Formatted,
                        ancombcLibCut = 100,
                        qvalCutoff = 0.05, showTopX = 3,
                        fileString = "WIS_bacteria_nonzero_ancombc_PT_1vsAll_",
                        taxaPlotLabel = "genus")

#----------------------------------------------------------#
# ANCOM-BC tumor stage I vs stage IV for fungi
# - Subset data to individual seq center and data type
#----------------------------------------------------------#
# Remove unclear or non-useful stages
metaQiitaCombined_Nonzero_DecontamV2_Path <- droplevels(metaQiitaCombined_Nonzero_DecontamV2[! (metaQiitaCombined_Nonzero_DecontamV2$pathologic_stage_label == "Not available" |
                                                                                                  metaQiitaCombined_Nonzero_DecontamV2$pathologic_stage_label == "I or II NOS" |
                                                                                                  metaQiitaCombined_Nonzero_DecontamV2$pathologic_stage_label == "Stage 0" |
                                                                                                  metaQiitaCombined_Nonzero_DecontamV2$pathologic_stage_label == "Stage IS" |
                                                                                                  metaQiitaCombined_Nonzero_DecontamV2$pathologic_stage_label == "Stage Tis" |
                                                                                                  metaQiitaCombined_Nonzero_DecontamV2$pathologic_stage_label == "Stage X"),])

tumorStageVector <- factor(metaQiitaCombined_Nonzero_DecontamV2_Path$pathologic_stage_label)
levels(tumorStageVector) <- list(StageI = c("Stage I", "Stage IA", "Stage IB", "Stage IC"),
                                 StageII = c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC"),
                                 StageIII = c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),
                                 StageIV = c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC"))
metaQiitaCombined_Nonzero_DecontamV2_Path$pathologic_stage_label_binned <- tumorStageVector
table(metaQiitaCombined_Nonzero_DecontamV2_Path$pathologic_stage_label_binned)

metaQiitaCombined_Nonzero_DecontamV2_Path_PT <- metaQiitaCombined_Nonzero_DecontamV2_Path %>% filter(sample_type == "Primary Tumor") %>% droplevels()
rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_Path_PT <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_Path_PT),]

# Run for decontaminated data v2
source("00-Functions.R") # for runAncomBC_Stage_Fungi() function
# Run Stage I vs Stage IV
runAncomBC_Stage_Fungi(countData=rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_Path_PT,
                       metaDataPath = metaQiitaCombined_Nonzero_DecontamV2_Path_PT,
                       taxTable = rep200TaxSplit_Fungi_Paired_to_Weizmann,
                       stagesCmp = c("Stage I","Stage IV"),
                       stageAllFlag = TRUE,
                       qvalCutoff = 0.05, showTopX = 3, decontamResV2 = decontamResultsV2,
                       ancombcLibCutWGS = 100,
                       fileString = "stage_decontamV2_fungi_ancombc_PT_StageI_vs_StageIV_",
                       taxaPlotLabel = "genus")
# Run Stage I-II ("early stage") vs Stage III-IV ("late stage")
runAncomBC_Stage_Fungi(countData=rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_Path_PT,
                        metaDataPath = metaQiitaCombined_Nonzero_DecontamV2_Path_PT,
                        taxTable = rep200TaxSplit_Fungi_Paired_to_Weizmann,
                        stagesCmp = c("Stage I","Stage II", "Stage III", "Stage IV"),
                        stageAllFlag = TRUE,
                        qvalCutoff = 0.05, showTopX = 3, decontamResV2 = decontamResultsV2,
                        ancombcLibCutWGS = 100,
                        fileString = "stage_decontamV2_fungi_ancombc_PT_Early_vs_Late_",
                        taxaPlotLabel = "genus")

#----------------------------------------------------------#
# ANCOM-BC tumor stage I vs stage IV for bacteria
# - Subset data to individual seq center and data type
#----------------------------------------------------------#

load("Interim_data/shared_bacterial_features_at_each_taxa_level_15Oct21.RData")

# Sanity check
rep200Data_Matched2ImmunePT_Bacteria_WIS <- rep200Data_Matched2ImmunePT_Bacteria[,colnames(rep200Data_Matched2ImmunePT_Bacteria) %in%
                                                                                   sharedSpeciesBacteria$intersectedOGUs]
all(rownames(rep200Data_Matched2ImmunePT_Bacteria_WIS)==rownames(metaImmunePT)) # TRUE
rep200Data_Matched2ImmunePT_Bacteria_WIS_Nonzero <- rep200Data_Matched2ImmunePT_Bacteria_WIS[!rowSums(rep200Data_Matched2ImmunePT_Bacteria_WIS)==0,]
metaImmunePT_Bacteria_WIS_Nonzero <- droplevels(metaImmunePT[!rowSums(rep200Data_Matched2ImmunePT_Bacteria_WIS)==0,])

metaImmunePT_Bacteria_WIS_Nonzero_Path <- droplevels(metaImmunePT_Bacteria_WIS_Nonzero[! (metaImmunePT_Bacteria_WIS_Nonzero$pathologic_stage_label == "Not available" |
                                                                                            metaImmunePT_Bacteria_WIS_Nonzero$pathologic_stage_label == "I or II NOS" |
                                                                                            metaImmunePT_Bacteria_WIS_Nonzero$pathologic_stage_label == "Stage 0" |
                                                                                            metaImmunePT_Bacteria_WIS_Nonzero$pathologic_stage_label == "Stage IS" |
                                                                                            metaImmunePT_Bacteria_WIS_Nonzero$pathologic_stage_label == "Stage Tis" |
                                                                                            metaImmunePT_Bacteria_WIS_Nonzero$pathologic_stage_label == "Stage X"),])

tumorStageVectorBacteria <- factor(metaImmunePT_Bacteria_WIS_Nonzero_Path$pathologic_stage_label)
levels(tumorStageVectorBacteria) <- list(StageI = c("Stage I", "Stage IA", "Stage IB", "Stage IC"),
                                 StageII = c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC"),
                                 StageIII = c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),
                                 StageIV = c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC"))
metaImmunePT_Bacteria_WIS_Nonzero_Path$pathologic_stage_label_binned <- tumorStageVectorBacteria
table(metaImmunePT_Bacteria_WIS_Nonzero_Path$pathologic_stage_label_binned)

metaImmunePT_Bacteria_WIS_Nonzero_Path_PT <- metaImmunePT_Bacteria_WIS_Nonzero_Path %>% filter(sample_type == "Primary Tumor") %>% droplevels()
rep200Data_Matched2ImmunePT_Bacteria_WIS_Nonzero_Path_PT <- rep200Data_Matched2ImmunePT_Bacteria_WIS_Nonzero[rownames(metaImmunePT_Bacteria_WIS_Nonzero_Path_PT),]

# For WIS intersected features
source("00-Functions.R") # for runAncomBC_Stage_Bacteria() function
# Run Stage I vs Stage IV
runAncomBC_Stage_Bacteria(countData=rep200Data_Matched2ImmunePT_Bacteria_WIS_Nonzero_Path_PT,
                           metaData = metaImmunePT_Bacteria_WIS_Nonzero_Path_PT,
                           taxTable = rep200TaxSplit_Bacteria_Formatted,
                           stagesCmp = c("Stage I","Stage IV"),
                           ancombcLibCutWGS = 100,
                           qvalCutoff = 0.05, showTopX = 3,
                           fileString = "stage_WIS_bacteria_nonzero_ancombc_PT_StageI_vs_StageIV_",
                           taxaPlotLabel = "genus")
# Run Stage I-II ("early stage") vs Stage III-IV ("late stage")
runAncomBC_Stage_Bacteria(countData=rep200Data_Matched2ImmunePT_Bacteria_WIS_Nonzero_Path_PT,
                          metaData = metaImmunePT_Bacteria_WIS_Nonzero_Path_PT,
                          taxTable = rep200TaxSplit_Bacteria_Formatted,
                          stagesCmp = c("Stage I","Stage II", "Stage III", "Stage IV"),
                          stageAllFlag = TRUE,
                          ancombcLibCutWGS = 100,
                          qvalCutoff = 0.05, showTopX = 3,
                          fileString = "stage_WIS_bacteria_nonzero_ancombc_PT_StageI_vs_StageIV_",
                          taxaPlotLabel = "genus")

#----------------------------------------------------------#
# PVCA of intersected TCGA and WIS data
#----------------------------------------------------------#
## Subset to WGS only data, as MMvec will only be done on WGS data
# DecontamV2 fungi data
metaQiitaCombined_Nonzero_DecontamV2_PT_WGS <- metaQiitaCombined_Nonzero_DecontamV2_PT %>% filter(experimental_strategy == "WGS") %>% droplevels()
rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_PT_WGS <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_PT[rownames(metaQiitaCombined_Nonzero_DecontamV2_PT_WGS),]
# Shared WIS fungi data
rep200Data_Matched2ImmunePT_WGS_Fungi_species <- rep200Data_Matched2ImmunePT_Fungi_species[rownames(metaImmunePT_WGS),]
# Shared WIS bacterial data
metaImmunePT_WGS_Bacteria_WIS_Nonzero <- metaImmunePT_Bacteria_WIS_Nonzero %>% filter(experimental_strategy == "WGS") %>% droplevels()
rep200Data_Matched2ImmunePT_WGS_Bacteria_WIS_Nonzero <- rep200Data_Matched2ImmunePT_Bacteria_WIS_Nonzero[rownames(metaImmunePT_WGS_Bacteria_WIS_Nonzero),]
rep200Data_Matched2ImmunePT_Bacteria_MMvec <- rep200Data_Matched2ImmunePT_Bacteria[rownames(metaImmunePT_WGS_Bacteria_WIS_Nonzero),]

dim(metaQiitaCombined_Nonzero_DecontamV2_PT_WGS) # 2005 41
dim(metaImmunePT_WGS) # 1916 41
dim(metaImmunePT_WGS_Bacteria_WIS_Nonzero) # 1910 41

save(rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_PT_WGS,
     metaQiitaCombined_Nonzero_DecontamV2_PT_WGS,
     rep200Data_Matched2ImmunePT_WGS_Fungi_species,
     metaImmunePT_WGS,
     rep200Data_Matched2ImmunePT_WGS_Bacteria_WIS_Nonzero,
     metaImmunePT_WGS_Bacteria_WIS_Nonzero,
     rep200Data_Matched2ImmunePT_Bacteria_MMvec,
     file = "Interim_data/data_for_pvca_mmvec_fungi_and_bacteria_WGS_19Oct21.RData")

# Scripts: S15

# metaQiitaCombined_Nonzero_DecontamV2_PT <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(sample_type == "Primary Tumor") %>% droplevels()
# rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_PT <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_PT),]

## Import and run pvca script locally
source("Supporting_scripts/S02-pvca-function.R")
require(lme4) # for linear modelling calculations

metaFiltered <- metaImmunePT_WGS[,c("disease_type",
                              "data_submitting_center_label")]

pvcaThreshold <- 0.8
pvcaFungiDecontamV2_WGS <- PVCA(counts = t(rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_PT_WGS),
                           meta = metaQiitaCombined_Nonzero_DecontamV2_PT_WGS[,c("disease_type",
                                                      "data_submitting_center_label")],
                           threshold = pvcaThreshold,
                           inter = FALSE)
round(pvcaFungiDecontamV2_WGS, 3)
pvcaFungiSharedWIS_WGS <- PVCA(counts = t(rep200Data_Matched2ImmunePT_WGS_Fungi_species),
                            meta = metaImmunePT_WGS[,c("disease_type",
                                                       "data_submitting_center_label")],
                            threshold = pvcaThreshold,
                            inter = FALSE)
round(pvcaFungiSharedWIS_WGS, 3)
pvcaBacteriaSharedWIS_WGS <- PVCA(counts = t(rep200Data_Matched2ImmunePT_WGS_Bacteria_WIS_Nonzero),
                               meta = metaImmunePT_WGS_Bacteria_WIS_Nonzero[,c("disease_type",
                                                          "data_submitting_center_label")],
                               threshold = pvcaThreshold,
                               inter = FALSE)
round(pvcaBacteriaSharedWIS_WGS, 3)





