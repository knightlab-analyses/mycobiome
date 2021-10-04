#-----------------------------------------------------------------------------
# 07-TCGA-RA-barplot-and-compare-tumor-vs-normal-bray-curtis.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Plot relative abundance barplot at Class level of batch corrected data
# - Calculate average relative abundances on TCGA data across groups (tumor vs. normal)
# - Use bray curtis to calculate distances between groups (tumor vs. normal)
# - Plot results in 2D (using ggplot) and 3D (using plotly)
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
require(reshape2)
require(ggpubr)
require(ggsci)
require(ggforce)
require(concaveman)

numCores <- detectCores()
registerDoMC(cores=numCores)

#----------------------------------------------------------#
# Import data: cancers matched to Weizmann
#----------------------------------------------------------#

# save(metaQiitaCombined_Nonzero_8cancer,
#      rep200FungiPhylum_8cancer,
#      rep200FungiClass_8cancer,
#      rep200FungiOrder_8cancer,
#      rep200FungiFamily_8cancer,
#      rep200FungiGenus_8cancer,
#      rep200FungiSpecies_8cancer,
#      file = "Interim_data/raw_data_tcga_8_cancers_16Sep21.RData")

load("Interim_data/raw_data_tcga_8_cancers_16Sep21.RData")

# save(metaQiitaCombined_Nonzero_8cancer,
#      rep200FungiOrder_8cancer_VSNM,
#      rep200FungiFamily_8cancer_VSNM,
#      rep200FungiGenus_8cancer_VSNM,
#      rep200FungiSpecies_8cancer_VSNM,
#      rep200FungiFamily_8cancer_VSNM_CT,
#      rep200FungiGenus_8cancer_VSNM_CT,
#      rep200FungiSpecies_8cancer_VSNM_CT,
#      file = "Interim_data/data_for_ml_8_cancers_13Sep21.RData")

load("Interim_data/data_for_ml_8_cancers_13Sep21.RData")

# save(metaQiitaCombined_Nonzero,
#      rep200FungiPhylum,
#      rep200FungiClass,
#      rep200FungiOrder,
#      rep200FungiFamily,
#      rep200FungiGenus,
#      rep200FungiSpecies,
#      file = "Interim_data/raw_data_tcga_all_cancers_16Sep21.RData")

load("Interim_data/raw_data_tcga_all_cancers_16Sep21.RData")

# save(metaQiitaCombined_Nonzero,
#      rep200FungiOrderVSNM,
#      rep200FungiFamilyVSNM,
#      rep200FungiGenusVSNM,
#      rep200FungiSpeciesVSNM,
#      rep200FungiGenusVSNM_CT,
#      rep200FungiSpeciesVSNM_CT,
#      snmDataOGUFungi,
#      file = "Interim_data/data_for_ml_tcga_13Sep21.RData")

load("Interim_data/data_for_ml_tcga_13Sep21.RData")

# # snmDataOGUFungi,
# # vdge_dataE,
# # rep200Data_WGS_RNA_HiSeq_Fungi_Decontam_Nonzero,
# # metaQiitaCombined_Nonzero,
load("Interim_data/snmDataFungi_13Sep21.RData") # To load the metadata and raw data objects
# load("Interim_data/data_for_ml_tcga_by_seq_center_and_experimental_strategy_16Sep21.RData")

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
# Define phyloseq objects for batch corrected TCGA data
# based on scaled counts
#----------------------------------------------------------#

snmDataOGUFungi_ScaledCounts <- floor(((2^snmDataOGUFungi)/rowSums(2^snmDataOGUFungi))*10^6)

psFungiVSNM <- phyloseq(otu_table(snmDataOGUFungi_ScaledCounts, taxa_are_rows = FALSE), 
                        tax_table(as.matrix(rep200TaxSplit_Fungi)), sample_data(metaQiitaCombined_Nonzero))

# Aggregate counts
psFungiVSNM_phylum = aggregate_taxa(psFungiVSNM, "Phylum")
psFungiVSNM_class = aggregate_taxa(psFungiVSNM, "Class")
psFungiVSNM_order = aggregate_taxa(psFungiVSNM, "Order")
psFungiVSNM_family = aggregate_taxa(psFungiVSNM, "Family")
psFungiVSNM_genus = aggregate_taxa(psFungiVSNM, "Genus")
psFungiVSNM_species = aggregate_taxa(psFungiVSNM, "Species")

# Define for WGS
psFungiVSNM_WGS <- subset_samples(psFungiVSNM, experimental_strategy == "WGS")
psFungiVSNM_RNA <- subset_samples(psFungiVSNM, experimental_strategy == "RNA-Seq")
psFungiVSNM_PT <- subset_samples(psFungiVSNM, sample_type == "Primary Tumor")
psFungiVSNM_PT_WGS <- subset_samples(psFungiVSNM_PT, experimental_strategy == "WGS")
psFungiVSNM_PT_RNA <- subset_samples(psFungiVSNM_PT, experimental_strategy == "RNA-Seq")

#----------------------------------------------------------#
# Calculate batch corrected relative abundance barplots at Class level
#----------------------------------------------------------#

## VSNM corrected all samples
ps <- tax_glom(psFungiVSNM, "Class")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "investigation")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
ps2 %>% psmelt() %>% ggplot(aes(Sample,Abundance, fill = Class, group = Class)) + geom_bar(stat = "identity") +
  theme_bw() + rotate_x_text(90) + xlab("TCGA Investigation") + ylab("Relative Abundance") +
  ggtitle("Voom-SNM Corrected Fungi Relative Abundance Bar Plot at Class Level for all TCGA samples") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggsave(filename = "Figures/Figure_4/rel_abund_barplot_fungi_class_TCGA_all_samples_VSNM.jpeg",
         dpi = "retina", units = "in", width = 12, height = 6.5)

## VSNM corrected PT only
ps <- tax_glom(psFungiVSNM_PT, "Class")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "investigation")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
ps2 %>% psmelt() %>% ggplot(aes(Sample,Abundance, fill = Class, group = Class)) + geom_bar(stat = "identity") +
  theme_bw() + rotate_x_text(90) + xlab("TCGA Investigation") + ylab("Relative Abundance") +
  ggtitle("Batch Corrected Fungi Relative Abundance Bar Plot at Class Level for TCGA primary tumor samples") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggsave(filename = "Figures/Figure_4/rel_abund_barplot_fungi_class_TCGA_primary_tumor_VSNM.jpeg",
         dpi = "retina", units = "in", width = 12, height = 6.5)
ps2 %>% psmelt() %>% select(OTU, Sample, Abundance, Domain, Phylum, Class) %>%
  write.csv("Figures/Figure_4/rel_abund_barplot_fungi_class_TCGA_primary_tumor_VSNM.csv")

## VSNM corrected PT WGS only
ps <- tax_glom(psFungiVSNM_PT_WGS, "Class")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "investigation")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
ps2 %>% psmelt() %>% ggplot(aes(Sample,Abundance, fill = Class, group = Class)) + geom_bar(stat = "identity") +
  theme_bw() + rotate_x_text(90) + xlab("TCGA Investigation") + ylab("Relative Abundance") +
  ggtitle("Batch Corrected Fungi Relative Abundance Bar Plot at Class Level for TCGA primary tumor WGS samples") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggsave(filename = "Figures/Figure_4/rel_abund_barplot_fungi_class_TCGA_primary_tumor_wgs_VSNM.jpeg",
         dpi = "retina", units = "in", width = 12, height = 6.5)

## VSNM corrected PT RNA only
ps <- tax_glom(psFungiVSNM_PT_RNA, "Class")
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
ps1 <- merge_samples(ps0, "investigation")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
ps2 %>% psmelt() %>% ggplot(aes(Sample,Abundance, fill = Class, group = Class)) + geom_bar(stat = "identity") +
  theme_bw() + rotate_x_text(90) + xlab("TCGA Investigation") + ylab("Relative Abundance") +
  ggtitle("Batch Corrected Fungi Relative Abundance Bar Plot at Class Level for TCGA primary tumor RNA-Seq samples") + theme(plot.title = element_text(hjust = 0.5)) +
  ggsave(filename = "Figures/Figure_4/rel_abund_barplot_fungi_class_TCGA_primary_tumor_rna_batch_corrected.jpeg",
         dpi = "retina", units = "in", width = 12, height = 6.5)

#----------------------------------------------------------#
# Call Bray Curtis function
#----------------------------------------------------------#

all(rownames(rep200FungiSpecies_8cancer_VSNM) == rownames(metaQiitaCombined_Nonzero_8cancer)) # TRUE

source("00-Functions.R") # for compareTvsNATbrayCurtis() function

bcPlot_SpeciesVSNM <- compareTvsNATbrayCurtis(rep200FungiSpeciesVSNM, metaQiitaCombined_Nonzero, axes=c(1,2),
                                grepFlag = TRUE, 
                                greplString = "Lung|Breast|Colon|Rectum|Ovarian", 
                                scalar = 10^4)
bcPlot_SpeciesVSNM$brayPlot2D +
  ggsave("Figures/Figure_4/figure_4_I__tcga_2D_bray_curtis_species_vsnm_scaled.jpeg",
         dpi = "retina", units = "in", width = 8, height = 6)
bcPlot_SpeciesVSNM$brayData %>%
  write.csv("Figures_data/Figure_4/figure_4_I__tcga_2D_bray_curtis_species_vsnm_scaled.csv")


