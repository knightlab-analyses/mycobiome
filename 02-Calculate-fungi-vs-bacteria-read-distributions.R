#-----------------------------------------------------------------------------
# 02-Calculate-fungi-vs-bacteria-read-distributions.R
# Copyright (c) 2021--, Greg Poore
# Purposes:
# - Calculate relative abundances of fungi vs bacteria with and without genome size correction
# - Calculate fungi vs bacteria read %s among all TCGA cancer types and samples
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

load("Interim_data/snmDataFungi_13Sep21.RData") # To load the metaQiitaCombined_Nonzero object
idxstats <- read.csv("Input_data/cgc_idxstats_mycobiome_all_total_reads_gdp.csv", stringsAsFactors = FALSE, row.names = 1)

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
# Remerge WGS and RNA-Seq data while retaining bacteria and fungal data
#----------------------------------------------------------#

#-----------------------Import WGS rep200 data-----------------------#

## Import metadata and read count data
rep200Data_WGS_BIOM <- read_biom(biom_file = "Input_data/Qiita_results/rep200_TCGA_WGS_OGU.biom")
rep200Data_WGS <- t(as(biom_data(rep200Data_WGS_BIOM), "matrix"))
dim(rep200Data_WGS) # 4736 11585

# Check rowname overlap and subset metadata
sum(rownames(rep200Data_WGS) %in% rownames(metaQiitaCombined_Nonzero)) # 4379
sum(rownames(metaQiitaCombined_Nonzero) %in% rownames(rep200Data_WGS)) # 4379

metaQiitaCombined_Nonzero_WGS <- metaQiitaCombined_Nonzero %>% filter(experimental_strategy == "WGS") %>% droplevels()
rep200Data_WGS_Matched <- rep200Data_WGS[rownames(metaQiitaCombined_Nonzero_WGS),]
dim(rep200Data_WGS_Matched) # 4379 11585

#-----------------------Import RNA rep200 data-----------------------#

rep200Data_RNA_BIOM <- read_biom(biom_file = "Input_data/Qiita_results/rep200_TCGA_RNA_OGU_Qiita_analysis_47017_updated_13Sep21.biom")
rep200Data_RNA <- t(as(biom_data(rep200Data_RNA_BIOM), "matrix"))
rownames(rep200Data_RNA) <- gsub("^11[0-9]+\\.","",rownames(rep200Data_RNA)) # Qiita IDs get appended to name; this removes them
dim(rep200Data_RNA)

# Check rowname overlap and subset metadata
sum(rownames(rep200Data_RNA) %in% rownames(metaQiitaCombined_Nonzero)) # 10167
sum(rownames(metaQiitaCombined_Nonzero) %in% rownames(rep200Data_RNA)) # 10167

metaQiitaCombined_Nonzero_RNA <- metaQiitaCombined_Nonzero %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()
rep200Data_RNA_Matched <- rep200Data_RNA[rownames(metaQiitaCombined_Nonzero_RNA),]
dim(rep200Data_RNA_Matched) # 10167 11735
#-----------------------Combine WGS and RNA rep200 data-----------------------#

sum(colnames(rep200Data_WGS_Matched) %in% colnames(rep200Data_RNA_Matched)) # 11526

rep200Data_WGS_MatchedDf <- rep200Data_WGS_Matched %>% data.frame() %>% rownames_to_column("sampleID")
rep200Data_RNA_MatchedDf <- rep200Data_RNA_Matched %>% data.frame() %>% rownames_to_column("sampleID")

rep200Data_WGS_RNA_Matched <- plyr::rbind.fill(rep200Data_WGS_MatchedDf, rep200Data_RNA_MatchedDf) %>% column_to_rownames("sampleID")
rep200Data_WGS_RNA_Matched[is.na(rep200Data_WGS_RNA_Matched)] <- 0 # rbind.fill places NAs for missing entries; replace them with 0
dim(rep200Data_WGS_RNA_Matched) # 14546 11794

## Save "rep200Data_WGS_RNA_Matched" as phyloseq table for later use
psMetaQiitaCombined_Nonzero <- metaQiitaCombined_Nonzero
psMetaQiitaCombined_Nonzero$library_size <- rowSums(rep200Data_WGS_RNA_Matched)
psMetaQiitaCombined_Nonzero$library_size_log10 <- log10(psMetaQiitaCombined_Nonzero$library_size)

psRep200All <- phyloseq(otu_table(rep200Data_WGS_RNA_Matched, taxa_are_rows = FALSE), 
                                  tax_table(as.matrix(rep200TaxSplit)), sample_data(psMetaQiitaCombined_Nonzero))
save(psRep200All, fungiOGUs, bacteriaOGUs,
     file = "Interim_data/phyloseq_tcga_rep200_all_OGUs_16Sep21.RData")

## Find bacterial and fungal OGUs
bacteria_fungi_OGUs <- c(rownames(rep200TaxSplit_Fungi), rownames(rep200TaxSplit_Bacteria))

rep200Data_WGS_RNA_Matched_Filt <- rep200Data_WGS_RNA_Matched[,colnames(rep200Data_WGS_RNA_Matched) %in% bacteria_fungi_OGUs]
dim(rep200Data_WGS_RNA_Matched_Filt) # 14546 11390

rep200Data_WGS_RNA_Matched_Bacteria <- rep200Data_WGS_RNA_Matched[,colnames(rep200Data_WGS_RNA_Matched) %in% rownames(rep200TaxSplit_Bacteria)]
dim(rep200Data_WGS_RNA_Matched_Bacteria) # 14546 11071

rep200Data_WGS_RNA_Matched_Fungi <- rep200Data_WGS_RNA_Matched[,colnames(rep200Data_WGS_RNA_Matched) %in% rownames(rep200TaxSplit_Fungi)]
dim(rep200Data_WGS_RNA_Matched_Fungi) # 14546 319

#----------------------------------------------------------#
# Construct phyloseq object and summarize counts to domain level
#----------------------------------------------------------#
psRep200BacteriaFungi <- phyloseq(otu_table(rep200Data_WGS_RNA_Matched_Filt, taxa_are_rows = FALSE), 
                            tax_table(as.matrix(rep200TaxSplit)), sample_data(metaQiitaCombined_Nonzero))
# Separate into various subsets
psRep200BacteriaFungi_PT <- subset_samples(psRep200BacteriaFungi, sample_type == "Primary Tumor")
psRep200BacteriaFungi_PT_WGS <- subset_samples(psRep200BacteriaFungi_PT, experimental_strategy == "WGS")
psRep200BacteriaFungi_PT_RNA <- subset_samples(psRep200BacteriaFungi_PT, experimental_strategy == "RNA-Seq")

#-----------------------Aggregate counts to domain level-----------------------#
# Note that k__Bacteria = bacteria and k__Eukaryota = fungi
psRep200BacteriaFungi_domain = aggregate_taxa(psRep200BacteriaFungi, "Domain")
psRep200BacteriaFungi_PT_domain = aggregate_taxa(psRep200BacteriaFungi_PT, "Domain")
psRep200BacteriaFungi_PT_WGS_domain = aggregate_taxa(psRep200BacteriaFungi_PT_WGS, "Domain")
psRep200BacteriaFungi_PT_RNA_domain = aggregate_taxa(psRep200BacteriaFungi_PT_RNA, "Domain")

#----------------------Plot grouped taxa bar plot----------------------#
source("00-Functions.R") # contains avgRAbarplot() function
ra_WGS_RNA <- avgRAbarplot(psRep200BacteriaFungi_PT_domain, "WGS & RNA", "avg_rel_abun_TCGA_primary_tumor_wgs_rna")
ra_WGS <- avgRAbarplot(psRep200BacteriaFungi_PT_WGS_domain, "WGS only", "avg_rel_abun_TCGA_primary_tumor_wgs")
ra_RNA <- avgRAbarplot(psRep200BacteriaFungi_PT_RNA_domain, "RNA only", "avg_rel_abun_TCGA_primary_tumor_rna")

#----------------------------------------------------------#
# Aggregate across groups and normalize by genome size
#----------------------------------------------------------#

genomeSizes <- read.csv("Supporting_data/rep200_genome_lengths.csv", stringsAsFactors = FALSE, row.names = 1)
genomeSizesOrdered <- genomeSizes[colnames(rep200Data_WGS_RNA_Matched_Filt),]

psRep200DataBacteriaFungiNorm <- sweep(rep200Data_WGS_RNA_Matched_Filt, 2, genomeSizesOrdered, FUN = '/')
psRep200DataBacteriaFungiNorm[1:3,1:3] # entry [3,1] should be 4.18806e-07 | entry [3,3] should be 1.425966e-05

psRep200BacteriaFungiNorm <- phyloseq(otu_table(psRep200DataBacteriaFungiNorm, taxa_are_rows = FALSE), 
                                      tax_table(as.matrix(rep200TaxSplit)), sample_data(metaQiitaCombined_Nonzero))
# Separate into various subsets
psRep200BacteriaFungiNorm_PT <- subset_samples(psRep200BacteriaFungiNorm, sample_type == "Primary Tumor")
psRep200BacteriaFungiNorm_PT_WGS <- subset_samples(psRep200BacteriaFungiNorm_PT, experimental_strategy == "WGS")
psRep200BacteriaFungiNorm_PT_RNA <- subset_samples(psRep200BacteriaFungiNorm_PT, experimental_strategy == "RNA-Seq")

#-----------------------Aggregate counts to domain level-----------------------#
# Note that k__Bacteria = bacteria and k__Eukaryota = fungi
psRep200BacteriaFungiNorm_PT_domain = aggregate_taxa(psRep200BacteriaFungiNorm_PT, "Domain")
psRep200BacteriaFungiNorm_PT_WGS_domain = aggregate_taxa(psRep200BacteriaFungiNorm_PT_WGS, "Domain")
psRep200BacteriaFungiNorm_PT_RNA_domain = aggregate_taxa(psRep200BacteriaFungiNorm_PT_RNA, "Domain")

#----------------------Plot grouped taxa bar plot----------------------#
source("00-Functions.R") # contains avgRAbarplot() function
raNorm_WGS_RNA <- avgRAbarplot(psRep200BacteriaFungiNorm_PT_domain, "WGS & RNA", filename = "normalized_by_genome_size__avg_rel_abun_TCGA_primary_tumor_wgs_rna",
                               yAxisLab="Normalized Mean Relative Abundance", title="Normalized (by genome size) average relative abundances: TCGA bacteria vs. fungi")
raNorm_WGS <- avgRAbarplot(psRep200BacteriaFungiNorm_PT_WGS_domain, "WGS only", filename = "normalized_by_genome_size__avg_rel_abun_TCGA_primary_tumor_wgs",
                           yAxisLab="Normalized Mean Relative Abundance", title="Normalized (by genome size) average relative abundances: TCGA bacteria vs. fungi")
raNorm_RNA <- avgRAbarplot(psRep200BacteriaFungiNorm_PT_RNA_domain, "RNA only", filename = "normalized_by_genome_size_avg_rel_abun_TCGA_primary_tumor_rna",
                           yAxisLab="Normalized Mean Relative Abundance", title="Normalized (by genome size) average relative abundances: TCGA bacteria vs. fungi")

#----------------------------------------------------------#
# Calculate read %s
#----------------------------------------------------------#

sum(metaQiitaCombined_Nonzero$filename %in% rownames(idxstats)) # 14546
sum(rownames(idxstats) %in% metaQiitaCombined_Nonzero$filename) # 14546

metaQiitaCombined_Nonzero_WithBamcounts <- metaQiitaCombined_Nonzero
metaQiitaCombined_Nonzero_WithBamcounts$bam_total_reads <- idxstats[metaQiitaCombined_Nonzero_WithBamcounts$filename, "total"]
metaQiitaCombined_Nonzero_WithBamcounts$bam_mapped_reads <- idxstats[metaQiitaCombined_Nonzero_WithBamcounts$filename, "mapped"]
metaQiitaCombined_Nonzero_WithBamcounts$bam_unmapped_reads <- idxstats[metaQiitaCombined_Nonzero_WithBamcounts$filename, "unmapped"]
metaQiitaCombined_Nonzero_WithBamcounts$bam_ratio_unmapped <- idxstats[metaQiitaCombined_Nonzero_WithBamcounts$filename, "ratio_unmapped"]

save(metaQiitaCombined_Nonzero_WithBamcounts, 
     rep200Data_WGS_RNA_Matched_Bacteria,
     rep200Data_WGS_RNA_Matched_Fungi,
     file = "Interim_data/metaQiitaCombined_Nonzero_WithBamcounts_and_Data_13Sep21.RData")

# Sanity check
all(rownames(metaQiitaCombined_Nonzero_WithBamcounts) == rownames(rep200Data_WGS_RNA_Matched_Bacteria)) # TRUE

#----------------------Raw read %----------------------#
readSumRep200 <- rowSums(rep200Data_WGS_RNA_Matched)
readSumBacteria <- rowSums(rep200Data_WGS_RNA_Matched_Bacteria)
readSumFungi <- rowSums(rep200Data_WGS_RNA_Matched_Fungi)
# Sanity check
all(names(readSumRep200) == rownames(metaQiitaCombined_Nonzero_WithBamcounts)) # TRUE
metaQiitaCombined_Nonzero_WithBamcounts$sum_rep200 <- unname(readSumRep200)
metaQiitaCombined_Nonzero_WithBamcounts$sum_bacteria <- unname(readSumBacteria)
metaQiitaCombined_Nonzero_WithBamcounts$sum_fungi <- unname(readSumFungi)

metaQiitaCombined_Nonzero_WithBamcounts$ratio_rep200_total <- metaQiitaCombined_Nonzero_WithBamcounts$sum_rep200/metaQiitaCombined_Nonzero_WithBamcounts$bam_total_reads
metaQiitaCombined_Nonzero_WithBamcounts$ratio_rep200_unmapped <- metaQiitaCombined_Nonzero_WithBamcounts$sum_rep200/metaQiitaCombined_Nonzero_WithBamcounts$bam_unmapped_reads
metaQiitaCombined_Nonzero_WithBamcounts$ratio_bacteria_total <- metaQiitaCombined_Nonzero_WithBamcounts$sum_bacteria/metaQiitaCombined_Nonzero_WithBamcounts$bam_total_reads
metaQiitaCombined_Nonzero_WithBamcounts$ratio_bacteria_unmapped <- metaQiitaCombined_Nonzero_WithBamcounts$sum_bacteria/metaQiitaCombined_Nonzero_WithBamcounts$bam_unmapped
metaQiitaCombined_Nonzero_WithBamcounts$ratio_fungi_total <- metaQiitaCombined_Nonzero_WithBamcounts$sum_fungi/metaQiitaCombined_Nonzero_WithBamcounts$bam_total_reads
metaQiitaCombined_Nonzero_WithBamcounts$ratio_fungi_unmapped <- metaQiitaCombined_Nonzero_WithBamcounts$sum_fungi/metaQiitaCombined_Nonzero_WithBamcounts$bam_unmapped

metaQiitaCombined_Nonzero_WithBamcounts$percent_rep200_total <- 100*metaQiitaCombined_Nonzero_WithBamcounts$ratio_rep200_total
metaQiitaCombined_Nonzero_WithBamcounts$percent_rep200_unmapped <- 100*metaQiitaCombined_Nonzero_WithBamcounts$ratio_rep200_unmapped
metaQiitaCombined_Nonzero_WithBamcounts$percent_bacteria_total <- 100*metaQiitaCombined_Nonzero_WithBamcounts$ratio_bacteria_total
metaQiitaCombined_Nonzero_WithBamcounts$percent_bacteria_unmapped <- 100*metaQiitaCombined_Nonzero_WithBamcounts$ratio_bacteria_unmapped
metaQiitaCombined_Nonzero_WithBamcounts$percent_fungi_total <- 100*metaQiitaCombined_Nonzero_WithBamcounts$ratio_fungi_total
metaQiitaCombined_Nonzero_WithBamcounts$percent_fungi_unmapped <- 100*metaQiitaCombined_Nonzero_WithBamcounts$ratio_fungi_unmapped

cols2Keep <- c("investigation","sample_type","data_submitting_center_label","experimental_strategy",
               "percent_rep200_total", "percent_rep200_unmapped",
               "percent_bacteria_total", "percent_bacteria_unmapped",
               "percent_fungi_total", "percent_fungi_unmapped")
bamCountWithMicrobes <- droplevels(metaQiitaCombined_Nonzero_WithBamcounts[,cols2Keep])
bamCountWithMicrobes.melted <- bamCountWithMicrobes %>%
  rownames_to_column("sampleid") %>%
  mutate(investigation_short = gsub("TCGA-","",investigation)) %>%
  reshape2::melt(id.vars = c("sampleid","investigation_short","investigation","sample_type","data_submitting_center_label","experimental_strategy"))
bamCountWithMicrobes.melted$variable <- factor(bamCountWithMicrobes.melted$variable,
                                               levels = c("percent_bacteria_unmapped", "percent_fungi_unmapped", "percent_rep200_unmapped",
                                                          "percent_bacteria_total", "percent_fungi_total", "percent_rep200_total"))

save(metaQiitaCombined_Nonzero_WithBamcounts,
     bamCountWithMicrobes,
     bamCountWithMicrobes.melted,
     file = "Interim_data/data_for_read_percentage_plots_13Sep21.RData")
#-----------------------Plot PT fungi %s-----------------------#
filePath <- "Figures/Figure_1/"
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_fungi_unmapped","percent_fungi_total")) %>%
  filter(sample_type == "Primary Tumor") %>%
  ggplot(aes(reorder(investigation_short, value, FUN=median),value, fill=variable)) +
  geom_boxplot(position = "dodge") + theme_pubr() +
  ylab("Percentage of reads (%)") + xlab("TCGA Cancer Type") +
  ggtitle("Percentage of reads in TCGA primary tumors mapped to fungal genomes in rep200") +
  rotate_x_text(90) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
  scale_fill_aaas(labels = c("percent_fungi_unmapped" = "Percentage of unmapped reads\nclassified as fungal (%)",
                             "percent_fungi_total" = "Percentage of total reads\nclassified as fungal (%)")) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5)) +
  ggsave(filename = paste0(filePath,"percentage_fungi_unmapped_and_total_primary_tumors.jpeg"),
         dpi = "retina", units = "in", width = 12)

#-----------------------Plot PT bacteria %s-----------------------#
filePath <- "Figures/Supplementary_Figures/"
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_bacteria_unmapped","percent_bacteria_total")) %>%
  filter(sample_type == "Primary Tumor") %>%
  ggplot(aes(reorder(investigation_short, value, FUN=median),value, fill=variable)) +
  geom_boxplot(position = "dodge") + theme_pubr() +
  ylab("Percentage of reads (%)") + xlab("TCGA Cancer Type") +
  ggtitle("Percentage of reads in TCGA primary tumors mapped to bacterial genomes in rep200") +
  rotate_x_text(90) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
  scale_fill_aaas(labels = c("percent_bacteria_unmapped" = "Percentage of unmapped reads\nclassified as bacterial (%)",
                             "percent_bacteria_total" = "Percentage of total reads\nclassified as bacterial (%)")) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5)) +
  ggsave(filename = paste0(filePath,"percentage_bacteria_unmapped_and_total_primary_tumors.jpeg"),
         dpi = "retina", units = "in", width = 12)

#-----------------------Plot PT rep200 %s-----------------------#
filePath <- "Figures/Supplementary_Figures/"
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_rep200_unmapped","percent_rep200_total")) %>%
  filter(sample_type == "Primary Tumor") %>%
  ggplot(aes(reorder(investigation_short, value, FUN=median),value, fill=variable)) +
  geom_boxplot(position = "dodge") + theme_pubr() +
  ylab("Percentage of reads (%)") + xlab("TCGA Cancer Type") +
  ggtitle("Percentage of reads in TCGA primary tumors mapped to all microbial genomes in rep200") +
  rotate_x_text(90) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
  scale_fill_aaas(labels = c("percent_rep200_unmapped" = "Percentage of unmapped reads\nclassified as microbial (%)",
                             "percent_rep200_total" = "Percentage of total reads\nclassified as microbial (%)")) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5)) +
  ggsave(filename = paste0(filePath,"percentage_rep200_unmapped_and_total_primary_tumors.jpeg"),
         dpi = "retina", units = "in", width = 12)

#-----------------------Plot all fungi %s-----------------------#
filePath <- "Figures/Supplementary_Figures/"
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_fungi_unmapped","percent_fungi_total")) %>%
  ggplot(aes(reorder(investigation_short, value, FUN=median),value, fill=variable)) +
  geom_boxplot(position = "dodge") + theme_pubr() +
  ylab("Percentage of reads (%)") + xlab("TCGA Cancer Type") +
  ggtitle("Percentage of reads in TCGA across all sample types mapped to fungal genomes in rep200") +
  rotate_x_text(90) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
  scale_fill_aaas(labels = c("percent_fungi_unmapped" = "Percentage of unmapped reads\nclassified as fungal (%)",
                             "percent_fungi_total" = "Percentage of total reads\nclassified as fungal (%)")) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5)) +
  ggsave(filename = paste0(filePath,"percentage_fungi_unmapped_and_total_all_sample_types.jpeg"),
         dpi = "retina", units = "in", width = 12)

#-----------------------Plot all bacteria %s-----------------------#
filePath <- "Figures/Supplementary_Figures/"
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_bacteria_unmapped","percent_bacteria_total")) %>%
  ggplot(aes(reorder(investigation_short, value, FUN=median),value, fill=variable)) +
  geom_boxplot(position = "dodge") + theme_pubr() +
  ylab("Percentage of reads (%)") + xlab("TCGA Cancer Type") +
  ggtitle("Percentage of reads in TCGA across all sample types mapped to bacterial genomes in rep200") +
  rotate_x_text(90) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
  scale_fill_aaas(labels = c("percent_bacteria_unmapped" = "Percentage of unmapped reads\nclassified as bacterial (%)",
                             "percent_bacteria_total" = "Percentage of total reads\nclassified as bacterial (%)")) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5)) +
  ggsave(filename = paste0(filePath,"percentage_bacteria_unmapped_and_total_all_sample_types.jpeg"),
         dpi = "retina", units = "in", width = 12)

#-----------------------Plot all rep200 %s-----------------------#
filePath <- "Figures/Supplementary_Figures/"
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_rep200_unmapped","percent_rep200_total")) %>%
  ggplot(aes(reorder(investigation_short, value, FUN=median),value, fill=variable)) +
  geom_boxplot(position = "dodge") + theme_pubr() +
  ylab("Percentage of reads (%)") + xlab("TCGA Cancer Type") +
  ggtitle("Percentage of reads in TCGA across all sample types mapped to all microbial genomes in rep200") +
  rotate_x_text(90) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
  scale_fill_aaas(labels = c("percent_rep200_unmapped" = "Percentage of unmapped reads\nclassified as microbial (%)",
                             "percent_rep200_total" = "Percentage of total reads\nclassified as microbial (%)")) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5)) +
  ggsave(filename = paste0(filePath,"percentage_rep200_unmapped_and_total_all_sample_types.jpeg"),
         dpi = "retina", units = "in", width = 12)

#-----------------------Plot all fungal vs bacteria %s-----------------------#
filePath <- "Figures/Figure_1/"
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_fungi_total","percent_bacteria_total")) %>%
  ggplot(aes(reorder(investigation_short, value, FUN=median),value, fill=variable)) +
  geom_boxplot(position = "dodge") + theme_pubr() +
  ylab("Percentage of reads (%)") + xlab("TCGA Cancer Type") +
  ggtitle("Comparing percentage of reads in TCGA across all sample types\nmapped to bacterial vs. fungal genomes in rep200 (no correction for genome size)") +
  rotate_x_text(90) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
  scale_fill_aaas(labels = c("percent_fungi_total" = "Percentage of total reads\nclassified as fungal (%)",
                             "percent_bacteria_total" = "Percentage of total reads\nclassified as bacterial (%)")) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5)) +
  ggsave(filename = paste0(filePath,"percentage_of_total_fungal_vs_bacterial_all_sample_types.jpeg"),
         dpi = "retina", units = "in", width = 12)

#-----------------------Plot PT fungal vs bacteria %s-----------------------#
filePath <- "Figures/Figure_1/"
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_fungi_total","percent_bacteria_total")) %>%
  filter(sample_type == "Primary Tumor") %>%
  ggplot(aes(reorder(investigation_short, value, FUN=median),value, fill=variable)) +
  geom_boxplot(position = "dodge") + theme_pubr() +
  ylab("Percentage of reads (%)") + xlab("TCGA Cancer Type") +
  ggtitle("Comparing percentage of reads in TCGA across primary tumors\nmapped to bacterial vs. fungal genomes in rep200 (no correction for genome size)") +
  rotate_x_text(90) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
  scale_fill_aaas(labels = c("percent_fungi_total" = "Percentage of total reads\nclassified as fungal (%)",
                             "percent_bacteria_total" = "Percentage of total reads\nclassified as bacterial (%)")) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5)) +
  ggsave(filename = paste0(filePath,"percentage_of_total_fungal_vs_bacterial_primary_tumor.jpeg"),
         dpi = "retina", units = "in", width = 12)

#----------------------------------------------------------#
# Calculate norm read %s
#----------------------------------------------------------#
# Look across all genomes in rep200
genomeSizesOrderedNormReads <- genomeSizes[colnames(rep200Data_WGS_RNA_Matched),]

rep200Data_WGS_RNA_Matched_GenomeNorm <- sweep(rep200Data_WGS_RNA_Matched, 2, genomeSizesOrderedNormReads, FUN = '/')
rep200Data_WGS_RNA_Matched_GenomeNorm_Bacteria <- rep200Data_WGS_RNA_Matched_GenomeNorm[,colnames(rep200Data_WGS_RNA_Matched_GenomeNorm) %in% rownames(rep200TaxSplit_Bacteria)]
rep200Data_WGS_RNA_Matched_GenomeNorm_Fungi <- rep200Data_WGS_RNA_Matched_GenomeNorm[,colnames(rep200Data_WGS_RNA_Matched_GenomeNorm) %in% rownames(rep200TaxSplit_Fungi)]

#----------------------Norm read %----------------------#

readSumNormRep200 <- rowSums(rep200Data_WGS_RNA_Matched_GenomeNorm)
readSumNormBacteria <- rowSums(rep200Data_WGS_RNA_Matched_GenomeNorm_Bacteria)
readSumNormFungi <- rowSums(rep200Data_WGS_RNA_Matched_GenomeNorm_Fungi)

bamNormCountWithMicrobes <- metaQiitaCombined_Nonzero_WithBamcounts
bamNormCountWithMicrobes$investigation_short <- gsub("TCGA-","",bamNormCountWithMicrobes$investigation)
bamNormCountWithMicrobes$norm_rep200 <- unname(readSumNormRep200[rownames(bamNormCountWithMicrobes)])
bamNormCountWithMicrobes$norm_bacteria <- unname(readSumNormBacteria[rownames(bamNormCountWithMicrobes)])
bamNormCountWithMicrobes$norm_fungi <- unname(readSumNormFungi[rownames(bamNormCountWithMicrobes)])
bamNormCountWithMicrobes$norm_ratio_rep200_total <- bamNormCountWithMicrobes$norm_rep200/bamNormCountWithMicrobes$bam_total_reads
bamNormCountWithMicrobes$norm_ratio_bacteria_total <- bamNormCountWithMicrobes$norm_bacteria/bamNormCountWithMicrobes$bam_total_reads
bamNormCountWithMicrobes$norm_ratio_fungi_total <- bamNormCountWithMicrobes$norm_fungi/bamNormCountWithMicrobes$bam_total_reads
bamNormCountWithMicrobes$norm_ratio_rep200_unmapped <- bamNormCountWithMicrobes$norm_rep200/bamNormCountWithMicrobes$bam_unmapped_reads
bamNormCountWithMicrobes$norm_ratio_bacteria_unmapped <- bamNormCountWithMicrobes$norm_bacteria/bamNormCountWithMicrobes$bam_unmapped_reads
bamNormCountWithMicrobes$norm_ratio_fungi_unmapped <- bamNormCountWithMicrobes$norm_fungi/bamNormCountWithMicrobes$bam_unmapped_reads

cols2KeepNorm <- c("investigation_short","sample_type","data_submitting_center_label","experimental_strategy",
                   "norm_ratio_rep200_total", "norm_ratio_bacteria_total",
                   "norm_ratio_fungi_total", "norm_ratio_rep200_unmapped",
                   "norm_ratio_bacteria_unmapped", "norm_ratio_fungi_unmapped")
bamNormCountWithMicrobes_Filt <- droplevels(bamNormCountWithMicrobes[,cols2KeepNorm])
bamNormCountWithMicrobes_Filt_melted <- bamNormCountWithMicrobes_Filt %>%
  rownames_to_column("sampleid") %>%
  reshape2::melt(id.vars = c("sampleid","investigation_short","sample_type","data_submitting_center_label","experimental_strategy"))
bamNormCountWithMicrobes_Filt_melted$variable <- factor(bamNormCountWithMicrobes_Filt_melted$variable,
                                                        levels = c("norm_ratio_bacteria_unmapped", "norm_ratio_fungi_unmapped", "norm_ratio_rep200_unmapped",
                                                                   "norm_ratio_bacteria_total", "norm_ratio_fungi_total", "norm_ratio_rep200_total"))
#-----------------------Plot norm all fungal vs bacteria %s-----------------------#
filePath <- "Figures/Figure_1/"
bamNormCountWithMicrobes_Filt_melted %>%
  filter(variable %in% c("norm_ratio_fungi_total","norm_ratio_bacteria_total")) %>%
  ggplot(aes(reorder(investigation_short, value, FUN=median),value, fill=variable)) +
  geom_boxplot(position = "dodge") + theme_pubr() +
  ylab("Ratio of normalized read count to total reads") + xlab("TCGA Cancer Type") +
  ggtitle("Ratio of genome-size normalized read counts to total reads in TCGA across all sample types\nmapped to bacterial vs. fungal genomes in rep200") +
  rotate_x_text(90) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
  scale_fill_aaas(labels = c("norm_ratio_fungi_total" = "Ratio of (fungal reads/genome size):total read count",
                             "norm_ratio_bacteria_total" = "Ratio of (bacterial reads/genome size):total read count")) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5)) +
  ggsave(filename = paste0(filePath,"normalized_ratio_fungal_vs_bacterial_reads_to_total_reads_all_sample_types.jpeg"),
         dpi = "retina", units = "in", width = 12)

#-----------------------Plot norm PT fungal vs bacteria %s-----------------------#
filePath <- "Figures/Figure_1/"
bamNormCountWithMicrobes_Filt_melted %>%
  filter(variable %in% c("norm_ratio_fungi_total","norm_ratio_bacteria_total")) %>%
  filter(sample_type == "Primary Tumor") %>%
  ggplot(aes(reorder(investigation_short, value, FUN=median),value, fill=variable)) +
  geom_boxplot(position = "dodge") + theme_pubr() +
  ylab("Ratio of normalized read count to total reads") + xlab("TCGA Cancer Type") +
  ggtitle("Ratio of genome-size normalized read counts to total reads in TCGA across primary tumors\nmapped to bacterial vs. fungal genomes in rep200") +
  rotate_x_text(90) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
  scale_fill_aaas(labels = c("norm_ratio_fungi_total" = "Ratio of (fungal reads/genome size):total read count",
                             "norm_ratio_bacteria_total" = "Ratio of (bacterial reads/genome size):total read count")) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5)) +
  ggsave(filename = paste0(filePath,"normalized_ratio_fungal_vs_bacterial_reads_to_total_reads_primary_tumor.jpeg"),
         dpi = "retina", units = "in", width = 12)

#-----------------------Save data for plots-----------------------#
cols2KeepAll <- c("investigation_short","sample_type","data_submitting_center_label","experimental_strategy",
                  "norm_ratio_rep200_total", "norm_ratio_bacteria_total",
                  "norm_ratio_fungi_total", "norm_ratio_rep200_unmapped",
                  "norm_ratio_bacteria_unmapped", "norm_ratio_fungi_unmapped",
                  "percent_rep200_total", "percent_rep200_unmapped",
                  "percent_bacteria_total", "percent_bacteria_unmapped",
                  "percent_fungi_total", "percent_fungi_unmapped")
bamNormCountWithMicrobes_All <- droplevels(bamNormCountWithMicrobes[,cols2KeepAll])

filePathCSV <- "Figures_data/Figure_1/"
write.csv(bamNormCountWithMicrobes_All, file = paste0(filePathCSV,"data_for_plotting_read_percentages_and_ratios.csv"))
