#-----------------------------------------------------------------------------
# 02B-Calculate-fungi-vs-bacteria-read-distributions.R
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

# load("Interim_data/snmDataFungi_13Sep21.RData") # To load the metaQiitaWGS_RNA_AllSeqPlatforms_Joined object
load("Interim_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_29Sep21.RData")
idxstats <- read.csv("Input_data/cgc_idxstats_mycobiome_all_total_reads_gdp_29Sep21.csv", stringsAsFactors = FALSE, row.names = 1)

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
sum(rownames(rep200Data_WGS) %in% rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined)) # 4736
sum(rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined) %in% rownames(rep200Data_WGS)) # 4736

metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WGS <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined %>% filter(experimental_strategy == "WGS") %>% droplevels()
rep200Data_WGS_Matched <- rep200Data_WGS[rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WGS),]
dim(rep200Data_WGS_Matched) # 4736 11585

#-----------------------Import RNA rep200 data-----------------------#

rep200Data_RNA_BIOM <- read_biom(biom_file = "Input_data/Qiita_results/rep200_TCGA_RNA_OGU_Qiita_analysis_47017_updated_13Sep21.biom")
rep200Data_RNA <- t(as(biom_data(rep200Data_RNA_BIOM), "matrix"))
rownames(rep200Data_RNA) <- gsub("^11[0-9]+\\.","",rownames(rep200Data_RNA)) # Qiita IDs get appended to name; this removes them
dim(rep200Data_RNA)

# Check rowname overlap and subset metadata
sum(rownames(rep200Data_RNA) %in% rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined)) # 10776
sum(rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined) %in% rownames(rep200Data_RNA)) # 10776

metaQiitaWGS_RNA_AllSeqPlatforms_Joined_RNA <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()
rep200Data_RNA_Matched <- rep200Data_RNA[rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_RNA),]
dim(rep200Data_RNA_Matched) # 10776 11735
#-----------------------Combine WGS and RNA rep200 data-----------------------#

sum(colnames(rep200Data_WGS_Matched) %in% colnames(rep200Data_RNA_Matched)) # 11526

rep200Data_WGS_MatchedDf <- rep200Data_WGS_Matched %>% data.frame() %>% rownames_to_column("sampleID")
rep200Data_RNA_MatchedDf <- rep200Data_RNA_Matched %>% data.frame() %>% rownames_to_column("sampleID")

rep200Data_WGS_RNA_Matched <- plyr::rbind.fill(rep200Data_WGS_MatchedDf, rep200Data_RNA_MatchedDf) %>% column_to_rownames("sampleID")
rep200Data_WGS_RNA_Matched[is.na(rep200Data_WGS_RNA_Matched)] <- 0 # rbind.fill places NAs for missing entries; replace them with 0
dim(rep200Data_WGS_RNA_Matched) # 15512 11794

## Save "rep200Data_WGS_RNA_Matched" as phyloseq table for later use
psmetaQiitaWGS_RNA_AllSeqPlatforms_Joined <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined
psmetaQiitaWGS_RNA_AllSeqPlatforms_Joined$library_size <- rowSums(rep200Data_WGS_RNA_Matched)
psmetaQiitaWGS_RNA_AllSeqPlatforms_Joined$library_size_log10 <- log10(psmetaQiitaWGS_RNA_AllSeqPlatforms_Joined$library_size)

psRep200All <- phyloseq(otu_table(rep200Data_WGS_RNA_Matched, taxa_are_rows = FALSE), 
                                  tax_table(as.matrix(rep200TaxSplit)), sample_data(psmetaQiitaWGS_RNA_AllSeqPlatforms_Joined))
# save(psRep200All, fungiOGUs, bacteriaOGUs,
#      file = "Interim_data/phyloseq_tcga_rep200_all_OGUs_16Sep21.RData")

## Find bacterial and fungal OGUs
bacteria_fungi_OGUs <- c(rownames(rep200TaxSplit_Fungi), rownames(rep200TaxSplit_Bacteria))

rep200Data_WGS_RNA_Matched_Filt <- rep200Data_WGS_RNA_Matched[,colnames(rep200Data_WGS_RNA_Matched) %in% bacteria_fungi_OGUs]
dim(rep200Data_WGS_RNA_Matched_Filt) # 15512 11390

rep200Data_WGS_RNA_Matched_Bacteria <- rep200Data_WGS_RNA_Matched[,colnames(rep200Data_WGS_RNA_Matched) %in% rownames(rep200TaxSplit_Bacteria)]
dim(rep200Data_WGS_RNA_Matched_Bacteria) # 15512 11071

rep200Data_WGS_RNA_Matched_Fungi <- rep200Data_WGS_RNA_Matched[,colnames(rep200Data_WGS_RNA_Matched) %in% rownames(rep200TaxSplit_Fungi)]
dim(rep200Data_WGS_RNA_Matched_Fungi) # 15512 319

#------------------------------------------------------------------#
# WGS vs. RNA-Seq norm read counts
#------------------------------------------------------------------#

metaQiitaWGS_RNA_AllSeqPlatforms_Joined_DNAonly <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined %>% filter(experimental_strategy == "WGS") %>% droplevels()
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_RNAonly <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()

wgsData <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_DNAonly),]
rnaData <- rep200Data_WGS_RNA_Matched_Fungi[rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_RNAonly),]

wgsSampleCounts <- log10(unname(rowSums(wgsData)))
rnaSampleCounts <- log10(unname(rowSums(rnaData)))
summary(wgsSampleCounts)
summary(rnaSampleCounts)

combinedSampleCountsJoined <- data.frame(sample_counts = c(wgsSampleCounts,rnaSampleCounts),
                                   data_type = c(rep("WGS",length(wgsSampleCounts)), rep("RNA-Seq",length(rnaSampleCounts))))
require(EnvStats)
combinedSampleCountsJoined %>%
  ggviolin(x = "data_type",
           y = "sample_counts",
           fill = "data_type",
           palette = "nejm",
           legend = "none",
           draw_quantiles = c(0.25,0.50,0.75),
           xlab = "Experimental strategy",
           ylab = "log10(sample fungi read counts)",
           add = "mean",
           add.params = list(color="white",size=1)) +
  stat_compare_means(label.x.npc = 0.1, label.y = 7) +
  stat_n_text(y.pos = -1) + ylim(c(-1.1,7))
ggsave(filename = "Figures/Supplementary_Figures/tcga_read_count_wgs_vs_rna_allSeqPlatforms_04Oct21.svg",
         dpi = "retina",
         height = 5,
         width = 3,
         units = "in")

#----------------------------------------------------------#
# Construct phyloseq object and summarize counts to domain level
#----------------------------------------------------------#
psRep200BacteriaFungi <- phyloseq(otu_table(rep200Data_WGS_RNA_Matched_Filt, taxa_are_rows = FALSE), 
                            tax_table(as.matrix(rep200TaxSplit)), sample_data(metaQiitaWGS_RNA_AllSeqPlatforms_Joined))
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
ra_WGS_RNA <- avgRAbarplot(psRep200BacteriaFungi_PT_domain, "WGS & RNA", "avg_rel_abun_TCGA_primary_tumor_wgs_rna_allSeqPlatforms")
ra_WGS <- avgRAbarplot(psRep200BacteriaFungi_PT_WGS_domain, "WGS only", "avg_rel_abun_TCGA_primary_tumor_wgs_allSeqPlatforms")
ra_RNA <- avgRAbarplot(psRep200BacteriaFungi_PT_RNA_domain, "RNA only", "avg_rel_abun_TCGA_primary_tumor_rna_allSeqPlatforms")

#----------------------------------------------------------#
# Aggregate across groups and normalize by genome size
#----------------------------------------------------------#

genomeSizes <- read.csv("Supporting_data/rep200_genome_lengths.csv", stringsAsFactors = FALSE, row.names = 1)
genomeSizesOrdered <- genomeSizes[colnames(rep200Data_WGS_RNA_Matched_Filt),]

psRep200DataBacteriaFungiNorm <- sweep(rep200Data_WGS_RNA_Matched_Filt, 2, genomeSizesOrdered, FUN = '/')

psRep200BacteriaFungiNorm <- phyloseq(otu_table(psRep200DataBacteriaFungiNorm, taxa_are_rows = FALSE), 
                                      tax_table(as.matrix(rep200TaxSplit)), sample_data(metaQiitaWGS_RNA_AllSeqPlatforms_Joined))
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
raNorm_WGS_RNA <- avgRAbarplot(psRep200BacteriaFungiNorm_PT_domain, "WGS & RNA", filename = "normalized_by_genome_size__avg_rel_abun_TCGA_primary_tumor_wgs_rna_allSeqPlatforms",
                               yAxisLab="Normalized Mean Relative Abundance", title="Normalized (by genome size) average relative abundances: TCGA bacteria vs. fungi")
raNorm_WGS <- avgRAbarplot(psRep200BacteriaFungiNorm_PT_WGS_domain, "WGS only", filename = "normalized_by_genome_size__avg_rel_abun_TCGA_primary_tumor_wgs_allSeqPlatforms",
                           yAxisLab="Normalized Mean Relative Abundance", title="Normalized (by genome size) average relative abundances: TCGA bacteria vs. fungi")
raNorm_RNA <- avgRAbarplot(psRep200BacteriaFungiNorm_PT_RNA_domain, "RNA only", filename = "normalized_by_genome_size_avg_rel_abun_TCGA_primary_tumor_rna_allSeqPlatforms",
                           yAxisLab="Normalized Mean Relative Abundance", title="Normalized (by genome size) average relative abundances: TCGA bacteria vs. fungi")

#----------------------------------------------------------#
# Calculate read %s
#----------------------------------------------------------#

sum(metaQiitaWGS_RNA_AllSeqPlatforms_Joined$cgc_filename %in% rownames(idxstats)) # 15512
sum(rownames(idxstats) %in% metaQiitaWGS_RNA_AllSeqPlatforms_Joined$cgc_filename) # 15512

metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_total_reads <- idxstats[metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$cgc_filename, "total"]
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_mapped_reads <- idxstats[metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$cgc_filename, "mapped"]
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_unmapped_reads <- idxstats[metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$cgc_filename, "unmapped"]
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_ratio_unmapped <- idxstats[metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$cgc_filename, "ratio_unmapped"]

save(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts, 
     rep200Data_WGS_RNA_Matched,
     rep200Data_WGS_RNA_Matched_Bacteria,
     rep200Data_WGS_RNA_Matched_Fungi,
     file = "Interim_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_and_Data_29Sep21.RData")

# Sanity check
all(rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts) == rownames(rep200Data_WGS_RNA_Matched_Bacteria)) # TRUE

#----------------------Raw read %----------------------#
readSumRep200 <- rowSums(rep200Data_WGS_RNA_Matched)
readSumBacteria <- rowSums(rep200Data_WGS_RNA_Matched_Bacteria)
readSumFungi <- rowSums(rep200Data_WGS_RNA_Matched_Fungi)
# Sanity check
all(names(readSumRep200) == rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts)) # TRUE
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_rep200 <- unname(readSumRep200)
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_bacteria <- unname(readSumBacteria)
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_fungi <- unname(readSumFungi)

metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_rep200_total <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_rep200/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_total_reads
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_rep200_unmapped <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_rep200/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_unmapped_reads
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_bacteria_total <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_bacteria/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_total_reads
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_bacteria_unmapped <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_bacteria/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_unmapped
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_fungi_total <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_fungi/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_total_reads
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_fungi_unmapped <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_fungi/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$bam_unmapped

metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$percent_rep200_total <- 100*metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_rep200_total
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$percent_rep200_unmapped <- 100*metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_rep200_unmapped
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$percent_bacteria_total <- 100*metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_bacteria_total
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$percent_bacteria_unmapped <- 100*metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_bacteria_unmapped
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$percent_fungi_total <- 100*metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_fungi_total
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$percent_fungi_unmapped <- 100*metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_fungi_unmapped

# 447 samples have 0 fungal reads although all have >0 bacterial reads
# Running the following line of code will create a version of the metadata
# only with nonzero counts. However, since these 0 counts change the distribution
# of % fungal reads (including their cancer type median rankings), they will be maintained for now
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts_Nonzero <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>%
  filter(readSumFungi > 0) %>% droplevels()
  
cols2Keep <- c("investigation","sample_type","data_submitting_center_label","experimental_strategy",
               "percent_rep200_total", "percent_rep200_unmapped",
               "percent_bacteria_total", "percent_bacteria_unmapped",
               "percent_fungi_total", "percent_fungi_unmapped")
bamCountWithMicrobes <- droplevels(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts[,cols2Keep])
bamCountWithMicrobes.melted <- bamCountWithMicrobes %>%
  rownames_to_column("sampleid") %>%
  mutate(investigation_short = gsub("TCGA-","",investigation)) %>%
  reshape2::melt(id.vars = c("sampleid","investigation_short","investigation","sample_type","data_submitting_center_label","experimental_strategy"))
bamCountWithMicrobes.melted$variable <- factor(bamCountWithMicrobes.melted$variable,
                                               levels = c("percent_bacteria_unmapped", "percent_fungi_unmapped", "percent_rep200_unmapped",
                                                          "percent_bacteria_total", "percent_fungi_total", "percent_rep200_total"))

# save(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts,
#      bamCountWithMicrobes,
#      bamCountWithMicrobes.melted,
#      file = "Interim_data/data_for_read_percentage_plots_AllSeqPlatforms_29Sep21.RData")
#-----------------------Plot PT fungi %s-----------------------#
filePath <- "Figures/Figure_1/"
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_fungi_unmapped","percent_fungi_total")) %>%
  # filter(variable %in% c("percent_fungi_total")) %>%
  filter(sample_type == "Primary Tumor") %>%
  ggplot(aes(reorder(investigation_short, value, FUN=median),value, fill=variable)) +
  geom_boxplot(position = "dodge") + theme_pubr() +
  ylab("Percentage of reads (%)") + xlab("TCGA Cancer Type") +
  ggtitle("Percentage of reads in TCGA primary tumor samples mapped to fungal genomes in rep200") +
  rotate_x_text(90) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
  scale_fill_aaas(labels = c("percent_fungi_unmapped" = "Percentage of unmapped reads\nclassified as fungal (%)",
                             "percent_fungi_total" = "Percentage of total reads\nclassified as fungal (%)")) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5)) +
  # geom_label(label = "Kruskal-Wallis, p < 2.2e-16", x = 0.1) +
  # stat_compare_means(label.x.npc = 0.1, label.y.npc = 0.8) +
  stat_summary(geom = "text", angle = 0,
               fun.data = function(x){c(y = -8, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0))
ggsave(filename = paste0(filePath,"percentage_fungi_unmapped_and_total_primary_tumors_AllSeqPlatforms_29Sep21.pdf"),
         dpi = "retina", units = "in", width = 12, height = 4)

require(rstatix)
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_fungi_unmapped","percent_fungi_total")) %>%
  filter(sample_type == "Primary Tumor") %>%
  group_by(variable) %>%
  anova_test(value ~ investigation_short) %>% data.frame()

#-----------------------Plot PT bacteria %s-----------------------#
filePath <- "Figures/Supplementary_Figures/"
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_bacteria_unmapped","percent_bacteria_total")) %>%
  # filter(variable %in% c("percent_bacteria_unmapped")) %>%
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
  # stat_compare_means(label.x.npc = 0.1, label.y.npc = 0.8) +
  stat_summary(geom = "text", angle = 0,
               fun.data = function(x){c(y = -6, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0))
ggsave(filename = paste0(filePath,"percentage_bacteria_unmapped_and_total_primary_tumors.pdf"),
         dpi = "retina", units = "in", width = 12, height = 4)

require(rstatix)
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_bacteria_unmapped","percent_bacteria_total")) %>%
  filter(sample_type == "Primary Tumor") %>%
  group_by(variable) %>%
  anova_test(value ~ investigation_short) %>% data.frame()

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
  stat_summary(geom = "text", angle = 0,
               fun.data = function(x){c(y = -6, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0))
ggsave(filename = paste0(filePath,"percentage_rep200_unmapped_and_total_primary_tumors.pdf"),
         dpi = "retina", units = "in", width = 12, height = 4)

require(rstatix)
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_rep200_unmapped","percent_rep200_total")) %>%
  filter(sample_type == "Primary Tumor") %>%
  group_by(variable) %>%
  anova_test(value ~ investigation_short) %>% data.frame()

#-----------------------Plot all fungi %s-----------------------#
filePath <- "Figures/Figure_1/"
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
  stat_summary(geom = "text", angle = 0,
               fun.data = function(x){c(y = -8, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0))
ggsave(filename = paste0(filePath,"percentage_fungi_unmapped_and_total_all_sample_types.pdf"),
         dpi = "retina", units = "in", width = 12, height = 4)

require(rstatix)
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_fungi_unmapped","percent_fungi_total")) %>%
  group_by(variable) %>%
  anova_test(value ~ investigation_short) %>% data.frame()

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
  stat_summary(geom = "text", angle = 0,
               fun.data = function(x){c(y = -6, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0))
ggsave(filename = paste0(filePath,"percentage_bacteria_unmapped_and_total_all_sample_types.pdf"),
         dpi = "retina", units = "in", width = 12, height = 4)

require(rstatix)
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_bacteria_unmapped","percent_bacteria_total")) %>%
  group_by(variable) %>%
  anova_test(value ~ investigation_short) %>% data.frame()

#-----------------------Plot all rep200 %s-----------------------#
filePath <- "Figures/Supplementary_Figures/"
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_rep200_unmapped","percent_rep200_total")) %>%
  # filter(variable %in% c("percent_rep200_unmapped")) %>%
  ggplot(aes(reorder(investigation_short, value, FUN=median),value, fill=variable)) +
  geom_boxplot(position = "dodge") + theme_pubr() +
  ylab("Percentage of reads (%)") + xlab("TCGA Cancer Type") +
  ggtitle("Percentage of reads in TCGA across all sample types mapped to all microbial genomes in rep200") +
  rotate_x_text(90) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) + 
  scale_fill_aaas(labels = c("percent_rep200_unmapped" = "Percentage of unmapped reads\nclassified as microbial (%)",
                             "percent_rep200_total" = "Percentage of total reads\nclassified as microbial (%)")) +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5)) +
  # stat_compare_means() +
  stat_summary(geom = "text", angle = 0,
               fun.data = function(x){c(y = -6, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0))
ggsave(filename = paste0(filePath,"percentage_rep200_unmapped_and_total_all_sample_types.pdf"),
         dpi = "retina", units = "in", width = 12, height = 4)

require(rstatix)
bamCountWithMicrobes.melted %>%
  filter(variable %in% c("percent_rep200_unmapped","percent_rep200_total")) %>%
  group_by(variable) %>%
  anova_test(value ~ investigation_short) %>% data.frame()

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
  stat_summary(geom = "text", angle = 90,
               fun.data = function(x){c(y = -7.5, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0.75))
ggsave(filename = paste0(filePath,"percentage_of_total_fungal_vs_bacterial_all_sample_types.pdf"),
         dpi = "retina", units = "in", width = 12, height = 4)

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
  # stat_compare_means(label = "p.signif", method = "wilcox.test", method.args = list(alternative = "less")) +
  stat_compare_means(label = "p.signif", method = "wilcox.test") +
  stat_summary(geom = "text", angle = 90,
               fun.data = function(x){c(y = -7.5, label = length(x) )}, 
               colour = "blue", size = 3,
               position = position_dodge(width = 0.75))
ggsave(filename = paste0(filePath,"percentage_of_total_fungal_vs_bacterial_primary_tumor.pdf"),
         dpi = "retina", units = "in", width = 12, height = 5.5)

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

bamNormCountWithMicrobes <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts
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
  stat_summary(geom = "text", angle = 90,
               fun.data = function(x){c(y = -18, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0.75))
ggsave(filename = paste0(filePath,"normalized_ratio_fungal_vs_bacterial_reads_to_total_reads_all_sample_types.pdf"),
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
  stat_compare_means(label = "p.signif", method = "wilcox.test") +
  stat_summary(geom = "text", angle = 90,
               fun.data = function(x){c(y = -18, label = length(x) )}, 
               colour = "blue",
               position = position_dodge(width = 0.75))
ggsave(filename = paste0(filePath,"normalized_ratio_fungal_vs_bacterial_reads_to_total_reads_primary_tumor.pdf"),
         dpi = "retina", units = "in", width = 12, height = 5.5)

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
write.csv(bamNormCountWithMicrobes_All, file = paste0(filePathCSV,"data_for_plotting_read_percentages_and_ratios_allSeqPlatforms_30Sep21.csv"))

#----------------------------------------------------------#
# Calculate read statistics for main text
#----------------------------------------------------------#

# Percent bacteria and fungi of total reads
require(Rmisc)
formatC(CI(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_bacteria_total)*100)
formatC(CI(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$ratio_fungi_total)*100)

# Sum of total, bacterial, and fungal reads
formatC(sum(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_rep200))
formatC(sum(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_bacteria))
formatC(sum(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_fungi))

# Ratio of bacterial and fungal to rep200
ratio_bacteria_rep200 <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_bacteria/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_rep200
formatC(CI(ratio_bacteria_rep200)*100)
ratio_fungi_rep200 <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_fungi/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts$sum_rep200
formatC(CI(ratio_fungi_rep200)*100)

#----------------------------------------------------------#
# Calculate fungal vs. bacterial % read correlation
#----------------------------------------------------------#
# Color by sample type -- limit to PT/BDN/NAT
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>%
  filter(sample_type %in% c("Primary Tumor","Blood Derived Normal","Solid Tissue Normal")) %>%
  mutate(sample_type = factor(sample_type, levels = c("Primary Tumor","Blood Derived Normal","Solid Tissue Normal"))) %>%
  ggplot(aes(x = percent_bacteria_total, y = percent_fungi_total, color=sample_type)) +
  geom_point(alpha = 0.2) + geom_smooth(method='lm') +
  scale_x_log10() + scale_y_log10() +
  scale_color_aaas() + labs(x = "Percent bacterial reads of total reads (%)", 
                            y = "Percent fungal reads of total reads (%)",
                            color = "Sample Type") +
  theme_pubr() + theme(aspect.ratio=1, legend.position = "right") + coord_fixed()
ggsave(filename = "Figures/Figure_1/corr_log10_fungal_and_bacterial_read_percentages_sample_type_PT_BDN_NAT.pdf",
         dpi = "retina", units = "in", width = 8, height = 8)

metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>%
  filter(sample_type %in% c("Primary Tumor","Blood Derived Normal","Solid Tissue Normal")) %>%
  group_by(sample_type) %>%
  cor_test(percent_fungi_total, percent_bacteria_total, method = "spearman") %>%
  data.frame()

# Color by sample type -- others
# NOTE: "Additional Metastatic" contains only 1 sample and "Additional - New Primary"
# contains only 11 samples, so both of them are removed
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>%
  filter(!(sample_type %in% c("Primary Tumor","Blood Derived Normal","Solid Tissue Normal", "Additional Metastatic", "Additional - New Primary"))) %>%
  ggplot(aes(x = percent_bacteria_total, y = percent_fungi_total, color=sample_type)) +
  geom_point(alpha = 0.2) + 
  scale_x_log10() + scale_y_log10() +
  geom_smooth(method='lm', alpha = 0.4) +
  scale_color_aaas(position = "right") + labs(x = "Percent bacterial reads of total reads (%)", 
                                              y = "Percent fungal reads of total reads (%)",
                                              color = "Sample Type") +
  theme_pubr() + theme(aspect.ratio=1, legend.position = "right") + coord_fixed()
ggsave(filename = "Figures/Figure_1/corr_log10_fungal_and_bacterial_read_percentages_sample_type_others.pdf",
       dpi = "retina", units = "in", width = 8, height = 8)

metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>%
  filter(!(sample_type %in% c("Primary Tumor","Blood Derived Normal","Solid Tissue Normal", "Additional Metastatic", "Additional - New Primary"))) %>%
  group_by(sample_type) %>%
  cor_test(percent_fungi_total, percent_bacteria_total, method = "spearman") %>%
  data.frame()

# Color by experimental strategy
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>%
  filter(!(sample_type %in% c("Additional Metastatic","Additional - New Primary"))) %>%
  ggplot(aes(x = percent_bacteria_total, y = percent_fungi_total, color=experimental_strategy)) +
  geom_point(alpha = 0.2) + geom_smooth(method='lm') + 
  # stat_cor(method = "spearman", cor.coef.name = "rho", show.legend = FALSE) + 
  scale_x_log10() + scale_y_log10() +
  scale_color_aaas() + labs(x = "Percent bacterial reads of total reads (%)", 
                            y = "Percent fungal reads of total reads (%)",
                            color = "Experimental Strategy") +
  theme_bw() + theme(aspect.ratio=1) + coord_fixed()
ggsave(filename = "Figures/Figure_1/corr_log10_fungal_and_bacterial_read_percentages_experimental_strategy.pdf",
         dpi = "retina", units = "in", width = 7, height = 7)

metaQiitaWGS_RNA_AllSeqPlatforms_Joined_WithBamcounts %>%
  filter(!(sample_type %in% c("Additional Metastatic","Additional - New Primary"))) %>%
  group_by(experimental_strategy) %>%
  cor_test(percent_fungi_total, percent_bacteria_total, method = "spearman") %>%
  data.frame()

#----------------------------------------------------------#
# TCGA metadata for table S5
#----------------------------------------------------------#
dataTableList <- list()
for(ii in seq_along(names(table(metaQiitaWGS_RNA_AllSeqPlatforms_Joined$investigation)))){
  metaData <- metaQiitaWGS_RNA_AllSeqPlatforms_Joined
  metaData$investigation <- factor(metaData$investigation)
  study <- levels(metaData$investigation)[ii]
  subsetMetadata <- metaData %>% filter(investigation == study) %>% droplevels()
  seqCenterNum <- length(table(subsetMetadata$data_submitting_center_label))
  subsetMetadataPT <- subsetMetadata %>% filter(sample_type == "Primary Tumor") %>% droplevels()
  subsetMetadataNAT <- subsetMetadata %>% filter(sample_type == "Solid Tissue Normal") %>% droplevels()
  subsetMetadataBDN <- subsetMetadata %>% filter(sample_type == "Blood Derived Normal") %>% droplevels()
  subsetMetadataMet <- subsetMetadata %>% filter(sample_type == "Metastatic") %>% droplevels()
  subsetMetadataOther <- subsetMetadata %>% filter(!(sample_type %in% c("Primary Tumor",
                                                                        "Solid Tissue Normal",
                                                                        "Blood Derived Normal",
                                                                        "Metastatic"))) %>% droplevels()
  subsetMetadataRNA <- subsetMetadata %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()
  subsetMetadataWGS <- subsetMetadata %>% filter(experimental_strategy == "WGS") %>% droplevels()
  dfRes <- data.frame(Study = study, 
                      Centers = seqCenterNum,
                      PercRNA = round(nrow(subsetMetadataRNA)/(nrow(subsetMetadataRNA)+nrow(subsetMetadataWGS)),4),
                      Total = nrow(subsetMetadata),
                      NAT = nrow(subsetMetadataNAT),
                      PT = nrow(subsetMetadataPT),
                      Met = nrow(subsetMetadataMet),
                      BDN = nrow(subsetMetadataBDN),
                      Other = nrow(subsetMetadataOther))
  
  dataTableList[[ii]] <- dfRes}
dataTableRes <- do.call(rbind, dataTableList)
dataTableRes %>% write.csv(file = "Figures/Supplementary_Figures/tcga_data_table_AllSeqPlatforms_21Nov21.csv",
                           row.names = FALSE)

