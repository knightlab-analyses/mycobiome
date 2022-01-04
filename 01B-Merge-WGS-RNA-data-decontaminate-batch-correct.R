#-----------------------------------------------------------------------------
# 01-Merge-WGS-RNA-data-decontaminate-batch-correct.R
# Copyright (c) 2021--, Greg Poore
# Purposes: 
# - Merge WGS and RNA-Seq datasets
# - Decontaminate using the plate-center method
# - Measure batch effects using PVCA
# - Batch correct using Voom-SNM
#-----------------------------------------------------------------------------

#----------------------------------------------------------#
# Load environments
#----------------------------------------------------------#
require(doMC)
require(plyr)
require(dplyr)
require(tibble)
require(biomformat)
require(rhdf5)
require(ggpubr)
require(ggsci)

numCores <- detectCores()
registerDoMC(cores=numCores)
#----------------------------------------------------------#
# Rep200 fungal species identification
#----------------------------------------------------------#

rep200TaxSplit <- read.csv("Supporting_data/rep200_lineage_map_split.csv", stringsAsFactors = FALSE, row.names = 1)
rep200Kingdoms <- read.csv("Supporting_data/rep200_gOTU_kingdom_mapping.csv", stringsAsFactors = FALSE)
rep200Kingdoms_Fungi <- rep200Kingdoms[which(rep200Kingdoms$kingdom == "fungi"),]
rep200TaxSplit_Fungi <- rep200TaxSplit[rownames(rep200TaxSplit) %in% rep200Kingdoms_Fungi$genomeID,,drop=FALSE]
dim(rep200TaxSplit_Fungi) # 320   7
#----------------------------------------------------------#
# Microbial data import - WGS
#----------------------------------------------------------#

## Import metadata and read count data
metaQiita <- read.csv("Input_data/tcga_wgs_reprocess_qiita_metadata_2Apr21.csv", stringsAsFactors = FALSE, row.names = 1)
rep200Data_WGS_BIOM <- read_biom(biom_file = "Input_data/Qiita_results/rep200_TCGA_WGS_OGU.biom")
rep200Data_WGS <- t(as(biom_data(rep200Data_WGS_BIOM), "matrix"))
dim(rep200Data_WGS) # 4736 11585

# Check rowname overlap and subset metadata
sum(rownames(rep200Data_WGS) %in% rownames(metaQiita)) # 4736
metaQiitaWGS <- droplevels(metaQiita[rownames(rep200Data_WGS),])

## Extract only fungal features
rep200Data_WGS_Fungi <- rep200Data_WGS[,colnames(rep200Data_WGS) %in% rownames(rep200TaxSplit_Fungi)]
dim(rep200Data_WGS_Fungi) # 4736  318

# Extract HiSeq samples
metaQiitaWGS_HiSeq <- metaQiitaWGS %>% filter(platform_tcga == "Illumina HiSeq") %>% droplevels()
dim(metaQiitaWGS_HiSeq) # 4387   39

wgsSampleIDintersectHiSeq <- intersect(rownames(metaQiitaWGS_HiSeq), rownames(rep200Data_WGS))
metaQiitaWGS_HiSeq_Filt <- droplevels(metaQiitaWGS_HiSeq[wgsSampleIDintersectHiSeq,])
rep200Data_WGS_HiSeq <- rep200Data_WGS[wgsSampleIDintersectHiSeq,]

rep200Data_WGS_HiSeq_Fungi <- rep200Data_WGS_HiSeq[,colnames(rep200Data_WGS_HiSeq) %in% rownames(rep200TaxSplit_Fungi)]
# rep200Data_WGS_HiSeq_Fungi <- rep200Data_WGS_HiSeq[,grep("mycota", colnames(rep200Data_WGS_HiSeq))]
dim(rep200Data_WGS_HiSeq_Fungi) # 4387  318

#----------------------------------------------------------#
# Microbial data import - RNA
#----------------------------------------------------------#
## Import metadata and read count data
metaQiitaRNA <- read.csv("Input_data/qiita_metadata_tcga_rna_reprocess_28Apr21.csv", stringsAsFactors = FALSE)
rownames(metaQiitaRNA) <- paste0("13767.",metaQiitaRNA$sample_name)
rep200Data_RNA_BIOM <- read_biom(biom_file = "Input_data/Qiita_results/rep200_TCGA_RNA_OGU_Qiita_analysis_47017_updated_13Sep21.biom")
rep200Data_RNA <- t(as(biom_data(rep200Data_RNA_BIOM), "matrix"))
rownames(rep200Data_RNA) <- gsub("^11[0-9]+\\.","",rownames(rep200Data_RNA)) # Qiita IDs get appended to name; this removes them
dim(rep200Data_RNA) # 10776 11735

# Extract HiSeq samples
metaQiitaRNA_HiSeq <- metaQiitaRNA %>% filter(cgc_platform == "Illumina HiSeq") %>% droplevels()
dim(metaQiitaRNA_HiSeq) # 10701    41

rnaSampleIDintersectHiSeq <- intersect(rownames(metaQiitaRNA_HiSeq),
                                       rownames(rep200Data_RNA))
metaQiitaRNA_HiSeq_Filt <- droplevels(metaQiitaRNA_HiSeq[rnaSampleIDintersectHiSeq,])
rep200Data_RNA_HiSeq <- rep200Data_RNA[rnaSampleIDintersectHiSeq,]

rep200Data_RNA_HiSeq_Fungi <- rep200Data_RNA_HiSeq[,colnames(rep200Data_RNA_HiSeq) %in% rownames(rep200TaxSplit_Fungi)]
dim(rep200Data_RNA_HiSeq_Fungi) # 10701   319

#-----------------------------------------------#
# Combine metadata data - all sequencing platforms
# To be used in "02-Calculate-fungi-vs-bacteria-read-distributions.R"
#-----------------------------------------------#
# metaQiitaWGS, metaQiitaRNA
sum(colnames(metaQiitaWGS) %in% colnames(metaQiitaRNA)) # 28
colnames(metaQiitaRNA)[!(colnames(metaQiitaRNA) %in% colnames(metaQiitaWGS))]
colnames(metaQiitaWGS)[!(colnames(metaQiitaWGS) %in% colnames(metaQiitaRNA))]

colnames(metaQiitaWGS)[which(colnames(metaQiitaWGS) == "id_cgc")] <- "cgc_id"
colnames(metaQiitaWGS)[which(colnames(metaQiitaWGS) == "filename_cgc")] <- "cgc_filename"
colnames(metaQiitaWGS)[which(colnames(metaQiitaWGS) == "sample_id_tcga")] <- "tcga_sample_id"
colnames(metaQiitaWGS)[which(colnames(metaQiitaWGS) == "case_uuid_cgc")] <- "cgc_case_uuid"
colnames(metaQiitaWGS)[which(colnames(metaQiitaWGS) == "platform_tcga")] <- "cgc_platform"
colnames(metaQiitaWGS)[which(colnames(metaQiitaWGS) == "gdc_uuid")] <- "gdc_file_uuid"
colnames(metaQiitaWGS)[which(colnames(metaQiitaWGS) == "sample_uuid_cgc")] <- "cgc_sample_uuid"
colnames(metaQiitaWGS)[which(colnames(metaQiitaWGS) == "aliquot_uuid_cgc")] <- "cgc_aliquot_uuid"
colnames(metaQiitaWGS)[which(colnames(metaQiitaWGS) == "aliquot_id_tcga")] <- "tcga_aliquot_id"
colnames(metaQiitaWGS)[which(colnames(metaQiitaWGS) == "case_id_tcga")] <- "tcga_case_id"
colnames(metaQiitaWGS)[which(colnames(metaQiitaWGS) == "base_name")] <- "cgc_base_name"

# Sanity check
sum(colnames(metaQiitaWGS) %in% colnames(metaQiitaRNA)) # 39
# Intersect overlapping metadata columns and rbind
intersectingMetadataColumns <- intersect(colnames(metaQiitaWGS),
                                         colnames(metaQiitaRNA))
metaQiitaWGS_RNA_AllSeqPlatforms_Joined <- rbind(metaQiitaWGS[,intersectingMetadataColumns],
                                                metaQiitaRNA[,intersectingMetadataColumns])
dim(metaQiitaWGS_RNA_AllSeqPlatforms_Joined) # 15512    39

metaQiitaWGS_RNA_AllSeqPlatforms_Joined %>% count(data_submitting_center_label)
missingSeqCenterAllSeqPlatforms <- rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined)[which(is.na(metaQiitaWGS_RNA_AllSeqPlatforms_Joined$data_submitting_center_label))]
metaQiitaWGS_RNA_AllSeqPlatforms_Joined_Filt <- droplevels(metaQiitaWGS_RNA_AllSeqPlatforms_Joined[!(rownames(metaQiitaWGS_RNA_AllSeqPlatforms_Joined) %in% missingSeqCenterAllSeqPlatforms),])
dim(metaQiitaWGS_RNA_AllSeqPlatforms_Joined_Filt) # 15484    39

# save(metaQiitaWGS_RNA_AllSeqPlatforms_Joined,
#      metaQiitaWGS_RNA_AllSeqPlatforms_Joined_Filt,
#      file = "Interim_data/metaQiitaWGS_RNA_AllSeqPlatforms_Joined_29Sep21.RData")

#-----------------------------------------------#
# Combine WGS and RNA data - HiSeq only
#-----------------------------------------------#

dim(rep200Data_WGS_HiSeq_Fungi) # 4387  318
dim(rep200Data_RNA_HiSeq_Fungi) # 10701   319
dim(metaQiitaWGS_HiSeq_Filt) # 4387   39
dim(metaQiitaRNA_HiSeq_Filt) # 10701    41

sum(colnames(rep200Data_WGS_HiSeq_Fungi) %in% colnames(rep200Data_RNA_HiSeq_Fungi)) # 318
sum(colnames(metaQiitaWGS_HiSeq_Filt) %in% colnames(metaQiitaRNA_HiSeq_Filt)) # 28
colnames(metaQiitaRNA_HiSeq_Filt)[!(colnames(metaQiitaRNA_HiSeq_Filt) %in% colnames(metaQiitaWGS_HiSeq_Filt))]
colnames(metaQiitaWGS_HiSeq_Filt)[!(colnames(metaQiitaWGS_HiSeq_Filt) %in% colnames(metaQiitaRNA_HiSeq_Filt))]

colnames(metaQiitaWGS_HiSeq_Filt)[which(colnames(metaQiitaWGS_HiSeq_Filt) == "id_cgc")] <- "cgc_id"
colnames(metaQiitaWGS_HiSeq_Filt)[which(colnames(metaQiitaWGS_HiSeq_Filt) == "filename_cgc")] <- "cgc_filename"
colnames(metaQiitaWGS_HiSeq_Filt)[which(colnames(metaQiitaWGS_HiSeq_Filt) == "sample_id_tcga")] <- "tcga_sample_id"
colnames(metaQiitaWGS_HiSeq_Filt)[which(colnames(metaQiitaWGS_HiSeq_Filt) == "case_uuid_cgc")] <- "cgc_case_uuid"
colnames(metaQiitaWGS_HiSeq_Filt)[which(colnames(metaQiitaWGS_HiSeq_Filt) == "platform_tcga")] <- "cgc_platform"
colnames(metaQiitaWGS_HiSeq_Filt)[which(colnames(metaQiitaWGS_HiSeq_Filt) == "gdc_uuid")] <- "gdc_file_uuid"
colnames(metaQiitaWGS_HiSeq_Filt)[which(colnames(metaQiitaWGS_HiSeq_Filt) == "sample_uuid_cgc")] <- "cgc_sample_uuid"
colnames(metaQiitaWGS_HiSeq_Filt)[which(colnames(metaQiitaWGS_HiSeq_Filt) == "aliquot_uuid_cgc")] <- "cgc_aliquot_uuid"
colnames(metaQiitaWGS_HiSeq_Filt)[which(colnames(metaQiitaWGS_HiSeq_Filt) == "aliquot_id_tcga")] <- "tcga_aliquot_id"
colnames(metaQiitaWGS_HiSeq_Filt)[which(colnames(metaQiitaWGS_HiSeq_Filt) == "case_id_tcga")] <- "tcga_case_id"
colnames(metaQiitaWGS_HiSeq_Filt)[which(colnames(metaQiitaWGS_HiSeq_Filt) == "base_name")] <- "cgc_base_name"

# Sanity check
sum(colnames(metaQiitaWGS_HiSeq_Filt) %in% colnames(metaQiitaRNA_HiSeq_Filt)) # 39
# Intersect overlapping metadata columns and rbind
intersectingMetadataColumns <- intersect(colnames(metaQiitaWGS_HiSeq_Filt),
                                         colnames(metaQiitaRNA_HiSeq_Filt))
metaQiitaWGS_RNA_HiSeq_Joined <- rbind(metaQiitaWGS_HiSeq_Filt[,intersectingMetadataColumns],
                                     metaQiitaRNA_HiSeq_Filt[,intersectingMetadataColumns])
dim(metaQiitaWGS_RNA_HiSeq_Joined) # 15088    39

# Rbind fungi data
require(gtools)
colnames(rep200Data_RNA_HiSeq_Fungi)[!(colnames(rep200Data_RNA_HiSeq_Fungi) %in% colnames(rep200Data_WGS_HiSeq_Fungi))] # G000277815
rep200Data_WGS_RNA_HiSeq_Fungi <- smartbind(cbind(rep200Data_WGS_HiSeq_Fungi, G000277815=0),
                                        rep200Data_RNA_HiSeq_Fungi)
rownames(rep200Data_WGS_RNA_HiSeq_Fungi) <- c(rownames(rep200Data_WGS_HiSeq_Fungi),
                                              rownames(rep200Data_RNA_HiSeq_Fungi))
dim(rep200Data_WGS_RNA_HiSeq_Fungi) # 15088   319
rep200Data_WGS_HiSeq_Fungi[1:3,1:3]

# Subset WGS and RNA data and save --> to be used for alpha rarefaction
metaQiitaWGS_RNA_HiSeq_Joined_WGS <- metaQiitaWGS_RNA_HiSeq_Joined %>% filter(experimental_strategy == "WGS") %>% droplevels()
metaQiitaWGS_RNA_HiSeq_Joined_RNA <- metaQiitaWGS_RNA_HiSeq_Joined %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()
rep200Data_WGS_RNA_HiSeq_Fungi_WGS <- rep200Data_WGS_RNA_HiSeq_Fungi[rownames(metaQiitaWGS_RNA_HiSeq_Joined_WGS),]
rep200Data_WGS_RNA_HiSeq_Fungi_RNA <- rep200Data_WGS_RNA_HiSeq_Fungi[rownames(metaQiitaWGS_RNA_HiSeq_Joined_RNA),]

# save(metaQiitaWGS_RNA_HiSeq_Joined,
#      rep200Data_WGS_RNA_HiSeq_Fungi,
#      metaQiitaWGS_RNA_HiSeq_Joined_WGS,
#      rep200Data_WGS_RNA_HiSeq_Fungi_WGS,
#      metaQiitaWGS_RNA_HiSeq_Joined_RNA,
#      rep200Data_WGS_RNA_HiSeq_Fungi_RNA,
#      file = "Interim_data/data_for_alpha_rarefaction_29Sep21.RData")

metaQiitaWGS_RNA_HiSeq_Joined %>% count(data_submitting_center_label)
missingSeqCenter <- rownames(metaQiitaWGS_RNA_HiSeq_Joined)[which(is.na(metaQiitaWGS_RNA_HiSeq_Joined$data_submitting_center_label))]
smallHopkins <- rownames(metaQiitaWGS_RNA_HiSeq_Joined)[metaQiitaWGS_RNA_HiSeq_Joined$data_submitting_center_label == "Johns Hopkins / University of Southern California"]
samples2Remove <- c(missingSeqCenter, na.omit(smallHopkins))
metaQiitaWGS_RNA_HiSeq_Filt <- droplevels(metaQiitaWGS_RNA_HiSeq_Joined[!(rownames(metaQiitaWGS_RNA_HiSeq_Joined) %in% samples2Remove),])
dim(metaQiitaWGS_RNA_HiSeq_Filt) # 15059    39
#---------------------------------------------#
# Using plate number for batch to decontaminate
#---------------------------------------------#

# Function for extracting last n characters from R string
# URL: https://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

# NB: metaQiitaWGS_RNA_HiSeq_Filt IS A METADATA FILE OF TCGA SAMPLES CONTAINING ALIQUOT IDS FOR ALL SAMPLES
# IN A COLUMN CALLED "tcga_aliquot_id"
tmp <- as.character(metaQiitaWGS_RNA_HiSeq_Filt$tcga_aliquot_id)
metaQiitaWGS_RNA_HiSeq_Filt$PlateCenter <- factor(substrRight(tmp, 7))
# NB: SINCE DECONTAM ESSENTIALLY PERFORMS LINEAR REGRESSION BETWEEN READ FRACTIONS AND 
# ANALYTE CONCENTRATIONS, AT LEAST 10 SAMPLES ARE REQUIRED PER PLATE-CENTER COMBINATION
# TO BE PROCESSED FOR IDENTIFYING PUTATIVE CONTAMINANTS. NOTE ALSO THAT ANY CONTAMINANT
# IDENTIFIED IN ANY ONE PLATE-CENTER BATCH WILL BE REMOVED FROM THE WHOLE DATASET
booleanPlateCenter <- as.logical(table(metaQiitaWGS_RNA_HiSeq_Filt$PlateCenter)>=10)
sufficientPlateCenter <- names(table(metaQiitaWGS_RNA_HiSeq_Filt$PlateCenter))[booleanPlateCenter]
length(sufficientPlateCenter) # 329
metaQiitaWGS_RNA_HiSeq_Filt$PlateCenterFlag <- (metaQiitaWGS_RNA_HiSeq_Filt$PlateCenter %in% sufficientPlateCenter)
metaQiitaWGS_RNA_HiSeq_Filt_PlateCenterSubset <- droplevels(metaQiitaWGS_RNA_HiSeq_Filt[metaQiitaWGS_RNA_HiSeq_Filt$PlateCenterFlag &
                                                                                          !is.na(metaQiitaWGS_RNA_HiSeq_Filt$aliquot_concentration),])
dim(metaQiitaWGS_RNA_HiSeq_Filt_PlateCenterSubset) # 14374    41

# NB: rep200Data_WGS_RNA_HiSeq_Fungi_PlateCenterSubset CONTAINS RAW TCGA FUNGI DATA FROM rep200
rep200Data_WGS_RNA_HiSeq_Fungi_PlateCenterSubset <- rep200Data_WGS_RNA_HiSeq_Fungi[rownames(metaQiitaWGS_RNA_HiSeq_Filt_PlateCenterSubset),]

# Decontam
require(decontam)
countDataPlateCenter <- rep200Data_WGS_RNA_HiSeq_Fungi_PlateCenterSubset
countMetadataPlateCenter <- metaQiitaWGS_RNA_HiSeq_Filt_PlateCenterSubset

contamdf.freq.fungi.plateCenter <- isContaminant(seqtab = as.matrix(countDataPlateCenter), 
                                                 conc = countMetadataPlateCenter$aliquot_concentration, 
                                                 method = "frequency", 
                                                 batch = countMetadataPlateCenter$PlateCenter,
                                                 threshold = 0.1) # DEFAULT VALUE IS 0.1
# save(contamdf.freq.fungi.plateCenter, file = "Interim_data/contamdf.freq.fungi.plateCenter_13Sep21.RData")
load("Interim_data/contamdf.freq.fungi.plateCenter_13Sep21.RData")
table(contamdf.freq.fungi.plateCenter$contaminant) # 71 TRUE | 248 FALSE
hist(contamdf.freq.fungi.plateCenter$p)
contaminants <- rownames(contamdf.freq.fungi.plateCenter)[contamdf.freq.fungi.plateCenter$contaminant]
contaminantsFungiPlateCenterTCGA <- rep200TaxSplit_Fungi[contaminants,]
# save(contaminantsFungiPlateCenterTCGA, file = "Interim_data/contaminantsFungiPlateCenterTCGA_13Sep21.RData")

notContamSumFreq <- colSums(as.matrix(rep200Data_WGS_RNA_HiSeq_Fungi)[,!contamdf.freq.fungi.plateCenter$contaminant])
contamSumFreq <- colSums(as.matrix(rep200Data_WGS_RNA_HiSeq_Fungi)[,contamdf.freq.fungi.plateCenter$contaminant])
sum(contamSumFreq)/sum(colSums(as.matrix(rep200Data_WGS_RNA_HiSeq_Fungi))) #--> 0.009430594 (13 Sep 21)

#------------------------Version 2 of decontamination------------------------#
## Save data for literature searching
# Added later. Goal is to cross-examine decontam results with
# biological plausibility (from the literature) and WIS results
load("Interim_data/shared_fungi_features_at_each_taxa_level_13Sep21.RData") ## Load shared features with Weizmann
rep200TaxSplit_Fungi_Paired_to_Weizmann <- read.csv("Supporting_data/rep200TaxSplit_Fungi_Paired_To_Weizmann_Final.csv", stringsAsFactors = FALSE, row.names = 1)
decontamResults <- contamdf.freq.fungi.plateCenter
decontamResults$species <- rep200TaxSplit_Fungi_Paired_to_Weizmann[rownames(decontamResults),"species"]
decontamResults$sharedWIS <- ifelse(decontamResults$species %in% sharedSpecies, yes = "YES", no = "NO")
# decontamResults %>% write.csv(file = "Interim_data/contaminantsFungiPlateCenterTCGA_updated_annotations_12Oct21.csv")

## Load data after literature searching is complete
decontamResultsV2 <- read.csv("Supporting_data/mycobiome_contaminant_analyses_updated_annotations_12Oct21.csv", row.names = 1, stringsAsFactors = FALSE)
# Extract PMIDs
require(stringr)
uniquePMIDs <- na.omit(unique(str_extract(string = decontamResultsV2$comments_and_literature, "[0-9]{7,8}")))
length(uniquePMIDs)
# Process rest of decontamV2
decontamResultsV2 <- decontamResultsV2[,!(colnames(decontamResultsV2) == "comments_and_literature")]
decontamResultsV2$decision <- ifelse(decontamResultsV2$shared_with_WIS == "YES" |
                                       decontamResultsV2$in_hmp_gut_mycobiome_metagenomic_data %in% c("YES","YES*") |
                                       decontamResultsV2$known_human_association_literature == "YES" |
                                       (decontamResultsV2$known_human_association_literature == "UNKNOWN" & 
                                          decontamResultsV2$decontam_predicted_contaminant == "FALSE"),
                                     yes = "KEEP", no = "DISCARD")
table(decontamResultsV2$decision) # KEEP 224 | DISCARD 95

contaminantsV2 <- rownames(decontamResultsV2)[which(decontamResultsV2$decision == "DISCARD")]

decontamResultsV2$reason <- decontamResultsV2$known_human_association_literature
decontamResultsV2$reason[decontamResultsV2$known_human_association_literature == "YES"] <- "Known human association"
decontamResultsV2$reason[decontamResultsV2$known_human_association_literature == "NO"] <- "No known human association"
decontamResultsV2$reason[decontamResultsV2$known_human_association_literature == "UNKNOWN" & decontamResultsV2$decontam_predicted_contaminant == "FALSE"] <- "Unknown human association but\nnot predicted contaminant"
decontamResultsV2$reason[decontamResultsV2$known_human_association_literature == "UNKNOWN" & decontamResultsV2$decontam_predicted_contaminant == "TRUE"] <- "Unknown human association and\npredicted contaminant"
decontamResultsV2$reason[decontamResultsV2$in_hmp_gut_mycobiome_metagenomic_data %in% c("YES","YES*")] <- "In HMP gut mycobiome data"
decontamResultsV2$reason[decontamResultsV2$shared_with_WIS == "YES"] <- "Shared with WIS"
decontamResultsV2$reason <- factor(decontamResultsV2$reason, levels = c("Shared with WIS","In HMP gut mycobiome data","Known human association",
                                                                        "Unknown human association but\nnot predicted contaminant",
                                                                        "Unknown human association and\npredicted contaminant",
                                                                        "No known human association"))

save(decontamResultsV2, contaminantsV2, file = "Interim_data/decontamResultsV2_13Oct21.RData")

contaminantsV2boolean <- ifelse(colnames(rep200Data_WGS_RNA_HiSeq_Fungi) %in% contaminantsV2,
                                yes = TRUE, no = FALSE)
notContamSumFreq <- colSums(as.matrix(rep200Data_WGS_RNA_HiSeq_Fungi)[,!contaminantsV2boolean])
contamSumFreq <- colSums(as.matrix(rep200Data_WGS_RNA_HiSeq_Fungi)[,contaminantsV2boolean])
sum(contamSumFreq)/sum(colSums(as.matrix(rep200Data_WGS_RNA_HiSeq_Fungi))) #--> 0.02234265 (14 Oct 21)



#-------------------------------------------------------------#

rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2 <- rep200Data_WGS_RNA_HiSeq_Fungi[,!(colnames(rep200Data_WGS_RNA_HiSeq_Fungi) %in% contaminantsV2)]
dim(rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2) # 15088   224

sum(rowSums(rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2) == 0) # 564
emptySamplesAfterDecontamV2 <- names(which(rowSums(rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2) == 0))
samples2RemoveWithDecontamV2 <- c(emptySamplesAfterDecontamV2, missingSeqCenter, na.omit(smallHopkins))

rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2[!(rownames(rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2) %in% samples2RemoveWithDecontamV2),]
metaQiitaWGS_RNA_HiSeq_Filt_Nonzero_DecontamV2 <- droplevels(metaQiitaWGS_RNA_HiSeq_Filt[!(rownames(metaQiitaWGS_RNA_HiSeq_Filt) %in% samples2RemoveWithDecontamV2),])

dim(rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero) # 14495   224
dim(metaQiitaWGS_RNA_HiSeq_Filt_Nonzero_DecontamV2) # 14495    41

## Convert strings to factors
metaQiitaWGS_RNA_HiSeq_Filt_Nonzero_DecontamV2$sample_type <- factor(metaQiitaWGS_RNA_HiSeq_Filt_Nonzero_DecontamV2$sample_type)
metaQiitaWGS_RNA_HiSeq_Filt_Nonzero_DecontamV2$disease_type <- factor(metaQiitaWGS_RNA_HiSeq_Filt_Nonzero_DecontamV2$disease_type)
metaQiitaWGS_RNA_HiSeq_Filt_Nonzero_DecontamV2$experimental_strategy <- factor(metaQiitaWGS_RNA_HiSeq_Filt_Nonzero_DecontamV2$experimental_strategy)
metaQiitaWGS_RNA_HiSeq_Filt_Nonzero_DecontamV2$data_submitting_center_label <- factor(metaQiitaWGS_RNA_HiSeq_Filt_Nonzero_DecontamV2$data_submitting_center_label)
metaQiitaCombined_Nonzero_DecontamV2 <- metaQiitaWGS_RNA_HiSeq_Filt_Nonzero_DecontamV2

save(rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2,
     file = "Interim_data/fungi_data_WGS_RNA_HiSeq_DecontamV2_13Oct21.RData")

#------------------------------------------------------------------#
# WGS vs. RNA-Seq norm read counts
#------------------------------------------------------------------#

metaQiitaCombined_Nonzero_DecontamV2_DNAonly <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(experimental_strategy == "WGS") %>% droplevels()
metaQiitaCombined_Nonzero_DecontamV2_RNAonly <- metaQiitaCombined_Nonzero_DecontamV2 %>% filter(experimental_strategy == "RNA-Seq") %>% droplevels()

wgsData <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_DNAonly),]
rnaData <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero[rownames(metaQiitaCombined_Nonzero_DecontamV2_RNAonly),]

wgsSampleCounts <- log10(unname(rowSums(wgsData)))
rnaSampleCounts <- log10(unname(rowSums(rnaData)))
summary(wgsSampleCounts)
summary(rnaSampleCounts)

combinedSampleCounts <- data.frame(sample_counts = c(wgsSampleCounts,rnaSampleCounts),
                                   data_type = c(rep("WGS",length(wgsSampleCounts)), rep("RNA-Seq",length(rnaSampleCounts))))
require(EnvStats)
combinedSampleCounts %>%
  ggviolin(x = "data_type",
           y = "sample_counts",
           fill = "data_type",
           palette = "nejm",
           legend = "none",
           draw_quantiles = c(0.25,0.50,0.75),
           xlab = "Experimental strategy",
           ylab = "log10(decontaminated sample fungi read counts)",
           add = "mean",
           add.params = list(color="white",size=1)) +
  stat_compare_means(label.x.npc = 0.1, label.y = 7) +
  stat_n_text() + 
  ggsave(filename = "Figures/Supplementary_Figures/tcga_read_count_wgs_vs_rna_decontamV2.jpeg",
         dpi = "retina",
         height = 5,
         width = 3,
         units = "in")

#-----------------------------------------------#
# Voom-SNM
#-----------------------------------------------#
require(limma)
require(edgeR)
require(dplyr)
require(snm)

qcMetadata <- metaQiitaCombined_Nonzero_DecontamV2 
qcData <- rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero 

qcMetadata %>% count(data_submitting_center_label)

# Set up design matrix
covDesignNorm <- model.matrix(~0 + disease_type + sample_type +
                                data_submitting_center_label +
                                experimental_strategy,
                              data = qcMetadata)
# Check row dimensions
dim(covDesignNorm)[1] == dim(qcData)[1]

print(colnames(covDesignNorm))
colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
print(colnames(covDesignNorm))

# Set up counts matrix
counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)

# Normalize using edgeR and then plug into voom
dge <- DGEList(counts = counts)
vdge_data <- voom(dge, design = covDesignNorm, plot = TRUE, save.plot = TRUE,
                  normalize.method="quantile")
vdge_dataE_DecontamV2 <- t(vdge_data$E)

# Apply
bio.var <- model.matrix(~disease_type + sample_type,
                        data=qcMetadata)
adj.var <- model.matrix(~data_submitting_center_label +
                          experimental_strategy,
                        data=qcMetadata)
colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
print(dim(adj.var))
print(dim(bio.var))
print(dim(t(vdge_data$E)))
print(dim(covDesignNorm))

snmDataObjOnly <- snm(raw.dat = vdge_data$E, 
                      bio.var = bio.var, 
                      adj.var = adj.var, 
                      rm.adj=TRUE,
                      verbose = TRUE,
                      diagnose = TRUE)
snmDataOGUFungiDecontamV2 <- t(snmDataObjOnly$norm.dat)

save(snmDataOGUFungiDecontamV2,
     vdge_dataE_DecontamV2,
     rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero,
     metaQiitaCombined_Nonzero_DecontamV2,
     file = "Interim_data/snmDataFungi_DecontamV2_13Oct21.RData")

#---------------------------------------------
# PCA
#---------------------------------------------

source("00-Functions.R") # contains pcaPlotting() function

voomPca_SeqCenter <- pcaPlotting(pcaObject = prcomp(vdge_dataE_DecontamV2),
                                pcChoices = c(1,2),
                                dataLabels = qcMetadata$data_submitting_center_label,
                                factorString = "Sequencing Center",
                                titleString = "PCA w/o Batch Correction") +
  ggsave(filename = "Figures/Supplementary_Figures/tcga_pca_voom_only_seq_center_DecontamV2.jpeg",
              dpi = "retina", width = 6, height = 6, units = "in")

voomPca_ExpStrategy <- pcaPlotting(pcaObject = prcomp(vdge_dataE_DecontamV2),
                    pcChoices = c(1,2),
                    dataLabels = qcMetadata$experimental_strategy,
                    factorString = "Experimental Strategy",
                    titleString = "PCA w/o Batch Correction") +
  ggsave(filename = "Figures/Supplementary_Figures/tcga_pca_voom_only_exp_strategy_DecontamV2.jpeg",
         dpi = "retina", width = 6, height = 6, units = "in")

vsnmPca_SeqCenter <- pcaPlotting(pcaObject = prcomp(snmDataOGUFungiDecontamV2),
                       pcChoices = c(1,2),
                       dataLabels = qcMetadata$data_submitting_center_label,
                       factorString = "Sequencing Center",
                       titleString = "PCA w/ Voom-SNM Correction") +
  ggsave(filename = "Figures/Supplementary_Figures/tcga_pca_vsnm_seq_center_DecontamV2.jpeg",
         dpi = "retina", width = 6, height = 6, units = "in")

vsnmPca_ExpStrategy <- pcaPlotting(pcaObject = prcomp(snmDataOGUFungiDecontamV2),
                       pcChoices = c(1,2),
                       dataLabels = qcMetadata$experimental_strategy,
                       factorString = "Sequencing Center",
                       titleString = "PCA w/ Voom-SNM Correction") +
  ggsave(filename = "Figures/Supplementary_Figures/tcga_pca_vsnm_exp_strategy_DecontamV2.jpeg",
         dpi = "retina", width = 6, height = 6, units = "in")

#------------------------------------------------------------------#
# PVCA
#------------------------------------------------------------------#

# The following scripts were run on a compute cluster to generate results:
# "Supplementary_scripts/S01B-Run-pvca-fungi.R"
# "Supplementary_scripts/S02-pvca-function.R"

# load("Interim_data/pvca_fungi_results_raw_Voom_VSNM_13Sep21.RData")
# 
# pvcaRawRound <- round(pvcaRaw,3)
# pvcaVoomRound <- round(pvcaVoom,3)
# pvcaVSNMRound <- round(pvcaVSNM,3)
# 
# pvcaRes <- data.frame('Sample Type' = c(pvcaRawRound[1], pvcaVoomRound[1], pvcaVSNMRound[1]),
#                       'Disease Type' = c(pvcaRawRound[2], pvcaVoomRound[2], pvcaVSNMRound[2]),
#                       'Sequencing Center' = c(pvcaRawRound[3], pvcaVoomRound[3], pvcaVSNMRound[3]),
#                       'Experimental Strategy' = c(pvcaRawRound[4], pvcaVoomRound[4], pvcaVSNMRound[4]),
#                       'Residual\n(not explained by\ntechnical variation)' = c(pvcaRawRound[5], pvcaVoomRound[5], pvcaVSNMRound[5]),
#                       data_type = factor(c("Raw count data","Voom Normalized Data","Voom Normalized & SNM Corrected Data"),
#                                          levels = c("Raw count data","Voom Normalized Data","Voom Normalized & SNM Corrected Data")),
#                       check.names = FALSE)
# 
# pvcaRes.melted <- reshape2::melt(pvcaRes, id.vars = "data_type")
# pvcaRes.melted %>%
#   ggbarplot(x = "variable",
#             y = "value",
#             fill = "data_type",
#             palette = "nejm",
#             legend = "top",
#             ylim = c(0,1),
#             xlab = "Biological Effects & Technical Effects",
#             ylab = "Weighted average proportion variance",
#             label = TRUE,
#             position = position_dodge(0.9)) +
#   labs(fill = "Data type") +
#   ggsave("Figures/Supplementary_Figures/pvca_plot_OGUs_13Sep21.jpeg",
#          dpi = "retina",
#          width = 12,
#          height = 3,
#          units = "in")
