#-----------------------------------------------------------------------------
# 00-Functions.R
# Copyright (c) 2021--, Greg Poore
# Purpose: List main functions for R scripts
#-----------------------------------------------------------------------------

pcaPlotting <- function(pcaObject,pcChoices, dataLabels, factorString, titleString){
  require(ggbiplot)
  theme_update(plot.title = element_text(hjust = 0.5))
  g <- ggbiplot(pcaObject,pcChoices, obs.scale = 1, var.scale = 1,
                groups = dataLabels, ellipse = TRUE,
                alpha = 0.05,
                circle = TRUE,var.axes=FALSE) + 
    # scale_color_nejm(name = factorString) +
    theme_bw() + 
    #theme(legend.direction = "horizontal", legend.position = "top") +
    ggtitle(titleString) + theme(plot.title = element_text(hjust = 0.5))
  
  print(g)
}

avgRAbarplot <- function(psObj, dataType, filename, yAxisLab="Mean Relative Abundance",
                         title="TCGA bacteria vs. fungi average relative abundances"){
  figFilepath <- "Figures/Supplementary_Figures/"
  dataFilepath <- "Figures_data/Supplementary_Figures/"
  ps <- tax_glom(psObj, "Domain")
  ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
  ps1 <- merge_samples(ps0, "investigation")
  ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
  ps2_MeltFilt <- psmelt(ps2) %>% 
    select(OTU, Sample, Abundance, Domain) %>% 
    mutate(label = ifelse(Domain == "k__Eukaryota",yes = round(Abundance, 3),no = NA)) %>%
    droplevels() %>% data.frame()
  plot_bar(ps2, fill="Domain") + scale_fill_manual(values = c("#0072B5FF","#BC3C29FF")) +
    ylab(yAxisLab) + xlab("TCGA Cancer Type") + 
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
    ggtitle(paste0(title," (Primary Tumor | ",dataType,")")) +
    geom_text(data = ps2_MeltFilt, aes(x=Sample, label=label), 
            color = "white", position=position_dodge(width=0.9), hjust=-0.25,vjust=-0.01,angle = 90) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0(figFilepath,filename,".pdf"),
         dpi = "retina", width = 12, height = 4, units = "in")
  ps2_otu <- data.frame(otu_table(ps2))
  ps2_otu %>% write.csv(file = paste0(dataFilepath,filename,".csv"))
  ps2Df <- data.frame(otu_table(ps2))
  print("Normalized relative abundance means (%):")
  print(round(100*colMeans(ps2Df),2))
  return(ps2)
}

vsnmFunctionTCGA <- function(qcData, qcMetadata=metaQiitaCombined_Nonzero_DecontamV2, cancerTypeFlag=FALSE, filename){
  require(limma)
  require(edgeR)
  require(dplyr)
  require(snm)
  # Set up design matrix
  if(cancerTypeFlag){
    covDesignNorm <- model.matrix(~0 + disease_type + sample_type +
                                  data_submitting_center_label +
                                  experimental_strategy,
                                data = qcMetadata)
    } else{
      covDesignNorm <- model.matrix(~0 + sample_type +
                                  data_submitting_center_label +
                                  experimental_strategy,
                                data = qcMetadata)
    }
  
  # Check row dimensions
  print(dim(covDesignNorm))
  print(dim(qcData))
  print(dim(covDesignNorm)[1] == dim(qcData)[1])
  colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
  
  # Set up counts matrix
  counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)
  
  # Normalize using edgeR and then plug into voom
  dge <- DGEList(counts = counts)
  vdge_data <- voom(dge, design = covDesignNorm, plot = TRUE, save.plot = TRUE,
                    normalize.method="quantile")
  vdge_dataE <- t(vdge_data$E)
  
  # Apply
  if(cancerTypeFlag){
    bio.var <- model.matrix(~disease_type + sample_type,
                          data=qcMetadata)
    } else{
      bio.var <- model.matrix(~sample_type,
                          data=qcMetadata)
    }
  adj.var <- model.matrix(~data_submitting_center_label +
                            experimental_strategy,
                          data=qcMetadata)
  colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
  colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
  
  snmDataObjOnly <- snm(raw.dat = vdge_data$E, 
                        bio.var = bio.var, 
                        adj.var = adj.var, 
                        rm.adj=TRUE,
                        verbose = TRUE,
                        diagnose = TRUE)
  snmData <- t(snmDataObjOnly$norm.dat)
  res <- list(qcData=qcData,
              vdge_dataE=vdge_dataE,
              snmData=snmData)
  return(res)
}

findSharedTaxa <- function(psTCGA, psWeizmann, taxLevel){
  if(taxLevel == "species"){
    tcgaFeat <- colnames(psFungiHiSeqFungi_Paired2Wz_species) # rownames(otu_table(psTCGA))
    wzFeat <- data.frame(tax_table(psWeizmann))[,taxLevel]
  } else{
    tcgaFeat <- data.frame(tax_table(psTCGA))[,taxLevel] # rownames(otu_table(psTCGA))
    wzFeat <- data.frame(tax_table(psWeizmann))[,taxLevel]
  }
  tcgaFeatFilt <- tcgaFeat[!grepl("other",tcgaFeat)]
  wzFeatFilt <- wzFeat[!grepl("other",wzFeat)]
  sharedFeatures <- intersect(tcgaFeatFilt, wzFeatFilt)
  print(sharedFeatures)
  return(sharedFeatures)
}

findSharedTaxaWISRep200 <- function(rep200TaxTableX, psWeizmann, taxLevel){
  wzFeat <- data.frame(tax_table(psWeizmann))[,taxLevel]
  wzFeatFilt <- wzFeat[!grepl("other",wzFeat)]
  rep200TaxTableXFilt <- rep200TaxTableX[!grepl("other",rep200TaxTableX)]
  sharedFeatures <- intersect(rep200TaxTableXFilt, wzFeatFilt)
  print(sharedFeatures)
  return(sharedFeatures)
}

findSharedTaxaBacteria <- function(psTCGA, psWeizmann, taxLevel){

  tcgaFeat <- data.frame(tax_table(psTCGA))
  wzFeat <- data.frame(tax_table(psWeizmann))

  tcgaFeat$formatted <- gsub('([[:punct:]])|\\s+','',tcgaFeat[,taxLevel])
  wzFeat$formatted <- gsub('([[:punct:]])|\\s+','',wzFeat[,taxLevel])

  tcgaFeatFilt1 <- tcgaFeat[!tcgaFeat[,taxLevel] == "other",]
  wzFeatFilt1 <- wzFeat[!wzFeat[,taxLevel] == "other",]

  sharedFeatures <- intersect(tcgaFeatFilt1$formatted, wzFeatFilt1$formatted)
  tcgaFeatFilt2 <- tcgaFeatFilt1[tcgaFeatFilt1$formatted %in% sharedFeatures,]
  wzFeatFilt2 <- wzFeatFilt1[wzFeatFilt1$formatted %in% sharedFeatures,]

  print(unique(tcgaFeatFilt2[,taxLevel]))
  print(sprintf("Number of uniquely shared features: %d",length(unique(tcgaFeatFilt2[,taxLevel]))))

  res <- list(intersectedOGUs = rownames(tcgaFeatFilt2),
    intersectedTaxa = tcgaFeatFilt2[,taxLevel])

  return(res)
}

findPrevGenus <- function(metaData=metaQiitaCombined_Nonzero_AllSeqPlatforms_8cancer_GenusShared, 
                     countData=rep200FungiAllSeqPlatformsGenusShared_8cancer_Nonzero, 
                     sampleType=c("Primary Tumor","Solid Tissue Normal"), 
                     cancerType=c("Breast Cancer"), seqCenter=c("Harvard Medical School"),
                     taxaFlag=FALSE, weizmannIntersectFlag=FALSE, weizmannTissue=c("breast"),
                     weizmannTissueType=c("tumor","nat"),
                     weizmannMeta=weizmannMetaGenusShared_Nonzero, weizmannData=weizmannGenusShared_Nonzero){
  metaData_seqCenter_Dz <- metaData %>% filter(data_submitting_center_label %in% seqCenter) %>%
    filter(sample_type %in% sampleType) %>% filter(disease_type %in% cancerType) %>% droplevels()
  countData_subset <- countData[rownames(metaData_seqCenter_Dz),]
  taxaBoolean <- (colSums(countData_subset) != 0)
  taxaWithCounts <- names(colSums(countData_subset))[taxaBoolean]
  taxaWithCountsNamed <- taxaWithCounts # rep200FungiTaxIntersect[taxaWithCounts, "Genus"]
  
  print(sprintf("TCGA taxa with counts: %d", sum(taxaBoolean)))
  # print(taxaWithCountsNamed)
  
  taxaWithCountsNamedFormatted <- taxaWithCountsNamed # gsub('^g__|([[:punct:]])','',taxaWithCountsNamed)
  if(taxaFlag){print(sort(taxaWithCountsNamedFormatted))}
  if(weizmannIntersectFlag){
    weizmannMeta_Dz <- weizmannMeta %>% filter(tissue %in% weizmannTissue) %>% filter(type.detail %in% weizmannTissueType) %>% droplevels()
    weizmannData_Dz <- weizmannData[rownames(weizmannMeta_Dz),]
    weizmannTaxaWithCounts <- names(colSums(weizmannData_Dz))[colSums(weizmannData_Dz) != 0]
    # print(weizmannTaxaWithCounts)
    weizmannTaxaWithCountsFormatted <- weizmannTaxaWithCounts # gsub('([[:punct:]])', ' ', gsub('^g__','',weizmannTaxaWithCounts))
    # weizmannTaxaWithCountsFormatted <- gsub('^s__', '', gsub('([[:punct:]])',' ',weizmannTaxaWithCounts))
    
    # print(weizmannTaxaWithCountsFormatted)
    taxaIntersection <- intersect(taxaWithCountsNamedFormatted, weizmannTaxaWithCountsFormatted)
    cat(sprintf("\nNumber of intersected taxa: %d/%d (Weizmann) | %d/%d (TCGA)\n",
                length(taxaIntersection), length(weizmannTaxaWithCountsFormatted),
                length(taxaIntersection), length(taxaWithCountsNamedFormatted)))
    print(taxaIntersection)
  }
}

findPrevSpecies <- function(metaData=metaQiitaCombined_Nonzero_AllSeqPlatforms_8cancer_SpeciesShared, 
                     countData=rep200FungiAllSeqPlatformsSpeciesShared_8cancer_Nonzero, 
                     sampleType=c("Primary Tumor","Solid Tissue Normal"), 
                     cancerType=c("Breast Cancer"), seqCenter=c("Harvard Medical School"),
                     taxaFlag=FALSE, weizmannIntersectFlag=FALSE, weizmannTissue=c("breast"),
                     weizmannTissueType=c("tumor","nat"),
                     weizmannMeta=weizmannMetaSpeciesShared_Nonzero, weizmannData=weizmannSpeciesShared_Nonzero){
  metaData_seqCenter_Dz <- metaData %>% filter(data_submitting_center_label %in% seqCenter) %>%
    filter(sample_type %in% sampleType) %>% filter(disease_type %in% cancerType) %>% droplevels()
  countData_subset <- countData[rownames(metaData_seqCenter_Dz),]
  taxaBoolean <- (colSums(countData_subset) != 0)
  taxaWithCounts <- names(colSums(countData_subset))[taxaBoolean]
  taxaWithCountsNamed <- taxaWithCounts # rep200FungiTaxIntersect[taxaWithCounts, "Species"]
  
  print(sprintf("TCGA taxa with counts: %d", sum(taxaBoolean)))
  
  taxaWithCountsNamedFormatted <- taxaWithCountsNamed # gsub('^s__|([[:punct:]])','',taxaWithCountsNamed)
  if(taxaFlag){print(sort(taxaWithCountsNamedFormatted))}
  if(weizmannIntersectFlag){
    weizmannMeta_Dz <- weizmannMeta %>% filter(tissue %in% weizmannTissue) %>% filter(type.detail %in% weizmannTissueType) %>% droplevels()
    weizmannData_Dz <- weizmannData[rownames(weizmannMeta_Dz),]
    weizmannTaxaWithCounts <- names(colSums(weizmannData_Dz))[colSums(weizmannData_Dz) != 0]
    weizmannTaxaWithCountsFormatted <- weizmannTaxaWithCounts # gsub('([[:punct:]])', ' ', gsub('^s__','',weizmannTaxaWithCounts))
    # weizmannTaxaWithCountsFormatted <- gsub('^s__', '', gsub('([[:punct:]])',' ',weizmannTaxaWithCounts))
    # print(weizmannTaxaWithCountsFormatted)
    taxaIntersection <- intersect(taxaWithCountsNamedFormatted, weizmannTaxaWithCountsFormatted)
    cat(sprintf("\nNumber of intersected taxa: %d/%d (Weizmann) | %d/%d (TCGA)\n", 
                  length(taxaIntersection), length(weizmannTaxaWithCountsFormatted),
                  length(taxaIntersection), length(taxaWithCountsNamedFormatted)))
    print(taxaIntersection)
  }
}

compareTvsNATbrayCurtis <- function(featureData, metaData = metaQiitaCombined_Nonzero_8cancer,
                                    # taxTable = psRep200FungiTax,
                                    scalar=10^6, snmFlag = TRUE, axes=c(1,2),
                                    grepFlag = TRUE, greplString = "Lung|Melanoma|Breast|Colorectal",
                                    plotlyAspectmode='cube', apiUploadFlag = FALSE, brayHullFlag = TRUE){
  require(tibble)
  require(stringr)
  require(plotly)
  # metaData processing
  metaDataFilt <- droplevels(metaData[metaData$sample_type %in% c("Primary Tumor","Solid Tissue Normal"),])
  metaDataFilt$sample_type <- ifelse(metaDataFilt$sample_type == "Primary Tumor", yes = "PT", no = "NAT")
  # metaDataFilt$sample_type[metaDataFilt$sample_type == "Primary Tumor"] <- "PT"
  # metaDataFilt$sample_type[metaDataFilt$sample_type == "Solid Tissue Normal"] <- "NAT"
  metaDataFilt$disease_and_sample_type <- paste(metaDataFilt$disease_type,metaDataFilt$sample_type)
  print(table(metaDataFilt$disease_and_sample_type))
  
  if(snmFlag){
    countData <- floor(((2^featureData)/rowSums(2^featureData))*scalar)
  } else{
    countData <- featureData
  }
  countDataRA <- countData/rowSums(countData)
  countDataRAfilt <- countDataRA[rownames(metaDataFilt),]
  
  # Calculate mean 
  groupedCountDataRA <- aggregate(countDataRAfilt,
                                  by = list(metaDataFilt$disease_and_sample_type),
                                  FUN = mean)
  # rownames(groupedCountDataRA) <- levels(metaDataFilt$disease_and_sample_type)
  groupedCountDataRAformatted <- groupedCountDataRA %>% column_to_rownames("Group.1")
  
  ptNatVec <- ifelse(grepl("NAT",rownames(groupedCountDataRAformatted)), yes="NAT",no="PT")
  dzSplitVec <- sapply(strsplit(rownames(groupedCountDataRAformatted)," NAT| PT"), `[`, 1) 
  groupedMetaDataFilt <- data.frame(disease_type = dzSplitVec,
                                    sample_type = ptNatVec,
                                    disease_and_sample_type = levels(factor(metaDataFilt$disease_and_sample_type)))
  rownames(groupedMetaDataFilt) <- levels(factor(metaDataFilt$disease_and_sample_type))
  
  # print(groupedMetaDataFilt)
  
  if(grepFlag){
    groupedMetaDataFiltGrep <- droplevels(groupedMetaDataFilt[grepl(greplString,groupedMetaDataFilt$disease_type),])
    groupedCountDataRAformattedGrep <- groupedCountDataRAformatted[rownames(groupedMetaDataFiltGrep),]
  } else{
    groupedMetaDataFiltGrep <- groupedMetaDataFilt
    groupedCountDataRAformattedGrep <- groupedCountDataRAformatted
  }
  
  combinedPSobj <- phyloseq(otu_table(groupedCountDataRAformattedGrep, taxa_are_rows = FALSE), 
                            sample_data(groupedMetaDataFiltGrep))

  bray_dist = phyloseq::distance(combinedPSobj, method="bray") # need to rarefy
  bray_ordination = ordinate(combinedPSobj, method="PCoA", distance=bray_dist)
  brayPlot2D <- plot_ordination(combinedPSobj, bray_ordination, color="disease_type", shape = "sample_type", axes = axes) +
    theme(aspect.ratio=1) + theme_pubr() + geom_point(size = 3) + scale_color_nejm() + coord_fixed() + theme(legend.position = "right") +
    labs(shape="Sample type", color="Disease type")

  if(brayHullFlag){
    brayPlot2D <- brayPlot2D + geom_mark_hull(aes(group = disease_type, color = disease_type), expand = 0)
  }
  
  brayData <- data.frame(bray_ordination$vectors, groupedMetaDataFiltGrep)
  brayPlot3D <- plot_ly(data = brayData, x= ~Axis.1, y= ~Axis.2, z= ~Axis.3) %>%
    add_markers(color= ~disease_type, symbol = ~sample_type, symbols = c(200, 18)) %>%
    layout(scene = list(xaxis = list(title = 'PCoA 1'),
                        yaxis = list(title = 'PCoA 2'),
                        zaxis = list(title = 'PCoA 3'),
                        aspectmode=plotlyAspectmode))
  if(apiUploadFlag){
    api_create(brayPlot3D)
  }
  
  print(brayPlot3D)
  
  res <- list(brayPlot2D=brayPlot2D,
              brayPlot3D=brayPlot3D,
              bray_dist=bray_dist,
              bray_ordination=bray_ordination,
              brayData=brayData,
              countDataRAfilt=countDataRAfilt,
              metaDataFilt=metaDataFilt,
              groupedCountDataRAformattedGrep=groupedCountDataRAformattedGrep,
              groupedMetaDataFiltGrep=groupedMetaDataFiltGrep,
              combinedPSobj=combinedPSobj)
  return(res)
}

calcCorrAlpha <- function(psAll, seqCenter="Harvard Medical School", sampleType="Primary Tumor",
                          fungi_gOTUs=fungiOGUs, bacteria_gOTUs=bacteriaOGUs,
                          dataType="WGS", rarefyNumber=NA, statCorMethod = "Pearson", libraryCorrFlag = FALSE,
                          filePath = "Figures/Supplementary_Figures/"){
  if(statCorMethod == "Pearson"){
    corCoefName <- "R"
  } else if(statCorMethod == "Spearman"){
    corCoefName <- "rho"
  }
  # print(seqCenter)
  if(seqCenter=="All Seq Centers"){
    filtBoolean <- ((sample_data(psAll)[["sample_type"]] == sampleType) &
                      (sample_data(psAll)[["experimental_strategy"]] == dataType))
    figWidthPerCT <- 18
    figRowsPerCT <- 4
  } else{
    filtBoolean <- ((sample_data(psAll)[["data_submitting_center_label"]] == seqCenter) &
                      (sample_data(psAll)[["sample_type"]] == sampleType) &
                      (sample_data(psAll)[["experimental_strategy"]] == dataType))
    figWidthPerCT <- 12
    figRowsPerCT <- 3
  }
  
  psAll_Filt <- prune_samples(filtBoolean, psAll)
  
  if(!is.na(rarefyNumber)){
    rarefyLevelInt <- rarefyNumber
  } else{
    print("Sample read distribution:")
    print(summary(sample_sums(psAll_Filt)))
    rarefyLevel <- readline(prompt="Enter rarefaction amount for sample data: ")
    rarefyLevelInt <- as.integer(rarefyLevel)
  }
  print("Ok. Rarefying sample data...")
  psAll_Filt_rarefied <- rarefy_even_depth(psAll_Filt, sample.size = rarefyLevelInt,
                                             rngseed = 42, replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
  rarefactionTitle <- sprintf("Rarefied %d reads",rarefyLevelInt)
  plotTitlePerCT <- paste(rarefactionTitle, seqCenter, sampleType, dataType, statCorMethod, sep = " | ")
  plotTitleAggregated <- paste("Aggregated data",rarefactionTitle, seqCenter, sampleType, dataType, statCorMethod, sep = " | ")
  
  taxFungiBoolean <- (rownames(tax_table(psAll_Filt_rarefied)) %in% fungi_gOTUs)
  taxBacteriaBoolean <- (rownames(tax_table(psAll_Filt_rarefied)) %in% bacteria_gOTUs)
  
  psFungi_Filt_rarefied <- prune_taxa(taxFungiBoolean, psAll_Filt_rarefied)
  psBacteria_Filt_rarefied <- prune_taxa(taxBacteriaBoolean, psAll_Filt_rarefied)
  
  meta2Keep <- sample_data(psFungi_Filt_rarefied)[,c("experimental_strategy","disease_type","sample_type",
                                            "cgc_platform","investigation","data_submitting_center_label",
                                            "pathologic_stage_label","library_size","library_size_log10")]
  print("Now calculating richness...")
  fungi_rarefied_alpha <- data.frame(estimate_richness(psFungi_Filt_rarefied, measures = c("Observed","Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")))
  colnames(fungi_rarefied_alpha) <- paste0("fungi_",colnames(fungi_rarefied_alpha))
  fungi_rarefied_alpha_formatted <- fungi_rarefied_alpha %>% rownames_to_column("sampleid")
  bacteria_rarefied_alpha <- data.frame(estimate_richness(psBacteria_Filt_rarefied, measures = c("Observed","Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")))
  colnames(bacteria_rarefied_alpha) <- paste0("bacteria_",colnames(bacteria_rarefied_alpha))
  bacteria_rarefied_alpha_formatted <- bacteria_rarefied_alpha %>% rownames_to_column("sampleid")
  
  fungi_rarefied_alpha_formatted %>%
    left_join(bacteria_rarefied_alpha_formatted, by = "sampleid") %>%
    column_to_rownames("sampleid") -> fungi_and_bacteria_alpha
  
  rownames(meta2Keep) <- paste0("X",rownames(meta2Keep))
  fungi_and_bacteria_alpha_with_meta <- cbind(fungi_and_bacteria_alpha, meta2Keep[rownames(fungi_and_bacteria_alpha),])
  
  maxFungiObserved <- max(fungi_and_bacteria_alpha_with_meta$fungi_Observed)
  fungi_and_bacteria_alpha_with_meta %>%
    ggplot(aes(bacteria_Observed,fungi_Observed)) +
    geom_point() +
    xlab("Bacterial richness") + ylab("Fungal richness") + theme_pubr() +
    facet_wrap(facets = vars(investigation), nrow = figRowsPerCT) +
    scale_y_continuous(limits = c(0, maxFungiObserved+10)) +
    geom_smooth(method='lm') +
    ggtitle(plotTitlePerCT) +
    theme(plot.title = element_text(hjust = 0.5)) +
    # scale_x_log10() +
    # rotate_x_text(30) +
    scale_fill_igv() +
    stat_cor(method = tolower(statCorMethod), cor.coef.name = corCoefName) -> richnessCorrPlotPerCT
  
  print(richnessCorrPlotPerCT)
  baseName <- paste(paste0("rarefied",rarefyLevelInt),
                    gsub('([[:punct:]])|\\s+','',seqCenter), 
                    gsub('([[:punct:]])|\\s+','',sampleType),
                    gsub('([[:punct:]])|\\s+','',dataType), 
                    statCorMethod, sep = "_")
  richnessCorrPlotPerCT + ggsave(filename = paste0("alpha_div_per_CT_",baseName,".png"), path = filePath,
                                 dpi = "retina", units = "in", width = figWidthPerCT)
  
  fungi_and_bacteria_alpha_with_meta %>%
    ggplot(aes(bacteria_Observed,fungi_Observed)) +
    geom_point() +
    xlab("Bacterial richness") + ylab("Fungal richness") + theme_pubr() +
    scale_y_continuous(limits = c(0, NA)) +
    geom_smooth(method='lm') +
    ggtitle(plotTitleAggregated) +
    theme(plot.title = element_text(hjust = 0.5)) +
    # scale_x_log10() +
    # rotate_x_text(30) +
    scale_fill_igv() +
    stat_cor(method = tolower(statCorMethod), cor.coef.name = corCoefName) -> richnessCorrPlotAggregated
  
  print(richnessCorrPlotAggregated)
  richnessCorrPlotAggregated + ggsave(filename = paste0("alpha_div_All_CT_",baseName,".png"), path = filePath,
                                      dpi = "retina", units = "in", width = 12)
  
  res <- list(fungi_and_bacteria_alpha_with_meta=fungi_and_bacteria_alpha_with_meta,
              richnessCorrPlotPerCT=richnessCorrPlotPerCT,
              richnessCorrPlotAggregated=richnessCorrPlotAggregated,
              meta2Keep=meta2Keep)
  
  if(libraryCorrFlag){
    
    plotTitlePerCT_Library_Fungi <- paste("Fungi vs. Library Size", rarefactionTitle, seqCenter, sampleType, dataType, statCorMethod, sep = " | ")
    plotTitleAggregated_Library_Fungi <- paste("Fungi vs. Library Size", "Aggregated data",rarefactionTitle, seqCenter, sampleType, dataType, statCorMethod, sep = " | ")
    plotTitlePerCT_Library_Bacteria <- paste("Bacteria vs. Library Size", rarefactionTitle, seqCenter, sampleType, dataType, statCorMethod, sep = " | ")
    plotTitleAggregated_Library_Bacteria <- paste("Bacteria vs. Library Size", "Aggregated data",rarefactionTitle, seqCenter, sampleType, dataType, statCorMethod, sep = " | ")
    
    ## Fungi vs. library size
    fungi_and_bacteria_alpha_with_meta %>%
      ggplot(aes(library_size_log10,fungi_Observed)) +
      geom_point(aes(colour = library_size_log10)) +
      scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
      xlab("Original microbial library Size (log10)") + ylab("Fungal richness") + theme_pubr() +
      # facet_wrap(facets = vars(investigation), scales = "free", nrow = 3) +
      facet_wrap(facets = vars(investigation), nrow = 3) +
      scale_y_continuous(limits = c(0, maxFungiObserved+10)) +
      geom_smooth(method='lm') +
      ggtitle(plotTitlePerCT_Library_Fungi) +
      theme(plot.title = element_text(hjust = 0.5)) +
      # scale_x_log10() +
      # rotate_x_text(30) +
      scale_fill_igv() +
      stat_cor(method = tolower(statCorMethod), cor.coef.name = corCoefName) -> libraryFungiCorrPlotPerCT
    print(libraryFungiCorrPlotPerCT)
    
    ## Bacteria vs. library size
    fungi_and_bacteria_alpha_with_meta %>%
      ggplot(aes(library_size_log10,bacteria_Observed)) +
      geom_point(aes(colour = library_size_log10)) +
      scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
      xlab("Original microbial library Size (log10)") + ylab("Bacterial richness") + theme_pubr() +
      # facet_wrap(facets = vars(investigation), scales = "free", nrow = 3) +
      facet_wrap(facets = vars(investigation), nrow = 3) +
      scale_y_continuous(limits = c(0, maxFungiObserved+10)) +
      geom_smooth(method='lm') +
      ggtitle(plotTitlePerCT_Library_Bacteria) +
      theme(plot.title = element_text(hjust = 0.5)) +
      # scale_x_log10() +
      # rotate_x_text(30) +
      scale_fill_igv() +
      stat_cor(method = tolower(statCorMethod), cor.coef.name = corCoefName) -> libraryBacteriaCorrPlotPerCT
    print(libraryBacteriaCorrPlotPerCT)
    
    ## Aggregated fungi vs. library size
    fungi_and_bacteria_alpha_with_meta %>%
      ggplot(aes(library_size_log10,fungi_Observed)) +
      geom_point(aes(colour = library_size_log10)) +
      scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
      xlab("Original microbial library Size (log10)") + ylab("Fungal richness") + theme_pubr() +
      scale_y_continuous(limits = c(0, NA)) +
      # facet_wrap(facets = vars(investigation), scales = "free", nrow = 3) +
      # facet_wrap(facets = vars(investigation), nrow = 3) +
      geom_smooth(method='lm') +
      ggtitle(paste(plotTitleAggregated_Library_Fungi)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      # scale_x_log10() +
      # rotate_x_text(30) +
      scale_fill_igv() +
      stat_cor(method = tolower(statCorMethod), cor.coef.name = corCoefName) -> libraryFungiCorrPlotAggregated
    print(libraryFungiCorrPlotAggregated)
    
    ## Aggregated bacteria vs. library size
    fungi_and_bacteria_alpha_with_meta %>%
      ggplot(aes(library_size_log10,bacteria_Observed)) +
      geom_point(aes(colour = library_size_log10)) +
      scale_colour_gradient(low = "yellow", high = "red", na.value = NA) +
      xlab("Original microbial library Size (log10)") + ylab("Bacterial richness") + theme_pubr() +
      scale_y_continuous(limits = c(0, NA)) +
      # facet_wrap(facets = vars(investigation), scales = "free", nrow = 3) +
      # facet_wrap(facets = vars(investigation), nrow = 3) +
      geom_smooth(method='lm') +
      ggtitle(plotTitleAggregated_Library_Bacteria) +
      theme(plot.title = element_text(hjust = 0.5)) +
      # scale_x_log10() +
      # rotate_x_text(30) +
      scale_fill_igv() +
      stat_cor(method = tolower(statCorMethod), cor.coef.name = corCoefName) -> libraryBacteriaCorrPlotAggregated
    print(libraryBacteriaCorrPlotAggregated)
    
    res <- list(libraryFungiCorrPlotPerCT=libraryFungiCorrPlotPerCT,
                libraryBacteriaCorrPlotPerCT=libraryBacteriaCorrPlotPerCT,
                libraryFungiCorrPlotAggregated=libraryFungiCorrPlotAggregated,
                libraryBacteriaCorrPlotAggregated=libraryBacteriaCorrPlotAggregated,
                fungi_and_bacteria_alpha_with_meta=fungi_and_bacteria_alpha_with_meta,
                richnessCorrPlotPerCT=richnessCorrPlotPerCT,
                richnessCorrPlotAggregated=richnessCorrPlotAggregated,
                meta2Keep=meta2Keep)
    
  }
  
  return(res)
  
}

mlCristiano <- function(metaData=metaCristianoTxNaiveFilt,
                       countData,
                       col2Predict = "cancer_status",
                       dzOfInterest = "Cancer vs. Healthy",
                       stageNum = "I",
                       modelType = "gbm",
                       dataString = "fungi_only",
                       # trainTestFlag = TRUE,
                       cutPoint = 0.5, 
                       numKFold = 10,
                       numResampleIter = 1,
                       varImpFlag = FALSE){
  require(caret) # for model building
  require(gbm) # for machine learning
  require(PRROC) # for precision-recall curves
  require(MLmetrics) # for multiclass ML
  require(e1071)
  require(ggpubr)
  
  if(col2Predict == "phenotype"){
    dzOfInterestFixed <- gsub('([[:punct:]])|\\s+','',dzOfInterest)
    metaDataFixed <- droplevels(metaData[metaData$cancer_status == "Cancer",])
    metaDataFixed$condition <- factor(ifelse(metaDataFixed[,col2Predict] == dzOfInterest,
                                             yes = dzOfInterestFixed, no = "Other"),
                                      levels = c(dzOfInterestFixed,"Other"))
    positiveClass <- dzOfInterestFixed
    negativeClass <- "Other"
  } else if(col2Predict == "cancer_status"){
    metaDataFixed <- metaData
    metaDataFixed$condition <- metaDataFixed$cancer_status
    positiveClass <- "Cancer"
    negativeClass <- "Healthy"
  } else if(col2Predict == "one_cancer_vs_healthy"){
    dzOfInterestFixed <- gsub('([[:punct:]])|\\s+','',dzOfInterest)
    metaDataFixed <- droplevels(metaData[metaData$phenotype %in% c("Healthy",dzOfInterest),])
    metaDataFixed$condition <- factor(ifelse(metaDataFixed[,"phenotype"] == dzOfInterest,
                                             yes = dzOfInterestFixed, no = "Healthy"),
                                      levels = c(dzOfInterestFixed,"Healthy"))
    positiveClass <- dzOfInterestFixed
    negativeClass <- "Healthy"
  } else if(col2Predict == "stage"){
    metaDataFixed <- metaData %>% filter(stage %in% c(stageNum, "Healthy")) %>% droplevels()
    metaDataFixed$condition <- factor(metaDataFixed$stage, levels = c(stageNum,"Healthy"))
    # print(summary(metaDataFixed$condition))
    positiveClass <- stageNum
    negativeClass <- "Healthy"
  }
  
  mlDataY <- metaDataFixed
  mlDataX <- countData[rownames(mlDataY),]
  seedNum <- 42
  
  set.seed(seedNum)
  trainX <- mlDataX
  trainY <- mlDataY[,"condition"]
  refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
  
  set.seed(seedNum)
  ctrl <- trainControl(method = "repeatedcv",
                       number = numKFold,
                       repeats = numResampleIter,
                       sampling = "up",
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,
                       verboseIter = TRUE,
                       savePredictions = TRUE,
                       allowParallel=TRUE)
  
  set.seed(seedNum)
  mlModel <- train(x = trainX,
                   y = refactoredTrainY,
                   method = modelType,
                   preProcess = c("zv"),
                   trControl = ctrl,
                   tuneGrid=data.frame(n.trees=150, interaction.depth=3,
                                       shrinkage=0.1,
                                       n.minobsinnode=1),
                   verbose = TRUE,
                   metric = "ROC")
  # print(mlModel)
  
  resPred <- mlModel$pred
  
  ## Calculate performance on concatenated fold predictions
  predProbs <- resPred
  multiClass <- resPred
  # print(head(multiClass))
  multiClass$pred <- relevel(multiClass$pred, positiveClass)
  multiClass$obs <- relevel(multiClass$obs, positiveClass)
  fg <- predProbs[resPred$obs == positiveClass,positiveClass]
  bg <- predProbs[resPred$obs == negativeClass,positiveClass]
  
  confusionMatrix <- confusionMatrix(multiClass$pred, multiClass$obs)
  prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
  
  ## Estimate AUROC CIs using cvAUC on concatenated fold predictions
  require(cvAUC)
  out95 <- ci.cvAUC(predictions = resPred[,positiveClass], labels = resPred$obs, label.ordering = c(negativeClass,positiveClass), folds = resPred$Resample, confidence = 0.95)
  out99 <- ci.cvAUC(predictions = resPred[,positiveClass], labels = resPred$obs, label.ordering = c(negativeClass,positiveClass), folds = resPred$Resample, confidence = 0.99)
  # resCvAUC <- data.frame(estimate = c(out95$cvAUC, out99$cvAUC), se = c(out95$se, out99$se), lowerCI = c(out95$ci[1], out99$ci[1]), upperCI = c(out95$ci[1], out99$ci[2]), levelCI = c(0.95,0.99))
  
  ## Split folds and calculate perf on each fold
  resPredSplit <- split(resPred, resPred$Resample)
  repX_perf <- list()
  for(zz in seq_along(resPredSplit)){
    resPredSingleRep <- resPredSplit[[zz]]
    predProbsSingleRep <- resPredSingleRep
    multiClassSingleRep <- resPredSingleRep
    multiClassSingleRep$pred <- relevel(multiClassSingleRep$pred, positiveClass)
    multiClassSingleRep$obs <- relevel(multiClassSingleRep$obs, positiveClass)
    fgSingleRep <- predProbsSingleRep[resPredSingleRep$obs == positiveClass,positiveClass]
    bgSingleRep <- predProbsSingleRep[resPredSingleRep$obs == negativeClass,positiveClass]
    
    rep_roc <- roc.curve(scores.class0 = fgSingleRep, scores.class1 = bgSingleRep, curve = T)
    rep_pr <- pr.curve(scores.class0 = fgSingleRep, scores.class1 = bgSingleRep, curve = T, rand.compute=T)
    
    repX_perf[[zz]] <- data.frame(auroc=rep_roc$auc,
                                  aupr=rep_pr$auc.integral,
                                  rep=paste0("Fold",zz), 
                                  dataString=dataString,
                                  minorityClassSize=sum(multiClass$obs==positiveClass),
                                  majorityClassSize=sum(multiClass$obs==negativeClass),
                                  diseaseType = dzOfInterest,
                                  col2Predict = col2Predict)
  }
  
  # SUMMARIZE MODEL PERFORMANCES
  rep_perf <- do.call(rbind, repX_perf)
  perf <- data.frame(auroc = prroc_roc$auc,
                     aupr = prroc_pr$auc.integral,
                     aucEstimate = out95$cvAUC, # either out95 or out99 work (same result)
                     aucSE95 = out95$se,
                     lowerCI95 = out95$ci[1],
                     upperCI95 = out95$ci[2],
                     aucSE99 = out99$se,
                     lowerCI99 = out99$ci[1],
                     upperCI99 = out99$ci[2],
                     diseaseType = dzOfInterest,
                     dataString=dataString,
                     col2Predict = col2Predict)
  
  # predProbs <- mlModel$pred
  # multiClass <- mlModel$pred
  # multiClass$pred <- relevel(multiClass$pred, positiveClass)
  # multiClass$obs <- relevel(multiClass$obs, positiveClass)
  # fg <- predProbs[mlModel$pred$obs == positiveClass,positiveClass]
  # bg <- predProbs[mlModel$pred$obs == negativeClass,positiveClass]
  # 
  # prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  # prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
  # 
  print(confusionMatrix(multiClass$pred, multiClass$obs))
  print(multiClassSummary(multiClass, lev = levels(multiClass$obs)))
  
  varImpBestModelDF <- as.data.frame(varImp(mlModel$finalModel, scale = FALSE))
  varImpBestModelDF2 <- rownames_to_column(varImpBestModelDF, "Taxa")
  varImpBestModelDF2Ordered <- varImpBestModelDF2[order(-varImpBestModelDF2$Overall),]
  colnames(varImpBestModelDF2Ordered)[2] <- "varImp"
  varImpBestModelDF2OrderedNonzero <- varImpBestModelDF2Ordered[varImpBestModelDF2Ordered$varImp != 0,]
  
  if(varImpFlag){
    print(head(varImpBestModelDF2OrderedNonzero))
  }
  
  plot(prroc_pr)
  plot(prroc_roc)
  
  res <- list(mlModel=mlModel,
              perf=perf,
              rep_perf=rep_perf,
              multiClass=multiClass,
              varImp=varImpBestModelDF2OrderedNonzero)
  return(res)
  
}

ml1VsAllCristiano10kRep1_Iterate <- function(metaData=metaCristianoTxNaiveFilt,
                                     listOfData = list(cristiano_rep200Data_Filt,
                                                        cristianoRep200BacteriaSpecies,
                                                        cristianoRep200FungiSpecies,
                                                        cristiano_rep200Data_Filt_Fungi_Decontam,
                                                        cristianoRep200FungiSpeciesShared,
                                                        cristianoRep200BacteriaSpeciesShared,
                                                        cbind(cristianoRep200FungiSpeciesShared, 
                                                          cristianoRep200BacteriaSpeciesShared),
                                                        topXHopkinsData),
                                     dataStringList = c("full_rep200",
                                                        "bacteria_only",
                                                        "fungi_only",
                                                        "fungi_decontam",
                                                        "fungi_intersected_with_Weizmann",
                                                        "bacteria_intersected_with_Weizmann",
                                                        "fungi_and_bacteria_intersected_with_Weizmann",
                                                        "topX_fungi"),
                                     col2Predict = "one_cancer_vs_healthy",
                                     dzOfInterestList = list("Bile Duct Cancer",
                                                              "Breast Cancer",
                                                              "Colorectal Cancer",
                                                              "Gastric cancer",
                                                              "Lung Cancer",
                                                              "Ovarian Cancer",
                                                              "Pancreatic Cancer"),
                                     stageNum = "stageI",
                                     scrambleFlag=FALSE,
                                     shuffleFlag=FALSE,
                                     modelType = "gbm",
                                     cutPoint = 0.5, 
                                     numKFold = 10,
                                     numResampleIter = 1,
                                     varImpFlag = FALSE){
  require(caret) # for model building
  require(gbm) # for machine learning
  require(PRROC) # for precision-recall curves
  require(MLmetrics) # for multiclass ML
  require(e1071)
  require(ggpubr)

  rep_perfTmp <- list()
  perfTmp <- list()
  rep_perfTmp2 <- list()
  perfTmp2 <- list()

  for(jj in seq_along(listOfData)){

    snmData <- listOfData[[jj]]
    dataString <- dataStringList[[jj]]
    print(sprintf("Working on dataset: %s",dataString))

    for(ii in seq_along(dzOfInterestList)){

    dzOfInterest <- dzOfInterestList[[ii]]
    print(sprintf("Working on cancer type: %s",dzOfInterest))

    if(col2Predict == "phenotype"){
    dzOfInterestFixed <- gsub('([[:punct:]])|\\s+','',dzOfInterest)
    metaDataFixed <- droplevels(metaData[metaData$cancer_status == "Cancer",])
    metaDataFixed$condition <- factor(ifelse(metaDataFixed[,col2Predict] == dzOfInterest,
                                             yes = dzOfInterestFixed, no = "Other"),
                                      levels = c(dzOfInterestFixed,"Other"))
    positiveClass <- dzOfInterestFixed
    negativeClass <- "Other"
  } else if(col2Predict == "cancer_status"){
    metaDataFixed <- metaData
    metaDataFixed$condition <- metaDataFixed$cancer_status
    positiveClass <- "Cancer"
    negativeClass <- "Healthy"
  } else if(col2Predict == "one_cancer_vs_healthy"){
    dzOfInterestFixed <- gsub('([[:punct:]])|\\s+','',dzOfInterest)
    metaDataFixed <- droplevels(metaData[metaData$phenotype %in% c("Healthy",dzOfInterest),])
    metaDataFixed$condition <- factor(ifelse(metaDataFixed[,"phenotype"] == dzOfInterest,
                                             yes = dzOfInterestFixed, no = "Healthy"),
                                      levels = c(dzOfInterestFixed,"Healthy"))
    positiveClass <- dzOfInterestFixed
    negativeClass <- "Healthy"
  } else if(col2Predict == "stage"){
    metaDataFixed <- metaData %>% filter(stage %in% c(stageNum, "Healthy")) %>% droplevels()
    metaDataFixed$condition <- factor(metaDataFixed$stage, levels = c(stageNum,"Healthy"))
    # print(summary(metaDataFixed$condition))
    positiveClass <- stageNum
    negativeClass <- "Healthy"
  }
  
  mlDataY <- metaDataFixed
  mlDataX <- snmData[rownames(mlDataY),]

  if(scrambleFlag){
    set.seed(42)
    mlDataY[,"condition"] <- sample(mlDataY[,"condition"])
  } else if(shuffleFlag){
    set.seed(42)
    mlDataX_shuffled <- mlDataX[sample(nrow(mlDataX)),]
    rownames(mlDataX_shuffled) <- rownames(mlDataX)
    mlDataX <- mlDataX_shuffled
  }
  
  set.seed(42)
  trainX <- mlDataX
  trainY <- mlDataY[,"condition"]
  refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
  
  set.seed(42)
  ctrl <- trainControl(method = "repeatedcv",
                       number = numKFold,
                       repeats = numResampleIter,
                       sampling = "up",
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,
                       verboseIter = TRUE,
                       savePredictions = TRUE,
                       allowParallel=TRUE)
  
  set.seed(42)
  mlModel <- train(x = trainX,
                   y = refactoredTrainY,
                   method = modelType,
                   preProcess = c("zv"),
                   trControl = ctrl,
                   tuneGrid=data.frame(n.trees=150, interaction.depth=3,
                                       shrinkage=0.1,
                                       n.minobsinnode=1),
                   verbose = TRUE,
                   metric = "ROC")
  # print(mlModel)
  
  resPred <- mlModel$pred
  
  ## Calculate performance on concatenated fold predictions
  predProbs <- resPred
  multiClass <- resPred
  # print(head(multiClass))
  multiClass$pred <- relevel(multiClass$pred, positiveClass)
  multiClass$obs <- relevel(multiClass$obs, positiveClass)
  fg <- predProbs[resPred$obs == positiveClass,positiveClass]
  bg <- predProbs[resPred$obs == negativeClass,positiveClass]
  
  confusionMatrix <- confusionMatrix(multiClass$pred, multiClass$obs)
  prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
  
  ## Estimate AUROC CIs using cvAUC on concatenated fold predictions
  require(cvAUC)
  out95 <- ci.cvAUC(predictions = resPred[,positiveClass], labels = resPred$obs, label.ordering = c(negativeClass,positiveClass), folds = resPred$Resample, confidence = 0.95)
  out99 <- ci.cvAUC(predictions = resPred[,positiveClass], labels = resPred$obs, label.ordering = c(negativeClass,positiveClass), folds = resPred$Resample, confidence = 0.99)
  # resCvAUC <- data.frame(estimate = c(out95$cvAUC, out99$cvAUC), se = c(out95$se, out99$se), lowerCI = c(out95$ci[1], out99$ci[1]), upperCI = c(out95$ci[1], out99$ci[2]), levelCI = c(0.95,0.99))
  
  ## Split folds and calculate perf on each fold
  resPredSplit <- split(resPred, resPred$Resample)
  repX_perf <- list()
  for(zz in seq_along(resPredSplit)){
    resPredSingleRep <- resPredSplit[[zz]]
    predProbsSingleRep <- resPredSingleRep
    multiClassSingleRep <- resPredSingleRep
    multiClassSingleRep$pred <- relevel(multiClassSingleRep$pred, positiveClass)
    multiClassSingleRep$obs <- relevel(multiClassSingleRep$obs, positiveClass)
    fgSingleRep <- predProbsSingleRep[resPredSingleRep$obs == positiveClass,positiveClass]
    bgSingleRep <- predProbsSingleRep[resPredSingleRep$obs == negativeClass,positiveClass]
    
    rep_roc <- roc.curve(scores.class0 = fgSingleRep, scores.class1 = bgSingleRep, curve = T)
    rep_pr <- pr.curve(scores.class0 = fgSingleRep, scores.class1 = bgSingleRep, curve = T, rand.compute=T)
    
    repX_perf[[zz]] <- data.frame(auroc=rep_roc$auc,
                                  aupr=rep_pr$auc.integral,
                                  rep=paste0("Fold",zz), 
                                  dataString=dataString,
                                  minorityClassSize=sum(multiClass$obs==positiveClass),
                                  majorityClassSize=sum(multiClass$obs==negativeClass),
                                  diseaseType = dzOfInterest,
                                  col2Predict = col2Predict)
  }
  
  # SUMMARIZE MODEL PERFORMANCES
  rep_perfTmp[[ii]] <- do.call(rbind, repX_perf)
  perfTmp[[ii]] <- data.frame(auroc = prroc_roc$auc,
                     aupr = prroc_pr$auc.integral,
                     aucEstimate = out95$cvAUC, # either out95 or out99 work (same result)
                     aucSE95 = out95$se,
                     lowerCI95 = out95$ci[1],
                     upperCI95 = out95$ci[2],
                     aucSE99 = out99$se,
                     lowerCI99 = out99$ci[1],
                     upperCI99 = out99$ci[2],
                     diseaseType = dzOfInterest,
                     dataString=dataString,
                     col2Predict = col2Predict)

  print(confusionMatrix(multiClass$pred, multiClass$obs))
  print(multiClassSummary(multiClass, lev = levels(multiClass$obs)))
  
  varImpBestModelDF <- as.data.frame(varImp(mlModel$finalModel, scale = FALSE))
  varImpBestModelDF2 <- rownames_to_column(varImpBestModelDF, "Taxa")
  varImpBestModelDF2Ordered <- varImpBestModelDF2[order(-varImpBestModelDF2$Overall),]
  colnames(varImpBestModelDF2Ordered)[2] <- "varImp"
  varImpBestModelDF2OrderedNonzero <- varImpBestModelDF2Ordered[varImpBestModelDF2Ordered$varImp != 0,]
  
  if(varImpFlag){print(head(varImpBestModelDF2OrderedNonzero))}
  
  plot(prroc_pr)
  plot(prroc_roc)
  
  }

  rep_perfTmp2[[jj]] <- do.call(rbind, rep_perfTmp)
  perfTmp2[[jj]] <- do.call(rbind, perfTmp)

  }

  rep_perf <- do.call(rbind, rep_perfTmp2)
  perf <- do.call(rbind, perfTmp2)

  res <- list(mlModel=mlModel,
              perf=perf,
              rep_perf=rep_perf,
              multiClass=multiClass,
              varImp=varImpBestModelDF2OrderedNonzero)
  return(res)
  
}

ml1VsAllCristiano10kRep1_Iterate_Stage <- function(metaData=metaCristianoTxNaiveFilt,
                                     listOfData = list(cristiano_rep200Data_Filt,
                                                        cristianoRep200BacteriaSpecies,
                                                        cristianoRep200FungiSpecies,
                                                        cristiano_rep200Data_Filt_Fungi_Decontam,
                                                        cristianoRep200FungiSpeciesShared,
                                                        cristianoRep200BacteriaSpeciesShared,
                                                        cbind(cristianoRep200FungiSpeciesShared, 
                                                          cristianoRep200BacteriaSpeciesShared),
                                                        topXHopkinsData),
                                     dataStringList = c("full_rep200",
                                                        "bacteria_only",
                                                        "fungi_only",
                                                        "fungi_decontam",
                                                        "fungi_intersected_with_Weizmann",
                                                        "bacteria_intersected_with_Weizmann",
                                                        "fungi_and_bacteria_intersected_with_Weizmann",
                                                        "topX_fungi"),
                                     col2Predict = "stage",
                                     dzOfInterestList = list("Bile Duct Cancer",
                                                              "Breast Cancer",
                                                              "Colorectal Cancer",
                                                              "Gastric cancer",
                                                              "Lung Cancer",
                                                              "Ovarian Cancer",
                                                              "Pancreatic Cancer"),
                                     stageNumList = list("I",
                                                        "II", 
                                                        "III", 
                                                        "IV"),
                                     modelType = "gbm",
                                     cutPoint = 0.5, 
                                     numKFold = 10,
                                     numResampleIter = 1,
                                     varImpFlag = FALSE){
  require(caret) # for model building
  require(gbm) # for machine learning
  require(PRROC) # for precision-recall curves
  require(MLmetrics) # for multiclass ML
  require(e1071)
  require(ggpubr)

  rep_perfTmp <- list()
  perfTmp <- list()
  rep_perfTmp2 <- list()
  perfTmp2 <- list()

  for(jj in seq_along(listOfData)){

    countData <- listOfData[[jj]]
    dataString <- dataStringList[[jj]]
    print(sprintf("Working on dataset: %s",dataString))

    for(ii in seq_along(stageNumList)){

    dzOfInterest <- dzOfInterestList[[ii]]
    stageNum <- stageNumList[[ii]]
    print(sprintf("Working on stage: %s",stageNum))

    if(col2Predict == "phenotype"){
    dzOfInterestFixed <- gsub('([[:punct:]])|\\s+','',dzOfInterest)
    metaDataFixed <- droplevels(metaData[metaData$cancer_status == "Cancer",])
    metaDataFixed$condition <- factor(ifelse(metaDataFixed[,col2Predict] == dzOfInterest,
                                             yes = dzOfInterestFixed, no = "Other"),
                                      levels = c(dzOfInterestFixed,"Other"))
    positiveClass <- dzOfInterestFixed
    negativeClass <- "Other"
  } else if(col2Predict == "cancer_status"){
    metaDataFixed <- metaData
    metaDataFixed$condition <- metaDataFixed$cancer_status
    positiveClass <- "Cancer"
    negativeClass <- "Healthy"
  } else if(col2Predict == "one_cancer_vs_healthy"){
    dzOfInterestFixed <- gsub('([[:punct:]])|\\s+','',dzOfInterest)
    metaDataFixed <- droplevels(metaData[metaData$phenotype %in% c("Healthy",dzOfInterest),])
    metaDataFixed$condition <- factor(ifelse(metaDataFixed[,"phenotype"] == dzOfInterest,
                                             yes = dzOfInterestFixed, no = "Healthy"),
                                      levels = c(dzOfInterestFixed,"Healthy"))
    positiveClass <- dzOfInterestFixed
    negativeClass <- "Healthy"
  } else if(col2Predict == "stage"){
    metaDataFixed <- metaData %>% filter(stage %in% c(stageNum, "Healthy")) %>% droplevels()
    metaDataFixed$condition <- factor(metaDataFixed$stage, levels = c(stageNum,"Healthy"))
    # print(summary(metaDataFixed$condition))
    positiveClass <- stageNum
    negativeClass <- "Healthy"
  }
  
  mlDataY <- metaDataFixed
  mlDataX <- countData[rownames(mlDataY),]
  
  set.seed(42)
  trainX <- mlDataX
  trainY <- mlDataY[,"condition"]
  refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
  
  set.seed(42)
  ctrl <- trainControl(method = "repeatedcv",
                       number = numKFold,
                       repeats = numResampleIter,
                       sampling = "up",
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,
                       verboseIter = TRUE,
                       savePredictions = TRUE,
                       allowParallel=TRUE)
  
  set.seed(42)
  mlModel <- train(x = trainX,
                   y = refactoredTrainY,
                   method = modelType,
                   preProcess = c("zv"),
                   trControl = ctrl,
                   tuneGrid=data.frame(n.trees=150, interaction.depth=3,
                                       shrinkage=0.1,
                                       n.minobsinnode=1),
                   verbose = TRUE,
                   metric = "ROC")
  # print(mlModel)
  
  resPred <- mlModel$pred
  
  ## Calculate performance on concatenated fold predictions
  predProbs <- resPred
  multiClass <- resPred
  # print(head(multiClass))
  multiClass$pred <- relevel(multiClass$pred, positiveClass)
  multiClass$obs <- relevel(multiClass$obs, positiveClass)
  fg <- predProbs[resPred$obs == positiveClass,positiveClass]
  bg <- predProbs[resPred$obs == negativeClass,positiveClass]
  
  confusionMatrix <- confusionMatrix(multiClass$pred, multiClass$obs)
  prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
  
  ## Estimate AUROC CIs using cvAUC on concatenated fold predictions
  require(cvAUC)
  out95 <- ci.cvAUC(predictions = resPred[,positiveClass], labels = resPred$obs, label.ordering = c(negativeClass,positiveClass), folds = resPred$Resample, confidence = 0.95)
  out99 <- ci.cvAUC(predictions = resPred[,positiveClass], labels = resPred$obs, label.ordering = c(negativeClass,positiveClass), folds = resPred$Resample, confidence = 0.99)
  # resCvAUC <- data.frame(estimate = c(out95$cvAUC, out99$cvAUC), se = c(out95$se, out99$se), lowerCI = c(out95$ci[1], out99$ci[1]), upperCI = c(out95$ci[1], out99$ci[2]), levelCI = c(0.95,0.99))
  
  ## Split folds and calculate perf on each fold
  resPredSplit <- split(resPred, resPred$Resample)
  repX_perf <- list()
  for(zz in seq_along(resPredSplit)){
    resPredSingleRep <- resPredSplit[[zz]]
    predProbsSingleRep <- resPredSingleRep
    multiClassSingleRep <- resPredSingleRep
    multiClassSingleRep$pred <- relevel(multiClassSingleRep$pred, positiveClass)
    multiClassSingleRep$obs <- relevel(multiClassSingleRep$obs, positiveClass)
    fgSingleRep <- predProbsSingleRep[resPredSingleRep$obs == positiveClass,positiveClass]
    bgSingleRep <- predProbsSingleRep[resPredSingleRep$obs == negativeClass,positiveClass]
    
    rep_roc <- roc.curve(scores.class0 = fgSingleRep, scores.class1 = bgSingleRep, curve = T)
    rep_pr <- pr.curve(scores.class0 = fgSingleRep, scores.class1 = bgSingleRep, curve = T, rand.compute=T)
    
    repX_perf[[zz]] <- data.frame(auroc=rep_roc$auc,
                                  aupr=rep_pr$auc.integral,
                                  rep=paste0("Fold",zz), 
                                  dataString=dataString,
                                  stageNum=stageNum,
                                  minorityClassSize=sum(multiClass$obs==positiveClass),
                                  majorityClassSize=sum(multiClass$obs==negativeClass),
                                  diseaseType = dzOfInterest,
                                  col2Predict = col2Predict)
  }
  
  # SUMMARIZE MODEL PERFORMANCES
  rep_perfTmp[[ii]] <- do.call(rbind, repX_perf)
  perfTmp[[ii]] <- data.frame(auroc = prroc_roc$auc,
                     aupr = prroc_pr$auc.integral,
                     aucEstimate = out95$cvAUC, # either out95 or out99 work (same result)
                     aucSE95 = out95$se,
                     lowerCI95 = out95$ci[1],
                     upperCI95 = out95$ci[2],
                     aucSE99 = out99$se,
                     lowerCI99 = out99$ci[1],
                     upperCI99 = out99$ci[2],
                     diseaseType = dzOfInterest,
                     dataString=dataString,
                     stageNum=stageNum,
                     col2Predict = col2Predict)

  print(confusionMatrix(multiClass$pred, multiClass$obs))
  print(multiClassSummary(multiClass, lev = levels(multiClass$obs)))
  
  varImpBestModelDF <- as.data.frame(varImp(mlModel$finalModel, scale = FALSE))
  varImpBestModelDF2 <- rownames_to_column(varImpBestModelDF, "Taxa")
  varImpBestModelDF2Ordered <- varImpBestModelDF2[order(-varImpBestModelDF2$Overall),]
  colnames(varImpBestModelDF2Ordered)[2] <- "varImp"
  varImpBestModelDF2OrderedNonzero <- varImpBestModelDF2Ordered[varImpBestModelDF2Ordered$varImp != 0,]
  
  if(varImpFlag){print(head(varImpBestModelDF2OrderedNonzero))}
  
  plot(prroc_pr)
  plot(prroc_roc)
  
  }

  rep_perfTmp2[[jj]] <- do.call(rbind, rep_perfTmp)
  perfTmp2[[jj]] <- do.call(rbind, perfTmp)

  }

  rep_perf <- do.call(rbind, rep_perfTmp2)
  perf <- do.call(rbind, perfTmp2)

  res <- list(mlModel=mlModel,
              perf=perf,
              rep_perf=rep_perf,
              multiClass=multiClass,
              varImp=varImpBestModelDF2OrderedNonzero)
  return(res)
  
}

plotMLWithCIs <- function(modelOutput, showRepCurves=FALSE, sizeRepCurves=0.5,
                          sizeAvgCurve=1, colorAvgCurve="blue",ciFillColor="orange",
                          colorRepCurves="lightgray",
                          interpXLength = 1e3, ciAlpha = 0.3, ciLevel = 0.95,
                          positiveClass="Cancer", negativeClass="Healthy"){
  require(gmodels)
  # TO DO: INPUT LIST OF MODELS, ITERATE OVER LIST, AND PLOT TOGETHER
  interpX <- seq(0, 1, length.out = interpXLength)
  
  tmpML_Res <- modelOutput[["multiClass"]]
  tmpML_Res$Rep <- gsub("^Fold[0-9]+\\.","",tmpML_Res$Resample)
  numRep <- length(names(table(tmpML_Res$Rep)))
  tmpML_Res_Split <- split(tmpML_Res, tmpML_Res$Rep)
  
  interpROCY <- list()
  interpPRY <- list()
  rocCurveData <- list()
  auroc <- vector()
  prCurveData <- list()
  aupr <- vector()
  for(ii in seq_along(tmpML_Res_Split)){
    mlPreds <- tmpML_Res_Split[[ii]]
    fg <- mlPreds[mlPreds$obs == positiveClass, positiveClass]
    bg <- mlPreds[mlPreds$obs == negativeClass, positiveClass]
    prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
    rocCurveData[[ii]] <- data.frame(prroc_roc$curve)
    colnames(rocCurveData[[ii]]) <- c("fpr","tpr","thresh")
    auroc[ii] <- prroc_roc$auc
    interpROCVals <- approx(x = rocCurveData[[ii]][,"fpr"], y = rocCurveData[[ii]][,"tpr"], xout = interpX)
    interpROCY[[ii]] <- interpROCVals$y
    
    prCurveData[[ii]] <- data.frame(prroc_pr$curve)
    aupr[ii] <- prroc_pr$auc.integral
    colnames(prCurveData[[ii]]) <- c("recall","precision","thresh")
    interpPRVals <- approx(x = prCurveData[[ii]][,"recall"], y = prCurveData[[ii]][,"precision"], xout = interpX)
    interpPRY[[ii]] <- interpPRVals$y
  }
  # ROC
  aurocCI <- ci(as.vector(auroc))
  interpROCYDf <- data.frame(interpROCY)
  colnames(interpROCYDf) <- paste0("Rep",1:numRep)
  interpROCYDf_CI <- data.frame(t(apply(as.matrix(interpROCYDf), 1, function(x) ci(x, confidence = ciLevel))))
  interpROCYDf_CI$xval <- interpX
  interpROCYDf_CI <- rbind(rep(0,dim(interpROCYDf_CI)[2]), # make sure plot begins at 0,0
                           interpROCYDf_CI,
                           rep(1,dim(interpROCYDf_CI)[2])) # make sure plot ends at 1,1
  # PR
  auprCI <- ci(as.vector(aupr))
  interpPRYDf <- data.frame(interpPRY)
  colnames(interpPRYDf) <- paste0("Rep",1:numRep)
  interpPRYDf_CI <- data.frame(t(apply(as.matrix(interpPRYDf), 1, function(x) ci(x, confidence = ciLevel))))
  interpPRYDf_CI$xval <- interpX
  interpPRYDf_CI <- rbind(c(1,0,0,0,0), # make sure plot begins at 0,1
                          interpPRYDf_CI,
                          c(0,0,0,0,1)) # make sure plot ends at 1,0
  
  if(showRepCurves){
    # ROC
    rocPlot <- ggplot(data = rocCurveData[[1]], aes(x = fpr, y = tpr)) + 
      geom_path(color = colorRepCurves, size = sizeRepCurves) + theme_minimal() +
      geom_abline(intercept = 0, slope = 1, color = "gray", size = 0.5) +
      coord_equal(ratio=1) + xlim(0, 1) + ylim(0,1) +
      labs(x = "False Positive Rate", y = "True Positive Rate")
    for(jj in 2:numRep){
      rocPlot <- rocPlot + geom_path(data = rocCurveData[[jj]], aes(x = fpr, y = tpr), color = colorRepCurves, size = sizeRepCurves)
    }
    rocPlot <- rocPlot + geom_path(data = interpROCYDf_CI, aes(x = xval, y = Estimate), color = colorAvgCurve, size = sizeAvgCurve) +
      geom_ribbon(data = interpROCYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
                  fill = ciFillColor, alpha = ciAlpha, inherit.aes = F)
    # PR
    prPlot <- ggplot(data = prCurveData[[1]], aes(x = recall, y = precision)) + 
      geom_path(color = colorRepCurves, size = sizeRepCurves) + theme_minimal() +
      coord_equal(ratio=1) + xlim(0, 1) + ylim(0,1) +
      labs(x = "Recall", y = "Precision")
    for(jj in 2:numRep){
      prPlot <- prPlot + geom_path(data = prCurveData[[jj]], aes(x = recall, y = precision), color = colorRepCurves, size = sizeRepCurves)
    }
    prPlot <- prPlot + geom_path(data = interpPRYDf_CI, aes(x = xval, y = Estimate), color = colorAvgCurve, size = sizeAvgCurve) +
      geom_ribbon(data = interpPRYDf_CI, aes(x = xval, ymin = CI.lower, ymax = CI.upper), 
                  fill = ciFillColor, alpha = ciAlpha, inherit.aes = F)
  } else{
    # ROC
    rocPlot <- ggplot(data = interpROCYDf_CI, aes(x = xval, y = Estimate)) + 
      geom_path(color = colorAvgCurve, size = sizeAvgCurve) + theme_minimal() +
      geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), fill = ciFillColor, alpha = 0.3) +
      geom_abline(intercept = 0, slope = 1, color = "gray", size = 0.5) +
      coord_equal(ratio=1) + xlim(0, 1) + ylim(0,1) +
      labs(x = "False Positive Rate", y = "True Positive Rate")
    # PR
    prPlot <- ggplot(data = interpPRYDf_CI, aes(x = xval, y = Estimate)) + 
      geom_path(color = colorAvgCurve, size = sizeAvgCurve) + theme_minimal() +
      geom_ribbon(aes(ymin = CI.lower, ymax = CI.upper), fill = ciFillColor, alpha = 0.3) +
      coord_equal(ratio=1) + xlim(0, 1) + ylim(0,1) +
      labs(x = "Recall", y = "Precision")
  }
  combinedPlots <- ggarrange(rocPlot, prPlot, ncol = 2)
  print(combinedPlots)
  res <- list(
    rocPlot=rocPlot,
    prPlot=prPlot,
    combinedPlots=combinedPlots,
    rocCurveData=rocCurveData,
    prCurveData=prCurveData,
    interpROCYDf_CI=interpROCYDf_CI,
    interpPRYDf_CI=interpPRYDf_CI,
    auroc=auroc,
    aurocCI=aurocCI,
    aupr=aupr,
    auprCI=auprCI
  )
  return(res)
}

vsnmFunctionUCSD <- function(qcData, qcMetadata){
  require(limma)
  require(edgeR)
  require(dplyr)
  require(snm)
  # Set up design matrix
  covDesignNorm <- model.matrix(~0 + disease_type_consol +
                                  host_age +
                                  sex,
                                data = qcMetadata)
  # Check row dimensions
  print(dim(covDesignNorm))
  print(dim(qcData))
  print(dim(covDesignNorm)[1] == dim(qcData)[1])
  colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
  
  # Set up counts matrix
  counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)
  
  # Normalize using edgeR and then plug into voom
  dge <- DGEList(counts = counts)
  vdge_data <- voom(dge, design = covDesignNorm, plot = TRUE, save.plot = TRUE,
                    normalize.method="quantile")
  vdge_dataE <- t(vdge_data$E)
  
  # Apply
  bio.var <- model.matrix(~disease_type_consol,
                          data=qcMetadata)
  adj.var <- model.matrix(~host_age +
                            sex,
                          data=qcMetadata)
  colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
  colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
  
  snmDataObjOnly <- snm(raw.dat = vdge_data$E, 
                        bio.var = bio.var, 
                        adj.var = adj.var, 
                        rm.adj=TRUE,
                        verbose = TRUE,
                        diagnose = TRUE)
  snmData <- t(snmDataObjOnly$norm.dat)
  return(snmData)
}

vsnmFunctionUCSD_IOR <- function(qcData, qcMetadata){
  require(limma)
  require(edgeR)
  require(dplyr)
  require(snm)
  # Set up design matrix
  covDesignNorm <- model.matrix(~0 + disease_type_consol +
                                  Responder +
                                  host_age +
                                  sex,
                                data = qcMetadata)
  # Check row dimensions
  print(dim(covDesignNorm))
  print(dim(qcData))
  print(dim(covDesignNorm)[1] == dim(qcData)[1])
  colnames(covDesignNorm) <- gsub('([[:punct:]])|\\s+','',colnames(covDesignNorm))
  
  # Set up counts matrix
  counts <- t(qcData) # DGEList object from a table of counts (rows=features, columns=samples)
  
  # Normalize using edgeR and then plug into voom
  dge <- DGEList(counts = counts)
  vdge_data <- voom(dge, design = covDesignNorm, plot = TRUE, save.plot = TRUE,
                    normalize.method="quantile")
  vdge_dataE <- t(vdge_data$E)
  
  # Apply
  bio.var <- model.matrix(~disease_type_consol + Responder,
                          data=qcMetadata)
  adj.var <- model.matrix(~host_age +
                            sex,
                          data=qcMetadata)
  colnames(bio.var) <- gsub('([[:punct:]])|\\s+','',colnames(bio.var))
  colnames(adj.var) <- gsub('([[:punct:]])|\\s+','',colnames(adj.var))
  
  snmDataObjOnly <- snm(raw.dat = vdge_data$E, 
                        bio.var = bio.var, 
                        adj.var = adj.var, 
                        rm.adj=TRUE,
                        verbose = TRUE,
                        diagnose = TRUE)
  snmData <- t(snmDataObjOnly$norm.dat)
  return(snmData)
}

ml1VsAllUCSD <- function(metaData,
                              snmData,
                              col2Predict = "HvsC",
                              dzOfInterest = "NSCLC",
                              modelType = "gbm",
                              dataString = "fungi_only",
                              # trainTestFlag = TRUE,
                              cutPoint = 0.5, 
                              numKFold = 10,
                              scrambleFlag=FALSE,
                              shuffleFlag=FALSE,
                              numResampleIter = 1,
                              varImpFlag = FALSE){
  require(caret) # for model building
  require(gbm) # for machine learning
  require(PRROC) # for precision-recall curves
  require(MLmetrics) # for multiclass ML
  require(e1071)
  require(ggpubr)
  
  if(col2Predict == "HvsC"){
    metaDataFixed <- metaData
    metaDataFixed$condition <- factor(metaDataFixed$HvsC, levels = c("Cancer","Control"))
    positiveClass <- "Cancer"
    negativeClass <- "Control"
  } else if(col2Predict == "one_cancer_vs_controls"){
    dzOfInterestFixed <- gsub('([[:punct:]])|\\s+','',dzOfInterest)
    metaDataFixed <- droplevels(metaData[metaData$disease_type_consol %in% c("Control",dzOfInterest),])
    metaDataFixed$condition <- factor(ifelse(metaDataFixed[,"disease_type_consol"] == dzOfInterest,
                                             yes = dzOfInterestFixed, no = "Control"),
                                      levels = c(dzOfInterestFixed,"Control"))
    positiveClass <- dzOfInterestFixed
    negativeClass <- "Control"
  } else if(col2Predict == "one_cancer_vs_others"){
    dzOfInterestFixed <- gsub('([[:punct:]])|\\s+','',dzOfInterest)
    metaDataFixed <- droplevels(metaData[metaData$disease_type_consol != "Control",])
    metaDataFixed$condition <- factor(ifelse(metaDataFixed[,"disease_type_consol"] == dzOfInterest,
                                             yes = dzOfInterestFixed, no = "Other"),
                                      levels = c(dzOfInterestFixed,"Other"))
    positiveClass <- dzOfInterestFixed
    negativeClass <- "Other"
    
  }
  
  mlDataY <- metaDataFixed
  mlDataX <- snmData[rownames(mlDataY),]

  if(scrambleFlag){
    set.seed(42)
    mlDataY[,"condition"] <- sample(mlDataY[,"condition"])
  } else if(shuffleFlag){
    set.seed(42)
    mlDataX_shuffled <- mlDataX[sample(nrow(mlDataX)),]
    rownames(mlDataX_shuffled) <- rownames(mlDataX)
    mlDataX <- mlDataX_shuffled
  }
  
  set.seed(42)
  trainX <- mlDataX
  trainY <- mlDataY[,"condition"]
  refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
  
  set.seed(42)
  ctrl <- trainControl(method = "repeatedcv",
                       number = numKFold,
                       repeats = numResampleIter,
                       sampling = "up",
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,
                       verboseIter = TRUE,
                       savePredictions = TRUE,
                       allowParallel=TRUE)
  
  set.seed(42)
  mlModel <- train(x = trainX,
                   y = refactoredTrainY,
                   method = modelType,
                   preProcess = c("nzv"),
                   trControl = ctrl,
                   tuneGrid=data.frame(n.trees=150, interaction.depth=3,
                                       shrinkage=0.1,
                                       n.minobsinnode=1),
                   verbose = TRUE,
                   metric = "ROC")
  # print(mlModel)
  
  resPred <- mlModel$pred
  
  ## Calculate performance on concatenated fold predictions
  predProbs <- resPred
  multiClass <- resPred
  multiClass$pred <- relevel(multiClass$pred, positiveClass)
  multiClass$obs <- relevel(multiClass$obs, positiveClass)
  fg <- predProbs[resPred$obs == positiveClass,positiveClass]
  bg <- predProbs[resPred$obs == negativeClass,positiveClass]
  
  confusionMatrix <- confusionMatrix(multiClass$pred, multiClass$obs)
  prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
  
  ## Estimate AUROC CIs using cvAUC on concatenated fold predictions
  require(cvAUC)
  out95 <- ci.cvAUC(predictions = resPred[,positiveClass], labels = resPred$obs, label.ordering = c(negativeClass,positiveClass), folds = resPred$Resample, confidence = 0.95)
  out99 <- ci.cvAUC(predictions = resPred[,positiveClass], labels = resPred$obs, label.ordering = c(negativeClass,positiveClass), folds = resPred$Resample, confidence = 0.99)
  # resCvAUC <- data.frame(estimate = c(out95$cvAUC, out99$cvAUC), se = c(out95$se, out99$se), lowerCI = c(out95$ci[1], out99$ci[1]), upperCI = c(out95$ci[1], out99$ci[2]), levelCI = c(0.95,0.99))
  
  ## Split folds and calculate perf on each fold
  resPredSplit <- split(resPred, resPred$Resample)
  repX_perf <- list()
  for(zz in seq_along(resPredSplit)){
    resPredSingleRep <- resPredSplit[[zz]]
    predProbsSingleRep <- resPredSingleRep
    multiClassSingleRep <- resPredSingleRep
    multiClassSingleRep$pred <- relevel(multiClassSingleRep$pred, positiveClass)
    multiClassSingleRep$obs <- relevel(multiClassSingleRep$obs, positiveClass)
    fgSingleRep <- predProbsSingleRep[resPredSingleRep$obs == positiveClass,positiveClass]
    bgSingleRep <- predProbsSingleRep[resPredSingleRep$obs == negativeClass,positiveClass]
    
    rep_roc <- roc.curve(scores.class0 = fgSingleRep, scores.class1 = bgSingleRep, curve = T)
    rep_pr <- pr.curve(scores.class0 = fgSingleRep, scores.class1 = bgSingleRep, curve = T, rand.compute=T)
    
    repX_perf[[zz]] <- data.frame(auroc=rep_roc$auc,
                                  aupr=rep_pr$auc.integral,
                                  rep=paste0("Fold",zz), 
                                  dataString=dataString,
                                  positiveClassSize=sum(multiClass$obs==positiveClass),
                                  negativeClassSize=sum(multiClass$obs==negativeClass),
                                  diseaseType = dzOfInterest,
                                  col2Predict = col2Predict)
  }
  
  # SUMMARIZE MODEL PERFORMANCES
  rep_perf <- do.call(rbind, repX_perf)
  perf <- data.frame(auroc = prroc_roc$auc,
                           aupr = prroc_pr$auc.integral,
                           aucEstimate = out95$cvAUC, # either out95 or out99 work (same result)
                           aucSE95 = out95$se,
                           lowerCI95 = out95$ci[1],
                           upperCI95 = out95$ci[2],
                           aucSE99 = out99$se,
                           lowerCI99 = out99$ci[1],
                           upperCI99 = out99$ci[2],
                           diseaseType = dzOfInterest,
                           dataString=dataString,
                           col2Predict = col2Predict)
   
  print(confusionMatrix(multiClass$pred, multiClass$obs))
  print(multiClassSummary(multiClass, lev = levels(multiClass$obs)))
  
  varImpBestModelDF <- as.data.frame(varImp(mlModel$finalModel, scale = FALSE))
  varImpBestModelDF2 <- rownames_to_column(varImpBestModelDF, "Taxa")
  varImpBestModelDF2Ordered <- varImpBestModelDF2[order(-varImpBestModelDF2$Overall),]
  colnames(varImpBestModelDF2Ordered)[2] <- "varImp"
  varImpBestModelDF2OrderedNonzero <- varImpBestModelDF2Ordered[varImpBestModelDF2Ordered$varImp != 0,]
  
  if(varImpFlag){
    print(head(varImpBestModelDF2OrderedNonzero))
  }
  
  plot(prroc_pr)
  plot(prroc_roc)
  
  res <- list(mlModel=mlModel,
              perf=perf,
              rep_perf=rep_perf,
              multiClass=multiClass,
              varImp=varImpBestModelDF2OrderedNonzero)
  return(res)
  
}

loocvIOR <- function(metaData, 
                     snmData, 
                     dzType = "SKCM",
                     caretTuneGrid = customGBMGrid,
                     savePlotFlag=FALSE,
                     dataStringForPlotFilename=NULL){
  
  defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                                 n.trees = floor((1:3) * 50),
                                 shrinkage = 0.1,
                                 n.minobsinnode = 5)
  customGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                                n.trees = floor((1:3) * 50),
                                shrinkage = 0.1,
                                n.minobsinnode = 1)
  
  if(dzType == "all"){
    mlDataY <- droplevels(metaData[metaData$Responder %in% c("Yes","No"),])
    mlDataX <- snmData[rownames(mlDataY),]
    print(table(mlDataY$Responder))
  } else if(dzType == "NSCLC"){
    mlDataY <- droplevels(metaData[metaData$Responder %in% c("Yes","No") &
                                                             metaData$disease_type_consol == "NSCLC",])
    mlDataX <- snmData[rownames(mlDataY),]
    print(table(mlDataY$Responder))
  } else if(dzType == "SKCM"){
    mlDataY <- droplevels(metaData[metaData$Responder %in% c("Yes","No") &
                                                             metaData$disease_type_consol == "SKCM",])
    mlDataX <- snmData[rownames(mlDataY),]
    print(table(mlDataY$Responder))
  }
  
  # Do LOOCV model building and testing
  
  multiClassSummaryStats <- list()
  multiClassSummaryStatsDist <- list()
  numKFold <- 4
  numResampleIter <- 1
  # metaData <- metaTmpX
  # snmData <- snmData # dataPSUniqueDecontamQC # 
  iterSize <- 1
  
  indexSuper <- 1:dim(mlDataY)[1]
  predProbs <- list()
  obsClass <- vector()
  predClass <- vector()
  # varImpBestModelDF2OrderedNonzeroList <- list()
  
  for(ii in 1:length(indexSuper)){
    print(sprintf("Iteration: %d/%d", ii, length(indexSuper)))
    index <- indexSuper[ii]
    # print(index)
    trainX <- mlDataX[-index,]
    trainY <- mlDataY[-index,]$Responder
    testX <- mlDataX[index,,drop=FALSE]
    testY <- mlDataY[index,,drop=FALSE]$Responder
    print(as.character(testY))
    
    refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
    refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
    
    obsClass[ii] <- as.character(refactoredTestY)
    
    set.seed(42)
    ctrl <- trainControl(method = "repeatedcv",
                         number = numKFold,
                         repeats = numResampleIter,
                         sampling = "up",
                         summaryFunction = twoClassSummary,
                         classProbs = TRUE,
                         verboseIter = FALSE,
                         savePredictions = TRUE,
                         allowParallel=TRUE)
    
    mlModel <- train(x = trainX,
                     y = refactoredTrainY,
                     method = "gbm",
                     preProcess = c("scale","center"),
                     trControl = ctrl,
                     verbose = FALSE,
                     tuneGrid = caretTuneGrid,
                     metric = "ROC")
    
    predProbs[ii] <- list(predict(mlModel, newdata = testX, type = "prob"))
    predClass[ii] <- as.character(predict(mlModel, newdata = testX, type = "raw"))
    
    rm(mlModel)
  }
  
  classes <- c("Yes","No")
  loocvPreds <- cbind(obs = factor(obsClass,
                                   levels = classes),
                      pred = factor(predClass,
                                    levels = classes),
                      do.call(rbind,predProbs))
  multiClassSummaryStats <- multiClassSummary(loocvPreds, lev = classes)
  print(multiClassSummaryStats)
  print(confusionMatrix(loocvPreds$pred, loocvPreds$obs))
  
  positiveClass <- "Yes"
  negativeClass <- "No"
  predProbs <- loocvPreds[,"Yes"]
  fg <- predProbs[loocvPreds$obs == positiveClass]
  bg <- predProbs[loocvPreds$obs == negativeClass]
  
  prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
  
  rocCurveData <- data.frame(prroc_roc$curve)
  colnames(rocCurveData) <- c("fpr","tpr","thresh")
  auroc <- prroc_roc$auc
  prCurveData <- data.frame(prroc_pr$curve)
  colnames(prCurveData) <- c("recall","precision","thresh")
  aupr <- prroc_pr$auc.integral
  
  rocPlot <- ggplot(data = rocCurveData, aes(x = fpr, y = tpr, color=thresh)) + 
    geom_path(size=1) + theme_minimal() +
    scale_color_gradientn(colors = rainbow(10), name = "Probability\ncutoff") +
    geom_abline(intercept = 0, slope = 1, color = "gray", size = 0.5) +
    coord_equal(ratio=1) + xlim(0, 1) + ylim(0,1) +
    labs(x = "False Positive Rate", y = "True Positive Rate") +
    annotate("text", x = 0.3, y = 0.15, 
             color = "black", 
             label = paste0("AUROC: ",
                            paste0(100*round(auroc,4),"%")), hjust = 0)
  
  prPlot <- ggplot(data = prCurveData, aes(x = recall, y = precision, color=thresh)) + 
    geom_path(size=1) + theme_minimal() +
    scale_color_gradientn(colors = rainbow(10), name = "Probability\ncutoff") +
    coord_equal(ratio=1) + xlim(0, 1) + ylim(0,1) +
    labs(x = "Recall", y = "Precision") +
    annotate("text", x = 0.3, y = 0.15, 
             color = "black", 
             label = paste0("AUPR: ",
                            paste0(100*round(aupr,4),"%")), hjust = 0)
  
  combinedPlots <- ggarrange(rocPlot, prPlot, ncol = 2, 
                             common.legend = TRUE, legend = "right")
  combinedPlotsAnnotated <- annotate_figure(combinedPlots, 
                                            top = text_grob("UCSD plasma cohort\nPredicting immunotherapy response (LOOCV)", 
                                                            color = "black", face = "bold", size = 14))
  plot(prroc_pr)
  plot(prroc_roc)
  print(combinedPlotsAnnotated)
  if(savePlotFlag){
    fileName <- paste0("Figures/Figure_5/ucsd_roc_pr_IOR_",
                       dzType,"_",dataStringForPlotFilename,".svg")
    ggsave(filename = fileName, units = "in", width = 8, height = 5)
  }
  
  res <- list(loocvPreds=loocvPreds,
              rocPlot=rocPlot,
              prPlot=prPlot,
              rocCurveData=rocCurveData,
              prCurveData=prCurveData,
              combinedPlots=combinedPlots)
  
  return(res)
}

runAncomBC_TvsNAT <- function(psSeqCenter, ancombcLibCut = 1000, showTopX = 3, qvalCutoff = 0.05, decontamResV2=decontamResultsV2){
  metaData <- data.frame(sample_data(psSeqCenter))
  SeqCenter <- as.character(metaData$data_submitting_center_label[1])
  SeqCenterFormatted <- gsub('([[:punct:]])|\\s+','',SeqCenter)
  expStrategy <- as.character(metaData$experimental_strategy[1])
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
      diff_abn = ifelse(unlist(ancombc_Fungi_TvsNAT_X$res$q_val)<=qvalCutoff, yes = TRUE, no = FALSE),
      genus = rep200TaxSplit_Fungi_Paired_to_Weizmann[row.names(ancombc_Fungi_TvsNAT_X$res$beta),"genus"],
      species = rep200TaxSplit_Fungi_Paired_to_Weizmann[row.names(ancombc_Fungi_TvsNAT_X$res$beta),"species"],
      OGUs = row.names(ancombc_Fungi_TvsNAT_X$res$beta),
      row.names = row.names(ancombc_Fungi_TvsNAT_X$res$beta))
    ancom_res_df_Fungi_X$diff_name_flag <- ancom_res_df_Fungi_X$diff_abn + 0
    ancom_res_df_Fungi_X_sorted <- ancom_res_df_Fungi_X[order(ancom_res_df_Fungi_X$q_val),]
    
    ancom_res_df_Fungi_X_sorted$diff_label <- ""
    ancom_res_df_Fungi_X_sorted$diff_label[1:showTopX] <- ifelse(ancom_res_df_Fungi_X_sorted$diff_abn, 
                                                          yes = paste0(ancom_res_df_Fungi_X_sorted$OGUs,"\n(",ancom_res_df_Fungi_X_sorted$genus,")"),
                                                          no = "")[1:showTopX]
    ancom_res_df_Fungi_X_sorted$origin <- decontamResV2[ancom_res_df_Fungi_X_sorted$OGUs, "reason"]
    
    print(sprintf("Number of differentially abundant fungi: %d", sum(ancom_res_df_Fungi_X_sorted$diff_name_flag)))
    cat("\n")
    
    ## Set up sub-axis tumor (right) vs. NAT (left) labels
    library(grid)
    text_high <- textGrob("Tumor", gp=gpar(fontsize=11, fontface="bold.italic"))
    text_low <- textGrob("NAT", gp=gpar(fontsize=11, fontface="bold.italic"))
    text_high_xpos <- ifelse(max(ancom_res_df_Fungi_X_sorted$beta)<=0.2, yes = 0.2, no = max(ancom_res_df_Fungi_X_sorted$beta))
    text_low_xpos <- ifelse(min(ancom_res_df_Fungi_X_sorted$beta)>=-0.2, yes = -0.2, no = min(ancom_res_df_Fungi_X_sorted$beta))
    text_ypos <- ifelse(1.1*max(-log10(ancom_res_df_Fungi_X_sorted$p_val))<=-log10(0.05), yes = 1.1*-log10(0.05), no = 1.1*max(-log10(ancom_res_df_Fungi_X_sorted$p_val)))
    
    plotFilePath <- "Figures/Supplementary_Figures/"
    ancom_res_df_Fungi_X_sorted %>%
      ggplot(aes(x = beta, y = -log10(p_val), color = diff_abn, label = diff_label, shape=origin)) + geom_point(size = 2) +
      theme_bw() + geom_hline(yintercept=-log10(0.05), col="black", linetype='dashed') +
      geom_vline(xintercept=c(-0.2, 0.2), col="black", linetype='dashed') + 
      theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + coord_fixed(ylim = c(0, NA), clip = "off") +
      scale_color_aaas() +
      labs(x = "log-fold change\n(ANCOM-BC beta)", y = "-Log10(p-value)", 
           color = paste0("Differentially\nabundant (q<=",qvalCutoff,")"), 
           title = paste(paste0(SeqCenter," (",expStrategy,")"), Dz, "", sep = "\n"),
           shape = "Source reference") +
      geom_label_repel(force = 20, size = 2, box.padding = 2, point.padding = 1e-06, label.size = 0.2, show.legend = FALSE, color = "black", max.overlaps = 10) +
      theme(plot.margin = unit(c(1,1,2,1), "lines")) +
      annotation_custom(text_high,xmin=text_high_xpos,xmax=text_high_xpos,ymin=text_ypos, ymax=text_ypos) + 
      annotation_custom(text_low,xmin=text_low_xpos,xmax=text_low_xpos,ymin=text_ypos, ymax=text_ypos) +
      ggsave(filename = paste0(plotFilePath, "ancombc_TvsNAT_", SeqCenterFormatted,"_",DzFormatted,".pdf"), dpi = "retina",
             width = 8, height = 6, units = "in")
    # Write data to file
    dataFilePath <- "Figures_data/Supplementary_Figures/"
    ancom_res_df_Fungi_X_sorted %>% write.csv(file = paste0(dataFilePath, "ancombc_TvsNAT_", SeqCenterFormatted,"_",DzFormatted,".csv"))
  }
}

runAncomBC_1VsAll_Fungi <- function(countData=rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_PT,
                              metaData = metaQiitaCombined_Nonzero_DecontamV2_PT,
                              taxTable = rep200TaxSplit_Fungi_Paired_to_Weizmann,
                              qvalCutoff = 0.05,
                              showTopX = 3,
                              ancombcLibCut = 100,
                              decontamResV2 = decontamResultsV2,
                              fileString = "ancombc_PT_1vsAll_",
                              taxaPlotLabel = "genus"){
  metaData$investigation_formatted <- gsub("^TCGA-","",metaData$investigation)
  
  SeqCenters <- as.character(unique(metaData$data_submitting_center_label))
  expStrategies <- as.character(unique(metaData$experimental_strategy))
  
  for(zz in seq_along(expStrategies)){
    
    expStrategy <- expStrategies[zz]
    # ancombcLibCut <- ifelse(expStrategy == "WGS", yes = ancombcLibCutWGS, no = 0)
    
    for(jj in seq_along(SeqCenters)){
      
      SeqCenter <- SeqCenters[jj]
      SeqCenterFormatted <- gsub('([[:punct:]])|\\s+','',SeqCenter)
      
      metaDataFilt <- metaData %>% filter(data_submitting_center_label == SeqCenter & experimental_strategy == expStrategy) %>% droplevels()
      countDataFilt <- countData[rownames(metaDataFilt),]
      
      cancerTypes <- as.character(unique(metaDataFilt$disease_type))
      
      for(ii in seq_along(cancerTypes)){
        Dz <- cancerTypes[ii]
        DzFormatted <- gsub('([[:punct:]])|\\s+','',Dz)
        
        metaDataFilt$predY <- factor(ifelse(metaDataFilt$disease_type == Dz,
                                            yes = Dz, no = "Other"),
                                     levels = c("Other",Dz))
        countDataFilt_Dz <- countDataFilt[which(metaDataFilt$disease_type == Dz),]
        countDataFilt_Other <- countDataFilt[which(metaDataFilt$disease_type != Dz),]
        investigationText <- metaDataFilt$investigation_formatted[which(metaDataFilt$disease_type == Dz)[1]]
        
        # SKIP CANCER TYPES THAT ONLY HAVE ONE CLASS OF A COMPARISON
        # in case a seq center only has 1 cancer type
        if(length(table(metaDataFilt$predY)) < 2){next}
        
        # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF 10 SAMPLES IN EITHER CLASS
        if(any(table(metaDataFilt$predY) < 10)){next}

        # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF ancombcLibCut READS/SAMPLE IN EITHER CLASS
        if(all(rowSums(countDataFilt_Dz) < ancombcLibCut)){next}
        if(all(rowSums(countDataFilt_Other) < ancombcLibCut)){next}
        
        # If sufficient samples, then:
        print(SeqCenter)
        print(Dz)
        print(sprintf("Number of samples (Dz|Other): %d | %d", 
                      unname(table(metaDataFilt$predY)[2]),
                      unname(table(metaDataFilt$predY)[1])))
        
        psFungiTCGADecontam_1vsAll <- phyloseq(otu_table(countDataFilt, taxa_are_rows = FALSE), 
                                               tax_table(as.matrix(taxTable)), 
                                               sample_data(metaDataFilt))
        # print(psFungiTCGADecontam_1vsAll)
        print(sprintf("Read count distribution:"))
        print(summary(rowSums(otu_table(psFungiTCGADecontam_1vsAll)))) # Sample read distribution
        print(sprintf("Now running ANCOM-BC..."))
        ancombc_Fungi_1vsAll_X <- ancombc(phyloseq = psFungiTCGADecontam_1vsAll, 
                                          formula = "predY", p_adj_method = "BH", zero_cut = 0.999, 
                                          lib_cut = ancombcLibCut, 
                                          # group = "sample_type", struc_zero = FALSE, neg_lb = FALSE,
                                          tol = 1e-5, max_iter = 100, conserve = FALSE, alpha = 0.05, global = FALSE)
        print(sprintf("Plotting and saving data..."))
        ancom_res_df_Fungi_1vsAll_X <- data.frame(
          beta = unlist(ancombc_Fungi_1vsAll_X$res$beta),
          se = unlist(ancombc_Fungi_1vsAll_X$res$se),
          W = unlist(ancombc_Fungi_1vsAll_X$res$W),
          p_val = unlist(ancombc_Fungi_1vsAll_X$res$p_val),
          q_val = unlist(ancombc_Fungi_1vsAll_X$res$q_val),
          diff_abn = ifelse(unlist(ancombc_Fungi_1vsAll_X$res$q_val)<=qvalCutoff, yes = TRUE, no = FALSE),
          genus = taxTable[row.names(ancombc_Fungi_1vsAll_X$res$beta),"genus"],
          species = taxTable[row.names(ancombc_Fungi_1vsAll_X$res$beta),"species"],
          OGUs = row.names(ancombc_Fungi_1vsAll_X$res$beta),
          row.names = row.names(ancombc_Fungi_1vsAll_X$res$beta))
        ancom_res_df_Fungi_1vsAll_X$diff_name_flag <- ancom_res_df_Fungi_1vsAll_X$diff_abn + 0
        ancom_res_df_Fungi_1vsAll_X_sorted <- ancom_res_df_Fungi_1vsAll_X[order(ancom_res_df_Fungi_1vsAll_X$q_val),]
        
        ancom_res_df_Fungi_1vsAll_X_sorted$diff_label <- ""
        ancom_res_df_Fungi_1vsAll_X_sorted$diff_label[1:showTopX] <- ifelse(ancom_res_df_Fungi_1vsAll_X_sorted$diff_abn, 
                                                                            yes = paste0(ancom_res_df_Fungi_1vsAll_X_sorted$OGUs,"\n(",ancom_res_df_Fungi_1vsAll_X_sorted[,taxaPlotLabel],")"),
                                                                            no = "")[1:showTopX]
        ancom_res_df_Fungi_1vsAll_X_sorted$origin <- decontamResV2[ancom_res_df_Fungi_1vsAll_X_sorted$OGUs, "reason"]
        
        print(sprintf("Number of differentially abundant fungi: %d", sum(ancom_res_df_Fungi_1vsAll_X_sorted$diff_name_flag)))
        cat("\n")
        
        ## Set up sub-axis CT (right) vs. Others (left) labels
        library(grid)
        text_high <- textGrob(investigationText, gp=gpar(fontsize=11, fontface="bold.italic"))
        text_low <- textGrob("Other", gp=gpar(fontsize=11, fontface="bold.italic"))
        text_high_xpos <- ifelse(max(ancom_res_df_Fungi_1vsAll_X_sorted$beta)<=0.2, yes = 0.2, no = max(ancom_res_df_Fungi_1vsAll_X_sorted$beta))
        text_low_xpos <- ifelse(min(ancom_res_df_Fungi_1vsAll_X_sorted$beta)>=-0.2, yes = -0.2, no = min(ancom_res_df_Fungi_1vsAll_X_sorted$beta))
        yval <- -log10(ancom_res_df_Fungi_1vsAll_X_sorted$p_val)
        text_ypos <- ifelse(1.1*max(yval)<=-log10(0.05), 
                            yes = 1.1*-log10(0.05), no = 1.1*max(yval))
        text_ypos <- ifelse(is.finite(text_ypos), yes = text_ypos, no = 1.1*max(yval[is.finite(yval)]))
        
        plotFilePath <- "Figures/Supplementary_Figures/"
        ancom_res_df_Fungi_1vsAll_X_sorted %>%
          ggplot(aes(x = beta, y = -log10(p_val), color = diff_abn, label = diff_label, shape=origin)) + geom_point(size = 2) +
          theme_bw() + geom_hline(yintercept=-log10(0.05), col="black", linetype='dashed') +
          geom_vline(xintercept=c(-0.2, 0.2), col="black", linetype='dashed') + 
          theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + coord_fixed(ylim = c(0, NA), clip = "off") +
          scale_color_aaas() +
          labs(x = "log-fold change\n(ANCOM-BC beta)", y = "-Log10(p-value)", 
               color = paste0("Differentially\nabundant\nfungi\n(q<=",qvalCutoff,")"), 
               title = paste(paste0(SeqCenter," (",expStrategy,")"), Dz, "", sep = "\n"),
               shape = "Source reference") +
          geom_label_repel(force = 50, size = 2, box.padding = 2, point.padding = 1e-06, label.size = 0.2, show.legend = FALSE, color = "black", max.overlaps = 10) +
          theme(plot.margin = unit(c(1,1,2,1), "lines")) +
          annotation_custom(text_high,xmin=text_high_xpos,xmax=text_high_xpos,ymin=text_ypos, ymax=text_ypos) + 
          annotation_custom(text_low,xmin=text_low_xpos,xmax=text_low_xpos,ymin=text_ypos, ymax=text_ypos) +
          ggsave(filename = paste0(plotFilePath, fileString, SeqCenterFormatted,"_",DzFormatted,".pdf"), dpi = "retina",
                 width = 8, height = 6, units = "in")
        # Write data to file
        dataFilePath <- "Figures_data/Supplementary_Figures/"
        ancom_res_df_Fungi_1vsAll_X_sorted %>% write.csv(file = paste0(dataFilePath, fileString, SeqCenterFormatted,"_",DzFormatted,".csv"))
      }
    }
  }
}

runAncomBC_1VsAll_Bacteria <- function(countData=rep200Data_Matched2ImmunePT_Bacteria_WIS,
                              metaData = metaImmunePT,
                              taxTable = rep200TaxSplit_Bacteria_Formatted,
                              qvalCutoff = 0.05,
                              showTopX = 3,
                              ancombcLibCut = 100,
                              fileString = "WIS_bacteria_ancombc_PT_1vsAll_",
                              taxaPlotLabel = "genus"){
  metaData$investigation_formatted <- gsub("^TCGA-","",metaData$investigation)
  
  SeqCenters <- as.character(unique(metaData$data_submitting_center_label))
  expStrategies <- as.character(unique(metaData$experimental_strategy))
  
  for(zz in seq_along(expStrategies)){
    
    expStrategy <- expStrategies[zz]
    # ancombcLibCut <- ifelse(expStrategy == "WGS", yes = 1000, no = 0)
    
    for(jj in seq_along(SeqCenters)){
      
      SeqCenter <- SeqCenters[jj]
      SeqCenterFormatted <- gsub('([[:punct:]])|\\s+','',SeqCenter)
      
      metaDataFilt <- metaData %>% filter(data_submitting_center_label == SeqCenter & experimental_strategy == expStrategy) %>% droplevels()
      countDataFilt <- countData[rownames(metaDataFilt),]
      
      cancerTypes <- as.character(unique(metaDataFilt$disease_type))
      
      for(ii in seq_along(cancerTypes)){
        Dz <- cancerTypes[ii]
        DzFormatted <- gsub('([[:punct:]])|\\s+','',Dz)
        
        metaDataFilt$predY <- factor(ifelse(metaDataFilt$disease_type == Dz,
                                            yes = Dz, no = "Other"),
                                     levels = c("Other",Dz))
        countDataFilt_Dz <- countDataFilt[which(metaDataFilt$disease_type == Dz),]
        countDataFilt_Other <- countDataFilt[which(metaDataFilt$disease_type != Dz),]
        investigationText <- metaDataFilt$investigation_formatted[which(metaDataFilt$disease_type == Dz)[1]]
        
        # SKIP CANCER TYPES THAT ONLY HAVE ONE CLASS OF A COMPARISON
        # in case a seq center only has 1 cancer type
        if(length(table(metaDataFilt$predY)) < 2){next}
        
        # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF 10 SAMPLES IN EITHER CLASS
        if(any(table(metaDataFilt$predY) < 10)){next}

        # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF ancombcLibCut READS/SAMPLE IN EITHER CLASS
        if(all(rowSums(countDataFilt_Dz) < ancombcLibCut)){next}
        if(all(rowSums(countDataFilt_Other) < ancombcLibCut)){next}
        
        # If sufficient samples, then:
        print(SeqCenter)
        print(Dz)
        print(sprintf("Number of samples (Dz|Other): %d | %d", 
                      unname(table(metaDataFilt$predY)[2]),
                      unname(table(metaDataFilt$predY)[1])))
        
        psFungiTCGADecontam_1vsAll <- phyloseq(otu_table(countDataFilt, taxa_are_rows = FALSE), 
                                               tax_table(as.matrix(taxTable)), 
                                               sample_data(metaDataFilt))
        # print(psFungiTCGADecontam_1vsAll)
        print(sprintf("Read count distribution:"))
        print(summary(rowSums(otu_table(psFungiTCGADecontam_1vsAll)))) # Sample read distribution
        print(sprintf("Now running ANCOM-BC..."))
        ancombc_Bacteria_1vsAll_X <- ancombc(phyloseq = psFungiTCGADecontam_1vsAll, 
                                          formula = "predY", p_adj_method = "BH", zero_cut = 0.999, 
                                          lib_cut = ancombcLibCut, 
                                          # group = "sample_type", struc_zero = FALSE, neg_lb = FALSE,
                                          tol = 1e-5, max_iter = 100, conserve = FALSE, alpha = 0.05, global = FALSE)
        print(sprintf("Plotting and saving data..."))
        ancom_res_df_Bacteria_1vsAll_X <- data.frame(
          beta = unlist(ancombc_Bacteria_1vsAll_X$res$beta),
          se = unlist(ancombc_Bacteria_1vsAll_X$res$se),
          W = unlist(ancombc_Bacteria_1vsAll_X$res$W),
          p_val = unlist(ancombc_Bacteria_1vsAll_X$res$p_val),
          q_val = unlist(ancombc_Bacteria_1vsAll_X$res$q_val),
          diff_abn = ifelse(unlist(ancombc_Bacteria_1vsAll_X$res$q_val)<=qvalCutoff, yes = TRUE, no = FALSE),
          genus = taxTable[row.names(ancombc_Bacteria_1vsAll_X$res$beta),"genus"],
          species = taxTable[row.names(ancombc_Bacteria_1vsAll_X$res$beta),"species"],
          OGUs = row.names(ancombc_Bacteria_1vsAll_X$res$beta),
          row.names = row.names(ancombc_Bacteria_1vsAll_X$res$beta))
        ancom_res_df_Bacteria_1vsAll_X$diff_name_flag <- ancom_res_df_Bacteria_1vsAll_X$diff_abn + 0
        ancom_res_df_Bacteria_1vsAll_X_sorted <- ancom_res_df_Bacteria_1vsAll_X[order(ancom_res_df_Bacteria_1vsAll_X$q_val),]
        
        ancom_res_df_Bacteria_1vsAll_X_sorted$diff_label <- ""
        ancom_res_df_Bacteria_1vsAll_X_sorted$diff_label[1:showTopX] <- ifelse(ancom_res_df_Bacteria_1vsAll_X_sorted$diff_abn, 
                                                                            yes = paste0(ancom_res_df_Bacteria_1vsAll_X_sorted$OGUs,"\n(",ancom_res_df_Bacteria_1vsAll_X_sorted[,taxaPlotLabel],")"),
                                                                            no = "")[1:showTopX]
        
        print(sprintf("Number of differentially abundant bacteria: %d", sum(ancom_res_df_Bacteria_1vsAll_X_sorted$diff_name_flag)))
        cat("\n")
        
        ## Set up sub-axis CT (right) vs. Others (left) labels
        library(grid)
        text_high <- textGrob(investigationText, gp=gpar(fontsize=11, fontface="bold.italic"))
        text_low <- textGrob("Other", gp=gpar(fontsize=11, fontface="bold.italic"))
        text_high_xpos <- ifelse(max(ancom_res_df_Bacteria_1vsAll_X_sorted$beta)<=0.2, yes = 0.2, no = max(ancom_res_df_Bacteria_1vsAll_X_sorted$beta))
        text_low_xpos <- ifelse(min(ancom_res_df_Bacteria_1vsAll_X_sorted$beta)>=-0.2, yes = -0.2, no = min(ancom_res_df_Bacteria_1vsAll_X_sorted$beta))
        yval <- -log10(ancom_res_df_Bacteria_1vsAll_X_sorted$p_val)
        text_ypos <- ifelse(1.1*max(yval)<=-log10(0.05), 
                            yes = 1.1*-log10(0.05), no = 1.1*max(yval))
        text_ypos <- ifelse(is.finite(text_ypos), yes = text_ypos, no = 1.1*max(yval[is.finite(yval)]))
        
        plotFilePath <- "Figures/Supplementary_Figures/"
        ancom_res_df_Bacteria_1vsAll_X_sorted %>%
          ggplot(aes(x = beta, y = -log10(p_val), color = diff_abn, label = diff_label)) + geom_point(size = 2) +
          theme_bw() + geom_hline(yintercept=-log10(0.05), col="black", linetype='dashed') +
          geom_vline(xintercept=c(-0.2, 0.2), col="black", linetype='dashed') + 
          theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + coord_fixed(ylim = c(0, NA), clip = "off") +
          scale_color_aaas() +
          labs(x = "log-fold change\n(ANCOM-BC beta)", y = "-Log10(p-value)", 
               color = paste0("Differentially\nabundant\nbacteria\n(q<=",qvalCutoff,")"), 
               title = paste(paste0(SeqCenter," (",expStrategy,")"), Dz, "", sep = "\n")) +
          geom_label_repel(force = 20, size = 2, box.padding = 2, point.padding = 1e-06, label.size = 0.2, show.legend = FALSE, color = "black", max.overlaps = 10) +
          theme(plot.margin = unit(c(1,1,2,1), "lines")) +
          annotation_custom(text_high,xmin=text_high_xpos,xmax=text_high_xpos,ymin=text_ypos, ymax=text_ypos) + 
          annotation_custom(text_low,xmin=text_low_xpos,xmax=text_low_xpos,ymin=text_ypos, ymax=text_ypos) +
          ggsave(filename = paste0(plotFilePath, fileString, SeqCenterFormatted,"_",DzFormatted,".pdf"), dpi = "retina",
                 width = 8, height = 6, units = "in")
        # Write data to file
        dataFilePath <- "Figures_data/Supplementary_Figures/"
        ancom_res_df_Bacteria_1vsAll_X_sorted %>% write.csv(file = paste0(dataFilePath, fileString, SeqCenterFormatted,"_",DzFormatted,".csv"))
      }
    }
  }
}

runAncomBC_Stage_Fungi <- function(countData=rep200Data_WGS_RNA_HiSeq_Fungi_DecontamV2_Nonzero_Path_PT,
                                    metaDataPath = metaQiitaCombined_Nonzero_DecontamV2_Path_PT,
                                    taxTable = rep200TaxSplit_Fungi_Paired_to_Weizmann,
                                    stagesCmp = c("Stage I","Stage IV"),
                                    stageAllFlag = FALSE,
                                    qvalCutoff = 0.05,
                                    showTopX = 3,
                                    ancombcLibCutWGS = 100,
                                    decontamResV2 = decontamResultsV2,
                                    fileString = "ancombc_PT_Stage_",
                                    taxaPlotLabel = "genus"){
  metaDataPath$investigation_formatted <- gsub("^TCGA-","",metaDataPath$investigation)
  stagesCmpFormatted <- gsub('([[:punct:]])|\\s+','',stagesCmp)
  
  SeqCenters <- as.character(unique(metaDataPath$data_submitting_center_label))
  expStrategies <- as.character(unique(metaDataPath$experimental_strategy))
  
  for(zz in seq_along(expStrategies)){
    
    expStrategy <- expStrategies[zz]
    ancombcLibCut <- ifelse(expStrategy == "WGS", yes = ancombcLibCutWGS, no = 0)
    
    for(jj in seq_along(SeqCenters)){
      
      SeqCenter <- SeqCenters[jj]
      SeqCenterFormatted <- gsub('([[:punct:]])|\\s+','',SeqCenter)
      
      metaDataPathFilt <- metaDataPath %>% filter(data_submitting_center_label == SeqCenter,
                                                  experimental_strategy == expStrategy,
                                                  pathologic_stage_label_binned %in% stagesCmpFormatted) %>% droplevels()
      countDataFilt <- countData[rownames(metaDataPathFilt),]
      
      cancerTypes <- as.character(unique(metaDataPathFilt$disease_type))
      
      for(ii in seq_along(cancerTypes)){
        Dz <- cancerTypes[ii]
        DzFormatted <- gsub('([[:punct:]])|\\s+','',Dz)
        
        metaDataPathFiltDz <- metaDataPathFilt %>% filter(disease_type == Dz) %>% droplevels()
        if(stageAllFlag){
          metaDataPathFiltDz$predY <- factor(ifelse(metaDataPathFiltDz$pathologic_stage_label_binned %in% c("StageI","StageII"),
                                             yes = "EarlyStage", no = "LateStage"), levels = c("EarlyStage","LateStage"))
          textRight <- "Late Stage"
          textLeft <- "Early Stage"
        } else{
          metaDataPathFiltDz$predY <- metaDataPathFiltDz$pathologic_stage_label_binned
          textRight <- stagesCmp[2]
          textLeft <- stagesCmp[1]
        }
        
        countDataFilt_Dz <- countDataFilt[rownames(metaDataPathFiltDz),]
        
        # SKIP CANCER TYPES THAT ONLY HAVE ONE CLASS OF A COMPARISON
        # in case a seq center only has 1 stage type
        if(length(table(metaDataPathFiltDz$predY)) < 2){next}
        
        # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF 10 SAMPLES PER STAGE IN EITHER CLASS
        if(any(table(metaDataPathFiltDz$predY) < 10)){next}
        
        # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF ancombcLibCut READS/SAMPLE IN EITHER CLASS
        if(all(rowSums(countDataFilt_Dz) < ancombcLibCut)){next}
        
        # If sufficient samples, then:
        print(SeqCenter)
        print(Dz)
        print(sprintf("Number of samples (%s|%s): %d | %d", 
                      names(table(metaDataPathFiltDz$predY)[1]), 
                      names(table(metaDataPathFiltDz$predY)[2]), 
                      unname(table(metaDataPathFiltDz$predY)[1]),
                      unname(table(metaDataPathFiltDz$predY)[2])))
        
        psFungiTCGADecontam_Stage <- phyloseq(otu_table(countDataFilt_Dz, taxa_are_rows = FALSE), 
                                               tax_table(as.matrix(taxTable)), 
                                               sample_data(metaDataPathFiltDz))
        # print(psFungiTCGADecontam_1vsAll)
        print(sprintf("Read count distribution:"))
        print(summary(rowSums(otu_table(psFungiTCGADecontam_Stage)))) # Sample read distribution
        print(sprintf("Now running ANCOM-BC..."))
        ancombc_Fungi_Stage_X <- ancombc(phyloseq = psFungiTCGADecontam_Stage, 
                                          formula = "predY", p_adj_method = "BH", zero_cut = 0.999, 
                                          lib_cut = ancombcLibCut, 
                                          # group = "sample_type", struc_zero = FALSE, neg_lb = FALSE,
                                          tol = 1e-5, max_iter = 100, conserve = FALSE, alpha = 0.05, global = FALSE)
        print(sprintf("Plotting and saving data..."))
        ancom_res_df_Fungi_Stage_X <- data.frame(
          beta = unlist(ancombc_Fungi_Stage_X$res$beta),
          se = unlist(ancombc_Fungi_Stage_X$res$se),
          W = unlist(ancombc_Fungi_Stage_X$res$W),
          p_val = unlist(ancombc_Fungi_Stage_X$res$p_val),
          q_val = unlist(ancombc_Fungi_Stage_X$res$q_val),
          diff_abn = ifelse(unlist(ancombc_Fungi_Stage_X$res$q_val)<=qvalCutoff, yes = TRUE, no = FALSE),
          genus = taxTable[row.names(ancombc_Fungi_Stage_X$res$beta),"genus"],
          species = taxTable[row.names(ancombc_Fungi_Stage_X$res$beta),"species"],
          OGUs = row.names(ancombc_Fungi_Stage_X$res$beta),
          row.names = row.names(ancombc_Fungi_Stage_X$res$beta))
        ancom_res_df_Fungi_Stage_X$diff_name_flag <- ancom_res_df_Fungi_Stage_X$diff_abn + 0
        ancom_res_df_Fungi_Stage_X_sorted <- ancom_res_df_Fungi_Stage_X[order(ancom_res_df_Fungi_Stage_X$q_val),]
        
        ancom_res_df_Fungi_Stage_X_sorted$diff_label <- ""
        ancom_res_df_Fungi_Stage_X_sorted$diff_label[1:showTopX] <- ifelse(ancom_res_df_Fungi_Stage_X_sorted$diff_abn, 
                                                                            yes = paste0(ancom_res_df_Fungi_Stage_X_sorted$OGUs,"\n(",ancom_res_df_Fungi_Stage_X_sorted[,taxaPlotLabel],")"),
                                                                            no = "")[1:showTopX]
        ancom_res_df_Fungi_Stage_X_sorted$origin <- decontamResV2[ancom_res_df_Fungi_Stage_X_sorted$OGUs, "reason"]
        
        print(sprintf("Number of differentially abundant fungi: %d", sum(ancom_res_df_Fungi_Stage_X_sorted$diff_name_flag)))
        cat("\n")
        
        ## Set up sub-axis CT (right) vs. Others (left) labels
        library(grid)
        text_high <- textGrob(textRight, gp=gpar(fontsize=11, fontface="bold.italic"))
        text_low <- textGrob(textLeft, gp=gpar(fontsize=11, fontface="bold.italic"))
        text_high_xpos <- ifelse(max(ancom_res_df_Fungi_Stage_X_sorted$beta)<=0.2, yes = 0.2, no = max(ancom_res_df_Fungi_Stage_X_sorted$beta))
        text_low_xpos <- ifelse(min(ancom_res_df_Fungi_Stage_X_sorted$beta)>=-0.2, yes = -0.2, no = min(ancom_res_df_Fungi_Stage_X_sorted$beta))
        yval <- -log10(ancom_res_df_Fungi_Stage_X_sorted$p_val)
        text_ypos <- ifelse(1.1*max(yval)<=-log10(0.05), 
                            yes = 1.1*-log10(0.05), no = 1.1*max(yval))
        text_ypos <- ifelse(is.finite(text_ypos), yes = text_ypos, no = 1.1*max(yval[is.finite(yval)]))
        
        plotFilePath <- "Figures/Supplementary_Figures/"
        ancom_res_df_Fungi_Stage_X_sorted %>%
          ggplot(aes(x = beta, y = -log10(p_val), color = diff_abn, label = diff_label, shape=origin)) + geom_point(size = 2) +
          theme_bw() + geom_hline(yintercept=-log10(0.05), col="black", linetype='dashed') +
          geom_vline(xintercept=c(-0.2, 0.2), col="black", linetype='dashed') + 
          theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + coord_fixed(ylim = c(0, NA), clip = "off") +
          scale_color_aaas() +
          labs(x = "log-fold change\n(ANCOM-BC beta)", y = "-Log10(p-value)", 
               color = paste0("Differentially\nabundant\nfungi\n(q<=",qvalCutoff,")"), 
               title = paste(paste0(SeqCenter," (",expStrategy,")"), Dz, "", sep = "\n"),
               shape = "Source reference") +
          geom_label_repel(force = 50, size = 2, box.padding = 2, point.padding = 1e-06, label.size = 0.2, show.legend = FALSE, color = "black", max.overlaps = 10) +
          theme(plot.margin = unit(c(1,1,2,1), "lines")) +
          annotation_custom(text_high,xmin=text_high_xpos,xmax=text_high_xpos,ymin=text_ypos, ymax=text_ypos) + 
          annotation_custom(text_low,xmin=text_low_xpos,xmax=text_low_xpos,ymin=text_ypos, ymax=text_ypos) +
          ggsave(filename = paste0(plotFilePath, fileString, SeqCenterFormatted,"_",DzFormatted,".pdf"), dpi = "retina",
                 width = 8, height = 6, units = "in")
        # Write data to file
        dataFilePath <- "Figures_data/Supplementary_Figures/"
        ancom_res_df_Fungi_Stage_X_sorted %>% write.csv(file = paste0(dataFilePath, fileString, SeqCenterFormatted,"_",DzFormatted,".csv"))
      }
    }
  }
}

runAncomBC_Stage_Bacteria <- function(countData=rep200Data_Matched2ImmunePT_Bacteria_WIS,
                              metaDataPath = metaImmunePT,
                              taxTable = rep200TaxSplit_Bacteria_Formatted,
                              stagesCmp = c("Stage I","Stage IV"),
                              stageAllFlag = FALSE,
                              qvalCutoff = 0.05,
                              showTopX = 3,
                              ancombcLibCutWGS = 100,
                              fileString = "WIS_bacteria_ancombc_PT_Stage_",
                              taxaPlotLabel = "genus"){
  metaDataPath$investigation_formatted <- gsub("^TCGA-","",metaDataPath$investigation)
  stagesCmpFormatted <- gsub('([[:punct:]])|\\s+','',stagesCmp)
  
  SeqCenters <- as.character(unique(metaDataPath$data_submitting_center_label))
  expStrategies <- as.character(unique(metaDataPath$experimental_strategy))
  
  for(zz in seq_along(expStrategies)){
    
    expStrategy <- expStrategies[zz]
    ancombcLibCut <- ifelse(expStrategy == "WGS", yes = ancombcLibCutWGS, no = 0)
    
    for(jj in seq_along(SeqCenters)){
      
      SeqCenter <- SeqCenters[jj]
      SeqCenterFormatted <- gsub('([[:punct:]])|\\s+','',SeqCenter)
      
      metaDataPathFilt <- metaDataPath %>% filter(data_submitting_center_label == SeqCenter,
                                                  experimental_strategy == expStrategy,
                                                  pathologic_stage_label_binned %in% stagesCmpFormatted) %>% droplevels()
      countDataFilt <- countData[rownames(metaDataPathFilt),]
      
      cancerTypes <- as.character(unique(metaDataPathFilt$disease_type))
      
      for(ii in seq_along(cancerTypes)){
        Dz <- cancerTypes[ii]
        DzFormatted <- gsub('([[:punct:]])|\\s+','',Dz)
        
        metaDataPathFiltDz <- metaDataPathFilt %>% filter(disease_type == Dz) %>% droplevels()
        if(stageAllFlag){
          metaDataPathFiltDz$predY <- factor(ifelse(metaDataPathFiltDz$pathologic_stage_label_binned %in% c("StageI","StageII"),
                                             yes = "EarlyStage", no = "LateStage"), levels = c("EarlyStage","LateStage"))
          textRight <- "Late Stage"
          textLeft <- "Early Stage"
        } else{
          metaDataPathFiltDz$predY <- metaDataPathFiltDz$pathologic_stage_label_binned
          textRight <- stagesCmp[2]
          textLeft <- stagesCmp[1]
        }

        countDataFilt_Dz <- countDataFilt[rownames(metaDataPathFiltDz),]
        
        # SKIP CANCER TYPES THAT ONLY HAVE ONE CLASS OF A COMPARISON
        # in case a seq center only has 1 cancer type
        if(length(table(metaDataPathFiltDz$predY)) < 2){next}
        
        # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF 10 SAMPLES IN EITHER CLASS
        if(any(table(metaDataPathFiltDz$predY) < 10)){next}

        # SKIP CANCER TYPES THAT DO NOT HAVE A MINIMUM OF ancombcLibCut READS/SAMPLE IN EITHER CLASS
        if(all(rowSums(countDataFilt_Dz) < ancombcLibCut)){next}
        
        # If sufficient samples, then:
        print(SeqCenter)
        print(Dz)
        print(sprintf("Number of samples (%s|%s): %d | %d", 
                      names(table(metaDataPathFiltDz$predY)[1]), 
                      names(table(metaDataPathFiltDz$predY)[2]), 
                      unname(table(metaDataPathFiltDz$predY)[1]),
                      unname(table(metaDataPathFiltDz$predY)[2])))
        
        psFungiTCGADecontam_Stage <- phyloseq(otu_table(countDataFilt_Dz, taxa_are_rows = FALSE), 
                                               tax_table(as.matrix(taxTable)), 
                                               sample_data(metaDataPathFiltDz))
        # print(psFungiTCGADecontam_Stage)
        print(sprintf("Read count distribution:"))
        print(summary(rowSums(otu_table(psFungiTCGADecontam_Stage)))) # Sample read distribution
        print(sprintf("Now running ANCOM-BC..."))
        ancombc_Bacteria_Stage_X <- ancombc(phyloseq = psFungiTCGADecontam_Stage, 
                                          formula = "predY", p_adj_method = "BH", zero_cut = 0.999, 
                                          lib_cut = ancombcLibCut, 
                                          tol = 1e-5, max_iter = 100, conserve = FALSE, alpha = 0.05, global = FALSE)
        print(sprintf("Plotting and saving data..."))
        ancom_res_df_Bacteria_Stage_X <- data.frame(
          beta = unlist(ancombc_Bacteria_Stage_X$res$beta),
          se = unlist(ancombc_Bacteria_Stage_X$res$se),
          W = unlist(ancombc_Bacteria_Stage_X$res$W),
          p_val = unlist(ancombc_Bacteria_Stage_X$res$p_val),
          q_val = unlist(ancombc_Bacteria_Stage_X$res$q_val),
          diff_abn = ifelse(unlist(ancombc_Bacteria_Stage_X$res$q_val)<=qvalCutoff, yes = TRUE, no = FALSE),
          genus = taxTable[row.names(ancombc_Bacteria_Stage_X$res$beta),"genus"],
          species = taxTable[row.names(ancombc_Bacteria_Stage_X$res$beta),"species"],
          OGUs = row.names(ancombc_Bacteria_Stage_X$res$beta),
          row.names = row.names(ancombc_Bacteria_Stage_X$res$beta))
        ancom_res_df_Bacteria_Stage_X$diff_name_flag <- ancom_res_df_Bacteria_Stage_X$diff_abn + 0
        ancom_res_df_Bacteria_Stage_X_sorted <- ancom_res_df_Bacteria_Stage_X[order(ancom_res_df_Bacteria_Stage_X$q_val),]
        
        ancom_res_df_Bacteria_Stage_X_sorted$diff_label <- ""
        ancom_res_df_Bacteria_Stage_X_sorted$diff_label[1:showTopX] <- ifelse(ancom_res_df_Bacteria_Stage_X_sorted$diff_abn, 
                                                                            yes = paste0(ancom_res_df_Bacteria_Stage_X_sorted$OGUs,"\n(",ancom_res_df_Bacteria_Stage_X_sorted[,taxaPlotLabel],")"),
                                                                            no = "")[1:showTopX]
        
        print(sprintf("Number of differentially abundant bacteria: %d", sum(ancom_res_df_Bacteria_Stage_X_sorted$diff_name_flag)))
        cat("\n")
        
        ## Set up sub-axis CT (right) vs. Others (left) labels
        library(grid)
        text_high <- textGrob(textRight, gp=gpar(fontsize=11, fontface="bold.italic"))
        text_low <- textGrob(textLeft, gp=gpar(fontsize=11, fontface="bold.italic"))
        text_high_xpos <- ifelse(max(ancom_res_df_Bacteria_Stage_X_sorted$beta)<=0.2, yes = 0.2, no = max(ancom_res_df_Bacteria_Stage_X_sorted$beta))
        text_low_xpos <- ifelse(min(ancom_res_df_Bacteria_Stage_X_sorted$beta)>=-0.2, yes = -0.2, no = min(ancom_res_df_Bacteria_Stage_X_sorted$beta))
        yval <- -log10(ancom_res_df_Bacteria_Stage_X_sorted$p_val)
        text_ypos <- ifelse(1.1*max(yval)<=-log10(0.05), 
                            yes = 1.1*-log10(0.05), no = 1.1*max(yval))
        text_ypos <- ifelse(is.finite(text_ypos), yes = text_ypos, no = 1.1*max(yval[is.finite(yval)]))
        
        plotFilePath <- "Figures/Supplementary_Figures/"
        ancom_res_df_Bacteria_Stage_X_sorted %>%
          ggplot(aes(x = beta, y = -log10(p_val), color = diff_abn, label = diff_label)) + geom_point(size = 2) +
          theme_bw() + geom_hline(yintercept=-log10(0.05), col="black", linetype='dashed') +
          geom_vline(xintercept=c(-0.2, 0.2), col="black", linetype='dashed') + 
          theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + coord_fixed(ylim = c(0, NA), clip = "off") +
          scale_color_aaas() +
          labs(x = "log-fold change\n(ANCOM-BC beta)", y = "-Log10(p-value)", 
               color = paste0("Differentially\nabundant\nbacteria\n(q<=",qvalCutoff,")"), 
               title = paste(paste0(SeqCenter," (",expStrategy,")"), Dz, "", sep = "\n")) +
          geom_label_repel(force = 20, size = 2, box.padding = 2, point.padding = 1e-06, label.size = 0.2, show.legend = FALSE, color = "black", max.overlaps = 10) +
          theme(plot.margin = unit(c(1,1,2,1), "lines")) +
          annotation_custom(text_high,xmin=text_high_xpos,xmax=text_high_xpos,ymin=text_ypos, ymax=text_ypos) + 
          annotation_custom(text_low,xmin=text_low_xpos,xmax=text_low_xpos,ymin=text_ypos, ymax=text_ypos) +
          ggsave(filename = paste0(plotFilePath, fileString, SeqCenterFormatted,"_",DzFormatted,".pdf"), dpi = "retina",
                 width = 8, height = 6, units = "in")
        # Write data to file
        dataFilePath <- "Figures_data/Supplementary_Figures/"
        ancom_res_df_Bacteria_Stage_X_sorted %>% write.csv(file = paste0(dataFilePath, fileString, SeqCenterFormatted,"_",DzFormatted,".csv"))
      }
    }
  }
}

wzML1VsAll10k <- function(metaData,
                             snmData,
                             col2Predict = "tissue",
                             dzOfInterest = "breast",
                             modelType = "gbm",
                             trainTestFlag = FALSE,
                             cutPoint = 0.5, 
                             numKFold = 10,
                             numResampleIter = 1,
                             tumorVsNormalFlag = FALSE,
                             tumorVsNormalClasses = c("tumor","nat"),
                             pcrBatchFlag = FALSE,
                             pcrBatchID = "T16",
                             varImpFlag = FALSE){
  require(caret) # for model building
  require(gbm) # for machine learning
  require(PRROC) # for precision-recall curves
  require(MLmetrics) # for multiclass ML
  require(e1071)
  require(ggpubr)
  require(tibble)
  
  if(pcrBatchFlag){
    metaData <- droplevels(metaData[metaData$pcr_batch == pcrBatchID,])
  }
  
  if(tumorVsNormalFlag){
    dzOfInterestFixed <- gsub('([[:punct:]])|\\s+','',dzOfInterest)
    dt <- dzOfInterestFixed
    metaData <- droplevels(metaData[metaData$tissue == dzOfInterestFixed &
                                      metaData$type.detail %in% tumorVsNormalClasses,])
    metaData$condition <- factor(metaData$type.detail,
                                 levels = tumorVsNormalClasses)
    print(table(metaData$condition))
    metaDataFixed <- metaData
    positiveClass <- tumorVsNormalClasses[1]
    negativeClass <- tumorVsNormalClasses[2]
    
  } else{
    dzOfInterestFixed <- gsub('([[:punct:]])|\\s+','',dzOfInterest)
    dt <- dzOfInterestFixed
    metaData$condition <- factor(ifelse(metaData[,col2Predict] == dzOfInterest,
                                        yes = dzOfInterestFixed, no = "Other"),
                                 levels = c(dzOfInterestFixed,"Other"))
    print(table(metaData$condition))
    metaDataFixed <- metaData
    positiveClass <- dzOfInterestFixed
    negativeClass <- "Other"
  }
  
  mlDataY <- metaDataFixed
  mlDataX <- snmData[rownames(mlDataY),]
  
  if(trainTestFlag){
    set.seed(42)
    index <- createDataPartition(mlDataY[,"condition"], p = 0.7, list = FALSE)
    trainX <- mlDataX[index,]
    trainY <- mlDataY[index,"condition"]
    testX <- mlDataX[-index,]
    testY <- mlDataY[-index,"condition"]
    # print(testY)
    
    refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
    refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
    
    set.seed(42)
    ctrl <- trainControl(method = "repeatedcv",
                         number = numKFold,
                         repeats = numResampleIter,
                         sampling = "up",
                         summaryFunction = twoClassSummary,
                         classProbs = TRUE,
                         verboseIter = TRUE,
                         savePredictions = TRUE,
                         allowParallel=TRUE)
    
    set.seed(42)
    mlModel <- train(x = trainX,
                     y = refactoredTrainY,
                     method = modelType,
                     # preProcess = c("nzv"),
                     trControl = ctrl,
                     # tuneGrid=data.frame(n.trees=150, interaction.depth=3,
                     #                     shrinkage=0.1,
                     #                     n.minobsinnode=10),
                     metric = "ROC",
                     verbose = TRUE)
    # print(mlModel)
    
    predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[,positiveClass])
    fg <- predProbs[refactoredTestY == positiveClass]
    bg <- predProbs[refactoredTestY == negativeClass]
    
    prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
    
    multiClass <- data.frame(obs = refactoredTestY,
                             pred = predict(mlModel, newdata = testX),
                             predict(mlModel, newdata = testX, type = "prob"))
    
    print(confusionMatrix(multiClass$pred, multiClass$obs))
    print(multiClassSummary(multiClass, lev = levels(multiClass$obs)))
    
    res <- list(mlModel=mlModel,
                multiClass=multiClass)
    
  } else{
    set.seed(42)
    trainX <- mlDataX[,]
    trainY <- mlDataY[,"condition"]
    refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
    
    set.seed(42)
    ctrl <- trainControl(method = "repeatedcv",
                         number = numKFold,
                         repeats = numResampleIter,
                         sampling = "up",
                         summaryFunction = twoClassSummary,
                         classProbs = TRUE,
                         verboseIter = TRUE,
                         savePredictions = TRUE,
                         allowParallel=TRUE)
    
    set.seed(42)
    mlModel <- train(x = trainX,
                     y = refactoredTrainY,
                     method = modelType,
                     # preProcess = c("nzv"),
                     trControl = ctrl,
                     tuneGrid=data.frame(n.trees=150, interaction.depth=3,
                                         shrinkage=0.1,
                                         n.minobsinnode=1),
                     metric = "ROC",
                     verbose = TRUE)
    # print(mlModel)
    
    resPred <- mlModel$pred
    
    ## Calculate performance on concatenated fold predictions
    predProbs <- resPred
    multiClass <- resPred
    multiClass$pred <- relevel(multiClass$pred, positiveClass)
    multiClass$obs <- relevel(multiClass$obs, positiveClass)
    fg <- predProbs[resPred$obs == positiveClass,positiveClass]
    bg <- predProbs[resPred$obs == negativeClass,positiveClass]
    
    confusionMatrix <- confusionMatrix(multiClass$pred, multiClass$obs)
    prroc_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    prroc_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
    rocCurveData <- cbind(as.data.frame(prroc_roc$curve), disease_type = dt)
    prCurveData <- cbind(as.data.frame(prroc_pr$curve), disease_type = dt)
    
    ## Estimate AUROC CIs using cvAUC on concatenated fold predictions
    require(cvAUC)
    out95 <- ci.cvAUC(predictions = resPred[,positiveClass], labels = resPred$obs, label.ordering = c(negativeClass,positiveClass), folds = resPred$Resample, confidence = 0.95)
    out99 <- ci.cvAUC(predictions = resPred[,positiveClass], labels = resPred$obs, label.ordering = c(negativeClass,positiveClass), folds = resPred$Resample, confidence = 0.99)
    # resCvAUC <- data.frame(estimate = c(out95$cvAUC, out99$cvAUC), se = c(out95$se, out99$se), lowerCI = c(out95$ci[1], out99$ci[1]), upperCI = c(out95$ci[1], out99$ci[2]), levelCI = c(0.95,0.99))
    
    ## Split folds and calculate perf on each fold
    resPredSplit <- split(resPred, resPred$Resample)
    repX_perf <- list()
    for(zz in seq_along(resPredSplit)){
      resPredSingleRep <- resPredSplit[[zz]]
      predProbs <- resPredSingleRep
      multiClass <- resPredSingleRep
      multiClass$pred <- relevel(multiClass$pred, positiveClass)
      multiClass$obs <- relevel(multiClass$obs, positiveClass)
      fg <- predProbs[resPredSingleRep$obs == positiveClass,positiveClass]
      bg <- predProbs[resPredSingleRep$obs == negativeClass,positiveClass]
      
      rep_roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
      rep_pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
      
      repX_perf[[zz]] <- data.frame(auroc=rep_roc$auc,
                                    aupr=rep_pr$auc.integral,
                                    rep=paste0("Fold",zz), 
                                    diseaseType = dt)
    }
    
    # SUMMARIZE MODEL PERFORMANCES
    rep_perf <- do.call(rbind, repX_perf)
    perf <- data.frame(auroc = prroc_roc$auc,
                       aupr = prroc_pr$auc.integral,
                       aucEstimate = out95$cvAUC, # either out95 or out99 work (same result)
                       aucSE95 = out95$se,
                       lowerCI95 = out95$ci[1],
                       upperCI95 = out95$ci[2],
                       aucSE99 = out99$se,
                       lowerCI99 = out99$ci[1],
                       upperCI99 = out99$ci[2],
                       diseaseType = dt)
    
    print(perf)
    print(confusionMatrix)
    
    res <- list(mlModel=mlModel,
                multiClass=multiClass,
                rep_perf=rep_perf,
                perf=perf)
  }
  
  if(varImpFlag){
    varImpBestModelDF <- as.data.frame(varImp(mlModel$finalModel, scale = FALSE))
    varImpBestModelDF2 <- rownames_to_column(varImpBestModelDF, "Taxa")
    varImpBestModelDF2Ordered <- varImpBestModelDF2[order(-varImpBestModelDF2$Overall),]
    colnames(varImpBestModelDF2Ordered)[2] <- "varImp"
    varImpBestModelDF2OrderedNonzero <- varImpBestModelDF2Ordered[varImpBestModelDF2Ordered$varImp != 0,]
    print(head(varImpBestModelDF2OrderedNonzero))
  }
  
  plot(prroc_pr)
  plot(prroc_roc)
  
  return(res)
  
}
