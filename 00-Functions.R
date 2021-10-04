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
  plot_bar(ps2, fill="Domain") + scale_fill_manual(values = c("#0072B5FF","#BC3C29FF")) +
    ylab(yAxisLab) + xlab("TCGA Cancer Type") + 
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
    ggtitle(paste0(title," (Primary Tumor | ",dataType,")")) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0(figFilepath,filename,".png"),
         dpi = "retina", width = 12, height = 4, units = "in")
  ps2_otu <- data.frame(otu_table(ps2))
  ps2_otu %>% write.csv(file = paste0(dataFilepath,filename,".csv"))
  ps2Df <- data.frame(otu_table(ps2))
  print("Normalized relative abundance means (%):")
  print(round(100*colMeans(ps2Df),2))
  return(ps2)
}

vsnmFunctionTCGA <- function(qcData, qcMetadata=metaQiitaCombined_Nonzero, cancerTypeFlag=FALSE, filename){
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

findPrevGenus <- function(metaData=metaQiitaCombined_Nonzero_8cancer_shared, 
                     countData=rep200FungiGenusShared_8cancer, 
                     sampleType=c("Primary Tumor","Solid Tissue Normal"), 
                     cancerType=c("Breast Cancer"), seqCenter=c("Harvard Medical School"),
                     taxaFlag=FALSE, weizmannIntersectFlag=FALSE, weizmannTissue=c("breast"),
                     weizmannTissueType=c("tumor","nat"),
                     weizmannMeta=weizmannMetaGenusShared, weizmannData=weizmannGenusShared){
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

findPrevSpecies <- function(metaData=metaQiitaCombined_Nonzero_8cancer_shared, 
                     countData=rep200FungiSpeciesShared_8cancer, 
                     sampleType=c("Primary Tumor","Solid Tissue Normal"), 
                     cancerType=c("Breast Cancer"), seqCenter=c("Harvard Medical School"),
                     taxaFlag=FALSE, weizmannIntersectFlag=FALSE, weizmannTissue=c("breast"),
                     weizmannTissueType=c("tumor","nat"),
                     weizmannMeta=weizmannMetaSpeciesShared, weizmannData=weizmannSpeciesShared){
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
              brayData=brayData,
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
                       stageNum = "stageI",
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
                                                        cristianoRep200FungiSpeciesShared),
                                     dataStringList = c("full_rep200",
                                                        "bacteria_only",
                                                        "fungi_only",
                                                        "fungi_decontam",
                                                        "fungi_intersected_with_Weizmann"),
                                     col2Predict = "one_cancer_vs_healthy",
                                     dzOfInterestList = list("Bile Duct Cancer",
                                                              "Breast Cancer",
                                                              "Colorectal Cancer",
                                                              "Gastric cancer",
                                                              "Lung Cancer",
                                                              "Ovarian Cancer",
                                                              "Pancreatic Cancer"),
                                     stageNum = "stageI",
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
                                                        cristianoRep200FungiSpeciesShared),
                                     dataStringList = c("full_rep200",
                                                        "bacteria_only",
                                                        "fungi_only",
                                                        "fungi_decontam",
                                                        "fungi_intersected_with_Weizmann"),
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

ml1VsAllUCSD <- function(metaData,
                              snmData,
                              col2Predict = "HvsC",
                              dzOfInterest = "NSCLC",
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
  }
  
  mlDataY <- metaDataFixed
  mlDataX <- snmData[rownames(mlDataY),]
  
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
