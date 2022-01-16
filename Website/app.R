# RShiny_Microbe_DA_app.R
# Gregory Poore
# Nov 28, 2021
# Purpose: Create an interactive website for fungal abundances and ML perf

# # Load dependencies
library(shiny)
library(shinyauthr)
library(shinyjs)
library(shinydashboard)
library(ggpubr)
library(ggsci)
library(plotly)
library(dplyr)
library(tibble)
library(glue)
library(DT)
library(markdown) # for displaying text files

## Load data
# Decontamination information
load("data/decontamResultsV2_13Oct21.RData")
decontamResultsV2$species <- gsub("[[:space:]]","_",decontamResultsV2$species)
decontamResultsV2FormattedSpecies <- decontamResultsV2 %>% rownames_to_column("OGU") %>%
  column_to_rownames("species") %>% select(reason, OGU)
decontamResultsV2FormattedOGU <- decontamResultsV2

# Normalized abundances
load("data/data_for_ml_tcga_decontamV2_14Oct21.RData")
mergedSnmDataAndMetadataQC_Genus <- cbind(rep200FungiDecontamV2GenusVSNM, metaQiitaCombined_Nonzero_DecontamV2)
mergedSnmDataAndMetadataQC_Species <- cbind(rep200FungiDecontamV2SpeciesVSNM, metaQiitaCombined_Nonzero_DecontamV2)

#------------------------------------------------------------#
# One-time operations
autocompleteNamesGenus <- sort(colnames(rep200FungiDecontamV2GenusVSNM)) # tail(names(mergedSnmDataAndMetadataQC),-41)
autocompleteNamesSpecies <- sort(colnames(rep200FungiDecontamV2SpeciesVSNM)) # tail(names(mergedSnmDataAndMetadataQC),-41)

abbreviationTable <- read.csv(file = "./data/abbreviationTable.csv",
                              header = TRUE,
                              strip.white = TRUE,
                              stringsAsFactors = FALSE)
colnames(abbreviationTable) <- c("TCGA Abbreviation", "Cancer Type")

dataTypeList <- c("Decontaminated Batch-Corrected Genera",
                  "Decontaminated Batch-Corrected Species",
                  "Harvard Decontaminated Raw Counts (Species)",
                  "Baylor Decontaminated Raw Counts (Species)",
                  "MD Anderson Decontaminated Raw Counts (Species)",
                  "WashU Decontaminated Raw Counts (Species)",
                  "Broad Decontaminated Raw Counts (Species)",
                  "UNC Decontaminated Raw Counts (Species)",
                  "CMS Decontaminated Raw Counts (Species)")

selectedSampleTypeOrContrastList <- c("Primary Tumor",
                                      "Blood Derived Normal",
                                      "Primary Tumor vs Solid Tissue Normal")
ancombcContrastList <- c("PT", # Primary Tumor
                        "BDN", # Blood Derived Normal
                        "TvsNAT") # Tumor versus NAT

ancombcSeqCenterList <- c("Harvard Medical School",
                          "MD Anderson Institute for Applied Cancer Science",
                          "Baylor College of Medicine",
                          "Washington University School of Medicine",
                          "Broad Institute of MIT and Harvard",
                          "University of North Carolina",
                          "Canadas Michael Smith Genome Sciences Centre")

# Order of colors for differential abundance plots
# Order is primary tumor (top) > solid tissue normal (middle) > blood (bottom)
# Ref: https://www.stat.ubc.ca/~jenny/STAT545A/block17_colorsGgplot2Qualitative.html
daColors <- c("#00468BFF", "#42B540FF", "#ED0000FF")
names(daColors) <- c("Primary Tumor", "Solid Tissue Normal", "Blood Derived Normal")
#------------------------------------------------------------#

user_base <- data.frame(
  user = c("reviewer"), # knightlab
  password = c("fungi"), # cancermicrobiome
  permissions = c("standard"),
  name = c("Reviewer"),
  stringsAsFactors = FALSE,
  row.names = NULL
)

ui <- dashboardPage(
  
  skin = "purple",
  
  dashboardHeader(title = "Pan-Cancer Mycobiome",
                  tags$li(class = "dropdown", style = "padding: 8px;",
                          shinyauthr::logoutUI("logout")),
                  titleWidth = 350 #,
                  # tags$li(class = "dropdown",
                  #         tags$a(icon("database"), 
                  #                href = "ftp://ftp.microbio.me/pub/cancer_microbiome_analysis/",
                  #                title = "FTP link to download microbe data")),
                  # tags$li(class = "dropdown",
                  #         tags$a(icon("github"), 
                  #                href = "https://github.com/biocore/tcga",
                  #                title = "See the code repository on Github")
                  # )
  ),
  
  dashboardSidebar(width = 350,
                   collapsed = FALSE, 
                   div(textOutput("welcome"), 
                       style = "padding: 5px",
                       imageOutput("klLogo", width = "100%", inline = TRUE),
                       # id = "sidebar", # id important for updateTabItems
                       sidebarMenu(
                         menuItem("Interactive differential abundances (ANCOM-BC)", tabName = "home", icon = icon("chart-bar")),
                         menuItem("TCGA ML model performance and feature lists", tabName = "tab2", icon = icon("table")),
                         menuItem("Batch-corrected fungal abundances (genus)", tabName = "tab3", icon = icon("chart-bar")),
                         menuItem("Batch-corrected fungal abundances (species)", tabName = "tab4", icon = icon("chart-bar"))),
                       tableOutput("abbrev"))
  ),
  
  dashboardBody(
    # must turn shinyjs on
    shinyjs::useShinyjs(),
    # add logout button UI 
    div(class = "pull-right", shinyauthr::logoutUI(id = "logout")),
    # add login panel UI function
    shinyauthr::loginUI(id = "login"),
    
    tags$style(type="text/css",".shiny-output-error { visibility: hidden; }",".shiny-output-error:before { visibility: hidden; }"),
    tabItems(
      tabItem("home", uiOutput("ancombc")),
      tabItem("tab2", uiOutput("mlPerfUI")), 
      tabItem("tab3", uiOutput("generaAbundances")),
      tabItem("tab4", uiOutput("speciesAbundances"))
    )
  )
  
  
  
)

# Define server logic 
server <- function(input, output,session) {
  
  # call the logout module with reactive trigger to hide/show
  logout_init <- callModule(shinyauthr::logout, 
                            id = "logout", 
                            active = reactive(credentials()$user_auth))
  
  # call login module supplying data frame, user and password cols
  # and reactive trigger
  credentials <- callModule(shinyauthr::login, 
                            id = "login", 
                            data = user_base,
                            user_col = user,
                            pwd_col = password,
                            log_out = reactive(logout_init()))
  
  observe({
    if(credentials()$user_auth) {
      shinyjs::removeClass(selector = "body", class = "sidebar-collapse")
    } else {
      shinyjs::addClass(selector = "body", class = "sidebar-collapse")
    }
  })
  
  # pulls out the user information returned from login module
  user_data <- reactive({credentials()$info})
  
  user_info <- reactive({credentials()$info})
  
  output$welcome <- renderText({
    req(credentials()$user_auth)
    
    glue("Welcome {user_info()$name}")
  })
  
  output$klLogo <- renderImage({
    req(credentials()$user_auth)
    list(src = "./logos/stacked_logos_v2.jpeg",
         width = 340)
    }, deleteFile = FALSE)
  
  output$abbrev <- renderTable({
    req(credentials()$user_auth)
    abbreviationTable
    })
    
  output$generaAbundances <- renderUI({
    req(credentials()$user_auth)
    
    fluidPage(
      titlePanel("Normalized Genus-Level Fungal Abundances in TCGA"),
      
      helpText("Select fungi of interest to display its
               normalized abundance distribution across 
               TCGA cancer types."),
      
      p("Note: Fungi shown include those that passed decontaminated in TCGA (see Methods of manuscript)", style = "color:red"),
      
      selectizeInput("kraken_microbe_name_genus", "Either (i) select from the drop-down list or
                    (ii) click on the name and hit backspace, then start typing the name of your fungi of interest
                     (it will autocomplete). After selecting your microbe of interest, please wait a second
                     for the plot to appear below. A basic Kruskal-Wallis test is performed for each sample type
                     to show if the fungi varies across cancer types. The genus Malassezia 
                     is pre-selected for you.",
                     choices = autocompleteNamesGenus,
                     selected = "Malassezia",
                     width = '100%'),
      
      mainPanel(
        
      renderPlot(width = 900, height = 750, expr = {
                   
                   microbeName <- input$kraken_microbe_name_genus
                   mergedSnmDataAndMetadataQC_Genus$sample_type <- factor(mergedSnmDataAndMetadataQC_Genus$sample_type, 
                                                                    levels = c("Primary Tumor", "Solid Tissue Normal", "Blood Derived Normal"))
                   mergedSnmDataAndMetadataQC_Genus$investigation <- factor(mergedSnmDataAndMetadataQC_Genus$investigation)
                   mergedSnmDataAndMetadataQC_Genus %>%
                     select(investigation, sample_type, microbeName) %>%
                     filter(sample_type %in% c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal")) %>%
                     ggboxplot(x = "investigation", y = microbeName, 
                               color = "sample_type",
                               palette = daColors,
                               facet.by = "sample_type",
                               nrow = 3,
                               add = "jitter",
                               add.params = list(alpha=0.2),
                               xlab = "Disease Type", ylab = "Normalized Abundance (log2-cpm)", 
                               title = paste(microbeName, "\n\nAbundances Across TCGA Cancer Types"),
                               legend = "top",
                               legend.title = "Sample Type") +
                     rotate_x_text(angle = 90) +
                     theme(plot.title = element_text(hjust = 0.5), strip.text.x = element_text(size = 16)) +
                     stat_compare_means(label.x = 5, label.y.npc = 0.9)
                   
                 })
        )
      )
    })
  
  output$mlPerfUI <- renderUI({
    req(credentials()$user_auth)
    
    fluidPage(
      
      titlePanel("Machine Learning Performances Using Fungal Abundances in TCGA"),
      
      fluidRow(
        
        column(4,
               selectInput("dataTypeSelected", 
                           label = "Select Data Used for Model Training/Testing:", 
                           selected = "Decontaminated Fungal Species",
                           choices = dataTypeList)
        ),
        
        column(4,
               selectInput("selectedCancerType",
                           label = "Select Cancer Type:", 
                           selected = "Ovarian Serous Cystadenocarcinoma",
                           choices = abbreviationTable$`Cancer Type`)
        ),
        
        column(4,
               selectInput("selectedSampleTypeOrContrast",
                           label = "Select Sample Type OR Comparison:",
                           selected = "Primary Tumor",
                           choices = selectedSampleTypeOrContrastList)
        )
        
      ),
      
      helpText("Color bars show probability cutoff thresholds. 
               Performance metrics are one-cancer-type-versus-all-others unless looking at
               tumor versus normal comparisons. A constant random number seed was used to ensure
               that the models are comparable in their performance. Ten-fold cross validation was used, and
               the ROC and PR curves were generated on the concatenated predictions across all ten holdout sets."),
      
      strong(helpText("NOTE: If no plots appear, then the comparison was not modeled. This is
               due to having <20 samples in the minority class.")),
      
      p("Scroll down to the very bottom of this page to see a list of all available machine learning comparisons on 
        the raw sequencing center data subsets.", style = "color:red"),
      
      fluidRow(class = "plotRow",
        
        column(6,
               renderImage({
                 validate(
                   need(file.exists(paste0("newimages/",
                                           input$dataTypeSelected,
                                           "/",
                                           input$selectedCancerType,
                                           " -- ",
                                           input$selectedSampleTypeOrContrast,
                                           " -- ROC.png")),
                        "Comparison not available. Nothing will be plotted."
                   ))
                 
                 list(src = paste0("newimages/",
                                   input$dataTypeSelected,
                                   "/",
                                   input$selectedCancerType,
                                   " -- ",
                                   input$selectedSampleTypeOrContrast,
                                   " -- ROC.png"),
                      width = 450)
               }, deleteFile = FALSE)
        ),
        
        column(6,
               renderImage({
                 validate(
                   need(file.exists(paste0("newimages/",
                                           input$dataTypeSelected,
                                           "/",
                                           input$selectedCancerType,
                                           " -- ",
                                           input$selectedSampleTypeOrContrast,
                                           " -- PR.png")),
                        "Comparison not available. Nothing will be plotted."
                   ))
                 
                 list(src = paste0("newimages/",
                                   input$dataTypeSelected,
                                   "/",
                                   input$selectedCancerType,
                                   " -- ",
                                   input$selectedSampleTypeOrContrast,
                                   " -- PR.png"),
                      width = 450)
               }, deleteFile = FALSE)
        ),
        tags$head(tags$style(".plotRow{margin-bottom:-400px;}"))
      ),
      
      fluidRow(
        column(12, 
               align = "center",
               strong("Confusion Matrix (using a 50% probability cutoff threshold on the 30% holdout test set):", style = "color:black"),
               p("Note: A 50% probability cutoff may *not* always be the best choice for class discrimination", style = "color:red"),
               pre(renderText({
                 
                 
                 
                 cmPath <- paste0("./newconfusionmatrices/",
                                  input$dataTypeSelected,
                                  "/",
                                  input$selectedCancerType,
                                  " -- ",
                                  input$selectedSampleTypeOrContrast,
                                  " -- CM.txt")
                 
                 includeText(cmPath)
                 }))
        )
      ),
      
      h3("Ranked list of model features:"),
      
      strong(helpText("NOTE: If using species, the ranked feature importances are cross-referenced
                      against the TCGA decontamination decision when displayed in the below table.")),
      
      helpText("The Variable Importance Score is calculated based on the R gbm package."),
      
      p("Note: Feature importance scores are *not* directly compatible between models, as their calculation depends on 
        how the model was tuned during training. It is thus a *relative* feature importance score. Moreover, a high feature 
        importance score for a given taxon does *not* guarantee or imply an overabundance of that taxon for this comparison.
        A high feature importance score only means that the taxon was important for making predictions; whether the taxon is 
        less or more abundant in certain samples requires statistical testing. Lastly, the models shown in this paper were 
        *not* regularized, meaning that they could theoretically use all taxa features available in the data to make predictions. 
        All features with non-zero feature importance scores are shown here.", style = "color:red"),
      
      renderDataTable({
        
        df <- read.csv(paste0("./newfeatures/",
                              input$dataTypeSelected,
                              "/",
                              input$selectedCancerType,
                              " -- ",
                              input$selectedSampleTypeOrContrast,
                              " -- Features.csv"))
        colnames(df) <- c("Taxonomy", "Variable Importance Score")
        if(grepl("Batch-Corrected Species", input$dataTypeSelected)){
          df$reason <- decontamResultsV2FormattedSpecies[df$Taxonomy,"reason"]
          df$OGU <- decontamResultsV2FormattedSpecies[df$Taxonomy,"OGU"]
          colnames(df)[3] <- "Reason passed decontamination"
        } else if(grepl("Raw Counts", input$dataTypeSelected)){
          df$species <- decontamResultsV2FormattedOGU[df$Taxonomy,"species"]
          df$reason <- decontamResultsV2FormattedOGU[df$Taxonomy,"reason"]
          colnames(df)[3] <- "Species"
          colnames(df)[4] <- "Reason passed decontamination"
        }
        
        DT::datatable(
          df,
          filter = 'top', extensions = c('Buttons', 'Scroller'),
          options = list(scrollY = 650,
                         # scrollX = 500,
                         deferRender = TRUE,
                         scroller = TRUE))
        
      }),
      
      br(),
      h3("Full list of available machine learning comparisons when examining sequencing center subsets:"),
      
      strong(helpText("Primary tumor comparisons")),
      helpText("Harvard has the following PT comparisons: BLCA, LGG, COAD, HNSC, LUAD, PRAD, READ, STAD, THCA, UCEC"), #
      helpText("Baylor has the following PT comparisons: COAD, HNSC, KICH, KIRC, KIRP, LIHC, OV"), #
      helpText("MD Anderson has the following PT comparisons: BLCA, CESC, ESCA, HNSC, THCA, UVM"), #
      helpText("WashU has the following PT comparisons: BRCA, SARC, UCEC"), #
      helpText("Broad Insitute has the following PT comparisons: BLCA, LGG, GBM, HNSC, LUAD, LUSC, OV, PRAD, STAD, THCA"), #
      helpText("UNC has the following PT comparisons: ACC, BLCA, LGG, BRCA, CESC, CHOL, COAD, KICH, KIRC, KIRP, LIHC, 
               LUAD, LUSC, DLBC, MESO, PAAD, PCPG, PRAD, READ, SKCM, TGCT, THYM, THCA, UCS, UCEC, UVM"), #
      helpText("CMS has the following PT comparisons: ESCA, OV, STAD"), #
      strong(helpText("Blood comparisons")),
      helpText("Harvard has the following blood comparisons: BLCA, LGG, COAD, HNSC, LUAD, PRAD, READ, SKCM, STAD, THCA, UCEC"), #
      helpText("Baylor has the following blood comparisons: COAD, HNSC, KIRP, LIHC, OV"), #
      helpText("MD Anderson has the following blood comparisons: BLCA, CESC, ESCA, HNSC, THCA, UVM"), #
      helpText("WashU has the following blood comparisons: BRCA, CESC, SARC, UCEC"), #
      helpText("Broad Insitute has the following blood comparisons: BLCA, LGG, GBM, HNSC, LUAD, SKCM, STAD, THCA"), #
      strong(helpText("Tumor versus NAT comparisons")),
      helpText("Harvard has the following TvsNAT comparisons: LUAD, STAD"), #
      helpText("Baylor has the following TvsNAT comparisons: KICH, KIRC, LIHC"), #
      helpText("MD Anderson did not have sufficient sample sizes to perform ML on TvsNAT comparisons"), #
      helpText("WashU did not have sufficient sample sizes to perform ML on TvsNAT comparisons"), #
      helpText("Broad Insitute has the following TvsNAT comparisons: LUSC"), #
      helpText("UNC has the following TvsNAT comparisons: BRCA, COAD, KICH, KIRC, KIRP, LIHC, LUAD, LUSC, PRAD, THCA, UCEC"), #
      helpText("CMS has the following TvsNAT comparisons: STAD"), #
      
    )
    
  })
  
  output$speciesAbundances <- renderUI({
    req(credentials()$user_auth)
    
    fluidPage(
      titlePanel("Normalized Species-Level Fungal Abundances in TCGA"),
      
      helpText("Select fungi of interest to display its
               normalized abundance distribution across 
               TCGA cancer types."),
      
      p("Note: Fungi shown include those that passed decontaminated in TCGA (see Methods of manuscript)", style = "color:red"),
      
      selectizeInput("kraken_microbe_name_species", "Either (i) select from the drop-down list or
                    (ii) click on the name and hit backspace, then start typing the name of your fungi of interest
                     (it will autocomplete). After selecting your microbe of interest, please wait a second
                     for the plot to appear below. A basic Kruskal-Wallis test is performed for each sample type
                     to show if the fungi varies across cancer types. The species Malassezia globosa 
                     is pre-selected for you.",
                     choices = autocompleteNamesSpecies,
                     selected = "Malassezia_globosa",
                     width = '100%'),
      
      mainPanel(
        
        renderPlot(width = 900, height = 750, expr = {
          
          microbeName <- input$kraken_microbe_name_species
          mergedSnmDataAndMetadataQC_Species$sample_type <- factor(mergedSnmDataAndMetadataQC_Species$sample_type, 
                                                           levels = c("Primary Tumor", "Solid Tissue Normal", "Blood Derived Normal"))
          mergedSnmDataAndMetadataQC_Species$investigation <- factor(mergedSnmDataAndMetadataQC_Species$investigation)
          mergedSnmDataAndMetadataQC_Species %>%
            select(investigation, sample_type, microbeName) %>%
            filter(sample_type %in% c("Solid Tissue Normal", "Primary Tumor", "Blood Derived Normal")) %>%
            ggboxplot(x = "investigation", y = microbeName, 
                      color = "sample_type",
                      palette = daColors,
                      facet.by = "sample_type",
                      nrow = 3,
                      add = "jitter",
                      add.params = list(alpha=0.2),
                      xlab = "Disease Type", ylab = "Normalized Abundance (log2-cpm)", 
                      title = paste(microbeName, "\n\nAbundances Across TCGA Cancer Types"),
                      legend = "top",
                      legend.title = "Sample Type") +
            rotate_x_text(angle = 90) +
            theme(plot.title = element_text(hjust = 0.5), strip.text.x = element_text(size = 16)) +
            stat_compare_means(label.x = 5, label.y.npc = 0.9)
          
        })
      )
    )
  })
  
  output$ancombc <- renderUI({
    req(credentials()$user_auth)
    
    fluidPage(
      
      titlePanel("Interactive Volcano Plots of Differentially Abundant Fungi (ANCOM-BC)"),
      
      fluidRow(
        
        column(4,
               selectInput("seqCenterSelected_ANCOMBC", 
                           label = "Select Sequencing Center Subset:", 
                           selected = "Harvard Medical School",
                           choices = ancombcSeqCenterList)
        ),
        
        column(4,
               selectInput("selectedCancerType_ANCOMBC",
                           label = "Select Cancer Type:", 
                           selected = "Colon Adenocarcinoma",
                           choices = abbreviationTable$`Cancer Type`)
        ),
        
        column(4,
               selectInput("selectedSampleTypeOrContrast_ANCOMBC",
                           label = "Select Sample Type OR Comparison:",
                           selected = "PT",
                           choices = ancombcContrastList)
        )
        
      ),
      
      strong(helpText("NOTE #1: In primary tumor (PT) and blood (BDN) comparisons, ANCOM-BC modeled 
                      one cancer type versus all others, and right-sided features in the plot are more associated with 
                      that particular cancer type of interest. In tumor versus NAT (TvsNAT) comparisons, 
                      the right-sided features in the plot are more associated with tumor tissue. 
                      In all cases, raw decontaminated count data were used.")),
      
      strong(helpText("NOTE #2: These plots are interactive. Hover over points for more information or use
                      the options on the top right of the plot to zoom in, take a photo, etc.")),
      
      strong(helpText("NOTE #3: If no plots appear, then the comparison was not modeled. This is
               due to having <10 samples in either class after enforcing minimum library size requirements (see Methods).")),
      
      p("Scroll down to the very bottom of this page to see a list of all available comparisons.", style = "color:red"),
      
      fluidRow(
          
        renderPlotly(expr = {
          
          domain <- "fungi"
          comparison <- input$selectedSampleTypeOrContrast_ANCOMBC
          seqCenter <- input$seqCenterSelected_ANCOMBC
          cancerType <- input$selectedCancerType_ANCOMBC
          seqCenterFormatted <- gsub("[[:space:]]","",seqCenter)
          cancerTypeFormatted <- gsub("[[:space:]]","",cancerType)
          comparisonFormatted <- ifelse(grepl("PT|BDN",comparison),
                                        yes = paste0(comparison,"_1vsAll"),
                                        no = comparison)
          ancombcRes <- read.csv(paste0("ancombc/",
                                        domain,
                                        "_ancombc_",
                                        comparisonFormatted,"_",
                                        seqCenterFormatted,"_",
                                        cancerTypeFormatted,".csv"),
                                 row.names = 1)
          
          ancombcPlot <- ggplot(data = ancombcRes, 
                                aes(x=beta, y=-log10(p_val), 
                                    color=diff_abn, label=diff_label, shape=origin,
                                    text = paste0(gsub("_"," ",species)," (",OGUs,")\n",
                                                  "origin: ",origin,"\n",
                                                  "beta: ",signif(beta,3),"\n",
                                                  "se: ",signif(se,3),"\n",
                                                  "pval: ",signif(p_val,3),"\n",
                                                  "qval: ",signif(q_val,3),"\n",
                                                  "diff abundant: ",diff_abn))) + 
            geom_point(size = 1.5) +
            theme_bw() + theme(plot.title = element_text(hjust = 0.5,
                                                         margin = c(30,0,0,0))) +
            labs(x = "log-fold change\n(ANCOM-BC beta)", y = "-Log10(p-value)", 
                 color = paste0("Differentially\nabundant (q<=",0.05,") &"), 
                 title = paste0(seqCenter," | ", cancerType," | ",comparison, "\n", 
                                "<----- Less associated | More associated ----->","\n"),
                 shape = "Source reference") +
            scale_color_aaas()
          
          # ancombcPlot
          
          ggplotly(ancombcPlot, tooltip = "text", height = 600, width = 1000) %>%
            layout(margin = list(t = 100))
          
        })
      ),
      br(), br(), br(),
      br(), br(), br(),
      br(), br(), br(),
      br(),
      h3("Ranked list of ANCOM-BC features:"),
      
      p("Note: These are the outputs of ANCOM-BC modeling, ranked by q-value.
        This is a wide table, so please scroll to the right for more information.", style = "color:red"),
      
      renderDataTable({
        
        domain <- "fungi"
        comparison <- input$selectedSampleTypeOrContrast_ANCOMBC
        seqCenter <- input$seqCenterSelected_ANCOMBC
        cancerType <- input$selectedCancerType_ANCOMBC
        seqCenterFormatted <- gsub("[[:space:]]","",seqCenter)
        cancerTypeFormatted <- gsub("[[:space:]]","",cancerType)
        comparisonFormatted <- ifelse(grepl("PT|BDN",comparison),
                                      yes = paste0(comparison,"_1vsAll"),
                                      no = comparison)
        ancombcRes <- read.csv(paste0("ancombc/",
                                      domain,
                                      "_ancombc_",
                                      comparisonFormatted,"_",
                                      seqCenterFormatted,"_",
                                      cancerTypeFormatted,".csv"),
                               row.names = 1)
        
        ancombcResFormatted <- subset(ancombcRes, 
                                      select = -c(diff_name_flag,
                                                  diff_label,
                                                  OGUs)) %>%
          rownames_to_column("OGUs") %>%
          relocate(genus, species, origin, .before = beta) %>%
          relocate(q_val, diff_abn, .before = genus)
        
        DT::datatable(
          ancombcResFormatted,
          filter = 'top', extensions = c('Buttons', 'Scroller'),
          options = list(scrollY = 650,
                         scrollX = 500,
                         deferRender = TRUE,
                         scroller = TRUE))
      }),
      
      br(),
      h3("Full list of available ANCOM-BC comparisons:"),
      
      strong(helpText("Primary tumor comparisons")),
      helpText("Harvard has the following PT comparisons: BLCA, LGG, BRCA, COAD, HNSC, LUAD, PRAD, READ, SKCM, STAD, THCA, UCEC"),
      helpText("Baylor has the following PT comparisons: COAD, HNSC, KICH, KIRC, KIRP, LIHC, OV, READ"),
      helpText("MD Anderson has the following PT comparisons: BLCA, CESC, ESCA, HNSC, THCA, UVM"),
      helpText("WashU has the following PT comparisons: BRCA, CESC, SARC, UCEC"),
      helpText("Broad Insitute has the following PT comparisons: BLCA, LGG, GBM, HNSC, LUAD, LUSC, OV, PRAD, STAD, THCA"),
      helpText("UNC has the following PT comparisons: BLCA, LGG, BRCA, CESC, COAD, KICH, KIRC, KIRP, LIHC, 
               LUAD, LUSC, DLBC, PAAD, PRAD, READ, SKCM, THYM, THCA, UCS, UCEC"),
      helpText("CMS has the following PT comparisons: ESCA, OV, STAD"),
      strong(helpText("Blood comparisons")),
      helpText("Harvard has the following blood comparisons: BLCA, LGG, BRCA, COAD, HNSC, LUAD, PRAD, READ, SKCM, STAD, THCA, UCEC"),
      helpText("Baylor has the following blood comparisons: COAD, HNSC, KIRP, LIHC, OV, READ"),
      helpText("MD Anderson has the following blood comparisons: BLCA, CESC, ESCA, HNSC, THCA, UVM"),
      helpText("WashU has the following blood comparisons: BRCA, CESC, SARC, UCEC"),
      helpText("Broad Insitute has the following blood comparisons: BLCA, LGG, GBM, HNSC, LUAD, LUSC, OV, PRAD, SKCM, STAD, THCA"),
      strong(helpText("Tumor versus NAT comparisons")),
      helpText("Harvard has the following TvsNAT comparisons: BLCA, COAD, HNSC, LUAD, PRAD, STAD"),
      helpText("Baylor has the following TvsNAT comparisons: COAD, KICH, KIRC, LIHC"),
      helpText("MD Anderson has the following TvsNAT comparisons: ESCA"),
      helpText("WashU has the following TvsNAT comparisons: BRCA"),
      helpText("Broad Insitute has the following TvsNAT comparisons: LUAD, LUSC"),
      helpText("UNC has the following TvsNAT comparisons: BLCA, BRCA, COAD, KICH, KIRC, KIRP, LIHC, LUAD, LUSC, PRAD, THCA, UCEC"),
      helpText("CMS has the following TvsNAT comparisons: ESCA, STAD"),
      
    )
    
  })
  
}

shinyApp(ui, server)