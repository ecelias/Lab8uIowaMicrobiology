#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(BiocManager)
library(shiny)
library(tidyverse)
library(ggplot2)
library(vegan)
library(data.table)
library(readr)
# installing phyloseq requires Bioconductor and installation of the "igraph" package beforehand
library(phyloseq)
library(phyloseqCompanion)
library(broom)
library(shinycssloaders)
library(bslib)

# define helper functions to be utilized in front and back end operations
# function to take the level 5 file and separate for future use
separate_taxa <- function(level5, meta){
  numSamples <<- ncol(level5) - 1
  # sets the first column of the "level 5" taxonomy file to the index
  colnames(level5)[1] <- "index" 
  # separates the data within the level 5 file using ";" as the delimiter
  level_5_sep <- separate(level5, "index", into=c("Kingdom", "Phylum", "Class", "Order", "Family"), sep=";")
  
  # dataset of samples with sample name, abundance and 1% threshold
  # uses "na.rm=TRUE to remove any null values"
  samples <- level5 %>% summarize_if(is.numeric, sum, na.rm = TRUE)
  samples <- pivot_longer(samples, 1:ncol(samples), names_to = "Sample", values_to = "Abundance")
  
  samples$threshold <- samples$Abundance * 0.01
  samples$threshold_rare <- min(samples$Abundance) * 0.01
  
  # dataset of treatments and how they relate to sample columns
  treatments <- as.data.frame(unique(meta[,2]))
  colnames(treatments) <- c("treatment")

  for (level5 in 1:nrow(treatments)) {
    treatment_rows <- which(meta$treatment == treatments[level5, 1])
    treatments[level5, 2] <- min(treatment_rows) + 1  
    treatments[level5, 3] <- max(treatment_rows) + 1
  }
  
  colnames(treatments) <- c("treatment", "min_col", "max_col")
  
  # create different datasets based on taxonomy
  phylum <- cbind(paste(level_5_sep$Kingdom, level_5_sep$Phylum), level_5_sep[,6:ncol(level_5_sep)])
  colnames(phylum)[1] <- "Phylum"
  
  class <- cbind(paste(level_5_sep$Kingdom, level_5_sep$Phylum, level_5_sep$Class), level_5_sep[,6:ncol(level_5_sep)])
  colnames(class)[1] <- "Class"
  
  order <- cbind(paste(level_5_sep$Kingdom, level_5_sep$Phylum, level_5_sep$Class, level_5_sep$Order), level_5_sep[,6:ncol(level_5_sep)])
  colnames(order)[1] <- "Order"
  
  family <- cbind(paste(level_5_sep$Kingdom, level_5_sep$Phylum, level_5_sep$Class, level_5_sep$Order, level_5_sep$Family), level_5_sep[,6:ncol(level_5_sep)])
  colnames(family)[1] <- "Family"
  
  # create list with separate taxa files, samples files, and treatment files
  # to access an item in this list: sep_taxa[[ITEM_NAME]] --> i.e sep_taxa[[phylum]]
  sep_taxa <- list(phylum=phylum, class=class, order=order, family=family, treatments=treatments, samples=samples)
  return(sep_taxa)
}

# function to create "other" category and reorder
create_other <- function(longdata, rarified) {
  
  # debugging statement
  # if(rarified == TRUE){
  #  print(paste("Dimensions of longdata:", nrow(longdata), "rows,", ncol(longdata), "columns"))
  # }
  
  # long data with "other" category
  otherdata <- longdata
  otherdata <- otherdata %>% mutate_if(is.factor, as.character)
  
  otherrows = nrow(longdata)
  
  # create "Other" category
  for(i in 1:numSamples)
    for(j in seq(i, otherrows, by=numSamples))
      # rarified data
      if(rarified == TRUE){
        if(otherdata$Abundance[j]<samples$threshold_rare[i])
          otherdata[j,1] <- "Other"
      }
  # raw data
  else {
    if(otherdata$Abundance[j]<samples$threshold[i])
      otherdata[j,1] <- "Other"
  }
  
  # convert taxa back to a factor
  otherdata <- otherdata %>% mutate_if(is.character, as.factor)
  otherdata <- as.data.frame(otherdata)
  
  # convert to wide format to sort taxa by overall abundance
  widedata <- otherdata
  widedata <- widedata %>% pivot_wider(names_from = "Sample", values_from = "Abundance", values_fill = list(Abundance = 0), values_fn = list(Abundance = sum))
  
  # calculate overall abundance and reorder taxa
  widedata$abundance <- rowSums(widedata[,-1], na.rm = TRUE)
  widedata[[1]] <- reorder(widedata[[1]], -widedata$abundance)
  
  # delete overall abundance and convert back to long
  widedata <- widedata[, -ncol(widedata)]
  newdata <- widedata %>% pivot_longer(cols = -1, names_to = "Sample", values_to = "Abundance")
  
  # add treatment data from meta data
  newdata <- merge(newdata, meta_global, by.x = "Sample", by.y = 1)
  names(newdata)[ncol(newdata)] <- "treatment"
  
  return(newdata)
}

# makes datasets to use in future analysis
create_datasets <- function(sep_taxa, taxa) {
  
  # select taxonomic dataset
  x <- as.data.frame(sep_taxa[[taxa]])
  
  # count the number of samples
  numSamples = nrow(sep_taxa$samples)
  
  numeric_cols <- sapply(x, is.numeric)
  
  # count the number of samples in columns (this will give the abundance of each sample)
  #columnData <- x %>% group_by_at(1) %>% summarize_if(is.numeric, sum, na.rm = TRUE)
  columnData <- x %>%
    group_by_at(1) %>%
    summarize(across(where(is.numeric), sum, na.rm = TRUE))
  #reorder column data by overall abundance
  columnData$abundance<-rowSums(columnData[,-1],na.rm = TRUE)
  columnData[[1]] <- factor(columnData[[1]], levels = columnData[[1]][order(columnData$abundance, decreasing = TRUE)])
  # columnData[[1]] <- reorder(columnData[[1]], columnData$abundance)
  columnData <- columnData[,-ncol(columnData)]
  
  # long data
  longData <- columnData %>% pivot_longer(cols = -1, names_to = "Sample", values_to = "Abundance")
  
  # long data with other (raw) category
  longDataOther <- create_other(longData, FALSE)
  
  # create datasets to use with phyloseq
  otuMatrix <- data.matrix(subset(columnData, select=c(2:ncol(columnData))))
  rownames(otuMatrix) <- paste0("sp", 1:nrow(otuMatrix))
  
  taxaMatrix <- as.matrix(subset(columnData, select=c(1)))
  rownames(taxaMatrix) <- paste0("sp", 1:nrow(taxaMatrix))
  
  OTU = otu_table(otuMatrix, taxa_are_rows=TRUE)
  TAX = tax_table(taxaMatrix)
  OTU_t <- t(OTU)
  physeqSamples <- meta_global
  rownames(physeqSamples) <- meta_global$sample
  physeq = phyloseq(OTU, TAX, sample_data(as.data.frame(physeqSamples)))
  
  # create rarefy datasets
  physeqRare <- rarefy_even_depth(physeq, rngseed = 10, sample.size = min(sample_sums(physeq)), replace = TRUE)
  otuRarefy <- otu_table(physeqRare)
  otuRarefy_t <- t(otuRarefy)
  taxaRarefy <- tax_table(physeqRare)
  
  # Align OTU identifiers between taxaRarefy and otuRarefy
  common_otus <- intersect(rownames(taxaRarefy), colnames(otuRarefy_t))

  if (length(common_otus) == 0) {
    stop("No common OTUs found between taxaRarefy and otuRarefy. Check your data.")
  }
  # Subset taxaRarefy and otuRarefy to include only common OTUs
  taxaRarefy <- taxaRarefy[common_otus, ]
  otuRarefy_t <- otuRarefy_t[, common_otus]
  
  # Convert taxaRarefy to a data frame
  taxaRarefy_df <- as.data.frame(taxaRarefy)
  rownames(taxaRarefy_df) <- common_otus
  
  # Convert otuRarefy_t to a data frame
  otuRarefy_t_df <- as.data.frame(otuRarefy_t)
  colnames(otuRarefy_t_df) <- common_otus
  
  # Transpose otuRarefy_t_df to have OTUs as rows and samples as columns
  otuRarefy_t_df <- t(otuRarefy_t_df)
  
  # Merge taxaRarefy_df and otuRarefy_t_df
  colDataRare <- merge(taxaRarefy_df, otuRarefy_t_df, by = "row.names")
  
  # Debugging: Print the result of the merge
  # print("Result of merge:")
  # print(colDataRare)
  
  # Delete row names
  colDataRare <- colDataRare[-1]
  
  # reorder column data by overall abundance
  colDataRare$abundance <- rowSums(colDataRare[,-1], na.rm = TRUE)
  colDataRare[[1]] <- reorder(colDataRare[[1]], colDataRare$abundance)
  colDataRare <- colDataRare[,-ncol(colDataRare)]
  longDataRare <- colDataRare %>% pivot_longer(cols = -1, names_to = "Sample", values_to = "Abundance")
  
  # long data (rarified) with other
  longDataRareOther <- create_other(longDataRare, TRUE)
  
  # create dataset with diversity indices
  richness <- specnumber(OTU_t)
  shannon <- diversity(OTU_t, index = "shannon")
  simpson <- diversity(OTU_t, index = "simpson")
  diversityResults <- cbind(meta_global[,2:ncol(meta_global)], richness, shannon, simpson)
  
  colnames(diversityResults) <- c("treatment", "richness", "shannon", "simpson")
  
  # create dataset with diversity indices based on rarefy data
  richnessRarefy <- specnumber(otuRarefy_t)
  shannonRarefy <- diversity(otuRarefy_t, index = "shannon")
  simpsonRarefy <- diversity(otuRarefy_t, index = "simpson")
  diversityResultsRarefy <- cbind(meta_global[,2:ncol(meta_global)], richnessRarefy, shannonRarefy, simpsonRarefy)
  
  # rename columns of rarefied diversity results to be consisted with raw diversity results
  colnames(diversityResultsRarefy) <- c("treatment", "richness", "shannon", "simpson")
  
  # create a list with datasets
  datasetList <- list(columnData = columnData, longData = longData, longDataOther = longDataOther, OTU = OTU, OTU_t = OTU_t, physeq = physeq,
                      otuRarefy = otuRarefy, otuRarefy_t = otuRarefy_t, physeqRare = physeqRare, colDataRare = colDataRare, 
                      longDataRare = longDataRare, longDataRareOther = longDataRareOther, diversityResults = diversityResults,
                      diversityResultsRarefy = diversityResultsRarefy, taxa = TAX)
  
  return(datasetList)
  
}

# function to find core taxa found in all treatments and samples 
# parameters: columnData
# returns "core" variable which contains core taxa
core_taxa <- function(colData) {
  core <- colData %>% filter_if(is.numeric, all_vars(.>0))
  core$abundance <- rowSums(core[,-1], na.rm = TRUE)
  core <- core[order(-core$abundance)]
  return(core)
}

# create a title for Core Taxa table using metadata, core, and columnData
# parameters: dataframe (columnData) and core variable
# returns: a string with the defined caption
core_caption <- function(df, core) {
  treatment <- paste(unique(meta_global$treatment), collapse = "/")
  numCore = nrow(core)
  total = nrow(df)
  caption = paste("<p><b>Core taxa in", treatment, "treatments:", numCore, "of", total, "taxa are shared between treatments</b><p>")
  return(toString(caption))
}

# function to find taxa unique to a treatment
# parameters: columnData, treatments table
# returns: list of unique treatments
unique_taxa <- function(level5, treatments){
  unique_list <- list()
  
  for(i in 1:nrow(treatments)) {
    unique_sample <- level5 %>% select(!(treatments[i, "min_col"]:treatments[i, "max_col"])) %>% filter_if(is.numeric, all_vars(.==0))
    treatment <- as_tibble(treatments)
    
    if(nrow(unique_sample) > 0){
      unique_list[[i]] <- cbind(treatment[i,1], unique_sample[,1])
      
      # add the names to the objects in the "uniqueList" to allow creation of other tables later
      name <- treatment[i, 1]
      names(unique_list)[i] <- name
    }
  }
  return(unique_list)
}

# function to graph rarefaction curve with ggplot2
graphRare <- function(x, ylab, bytype) {
  rareSample = list()
  
  # Write values for each sample from a rarecurve function into a list
  # Sample names obtained from "samples" dataframe
  for(i in 1:numSamples) {
    # rareSample[[i]] <- cbind(samples[i,1], as.data.frame(x[[i]]), attributes((x[[i]])$Subsample))
    rareSample[[i]] <- cbind(samples[i,1], as.data.frame(x[[i]]), attr(x[[i]], "Subsample"))
  }
  
  # Bind the data together from the "rareSample" list
  rareGraph <- bind_rows(rareSample)
  rareGraph <- merge(rareGraph, meta_global, by.x = "Sample", by.y = 1)
  colnames(rareGraph) <- c("Sample", "num_taxa", "num_samples", "Treatment")
  ggplot(rareGraph, aes(x = num_samples, y = num_taxa, group = Sample, color = .data[[bytype]])) +
    labs(x = "Sample Size", y = ylab) + theme_bw() + geom_line(size = 1) + guides(fill = guide_legend(title = bytype))
}

# create options for drop down menus used in UI
# drop down choices will be used to determine visualizations selected in certain tabs
taxachoices <- list("Phylum", "Class", "Order", "Family")
raw_or_rare <- list("Raw Data", "Rarified Data")
abs_or_rel <- list("Absolute Abundance", "Relative Abundance")
sample_or_treatment <- list("Sample", "Treatment")
diversity_choices <- list("Richness"="richness", "Shannon Diversity"="shannon", "Simpson Diversity"="simpson")
distance_choices <- list("Jaccard" = "jaccard", "Bray-Curtis" = "bray", "Morisita-Horn" = "horn")
ordination_choices <- list("NMDS", "PCoA")

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Bean Beetle Microbiome Analysis"), 
  tabsetPanel(id = "tabs", 
              tabPanel(value = "tab1", title = "Home", 
                       mainPanel(
                         h3("Welcome to the Bean Beetle Microbiome Analysis App"),
                         p(""),
                         h4("The Bean Beetle Microbiome Project is a research/teaching collaboration of institutions across the US that is studying the microbiome of", em("Callosobruchus maculatus"), "in research experiences (CUREs)."),
                         p("This app is designed to lead students through the community analysis of level-5 (family level) datasets produced in DNA subway"),
                         p("Before proceeding, students should ensure the level 5 file is formatted correctly such that all unidentified taxa, chloroplasts, and mitochondria have been removed. The level 5 file should also be formatted with the first column as taxa and the subsequent columns as samples with unique sample identifiers."),
                         p("Students should also prepare a metadata file in the first column as samples and the second column as treatments.", strong("Both files should be in .csv format.")), 
                         p("Additionally, this app should be capable of community analysis for any level 5 data"),
                         p(""),
                         
                         p("This app was reconfigured for the University of Iowa MICR 2158 course based on the source materials from Huang et al., 2022 by Elizabeth Elias, an undergraduate student at the University of Iowa under the guidance of Dr. Regina McGrane, Department of Microbiology and Immunology, University of Iowa"),
                         p(tags$a("The original app can be found by clicking here.", href = "https://beanbeetles.shinyapps.io/BeanBeetleMicrobiome/")),
                         p(tags$a("For more information on this CURE project, please click here", href = "https://www.beanbeetles.org/microbiome/the-bean-beetle-microbiome-project/")),
                         p(tags$a("Click here to find the GitHub repo for this project", href = "https://github.com/ecelias/Lab8uIowaMicrobiology"))
                       )
              ),
              
              tabPanel(value = "tab2", title = "Data Upload", 
                       sidebarLayout(
                         sidebarPanel(
                           fileInput("file1", "Choose level 5 CSV file", 
                                     accept = c(
                                       "text/csv",
                                       "text/comma-separated-values, text/plain",
                                       ".csv")),
                           tags$hr(),
                           fileInput("file2", "Choose metadata CSV file", 
                                     accept = c(
                                       "text/csv",
                                       "text/comma-separated-values, text/plain",
                                       ".csv")),
                           tags$hr(),
                           card(
                             "Click run app after selecting files",
                             actionButton("run", "Run App")
                           )
                         ), 
                         mainPanel(
                           layout_column_wrap(
                             width = 1,
                             height = 1000,
                             card(
                               width = 1/2, 
                               height = 300,
                               layout_columns(
                                 card(
                                   height = 350,
                                   card_header(
                                     class = "bg-primary mb-3",
                                     "Level 5 File"
                                   ),
                                   card_body(
                                     tags$i("Ensure the first column is the combined taxa (including kingdom) separated by semi-colons. The remaining columns should contain the abundance data for each sample. The 6 most abundant taxa will be displayed below."),
                                     tags$hr(),
                                     tableOutput("level_5_contents"),
                                     tags$hr()
                                   )
                                 ),
                                 card(
                                   height = 350,
                                   card_header(
                                     class = "bg-primary mb-3",
                                     "Metadata File"
                                   ),
                                   card_body(
                                     tags$i("Ensure the first column lists the samples and the second column lists the treatments for each sample. Your metadata is displayed below."),
                                     tags$hr(),
                                     tableOutput("metadata_contents"),
                                     tags$hr()
                                   )
                                 )
                               )
                             ),
                             card(
                               height = 650,
                               card_header(
                                 class = "bg-secondary mb-3",
                                 "Rank Abundance Curve"
                               ),
                               card_body(
                                 
                                 plotOutput("rankabundancecurve", width = "100%", height = "auto")
                                 
                               )
                             )
                           )
                         )
                       )
              ),
              
              # Tab to display Core Taxa present in samples
              tabPanel(value = "tab3", title = "Core Taxa",
                       sidebarLayout(
                         sidebarPanel(
                           card(
                             ## is there a way to make it so they can quickly click through all 4 taxa levels? instead of dropdown menu
                             selectInput("taxon_core", "Select a taxon", taxachoices), width = 2
                           )
                         ),
                         mainPanel(
                           card(
                             h2("Core Taxa"),
                             p("Core taxa are those taxa found in all samples."),
                             verticalLayout(htmlOutput("coreCaption"), tableOutput("coreTaxa")), 
                             width = 10
                           )
                         )
                       )
                       
              ), 
              tabPanel(value = "tab4", title = "Unique Taxa", 
                       sidebarLayout(
                         sidebarPanel(
                           card(
                             ## is there a way to make it so they can quickly click through all 4 taxa levels? instead of dropdown menu
                             selectInput("taxon_unique", "Select a taxon", taxachoices), width = 2
                           )
                         ),
                         mainPanel(
                           card(
                             h2("Unique Taxa"),
                             p("Unique taxa are those taxa found in a single treatment."),
                             textOutput('no_unique_captions'),
                             uiOutput('unique_tables')
                           )
                         )
                       )
                       
              ),
              tabPanel(value = "tab5", title = "Rarefaction", 
                       sidebarLayout(
                         sidebarPanel(
                           selectInput('taxon_rarefy', 'Select a taxonomic level', taxachoices),
                           selectInput('by_type', "Graph by:", sample_or_treatment)
                         ),
                         mainPanel(
                             card(
                               h2("Sample Rarefaction Curves"), 
                               h3("Raw Data"),
                               plotOutput('initial_rarefaction', width='100%', height='400px') %>%
                                 withSpinner(color='#0dc5c1')
                             ),
                             card(
                               h3('Even rarefaction to minimum number of sequences'),
                               plotOutput('even_rarefaction', width='100%', height='400px') %>%
                                 withSpinner(color='#0dc5c1')
                             )
                         )
                       )
              ),
              tabPanel(value = "tab6", title = "Taxonomy Bar Graphs", fluid = TRUE,
                       sidebarLayout(
                         sidebarPanel(
                           selectInput('taxon_bar', 'Select a taxonomic level', taxachoices),
                           selectInput('rawrare_bar', 'Select which data to use', raw_or_rare),
                           selectInput('abs_rel', 'Graph by:', abs_or_rel), width = 4
                         ),
                         mainPanel(
                           plotOutput('bargraph',width='100%', height='auto')
                         )
                       )
              ),
              tabPanel(value = "tab7", title = "Taxonomy Heat Map", 
                       sidebarLayout(
                         sidebarPanel(
                           selectInput('taxon_heatmap', 'Select a taxonomic level', taxachoices),
                           selectInput('rawrare_heatmap', 'Select which data to use', raw_or_rare),
                           selectInput('abs_rel', 'Graph by:', abs_or_rel), width = 4
                         ),
                         mainPanel(
                           plotOutput('heatmap', height='auto')
                         )
                       )
              ),
              tabPanel(value = "tab8", title = "Alpha Diversity", 
                       sidebarLayout(
                         sidebarPanel(
                           selectInput('taxon_alpha_test', 'Select a taxonomic level', taxachoices),
                           selectInput('rawrare_alpha', 'Select which data to use', raw_or_rare),
                           selectInput('div_measure', 'Select diversity measure', diversity_choices)
                         ),
                         mainPanel(
                           h3('Alpha Diversity'),
                           plotOutput('alpha_plots'), 
                           tags$hr(), 
                           htmlOutput('alpha_caption'),
                           textOutput('anova_caption'),
                           tableOutput('alpha_anova'),
                           textOutput('posthoc_caption'),
                           tableOutput('alpha_posthoc')
                         )
                       )
              ),
              tabPanel(value = "tab9", title = "Beta Diversity", 
                       sidebarLayout(
                         sidebarPanel(
                           selectInput('taxon_beta_test', 'Select a taxonomic level', taxachoices),
                           selectInput('rawrare_beta', 'Select which data to use', raw_or_rare),
                           selectInput('dist_measure', 'Select distance measure', distance_choices),
                           selectInput('ord_method', 'Select ordination method', ordination_choices),
                           selectInput('samptreat', 'Graph by:', sample_or_treatment)
                           
                         ),
                         mainPanel(
                           h3('Beta Diversity'),
                           plotOutput('ordination_plot'), 
                           textOutput('ordination_caption'),
                           tags$hr(), 
                           textOutput('permanova_caption'),
                           tableOutput('beta_permanova')
                         )
                       )
                       
              ),
              tabPanel(value = 'tab10', title='Helpful Info!')
  ),
  theme = bs_theme(preset = "slate")
  # close UI
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # variables for level_5 and metadata
  level_5 <- 0
  metadata <- 0
  
  observeEvent(input$tabs,{
    
    # if the user goes to tab 2
    # functions to accept the csv files and return a data frame for level 5 and metadata
    # displays the tables the user inputs
    # will also generate a rank abundance curve 
    if(input$tabs == "tab2"){
      l5_upload <- reactive({
        inFile1 <- input$file1
        if (is.null(inFile1))
          return(NULL)
        read_csv(inFile1$datapath)
      })
      meta_upload <- reactive({
        inFile2 <- input$file2
        if (is.null(inFile2))
          return(NULL)
        read_csv(inFile2$datapath)
      })
      # upload level 5 and metadata files
      # convert to table format for cards
      observeEvent(input$run, {
        level_5 <- l5_upload()
        level_5 %>% filter_all(any_vars(. != 0))
        metadata <- meta_upload()
        output$level_5_contents <- renderTable({head(level_5)})
        output$metadata_contents <- renderTable({metadata})
        
        # set column names in metadata to sample and treatment
        colnames(metadata) <- c("sample", "treatment")
        
        sep_taxa <- separate_taxa(level_5, metadata)
        
        # create a global dataset for treatments and samples, indicated by double arrows
        samples <<- sep_taxa$samples
        treatments <<- sep_taxa$treatments
        meta_global <<- as.data.frame(metadata)
        
        # create global datasets for each taxonomix level
        Phylum <<- create_datasets(sep_taxa, "phylum")
        Class <<- create_datasets(sep_taxa, "class")
        Order <<- create_datasets(sep_taxa, "order")
        Family <<- create_datasets(sep_taxa, "family")
      })
    }
    # server side functions for core taxa visualization
    else if(input$tabs == 'tab3'){
      if(exists('Phylum')){
        
        # creates core dataset based on taxonomic level selected for tab4
        core_data <- reactive ({
          taxa <- input$taxon_core
          my_list <- get(taxa)
          core <- core_taxa(my_list$columnData)
          return(core)
        })
        
        # creates the caption for the core dataset based on the taxonomic level selected
        core_data_caption <- reactive({
          taxa <- input$taxon_core
          my_list <- get(taxa)
          core <- core_taxa(my_list$columnData)
          caption <- core_caption(my_list$columnData,core)
          return(caption)
        })
        # returns caption and core taxa table to output
        output$coreCaption <- renderUI({HTML(core_data_caption())})
        output$coreTaxa <- renderTable({core_data()},digits=0)
      }
    }
    # server side functions for unique taxa visualization
    else if (input$tabs == 'tab4'){
      if(exists('Phylum')){
        
        # create unique taxa output for tab 4 
        # creates unique taxa datset based on taxonomic level selected
        # Create a list with different objects for each treatment
        unique_data <- reactive({
          taxa <- input$taxon_unique
          my_list <- get(taxa)
          unique_taxa_list <- list()
          unique_taxa_list <- unique_taxa(my_list$columnData, treatments)
          
          return(unique_taxa_list)
        })
        
        # creates a list of the number of output tables based on the number of treatments 
        # tables are labeled with treatment names
        output$unique_tables <-
          renderUI({
            table_output_list <- lapply(treatments$treatment, function(i){
              tablename <- paste0("table", i)
              htmlOutput(tablename)
            })
            # generates a list of HTML tags Shiny uses to categorize each table output
            # allows Shiny to properly render the unique tables
            tagList(table_output_list)
          })
        # render the table
        # use the "observe" command to allow reactive functions to run
        observe({
          unique_data_tables <- unique_data()
          
          if (!length(unique_data_tables))
          {
            output$no_unique_caption <- renderText({"There are no unique taxa."})
          }
          
          for (i in 1:nrow(treatments)) {
            local({
              my_i <- treatments[i,1]
              # create the caption
              num_unique <- length(unique_data_tables[[my_i]]$treatment)
              taxa <- input$taxon_unique
              my_list <- get(taxa)
              total <- nrow(my_list$columnData)
              caption_title <- paste('Unique taxa in', my_i, 'treatment:', num_unique, "of", total, "taxa are unique to this treatment")
              
              tablename <- paste0("table", my_i)
              output[[tablename]] <- renderTable(
                {
                  unique_data_tables[[my_i]]
                },
                caption = caption_title, caption.placement = getOption(
                  'xtable.caption.placement', 'top'
                )
              )
            })
          }
        })
      }
      
    }
    # server side functions for rarefaction visualization
    else if (input$tabs == 'tab5'){
      if(exists('Phylum')){
        # select taxa level data and create the initial rarefaction graph for tab 5
        which_initial_rarefaction <- reactive({
          taxa <- input$taxon_rarefy
          by_type <- input$by_type
          my_list <- get(taxa)
          my_data <- my_list$OTU_t
          
          # Coerce the transposed OTU table into a matrix. rarecurve() only 
          # accepts matrix-like objects and will not accept an OTU table
          my_data <- otu.matrix(my_data)
          rarecurveData <- rarecurve(my_data, step=20, sample=20, xlab="Sample Size", ylab="Species", label=FALSE, tidy=FALSE)
          graphRare(rarecurveData, taxa, by_type)
        })
        output$initial_rarefaction <- renderPlot({which_initial_rarefaction()})
        
        # select data and create even rarefaction graph
        which_even_rarefaction <- reactive ({
          taxa <- input$taxon_rarefy
          by_type <- input$by_type
          my_list <- get(taxa)
          my_data <- my_list$otuRarefy_t
          my_data <- otu.matrix(my_data)
          rarecurveData <- rarecurve(my_data, step=20, sample=20, xlab="Sample Size", ylab="Species", label=FALSE, tidy=FALSE)
          graphRare(rarecurveData, taxa, by_type)
        })
        output$even_rarefaction <- renderPlot({which_even_rarefaction()})
      }
      
    }
    # server side functions for bar graphs
    else if (input$tabs == 'tab6'){
      if(exists('Phylum')) {
        which_bar_graph <- reactive({
          taxa <- input$taxon_bar
          my_list <- get(taxa)
          if(input$rawrare_bar=='Raw Data'){
            my_data <- my_list$longDataOther
          }
          else {
            my_data <- my_list$longDataRareOther
          }
          if(input$abs_rel == 'Absolute Abundance'){ # absolute abundance
            ggplot(data=my_data, aes(x=Sample, y=Abundance)) +
              geom_bar(aes(fill=my_data[,2]), position='stack', stat='identity')+
              labs(fill=taxa, y='Absolute Abundance')+
              facet_grid(.~treatment, space='free_x', scales='free_x')+
              theme(legend.position='bottom')+
              guides(fill=guide_legend(ncol=2))
          }
          else{
            ggplot(data=my_data, aes(x=Sample, y=Abundance))+
              geom_bar(aes(fill=my_data[,2]), position='fill', stat='identity')+
              labs(fill=taxa, y='Relative Abundance')+
              facet_grid(.~treatment, space='free_x', scales='free_x')+
              theme(legend.position='bottom')+
              guides(fill=guide_legend(ncol=2))
          }
        })
        # determine height of the taxonomy bar graphs in tab3
        bar_graph_height <- reactive({
          taxa <- input$taxon_bar
          my_list <- get(taxa)
          my_data <- my_list$longDataOther
          num_taxa <- length(unique(my_data[,2]))
          height = 300 + num_taxa*10
          return(height)
        })
        # output the taxonomy bar graph 
        observe({output$bargraph <- renderPlot({which_bar_graph()}, height = bar_graph_height())})
      }
    }
    # server side functions for taxa heatmaps
    else if (input$tabs == 'tab7'){
      heatmap_height <- reactive({
        taxa <- input$taxon_heatmap
        my_list <- get(taxa)
        my_data <- my_list$longDataOther
        num_taxa <- nrow(unique(my_data[,1]))
        height <- max(c(225, num_taxa*15))
        return(height)
      })
      heatmap <- reactive({
        taxa <- input$taxon_heatmap
        my_list <- get(taxa)
        if(input$rawrare_heatmap=='Raw Data'){
          my_data <- my_list$longData
        }
        else{ 
          my_data <- my_list$longDataRare
        }
          ggplot(my_data, aes(x=Sample, y=.data[[taxa]], fill=Abundance))+
            geom_tile(color='gray')+
            theme(legend.justification='top', axis.text.x=element_text(angle=-90, hjust=0.5))+
            scale_x_discrete(position='top')
      })
      # obeserve the heatmap
      observe({output$heatmap <- renderPlot({heatmap()}, height=heatmap_height())})
    }
    # server side functions for alpha diversity visualization
    else if (input$tabs == 'tab8'){
      if(exists("Phylum")){
      # create a boxplot of alpha diversity depending on taxa, data type, and index
      which_alpha_plot <- reactive({
        taxa <- input$taxon_alpha_test
        my_list <- get(taxa)
        if(input$rawrare_alpha == 'Raw Data'){
          my_data <- my_list$diversityResults
        }
        else{
          my_data <- my_list$diversityResultsRarefy
        }
        
        which_div <- input$div_measure
        
        # coerce "my_data" into a dataframe so that the group_by function
        # is able to use it, group_by() only accepts tbl data type
        my_data = as.data.frame(my_data)
        
        
        data_median <- summarise(group_by(my_data, treatment), 
                                 MD = round(median(.data[[which_div]]),2))
        print(data_median)
        
        y_label <- names(diversity_choices)[grep(which_div, diversity_choices)]
        
        ggplot(my_data,aes(x=treatment,y=.data[[which_div]]))+
          geom_boxplot()+theme_bw()+labs(x="Treatment",y=y_label)+
          geom_text(data = data_median, aes(treatment, MD, label = MD), 
                    position = position_dodge(width = 0.8), size = 3, vjust = -0.5)
        })
        output$alpha_plots <- renderPlot({which_alpha_plot()})
        if(nrow(treatments)==2){
        #create caption from results of t-test
          alpha_caption <- reactive({
            taxa <- input$taxon_alpha_test
            my_list <- get(taxa)
            if(input$rawrare_alpha=="Raw Data"){
              my_data <- my_list$diversityResults
            }
            else{
              my_data <- my_list$diversityResultsRarefy
            }
          
            which_div <- input$div_measure
            
            result <- t.test(my_data[[which_div]]~my_data[,1])
            
            result_t <- round(result$statistic,digits=2)
            result_df <- round(result$parameter,digits=2)
            result_p <- round(result$p.value,digits=2)
            pvalue <- ""
            if(result_p<0.01){
              pvalue <- "p<0.01"
            }
            else{
              pvalue <- paste("p=",result_p)
            }
          
            #create caption from components of result from t-test
            caption <- HTML(paste("Welch's Two-sided T-test <br> t=",result_t,", df=",result_df,",",pvalue))
          })
        
        output$alpha_caption<-renderUI({alpha_caption()})
        }
      
        else{
          #create ANOVA table and post-hoc comparisons
          alpha_anova <- reactive({
            taxa <- input$taxon_alpha_test
            my_list <- get(taxa)
            
            if(input$rawrare_alpha=="Raw Data"){
              my_data <- my_list$diversityResults
            }
            else{
              my_data <- my_list$diversityResultsRarefy
            }
            
            anova_result <- 0
            posthoc_result <- 0
            anova_tables <- list()
            
            which_div <- input$div_measure
            
            anova_result <- anova(aov(my_data[[which_div]]~my_data[,1]))
            posthoc_result <- tidy(TukeyHSD(aov(my_data[[which_div]]~my_data[,1])))
            
            posthoc_result <- set(posthoc_result,,1,NULL)
            
            rownames(anova_result) <- c("Treatment","Residuals")
            anova_caption <-"ANOVA Table"
            posthoc_caption <-"Tukey's HSD"
            anova_tables <- list(anova_result,anova_caption,posthoc_result,posthoc_caption)
            
            #return list
            return(anova_tables)
          })
        
        #output ANOVA table and caption
        output$anova_caption <- renderText({alpha_anova()[[2]]})
        output$alpha_anova <- renderTable({alpha_anova()[[1]]},rownames=TRUE)
        
        #output post-hoc table and caption
        output$posthoc_caption <- renderText({alpha_anova()[[4]]})
        output$alpha_posthoc<-renderTable({alpha_anova()[[3]]})
      }
    }
      # close alpha diversity observe
    }
    # server side functions for beta diversity visualization
    else if (input$tabs == 'tab9'){
      # return the physeq object
      which_phy_seq <- reactive({
        taxa <- input$taxa_beta_test
        my_list <- get(taxa)
        if(input$rawrare_beta=='Raw Data'){
          my_physeq <- my_list$physeq
        }
        else{
          my_physeq <- my_list$physeq2
        }
        return(my_physeq)
      })
      
      # make ordination data 
      which_ordination_data <- reactive({
        my_ord_data <- ordinate(which_phy_seq(), method=input$ord_method, distance=input$dist_measure)
        return(my_ord_data)
      })
      
      # plot the ordination
      which_ordination_plot <- reactive ({
        if(input$samptreat=='Sample'){
          plot_ordination(which_phy_seq(), which_ordination_data, color = 'treatment') +
            stat_ellipse(type='t')+
            theme_bw()+
            coord_fixed()
        }
      })
      
      # make the ordination caption (NMDS stress)
      which_ordination_caption <- reactive({
        ord_data <- which_ordination_data()
        if(input$ord_method=='NMDS'){
          ord_cap <- paste("NMDS results: stress =", round(ord_data$stress, digits=4))
        }
        # may need to include an else statement here for other ordination methods
      })
      
      # plot ordination in UI
      output$ordination_plot <- renderPlot({which_ordination_plot()})
      output$ordination_caption <- renderText({which_ordination_caption})
      
      # return adonis results
      which_permanova <- reactive({
        taxa <- input$taxon_beta_test
        my_list <- get(taxa)
        if(input$rawrare_beta == 'Raw Data'){
          perm_method <- my_list$OTU_t
        }
        else{
          perm_method <- my_list$OTU_rarefy_t
        }
        
        perm_out <- adonis(perm_method~treatment, data=meta_global, method=input$dist_measure)
        return(as.data.frame(perm_out$aov.tab))
      })
      output$permanova_caption <- renderText({'PERMANOVA results'})
      output$beta_permanova <- renderTable({which_permanova()}, rownames=TRUE)
    }
    # close the observe event
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
