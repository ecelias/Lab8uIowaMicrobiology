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
library(shinyWidgets)

# define helper functions to be utilized in front and back end operations

# seperate the taxa by taxonomic level from the level 5 file
# parameters: dataframes from level 5 and metadata files
# returns: a list of all data separated by taxonomy
separateTaxa <- function(level5, meta){
  # create a global variable for the number of samples in the data
  numSamples <<- ncol(level5) - 1
  # sets the first column of the "level 5" taxonomy file to the index
  colnames(level5)[1] <- "index" 
  # separates the data within the level 5 file using ";" as the delimiter
  level5Sep <- separate(level5, "index", into=c("Kingdom", "Phylum", "Class", "Order", "Family"), sep=";")
  
  # dataset of samples with sample name, abundance and 1% threshold
  # uses "na.rm=TRUE to remove any null values"
  samples <- level5 %>% summarize_if(is.numeric, sum, na.rm = TRUE)
  samples <- pivot_longer(samples, 1:ncol(samples), names_to = "Sample", values_to = "Abundance")
  
  # set the threshold for the raw and rare data based on abundance
  samples$threshold <- samples$Abundance * 0.01
  samples$thresholdRare <- min(samples$Abundance) * 0.01
  
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
  phylum <- cbind(paste(level5Sep$Kingdom, level5Sep$Phylum), level5Sep[,6:ncol(level5Sep)])
  colnames(phylum)[1] <- "Phylum"
  
  class <- cbind(paste(level5Sep$Kingdom, level5Sep$Phylum, level5Sep$Class), level5Sep[,6:ncol(level5Sep)])
  colnames(class)[1] <- "Class"
  
  order <- cbind(paste(level5Sep$Kingdom, level5Sep$Phylum, level5Sep$Class, level5Sep$Order), level5Sep[,6:ncol(level5Sep)])
  colnames(order)[1] <- "Order"
  
  family <- cbind(paste(level5Sep$Kingdom, level5Sep$Phylum, level5Sep$Class, level5Sep$Order, level5Sep$Family), level5Sep[,6:ncol(level5Sep)])
  colnames(family)[1] <- "Family"
  
  # create list with separate taxa files, samples files, and treatment files
  # to access an item in this list: sep_taxa[[ITEM_NAME]] --> i.e sep_taxa[[phylum]]
  sepTaxa <- list(phylum=phylum, class=class, order=order, family=family, treatments=treatments, samples=samples)
  return(sepTaxa)
}

# function to create "other" category and reorder
# parameters: long data and rarified data
# returns: an updated dataframe sorted by abundance and merged with metadata
createOther <- function(longdata, rarified) {
  
  # debugging statement
  # if(rarified == TRUE){
  #  print(paste("Dimensions of longdata:", nrow(longdata), "rows,", ncol(longdata), "columns"))
  # }
  
  # long data with "other" category
  otherdata <- longdata
  otherdata <- otherdata %>% mutate_if(is.factor, as.character)
  
  otherrows = nrow(longdata)
  
  # create "Other" category
  for(i in 1:numSamples) {
    for(j in seq(i, otherrows, by=numSamples)) {
      # rarified data
      if(rarified == TRUE){
        if(otherdata$Abundance[j]<samples$thresholdRare[i]) {
          otherdata[j,1] <- "Other"
        }
      }
      # raw data
      else if (rarified == FALSE) {
        if(otherdata$Abundance[j]<samples$threshold[i]){
          otherdata[j,1] <- "Other"
        }
      }
    }
  }
  
  # convert taxa back to a factor
  otherdata <- otherdata %>% mutate_if(is.character, as.factor)
  otherdata <- as.data.frame(otherdata)
  
  # convert to wide format to sort taxa by overall abundance
  widedata <- otherdata
  widedata <- widedata %>% pivot_wider(names_from = "Sample", 
                                       values_from = "Abundance", 
                                       values_fill = list(Abundance = 0), 
                                       values_fn = list(Abundance = sum))
  
  # calculate overall abundance and reorder taxa
  widedata$abundance <- rowSums(widedata[,-1], na.rm = TRUE)
  widedata[[1]] <- reorder(widedata[[1]], -widedata$abundance)
  
  # delete overall abundance and convert back to long
  widedata <- widedata[, -ncol(widedata)]
  newdata <- widedata %>% pivot_longer(cols = -1, names_to = "Sample", values_to = "Abundance")
  
  # add treatment data from meta data
  newdata <- merge(newdata, metaGlobal, by.x = "Sample", by.y = 1)
  names(newdata)[ncol(newdata)] <- "treatment"
  
  return(newdata)
}

# creates all datasets needed in the server-side functions for bioinformatic analysis
# parameters: a list of all separated taxa, a taxonomic level
# returns: a list of all datasets needed for server-side functions
createDatasets <- function(sepTaxa, taxa) {
  
  # select taxonomic dataset
  x <- as.data.frame(sepTaxa[[taxa]])
  
  # count the number of samples
  numSamples = nrow(sepTaxa$samples)
  numeric_cols <- sapply(x, is.numeric)
  
  # count the number of samples in columns (this will give the abundance of each sample)
  # na.rm = TRUE will remove any null value
  columnData <- x %>%
    group_by_at(1) %>%
    summarize(across(where(is.numeric), sum, na.rm = TRUE))
  #reorder column data by overall abundance
  columnData$abundance<-rowSums(columnData[,-1], na.rm = TRUE)
  columnData[[1]] <- factor(columnData[[1]], levels = columnData[[1]][order(columnData$abundance, decreasing = TRUE)])
  # columnData[[1]] <- reorder(columnData[[1]], columnData$abundance)
  columnData <- columnData[,-ncol(columnData)]
  
  # long data
  longData <- columnData %>% pivot_longer(cols = -1, names_to = "Sample", values_to = "Abundance")
  
  # long data with other (raw) category
  longDataOther <- createOther(longData, FALSE)
  
  # create the matrices used to create the OTU and TAX tables with phyloseq
  # update the rownames with the column that contains species ("sp") identifiers
  otuMatrix <- data.matrix(subset(columnData, select=c(2:ncol(columnData))))
  rownames(otuMatrix) <- paste0("sp", 1:nrow(otuMatrix))
  taxaMatrix <- as.matrix(subset(columnData, select=c(1)))
  rownames(taxaMatrix) <- paste0("sp", 1:nrow(taxaMatrix))
  
  # create the tables for OTU and TAXA from the matrices
  OTU = otu_table(otuMatrix, taxa_are_rows=TRUE)
  TAX = tax_table(taxaMatrix)
  OTU_t <- t(OTU)
  
  # create a dataframe with the sample names as row names
  # then create a physeq object which will contain the OTU and TAX tables as 
  # well as a dataframe of all of the sample data
  physeqSamples <- metaGlobal
  rownames(physeqSamples) <- metaGlobal$sample
  physeq = phyloseq(OTU, TAX, sample_data(as.data.frame(physeqSamples)))
  
  # create rarefy datasets
  physeqRare <- rarefy_even_depth(physeq, rngseed = 10, sample.size = min(sample_sums(physeq)), replace = TRUE)
  otuRarefy <- otu_table(physeqRare)
  otuRarefy_t <- t(otuRarefy)
  taxaRarefy <- tax_table(physeqRare)
  
  # Align OTU identifiers between taxaRarefy and otuRarefy
  commonOTUs <- intersect(rownames(taxaRarefy), colnames(otuRarefy_t))
  
  # debugging statement to ensure that row names of TAXA and OTU tables are aligned
  if (length(commonOTUs) == 0) {
    stop("No common OTUs found between taxaRarefy and otuRarefy. Check your data.")
  }
  
  # Subset taxaRarefy and otuRarefy to include only common OTUs
  taxaRarefy <- taxaRarefy[commonOTUs, ]
  otuRarefy_t <- otuRarefy_t[, commonOTUs]
  
  # Convert taxaRarefy to a data frame
  taxaRarefy_df <- as.data.frame(taxaRarefy)
  rownames(taxaRarefy_df) <- commonOTUs
  
  # Convert otuRarefy_t to a data frame
  otuRarefy_t_df <- as.data.frame(otuRarefy_t)
  colnames(otuRarefy_t_df) <- commonOTUs
  
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
  longDataRare <- colDataRare %>% pivot_longer(cols = -1, 
                                               names_to = "Sample", 
                                               values_to = "Abundance")
  
  # long data (rarified) with other
  longDataRareOther <- createOther(longDataRare, TRUE)
  
  # create dataset with diversity indices
  richness <- specnumber(OTU_t)
  shannon <- diversity(OTU_t, index = "shannon")
  simpson <- diversity(OTU_t, index = "simpson")
  diversityResults <- cbind(metaGlobal[,2:ncol(metaGlobal)], 
                            richness, shannon, simpson)
  
  # update the column names of the diversity results 
  colnames(diversityResults) <- c("treatment", "richness", "shannon", "simpson")
  
  # create dataset with diversity indices based on rarefy data
  richnessRarefy <- specnumber(otuRarefy_t)
  shannonRarefy <- diversity(otuRarefy_t, index = "shannon")
  simpsonRarefy <- diversity(otuRarefy_t, index = "simpson")
  diversityResultsRarefy <- cbind(metaGlobal[,2:ncol(metaGlobal)], 
                                  richnessRarefy, shannonRarefy, simpsonRarefy)
  
  # rename columns of rarefied diversity results to be consisted with raw diversity results
  colnames(diversityResultsRarefy) <- c("treatment", "richness", "shannon", "simpson")
  
  # create a list with all of the datasets needed for each figure in the server
  datasetList <- list(columnData = columnData, longData = longData, 
                      longDataOther = longDataOther, OTU = OTU, OTU_t = OTU_t, 
                      physeq = physeq, otuRarefy = otuRarefy, taxa = TAX,
                      otuRarefy_t = otuRarefy_t, physeqRare = physeqRare, 
                      colDataRare = colDataRare, longDataRare = longDataRare, 
                      longDataRareOther = longDataRareOther, 
                      diversityResults = diversityResults,
                      diversityResultsRarefy = diversityResultsRarefy)
  
  return(datasetList)
  
}

# function to find core taxa found in all treatments and samples 
# parameters: columnData
# returns "core" variable which contains core taxa
coreTaxa <- function(colData) {
  core <- colData %>% filter_if(is.numeric, all_vars(.>0))
  core$abundance <- rowSums(core[,-1], na.rm = TRUE)
  print(core)
  #core <- core[order(-core$abundance)]
  return(core)
}

# create a title for Core Taxa table using metadata, core, and columnData
# parameters: dataframe (columnData) and core variable
# returns: a string with the defined caption
coreCaption <- function(df, core) {
  treatment <- paste(unique(metaGlobal$treatment), collapse = "/")
  numCore = nrow(core)
  total = nrow(df)
  caption = paste("<p class='text-light'>Core taxa in", treatment, "treatments: <b>", numCore, 
                  "of", total, "taxa are shared between treatments</b></p>")
  return(toString(caption))
}

# function to find taxa unique to a treatment
# parameters: columnData, treatments table
# returns: list of unique treatments
uniqueTaxa <- function(level5, treatments){
  uniqueList <- list()
  
  for(i in 1:nrow(treatments)) {
    uniqueSample <- level5 %>% select(!(treatments[i, "min_col"]:treatments[i, "max_col"])) %>% filter_if(is.numeric, all_vars(.==0))
    treatment <- as_tibble(treatments)
    
    if(nrow(uniqueSample) > 0){
      uniqueList[[i]] <- cbind(treatment[i,1], uniqueSample[,1])
      
      # add the names to the objects in the "uniqueList" to allow creation of other tables later
      name <- treatment[i, 1]
      names(uniqueList)[i] <- name
    }
  }
  return(uniqueList)
}

# function to graph rarefaction curve with ggplot2
graphRare <- function(x, ylab, bytype, graphTitle) {
  rareSample = list()
  
  # Write values for each sample from a rarecurve function into a list
  # Sample names obtained from "samples" dataframe
  for(i in 1:numSamples) {
    # rareSample[[i]] <- cbind(samples[i,1], as.data.frame(x[[i]]), attributes((x[[i]])$Subsample))
    rareSample[[i]] <- cbind(samples[i,1], as.data.frame(x[[i]]), attr(x[[i]], "Subsample"))
  }
  
  # Bind the data together from the "rareSample" list
  rareGraph <- bind_rows(rareSample)
  rareGraph <- merge(rareGraph, metaGlobal, by.x = "Sample", by.y = 1)
  colnames(rareGraph) <- c("Sample", "num_taxa", "num_samples", "Treatment")
  ggplot(rareGraph, 
         aes(x = num_samples, y = num_taxa, group = Sample, color = .data[[bytype]])) + 
    labs(x = "Sample Size", y = ylab) + theme_bw() + ggtitle(graphTitle) +
    geom_line(size = 1) + guides(fill = guide_legend(title = bytype)) +
    # adjust format of axis and legend text
    theme(
      axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
      axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
      axis.text.x = element_text(size = 14),   
      axis.text.y = element_text(size = 14),   
      legend.title = element_text(size = 16, , face="bold"),  
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 18, hjust = 0.5, face="bold")
    )
}

# create options for drop down menus used in UI
# drop down choices will be used to determine visualizations selected in certain tabs
taxachoices <- list("Phylum", "Class", "Order", "Family")
rawOrRare <- list("Raw"="Raw Data", "Rare"="Rarified Data")
absOrRel <- list("Absolute Abundance", "Relative Abundance")
sampleOrTreatment <- list("Sample", "Treatment")
diversityChoices <- list("Richness"="richness", "Shannon Diversity"="shannon", "Simpson Diversity"="simpson")
distanceChoices <- list("Jaccard" = "jaccard", "Bray-Curtis" = "bray", "Morisita-Horn" = "horn")
ordinationChoices <- list("NMDS", "PCoA")

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel(h1("Bean Beetle Microbiome Analysis", class="text-light")), 
  # create multiple pages in a tab bar
  tabsetPanel(id = "tabs", 
              # home page, only contains text
              tabPanel(value = "tab1", title = "Home", 
                       mainPanel(
                         h3("Welcome to the Bean Beetle Microbiome Analysis App", class="text-light"),
                         p(""),
                         h4("The Bean Beetle Microbiome Project is a research/teaching collaboration of institutions across the US that is studying the microbiome of", em("Callosobruchus maculatus"), "in research experiences (CUREs).", class="text-light"),
                         p("This app is designed to lead students through the community analysis of level-5 (family level) datasets produced in DNA subway.", class="text-light"),
                         p("Before proceeding, students should ensure the level 5 file is formatted correctly such that all unidentified taxa, chloroplasts, and mitochondria have been removed. 
                           The level 5 file should also be formatted with the first column as taxa and the subsequent columns as samples with unique sample identifiers.", class="text-light"),
                         p("Students should also prepare a metadata file in the first column as samples and the second column as treatments.", strong("Both files should be in .csv format."), class="text-light"), 
                         p("Additionally, this app should be capable of community analysis for any level 5 data.", class="text-light"),
                         p(""),
                         p("This app was reconfigured for the University of Iowa MICR 2158 course based on the source materials from Huang et al., 2022 by Elizabeth Elias, an undergraduate student at the University of Iowa under the guidance of Dr. Regina McGrane, Department of Microbiology and Immunology, University of Iowa.", class="text-light"),
                         p(tags$a("The original app can be found by clicking here.", href = "https://beanbeetles.shinyapps.io/BeanBeetleMicrobiome/")),
                         p(tags$a("For more information on this CURE project, please click here", href = "https://www.beanbeetles.org/microbiome/the-bean-beetle-microbiome-project/")),
                         p(tags$a("Click here to find the GitHub repo for this project", href = "https://github.com/ecelias/Lab8uIowaMicrobiology"))
                       )
              ),
              # page for users to upload data
              # displays a table with their level5 and metadata as well as a rankabundance curve
              tabPanel(value = "tab2", title = "Data Upload", 
                       sidebarLayout(
                         sidebarPanel(
                           fileInput("file1", p("Choose level 5 CSV file",class="text-light"), 
                                     accept = c( 
                                       "text/csv",
                                       "text/comma-separated-values, text/plain",
                                       ".csv")), # accept a file input
                           tags$hr(),
                           fileInput("file2", p("Choose metadata CSV file",class="text-light"), 
                                     accept = c(
                                       "text/csv",
                                       "text/comma-separated-values, text/plain",
                                       ".csv")),
                           tags$hr(),
                           card(
                             p("Click run app after selecting files",class="text-light"),
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
                                 card(
                                   height = 350,
                                   card_header(
                                     class = "bg-primary mb-3",
                                     "Level 5 File"
                                   ),
                                   # card to hold the level 5 data
                                   card_body(
                                     tags$i(p("Ensure the first column is the combined taxa (including kingdom) separated by semi-colons. 
                                              The remaining columns should contain the abundance data for each sample. 
                                              The 6 most abundant taxa will be displayed below.",class="text-light")),
                                     tags$hr(),
                                     tableOutput("level5Contents"),
                                     tags$hr()
                                   )
                                 ),
                                 # card to hold the metadata
                                 card(
                                   height = 350,
                                   card_header(
                                     class = "bg-primary mb-3",
                                     "Metadata File"
                                   ),
                                   card_body(
                                     tags$i(p("Ensure the first column lists the samples and the second column lists the treatments for each sample. 
                                              Your metadata is displayed below.",class="text-light")),
                                     tags$hr(),
                                     tableOutput("metadataContents"),
                                     tags$hr()
                                   )
                                 )
                               
                             )
                             # card for the rank abundance curve
                             # card(
                               # height = 650,
                               # card header defines the format of the card, additional
                               # card formats can be found on bootswatch
                               # card_header(class = "bg-secondary mb-3","Rank Abundance Curve"),
                               # card_body(tags$i(p("A rank abundance curve or Whittaker plot represents the relative abundance of species within a community.
                                          # by plotting species abundance against its rank order, providing a visualization of both species
                                          # richness and species evenness. A steeper curve indicates that a community is comprised of a few dominant species
                                          # while the others are relatively rare whereas a flatter curve suggests more more even distrubution of species",
                                          # class="text-light")), tags$hr(),plotOutput("rankAbundanceCurve")))
                           )
                         )
                       )
              ),
              
              # UI to display Core Taxa present in samples
              tabPanel(value = "tab3", title = "Core Taxa",
                       sidebarLayout(
                         fluid = TRUE,
                         sidebarPanel(
                           card(
                             # uses shinyWidgets package to create a vertical list of buttons
                             # to select which taxon they want to view data for
                             radioGroupButtons(
                               inputId = "taxonCore",
                               label = p("Select a taxonomic level:",class="text-light"),
                               choices = taxachoices,
                               direction = "vertical"
                             ),
                             ## Use this to create a dropdown menu
                             #selectInput("taxonCore", "Select a taxon", taxachoices), width = 2
                           ), 
                           # manually set the width to 3
                           width=3
                         ),
                         mainPanel(
                           # header text for the main panel
                           h2("Core Taxa",class="text-light"),
                             p("Core taxa are those taxa found in all samples.",class="text-light"),
                             # ensure the text is above the table
                             verticalLayout(htmlOutput("coreCaption"), tableOutput("coreTaxa")), 
                         )
                       )
                       
              ), 
              # UI to display Unique Taxa present in samples
              tabPanel(value = "tab4", title = "Unique Taxa", 
                       sidebarLayout(
                         sidebarPanel(
                           card(
                             radioGroupButtons(
                               inputId = "taxonUnique",
                               label = p("Select a taxonomic level:",class="text-light"),
                               choices = taxachoices,
                               direction = "vertical"
                             )
                           ), 
                           width=3
                         ),
                         mainPanel(
                           h2("Unique Taxa",class="text-light"),
                           card(
                             p("Unique taxa are those taxa found in a single treatment.",class="text-light"),
                             textOutput('noUniqueCaptions'),
                             uiOutput('uniqueTables')
                           )
                         )
                       )
                       
              ),
              # UI to display Rarefaction plots 
              tabPanel(value = "tab5", title = "Rarefaction", 
                       sidebarLayout(
                         sidebarPanel(
                           # provides option for user to select taxonomic levels
                           radioGroupButtons(
                             inputId = "taxonRarefy",
                             label = p("Select a taxonomic level:",class="text-light"),
                             choices = taxachoices,
                             direction = "vertical"
                           ),
                           # provides option for user to select if they want to graph
                           # by sample or treatment
                           radioGroupButtons(
                             inputId = "byType",
                             label = p("Graph By:", class="text-light"),
                             choices = sampleOrTreatment,
                             justified = TRUE, 
                             direction = "vertical"
                           ), 
                           width=3
                         ),
                         mainPanel(
                           h2("Sample Rarefaction Curves",class="text-light"), 
                             card(
                               card_header(
                                 class = "bg-primary mb-3",
                                 "Raw Data"
                               ),
                               plotOutput('initialRarefaction', width='100%', height='400px') %>%
                                 withSpinner(color='mediumaquamarine'), 
                               # provides option to expand plot into full screen mode
                               full_screen = TRUE, 
                             ),
                           # provides a button to download a png of the plot
                           downloadButton("downloadInitialRarefaction", "Download Plot", class="btn-sm"),
                           tags$hr(),
                             card(
                               card_header(
                                 class = "bg-primary mb-3",
                                 "Even Rarefaction to Minimum Number of Sequences"
                               ),
                               plotOutput('evenRarefaction', width='100%', height='400px') %>%
                                 withSpinner(color='darkseagreen'), 
                               full_screen = TRUE,
                             ), 
                           downloadButton("downloadEvenRarefaction", "Download Plot", class="btn-sm")
                         )
                       )
              ),
              tabPanel(value = "tab6", title = "Taxonomy Bar Graphs", fluid = TRUE,
                       sidebarLayout(
                         sidebarPanel(
                           radioGroupButtons(
                             inputId = "taxonBar",
                             label = p("Select a taxonomic level:",class="text-light"),
                             choices = taxachoices,
                             direction = "vertical"
                           ),
                           radioGroupButtons(
                             inputId = "rawrareBar",
                             label = p("Select which data to use:", class="text-light"),
                             choices = rawOrRare
                           ),
                           radioGroupButtons(
                             inputId = "absRel",
                             label = p("Graph By:", class="text-light"),
                             choices = absOrRel,
                             direction = "vertical"
                           ), 
                           # provides an option for the user to toggle the legend
                           # on and off, mainly available for saving the figure in a preferred format
                           p("Display Legend?", class="text-light"),
                           switchInput(
                             inputId = "hideLegend",
                             onLabel = "Show",
                             offLabel = "Hide"
                           ), 
                           width=3
                         ),
                         mainPanel(
                           h2("Taxonomy Bar Graph",class="text-light"),
                            plotOutput('bargraph',width='100%', height='auto') %>%
                              withSpinner(color='thistle'),
                          downloadButton("downloadBar", "Download Plot", class="btn-sm")
                         )
                       )
              ),
              tabPanel(value = "tab7", title = "Taxonomy Heat Map", 
                       sidebarLayout(
                         sidebarPanel(
                           radioGroupButtons(
                             inputId = "taxonHeatmap",
                             label = p("Select a taxonomic level:",class="text-light"),
                             choices = taxachoices,
                             direction = "vertical"
                           ),
                           radioGroupButtons(
                             inputId = "rawrareHeatmap",
                             label = p("Select which data to use:", class="text-light"),
                             choices = rawOrRare
                           ), 
                           width=3
                         ),
                         mainPanel(
                           h2("Taxonomy Heatmap",class="text-light"),
                           card(
                             plotOutput('heatmap', height="auto") %>%
                               withSpinner(color='lemonchiffon'), 
                             full_screen = TRUE),
                           downloadButton("downloadHeatmap", "Download Plot", class="btn-sm")
                         )
                       )
              ),
              tabPanel(value = "tab8", title = "Alpha Diversity", 
                       sidebarLayout(
                         sidebarPanel(
                           radioGroupButtons(
                             inputId = "taxonAlphaTest",
                             label = p("Select a taxonomic level:",class="text-light"),
                             choices = taxachoices,
                             direction = "vertical"
                           ),
                           radioGroupButtons(
                             inputId = "rawrareAlpha",
                             label = p("Select which data to use:", class="text-light"),
                             choices = rawOrRare,
                           ),
                           radioGroupButtons(
                             inputId = "divMeasure",
                             label = p("Select a diversity measure:", class="text-light"),
                             choices = diversityChoices,
                             direction="vertical"
                           ), 
                           width=3
                        ),
                         mainPanel(
                           h3('Alpha Diversity', class="text-light"),
                           # creates a panel so user's can easily switch between viewing
                           # a single taxonomic levele and viewing all taxonomic levels
                           # side by side for easy comparison
                           navset_card_underline(
                             title = h5("Visualizations", class="text-light"),
                             # panel with single plot
                             nav_panel("Scaled", 
                                       plotOutput("alphaPlots")  %>%
                                         withSpinner(color='lavender')
                             ), 
                             nav_panel("Side-by-Side", 
                                       card(
                                         layout_columns(
                                           card(
                                             h6("Phylum", class="text-light"), 
                                             plotOutput("phylumAlpha")  %>%
                                               withSpinner(color='plum')
                                           ), 
                                           card(
                                             h6("Class", class="text-light"), 
                                             plotOutput("classAlpha")  %>%
                                               withSpinner(color='hotpink')
                                           )
                                         ),
                                         layout_columns(
                                           card(
                                             h6("Order", class="text-light"), 
                                             plotOutput("orderAlpha")  %>%
                                               withSpinner(color='violet'), 
                                           ), 
                                           card(
                                             h6("Family", class="text-light"), 
                                             plotOutput("familyAlpha")  %>%
                                               withSpinner(color='palevioletred')
                                           )
                                         )
                                       )
                             )
                           ), 
                           downloadButton("downloadAlpha", "Download Plot", class="btn-sm"),
                           tags$hr(), 
                           # display statistics for alpha diversity as a table
                           h5("Welch's Two-Sided T-test:", class="text-light"),
                           tableOutput('alphaStats'),
                           # displays ANOVA and posthoc results 
                           # likely not applicable for UI students
                           textOutput('anovaCaption'),
                           tableOutput('alphaAnova'), 
                           textOutput('posthocCaption'),
                           tableOutput('alphaPosthoc')
                         )
                       )
              ),
              tabPanel(value = "tab9", title = "Beta Diversity", 
                       sidebarLayout(
                         sidebarPanel(
                           radioGroupButtons(
                             inputId = "taxonBetaTest",
                             label = p("Select a taxonomic level:",class="text-light"),
                             choices = taxachoices,
                             direction = "vertical"
                           ),
                           radioGroupButtons(
                             inputId = "rawrareBeta",
                             label = p("Select which data to use:", class="text-light"),
                             choices = rawOrRare,
                           ),
                           radioGroupButtons(
                             inputId = "distMeasure",
                             label = p("Select distance measure:", class="text-light"),
                             choices = distanceChoices,
                             direction="vertical"
                           ),
                           radioGroupButtons(
                             inputId = "ordMethod",
                             label = p("Select ordination method:", class="text-light"),
                             choices = ordinationChoices,
                           ),
                           radioGroupButtons(
                             inputId = "samptreat",
                             label = p("Graph by:", class="text-light"),
                             choices = sampleOrTreatment,
                           ), 
                           width=3
                         ),
                         mainPanel(
                           h3('Beta Diversity', class='text-light'),
                           # creates a panel so user's can easily switch between scaled
                           # and unscaled data in addition to viewing all taxonomic levels
                           # side by side for easy comparison
                           navset_card_underline(
                             title = h5("Visualizations", class="text-light"),
                             # panel with unscaled plots
                             nav_panel("Unscaled", 
                                       plotOutput("unscaledOrdPlot") %>%
                                         withSpinner(color='darkslateblue')
                             ),
                             # panel with scaled plots 
                             nav_panel("Scaled", 
                                       plotOutput("scaledOrdPlot")  %>%
                                         withSpinner(color='darkcyan')
                             ), 
                             nav_panel("Side-by-Side", 
                                       card(
                                         layout_columns(
                                           card(
                                             h6("Phylum", class="text-light"), 
                                             plotOutput("phylumBeta")  %>%
                                               withSpinner(color='cornflowerblue')
                                           ), 
                                           card(
                                             h6("Class", class="text-light"), 
                                             plotOutput("classBeta")  %>%
                                               withSpinner(color='cadetblue')
                                           )
                                         ),
                                         layout_columns(
                                           card(
                                             h6("Order", class="text-light"), 
                                             plotOutput("orderBeta")  %>%
                                               withSpinner(color='paleturquoise'), 
                                           ), 
                                           card(
                                             h6("Family", class="text-light"), 
                                             plotOutput("familyBeta")  %>%
                                               withSpinner(color='royalblue')
                                           )
                                         )
                                       )
                             )
                           ), 
                           tags$hr(),
                           htmlOutput('ordinationCaption'),
                           downloadButton("downloadBetaUnscaled", "Download Unscaled Plot", class="btn-sm"),
                           downloadButton("downloadBetaScaled", "Download Scaled Plot", class="btn-sm"),
                           tags$hr(), 
                           h5("PERMANOVA Result:", class="text-light"),
                           tableOutput('betaPermanova')
                         )
                       )
                       
              )
  ),
  theme = bs_theme(preset = "slate")
  # close UI
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # variables for level_5 and metadata
  level5 <- 0
  metadata <- 0
  
  observeEvent(input$tabs,{
    # if the user goes to tab 2
    # functions to accept the csv files and return a data frame for level 5 and metadata
    # displays the tables the user inputs
    # will also generate a rank abundance curve 
    if(input$tabs == "tab2"){
      l5Upload <- reactive({
        inFile1 <- input$file1
        if (is.null(inFile1))
          return(NULL)
        read_csv(inFile1$datapath)
      })
      metaUpload <- reactive({
        inFile2 <- input$file2
        if (is.null(inFile2))
          return(NULL)
        read_csv(inFile2$datapath)
      })
      
      # upload level 5 and metadata files
      # convert to table format for cards
      observeEvent(input$run, {
        level5 <- l5Upload()
        level5 %>% filter_all(any_vars(. != 0))
        metadata <- metaUpload()
        output$level5Contents <- renderTable({head(level5)})
        output$metadataContents <- renderTable({metadata})
        
        
        # set column names in metadata to sample and treatment
        colnames(metadata) <- c("sample", "treatment")
        
        sepTaxa <- separateTaxa(level5, metadata)
        
        # create a global dataset for treatments and samples, indicated by double arrows
        samples <<- sepTaxa$samples
        treatments <<- sepTaxa$treatments
        metaGlobal <<- as.data.frame(metadata)
        
        # Prepare data for rank abundance calculation
        #raData <- level5 %>%
          #gather(Sample, Abundance, -Sample) %>%  # Reshape data
          #group_by(Sample) %>%
          #arrange(desc(Abundance), .by_group = TRUE) %>%  # Ensure correct ranking order
          #mutate(Rank = row_number()) %>%  # Assign ranks correctly
          #ungroup()  # Remove grouping to avoid errors in conversion
        
        #raData$Abundance <- as.numeric(raData$Abundance)
        #raData$Abundance <- as.numeric(raData$Abundance)
        #raVector <- sort(raData$Abundance, decreasing = TRUE)
        
        #raPlot <- radplot(raVector, Whittaker=FALSE, size=1) +
          #xlab("Rank") + ylab("Abundance") +
          #ggtitle("") + theme_bw() +
          #theme(axis.text.x = element_blank())
        
        # output the rank abundance plot
        #output$rankAbundanceCurve <- renderPlot({raPlot})
        
        # create global datasets for each taxonomic level
        Phylum <<- createDatasets(sepTaxa, "phylum")
        Class <<- createDatasets(sepTaxa, "class")
        Order <<- createDatasets(sepTaxa, "order")
        Family <<- createDatasets(sepTaxa, "family")
      })
    }
    # server side functions for core taxa visualization
    else if(input$tabs == 'tab3'){
      # all tabs after this will check if data has been uploaded before
      # attempting to display any data to prevent errors from popping up
      if(exists('Phylum')){
        
        # creates core dataset based on taxonomic level selected for tab4
        coreData <- reactive ({
          # chooses which taxonomic level to obtain data for based on user input
          taxa <- input$taxonCore
          # gets the correct taxonomic level data from the global datasets created after data upload
          myList <- get(taxa)
          # selects the column data from the list of all taxa-specific datasets
          core <- coreTaxa(myList$columnData)
          print(core)
          return(core)
        })
        
        # creates the caption for the core dataset based on the taxonomic level selected
        coreDataCaption <- reactive({
          taxa <- input$taxonCore
          myList <- get(taxa)
          core <- coreTaxa(myList$columnData)
          caption <- coreCaption(myList$columnData, core)
          return(caption)
        })
        # returns caption and core taxa table to output
        output$coreCaption <- renderUI({HTML(coreDataCaption())})
        output$coreTaxa <- renderTable({coreData()},
                                       digits=0, striped=TRUE, bordered=TRUE)
      }
    }
    # server side functions for unique taxa visualization
    else if (input$tabs == 'tab4'){
      if(exists('Phylum')){
        
        # create unique taxa output for tab 4 
        # creates unique taxa datset based on taxonomic level selected
        # general function: create a list with different objects for each treatment
        uniqueData <- reactive({
          taxa <- input$taxonUnique
          myList <- get(taxa)
          uniqueTaxaList <- list()
          uniqueTaxaList <- uniqueTaxa(myList$columnData, treatments)
          
          return(uniqueTaxaList)
        })
        
        # creates a list of the number of output tables based on the number of treatments 
        # tables are labeled with treatment names
        output$uniqueTables <-
          renderUI({
            tableOutputList <- lapply(treatments$treatment, function(i){
              tablename <- paste0("table", i)
              htmlOutput(tablename)
            })
            # generates a list of HTML tags Shiny uses to categorize each table output
            # allows Shiny to properly render the unique tables
            tagList(tableOutputList)
          })
        # render the table
        # use the "observe" command to allow reactive functions to run
        observe({
          uniqueDataTables <- uniqueData()
          
          if (!length(uniqueDataTables))
          {
            output$noUniqueCaption <- renderText({"<p class='text-light'>There are no unique taxa.</p>"})
          }
          
          for (i in 1:nrow(treatments)) {
            local({
              myI <- treatments[i,1]
              # create the caption
              numUnique <- length(uniqueDataTables[[myI]]$treatment)
              taxa <- input$taxonUnique
              myList <- get(taxa)
              total <- nrow(myList$columnData)
              # change the text color by changing the class. 
              # text color options can be viewed on bootswatch.com
              captionTitle = paste("<p class='text-light'>Unique taxa in", myI, "treatment: <b>", numUnique, 
                              "of", total, "taxa are unique to this treatment</b></p>")
              # provide functionality to render the table and correctly place the caption
              tablename <- paste0("table", myI)
              output[[tablename]] <- renderTable(
                {
                  uniqueDataTables[[myI]]
                },
                caption = captionTitle, caption.placement = getOption(
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
        whichInitialRarefaction <- reactive({
          taxa <- input$taxonRarefy
          byType <- input$byType
          myList <- get(taxa)
          myData <- myList$OTU_t
          
          # Coerce the transposed OTU table into a matrix. rarecurve() only 
          # accepts matrix-like objects and will not accept an OTU table
          myData <- otu.matrix(myData)
          rarecurveData <- rarecurve(myData, step=20, sample=20, xlab="Sample Size", ylab="Species", label=FALSE, tidy=FALSE)
          graphRare(rarecurveData, taxa, byType, "Initial Rarefaction")
        })
        
        output$initialRarefaction <- renderPlot({whichInitialRarefaction()})
        
        # functionality to download the initial rarefaction plot
        output$downloadInitialRarefaction <- downloadHandler(
          filename = function() {
            paste("initial_rarefaction_plot.png", sep="")
          }, 
          content = function(file) {
            width <- 2000
            height <- 2000
            png(file=file, width=width, height=height, res=150)
            plot(whichInitialRarefaction())
            dev.off()
          }
        )
        
        # Select data and create even rarefaction graph
        # Uses the rarified OTU table which must be coerced into a matrix
        # prior to visualization using otu.matrix
        whichEvenRarefaction <- reactive ({
          taxa <- input$taxonRarefy
          byType <- input$byType
          myList <- get(taxa)
          myData <- myList$otuRarefy_t
          myData <- otu.matrix(myData)
          rarecurveData <- rarecurve(myData, step=20, sample=20, xlab="Sample Size", ylab="Species", label=FALSE, tidy=FALSE)
          graphRare(rarecurveData, taxa, byType, "Even Rarefaction")
        })
        
        output$evenRarefaction <- renderPlot({whichEvenRarefaction()})
        
        # Functionality for downloading the the even rarefaction plot
        output$downloadEvenRarefaction <- downloadHandler(
          filename = function() {
            paste("even_rarefaction_plot.png", sep="")
            }, 
          content = function(file) {
            width <- 2000
            height <- 2000 
            png(file=file, width=width, height=height, res=150)
            plot(whichEvenRarefaction())
            dev.off()
          }
        )
        
      }
    }
    
    # server side functions for bar graphs
    else if (input$tabs == 'tab6'){
      if(exists('Phylum')) {
        whichBarGraph <- reactive({
          taxa <- input$taxonBar
          myList <- get(taxa)
          if(input$rawrareBar=='Raw Data'){
            myData <- myList$longDataOther
          }
          else {
            myData <- myList$longDataRareOther
          }
          # allows user to hide legend if they prefer to download the plot that way
          if(!input$hideLegend){
            # if user selects absolute abundance to graph by
            if(input$absRel == 'Absolute Abundance'){ 
              ggplot(data=myData, aes(x=Sample, y=Abundance)) +
                geom_bar(aes(fill=myData[,2]), position='stack', stat='identity')+
                labs(fill=taxa, y='Absolute Abundance')+
                facet_grid(.~treatment, space='free_x', scales='free_x')+
                theme(legend.position='bottom')+
                guides(fill=guide_legend(ncol=2)) +
                ggtitle("")+
                # customizes text elements 
                theme(
                  axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                  axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                  axis.text.x = element_text(size = 14, angle=45, vjust=0.5),   
                  axis.text.y = element_text(size = 14),   
                  legend.title = element_text(size = 14, , face="bold"),  
                  legend.text = element_text(size = 10),
                  legend.title.position = "top",
                  strip.text = element_text(size = 16), 
                  plot.title = element_text(size = 18, hjust = 0.5, face="bold")
                )
            }
            else{
              ggplot(data=myData, aes(x=Sample, y=Abundance))+
                geom_bar(aes(fill=myData[,2]), position='fill', stat='identity')+
                labs(fill=taxa, y='Relative Abundance')+
                facet_grid(.~treatment, space='free_x', scales='free_x')+
                theme(legend.position='bottom')+
                guides(fill=guide_legend(ncol=2))+
                theme(
                  axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                  axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                  axis.text.x = element_text(size = 14, angle=45, vjust=0.5),   
                  axis.text.y = element_text(size = 14),   
                  legend.title = element_text(size = 14, , face="bold"),  
                  legend.text = element_text(size = 10),
                  strip.text = element_text(size = 16),
                  plot.title = element_text(size = 18, hjust = 0.5, face="bold")
                )
            }
          } else{
            if(input$absRel == 'Absolute Abundance'){ # absolute abundance
              ggplot(data=myData, aes(x=Sample, y=Abundance)) +
                geom_bar(aes(fill=myData[,2]), position='stack', stat='identity')+
                labs(fill=taxa, y='Absolute Abundance')+
                facet_grid(.~treatment, space='free_x', scales='free_x')+
                theme(legend.position='none')+
                guides(fill=guide_legend(ncol=2)) +
                ggtitle("")+
                theme(
                  axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                  axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                  axis.text.x = element_text(size = 14, angle=45, vjust=0.5),   
                  axis.text.y = element_text(size = 14),   
                  strip.text = element_text(size = 16), 
                  plot.title = element_text(size = 18, hjust = 0.5, face="bold")
                )
            }
            else{
              ggplot(data=myData, aes(x=Sample, y=Abundance))+
                geom_bar(aes(fill=myData[,2]), position='fill', stat='identity')+
                labs(fill=taxa, y='Relative Abundance')+
                facet_grid(.~treatment, space='free_x', scales='free_x')+
                theme(legend.position='none')+
                guides(fill=guide_legend(ncol=2))+
                theme(
                  axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                  axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                  axis.text.x = element_text(size = 14, angle=45, vjust=0.5),   
                  axis.text.y = element_text(size = 14),   
                  strip.text = element_text(size = 16),
                  plot.title = element_text(size = 18, hjust = 0.5, face="bold")
                )
            }            
          }
        })
        
        # determine height of the taxonomy bar graphs in tab3
        barGraphHeight <- reactive({
          taxa <- input$taxonBar
          myList <- get(taxa)
          myData <- myList$longDataOther
          numTaxa <- length(unique(myData[,2]))
          height = 500 + numTaxa*10
          return(height)
        })
        
        # output the taxonomy bar graph with the correct height
        observe({output$bargraph <- renderPlot({whichBarGraph()}, height = barGraphHeight())})
        
        # Functionality for downloading the the Bar Graph plot
        output$downloadBar <- downloadHandler(
          filename = function() {
            paste("taxa_bargraph_plot.png", sep="")
          }, 
          content = function(file) {
            width <- 5000
            height <- 5000
            png(file=file, width=width, height=height, res=300)
            plot(whichBarGraph(), height=barGraphHeight())
            dev.off()
          }
        )
      }
    }
    
    # server side functions for taxa heatmaps
    else if (input$tabs == 'tab7'){
      if(exists("Phylum")){
        heatmapHeight <- reactive({
          taxa <- input$taxonHeatmap
          myList <- get(taxa)
          myData <- myList$longDataOther
          numTaxa <- nrow(unique(myData[,1]))
          
          # adjusts the height of the plot based on which taxonomic level is selected
          if(taxa == "Phylum") {
            height <- max(c(400, numTaxa*15))
          } else if(taxa == "Class") {
            height <- max(c(750, numTaxa*15))
          } else if(taxa == "Order") {
            height <- max(c(1000, numTaxa*15))
          } else if(taxa == "Family") {
            height <- max(c(1500, numTaxa*15))
          }
          return(height)
        })
        
        heatmap <- reactive({
          taxa <- input$taxonHeatmap
          myList <- get(taxa)
          if(input$rawrareHeatmap=='Raw Data'){
            myData <- myList$longData
          }
          else{ 
            myData <- myList$longDataRare
          }
          
          # ensures that y-axis labels only consist of that samples taxonomic level
          # this shortens the y-axis labels significantly allowing more readbility of figures 
          if(taxa == "Phylum") {
            # Extract the second part of the species name (p__Phylum)
            yLabel <- sapply(as.character(myData[[taxa]]), 
                             function(x) strsplit(x, " ")[[1]][2])
            labelTitle <- "Phylum"
          } else if(taxa == "Class") {
            # Extract the third part of the species name (c__Class)
            yLabel <- sapply(as.character(myData[[taxa]]), 
                             function(x) strsplit(x, " ")[[1]][3])  
            labelTitle <- "Class"
          } else if(taxa == "Order") {
            # Extract the fourth part of the species name (o__Order)
            yLabel <- sapply(as.character(myData[[taxa]]), 
                             function(x) strsplit(x, " ")[[1]][4])  
            labelTitle <- "Order"
          } else if(taxa == "Family") {
            # Extract the fifth part of the species name (f__Family)
            yLabel <- sapply(as.character(myData[[taxa]]), 
                             function(x) strsplit(x, " ")[[1]][5])  
            labelTitle <- "Family"
          }
          
          # generate the heatmap
          ggplot(myData, aes(x=Sample, y=yLabel, fill=Abundance))+
              geom_tile(color='gray')+ labs(y=labelTitle)+
              theme(
                legend.justification='top',
                axis.text.x=element_text(angle=-90, size =10),
                axis.title.x = element_text(size = 14, margin = margin(t = 30)),
                axis.title.y = element_text(size = 14, margin = margin(r = 10)),
                axis.text.y=element_text(size=8),
                )+
              scale_x_discrete(position='top')  # places sample IDs on the top of the graph
        })
        
        # obeserve the heatmap
        observe({output$heatmap <- renderPlot({heatmap()}, height=heatmapHeight())})
        
        # Functionality for downloading the the Bar Graph plot
        output$downloadHeatmap <- downloadHandler(
          filename = function() {
            paste("taxa_heatmap_plot.png", sep="")
          }, 
          content = function(file) {
            width <- 2000
            height <- 2000 
            png(file=file, width=width, height=height, res=150)
            plot(heatmap(), height=heatmapHeight())
            dev.off()
          }
      )
      }
    }
    # server side functions for alpha diversity visualization
    else if (input$tabs == 'tab8'){
      if(exists("Phylum")){
        # create a boxplot of alpha diversity depending on taxa, data type, and index
        whichAlphaPlot <- reactive({
          taxa <- input$taxonAlphaTest
          myList <- get(taxa)
          if(input$rawrareAlpha == 'Raw Data'){
            myData <- myList$diversityResults
          }
          else{
            myData <- myList$diversityResultsRarefy
          }
  
          whichDiv <- input$divMeasure
          
          # coerce "myData" into a dataframe so that the group_by function
          # is able to use it, group_by() only accepts tbl data type
          # Additionally, ensure data is numeric 
          myData = as.data.frame(myData)
          myData[[whichDiv]] <- as.numeric(myData[[whichDiv]])
          
          # for each treatment, calculate the median of whatever diversity measure the user selects
          # then create a tbl with the treatments and the median of the div. measure
          dataMedian <- summarise(group_by(myData, treatment), 
                                  MD = round(median(as.numeric(.data[[whichDiv]])), 2))
          
          # selects the value for the y-axis based on the list of possible choices 
          # to correspond with user selection of the diversity measure they want to use
          yLabel <- names(diversityChoices)[grep(whichDiv, diversityChoices)]
          
          ggplot(myData,aes(x=treatment,y=.data[[whichDiv]], fill=treatment))+
            # "alpha 0.3" makes the fill color of the boxes transluscent
            geom_boxplot(alpha=0.3)+theme_bw()+labs(x="Treatment",y=yLabel)+
            geom_text(data = dataMedian, aes(treatment, MD, label = MD), 
                      position = position_dodge(width=0.8), # displays the median value of each boxplot inside the plot
                      size = 5, vjust = -0.5, hjust = 0.5)+
            theme(
              axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
              axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
              axis.text.x = element_text(size = 14),   
              axis.text.y = element_text(size = 14),   
              strip.text = element_text(size = 16), 
              legend.position="none"
              #plot.title = element_text(size = 18, hjust = 0.5, face="bold")
            ) +
            scale_fill_brewer(palette="Accent") # brewer color palette used to fill the box plots
          })
        
        # create a boxplot of alpha diversity from phylum data only
        phylumAlphaPlot <- reactive({
          myList <- get("Phylum")
          if(input$rawrareAlpha == 'Raw Data'){
            myData <- myList$diversityResults
          }
          else{
            myData <- myList$diversityResultsRarefy
          }
          
          whichDiv <- input$divMeasure
          myData = as.data.frame(myData)
          myData[[whichDiv]] <- as.numeric(myData[[whichDiv]])
          dataMedian <- summarise(group_by(myData, treatment), 
                                  MD = round(median(as.numeric(.data[[whichDiv]])), 2))
          yLabel <- names(diversityChoices)[grep(whichDiv, diversityChoices)]
          
          ggplot(myData,aes(x=treatment,y=.data[[whichDiv]], fill=treatment))+
            geom_boxplot(alpha=0.3)+theme_bw()+labs(x="Treatment",y=yLabel)+
            geom_text(data = dataMedian, aes(treatment, MD, label = MD), 
                      position = position_dodge(width=0.8),
                      size = 5, vjust = -0.5, hjust = 0.5)+
            theme(
              axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
              axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
              axis.text.x = element_text(size = 14),   
              axis.text.y = element_text(size = 14),   
              strip.text = element_text(size = 16), 
              legend.position="none"
            ) +
            scale_fill_brewer(palette="Accent") 
        })
        
        # create a boxplot of alpha diversity from class data only
        classAlphaPlot <- reactive({
          myList <- get("Class")
          if(input$rawrareAlpha == 'Raw Data'){
            myData <- myList$diversityResults
          }
          else{
            myData <- myList$diversityResultsRarefy
          }
          
          whichDiv <- input$divMeasure
          myData = as.data.frame(myData)
          myData[[whichDiv]] <- as.numeric(myData[[whichDiv]])
          dataMedian <- summarise(group_by(myData, treatment), 
                                  MD = round(median(as.numeric(.data[[whichDiv]])), 2))
          yLabel <- names(diversityChoices)[grep(whichDiv, diversityChoices)]
          
          ggplot(myData,aes(x=treatment,y=.data[[whichDiv]], fill=treatment))+
            geom_boxplot(alpha=0.3)+theme_bw()+labs(x="Treatment",y=yLabel)+
            geom_text(data = dataMedian, aes(treatment, MD, label = MD), 
                      position = position_dodge(width=0.8),
                      size = 5, vjust = -0.5, hjust = 0.5)+
            theme(
              axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
              axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
              axis.text.x = element_text(size = 14),   
              axis.text.y = element_text(size = 14),   
              strip.text = element_text(size = 16), 
              legend.position="none"
            ) +
            scale_fill_brewer(palette="Accent") 
        })
        
        # create a boxplot of alpha diversity from order data only
        orderAlphaPlot <- reactive({
          myList <- get("Order")
          if(input$rawrareAlpha == 'Raw Data'){
            myData <- myList$diversityResults
          }
          else{
            myData <- myList$diversityResultsRarefy
          }
          
          whichDiv <- input$divMeasure
          myData = as.data.frame(myData)
          myData[[whichDiv]] <- as.numeric(myData[[whichDiv]])
          dataMedian <- summarise(group_by(myData, treatment), 
                                  MD = round(median(as.numeric(.data[[whichDiv]])), 2))
          yLabel <- names(diversityChoices)[grep(whichDiv, diversityChoices)]
          
          ggplot(myData,aes(x=treatment,y=.data[[whichDiv]], fill=treatment))+
            geom_boxplot(alpha=0.3)+theme_bw()+labs(x="Treatment",y=yLabel)+
            geom_text(data = dataMedian, aes(treatment, MD, label = MD), 
                      position = position_dodge(width=0.8),
                      size = 5, vjust = -0.5, hjust = 0.5)+
            theme(
              axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
              axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
              axis.text.x = element_text(size = 14),   
              axis.text.y = element_text(size = 14),   
              strip.text = element_text(size = 16), 
              legend.position="none"
            ) +
            scale_fill_brewer(palette="Accent") 
        })
        
        # create a boxplot of alpha diversity from family data only
        familyAlphaPlot <- reactive({
          myList <- get("Family")
          if(input$rawrareAlpha == 'Raw Data'){
            myData <- myList$diversityResults
          }
          else{
            myData <- myList$diversityResultsRarefy
          }
          
          whichDiv <- input$divMeasure
          myData = as.data.frame(myData)
          myData[[whichDiv]] <- as.numeric(myData[[whichDiv]])
          dataMedian <- summarise(group_by(myData, treatment), 
                                  MD = round(median(as.numeric(.data[[whichDiv]])), 2))
          yLabel <- names(diversityChoices)[grep(whichDiv, diversityChoices)]
          
          ggplot(myData,aes(x=treatment,y=.data[[whichDiv]], fill=treatment))+
            geom_boxplot(alpha=0.3)+theme_bw()+labs(x="Treatment",y=yLabel)+
            geom_text(data = dataMedian, aes(treatment, MD, label = MD), 
                      position = position_dodge(width=0.8),
                      size = 5, vjust = -0.5, hjust = 0.5)+
            theme(
              axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
              axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
              axis.text.x = element_text(size = 14),   
              axis.text.y = element_text(size = 14),   
              strip.text = element_text(size = 16), 
              legend.position="none"
            ) +
            scale_fill_brewer(palette="Accent") 
        })
      
        output$alphaPlots <- renderPlot({whichAlphaPlot()})
        output$phylumAlpha <- renderPlot({phylumAlphaPlot()})
        output$classAlpha <- renderPlot({classAlphaPlot()})
        output$orderAlpha <- renderPlot({orderAlphaPlot()})
        output$familyAlpha <- renderPlot({familyAlphaPlot()})
        
        # Functionality for downloading the the Alpha Diversity plot
        output$downloadAlpha <- downloadHandler(
          filename = function() {
            paste("alpha_diversity_plot.png", sep="")
          }, 
          content = function(file) {
            width <- 2000
            height <- 2000 
            png(file=file, width=width, height=height, res=150)
            plot(whichAlphaPlot())
            dev.off()
          }
        )
        
        if(nrow(treatments)==2){
        #create caption from results of t-test
          alphaStats <- reactive({
            taxa <- input$taxonAlphaTest
            myList <- get(taxa)
            if(input$rawrareAlpha=="Raw Data"){
              myData <- myList$diversityResults
            }
            else{
              myData <- myList$diversityResultsRarefy
            }
          
            whichDiv <- input$divMeasure
            
            # coerce "myData" into a dataframe so that the group_by function
            # is able to use it, group_by() only accepts tbl data type
            # Additionally, ensure data is numeric 
            myData = as.data.frame(myData)
            myData[[whichDiv]] <- as.numeric(myData[[whichDiv]])
            
            # use the t.test function on the data
            result <- t.test(myData[[whichDiv]]~myData[,1])
            
            # extract the degrees of freedom, p-value, and t-test statistic from
            # the result of the t.test function
            result_t <- round(as.numeric(result$statistic,digits=2))
            result_df <- round(as.numeric(result$parameter,digits=2))
            result_p <- round(result$p.value,digits=2)
            
            # format the p-value string
            pvalue <- ""
            if(result_p<0.01){
              pvalue <- "p < 0.01"
            }
            else{
              pvalue <- paste("p = ",result_p)
            }
            
            # format the t test stat and degrees of freedom string
            ttest <- sprintf("t = %d", result_t)
            degf <- sprintf("df = %d", result_df)
            
            # create a dataframe to show the results of Welch's two-sided t-test
            alphaStatsTable <- data.frame(
              Statistic = c("t-test statistic", "degrees of freedom", "p-value"), 
              Value = c(ttest, degf, pvalue)
            )
            return(alphaStatsTable)
          })
          
        # output the table to display Alpha diversity statistics in UI
        output$alphaStats <- renderTable({alphaStats()}, rownames=TRUE, 
                                         striped=TRUE, bordered=TRUE)
        
        }
        else{
          # this code will likely not be used by the University of Iowa because
          # the experiments run by undergraduates there only choose two 
          # experimental conditions
          
          #create ANOVA table and post-hoc comparisons
          alphaAnova <- reactive({
            taxa <- input$taxonAlphaTest
            myList <- get(taxa)
            
            if(input$rawrareAlpha=="Raw Data"){
              myData <- myList$diversityResults
            }
            else{
              myData <- myList$diversityResultsRarefy
            }
            
            anovaResult <- 0
            posthocResult <- 0
            anovaTables <- list()
            
            whichDiv <- input$divMeasure
            
            anovaResult <- anova(aov(myData[[whichDiv]]~myData[,1]))
            posthocResult <- tidy(TukeyHSD(aov(myData[[whichDiv]]~myData[,1])))
            
            posthocResult <- set(posthocResult,,1,NULL)
            
            rownames(anovaResult) <- c("Treatment","Residuals")
            anovaCaption <-"ANOVA Table"
            posthocCaption <-"Tukey's HSD"
            anovaTables <- list(anovaResult,anovaCaption,posthocResult,posthocCaption)
            
            #return list
            return(anovaTables)
          })
        
        #output ANOVA table and caption
        output$anovaCaption <- renderText({alphaAnova()[[2]]})
        output$alphaAnova <- renderTable({alphaAnova()[[1]]},rownames=TRUE, 
                                         striped=TRUE, bordered=TRUE)
        
        #output post-hoc table and caption
        output$posthocCaption <- renderText({alphaAnova()[[4]]})
        output$alphaPosthoc<-renderTable({alphaAnova()[[3]]})
      }
    }
      # close alpha diversity observe
    }
    # server side functions for beta diversity visualization
    else if (input$tabs == 'tab9'){
      if(exists("Phylum")){
        # get the taxonomic level the user selects
        whichTaxa <- reactive({
          taxa <- input$taxonBetaTest
          return(taxa)
        })
          
        # return the physeq object
        whichPhySeq <- reactive({
          taxa <- whichTaxa()
          myList <- get(taxa)
          if(input$rawrareBeta=='Raw Data'){
            myPhyseq <- myList$physeq
          }
          else{
            myPhyseq <- myList$physeqRare
          }
          return(myPhyseq)
        })
        
        # make ordination data using the ordinate() function from phyloseq
        whichOrdinationData <- reactive({
          myOrdData <- ordinate(whichPhySeq(), method=input$ordMethod, distance=input$distMeasure)
          return(myOrdData)
        })
        
        # plot the ordination, coord_fixed was removed to ensure that all axis
        # were of the same scaled
        whichScaledOrdinationPlot <- reactive ({
          if(input$samptreat=='Sample'){
            plot_ordination(whichPhySeq(), whichOrdinationData(), color = 'sample') +
              # change the axis title
              guides(color = guide_legend(title = "Sample"))+
              stat_ellipse(type='t')+
              theme_bw()+
              #coord_fixed()+
              theme(
                axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                axis.text.x = element_text(size = 14),   
                axis.text.y = element_text(size = 14),   
                legend.title = element_text(size = 16, , face="bold"),  
                legend.text = element_text(size = 14),
              )
          } else if (input$samptreat=='Treatment'){
            plot_ordination(whichPhySeq(), whichOrdinationData(), color = 'treatment') +
              guides(color = guide_legend(title = "Treatment"))+
              stat_ellipse(type='t')+
              theme_bw()+
              #coord_fixed()+
              theme(
                axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                axis.text.x = element_text(size = 14),   
                axis.text.y = element_text(size = 14),   
                legend.title = element_text(size = 16, , face="bold"),  
                legend.text = element_text(size = 14),
              )
          }
        })
        
        # plot the ordination for specifically phylum data
        phylumScaledOrdinationPlot <- reactive ({
          myList <- get("Phylum")
          if(input$rawrareBeta=='Raw Data'){
            myPhyseq <- myList$physeq
          }
          else{
            myPhyseq <- myList$physeqRare
          }
          
          myOrdData <- ordinate(myPhyseq, method=input$ordMethod, distance=input$distMeasure)
          
          if(input$samptreat=='Sample'){
            plot_ordination(myPhyseq, myOrdData, color = 'sample') +
              # change the axis title
              guides(color = guide_legend(title = "Sample"))+
              stat_ellipse(type='t')+
              theme_bw()+
              #coord_fixed()+
              theme(
                axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                axis.text.x = element_text(size = 14),   
                axis.text.y = element_text(size = 14),   
                legend.title = element_text(size = 16, , face="bold"),  
                legend.text = element_text(size = 14),
              )
          } else if (input$samptreat=='Treatment'){
            plot_ordination(myPhyseq, myOrdData, color = 'treatment') +
              guides(color = guide_legend(title = "Treatment"))+
              stat_ellipse(type='t')+
              theme_bw()+
              #coord_fixed()+
              theme(
                axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                axis.text.x = element_text(size = 14),   
                axis.text.y = element_text(size = 14),   
                legend.title = element_text(size = 16, , face="bold"),  
                legend.text = element_text(size = 14),
              )
          }
        })
        # plot the ordination for specifically class data
        classScaledOrdinationPlot <- reactive ({
          myList <- get("Class")
          if(input$rawrareBeta=='Raw Data'){
            myPhyseq <- myList$physeq
          }
          else{
            myPhyseq <- myList$physeqRare
          }
          
          myOrdData <- ordinate(myPhyseq, method=input$ordMethod, distance=input$distMeasure)
          
          if(input$samptreat=='Sample'){
            plot_ordination(myPhyseq, myOrdData, color = 'sample') +
              # change the axis title
              guides(color = guide_legend(title = "Sample"))+
              stat_ellipse(type='t')+
              theme_bw()+
              #coord_fixed()+
              theme(
                axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                axis.text.x = element_text(size = 14),   
                axis.text.y = element_text(size = 14),   
                legend.title = element_text(size = 16, , face="bold"),  
                legend.text = element_text(size = 14),
              )
          } else if (input$samptreat=='Treatment'){
            plot_ordination(myPhyseq, myOrdData, color = 'treatment') +
              guides(color = guide_legend(title = "Treatment"))+
              stat_ellipse(type='t')+
              theme_bw()+
              #coord_fixed()+
              theme(
                axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                axis.text.x = element_text(size = 14),   
                axis.text.y = element_text(size = 14),   
                legend.title = element_text(size = 16, , face="bold"),  
                legend.text = element_text(size = 14),
              )
          }
        })
        # plot the ordination for specifically phylum data
        orderScaledOrdinationPlot <- reactive ({
          myList <- get("Order")
          if(input$rawrareBeta=='Raw Data'){
            myPhyseq <- myList$physeq
          }
          else{
            myPhyseq <- myList$physeqRare
          }
          
          myOrdData <- ordinate(myPhyseq, method=input$ordMethod, distance=input$distMeasure)
          
          if(input$samptreat=='Sample'){
            plot_ordination(myPhyseq, myOrdData, color = 'sample') +
              # change the axis title
              guides(color = guide_legend(title = "Sample"))+
              stat_ellipse(type='t')+
              theme_bw()+
              #coord_fixed()+
              theme(
                axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                axis.text.x = element_text(size = 14),   
                axis.text.y = element_text(size = 14),   
                legend.title = element_text(size = 16, , face="bold"),  
                legend.text = element_text(size = 14),
              )
          } else if (input$samptreat=='Treatment'){
            plot_ordination(myPhyseq, myOrdData, color = 'treatment') +
              guides(color = guide_legend(title = "Treatment"))+
              stat_ellipse(type='t')+
              theme_bw()+
              #coord_fixed()+
              theme(
                axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                axis.text.x = element_text(size = 14),   
                axis.text.y = element_text(size = 14),   
                legend.title = element_text(size = 16, , face="bold"),  
                legend.text = element_text(size = 14),
              )
          }
        })
        
        # plot the ordination for specifically phylum data
        familyScaledOrdinationPlot <- reactive ({
          myList <- get("Family")
          if(input$rawrareBeta=='Raw Data'){
            myPhyseq <- myList$physeq
          }
          else{
            myPhyseq <- myList$physeqRare
          }
          
          myOrdData <- ordinate(myPhyseq, method=input$ordMethod, distance=input$distMeasure)
          
          if(input$samptreat=='Sample'){
            plot_ordination(myPhyseq, myOrdData, color = 'sample') +
              # change the axis title
              guides(color = guide_legend(title = "Sample"))+
              stat_ellipse(type='t')+
              theme_bw()+
              #coord_fixed()+
              theme(
                axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                axis.text.x = element_text(size = 14),   
                axis.text.y = element_text(size = 14),   
                legend.title = element_text(size = 16, , face="bold"),  
                legend.text = element_text(size = 14),
              )
          } else if (input$samptreat=='Treatment'){
            plot_ordination(myPhyseq, myOrdData, color = 'treatment') +
              guides(color = guide_legend(title = "Treatment"))+
              stat_ellipse(type='t')+
              theme_bw()+
              #coord_fixed()+
              theme(
                axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                axis.text.x = element_text(size = 14),   
                axis.text.y = element_text(size = 14),   
                legend.title = element_text(size = 16, , face="bold"),  
                legend.text = element_text(size = 14),
              )
          }
        })
        
        # plot the ordination without scaling all axes to be equal
        # this may cause some oddly sized plots based on the input data
        whichUnscaledOrdinationPlot <- reactive ({
          if(input$samptreat=='Sample'){
            plot_ordination(whichPhySeq(), whichOrdinationData(), color = 'sample') +
              # change the axis title
              guides(color = guide_legend(title = "Sample"))+
              stat_ellipse(type='t')+
              theme_bw()+
              coord_fixed()+ # coord-fix may cause some oddly size plots 
              theme(
                axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                axis.text.x = element_text(size = 14),   
                axis.text.y = element_text(size = 14),   
                legend.title = element_text(size = 16, , face="bold"),  
                legend.text = element_text(size = 14),
              )
          } else if (input$samptreat=='Treatment'){
            plot_ordination(whichPhySeq(), whichOrdinationData(), color = 'treatment') +
              guides(color = guide_legend(title = "Treatment"))+
              stat_ellipse(type='t')+
              theme_bw()+
              coord_fixed()+
              theme(
                axis.title.x = element_text(size = 16, margin = margin(t = 10), face="bold"),  
                axis.title.y = element_text(size = 16, margin = margin(r = 10), face="bold"),  
                axis.text.x = element_text(size = 14),   
                axis.text.y = element_text(size = 14),   
                legend.title = element_text(size = 16, , face="bold"),  
                legend.text = element_text(size = 14),
              )
          }
        })
        
        
        # make the ordination caption which displays the NMDS stress value
        whichOrdinationCaption <- reactive({
          ordData <- whichOrdinationData()
          if(input$ordMethod=='NMDS'){
            ordCap <- paste("<p class='text-light'>NMDS results: stress =", round(ordData$stress, digits=4), "</p>")
            ordCap <- toString(ordCap)
          }
          # may need to include an else statement here for other ordination methods
        })
        
        # plot the scaled and unscaled ordination plots in the UI with the caption
        output$scaledOrdPlot <- renderPlot({whichScaledOrdinationPlot()})
        output$unscaledOrdPlot <- renderPlot({whichUnscaledOrdinationPlot()})
        output$phylumBeta <- renderPlot({phylumScaledOrdinationPlot()})
        output$classBeta <- renderPlot({classScaledOrdinationPlot()})
        output$orderBeta <- renderPlot({orderScaledOrdinationPlot()})
        output$familyBeta <- renderPlot({familyScaledOrdinationPlot()})
        output$ordinationCaption <-renderUI({HTML(whichOrdinationCaption())})
        
        # Functionality for downloading the scaled Beta Diversity plot
        output$downloadBetaScaled <- downloadHandler(
          filename = function() {
            paste("scaled_beta_diversity_plot.png", sep="")
          }, 
          content = function(file) {
            width <- 2000
            height <- 2000 
            png(file=file, width=width, height=height, res=150)
            plot(whichScaledOrdinationPlot())
            dev.off()
          }
        )
        
        # Functionality for downloading the unscaled Beta Diversity plot
        output$downloadBetaUnscaled <- downloadHandler(
          filename = function() {
            paste("unscaled_beta_diversity_plot.png", sep="")
          }, 
          content = function(file) {
            width <- 2000
            height <- 2000 
            png(file=file, width=width, height=height, res=150)
            plot(whichUnscaledOrdinationPlot())
            dev.off()
          }
        )
        
        # function to create permutational multivariate analysis of variance (permanova) 
        # data for the beta diversity tab. 
        # Returns: a table containing adonis results
        whichPermanova <- reactive({
          taxa <- input$taxonBetaTest
          myList <- get(taxa)
          if(input$rawrareBeta == 'Raw Data'){
            permMethod <- myList$OTU_t
          }
          else{
            permMethod <- myList$otuRarefy_t
          }
          
          whichDist <- input$distMeasure
          # use adonis2 function to analyze the variance among distance matrices
          permaOut <- adonis2(permMethod ~ treatment, data = metaGlobal, method=whichDist)
          
          # format the permanova results as a dataframe
          permanovaResults <- as.data.frame(permaOut)
          # update the first rowname
          rownames(permanovaResults)[1] <- "Treatment"
          return(permanovaResults)
        })
        output$betaPermanova <- renderTable({whichPermanova()}, rownames=TRUE,
                                            striped=TRUE, bordered=TRUE)
        }
    }
    # close the observe event
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
