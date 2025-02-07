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
library(BiodiversityR)

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
                               layout_columns(
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
                             ),
                             # card for the rank abundance curve
                             card(
                               height = 650,
                               # card header defines the format of the card, additional
                               # card formats can be found on bootswatch
                               card_header(
                                 class = "bg-secondary mb-3",
                                 "Rank Abundance Curve"
                               ),
                               card_body(
                                 tags$i(p("A rank abundance curve or Whittaker plot represents the relative abundance of species within a community.
                                          by plotting species abundance against its rank order, providing a visualization of both species
                                          richness and species evenness. A steeper curve indicates that a community is comprised of a few dominant species
                                          while the others are relatively rare whereas a flatter curve suggests more more even distrubution of species",
                                          class="text-light")), 
                                 tags$hr(),
                                 plotOutput("rankabundancecurve", width = "100%", height = "auto")
                                 
                               )
                             )
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
                           card(
                             p("Core taxa are those taxa found in all samples.",class="text-light"),
                             # ensure the text is above the table
                             verticalLayout(htmlOutput("coreCaption"), tableOutput("coreTaxa")), 
                           )
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
                               withSpinner(color='#0dc5c1'), 
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
                               withSpinner(color='#0dc5c1'), 
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
                           card(
                             plotOutput('bargraph',width='100%', height='auto') %>%
                               withSpinner(color='#0dc5c1'),
                             full_screen = TRUE
                           ), 
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
                               withSpinner(color='#0dc5c1'), 
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
                                         withSpinner(color='#0dc5c1')
                             ), 
                             nav_panel("Side-by-Side", 
                                       card(
                                         layout_columns(
                                           card(
                                             h6("Phylum", class="text-light"), 
                                             plotOutput("phylumAlpha")  %>%
                                               withSpinner(color='#0dc5c1')
                                           ), 
                                           card(
                                             h6("Class", class="text-light"), 
                                             plotOutput("classAlpha")  %>%
                                               withSpinner(color='#0dc5c1')
                                           )
                                         ),
                                         layout_columns(
                                           card(
                                             h6("Order", class="text-light"), 
                                             plotOutput("orderAlpha")  %>%
                                               withSpinner(color='#0dc5c1'), 
                                           ), 
                                           card(
                                             h6("Family", class="text-light"), 
                                             plotOutput("familyAlpha")  %>%
                                               withSpinner(color='#0dc5c1')
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
                                         withSpinner(color='#0dc5c1')
                             ),
                             # panel with scaled plots 
                             nav_panel("Scaled", 
                                       plotOutput("scaledOrdPlot")  %>%
                                         withSpinner(color='#0dc5c1')
                             ), 
                             nav_panel("Side-by-Side", 
                                       card(
                                         layout_columns(
                                           card(
                                             h6("Phylum", class="text-light"), 
                                             plotOutput("phylumBeta")  %>%
                                               withSpinner(color='#0dc5c1')
                                           ), 
                                           card(
                                             h6("Class", class="text-light"), 
                                             plotOutput("classBeta")  %>%
                                               withSpinner(color='#0dc5c1')
                                           )
                                         ),
                                         layout_columns(
                                           card(
                                             h6("Order", class="text-light"), 
                                             plotOutput("orderBeta")  %>%
                                               withSpinner(color='#0dc5c1'), 
                                           ), 
                                           card(
                                             h6("Family", class="text-light"), 
                                             plotOutput("familyBeta")  %>%
                                               withSpinner(color='#0dc5c1')
                                           )
                                         )
                                       )
                             )
                           ), 
                           tags$hr(),
                           htmlOutput('ordinationCaption'),
                           downloadButton("downloadBetaScaled", "Download Unscaled Plot", class="btn-sm"),
                           downloadButton("downloadBetaUnscaled", "Download Scaled Plot", class="btn-sm"),
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