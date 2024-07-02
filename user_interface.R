
library(BiocManager)
library(shiny)
library(dyplr)
library(tidyverse)
library(ggplot2)
library(vegan)
library(data.table)
library(readr)
library(phyloseq)
library(broom)
library(shinycssloaders)
library(tibble)

# user interface

ui <- fluidPage(
  titlePanel("Bean Beetle Microbiome Analysis"), 
  tabsetPanel(id = "tabs", 
              tabPanel(value = "tab1", title = "Home", 
                       mainPanel(
                         h2("Welcome to the Bean Beetle Microbiome Analysis App"),
                         h4("The Bean Beetle Microbiome Project is a research/teachinf collaboration of institutions across the US that is studying the microbiome of", em("Callosobruchus maculatus"), "in research experiences (CUREs)."),
                         p("This app is designed to lead students through the community analysis of level-5 (family leve) datasets produced in DNA subway")
                         p("Before proceeding, students should ensure the level 5 file is formatted correctly such that all unidentified taxa, chloroplasts, and mitochondria have been removed. The level 5 file should also be formatted with the first column as taxa and the subsequent columns as samples with unique sample identifiers.")
                         p("Students should also prepare a metadata file in the first column as samples and the second column as treatments.", strong("Both files should be in .csv format.")), 
                         p("Additionally, this app should be capable for community analysis of any level 5 data")
                         p("")
                         
                         p("This app was reconfigured for the University of Iowa MICR 2158 course based on the source materials from Huang et al., 2022 by Elizabeth Elias, an undergraduate student at the University of Iowa under the guidance of Dr. Regina Mcgrane, Department of Microbiology and Immunology, University of Iowa")
                         p("The original app can be found using this link:", tags$a(href = https://beanbeetles.shinyapps.io/BeanBeetleMicrobiome/))
                         p("For more information on this CURE project, please visit:", tags$a(href = https://www.beanbeetles.org/microbiome/the-bean-beetle-microbiome-project/))
                        )
                    ),
              tabPanel(value = "tab2", title = "Home"
              )
  
  
  theme = bs_theme(preset = "slate")
  # close UI
)