
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
  core <- core[order(-core$abundance)]
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