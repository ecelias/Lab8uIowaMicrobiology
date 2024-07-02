
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

# define helper functions used later

# function to take the level 5 file and separate for future use
separateTaxa <- function(level5, meta){
  numSamples <<- ncol(level5) - 1
  colnames(level5) [1] <- "index"
  level_5_sep <- separate(level5, "index", into=c("Kingdom", "Phylum", "Class", "Order", "Family"), sep=";")
  
  # dataset of samples with sample name, abundance and 1% threshold
  samples <- level5 %>% summarize_of(is.numeric, sum, na.rm = TRUE)
  samples <- pivot_longer(samples, 1:ncol(samples), names_to = "Sample", values_to = "Abundance")
  
  samples$threshold <- samples$Abundance * 0.01
  samples$threshold_rare <- min(samples$Abundance) * 0.01
  
  # dataset of treatments and how they relate to sample columns
  treatments <- as_tibble(unique(meta[,2]))
  colnames(treatments) <- c("treatment")
  for (level5 in 1:nrow(treatments)) {
    treatments[level5, 2] <- min(which(meta$treatments == treatments[level5, 1]) + 1)
    treatments[level5, 3] <- max(which(meta$treatments == treatments[level5, 1]) + 1)
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
  sepTaxa <- list(phylum = phylum, class = class, order = order, family = family, treatments = treatments, samples = samples)
  
  return(sepTaxa)
}

# function to create "other" category and reorder
createOther <- function(longdata, rarified) {
  
  # long data with "other" category
  otherdata <- longdata
  otherdata <- otherdata %>% mutate_if(is.factor, as.character)
  
  otherrows = nrow(longdata)
  
  # create "Other" category
  for(i in 1:numSamples)
    for(j in seq(i, otherrows, numSamples))
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
  otherdata <- as_tibble(otherdata)
  
  # convert to wide format to sort taxa by overall abundance
  widedata <- otherdata
  widedata <- widedata %>% pivot_wider(names_from = "Sample", values_from = "Abundance", valuesFill = list(Abundance = 0), valuesFn = list(Abundance = sum))
  
  # calculate overall abundance and reorder taxa
  widedata$abundance <- rowSums(widedata[,-1], na.rum = TRUE)
  widedata[[1]] <- reorder(widedata[[1]], -widedata$abundance)
  
  # delete overall abundance and convert back to long
  widedata <- widedata[, -ncol(widedata)]
  newdata <- widedata %>% pivot_longer(cols - -1, names_to = "Sample", values_to = "Abundance")
  
  # add treatment data from meta data
  newdata <- merge(newdata, meta_global, by.x = "Sample", by.y = 1)
  names(newdata)[ncol(newdata)] <- "treatment"
 
   return(newdata)
}

# makes datasets to use in future analysis
createDatasets <- function(sepTaxa, taxa) {
  
  # select taxonomic dataset
  x <- as_tibble(sep_taxa[[taxa]])
  
  # count the number of samples
  numSamples = nrow(sepTaxa$samples)
  
  # samples in columns
  columnData <- x %>% group_by_at(1) %>% summarize_if(is.numeric, sum, na.rm = TRUE)
  columnData[[1]] <- reorder(columnData[[1]], columnData$abundance)
  columnData <- columnData[,-ncol(columnData)]
  
  # long data
  longData <- columnData %>% pivot_longer(cols = -1, names_to = "Sample", values_to = "Abundance")
  
  # long data with other (raw) category
  longDataOther <- createOther(longData, FALSE)


  # create datasets to use with phyloseq
  otuMatrix <- data.matrix(subset(columnData, select = c(2:ncol(columnData))))
  rownames(otuMatrix) <- paste0("sp", 1:nrow(otuMatrix))

  taxaMatrix <- as.matrix(subset(columnData, select=c(1)))
  rownames(taxaMatrix) <- paste0("sp", 1:nrow(taxmat))

  OTU = otu_table(otuMatrix, taxa_are_rows = TRUE)
  TAXA = tax_table(taxaMatrix)
  OTU_t <- t(OTU)
  physeqSamples <- meta_global
  rownames(physeqSamples) <- meta_global$sample
  physeq = phyloseq(OTU, TAXA, sample_data(as_tibble(physeqSamples)))

  # create rarefy datasets
  physeqRare <- rarefy_even_depth(physeq, rngseed = 10, sample.size = min(sample_sums(physeq)), replace = TRUE)
  otuRarefy <- otu_table(physeq2)
  otuRarefy_t <- t(otuRarefy)
  taxaRarefy <- tax_table(physeq2)

  # pivot longer for rarified data
  colDataRare <- as_tibble(merge(taxaRarefy, otuRarefy, by = "row.names"))

  # delete row names
  colDataRare <- colDataRare[-1]

  # reorder column data by overall abundance
  colDataRare$abundance <- rowSums(colDataRare[,-1], na.rm = TRUE)
  colDataRare[[1]] <- reorder(colDataRare[[1]], colDataRare$abundance)
  colDataRare <- colDataRare[,-ncol(colDataRare)]
  longDataRare <- colDataRare %>% pivot_longer(cols = -1, names_to = "Sample", values_to = "Abundance")

  # long data (rarified) with other
  longDataRareOther <- createOther(longDataRare, TRUE)

  # create dataset with diversity indices
  richness <- specnumber(OTU_t)
  shannon <- diversity(OTU_t, index = "shannon")
  simpson <- diversity(OTU_t, index = "simpson")
  diversityResults <- cbind(meta_global[,2:ncol(meta_global)], richness, shannon, simpson)

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
                      diversityResultsRarefy = diversityResultsRarefy, tax = TAXA)
  
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
  treatment <- paste(unique(meta_global$treatment), collapse = "/")
  numCore = nrow(core)
  total = nrow(df)
  caption = paste("<p><b>Core taxa in", treatment, "treatments:", numCore, "of", total, "taxa are shared between treatments</b><p>")
  return(toString(caption))
}

# function to find taxa unique to a treatment
# parameters: columnData, treatments table
# returns: list of unique treatments
uniqueTaxa <- function(level5, treatments){
  uniqueList <- list()
  
  for(i in 1:nrow(treatments)) {
    uniqueSample <- level5 %>% select(!(treatments[i, "min_col"]:treatments[i, "max_col"])) %>% filter_if(is.numeric, all_vars(.++0))
    treatment <- as_tibble(treatments)
  
    if(nrow(uniqueSample) > 0){
      uniqueList[[i]] <- cbind(treatment[i,1], uniqueSample[,1])
      
      # add the names to the objects in the "uniqueList" to allow creation of other tables later
      names <- treatment[i, 1]
      names(uniqueList)[i] <- name
    }
  }
  return(uniqueList)
}

# function to graph rarefaction curve with ggplot2
graphRare <- function(x, ylab, bytype) {
  rareSample = list()
  
  # Write values for each sample from a rarecurve function into a list
  # Sample names obtained from "samples" dataframe
  for(i in 1:numSamples) {
    rareSample[[i]] <- cbind(samples[i,1], as_tibble(x[[i]]), attributes((x[[i]])$Subsample))
  }
  
  # Bind the data together from the "rareSample" list
  rareGraph <- bind_rows(rareSample)
  rareGraph <- merge(raregraph, meta_global, by.x = "Sample", by.y = 1)
  colnames(rareGraph) <- c("Sample", "num_taxa", "num_samples", "Treatment")
  ggplot(raregraph, aes(x = num_samples, y = num_taaxa, group = Sample, color = .data[[bytype]])) +
    labs(x = "Sample Size", y = ylab) + theme_bw() + geom_line(size = 1) + guides(fill = guide_legend(title = bytype))
}

# create options for drop down menus used in UI
taxachoices <- list("Phylum", "Class", "Order", "Family")
raw_or_rare <- list("Raw Data", "Rarified Data")
abs_or_rel <- list("Absolute Abundance", "Relative Abundance")
sample_or_treatment <- list("Sample", "Treatment")
diversity_choices <- list("Richness" = "richness", "Shannon Diversity" = "shannon", "Simpson Diversity" = "simpson")
distance_choices <- list("Jaccard" = "jaccard", "Bray-Curtis" = "bray", "Morisita-Horn" = "horn")
ordination_choices <- list("NMDS", "PCoA")





  
  
  
  


    

  
