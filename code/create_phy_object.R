# Generate phyloseq object before moving forward with analysis

# Clear the Environment
rm(list = ls())

# 4.1. Load Packages

library(DECIPHER) # This package will help in importing, maintaining, analyzing, manipulating, and exporting a massive amount of sequences.

library(ape) # Analyses of Phylogenetics and Evolution package. Required for tree calculations to be used with phyloseq

library(DESeq2) # This package will help analyze "differential expression" in the microbiota alongside phyloseq

library(ggplot2) # Graphing package used in phyloseq. To edit the default setting of a plot, you need to use functions in this package.

library(phyloseq) # The phyloseq package seeks to address issues with multiple microbiome analysis packages by providing a set of functions that internally manage the organizing, linking, storing, and analyzing of phylogenetic sequencing data. In general, this package is used for UniFrac analyses.

library(plotly) # A package to create interactive web graphics of use in 3D plots

library(vegan) # The vegan package provides tools for descriptive community ecology. It has most basic functions of diversity analysis, community ordination and dissimilarity analysis. In general, this package is used for Bray-Curtis and Jaccard analyses.

library(philr) # This package provides functions for the analysis of compositional data 

library(tidyverse) # This package is designed to make it easy to install and load multiple 'tidyverse' packages in a single step

library(adespatial) # Tools for the multiscale spatial analysis of multivariate data

library(devtools) # Make package development easier by providing R functions that simplify and expedite common tasks

library(qiime2R) # A package for importing qiime artifacts into an R session

library(MicrobeR) # Data visualization

library(microbiome) # Data analysis and visualization

library(microbiomeSeq) # Data analysis and visualization

library("pander") # provide a minimal and easy tool for rendering R objects into Pandoc's markdown

library(ranacapa) # Data analysis 

library(grid) # support data visualization

library(gridExtra)  # support data visualization

library(knitr) # Provides a general-purpose tool for dynamic report generation in R using Literate Programming techniques.

library(png) # Figure download

library("ggdendro") #set of tools for dendrograms and tree plots using 'ggplot2'

library(ggpubr) # publication quality figures, based on ggplot2

library(RColorBrewer) # nice color options

library(microbiomeutilities) # some utility tools 

#4.2 Load Data
#Convert qiime artifacts directly to phyloseq

# Importing ASVs abundance file
ASVs <- read_qza("data/process/bs.table.qza")

#Importing metadata
metadata <- read.table("data/process/bs.metadata.tsv", sep='\t', header=T, row.names=1, comment="")
metadata <- metadata[-1,] # remove the second line that specifies the data type
metadata$Month <- fct_relevel(metadata$Month, c("1","2","3","4","5","6","8","9","10")) # Set month order for down stream plotting

# Importing tree
tree <- read_qza("data/process/rooted-tree.qza")
# Importing taxonomy
taxonomy <- read_qza("data/process/taxonomy.qza")
tax_table <- do.call(rbind, strsplit(as.character(taxonomy$data$Taxon), "; "))
colnames(tax_table) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
rownames(tax_table) <- taxonomy$data$Feature.ID

# Creating phyloseq object
phy <- phyloseq(
  otu_table(ASVs$data, taxa_are_rows = T),
  phy_tree(tree$data),
  tax_table(tax_table),
  sample_data(metadata))