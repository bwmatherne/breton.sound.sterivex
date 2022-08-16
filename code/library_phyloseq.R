################################################################################
#
# phyloseq_dependencies.R
#
# R script to help pre-install any difficult dependencies used with phyloseq
#
#									
#									
#
# Reference: https://uw-madison-microbiome-hub.github.io/Microbiome_analysis_in-_R/
#
################################################################################

# Clear the Environment
rm(list = ls())

# Install Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


library(devtools) # Make package development easier by providing R functions that simplify and expedite common tasks

install.packages("remotes")

library(phyloseq) # The phyloseq package seeks to address issues with multiple microbiome analysis packages by providing a set of functions that internally manage the organizing, linking, storing, and analyzing of phylogenetic sequencing data. In general, this package is used for UniFrac analyses.



# Load Packages
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

library(qiime2R) # A package for importing qiime artifacts into an R session

library(MicrobeR) # Data visualization

library(microbiome) # Data analysis and visualization

library(microbiomeSeq) # Data analysis and visualization ---> Not actively maintained and used KMDA dependency that must be installed using remotes::install_github("cran/KMDA")

library(pander) # provide a minimal and easy tool for rendering R objects into Pandoc's markdown

library(ranacapa) # Data analysis 

library(grid) # support data visualization

library(gridExtra)  # support data visualization

library(knitr) # Provides a general-purpose tool for dynamic report generation in R using Literate Programming techniques.

library(png) # Figure download

library(ggdendro) #set of tools for dendrograms and tree plots using 'ggplot2'

library(ggpubr) # publication quality figures, based on ggplot2

library(RColorBrewer) # nice color options

library(microbiomeutilities) # some utility tools 