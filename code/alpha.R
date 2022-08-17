# Alpha Diversity Graphics

# Reference: https://joey711.github.io/phyloseq/plot_richness-examples.html




library(phyloseq)
library(cluster)
library(ggplot2)
library(plyr)
library(tidyverse)

#ggplot theming
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)}

# Remove any OTUs that are not represented in the dataset

GP <- prune_species(speciesSums(phy) > 0, phy)

# Plot richness

plot_richness(GP)

# Plot just richness figures of interest
plot_richness(GP, measures=c("Chao1", "Shannon"))

# Group X axis

plot_richness(GP, x="Year", measures=c("Chao1", "Shannon"))
