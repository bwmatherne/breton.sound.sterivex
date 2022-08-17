# Plot Tree

# Reference: https://joey711.github.io/phyloseq/plot_tree-examples.html




library(phyloseq)
library(cluster)
library(ggplot2)
library(plyr)
library(tidyverse)




# Vibrio Only Tree
GP.vib <- subset_taxa(phy, Family=="f__Vibrionaceae")
vib<- plot_tree(GP.vib, color="Location", shape="Genus", label.tips="Species", size="abundance", plot.margin=0.6, ladderize=TRUE)

ggsave("data/phyloseq/vibrio_tree.png", vib,  width = 14, height = 10, dpi = 300)
