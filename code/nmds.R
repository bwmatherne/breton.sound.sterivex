# NMDS
# Reference https://web.stanford.edu/class/bios221/Pune/Labs/Lab_phyloseq/Phyloseq_Lab.html#documentation

#Load packages
require("tidyverse")
library("ggplot2")
library("dplyr")
library("phyloseq")

theme_set(theme_bw())

phy<-readRDS('phy_data.RDS')

ordinate(phy, "NMDS", "wUniFrac") %>% 
  plot_ordination(phy, ., color = "Month", title = "weighted-UniFrac")


library(vegan)
library(mgcv)
data(dune)
summarise(dune)
names(dune)
data(data.env)

ord <- ordinate(phy, "NMDS", "wUniFrac")

ord.fit <- envfit(ord ~ Temp + Month, data=metadata, perm=999)

ord.fit
plot(ord.fit)

plot(ord, dis="site")

ordisurf(ord, Salinity, add=TRUE)
