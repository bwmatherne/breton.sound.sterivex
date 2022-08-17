# Gap Statistic Tutorial
# How many clusters are there in an MDS plot?
# Reference: https://joey711.github.io/phyloseq/gap-statistic.html
# Preliminary filtering or preprocessing may be required



library(phyloseq)
library(cluster)
library(ggplot2)
library(plyr)
library(tidyverse)

theme_set(theme_bw())


#Load data
phy <- readRDS("phy_data.rds")

# ordination
exord = ordinate(phy, method="MDS", distance="jsd")

# Compute Gap Statistic

pam1 = function(x, k){list(cluster = pam(x,k, cluster.only=TRUE))}
x = phyloseq:::scores.pcoa(exord, display="sites")
# gskmn = clusGap(x[, 1:2], FUN=kmeans, nstart=20, K.max = 6, B = 500)
gskmn = clusGap(x[, 1:2], FUN=pam1, K.max = 6, B = 50)
gskmn


#Wrapper-function add-on
gap_statistic_ordination = function(ord, FUNcluster, type="sites", K.max=6, axes=c(1:2), B=500, verbose=interactive(), ...){
  require("cluster")
  #   If "pam1" was chosen, use this internally defined call to pam
  if(FUNcluster == "pam1"){
    FUNcluster = function(x,k) list(cluster = pam(x, k, cluster.only=TRUE))     
  }
  # Use the scores function to get the ordination coordinates
  x = phyloseq:::scores.pcoa(ord, display=type)
  #   If axes not explicitly defined (NULL), then use all of them
  if(is.null(axes)){axes = 1:ncol(x)}
  #   Finally, perform, and return, the gap statistic calculation using cluster::clusGap  
  clusGap(x[, axes], FUN=FUNcluster, K.max=K.max, B=B, verbose=verbose, ...)
}
# Plot Results

plot_clusgap = function(clusgap, title="Gap Statistic calculation results"){
  require("ggplot2")
  gstab = data.frame(clusgap$Tab, k=1:nrow(clusgap$Tab))
  p = ggplot(gstab, aes(k, gap)) + geom_line() + geom_point(size=5)
  p = p + geom_errorbar(aes(ymax=gap+SE.sim, ymin=gap-SE.sim))
  p = p + ggtitle(title)
  return(p)
}

gs = gap_statistic_ordination(exord, "pam1", B=50, verbose=FALSE)
print(gs, method="Tibs2001SEmax")
plot_clusgap(gs)

plot(gs, main = "Gap statistic for the 'Breton Sound' data")
mtext("Looks like 6 clusters is best")  
