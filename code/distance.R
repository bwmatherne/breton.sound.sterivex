# Distance Function Script

# Preliminary filtering or preprocessing may be required

library(phyloseq)
library(ggplot2)
library(plyr)

theme_set(theme_bw())


# Remove the OTUs that included all unassigned sequences ("-1")
phy2 <- subset_taxa(phy, Genus != "-1")

# The available distance methods coded in distance

dist_methods <- unlist(distanceMethodList)
print(dist_methods)

# Loop through each distance method, save each plot to a list, called plist
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(phy2, method=i)
  # Calculate ordination
  iMDS  <- ordinate(phy2, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(phy2, iMDS, color="Month", shape="Location")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}
