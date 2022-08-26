# Distance Function Script
# Multidimensional Scaling Plots (MDS)
# Reference: https://joey711.github.io/phyloseq/distance.html
# Preliminary filtering or preprocessing may be required
# WARNING the Loop step takes several hours to run without parallel computing


library(phyloseq)
library(ggplot2)
library(plyr)
library(tidyverse)

theme_set(theme_bw())

# Importing ASVs abundance file
ASVs <- read_qza("data/process/bs.table.qza")

#Importing metadata
metadata <- read_tsv("data/process/bin.metadata.tsv")
# Importing tree
tree <- read_qza("data/process/rooted-tree.qza")
# Importing taxonomy
taxonomy <- read_qza("data/process/taxonomy.qza")
tax_table <- parse_taxonomy(taxonomy$data)
tax_table <- tax_names(taxonomy)



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


# Remove the OTUs that included all unassigned sequences ("-1")
phy2 <- subset_taxa(phy, Genus != "-1")

# The available distance methods coded in distance

dist_methods <- unlist(distanceMethodList)
print(dist_methods)

# This is the user-defined method:
dist_methods["designdist"]

# Remove the user-defined distance
dist_methods = dist_methods[-which(dist_methods=="ANY")]

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

save(plist, file = "data/phyloseq/distance_results.RData")


#Combine results
# Shade according to Month of sample collection

df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"

# Properly order months

df$Month <- fct_relevel(df$Month, c("1","2","3","4","5","6","8","9","10"))

#Plot the figure
p = ggplot(df, aes(Axis.1, Axis.2, color=Month, shape=Location))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for Breton Sound dataset")
p

ggsave("data/phyloseq/MDS_Color_Month.png", p,  width = 14, height = 10, dpi = 300)

# Compare results
plist$wunifrac$data$Month <- fct_relevel(plist$wunifrac$data$Month, c("1","2","3","4","5","6","8","9","10"))
wunifrac <- print(plist[["wunifrac"]])
ggsave("data/phyloseq/wunifrac.png", wunifrac,  width = 14, height = 10, dpi = 300)


plist$binomial$data$Month <- fct_relevel(plist$binomial$data$Month, c("1","2","3","4","5","6","8","9","10"))
binomial <- print(plist[["binomial"]])
ggsave("data/phyloseq/binomial.png", binomial,  width = 14, height = 10, dpi = 300)

plist$dpcoa$data$Month <- fct_relevel(plist$dpcoa$data$Month, c("1","2","3","4","5","6","8","9","10"))
dpcoa <- print(plist[["dpcoa"]])
ggsave("data/phyloseq/dpcoa.png", dpcoa,  width = 14, height = 10, dpi = 300)

plist$altGower$data$Month <- fct_relevel(plist$altGower$data$Month, c("1","2","3","4","5","6","8","9","10"))
altGower <- print(plist[["altGower"]])
ggsave("data/phyloseq/altGower.png", altGower,  width = 14, height = 10, dpi = 300)

plist$unifrac$data$Month <- fct_relevel(plist$unifrac$data$Month, c("1","2","3","4","5","6","8","9","10"))
unifrac <- print(plist[["unifrac"]])
ggsave("data/phyloseq/unifrac.png", unifrac,  width = 14, height = 10, dpi = 300)
