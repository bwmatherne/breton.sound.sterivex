# DESeq2

# References: 
#   https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html 
#   https://joey711.github.io/phyloseq-extensions/DESeq2.html


rm(list = ls())

library(phyloseq)
library(cluster)
library(ggplot2)
library(plyr)
library(tidyverse)
library(DESeq2)

phy<-readRDS('data/phyloseq/phyloseq.RDS')

phy


# Convert data to DESeq2 format

#   Remove any samples with less than 500 reads

phy <- prune_samples(sample_sums(phy) > 500, phy)

phy # Shows that all samples were retained.
    # Output:
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 13369 taxa and 99 samples ]
## sample_data() Sample Data:       [ 99 samples by 15 sample variables ]
## tax_table()   Taxonomy Table:    [ 13369 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 13369 tips and 13358 internal nodes ]

salseq = phyloseq_to_deseq2(phy, ~ bin_sal)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(salseq), 1, gm_mean)
salseq = estimateSizeFactors(salseq, geoMeans = geoMeans)
salseq = DESeq(salseq, fitType="local")

    # Note: The default multiple-inference correction is Benjamini-Hochberg, and occurs within the DESeq function.

# Investigate test results table

res = results(salseq)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Let’s look at just the OTUs that were significantly present. First, cleaning up the table a little for legibility.
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
head(posigtab)

# Plot Results
# Here is a bar plot showing the log2-fold-change, showing Genus and Phylum. Uses some ggplot2 commands.

library("ggplot2")
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
plot01 <- ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave("data/phyloseq/deseq2_log2fold_01.png", plot01,  width = 14, height = 10, dpi = 300)


# Run analyses using settings from joey711 example

salseq2 = DESeq(salseq, test="Wald", fitType="parametric")
res2 = results(salseq2, cooksCutoff = FALSE)
alpha = 0.01
sigtab2 = res2[which(res2$padj < alpha),]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(phy)[rownames(sigtab2), ], "matrix"))
head(sigtab2)
dim(sigtab2)

# Let's look at the OTUs that were significantly different between the locations. The following makes a nice ggplot2 summary of the results.

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab2$log2FoldChange, sigtab2$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Phylum = factor(as.character(sigtab2$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab2$log2FoldChange, sigtab2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab2$Genus = factor(as.character(sigtab2$Genus), levels=names(x))
plot02 <- ggplot(sigtab2, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave("data/phyloseq/deseq2_log2fold_02.png", plot02,  width = 14, height = 10, dpi = 300)


# NOTE: plot01 is much easier to read and also removes the NA Genus values

# Rerun first example on "Location" since binning created three salinity groups

rm(list = ls())

phy<-readRDS('data/phyloseq/phyloseq.RDS')

phy


# Convert data to DESeq2 format

#   Remove any samples with less than 500 reads

phy <- prune_samples(sample_sums(phy) > 500, phy)

phy # Shows that all samples were retained.
# Output:
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 13369 taxa and 99 samples ]
## sample_data() Sample Data:       [ 99 samples by 15 sample variables ]
## tax_table()   Taxonomy Table:    [ 13369 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 13369 tips and 13358 internal nodes ]



salseq = phyloseq_to_deseq2(phy, ~ Location)

# Need to set a reference level in Location
salseq$Location <- relevel(salseq$Location, ref = "Inshore")

# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(salseq), 1, gm_mean)
salseq = estimateSizeFactors(salseq, geoMeans = geoMeans)
salseq = DESeq(salseq, fitType="local")

# Note: The default multiple-inference correction is Benjamini-Hochberg, and occurs within the DESeq function.

# Investigate test results table

res = results(salseq)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Let’s look at just the OTUs that were significantly present. First, cleaning up the table a little for legibility.
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
head(posigtab)

# Plot Results
# Here is a bar plot showing the log2-fold-change, showing Genus and Phylum. Uses some ggplot2 commands.

library("ggplot2")
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
plot03 <- ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave("data/phyloseq/deseq2_log2fold_03.png", plot03,  width = 14, height = 10, dpi = 300)


# CONTRASTS - https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts

