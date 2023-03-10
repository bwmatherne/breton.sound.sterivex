# Negative Binomial in Microbiome Differential Abundance Testing

# References: 
# https://bioconductor.riken.jp/packages/2.14/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html

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



salseq = phyloseq_to_deseq2(phy, ~ Location)

# Need to set a reference level in Location
salseq$Location <- relevel(salseq$Location, ref = "Inshore")




# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(salseq), 1, gm_mean)
salseq = estimateSizeFactors(salseq, geoMeans = geoMeans)

salseq = DESeq(salseq, test="Wald", fitType="parametric")
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
plot04 <- ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
plot04
ggsave("data/phyloseq/deseq2_log2fold_04.png", plot04,  width = 14, height = 10, dpi = 300)

write.csv(as.data.frame(posigtab),
          file="data/phyloseq/negative_binomial_results.csv")

# CONTRASTS - https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts
# Contrasts 2: https://github.com/tavareshugo/tutorial_DESeq2_contrasts
    # 1 factor (bin_sal) and 3 levels (inshore, coast.hi, coast.low)
# Multifactor 

    # We create a copy of the DESeqDataSet, so that we can rerun the analysis using a multi-factor design.

ddsMF <- salseq
levels(ddsMF$bin_sal)
colData(ddsMF)

# We can account for the different salinity bins, and get a clearer picture of the differences attributable to the location As location is the variable of interest, we put it at the end of the formula. Thus the results function will by default pull the locations results unless contrast or name arguments are specified.

# Then we can re-run DESeq:
  
design(ddsMF) <- formula(~ 1 + bin_sal)
ddsMF <- DESeq(ddsMF)

resMF <- results(ddsMF)
head(resMF)
resultsNames(ddsMF)

## [1] "Intercept"                     "bin_sal_coast.low_vs_coast.hi" "bin_sal_inshore_vs_coast.hi" 
    # NOTE: Intercept = coast.hi as the "reference"
    # Note: Consider a relevel to set another group as the reference

# First method 

res1_lowVhi <- results(ddsMF, contrast = list("bin_sal_coast.low_vs_coast.hi"))
res1_inshoreVhi <- results(ddsMF, contrast = list("bin_sal_inshore_vs_coast.hi"))
res1_lowVinshore <- results(ddsMF, contrast = list("bin_sal_coast.low_vs_coast.hi",
                                                   "bin_sal_inshore_vs_coast.hi"))

# Alternate method - better control with downstream outputs
  # define the model matrix
mod_mat <- model.matrix(design(ddsMF), colData(ddsMF))
mod_mat

  # calculate coefficient vectors for each group
coast.hi <- colMeans(mod_mat[ddsMF$bin_sal == "coast.hi", ])
coast.low <- colMeans(mod_mat[ddsMF$bin_sal == "coast.low", ])
inshore <- colMeans(mod_mat[ddsMF$bin_sal == "inshore", ])

  # Define any contrast wanted
  # Obtain results for each pairwise contrast
res2_lowVhi <- results(ddsMF, contrast = coast.low - coast.hi)
res2_inshoreVhi <- results(ddsMF, contrast = inshore - coast.hi)
res2_lowVinshore <- results(ddsMF, contrast = coast.low - inshore)
  
# plot the results from the two approaches to check that they are identical
    # This step is just an example to show that both methods give the same type of results
plot(res1_inshoreVhi$log2FoldChange, res2_inshoreVhi$log2FoldChange)
plot(res1_lowVhi$log2FoldChange, res2_lowVhi$log2FoldChange)
plot(res1_lowVinshore$log2FoldChange, res2_lowVinshore$log2FoldChange)

# Look for species that differ in coast.hi as compared to coast.low and inshore

combo.low <- colMeans(mod_mat[ddsMF$bin_sal %in% c("inshore", "coast.low"),])
# contrast of interest
res2_combo.low <- results(ddsMF, contrast = combo.low - coast.hi)

################################3
# Investigate test results table
  # Output = deseq2_combo.low.png
res2_combo.low = res2_combo.low[order(res2_combo.low$padj, na.last=NA), ]
alpha = 0.01
sigtab1 = res2_combo.low[(res2_combo.low$padj < alpha), ]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(phy)[rownames(sigtab1), ], "matrix"))
head(sigtab1)

# Let’s look at just the OTUs that were significantly present. First, cleaning up the table a little for legibility.
posigtab1 = sigtab1[sigtab1[, "log2FoldChange"] > 0, ]
posigtab1 = posigtab1[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
head(posigtab1)

# Plot Results
# Here is a bar plot showing the log2-fold-change, showing Genus and Phylum. Uses some ggplot2 commands.

library("ggplot2")
theme_set(theme_bw())
sigtabgen1 = subset(sigtab1, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen1$log2FoldChange, sigtabgen1$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen1$Phylum = factor(as.character(sigtabgen1$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen1$log2FoldChange, sigtabgen1$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen1$Genus = factor(as.character(sigtabgen1$Genus), levels=names(x))
plot.low <- ggplot(sigtabgen1, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
plot.low
ggsave("data/phyloseq/deseq2_combo.low.png", plot.low,  width = 14, height = 10, dpi = 300)

##############
# Output = deseq2_coastal.low.hi.png
res2_lowVhi= res2_lowVhi[order(res2_lowVhi$padj, na.last=NA), ]
alpha = 0.01
sigtab2 = res2_lowVhi[(res2_lowVhi$padj < alpha), ]
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(phy)[rownames(sigtab2), ], "matrix"))
head(sigtab2)

# Let’s look at just the OTUs that were significantly present. First, cleaning up the table a little for legibility.
posigtab2 = sigtab2[sigtab2[, "log2FoldChange"] > 0, ]
posigtab2 = posigtab2[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
head(posigtab2)

# Plot Results
# Here is a bar plot showing the log2-fold-change, showing Genus and Phylum. Uses some ggplot2 commands.

library("ggplot2")
theme_set(theme_bw())
sigtabgen2 = subset(sigtab2, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen2$log2FoldChange, sigtabgen2$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen2$Phylum = factor(as.character(sigtabgen2$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen2$log2FoldChange, sigtabgen2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen2$Genus = factor(as.character(sigtabgen2$Genus), levels=names(x))
plot.low.hi <- ggplot(sigtabgen2, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
plot.low.hi
ggsave("data/phyloseq/deseq2_coastal.low.hi.png", plot.low.hi,  width = 14, height = 10, dpi = 300)



##############
# Output = deseq2_inshore.coast.hi.png
res2_inshoreVhi= res2_inshoreVhi[order(res2_inshoreVhi$padj, na.last=NA), ]
alpha = 0.01
sigtab3 = res2_inshoreVhi[(res2_inshoreVhi$padj < alpha), ]
sigtab3 = cbind(as(sigtab3, "data.frame"), as(tax_table(phy)[rownames(sigtab3), ], "matrix"))
head(sigtab3)

# Let’s look at just the OTUs that were significantly present. First, cleaning up the table a little for legibility.
posigtab3 = sigtab3[sigtab3[, "log2FoldChange"] > 0, ]
posigtab3 = posigtab3[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
head(posigtab3)

# Plot Results
# Here is a bar plot showing the log2-fold-change, showing Genus and Phylum. Uses some ggplot2 commands.

library("ggplot2")
theme_set(theme_bw())
sigtabgen3 = subset(sigtab3, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen3$log2FoldChange, sigtabgen3$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen3$Phylum = factor(as.character(sigtabgen3$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen3$log2FoldChange, sigtabgen3$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen3$Genus = factor(as.character(sigtabgen3$Genus), levels=names(x))
plot.inshore.ch <- ggplot(sigtabgen3, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
plot.inshore.ch
ggsave("data/phyloseq/deseq2_inshore.coast.hi.png", plot.inshore.ch,  width = 14, height = 10, dpi = 300)

##############
# Output = deseq2_coast.low.inshore.png
res2_lowVinshore = res2_lowVinshore[order(res2_lowVinshore$padj, na.last=NA), ]
alpha = 0.01
sigtab4 = res2_lowVinshore[(res2_lowVinshore$padj < alpha), ]
sigtab4 = cbind(as(sigtab4, "data.frame"), as(tax_table(phy)[rownames(sigtab4), ], "matrix"))
head(sigtab4)

# Let’s look at just the OTUs that were significantly present. First, cleaning up the table a little for legibility.
posigtab4 = sigtab4[sigtab4[, "log2FoldChange"] > 0, ]
posigtab4 = posigtab4[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
head(posigtab4)

# Plot Results
# Here is a bar plot showing the log2-fold-change, showing Genus and Phylum. Uses some ggplot2 commands.

library("ggplot2")
theme_set(theme_bw())
sigtabgen4 = subset(sigtab4, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen4$log2FoldChange, sigtabgen4$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen4$Phylum = factor(as.character(sigtabgen4$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen4$log2FoldChange, sigtabgen4$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen4$Genus = factor(as.character(sigtabgen4$Genus), levels=names(x))
plot.cl.inshore <- ggplot(sigtabgen4, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
plot.cl.inshore
ggsave("data/phyloseq/deseq2_coast.low.inshore.png", plot.cl.inshore,  width = 14, height = 10, dpi = 300)


##
##
##

# continuous covariate

colData(ddsMF)

# We can account for the different salinity bins, and get a clearer picture of the differences attributable to the location As location is the variable of interest, we put it at the end of the formula. Thus the results function will by default pull the locations results unless contrast or name arguments are specified.

# Then we can re-run DESeq:
ddsSAL <- salseq
design(ddsSAL) <- formula(~ 1 + Salinity)
ddsSAL <- DESeq(ddsSAL)

resSAL <- results(ddsSAL)
head(resSAL)
resultsNames(ddsSAL) # Need to continue making figure

# Investigate test results table


resSAL = resSAL[order(resSAL$padj, na.last=NA), ]
alpha = 0.01
sigtab5 = resSAL[(resSAL$padj < alpha), ]
sigtab5 = cbind(as(sigtab5, "data.frame"), as(tax_table(phy)[rownames(sigtab5), ], "matrix"))
head(sigtab5)

# Let’s look at just the OTUs that were significantly present. First, cleaning up the table a little for legibility.
posigtab5 = sigtab5[sigtab5[, "log2FoldChange"] > 0, ]
posigtab5 = posigtab5[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
head(posigtab5)

write.csv(as.data.frame(posigtab5), 
          file="data/phyloseq/desq2_Salinity.csv")
