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

salseq = phyloseq_to_deseq2(phy, ~ bin_sal)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(salseq), 1, gm_mean)
salseq = estimateSizeFactors(salseq, geoMeans = geoMeans)

dds <- salseq

design(dds) <- formula(~ 1 + Salinity + Temp + Salinity:Temp)
dds <- DESeq(dds)

res <- results(dds)
head(res)
resultsNames(dds) 
## [1] "Intercept"     "Salinity"      "Temp"          "Salinity.Temp"

#get the model matrix

mod_mat <- model.matrix(design(dds), colData(dds))

# define coefficient vectors for each group





# Investigate test results table


res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phy)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Let’s look at just the OTUs that were significantly present. First, cleaning up the table a little for legibility.
posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
head(posigtab)

write.csv(as.data.frame(posigtab), 
          file="data/phyloseq/desq2_Temp.csv")
