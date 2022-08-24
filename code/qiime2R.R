# qiime2R

# Reference https://rdrr.io/github/jbisanz/qiime2R/f/README.md 

#Load dependencies

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")


library(qiime2R)
library(phyloseq)


# Create phyloseq object using wrapper
physeq<-qza_to_phyloseq(
  features="data/process/bs.table.qza",
  tree="data/process/rooted-tree.qza",
  taxonomy="data/process/taxonomy.qza",
  metadata = "data/process/sal.metadata.tsv"
)
physeq #View summary of object

saveRDS(physeq, "data/phyloseq/phyloseq.RDS") #Save the object for use later


