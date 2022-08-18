# qiime2R

# Reference https://rdrr.io/github/jbisanz/qiime2R/f/README.md 

#Load dependencies

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")


library(qiime2R)

# Importing ASVs abundance file
ASVs <- read_qza("data/process/bs.table.qza")

#Importing metadata
metadata <- read_q2metadata("data/process/bs.metadata.tsv")
metadata$Month  <- as.factor(metadata$Month)
metadata$Month <- fct_relevel(metadata$Month, c("1","2","3","4","5","6","8","9","10")) # Set month order for down stream plotting

# Importing tree
tree <- read_qza("data/process/rooted-tree.qza")
# Importing taxonomy
taxonomy <- read_qza("data/process/taxonomy.qza")
taxonomy <- parse_taxonomy(taxonomy$data)

# Create phyloseq object using wrapper
physeq<-qza_to_phyloseq(
  features="data/process/bs.table.qza",
  tree="data/process/rooted-tree.qza",
  taxonomy="data/process/taxonomy.qza",
  metadata = "data/process/bs.metadata.tsv"
)
physeq

# Creating phyloseq object based on Madison example
phy <- phyloseq(
  otu_table(ASVs$data, taxa_are_rows = T),
  phy_tree(tree$data),
  tax_table(taxonomy),
  sample_data(metadata))