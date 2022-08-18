
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")


library(qiime2R)
library(tidyverse)
if (!requireNamespace("gplots", quietly = TRUE)){install.packages("gplots")} #Needed for vinn diagram

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

################################
# Alpha diversity plot

# Load shannon data
shannon <- read.csv(file='data/phyloseq/alpha_div.csv')

# Add SampleID column name and remove Simpson column
shannon$SampleID <- shannon$X
shannon <- shannon %>%
  select(2:5) %>%
  relocate(SampleID, .before = Observed) %>%
  select(-Simpson)

head(shannon)

gplots::venn(list(metadata=metadata$SampleID, shannon=shannon$SampleID))

#Join shannon data with metadata file
metadata<-
  metadata %>% 
  left_join(shannon)
head(metadata)

#plot

p<- metadata %>%
  filter(!is.na(Shannon)) %>%
  ggplot(aes(x=Location, y=Shannon, fill=`Location`)) +
  stat_summary(geom="bar", fun.data=mean_se, color="black") + #here black is the outline for the bars
  geom_jitter(shape=21, width=0.2, height=0) +
  coord_cartesian(ylim=c(2,7)) + # adjust y-axis
  facet_grid(~`Month`) + # create a panel for each body site
  xlab("Sampling Location") +
  ylab("Shannon Diversity") +
  theme_q2r() +
  scale_fill_manual(values=c("cornflowerblue","indianred")) + #specify custom colors
  theme(legend.position="none") #remove the legend as it isn't needed

comps <- make_pairs(sample_data(metadata)$Location)
print(comps)

p <- p + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test") 
print(p)

ggsave("data/phyloseq/Shannon_by_Location.png", p,  width = 14, height = 10, dpi = 300)

