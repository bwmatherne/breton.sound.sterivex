################################################################################
#
# phyloseq_driver.R
#
# R script for taking QIIME2 output and conducting in-depth statistics
#
# Dependencies:   metadata, taxonomy, and feature table files from qiime2
#									
#									
# Produces:       results/figures/?????
#
# Reference: https://uw-madison-microbiome-hub.github.io/Microbiome_analysis_in-_R/
#
################################################################################

# Clear the Environment
rm(list = ls())

# 4.1. Load Packages

library(DECIPHER) # This package will help in importing, maintaining, analyzing, manipulating, and exporting a massive amount of sequences.

library(ape) # Analyses of Phylogenetics and Evolution package. Required for tree calculations to be used with phyloseq

library(DESeq2) # This package will help analyze "differential expression" in the microbiota alongside phyloseq

library(ggplot2) # Graphing package used in phyloseq. To edit the default setting of a plot, you need to use functions in this package.

library(phyloseq) # The phyloseq package seeks to address issues with multiple microbiome analysis packages by providing a set of functions that internally manage the organizing, linking, storing, and analyzing of phylogenetic sequencing data. In general, this package is used for UniFrac analyses.

library(plotly) # A package to create interactive web graphics of use in 3D plots

library(vegan) # The vegan package provides tools for descriptive community ecology. It has most basic functions of diversity analysis, community ordination and dissimilarity analysis. In general, this package is used for Bray-Curtis and Jaccard analyses.

library(philr) # This package provides functions for the analysis of compositional data 

library(tidyverse) # This package is designed to make it easy to install and load multiple 'tidyverse' packages in a single step

library(adespatial) # Tools for the multiscale spatial analysis of multivariate data

library(devtools) # Make package development easier by providing R functions that simplify and expedite common tasks

library(qiime2R) # A package for importing qiime artifacts into an R session

library(MicrobeR) # Data visualization

library(microbiome) # Data analysis and visualization

library(microbiomeSeq) # Data analysis and visualization

library("pander") # provide a minimal and easy tool for rendering R objects into Pandoc's markdown

library(ranacapa) # Data analysis 

library(grid) # support data visualization

library(gridExtra)  # support data visualization

library(knitr) # Provides a general-purpose tool for dynamic report generation in R using Literate Programming techniques.

library(png) # Figure download

library("ggdendro") #set of tools for dendrograms and tree plots using 'ggplot2'

library(ggpubr) # publication quality figures, based on ggplot2

library(RColorBrewer) # nice color options

library(microbiomeutilities) # some utility tools 

#4.2 Load Data
#Convert qiime artifacts directly to phyloseq

# Importing ASVs abundance file
ASVs <- read_qza("data/process/table-with-phyla-no-mitochondria-no-chloroplast.qza")

#Importing metadata
metadata <- read.table("data/raw/combined.metadata.tsv", sep='\t', header=T, row.names=1, comment="")
metadata <- metadata[-1,] # remove the second line that specifies the data type

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
  sample_data(metadata)
)

# check for features of data  
summarize_phyloseq(phy)
print_ps(phy)
summary(sample_sums(phy))

# Accessing the phyloseq object

ntaxa(phy)

nsamples(phy)

sample_names(phy)[1:5]  

rank_names(phy)  

sample_variables(phy)  

otu_table(phy)[1:5, 1:5]  

tax_table(phy)[1:5, 1:4]

#Filter data by Site

physeq <- subset_samples(phy, Site == "BS")


# Distribution of reads

plot_read_distribution(physeq, groups = "Location", plot.type = "density") + theme_biome_utils()

#Rarefy the phyloseq object to even depth prior various analysis

physeq_rarefy <- rarefy_even_depth(physeq, rngseed=1, sample.size=0.9*min(sample_sums(physeq)), replace=F)

# Check the taxa prevalence at Phylum level

plot_taxa_prevalence(physeq, "Phylum")

# 5. Composition plots
# Barplots
physeq_fam <- microbiome::aggregate_rare(physeq, level = "Family", detection = 50/100, prevalence = 70/100)

physeq.fam.rel <- microbiome::transform(physeq_fam, "compositional")

physeq.fam.rel <- physeq %>%
  aggregate_rare(level = "Family", detection = 50/100, prevalence = 70/100) %>%
  microbiome::transform(transform = "compositional")

plot_composition(physeq.fam.rel,sample.sort = "Location", x.label = "SampleID") + theme(legend.position = "bottom") + scale_fill_brewer("Family", palette = "Paired") + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))

# Barplot Option 2

taxa_barplot(Summarize.Taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "Site")

  # To make it interactive
ggplotly(taxa_barplot(Summarize.Taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "Site"))

  # save the plot
b.plot <- taxa_barplot(Summarize.Taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "Site")

ggsave("results/figures/barplot_family.png", b.plot,  width = 14, height = 10, dpi = 300)

# Heatmap

taxa_heatmap(Summarize.Taxa(ASVs$data, as.data.frame(tax_table))$Family, metadata, "Site")

# Heatmap 2

library(pheatmap)

p <- plot_taxa_heatmap(physeq,
                       subset.top = 25,
                       VariableA = c("Location","Month"),
                       transformation = "log10",
                       cluster_rows = T,
                       cluster_cols = F,
                       show_colnames = F,
                       heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
)
#the plot is stored here
p$plot

# table used for plot is here
p$tax_tab[1:3,1:3]


# Heatmap 3

h.map <- plot_heatmap(physeq.fam.rel, method="PCoA", distance="bray", taxa.label = "Family", sample.order = unique(sample_names(physeq))) + facet_grid(~Location, scales = "free_x", drop = TRUE) + theme_bw() + theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1)) + theme(legend.key = element_blank(),strip.background = element_rect(colour="black", fill="white"))

# Make bacterial names italics
h.map <- h.map + theme(axis.text.y = element_text(colour = 'black', size = 10, face = 'italic'))

# Change the color palette
h.map <- h.map + scale_fill_distiller("Abundance", palette = "RdYlBu")

# clean the x-axis
h.map <- h.map + rremove("x.text")

print(h.map)

# Saving the plot
ggsave("results/figures/heatmap_family.png", h.map,  width = 14, height = 10, dpi = 300)

# Boxplot

physeq_df <- microbiomeutilities::phy_to_ldf(physeq_fam, 
                                             transform.counts = "compositional")

# An additonal column Sam_rep with sample names is created for reference purpose
colnames(physeq_df)

# Box plot at Family level

ggstripchart(physeq_df, "Location", "Abundance", 
             facet.by = "Family", color = "Location",
             palette = "jco") + rremove("x.text")

# Plot relative abundance of top taxa

mycols <- c("coral", "steelblue2", "slategray2", "olivedrab")

t.plot <- plot_taxa_boxplot(physeq,
                  taxonomic.level = "Family",
                  top.otu = 6, 
                  group = "Location",
                  add.violin= FALSE,
                  title = "Top Six Families", 
                  keep.other = FALSE,
                  group.order = c("Coastal","Inshore"),
                  group.colors = mycols,
                  dot.size = 2) + theme_biome_utils()
print(t.plot)
# Saving the plot
ggsave("results/figures/top6families_plot.png", t.plot,  width = 14, height = 10, dpi = 300)

# plotting taxa specified by the user

physeq.f <- format_to_besthit(physeq)

top_taxa(physeq.f, 5)

select.taxa <- c("6247ad94b599d7cadc5c5feb32b1ca05:p__Proteobacteria", "096e9e07cd5e6ad29c90f2de5b9731ba:d__Bacteria")

p <- plot_listed_taxa(physeq.f, select.taxa, 
                      group= "Location",
                      group.order = c("Coastal","Inshore"),
                      group.colors = mycols,
                      add.violin = F,
                      dot.opacity = 0.25,
                      box.opacity = 0.25,
                      panel.arrange= "grid") + ylab("Relative abundance") + scale_y_continuous(labels = scales::percent)
# Adding statistical test with ggpubr::stat_compare_means()

# If more than two variables
comps <- make_pairs(sample_data(physeq.f)$Location)
print(comps)

p <- p + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test") 

p + scale_y_continuous(labels = scales::percent)

# Plot top four genera

physeq.genus <- aggregate_taxa(physeq, "Genus")
top_four <- top_taxa(physeq.genus, 4)
top_four

top_genera <- plot_listed_taxa(physeq.genus, top_four, 
                               group= "Location",
                               group.order = c("Coastal","Inshore"),
                               group.colors = mycols,
                               add.violin = F,
                               dot.opacity = 0.25,
                               box.opacity = 0.25,
                               panel.arrange= "wrap")
top_genera

tg.plot <- top_genera + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test")

print(tg.plot)
# Saving the plot
ggsave("results/figures/top4genera_plot.png", tg.plot,  width = 14, height = 10, dpi = 300)


#Dominant Taxa

physeq.gen <- aggregate_taxa(physeq,"Genus")

dom.tax <- dominant_taxa(physeq,level = "Genus", group="Location")
head(dom.tax$dominant_overview)

# Taxa summary - entire dataset

taxa_summary(physeq, "Phylum")

#Taxa summary by groups
    # For group specific abundances of taxa

grp_abund <- get_group_abundances(physeq, 
                                  level = "Phylum", 
                                  group="Location",
                                  transform = "compositional")

# clean names 
grp_abund$OTUID <- gsub("p__", "",grp_abund$OTUID)
grp_abund$OTUID <- ifelse(grp_abund$OTUID == "", 
                          "Unclassified", grp_abund$OTUID)

mean.plot <- grp_abund %>% # input data
  ggplot(aes(x= reorder(OTUID, mean_abundance), # rerorder based on mean abundance
             y= mean_abundance,
             fill=Location)) + # x and y axis 
  geom_bar(     stat = "identity", 
                position=position_dodge()) + 
  scale_fill_manual("Location", values=mycols) + # manually specify colors
  theme_bw() + # add a widely used ggplot2 theme
  ylab("Mean Relative Abundance") + # label y axis
  xlab("Phylum") + # label x axis
  coord_flip() # rotate plot 

mean.plot

# Saving the plot
ggsave("results/figures/mean_plot.png", mean.plot,  width = 14, height = 10, dpi = 300)

#Find samples dominated by specific taxa
bact_dom <- find_samples_taxa(physeq.gen, taxa = "g__Acinetobacter")

#Find  samples dominated by "g__Acinetobacter"
ps.sub <- prune_samples(sample_names(physeq.gen) %in% bact_dom, physeq.gen)
print(ps.sub) # Needs more steps to be useful


#Alpha Diversities
    # Plot rarefraction curve
ggrare(physeq, step = 50, color="Location", label = "Sample", se = TRUE)

#phyloseq also allows you to easily plot alpha diversity, both by sample and by group

plot_richness(physeq_rarefy, measures="Shannon")

a.div <- plot_richness(physeq_rarefy, x="Location", measures=c("Shannon", "simpson", "Observed"), color = "Location") + geom_boxplot() + theme_bw()

# adding statistical support
a.div + stat_compare_means(
  comparisons = comps,
  label = "p.signif",
  tip.length = 0.05,
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
    symbols = c("xxxx", "***", "**", "*", "ns")
  ),
  method = "wilcox.test")

# generate a csv file of the richness estimates using

richness <- estimate_richness(physeq_rarefy,measures=c("Shannon", "simpson", "Observed"))
write.csv(richness, file = "alpha_div.csv")

#Creating a plot with one index and stats

plot_diversity_stats(physeq, group = "Location", 
                     index = "diversity_shannon", 
                     group.order = c("Inshore","Coastal"),                      
                     group.colors = mycols,
                     label.format="p.format",
                     stats = TRUE) + ylab("Shannon Diversity") + xlab("")
# 9. Beta diversity metrices
# 9.1 Non - Phylogenetic beta diversity metrics
physeq.ord <- ordinate(physeq_rarefy, "PCoA", "bray")
b.div.bray <- plot_ordination(physeq_rarefy, physeq.ord, type= "samples", color= "Location") + geom_point(size=3)
b.div.bray <- b.div.bray + stat_ellipse() + ggtitle("Bray Curtis")  + theme_classic() + scale_color_brewer("Location", palette = "Set2")
print(b.div.bray)

#9.2 Phylogenetic beta diversity metrics
#Weighted Unifrac will consider the abundances of different taxa.
# convert ot relative abundance
physeq_rel <- microbiome::transform(physeq, "compositional")
physeq.ord.wuni <- ordinate(physeq_rel, "PCoA", "unifrac", weighted=T)
b.div.wuni <- plot_ordination(physeq_rel, physeq.ord.wuni, type= "samples", color= "Location") + geom_point(size=3)
b.div.wuni <- b.div.wuni + stat_ellipse() + ggtitle("Weighted Unifrac")  + theme_classic() + scale_color_brewer("Location", palette = "Set2")
print(b.div.wuni)


#9.3 Permanova
#Permutational multivariate analysis of variance further reading

otu <- abundances(physeq_rel)
meta <- meta(physeq_rel)


#Statistics - Bdiv
permanova <- adonis(t(otu) ~ Location, data = meta, permutations=99, method = "bray")

#P-value
print(as.data.frame(permanova$aov.tab)["Location", "Pr(>F)"])

#9.4 Checking the homogeneity condition
#More infromation can be found by typing ?betadisper

#Pair - wise stats
dist <- vegdist(t(otu), "bray")
anova(betadisper(dist, meta$Location))

permutest(betadisper(dist, meta$Location), pairwise = TRUE)

#10. Core microbiota
#Subset the data to keep only Inshore samples.

physeq.inshore <- subset_samples(physeq, Location == "Inshore")

# convert to relative abundance  
physeq.inshore.rel <- microbiome::transform(physeq.inshore, "compositional")

physeq.inshore.rel2 <- prune_taxa(taxa_sums(physeq.inshore.rel) > 0, physeq.inshore.rel)

#Check for the core ASVs

core.taxa.standard <- core_members(physeq.inshore.rel2, detection = 0.001, prevalence = 20/100)
print(core.taxa.standard)
#we only see IDs, not very informative. We can get the classification of these as below.

# Extract the taxonomy table
taxonomy_core <- as.data.frame(tax_table(physeq.inshore.rel2))

# Subset this taxonomy table to include only core OTUs
core_taxa_id <- subset(taxonomy_core, rownames(taxonomy_core) %in% core.taxa.standard)

DT::datatable(core_taxa_id)

###
#Subset the data to keep only Coastal samples.

physeq.coastal <- subset_samples(physeq, Location == "Coastal")

# convert to relative abundance  
physeq.coastal.rel <- microbiome::transform(physeq.coastal, "compositional")

physeq.coastal.rel2 <- prune_taxa(taxa_sums(physeq.coastal.rel) > 0, physeq.coastal.rel)

#Check for the core ASVs

core.taxa.standard.coast <- core_members(physeq.coastal.rel2, detection = 0.001, prevalence = 20/100)
print(core.taxa.standard.coast)
#we only see IDs, not very informative. We can get the classification of these as below.

# Extract the taxonomy table
taxonomy_core_coast <- as.data.frame(tax_table(physeq.coastal.rel2))

# Subset this taxonomy table to include only core OTUs
core_taxa_id_coast <- subset(taxonomy_core_coast, rownames(taxonomy_core_coast) %in% core.taxa.standard.coast)

DT::datatable(core_taxa_id_coast)

#10.1 Core abundance and diversity
#Total core abundance in each sample (sum of abundances of the core members):
  
  core.abundance <- sample_sums(core(physeq.inshore.rel2, detection = 0.001, prevalence = 20/100))

DT::datatable(as.data.frame(core.abundance))

################################# Skipped to section 12
#12. Microbiome network
#You can plot the distances between ASVs as a network.

plot_net(physeq_rel, maxdist = 0.8, color = "Location")
#change distance to Jaccard
plot_net(physeq_rel, maxdist = 0.8, color = "Location", distance="jaccard")

#12.1 igraph-based network
ig <- make_network(physeq_rel, max.dist=0.8)
plot_network(ig, physeq_rel)

# Add color label 
plot_network(ig, physeq_rel, color="Location", line_weight=0.4, label=NULL)

#replace the Jaccard (default) distance method with Bray-Curtis
ig <- make_network(physeq_rel, dist.fun="bray", max.dist=0.8)
plot_network(ig, physeq_rel, color="Location", line_weight=0.4, label=NULL)
      #####Note: For co-occurrence networks of OTUs, I suggest trying Gephi or Cytoscape

#13. Differential abundance testing
# DeSeq2 to test for differential abundance between categories

#Convert phyloseq object ot DeSeq
bsdds <- phyloseq_to_deseq2(physeq_rarefy, ~ Location)
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(bsdds), 1, gm_mean)
bsdds <- estimateSizeFactors(bsdds, geoMeans = geoMeans)
# DeSeq function tests for differential abundance 
bsdds <- DESeq(bsdds, test="Wald", fitType="parametric")

# Results function call creates a table of the results of the tests
res <- results(bsdds, cooksCutoff = FALSE)
alpha <- 0.01
sigtab <- res[which(res$padj < alpha), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(physeq_rarefy)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Cleaning up the table a little for legibility
posigtab <- sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab <- posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

# Bar plot showing the log2-fold-change, showing Genus and Phylum. Uses some ggplot2 commands
sigtabgen <- subset(sigtab, !is.na(Genus))
# Phylum order
x <- tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x <- sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x <- tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
