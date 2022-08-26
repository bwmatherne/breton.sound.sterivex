#Fast Parallel UniFrac (in R)

# Reference: https://joey711.github.io/phyloseq-demo/unifrac.html

# Two types of Unifrac calculations:

    #Weighted UniFrac - which does take into account differences in abundance of taxa between samples, but takes longer to calculate; and

    #Unweighted UniFrac - which only considers the presence/absence of taxa between sample pairs.



# Clear the Environment
rm(list = ls())

# Load Phyloseq

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("doParallel"); packageVersion("doParallel")
library("foreach"); packageVersion("foreach")

#Set random number generator (RNG) seed, for perfect reproducibility.

set.seed(629)


# Define a default theme for ggplot graphics.

theme_set(theme_bw())
fontsize = 18L
theme_update(axis.title.x = element_text(size=fontsize))
theme_update(axis.title.y = element_text(size=fontsize))
theme_update(plot.title = element_text(size=fontsize+2))

phy <- read_rds("data/phyloseq/phy.RDS")

wuni <- UniFrac(phy, weighted=TRUE, normalized=TRUE, parallel=TRUE, fast = TRUE)

time.frame.plot <- data.frame(time.frame, algorithm=time.frame[, "weighted"] )

ordinate(phy, "PCoA", "wUniFrac") %>% 
  plot_ordination(phy, ., color = "bin_sal", title = "weighted-UniFrac")


ordinate(phy, "NMDS", "wUniFrac") %>% 
  plot_ordination(phy, ., color = "bin_sal", title = "weighted-UniFrac")


plot_heatmap(physeq = phy, method = "MDS", distance = "wUniFrac", 
             title = "weighted-UniFrac", taxa.label = FALSE)


######

# Plot Ordination
library(plyr)

GP <- read_rds("data/phyloseq/phyloseq.RDS") # Sticking with "GP" to make it easier to follow tutorial
wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP)) # Remove OTUs that do not appear more than five times in more than half the samples
GP1 = prune_taxa(wh0, GP)
