#Preprocessing

PHYr  = transform_sample_counts(phy, function(x) x / sum(x) )
PHYfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)

# Actinobacteriota preprocessing
PHY.act = subset_taxa(phy, Phylum=="p__Actinobacteriota")
PHY.act = prune_samples(sample_sums(PHY.act)>=20, PHY.act)

PHY.act.merged = merge_taxa(PHY.act, taxa_names(PHY.act)[1:5])

transform_sample_counts(PHY.act, function(OTU) OTU/sum(OTU) )

#Proteobacteria preprocessing
PHY.gam = subset_taxa(phy, Class=="c__Gammaproteobacteria")
PHY.gam = prune_samples(sample_sums(PHY.gam)>=20, PHY.gam)
PHY.gam.merged = merge_taxa(PHY.gam, taxa_names(PHY.gam)[1:5])
gamma = tax_glom(PHY.gam, "Family")
plot_tree(gamma, color="Location", shape="Class", size="abundance")
transform_sample_counts(gamma, function(OTU) OTU/sum(OTU) )
transform_sample_counts(PHY.gam.merged, function(OTU) OTU/sum(OTU) )

g.tree <- plot_tree(gamma, color="Location", shape="Class", size="abundance")
ggsave("data/phyloseq/gamma_tree.png", g.tree,  width = 14, height = 10, dpi = 300)