# Metagenomics analysis
# diversity + differential abundance
# Assignment 3 - Shotgun Metagenomics

library(phyloseq)
library(biomformat)
library(vegan)
library(ggplot2)
library(ANCOMBC)
library(dplyr)
library(tidyr)
library(tibble)

#####################################
# Load Data
#####################################

biom_data <- read_biom("table.biom")
physeq <- import_biom(biom_data)

# taxonomy column names
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# remove prefixes
tax_table(physeq) <- gsub("^[a-z]__", "", tax_table(physeq))

# Clean sample names
sample_names(physeq) <- gsub("_bracken$", "", sample_names(physeq))

# Add metadata 
metadata <- data.frame(diet = c("omnivore", "omnivore", "omnivore", "vegan", "vegan", "vegan"), row.names = sample_names(physeq))
sample_data(physeq) <- sample_data(metadata)

# Verify
print(physeq)
print(as.data.frame(sample_data(physeq)))

#####################################
# Rarefraction curve 
#####################################
otu_mat <- as.data.frame(t(otu_table(physeq)))
rarecurve(otu_mat, step=500, col=c(rep("maroon", 3), rep("navy", 3)), label=TRUE, main="Rarefaction Curves")

#####################################
# Taxonomic Abundance
#####################################

# Convert to relative abundance
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# function to aggregation + plot
plot_abundance <- function(physeq_rel, taxrank, top_n=10) {
  physeq_agg <- tax_glom(physeq_rel, taxrank=taxrank)
  df <- psmelt(physeq_agg)
  
  top_taxa <- df %>%
    group_by(.data[[taxrank]]) %>%
    summarise(Mean=mean(Abundance), .groups="drop") %>%
    slice_max(Mean, n=top_n) %>%
    pull(.data[[taxrank]])
  
  df <- df %>%
    mutate(!!taxrank := ifelse(.data[[taxrank]] %in% top_taxa, .data[[taxrank]], "Other"))
  
  ggplot(df, aes(x=Sample, y=Abundance, fill=.data[[taxrank]])) +
    geom_bar(stat="identity", position="stack") +
    facet_wrap(~diet, scales="free_x") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="right") +
    labs(title=paste0(taxrank, "-Level Relative Abundance"), y="Relative Abundance", x="Sample", fill=taxrank, caption=paste0("Top ", top_n, " shown; remaining grouped as Other"))
}

# Generate abundance plots at each taxonomic level
plot_abundance(physeq_rel, "Phylum",  top_n=10)

plot_abundance(physeq_rel, "Class",   top_n=10)

plot_abundance(physeq_rel, "Order",   top_n=10)

plot_abundance(physeq_rel, "Genus",   top_n=10)

plot_abundance(physeq_rel, "Species", top_n=20)



#####################################
# Alpha diversity
#####################################

alpha_div <- estimate_richness(physeq, measures=c("Observed", "Shannon", "Simpson"))

# add diet
alpha_div$diet   <- sample_data(physeq)$diet
alpha_div$sample <- rownames(alpha_div)

# Verify
print(alpha_div)

# Summary table
alpha_div %>%
  group_by(diet) %>%
  summarise(median_observed = median(Observed), range_observed  = paste(min(Observed), max(Observed), sep="-"), median_shannon  = median(Shannon), range_shannon   = paste(round(min(Shannon),2), round(max(Shannon),2), sep="-"), median_simpson  = median(Simpson), .groups="drop") %>%
  print()

# Plot
alpha_long <- alpha_div %>%
  pivot_longer(cols=c("Observed", "Shannon", "Simpson"), names_to="Metric", values_to="Value")

p_alpha <- ggplot(alpha_long, aes(x=diet, y=Value, fill=diet, color=diet)) +
  geom_boxplot(alpha=0.5, outlier.shape=NA) +
  geom_jitter(width=0.15, size=2) +
  facet_wrap(~Metric, scales="free_y") +
  scale_fill_manual(values=c("vegan"="steelblue", "omnivore"="firebrick")) +
  scale_color_manual(values=c("vegan"="steelblue", "omnivore"="firebrick")) +
  theme_bw() +
  labs(title="Alpha Diversity by Diet Group", y="Diversity Value", x="Diet", caption="Chao1 omitted (no singletons in Bracken output). Each point = one sample.")

ggsave("fig6_alpha_diversity.png", p_alpha, width=10, height=8, dpi=300)
cat("Saved: fig6_alpha_diversity.png\n")

# Wilcoxon tests
wilcox.test(Observed ~ diet, data=alpha_div)
wilcox.test(Shannon  ~ diet, data=alpha_div)
wilcox.test(Simpson  ~ diet, data=alpha_div)


#####################################
# Alpha diversity
#####################################

# Bray-Curtis PCoA
ord.pcoa.bray <- ordinate(physeq_rel, method="PCoA", distance="bray")
                                      
# Extract coordinates for ggplot
pcoa_df <- as.data.frame(ord.pcoa.bray$vectors[, 1:2])
pcoa_df$diet <- sample_data(physeq_rel)$diet

# % variance for axis labels
eig <- ord.pcoa.bray$values$Eigenvalues
pct <- round(eig / sum(eig) * 100, 1)

ggplot(pcoa_df, aes(x=Axis.1, y=Axis.2, color=diet)) +
  geom_point(size=4) +
  scale_color_manual(values=c("omnivore"="maroon", "vegan"="navy")) +
  theme_classic() +
  labs(title="Beta Diversity: Bray-Curtis PCoA", x=paste0("PC1 [", pct[1], "%]"), y=paste0("PC2 [", pct[2], "%]"), color="Diet", caption="Each point = one sample")



#####################################
# PERMANOVA
#####################################
metadata <- as(sample_data(physeq_rel), "data.frame")
                                      
# bray
permanova_bray <- adonis2(phyloseq::distance(physeq_rel, method="bray") ~ diet, data=metadata, permutations=999)
print(permanova_bray)

# Jaccard
permanova_jaccard <- adonis2(phyloseq::distance(physeq_rel, method="jaccard") ~ diet, data=metadata, permutations=999)
print(permanova_jaccard)
