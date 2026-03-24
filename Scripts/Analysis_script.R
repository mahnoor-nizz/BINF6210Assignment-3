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
metadata <- data.frame(
  diet = c("omnivore", "omnivore", "omnivore", "vegan", "vegan", "vegan"),
  row.names = sample_names(physeq)
)
sample_data(physeq) <- sample_data(metadata)

# Verify
print(physeq)
print(as.data.frame(sample_data(physeq)))
