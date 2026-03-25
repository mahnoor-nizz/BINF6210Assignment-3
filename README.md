# BINF6210Assignment-3

## Introduction

The human gut microbiome is a complex community of microorganisms, consisting primarily of bacteria, that plays an essential role in the host metabolism, immune function, and protection against pathogens (Grice & Segre, 2012). The composition of the gut microbiome is highly variable between individuals and is shaped by a range of factors including genetics, age, geographic origin, antibiotic use, and notably, diet (Human Microbiome Project Consortium, 2012).

Diet is known to be one of the most powerful modulators of gut microbial composition. Plant based diets are rich in dietary fibre and complex carbohydrates which serve as substrates for microbial fermentation, driving the production of short-chain fatty acids (SCFAs) (Dahl et al., 2023). While omnivorous diets are characterized by higher intakes of animal protein and saturated fat, which have been associated with enrichment of bile tolerant species and reduced SCFA-producing bacteria (Sonnenburg & Bäckhed, 2016). Understanding differences associated with diet has implications for human health, as specific microbes and their metabolites have been linked to inflammation, metabolic syndrome, and colorectal cancer risk.

Shotgun metagenomics is an approach used for characterizing microbial communities without relying on a single marker gene. Unlike 16S rRNA amplicon sequencing (metabarcoding), which amplifies only a specific gene region and is limited by primer bias and database completeness, shotgun metagenomics sequences all DNA present in a sample. This allows for taxonomic profiling at a greater resolution, direct detection of functional pathways, and the possibility of metagenome-assembled genome (MAG) reconstruction (Quince et al., 2017). Shotgun approaches do require greater sequencing depth and more sophisticated bioinformatic pipelines but offer substantially higher sensitivity for rare taxa and better accuracy with species level identification (Durazzi et al., 2021).

Here, shotgun metagenomics was used to compare the gut microbiome composition of vegan and omnivore individuals using publicly available data from Dahl et al. (2019), accessed via NCBI SRA (accession SRP126540). Taxonomic classification was performed using Kraken2 and Bracken, alpha and beta diversity were assessed, and differential abundance between diet groups was tested for using ANCOM-BC2.

---

## Methods

### Data Source

Raw shotgun metagenomic sequencing reads were downloaded from the NCBI Sequence Read Archive (SRA) under BioProject accession SRP126540, corresponding to the study by Dahl et al. (2019). Six samples were selected for analysis: three from vegan participants (SRR8146944, SRR8146963, SRR8146968) and three from omnivore participants (SRR8146935, SRR8146936, SRR8146938). Reads were downloaded using SRA Toolkit v3.0.9 (`prefetch` and `fasterq-dump`) from the NCBI.
