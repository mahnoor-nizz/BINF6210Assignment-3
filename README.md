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

### Quality Control

Raw paired-end FASTQ reads were filtered using fastp v0.23.4 (Chen et al., 2018) with the following parameters. A minimum base quality score of 20 was used, corresponding to a 99% base call accuracy which is a widely used standard cutoff that removes unreliable bases while keeping a majority of high quality reads (Chen et al., 2018). A minimum read length of 50 bp was set after trimming, as shorter reads map with insufficient specificity during kmer-based classification and are likely to generate false positives in Kraken2 (Wood et al., 2019). Automatic adapter detection for paired-end reads (`--detect_adapter_for_pe`) was enabled rather than specifying adapter sequences manually, because fastp's overlap-based detection is more robust when adapter sequences are unknown or variable across library preparations (Chen et al., 2018). Quality reports were generated HTML format for each sample to verify filtering accuracy.

### Taxonomic Classification

Quality-filtered reads were classified using Kraken2 v2.1.6 (Wood et al., 2019) against the Standard 8 GB database (k2_standard_08gb_20241228). The 8 GB database was chosen over the full Standard database (which exceeds 60 GB) to reduce memory requirements while still getting coverage of the most common gut microbial taxa. A confidence threshold of 0.15 was applied (`--confidence 0.15`), which requires that at least 15% of the kmer evidence for a read must support its assigned taxon. This reduces false positives and is commonly recommended for gut microbiome studies (Wood et al., 2019). The `--memory-mapping` flag was omitted as Compute Canada nodes have sufficient RAM and loading the database directly into memory was significantly faster than memory-mapping.

Species-level abundance was re-estimated using Bracken v2.8 (Lu et al., 2017). Kraken2 classifies many reads at genus or higher taxonomic levels because closely related species share large fractions of their genomes, resulting in ambiguous kmer assignments. Bracken fixes this issue by redistributing higher level reads to the species level, producing more accurate species level abundance estimates (Lu et al., 2017). A minimum read threshold of 10 reads per taxon (`-t 10`) was applied, removing taxa supported by very few reads, which are likely to be noise or low-confidence assignments rather than real species. The read length parameter was set to 150 bp (`-r 150`) to match the sequencing data. All Bracken species reports were merged into a single BIOM-format file using kraken-biom v1.2.0 for import into R.

