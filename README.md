# BINF6210Assignment-3

## Introduction

The human gut microbiome is a complex community of microorganisms, consisting primarily of bacteria, that plays an essential role in the host metabolism, immune function, and protection against pathogens (Grice & Segre, 2012). The composition of the gut microbiome is highly variable between individuals and is shaped by a range of factors including genetics, age, geographic origin, antibiotic use, and notably, diet (Human Microbiome Project Consortium, 2012).

Diet is known to be one of the most powerful modulators of gut microbial composition. Plant based diets are rich in dietary fibre and complex carbohydrates which serve as substrates for microbial fermentation, driving the production of short-chain fatty acids (SCFAs) (Dahl et al., 2023). While omnivorous diets are characterized by higher intakes of animal protein and saturated fat, which have been associated with enrichment of bile tolerant species and reduced SCFA-producing bacteria (Sonnenburg & Bäckhed, 2016). Understanding differences associated with diet has implications for human health, as specific microbes and their metabolites have been linked to inflammation, metabolic syndrome, and colorectal cancer risk.

Shotgun metagenomics is an approach used for characterizing microbial communities without relying on a single marker gene. Unlike 16S rRNA amplicon sequencing (metabarcoding), which amplifies only a specific gene region and is limited by primer bias and database completeness, shotgun metagenomics sequences all DNA present in a sample. This allows for taxonomic profiling at a greater resolution, direct detection of functional pathways, and the possibility of metagenome-assembled genome (MAG) reconstruction (Quince et al., 2017). Shotgun approaches do require greater sequencing depth and more sophisticated bioinformatic pipelines but offer substantially higher sensitivity for rare taxa and better accuracy with species level identification (Durazzi et al., 2021).

Here, shotgun metagenomics was used to compare the gut microbiome composition of vegan and omnivore individuals using publicly available data from Fillipis et al. (2019), accessed via NCBI SRA (accession SRP126540). Taxonomic classification was performed using Kraken2 and Bracken, alpha and beta diversity were assessed, and differential abundance between diet groups was tested for using ANCOM-BC2.

---

## Methods

### Data Source

Raw shotgun metagenomic sequencing reads were downloaded from the NCBI Sequence Read Archive (SRA) under BioProject accession SRP126540, corresponding to the study by Fillipis et al. (2019). Six samples were selected for analysis: three from vegan participants (SRR8146944, SRR8146963, SRR8146968) and three from omnivore participants (SRR8146935, SRR8146936, SRR8146938). Reads were downloaded using SRA Toolkit v3.0.9 (`prefetch` and `fasterq-dump`) from the NCBI.

### Quality Control

Raw paired-end FASTQ reads were filtered using fastp v0.23.4 (Chen et al., 2018) with the following parameters. A minimum base quality score of 20 was used, corresponding to a 99% base call accuracy which is a widely used standard cutoff that removes unreliable bases while keeping a majority of high quality reads (Chen et al., 2018). A minimum read length of 50 bp was set after trimming, as shorter reads map with insufficient specificity during kmer-based classification and are likely to generate false positives in Kraken2 (Wood et al., 2019). Automatic adapter detection for paired-end reads (`--detect_adapter_for_pe`) was enabled rather than specifying adapter sequences manually, because fastp's overlap-based detection is more robust when adapter sequences are unknown or variable across library preparations (Chen et al., 2018). Quality reports were generated HTML format for each sample to verify filtering accuracy.

### Taxonomic Classification

Quality-filtered reads were classified using Kraken2 v2.1.6 (Wood et al., 2019) against the Standard 8 GB database (k2_standard_08gb_20241228). The 8 GB database was chosen over the full Standard database (which exceeds 60 GB) to reduce memory requirements while still getting coverage of the most common gut microbial taxa. A confidence threshold of 0.15 was applied (`--confidence 0.15`), which requires that at least 15% of the kmer evidence for a read must support its assigned taxon. This reduces false positives and is commonly recommended for gut microbiome studies (Wood et al., 2019). The `--memory-mapping` flag was omitted as Compute Canada nodes have sufficient RAM and loading the database directly into memory was significantly faster than memory-mapping.

Species-level abundance was re-estimated using Bracken v2.8 (Lu et al., 2017). Kraken2 classifies many reads at genus or higher taxonomic levels because closely related species share large fractions of their genomes, resulting in ambiguous kmer assignments. Bracken fixes this issue by redistributing higher level reads to the species level, producing more accurate species level abundance estimates (Lu et al., 2017). A minimum read threshold of 10 reads per taxon (`-t 10`) was applied, removing taxa supported by very few reads, which are likely to be noise or low-confidence assignments rather than real species. The read length parameter was set to 150 bp (`-r 150`) to match the sequencing data. All Bracken species reports were merged into a single BIOM-format file using kraken-biom v1.2.0 for import into R.

### Diversity Analysis and Visualization

All downstream analyses were performed in R v4.5.2. The BIOM file was imported into a phyloseq object (McMurdie & Holmes, 2013) using the `biomformat` v1.38.0 and `phyloseq`v1.54.2 packages. Sample metadata including diet group assignment was added to the phyloseq object.

**Rarefaction** curves were generated using the `vegan` v2.7.2 package (Oksanen et al., 2022) to assess whether sequencing depth was sufficient to capture community diversity. A plateau in the rarefaction curve indicates that increasing sequencing depth is unlikely to reveal substantially more species, providing confidence that observed richness reflects true community composition rather than sampling limitation.

**Taxonomic abundance** was visualized at the phylum, class, order, genus, and species levels using stacked bar charts after transformation to relative abundance (`transform_sample_counts`). Relative abundance transformation was used rather than raw counts because samples differ in total classified read depth, and raw counts would mix up true biological differences with sequencing depth differences. Taxa were aggregated using `tax_glom()` and the top 10 (20 for species) taxa by mean abundance were retained, with the remainder grouped as "Other" to produce legible figures.

**Alpha diversity** was estimated using `estimate_richness()` from phyloseq, including observed species richness, Shannon diversity index, and Simpson dominance index. Chao1 richness estimation was excluded because Bracken's minimum read threshold of 10 reads removes low-abundance taxa, eliminating singletons from the data. Chao1 relies on singleton counts to estimate unobserved species richness, without singletons it produces unreliable estimates (Cassol et al., 2025). Shannon entropy and Simpson's index were used instead, as both capture richness and evenness and neither requires singleton counts. Shannon is more sensitive to rare taxa while Simpson is more sensitive to dominant taxa. (Cassol et al., 2025). Differences in alpha diversity between diet groups were assessed using the Wilcoxon rank-sum test, a non-parametric test appropriate for small sample sizes where normality cannot be assumed.

**Beta diversity** was assessed using two distance metrics: Bray-Curtis and Jaccard. Bray-Curtis was chosen as the primary metric because it accounts for differences in taxon abundance between communities and is the most widely used beta diversity metric in gut microbiome studies (Lozupone et al., 2011). Jaccard distance was included to see whether observed differences between groups reflect changes in which taxa are present versus changes in their relative proportions, if both metrics produce similar group separation, the communities differ in composition, not just abundance (Lozupone et al., 2011). Principal Coordinates Analysis (PCoA) was performed on both distance matrices using `ordinate()`. Bray-Curtis NMDS ordination was also computed. NMDS optimizes rank-order preservation of distances and does not assume linearity, providing a robustness check on the PCoA results. Statistical significance of group separation was tested using PERMANOVA (`adonis2()`)which tests whether the centroids and dispersion of groups differ more than expected by chance using permutation of the distance matrix (Anderson, 2001).

**Differential abundance** analysis was done using ANCOM-BC2 (Lin & Peddada, 2023) via the `ANCOMBC` Bioconductor package (v2.12.0). ANCOM-BC2 was selected over DESeq2 or Wilcoxon tests because metagenomics count data is compositional and observed abundances are relative rather than absolute, which violates the independence assumptions of standard statistical tests (Gloor et al., 2017). ANCOM-BC2 models and corrects for compositional bias and has been shown to have lower false positive rates than other methods in benchmarking studies (Nearing et al., 2022). The Holm multiple testing correction was applied rather than Benjamini-Hochberg (BH/FDR) because Holm provides strong family-wise error rate control, which is preferred in a small-sample (Lin & Peddada, 2023). Structural zero detection was enabled (`struc_zero = TRUE`) to distinguish taxa that are truly absent in a group. A minimum library size of 1000 reads (`lib_cut = 1000`) was set to exclude any samples with very low coverage. Taxa present in fewer than 10% of samples were excluded (`prv_cut = 0.10`) to remove rare taxa with insufficient data for reliable inference, a standard filtering step recommended by ANCOM-BC2 developers.

All analyses were conducted on Compute Canada's Narval cluster.

---

## Results

### Sequencing Depth and Quality

After quality filtering with fastp, all six samples retained high read counts suitable for metagenomics analysis. Rarefaction curves are shown in Figure 1. All samples approached plateau, indicating that sequencing depth was sufficient to capture the majority of detectable species in these samples given the database used.

---

![Rarefaction Curves](https://github.com/mahnoor-nizz/BINF6210Assignment-3/blob/main/Figures/fig0_rarefaction.png)

**Figure 1.** Rarefaction curves for all six samples. Red = omnivore samples (SRR8146935, SRR8146936, SRR8146938); blue = vegan samples (SRR8146944, SRR8146963, SRR8146968). Curves approaching plateau indicate adequate sequencing depth.

---

### Taxonomic Composition

Taxonomic classification with Kraken2 and Bracken identified between 68 and 123 species per sample (Table 1). Total classified reads per sample ranged from 348,713 to 642,499. The gut microbiomes of all samples were dominated by the phyla Bacillota and Bacteroidota, consistent with the known composition of the healthy human gut microbiome (Human Microbiome Project Consortium, 2012).

---

**Table 1.** Summary of classified reads and species richness per sample.

| Sample | Diet | Classified Reads | Species Detected |
|--------|------|-----------------|-----------------|
| SRR8146935 | Omnivore | 399,593 | 111 |
| SRR8146936 | Omnivore | 463,048 | 68 |
| SRR8146938 | Omnivore | 348,713 | 123 |
| SRR8146944 | Vegan | 417,665 | 81 |
| SRR8146963 | Vegan | 393,352 | 87 |
| SRR8146968 | Vegan | 642,499 | 93 |

---


Phylum-level relative abundance is shown in Figure 2. Both diet groups had Bacillota and Bacteroidota in highest number, with the ratio between the two varying between individuals. Class, order, genus, and species-level composition plots were created, showing the high degree of inter-individual variability characteristic of gut microbiome studies.

---

![Phylum-Level Relative Abundance](https://github.com/mahnoor-nizz/BINF6210Assignment-3/blob/main/Figures/fig1_abundance_phylum.png)

**Figure 2.** Phylum-level relative abundance for omnivore (left) and vegan (right) samples. Each bar represents one sample. The top 10 phyla by mean abundance are shown; remaining phyla are grouped as "Other". Both groups are dominated by Bacillota (Firmicutes) and Bacteroidota (Bacteroidetes).

---

![Genus-Level Relative Abundance](https://github.com/mahnoor-nizz/BINF6210Assignment-3/blob/main/Figures/fig4_abundance_genus.png)

**Figure 3.** Genus-level relative abundance for omnivore (left) and vegan (right) samples. Top 10 genera by mean abundance are shown. *Segatella*, *Bacteroides*, and *Alistipes* are among the most prominent genera, with *Segatella* notably dominant in several vegan samples, consistent with its known association with plant-based diets.

---

![Species-Level Relative Abundance](https://github.com/mahnoor-nizz/BINF6210Assignment-3/blob/main/Figures/fig5_abundance_species.png)

**Figure 4.** Species-level relative abundance showing the top 20 species by mean abundance. *Segatella copri*, *Alistipes putredinis*, and *Bacteroides uniformis* were among the most abundant species detected. *Segatella copri* showed particularly high relative abundance in vegan samples, consistent with its well-documented enrichment in plant-based gut microbiomes (Dahl et al., 2019).

---

### Alpha Diversity

Alpha diversity estimates were summarized in Table 2 and visualized in Figure 5. Omnivore samples had a median observed richness of 111 species, compared to 87 in vegan samples. Shannon diversity was higher in omnivore samples compared to vegan samples. Simpson dominance was also higher in omnivore samples. None of the alpha diversity differences between groups reached statistical significance by Wilcoxon rank-sum test (Observed: p = 0.70; Shannon: p = 0.40; Simpson: p = 0.70), likely reflecting the small sample size of n = 3 per group.

**Table 2.** Alpha diversity summary by diet group.

| Diet | Median Observed | Range Observed | Median Shannon | Range Shannon | Median Simpson |
|------|----------------|---------------|---------------|--------------|---------------|
| Omnivore | 111 | 68–123 | 3.10 | 1.89–3.42 | 0.914 | 
| Vegan | 87 | 81–93 | 3.08 | 1.70–2.30 | 0.930 |

---

![Alpha Diversity](https://github.com/mahnoor-nizz/BINF6210Assignment-3/blob/main/Figures/fig6_alpha_diversity.png)

**Figure 5.** Alpha diversity by diet group. Observed species richness, Shannon diversity, and Simpson dominance are shown for omnivore (red) and vegan (blue) samples. Each point represents one sample. No significant differences were detected between groups (Wilcoxon rank-sum test, all p > 0.05). Chao1 was excluded as Bracken's minimum read threshold eliminates singletons required for it.

---

### Beta Diversity

Beta diversity was assessed using Bray Curtis dissimilarity and Jaccard distance. PCoA plots are shown in Figures 6 and 7. For Bray Curtis the PCoA suggests partial separation between diet groups along PC1, which explained 60.5% of variance, with PC2 explaining a further 16.3%. However, because of the small sample size there was overlap between groups. The Jaccard PCoA was similar with PC1 explaining 47.7% and PC2 explaining 17.3% of variance. PERMANOVA confirmed that diet group did not explain a statistically significant proportion of community variance (Table 3). 

---

![Bray-Curtis PCoA](https://github.com/mahnoor-nizz/BINF6210Assignment-3/blob/main/Figures/fig7a_beta_braycurtis_pcoa.png)

**Figure 6.** Principal Coordinates Analysis (PCoA) of Bray-Curtis dissimilarity. Each point represents one sample, coloured by diet group (red = omnivore, blue = vegan). PC1 and PC2 explain 60.5%  and 16.3% of total variance, respectively. Partial separation between groups is visible along PC1, though PERMANOVA did not detect a statistically significant difference (R² = 0.143, p = 0.70).

---

![Bray-Curtis NMDS](https://github.com/mahnoor-nizz/BINF6210Assignment-3/blob/main/Figures/fig7b_beta_braycurtis_nmds.png)

**Figure 7.** Non-metric Multidimensional Scaling (NMDS) of Bray-Curtis dissimilarity. The near-zero stress value (< 0.001) indicates the distance structure is well-represented in two dimensions, though this also reflects the very small dataset size (n = 6). Grouping pattern is consistent with the PCoA.

---

**Table 3.** PERMANOVA results for beta diversity.

| Distance Metric | R² | F-statistic | p-value | Permutations |
|----------------|-----|-------------|---------|-------------|
| Bray-Curtis | 0.143 | 0.666 | 0.70 | 720 |
| Jaccard | 0.161 | 0.765 | 0.70 | 720 |

---

### Differential Abundance

Differential abundance analysis was performed using ANCOM-BC2 at the genus level. After filtering taxa present in fewer than 10% of samples, 38 genera were retained for analysis. No genera reached statistical significance after Holm multiple testing correction (q < 0.05; Figure 8). ANCOM-BC2 explicitly warned that variance estimation may be unstable given the sample size of n = 3 per group, which is below the recommended minimum of 5 per group.
Despite the lack of statistically significant results, the log fold change estimates were able to reveal meaningful biological trends. Genera with the largest positive fold changes in vegans included Segatella (LFC = +2.58) and Ruminococcus (LFC = +2.52), both consistent with plant-based dietary patterns and fibre fermentation. Genera with the largest negative fold changes in vegans (i.e., higher in omnivores) included Mediterraneibacter (LFC = -2.25) and Lachnospira (LFC = -2.14). These trends are biologically plausible but should be interpreted cautiously given the low statistical power.
Structural zero analysis identified 17 genera completely absent in one diet group but present in the other (Table 4).

---

![Differential Abundance Dot Plot](https://github.com/mahnoor-nizz/BINF6210Assignment-3/blob/main/Figures/fig8_differential_abundancev2.png)

**Figure 8.** ANCOM-BC2 differential abundance results at the genus level. Each point represents one genus; the x-axis shows the log fold change (vegan vs. omnivore). Error bars indicate the standard error of the log fold change estimate. The red dashed line indicates no difference (LFC = 0). Points are coloured by direction: blue = higher in vegan, red = higher in omnivore. No genera reached statistical significance after Holm correction (n = 3 per group).

---

**Table 4.** Selected taxa with structural zeroes in one diet group only (completely absent in one group, present in the other).

| Taxon | Absent in Omnivore | Absent in Vegan |
|-------|-------------------|----------------|
| *Rothia* | | **X** |
| *Cutibacterium* | | **X** |
| *Eggerthella* | **X** | |
| *Gordonibacter* | **X** | |
| *Christensenella* | | **X** |
| *Clostridium* | | **X** |
| *Anaerotruncus* | **X** | |
| *Pseudoflavonifractor* | | **X** |
| *Hungatella* | **X** | |
| *Wansuia* | **X** | |
| *Wujia* | | **X** |
| *Holdemania* | | **X** |
| *Phascolarctobacterium* | | **X** |
| *Veillonella* | | **X** |
| *Shutterella* | **X** | |
| *Haemophilus* | | **X** |
| *Endlipuvirus* | | **X** |

---


