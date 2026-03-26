# BINF6210Assignment-3

## Table of Contents
1. [Introduction](#introduction)
2. [Methods](#methods)
3. [Results](#results)
4. [Discussion](#discussion)
5. [References](#references)


---

## Introduction

The human gut microbiome contains a complex community of microorganisms that play an essential role in the host’s metabolism and immune function. [^8] The composition of the gut microbiome is affected by a range of factors including genetics, age, geographic origin, antibiotic use, and diet, making it highly variable between individuals. [^9]

Diet is known to be one of the strongest factors effecting gut microbial composition. Plant based diets are rich in dietary fibre and complex carbohydrates, serving as substrates for microbial fermentation and driving the production of short-chain fatty acids (SCFAs). [^?21] While omnivorous diets are characterized by higher intakes of animal protein and saturated fat, which have been associated with enrichment of bile tolerant species and reduced SCFA-producing bacteria. [^19] Understanding differences associated with diet has implications for human health, as specific microbes and their metabolites have been linked to inflammation, metabolic syndrome, and colorectal cancer risk.

Shotgun metagenomics is an approach used for characterizing microbial communities without relying on a single marker gene. Unlike 16S rRNA amplicon sequencing (metabarcoding), which amplifies only a specific gene region and is limited by primer bias and database completeness, shotgun metagenomics sequences all DNA present in a sample. This allows for taxonomic profiling at a greater resolution, direct detection of functional pathways, and the possibility of metagenome-assembled genome (MAG) reconstruction. [^18] Shotgun approaches do require greater sequencing depth and more sophisticated bioinformatic pipelines but offer substantially higher sensitivity for rare taxa and better accuracy with species level identification [^5].

Here, shotgun metagenomics was used to compare the gut microbiome composition of vegan and omnivore individuals using publicly available data from Filippis et al, (2019) accessed via NCBI (accession SRP126540) Taxonomic classification was performed using Kraken2 and Bracken, alpha and beta diversity were assessed, and differential abundance between diet groups was tested for using ANCOM-BC2.

---

## Methods

### Data Source

Raw shotgun metagenomic sequencing reads were downloaded from the NCBI Sequence Read Archive (SRA) under BioProject accession SRP126540, corresponding to the study by Filippis et al, (2019). Six samples were selected for analysis: three from vegan participants (SRR8146944, SRR8146963, SRR8146968) and three from omnivore participants (SRR8146935, SRR8146936, SRR8146938). Reads were downloaded using SRA Toolkit v3.0.9 (`prefetch` and `fasterq-dump`) from the NCBI.

### Quality Control

Raw paired-end FASTQ reads were filtered using fastp v0.23.4 [^3] with the following parameters. A minimum base quality score of 20 was used, corresponding to a 99% base call accuracy which is a widely used standard cutoff that removes unreliable bases while keeping a majority of high quality reads [^3]. A minimum read length of 50 bp was set after trimming, as shorter reads map with insufficient specificity during kmer-based classification and are likely to generate false positives in Kraken2 [^20]. Automatic adapter detection for paired-end reads (`--detect_adapter_for_pe`) was enabled rather than specifying adapter sequences manually, because fastp's overlap-based detection is more robust when adapter sequences are unknown or variable across library preparations.[^3] Quality reports were generated HTML format for each sample to verify filtering accuracy.

### Taxonomic Classification

Quality-filtered reads were classified using Kraken2 v2.1.6 [^20] against the Standard 8 GB database (k2_standard_08gb_20241228). The 8 GB database was chosen over the full Standard database (which exceeds 60 GB) to reduce memory requirements while still getting coverage of the most common gut microbial taxa. A confidence threshold of 0.15 was applied (`--confidence 0.15`), which requires that at least 15% of the kmer evidence for a read must support its assigned taxon. This reduces false positives and is commonly recommended for gut microbiome studies. [^20] The `--memory-mapping` flag was omitted as Compute Canada nodes have sufficient RAM and loading the database directly into memory was significantly faster than memory-mapping.

Species-level abundance was re-estimated using Bracken v2.8. [^14] Kraken2 classifies many reads at genus or higher taxonomic levels because closely related species share large fractions of their genomes, resulting in ambiguous kmer assignments. Bracken fixes this issue by redistributing higher level reads to the species level, producing more accurate species level abundance estimates. [^14] A minimum read threshold of 10 reads per taxon (`-t 10`) was applied, removing taxa supported by very few reads, which are likely to be noise or low-confidence assignments rather than real species. The read length parameter was set to 150 bp (`-r 150`) to match the sequencing data. All Bracken species reports were merged into a single BIOM format file using krakenbiom v1.2.0 for import into R.

### Diversity Analysis and Visualization

All of the following analyses were performed in R v4.5.2. The BIOM file was imported as a phyloseq object[^15] using the `biomformat` v1.38.0 and `phyloseq` v1.54.2 packages. Sample metadata (diet group) was added to the phyloseq object.

**Rarefaction** curves were generated using the `vegan` v2.7.2 package[^17] to check whether sequencing depth was sufficient enough to capture the community diversity. A plateau in the rarefaction curve indicated that increasing sequencing depth was unlikely to reveal more species and observed richness indeed reflected true community composition.

**Taxonomic abundance** was visualized at the phylum, class, order, genus, and species levels using stacked bar charts after transformation to relative abundance. Relative abundance was used rather than raw counts because samples differed in read depth, and raw counts would show differences in sequencing depth rather than the real biological differences in abundance. Taxa were aggregated using `tax_glom()` and the top 10 (20 for species) taxa by mean abundance were kept, with the rest grouped as "Other" to produce figures.

**Alpha diversity** was estimated using `estimate_richness()` from phyloseq, including observed species richness, Shannon diversity index, and Simpson dominance index. Chao1 richness estimation was not included because Bracken's minimum read threshold removed low-abundance taxa and eliminated singletons from the data. Shannon entropy and Simpson's index were used instead, as both captured richness and evenness andante require singletons. Shannon is more sensitive to rare taxa while Simpson is more sensitive to dominant taxa. [^2] Differences in alpha diversity between diet groups were assessed using the Wilcoxon rank-sum test, as it is a non-parametric test that is useable with small sample sizes where a normal distribution is not assumed.

**Beta diversity** was assessed using Bray-Curtis and Jaccard. Bray-Curtis was chosen as the primary metric because it accounts for differences in taxon abundance between communities and is the most widely used beta diversity metric in gut microbiome studies [^13]. Jaccard distance was also included to try and see whether observed differences between groups reflected changes in which taxa are present versus changes in their relative proportions. So, if both metrics produce similar group separation, the communities differ in composition, not just abundance. [^13] Principal Coordinates Analysis (PCoA) was performed on both distance matrices using `ordinate()`. Bray-Curtis NMDS ordination was also computed. NMDS optimizes rank-order preservation of distances and does not assume linearity, providing a robustness check on the PCoA results. Statistical significance of group separation was tested using PERMANOVA (`adonis2()`)which tests whether the centroids and dispersion of groups differ more than expected by chance using permutation of the distance matrix [^1].

**Differential abundance** analysis was done using ANCOM-BC2 [^11] via the `ANCOMBC` Bioconductor package (v2.12.0). ANCOM-BC2 was selected over DESeq2 or Wilcoxon tests because metagenomics count data is compositional and observed abundances are relative rather than absolute, which violates the independence assumptions of standard statistical tests. [^7]ANCOM-BC2 models and corrects for compositional bias and has been shown to have lower false positive rates than other methods in benchmarking studies. [^16] The Holm multiple testing correction was applied rather than Benjamini-Hochberg (BH/FDR) because Holm provides strong family-wise error rate control, which is preferred in a small-sample.[^11] Structural zero detection was enabled (`struc_zero = TRUE`) to distinguish taxa that are truly absent in a group. A minimum library size of 1000 reads (`lib_cut = 1000`) was set to exclude any samples with very low coverage. Taxa present in fewer than 10% of samples were excluded (`prv_cut = 0.10`) to remove rare taxa with insufficient data for reliable inference, a standard filtering step recommended by ANCOM-BC2 developers.

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

Taxonomic classification with Kraken2 and Bracken identified between 68 and 123 species per sample (Table 1). Total classified reads per sample ranged from 348,713 to 642,499. The gut microbiomes of all samples were dominated by the phyla Bacillota and Bacteroidota, consistent with the known composition of the healthy human gut microbiome. [^9]

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

**Figure 4.** Species-level relative abundance showing the top 20 species by mean abundance. *Segatella copri*, *Alistipes putredinis*, and *Bacteroides uniformis* were among the most abundant species detected. *Segatella copri* showed particularly high relative abundance in vegan samples, consistent with its well-documented enrichment in plant-based gut microbiomes. [^28]

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

## Discussion

This study compared the gut microbiome composition of vegan and omnivore individuals using shotgun metagenomics, taxonomic classification with Kraken2 and Bracken, and downstream diversity and differential abundance analyses. While no statistically significant differences were found between diet groups in alpha diversity, beta diversity, or differential abundance, the results revealed trends that were biologically meaningful and consistent with the existing literature on diet microbiome interactions.

Bacillota (Firmicutes) and Bacteroidota (Bacteroidetes) had the highest rates in all samples at the phylum level, which is consistent with the known composition of the healthy human adult gut microbiome. [^9] The noticeable inter-individual variability observed across samples is also consistent with large-scale microbiome studies, which have repeatedly shown that personal microbiome composition varies more between individuals than between dietary groups, particularly in small cohorts.[^9] This variability likely contributed to the lack of statistically significant findings in this study and highlights the need for larger sample sizes in diet-microbiome research.

Omnivore samples showed higher median observed richness (111 species) and Shannon diversity (3.10) compared to vegan samples (87 species; Shannon 3.08), though these differences were not statistically significant. This was surprising as some studies show that plant-based diets may be associated with higher microbial diversity due to increased dietary fibres, which support a broader range of fermentative species. [^19] This difference varies across studies however, with some studies reporting higher diversity in omnivores, possibly a reflection of the broader range of substrates available from a mixed diet. [^7] Others have found that there was no significant difference between diet groups at all, with inter individual variation dominating over dietary effects. [^9] Since the sample size was so small (n = 3 per group), these results were not sufficient enough to be able to draw conclusions about the true relationship between diet and alpha diversity in this population.

The Bray-Curtis PCoA showed partial separation between diet groups along PC1, which explained 60.5% of variance, suggesting that community composition may differ between vegans and omnivores. The Jaccard PCoA showed a similar pattern, with PC1 explaining 47.7% of variance. The fact that both abundance-weighted (Bray-Curtis) and presence/absence (Jaccard) metrics showed similar patterns suggests that the dietary groups may differ in both which taxa are present and their relative proportions, rather than just one or the other. However, PERMANOVA confirmed that neither result was statistically significant, consistent with the limited statistical power at this sample size. 

For differential abundance trends, although no genera reached statistical significance after Holm correction, the log fold change estimates from ANCOM-BC2 reveal some trends that are well supported by literature. Segatella (LFC = +2.58) showed the largest positive fold change in vegans. Segatella copri is one of the most consistently reported diet associated bacteria, with many studies documenting its enrichment in individuals on plant rich, high fibre diets. [^19] It is able to ferment complex plant polysaccharides and produces propionate and succinate as fermentation end products, adding beneficial short chain fatty acids (SCFAs) to the gut. [^28]

Ruminococcus (LFC = +2.52) also showed a large positive fold change in vegans. Members of this genus are some of the most important degraders of resistant starch and plant cell wall polysaccharides in the human gut.[^6] Their enrichment in vegan samples is consistent with the higher dietary fibre intake associated with plant-based diets.

In contrast, Mediterraneibacter (LFC = -2.25) and Lachnospira (LFC = -2.14) were higher in omnivore samples. Both belong to the family Lachnospiraceae, which contains many butyrate-producing species. Lachnospira has been associated with omnivorous dietary patterns and has been shown to respond positively to animal protein intake.[^28] The higher abundance of these genera in omnivores in our dataset may reflect the different fermentable substrates available in a mixed diet, including animal-derived oligosaccharides and mucins that some Lachnospiraceae members can degrade.

Limitations and future directions. The primary limitation of this study is the small sample size of n = 3 per group, which limits statistical power. A minimum of 10–20 samples per group is recommended for gut microbiome studies in order to obtain reliable differential abundance results.[^10] [^16] The use of the 8 GB Standard Kraken2 database rather than a more expansive/comprehensive reference may also have resulted in some reads remaining unclassified, particularly for less common taxa. [^20] Future work should expand to include all available samples from SRP126540 and possibly include functional profiling to assess whether observed taxonomic differences translate into differences in metabolic pathway abundance.

In conclusion, while no statistically significant differences between vegan and omnivore gut microbiomes were detected in this study, the trends observed are biologically consistent with the known effects of plant-based diets on gut microbiome composition. The enrichment of Segatella and Ruminococcus in vegan samples and Mediterraneibacter and Lachnospira in omnivore samples aligns well with the existing literature and supports the hypothesis that diet is a meaningful driver of gut microbiome composition. Larger, well-powered studies are needed to confirm these observations.SRP


---

## References

[^1]: Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecology, 26(1), 32–46.
[^2]: Cassol, M., et al. (2025). A comprehensive guide to alpha diversity metrics in microbiome research. Scientific Reports, 14, 27864.
[^3]: Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884–i890.
[^4]: Dahl, W. J., et al. (2019). Diet, nutrients, and the microbiome. Cell Host & Microbe, 25(3), 332–346.
[^5]: Durazzi, F., et al. (2021). Comparison between 16S rRNA and shotgun sequencing data for the taxonomic characterization of the gut microbiota. Scientific Reports, 11, 3585.
[^6]: Flint, H. J., Scott, K. P., Duncan, S. H., Louis, P., & Forano, E. (2012). Microbial degradation of complex carbohydrates in the gut. Gut Microbes, 3(4), 289–306.
[^7]: Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J. (2017). Microbiome datasets are compositional: and this is not optional. Frontiers in Microbiology, 8, 2224.
[^8]: Grice, E. A., & Segre, J. A. (2012). The human microbiome: our second genome. Annual Review of Genomics and Human Genetics, 13, 151–170.
[^9]: Human Microbiome Project Consortium. (2012). Structure, function and diversity of the healthy human microbiome. Nature, 486, 207–214.
[^10]: Kelly, B. J., et al. (2015). Power and sample-size estimation for microbiome studies using pairwise distances and PERMANOVA. Bioinformatics, 31(15), 2461–2468.
[^11]: Lin, H., & Peddada, S. D. (2023). Analysis of compositions of microbiomes with bias correction 2 (ANCOM-BC2). Nature Communications, 14, 1918.
[^12]: Lopez-Siles, M., et al. (2017). Faecalibacterium prausnitzii: from microbiology to diagnostics and prognostics. ISME Journal, 11, 841–852.
[^13]: Lozupone, C. A., et al. (2011). Diversity, stability and resilience of the human gut microbiota. Nature, 489, 220–230.
[^14]: Lu, J., et al. (2017). Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science, 3, e104.
[^15]: McMurdie, P. J., & Holmes, S. (2013). phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE, 8(4), e61217.
[^16]: Nearing, J. T., et al. (2022). Microbiome differential abundance methods produce different results across 38 datasets. Nature Communications, 13, 342.
[^17]: Oksanen, J., et al. (2022). vegan: Community Ecology Package. R package version 2.6-4.
[^18]: Quince, C., et al. (2017). Shotgun metagenomics, from sampling to analysis. Nature Biotechnology, 35, 833–844.
[^19]: Sonnenburg, J. L., & Bäckhed, F. (2016). Diet–microbiota interactions as moderators of human metabolism. Nature, 535, 56–64.
[^20]: Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20, 257.
[^21]: Dahl, W. J., et al. (2021). Diet and the intestinal virome are interrelated and associated with host genetics. Cell Host & Microbe.
[^22]: Ghozlane, A., et al. (2025). Meteor2: improved metagenomic profiling of the human gut microbiome. Bioinformatics.
[^23]: Marić, J., et al. (2024). False positive reduction in Kraken2 with long reads. Genome Biology.
[^24]: Plovier, H., et al. (2017). A purified membrane protein from Akkermansia muciniphila or the pasteurized bacterium improves metabolism in obese and diabetic mice. Nature Medicine, 23, 107–113.
[^25]: Sokol, H., et al. (2009). Faecalibacterium prausnitzii is an anti-inflammatory commensal bacterium identified by gut microbiota analysis of Crohn disease patients. PNAS, 105(43), 16731–16736.
[^26]: Takeuchi, T., et al. (2021). Lactic acid bacteria and Phocaeicola vulgatus in gut microbiome shifts. Cell Host & Microbe.
[^27]: Waters, J. L., & Ley, R. E. (2019). The human gut bacteria Christensenellaceae are widespread, heritable, and associated with health. BMC Biology, 17, 83.
[^28}: De Filippis, Francesca et al. “Distinct Genetic and Functional Traits of Human Intestinal Prevotella copri Strains Are Associated with Different Habitual Diets.” Cell host & microbe vol. 25,3 (2019): 444-453.e3. doi:10.1016/j.chom.2019.01.004
