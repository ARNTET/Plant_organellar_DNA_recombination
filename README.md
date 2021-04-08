# Pipeline of analysis for recombination of organellar genomes

## Introduction
This pipeline was initialy built for Chevigny et al., 2020 and aims to analyze recombination of organellar genomes from Illumina paired-end data (2 × 150 bp). To have an overview of this pipeline, please refer to Graphical_summary_PlantOrgRec.pdf 

## Methods
### Data pre-processing
Illumina paired-end reads (2 × 150 bp) were trimmed with trimmomatic and overlapping reads joined with fastq-join. Then reads were mapped on the three genomes: nuclear, plastid and mitochondrial and only reads mapped in pair were kept using BWA mem. 

### Organellar genome coverage analysis
This part is mandatory because the recombination part of the study require coverage file "mean_cov.csv".From the bam file containing reads mapped in proper pair, coverage data were extracted with Bedtools as well as the number of reads mapped using Samtools view. Then, using JGAF-Plant_organellar_coverage.R coverage data are normalized to 1'000'000 mapped per genomes, and plot of the means coverage per kb are drawned. 

### Organellar genome recombination analysis
Reads properly pairing to the cpDNA and mtDNA were filtered to only keep those showing a short-clipping sequence (threshold 20 nucleotides) and no indel (looking for the presence of an S in the CIGAR string without any I, D and H). The short-clipping sequences were then extracted with SE-MEI/extractSoftclipped (github.com/dpryan79/SE-MEI) and aligned using bowtie2 against the corresponding genome. The positions of the short-clipping sequences mapping the cpDNA and of their relative read were rounded down to the upper kb to analyze the localization of the rearrangement (see JGAF-Plant_organellar_recombination.R). Those corresponding to the isomerization that results from the recombination involving the large inverted repeats were filtered out.
