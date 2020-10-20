# Study recombination of the Arabidopsis thaliana radA mutant on cpDNA

Materials

Raw reads fastq files are available on ebi.ac.uk/ena under the accession numbers:

WT1: ERR4686393

WT2: ERR4686394

radA-1: ERR4686477

radA-2: ERR4686478

radA-3: ERR4686479

radA-4: ERR4686480

Mapping bam files containing only reads with a short-clipped sequences are available on ebi.ac.uk/ena under the accession numbers:

WT1: ERR4695569

WT2: ERR4695570

radA-1: ERR4695869

radA-2: ERR4695870

radA-3: ERR4695871

radA-4: ERR4695872

Short-clipped fastq sequences are available on ebi.ac.uk/ena under the accession numbers:

WT1: ERR4695865

WT2: ERR4695866

radA-1: ERR4695861

radA-2: ERR4695862

radA-3: ERR4695863

radA-4: ERR4695864

Methods

Illumina sequencing

In this study Illumina paired-end short reads sequencing strategy was used to determine the impact of the Arabidopsis thaliana radA mutant on the cpDNA. To do so, total leaf DNA of WT and radA plants was quantified with a QuBit Fluorometer (Life Technologies) and libraries were prepared with the Nextera Flex Library kit, according to manufacturer’s recommendations (Illumina) using 100 ng of each DNA sample. Final libraries were quantified, checked on a Bioanalyzer 2100 (Agilent) and sequenced on an Illumina Miseq system (2 × 150 paired-end reads).

Coverage analysis

For Illumina sequence analysis of the cpDNA, reads were aligned against the Arabidopsis reference genomes using BWA and filtered to keep only those mapping to the cpDNA. For analysis of rearranged sequences, reads properly pairing to the cpDNA were filtered to only keep those showing a short-clipping sequence (threshold 20 nucleotides) and no indel (looking for the presence of an S in the CIGAR string without any I, D and H). The short-clipping sequences were then extracted with SE-MEI/extractSoftclipped (github.com/dpryan79/SE-MEI) and aligned using bowtie2 against the cpDNA (see JGAF-Atha_cpDNA_SoftClipped.sh). The positions of the short-clipping sequences mapping the cpDNA and of their relative read were rounded down to the upper kb to analyze the localization of the rearrangement (see JGAF-Atha_cpDNA_SoftClipped.R). Those corresponding to the isomerization that results from the recombination involving the large inverted repeats were filtered out.
