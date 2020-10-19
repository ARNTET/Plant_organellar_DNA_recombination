#!/bin/bash
#SBATCH --partition=fast
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=recombination
#SBATCH --output=JGAF-recombination-%j.out

# Modules needed ###############################################################
module load bwa
module load samtools
module load bowtie2

# WT-1 #####################################################################
bwa mem -t 16 \
Atha_cpDNA.fasta \
JGAF-WT-1_R1.fq \
JGAF-WT-1_R2.fq \
> JGAF-WT-1_cpDNA_mapped.sam
samtools view -b -F4 \
JGAF-WT-1_cpDNA_mapped.sam \
> JGAF-WT-1_cpDNA_mapped.unsorted.bam
samtools view \
-h JGAF-WT-1_cpDNA_mapped.unsorted.bam \
| awk '{if($0 ~ /^@/ || $6 ~ /S/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /I/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /D/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /H/) {print $0}}' \
| samtools view -Sb - > JGAF-WT-1_cpDNA_mapped.SC.unsorted.bam
samtools sort \
JGAF-WT-1_cpDNA_mapped.SC.unsorted.bam \
> JGAF-WT-1_cpDNA_mapped.SC.bam
SE-MEI/extractSoftclipped \
JGAF-WT-1_cpDNA_mapped.SC.bam \
> JGAF-WT-1.SC.fq.gz
gunzip JGAF-WT-1.SC.fq.gz
awk \
'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 20) {print header, seq, qheader, qseq}}' \
< ARA_RADA_SOFT/SOFT_CLIPPED/JGAF-WT-1.SC.fq > ARA_RADA_SOFT/SOFT_CLIPPED/JGAF-WT-1.SC.flt.fq
bowtie2 \
--no-mixed \
-x Atha_cpDNA.bt2 \
-U JGAF-WT-1.SC.flt.fq \
-S JGAF-WT-1.SC.sam
samtools view -b \
JGAF-WT-1.SC.sam \
> JGAF-WT-1.SC.unsorted.bam
samtools sort  \
JGAF-WT-1.SC.unsorted.bam \
> JGAF-WT-1.SC.bam
samtools view -b \
JGAF-WT-1.SC.bam \
> JGAF-WT-1.SC.sam

# WT-2 #####################################################################
bwa mem -t 16 \
Atha_cpDNA.fasta \
JGAF-WT-2_R1.fq \
JGAF-WT-2_R2.fq \
> JGAF-WT-2_cpDNA_mapped.sam
samtools view -b -F4 \
JGAF-WT-2_cpDNA_mapped.sam \
> JGAF-WT-2_cpDNA_mapped.unsorted.bam
samtools view \
-h JGAF-WT-2_cpDNA_mapped.unsorted.bam \
| awk '{if($0 ~ /^@/ || $6 ~ /S/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /I/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /D/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /H/) {print $0}}' \
| samtools view -Sb - > JGAF-WT-2_cpDNA_mapped.SC.unsorted.bam
samtools sort \
JGAF-WT-2_cpDNA_mapped.SC.unsorted.bam \
> JGAF-WT-2_cpDNA_mapped.SC.bam
SE-MEI/extractSoftclipped \
JGAF-WT-2_cpDNA_mapped.SC.bam \
> JGAF-WT-2.SC.fq.gz
gunzip JGAF-WT-2.SC.fq.gz
awk \
'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 20) {print header, seq, qheader, qseq}}' \
< ARA_RADA_SOFT/SOFT_CLIPPED/JGAF-WT-2.SC.fq > ARA_RADA_SOFT/SOFT_CLIPPED/JGAF-WT-2.SC.flt.fq
bowtie2 \
--no-mixed \
-x Atha_cpDNA.bt2 \
-U JGAF-WT-2.SC.flt.fq \
-S JGAF-WT-2.SC.sam
samtools view -b \
JGAF-WT-2.SC.sam \
> JGAF-WT-2.SC.unsorted.bam
samtools sort  \
JGAF-WT-2.SC.unsorted.bam \
> JGAF-WT-2.SC.bam
samtools view -b \
JGAF-WT-2.SC.bam \
> JGAF-WT-2.SC.sam

# radA-1-1 #####################################################################
bwa mem -t 16 \
Atha_cpDNA.fasta \
JGAF-radA-1-1_R1.fq \
JGAF-radA-1-1_R2.fq \
> JGAF-radA-1-1_cpDNA_mapped.sam
samtools view -b -F4 \
JGAF-radA-1-1_cpDNA_mapped.sam \
> JGAF-radA-1-1_cpDNA_mapped.unsorted.bam
samtools view \
-h JGAF-radA-1-1_cpDNA_mapped.unsorted.bam \
| awk '{if($0 ~ /^@/ || $6 ~ /S/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /I/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /D/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /H/) {print $0}}' \
| samtools view -Sb - > JGAF-radA-1-1_cpDNA_mapped.SC.unsorted.bam
samtools sort \
JGAF-radA-1-1_cpDNA_mapped.SC.unsorted.bam \
> JGAF-radA-1-1_cpDNA_mapped.SC.bam
SE-MEI/extractSoftclipped \
JGAF-radA-1-1_cpDNA_mapped.SC.bam \
> JGAF-radA-1-1.SC.fq.gz
gunzip JGAF-radA-1-1.SC.fq.gz
awk \
'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 20) {print header, seq, qheader, qseq}}' \
< ARA_RADA_SOFT/SOFT_CLIPPED/JGAF-radA-1-1.SC.fq > ARA_RADA_SOFT/SOFT_CLIPPED/JGAF-radA-1-1.SC.flt.fq
bowtie2 \
--no-mixed \
-x Atha_cpDNA.bt2 \
-U JGAF-radA-1-1.SC.flt.fq \
-S JGAF-radA-1-1.SC.sam
samtools view -b \
JGAF-radA-1-1.SC.sam \
> JGAF-radA-1-1.SC.unsorted.bam
samtools sort  \
JGAF-radA-1-1.SC.unsorted.bam \
> JGAF-radA-1-1.SC.bam
samtools view -b \
JGAF-radA-1-1.SC.bam \
> JGAF-radA-1-1.SC.sam

# radA-1-2 #####################################################################
bwa mem -t 16 \
Atha_cpDNA.fasta \
JGAF-radA-1-2_R1.fq \
JGAF-radA-1-2_R2.fq \
> JGAF-radA-1-2_cpDNA_mapped.sam
samtools view -b -F4 \
JGAF-radA-1-2_cpDNA_mapped.sam \
> JGAF-radA-1-2_cpDNA_mapped.unsorted.bam
samtools view \
-h JGAF-radA-1-2_cpDNA_mapped.unsorted.bam \
| awk '{if($0 ~ /^@/ || $6 ~ /S/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /I/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /D/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /H/) {print $0}}' \
| samtools view -Sb - > JGAF-radA-1-2_cpDNA_mapped.SC.unsorted.bam
samtools sort \
JGAF-radA-1-2_cpDNA_mapped.SC.unsorted.bam \
> JGAF-radA-1-2_cpDNA_mapped.SC.bam
SE-MEI/extractSoftclipped \
JGAF-radA-1-2_cpDNA_mapped.SC.bam \
> JGAF-radA-1-2.SC.fq.gz
gunzip JGAF-radA-1-2.SC.fq.gz
awk \
'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 20) {print header, seq, qheader, qseq}}' \
< ARA_RADA_SOFT/SOFT_CLIPPED/JGAF-radA-1-2.SC.fq > ARA_RADA_SOFT/SOFT_CLIPPED/JGAF-radA-1-2.SC.flt.fq
bowtie2 \
--no-mixed \
-x Atha_cpDNA.bt2 \
-U JGAF-radA-1-2.SC.flt.fq \
-S JGAF-radA-1-2.SC.sam
samtools view -b \
JGAF-radA-1-2.SC.sam \
> JGAF-radA-1-2.SC.unsorted.bam
samtools sort  \
JGAF-radA-1-2.SC.unsorted.bam \
> JGAF-radA-1-2.SC.bam
samtools view -b \
JGAF-radA-1-2.SC.bam \
> JGAF-radA-1-2.SC.sam

# radA-1-3 #####################################################################
bwa mem -t 16 \
Atha_cpDNA.fasta \
JGAF-radA-1-3_R1.fq \
JGAF-radA-1-3_R2.fq \
> JGAF-radA-1-3_cpDNA_mapped.sam
samtools view -b -F4 \
JGAF-radA-1-3_cpDNA_mapped.sam \
> JGAF-radA-1-3_cpDNA_mapped.unsorted.bam
samtools view \
-h JGAF-radA-1-3_cpDNA_mapped.unsorted.bam \
| awk '{if($0 ~ /^@/ || $6 ~ /S/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /I/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /D/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /H/) {print $0}}' \
| samtools view -Sb - > JGAF-radA-1-3_cpDNA_mapped.SC.unsorted.bam
samtools sort \
JGAF-radA-1-3_cpDNA_mapped.SC.unsorted.bam \
> JGAF-radA-1-3_cpDNA_mapped.SC.bam
SE-MEI/extractSoftclipped \
JGAF-radA-1-3_cpDNA_mapped.SC.bam \
> JGAF-radA-1-3.SC.fq.gz
gunzip JGAF-radA-1-3.SC.fq.gz
awk \
'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 20) {print header, seq, qheader, qseq}}' \
< ARA_RADA_SOFT/SOFT_CLIPPED/JGAF-radA-1-3.SC.fq > ARA_RADA_SOFT/SOFT_CLIPPED/JGAF-radA-1-3.SC.flt.fq
bowtie2 \
--no-mixed \
-x Atha_cpDNA.bt2 \
-U JGAF-radA-1-3.SC.flt.fq \
-S JGAF-radA-1-3.SC.sam
samtools view -b \
JGAF-radA-1-3.SC.sam \
> JGAF-radA-1-3.SC.unsorted.bam
samtools sort  \
JGAF-radA-1-3.SC.unsorted.bam \
> JGAF-radA-1-3.SC.bam
samtools view -b \
JGAF-radA-1-3.SC.bam \
> JGAF-radA-1-3.SC.sam

# radA-1-4 #####################################################################
bwa mem -t 16 \
Atha_cpDNA.fasta \
JGAF-radA-1-4_R1.fq \
JGAF-radA-1-4_R2.fq \
> JGAF-radA-1-4_cpDNA_mapped.sam
samtools view -b -F4 \
JGAF-radA-1-4_cpDNA_mapped.sam \
> JGAF-radA-1-4_cpDNA_mapped.unsorted.bam
samtools view \
-h JGAF-radA-1-4_cpDNA_mapped.unsorted.bam \
| awk '{if($0 ~ /^@/ || $6 ~ /S/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /I/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /D/) {print $0}}' \
| awk '{if($0 ~ /^@/ || $6 !~ /H/) {print $0}}' \
| samtools view -Sb - > JGAF-radA-1-4_cpDNA_mapped.SC.unsorted.bam
samtools sort \
JGAF-radA-1-4_cpDNA_mapped.SC.unsorted.bam \
> JGAF-radA-1-4_cpDNA_mapped.SC.bam
SE-MEI/extractSoftclipped \
JGAF-radA-1-4_cpDNA_mapped.SC.bam \
> JGAF-radA-1-4.SC.fq.gz
gunzip JGAF-radA-1-4.SC.fq.gz
awk \
'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 20) {print header, seq, qheader, qseq}}' \
< ARA_RADA_SOFT/SOFT_CLIPPED/JGAF-radA-1-4.SC.fq > ARA_RADA_SOFT/SOFT_CLIPPED/JGAF-radA-1-4.SC.flt.fq
bowtie2 \
--no-mixed \
-x Atha_cpDNA.bt2 \
-U JGAF-radA-1-4.SC.flt.fq \
-S JGAF-radA-1-4.SC.sam
samtools view -b \
JGAF-radA-1-4.SC.sam \
> JGAF-radA-1-4.SC.unsorted.bam
samtools sort  \
JGAF-radA-1-4.SC.unsorted.bam \
> JGAF-radA-1-4.SC.bam
samtools view -b \
JGAF-radA-1-4.SC.bam \
> JGAF-radA-1-4.SC.sam

# Extract Soft clips positions and sequences####################################
for f in *.sam; do
    mv -- "$f" "${f%.sam}.txt"
done
