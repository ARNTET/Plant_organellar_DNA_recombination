#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000

# Modules needed ###############################################################
module load bwa
module load samtools
module load bedtools
module load extractSoftclipped
module load bowtie2

# /!\ To fill /!\ ##############################################################
samples=()
genomes=() #ex: mtDNA cpDNA
lengthS=${#samples[@]}
lengthG=${#genomes[@]}
adapters=

# Genomes indexing #############################################################
G=0
while (($G<$lengthG)); do
  bwa index ${genomes[G]}".fasta"
  bowtie2-build ${genomes[G]}".fasta" ${genomes[G]}".bt2"
  let "G=G+1"
done

# Step one #####################################################################
S=0
mkdir coverage
mkdir recombination
mkdir a_trimmed
mkdir b_joined
mkdir c_mappingB
mkdir d_SC
mkdir e_SCfiltering
mkdir g_mappingBt
while (($S<$lengthS)); do
  G=0
  #A# Trimming
  trimmomatic \
  PE \
  ${samples[a]}"_R1.fastq.gz" \
  ${samples[a]}"_R2.fastq.gz" \
  a_trimmed/${samples[a]}".R1.fq.gz" \
  a_trimmed/${samples[a]}".R1U.fq.gz" \
  a_trimmed/${samples[a]}".R2.fq.gz" \
  a_trimmed/${samples[a]}".R2U.fq.gz" \
  ILLUMINACLIP:$adapters:2:30:10 SLIDINGWINDOW:4:15 MINLEN:75
  #B# Merging overlaping reads
  fastq-join \
  ${samples[S]}"_R1.fastq" \
  ${samples[S]}"_R2.fastq" \
  -o b_joined/${samples[S]}"_R1.fj.fastq" \
  -o b_joined/${samples[S]}"_R2.fj.fastq" \
  -o b_joined/${samples[S]}"_joined.fj.fastq"
  gzip b_joined/${samples[S]}"_R1.fj.fastq"
  gzip b_joined/${samples[S]}"_R2.fj.fastq"
  gzip b_joined/${samples[S]}"_joined.fj.fastq"
  while (($G<$lengthG)); do
    #C# Mapping reads on reference genome
    bwa mem -t 16 \
    ${genomes[G]}".fasta" \
    b_joined/${samples[S]}"_R1.fj.fastq.gz" \
    b_joined/${samples[S]}"_R2.fj.fastq.gz" \
    > c_mappingB/${samples[S]}"_"${genomes[G]}"_mapped.sam"
    samtools view -b -f2\
    c_mappingB/${samples[S]}"_"${genomes[G]}"_mapped.sam" \
    > c_mappingB/${samples[S]}"_"${genomes[G]}"_mapped.bam"
    samtools sort \
    c_mappingB/${samples[S]}"_"${genomes[G]}"_mapped.bam" \
    > c_mappingB/${samples[S]}"_"${genomes[G]}"_mapped.cov.bam"
    #D.1# Exctracting coverage data
    bedtools genomecov \
    -ibam c_mappingB/${samples[S]}"_"${genomes[G]}"_mapped.cov.bam" \
    -bg \
    > coverage/${samples[S]}"_"${genomes[G]}"_mapped.cov.txt"
    #D.2# Filtering bam file to keep reads with shortclips but without indel and hardcliping
    samtools view \
    -h c_mappingB/${samples[S]}"_"${genomes[G]}"_mapped.cov.bam" \
    | awk '{if($0 ~ /^@/ || $6 ~ /S/) {print $0}}' \
    | awk '{if($0 ~ /^@/ || $6 !~ /I/) {print $0}}' \
    | awk '{if($0 ~ /^@/ || $6 !~ /D/) {print $0}}' \
    | awk '{if($0 ~ /^@/ || $6 !~ /H/) {print $0}}' \
    | samtools view -Sb - > d_SC/${samples[S]}"_"${genomes[G]}"_mapped.US.bam"
    samtools sort \
    d_SC/${samples[S]}"_"${genomes[G]}"_mapped.US.bam" \
    > d_SC/${samples[S]}"_"${genomes[G]}"_mapped.SC.bam"
    #E# Exctracting softclip sequences
    ~/SE-MEI/extractSoftclipped \
    d_SC/${samples[S]}"_"${genomes[G]}"_mapped.SC.bam" \
    > e_SCfiltering/${samples[S]}"_"${genomes[G]}"_SC.fq.gz"
    gunzip e_SCfiltering/${samples[S]}"_"${genomes[G]}"_SC.fq.gz"
    #F# Filtering softclip sequences (SC>20 bases)
    awk \
    'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 20) {print header, seq, qheader, qseq}}' \
    < e_SCfiltering/${samples[S]}"_"${genomes[G]}"_SC.fq" \
    > e_SCfiltering/${samples[S]}"_"${genomes[G]}"_SC.flt.fq"
    #G# Mapping softclip sequences on reference genome
    bowtie2 \
    --no-mixed \
    -x "Atha_"${genomes[G]}".bt2" \
    -U e_SCfiltering/${samples[S]}"_"${genomes[G]}"_SC.flt.fq" \
    -S g_mappingBt/${samples[S]}"_"${genomes[G]}".SC.sam"
    samtools view -b \
    g_mappingBt/${samples[S]}"_"${genomes[G]}".SC.sam" \
    > g_mappingBt/${samples[S]}"_"${genomes[G]}".SC.US.bam"
    samtools sort  \
    g_mappingBt/${samples[S]}"_"${genomes[G]}".SC.US.bam" \
    > g_mappingBt/${samples[S]}"_"${genomes[G]}".SC.bam"
    #H# Bam file convertion to .txt to study it with R
    samtools view \
    g_mappingBt/${samples[S]}"_"${genomes[G]}".SC.bam" \
    > recombination/${samples[S]}"_"${genomes[G]}".SC.txt"
    let "G=G+1"
  done
  let "S=S+1"
done
