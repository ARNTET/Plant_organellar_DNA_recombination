library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(scales)
library(taRifx)


# /!\ To fill /!\ ###############################################################
setwd("")
samples <- c()
genomes <- c("mtDNA","cpDNA")


## Step one #####################################################################
# Normalization and mean coverage per kb

S <- 1
cov <- data.table()
while (S <= length(samples)) {
  covG <- data.table()
  G <- 1
  Sx <- 0
  #A# Determination of number of mapped reads per sample on all genomes
  while (G <= length(genomes)) {
    Sa <- paste(samples[S],genomes[G], sep="_")
    Sc <- paste(Sa, "_readsCount.txt", sep="")
    Sd <- fread(Sc, fill = TRUE)
    Sx <- Sx + Sd$V1[1]}
  G <- 1
  while (G <= length(genomes)) {
    Sa <- paste(samples[S],genomes[G], sep="_")
    Sb <- paste(Sa, '_mapped.cov.txt', sep="")
    Se <- fread(Sb, fill = TRUE)
    #B# Normalization to 1'000'000 reads mapped
    Se$V5 <- ceiling(Se$V4*(1000000/Sx))
    #C# Defining limits for according to genomes  
    if (G %% 2 == 0) {
      genomeLimit <- 129}
    else{
      genomeLimit <- 370}
    #D# Mean coverage of each 1 kb-range
    a = 0
    b = 1000
    Se_kbCov <- data.table()
    while (a < genomeLimit) {
      a_b <- Se[V3 >= a & V3 <= b]
      a_b_mean <- data.table((b/1000),ceiling(mean(a_b$V5)))
      Se_kbCov <- rbind(Se_kbCov, a_b_mean)
      a = a + 1000
      b = b + 1000}
    write.table(Se_kbCov,paste(Sa, '.kbNrmlzCov.csv', sep=""), row.names = TRUE, sep = ";")
    #E# Building coverage table needed in JGAF-Plant_organellar_recombination.R
    Sf <- data.table(mean(Se$V5))
    covG <- rbind(covG,Sf)
    G <- G+1}
  names(covG)[names(covG) == "V1"] = paste(samples[S])
  cov <- cbind(cov,covG)
  S <- S+1}
write.table(cov,"mean_cov.csv", row.names = TRUE, sep = ";")


## Step two #####################################################################
# Assemble of coverage data into a single table 

G <- 1
while (G <= length(genomes)) {
  St <- data.table()
  S <- 1
  while (S <= length(samples)) {
    Sa <- paste(samples[S],genomes[G], sep="_")
    Sb <- paste(Sa, '.kbNrmlzCov.csv', sep="")
    Sc <- fread(Sb, fill = TRUE)
    St <- cbind(St, Sc$V2)
    names(St)[names(St) == "V2"] = paste(samples[S])
    S <- S+1}
  St$Position <- 1:nrow(St)
  write.table(St,paste(genomes[G], "totCov.csv", sep="_"), row.names = TRUE, sep = ";")
  G <- G+1}


## Step Three ###################################################################
# Plots

G <- 1
while (G <= length(genomes)) {
  Sa <- paste(genomes[G], "totCov.csv", sep="_")
  Sb <- fread(Sa, fill = TRUE)
  plotCOV <- ggplot(data = Sb, mapping = aes(Position)) +
    fill_palette("jco") +
    ylim(0,50) +
    ylab("Coverage depth") +
    xlab("mtDNA coordinates (kbp)")
  S <- 1
  while (S <= length(samples)) {
    plotCOV <- plotCOV + geom_line(aes(y = paste(samples[S]), colour = paste(samples[S])))
    S <- S+1}
  Sk <- paste(Sa, '.CovPlot.pdf', sep="")
  ggsave(Sk, plotCOV)
  G <- G+1}
