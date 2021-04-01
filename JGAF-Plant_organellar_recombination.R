library(readr)
library(data.table)
library(ggplot2)
library(scales)
library(taRifx)
library(hrbrthemes)
library(BBmisc)
library(gridExtra)
library(magrittr)
library(dplyr)
library(plyr)

# /!\ To fill /!\ ##############################################################
setwd ("")
masterSamples <- c("") #
samples <- c("")
genomes <- c("")
coverages <- fread("cov.csv", fill = TRUE) # see JGAF-Plant_organellar_coverage.R
WildTypes <- c("")

## Step zero (to do the first time only) ########################################
# Simplification of data table to make them usable

S <- 1
while (S <= length(samples)) {
  G <- 1
  while (G <= length(genomes)) {
    Sa <- paste(samples[S],genomes[G], sep="_")
    Sb <- paste(Sa, '.tot', sep="")
    Sc <- paste(Sa, '.SC.txt', sep="")
    Sd <- fread(Sc, fill = TRUE)
    C1 <- data.table(Sd$V1)
    names(C1)[names(C1) == "V1"] = "V1"
    C4 <- data.table(Sd$V4)
    names(C4)[names(C4) == "V1"] = "V2"
    CSa <- cbind(C1,C4)
    Sb <- apply(CSa,2,as.character)
    Sc <- paste(Sa, '.csv', sep="")
    write.table(Sb,Sc, row.names = FALSE, sep = ";")
    G <- G+1}
  rm (Sa,Sb,Sd,Sc,C1,C4,Ca)
  S <- S+1}

print('Now you need to clean .csv files to remove the name of the genome and to keep only read and shortclip mapping positions')

## Step one #####################################################################
# Building individual plot per sample
S <- 1
while (S <= length(samples)) {
  G <- 1
  while (G <= length(genomes)) {
    Sa <- paste(samples[S],genomes[G], sep="_")
    Sb <- paste(Sa,".csv", sep="")
    Sc <- fread(Sb, fill = TRUE)
    names(Sc)[names(Sc) == "V1"] = "Pos RD" #Stands for read mapping positions
    names(Sc)[names(Sc) == "V2"] = "Pos SC" #Stands for shortclip mapping positions
    Sd <- paste(Sa, ".flt.csv", sep="")
    write.table(Sc,Sd, row.names = FALSE, sep = ";")
    #A# Position are rounded to the upper kb position
    Se <- data.table(ceiling(Sc$`Pos RD`/1000))
    Sf <- data.table(ceiling(Sc$`Pos SC`/1000))
    names(Se)[names(Se) == "V1"] = "Pos RD"
    names(Sf)[names(Sf) == "V1"] = "Pos SC"
    Sg <- cbind(Se,Sf)
    #B# Counting events per pairs of kb
    Sg$V3 <- 1 
    names(Sg)[names(Sg) == "V3"] = "Count"
    Si <- data.table(Sg)
    Si <- Si[, lapply(.SD, sum), by=list(`Pos RD`, `Pos SC`)]
    Sh <- sum(Sg$Count)
    assign(paste(Sa,".flt", sep=""), data.table(Si))
    assign(paste(Sa,".events", sep=""), Sh)
    Sj <- paste(Sa, '.sum.csv', sep="")
    write.table(Si,Sj, row.names = FALSE, sep = ";")
    #C# Plot
    ggMUT <- ggplot(Si, aes(`Pos RD`, `Pos SC`, fill= Count)) + 
      geom_tile() +
      ggtitle(Sa) +
      xlab("Read mapping position (kbp)") +
      ylab("Short-clip mapping position (kbp)") +
      scale_fill_gradient(low="#015EA9", high="#F4B826") +
      coord_fixed()
    Sk <- paste(Sa, '.plot.pdf', sep="")
    ggsave(Sk, ggMUT)
    G <-G+1 }
  rm (Sa,Sb,Sd,Sh,Sj,Sk,Sc,Se,Sf,Sg,Si)
  S <- S+1 }

## Step two #####################################################################
# Adjusting count by coverage to make data comparable between samples

S <- 1
while (S <= length(samples)) {
  G <- 1
  while (G <= length(genomes)) {
    Sa <- paste(samples[S],genomes[G], sep="_")
    Sb <- paste(Sa,".sum.csv", sep="")
    Sc <- fread(Sb, fill = TRUE)
    if (G %% 2 == 0) {
      Sc$V4 <- (Sc$Count / coverages[S,mtDNA])}
    else{
      Sc$V4 <- (Sc$Count / coverages[S,cpDNA])}
    names(Sc)[names(Sc) == "V4"] = "Adj. Count" #Stands for count adjusted with respect to the coverage
    write.table(Sc,Sb, row.names = FALSE, sep = ";")
    G <-G+1 }
  rm (Sa,Sb,Sc,G)
  S <- S+1 }

## Step three ###################################################################
# Adjusting count by coverage to make data comparable between samples

G <- 1
while (G <= length(genomes)) {
  W <- 1
  #A# Defining thresholds regarding genomes
  if (G %% 2 == 0){
    trS <- 0.1}
  else{
    trS <- 0.015}
  #B# Merging Wild Types
  while (W <= length(WildTypes)) {
    if (W %% 2 == 0) {
      Sa <- paste(WildTypes[W],genomes[G], sep="_")
      Sb <- paste(Sa,".sum.csv", sep="")
      Sc <- fread(Sb, fill = TRUE)
      setkeyv(Sc, c("Pos RD", "Pos SC"))}
    else{
      Sa <- paste(WildTypes[W],genomes[G], sep="_")
      Sb <- paste(Sa,".sum.csv", sep="")
      Sd <- fread(Sb, fill = TRUE)
      setkeyv(Sd, c("Pos RD", "Pos SC"))}
    W <- W +1}
  Se <- merge(Sd, Sc, all=TRUE,allow.cartesian=TRUE)
  Se[is.na(Se)] <- 0
  Se$V6 <- Se$`Adj. Count.x` + Se$`Adj. Count.y`
  #C# Filtering according to the threshold for adjusted count
  Se <- Se[V6 > trS]
  Sf <- Se$`Pos RD`
  Sg <- Se$`Pos SC`
  Sh <- data.table(cbind(Sf,Sg))
  names(Sh)[names(Sh) == "Sf"] = "Pos RD"
  names(Sh)[names(Sh) == "Sg"] = "Pos SC"
  Sb <- paste(genomes[G],".out", sep="") 
  assign(Sb, Sh)
  G <- G + 1}

## Step four ####################################################################
# Building individual plot per sample without outlayers

S <- 1
outlayers <- c()
while (S <= length(samples)) {
  G <- 1
  #A# Defining limits for plots  
  while (G <= length(genomes)) {
    if (G %% 2 == 0) {
      outlayers <- mtDNA.out
      lIm <- 37 }
    else{
      outlayers <- cpDNA.out
      lIm <- 13}
    Sa <- paste(samples[S],genomes[G], sep="_")
    Sb <- paste(Sa,".sum.csv", sep="")
    Sc <- fread(Sb, fill = TRUE)
    #B# Defining and removing outlayer positions
    c <- 1
    while(c <= length(Sc$Count)){
      a <- 1
      b <- 2
      while(a <= length(outlayers$`Pos RD`)){
        if(Sc[c,1] == outlayers[a,`Pos RD`] && Sc[c,2] == outlayers[a,`Pos SC`]){
          Sc[c,4] = 0}
        else{}
        a <- a + 1}
      c <- c+1 }
    Sb <- paste(Sa, ".noOut.csv", sep="")
    write.table(Sc,Sb, row.names = FALSE, sep = ";")
    assign(paste(Sa,".noOut", sep=""), data.table(Sc))
    Sc <- Sc[`Adj. Count` >= 0.000000001]
    #C# Plot
    ggMUT <- ggplot(Sc, aes(`Pos RD`, `Pos SC`, fill= `Adj. Count`)) + 
      geom_tile() +
      ggtitle(Sa) +
      xlab("Read mapping position (kbp)") +
      xlim(0,lIm) +
      ylab("Short-clip mapping position (kbp)") +
      ylim(0,lIm) +
      scale_fill_gradient(low="#015EA9", high="#F4B826") +
      coord_fixed()
    assign(paste(Sa, '.plot.noOut', sep=""), ggMUT)
    Sk <- paste(Sa, '.plot.noOut.pdf', sep="")
    ggsave(Sk, ggMUT)
    G <- G+1}
  rm (Sa,Sb,Sc,ggMUT,Sk,a,b,c,G)
  S <- S+1}

## Step five ####################################################################
# Building global plot per mutant


#A# Merge all data in a single table 
G <- 1
while (G <= length(genomes)) {
  S <- 1
  Sd <- data.table()
  Sf <- data.table()
  while (S <= length(samples)) {
    Sa <- paste(samples[S],genomes[G], sep="_")
    Sb <- fread(paste(Sa,".noOut.csv", sep=""), fill = TRUE)
    Sc <- Sb[`Adj. Count` >= 0.000000001]
    Sa <- paste(samples[S],genomes[G], sep="_")
    Sb <- fread(paste(Sa,".noOut.csv", sep=""), fill = TRUE)
    Sc <- Sb[`Adj. Count` >= 0.000000001]
    Sc$V5 = paste(samples[S])
    names(Sc)[names(Sc) == "V5"] = "Sample"
    Sc$Count <- NULL
    #B# Adding a mark according to the mutant (here we volounteerly let our mutant names to help to understand the process)
    if (S == 1) { #
      Sd <- rbind(Sd,Sc)
      SdN <- paste("oex1",genomes[G],sep="_")}
    else if (S >= 2 && S <= 5) {
      Sf <- rbind(Sf,Sc)
      SfN <- paste("radA",genomes[G],sep="_")}
    else { #This step corresponds to the WT, so reported to both tables
      Sd <- rbind(Sd,Sc)
      Sf <- rbind(Sf,Sc)}
    S <- S+1}
  assign(paste(SdN), Sd)
  SdS <- paste(SdN, ".tot.csv",sep="")
  write.table(Sd,SdS, row.names = FALSE, sep = ";")
  assign(paste(SfN), Sf)
  SfS <- paste(SfN, ".tot.csv",sep="")
  write.table(Sf,SfS, row.names = FALSE, sep = ";")
  G <- G+1}

#C# Plot
M <- 1
while (M <= length(masterSamples)) {
  G <- 1
  while (G <= length(genomes)) {
    if (G %% 2 == 0) {
      lIm <- 37 }
    else{
      lIm <- 13}
    Sa <- paste(masterSamples[M],genomes[G], sep="_")
    Sb <- fread(paste(Sa,".tot.csv", sep=""), fill = TRUE)
    ggMUT <- ggplot(Sb, aes(`Pos RD`, `Pos SC`, fill= `Adj. Count`)) + 
      geom_tile() +
      xlab("Read mapping position (kbp)") +
      xlim(0,lIm) +
      ylab("Short-clip mapping position (kbp)") +
      ylim(0,lIm) +
      scale_fill_gradient(low="#015EA9", high="#F4B826") +
      coord_fixed() + 
      facet_wrap(~ Sample, ncol=3)
    Sk <- paste(Sa, '.totPlot.pdf', sep="")
    ggsave(Sk, ggMUT)
    G <- G+1 }
  M <- M+1 }
  
