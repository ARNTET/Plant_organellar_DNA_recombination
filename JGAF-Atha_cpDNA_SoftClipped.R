library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(scales)
library(taRifx)

setwd ("/Users/ArnaudFERTET/ARA_RADA_reMapping_20")

##### Loading files

WT1 <- fread("JGAF-WT-1.SC.txt", fill = TRUE)
WT2 <- fread("JGAF-WT-2.SC.txt", fill = TRUE)
radA_1 <- fread("JGAF-radA-1-1.SC.txt", fill = TRUE)
radA_2  <- fread("JGAF-radA-1-2.SC.txt", fill = TRUE)
radA_3  <- fread("JGAF-radA-1-3.SC.txt", fill = TRUE)
radA_4  <- fread("JGAF-radA-1-4.SC.txt", fill = TRUE)


###########################################################

### WT1

table.WT1 <- data.table()
table.WT1$V1 <- WT1$V4
table.WT1$V2 <- nchar(WT1$V10)
table.WT1$V3 <- WT1$V10
table.WT1 <- table.WT1[order(table.WT1$V1, table.WT1$V2),]
table.WT1_a <- table.WT1[V1 >= 73 & V1 <= 84170]
table.WT1_b <- table.WT1[V1 >= 84370 & V1 <= 110359]
table.WT1_c <- table.WT1[V1 >= 110434 & V1 <= 128128]
table.WT1_NOiso <- rbind(table.WT1_a, table.WT1_b, table.WT1_c)
table.WT1_NOiso$V4 <- 1
a = 0
b = 10000
covSC_WT1 <- data.table()
a_b_mean <- data.table()
while (a < 130000) {
  if (b < 130000) {
    a_b <- table.WT1_NOiso[V1 >= a & V1 <= b]
    a_b_mean <- data.table(b,(sum(a_b$V4)))
    covSC_WT1 <- rbind(covSC_WT1, a_b_mean)
    a = b
    b = b + 10000 }
  else  {
    b = 128128
    a_b <- table.WT1_NOiso[V1 >= a & V1 <= b]
    a_b_mean <- data.table(b,(sum(a_b$V4)))
    covSC_WT1 <- rbind(covSC_WT1, a_b_mean)
    a = 130000 }}

### WT2

table.WT2 <- data.table()
table.WT2$V1 <- WT2$V4
table.WT2$V2 <- nchar(WT2$V10)
table.WT2$V3 <- WT2$V10
table.WT2 <- table.WT2[order(table.WT2$V1, table.WT2$V2),]
table.WT2_a <- table.WT2[V1 >= 73 & V1 <= 84170]
table.WT2_b <- table.WT2[V1 >= 84370 & V1 <= 110359]
table.WT2_c <- table.WT2[V1 >= 110434 & V1 <= 128128]
table.WT2_NOiso <- rbind(table.WT2_a, table.WT2_b, table.WT2_c)
table.WT2_NOiso$V4 <- 1
a = 0
b = 10000
covSC_WT2 <- data.table()
a_b_mean <- data.table()
while (a < 130000) {
  if (b < 130000) {
    a_b <- table.WT2_NOiso[V1 >= a & V1 <= b]
    a_b_mean <- data.table(b,(sum(a_b$V4)))
    covSC_WT2 <- rbind(covSC_WT2, a_b_mean)
    a = b
    b = b + 10000 }
  else  {
    b = 128128
    a_b <- table.WT2_NOiso[V1 >= a & V1 <= b]
    a_b_mean <- data.table(b,(sum(a_b$V4)))
    covSC_WT2 <- rbind(covSC_WT2, a_b_mean)
    a = 130000 }}

### radA_1

table.radA_1 <- data.table()
table.radA_1$V1 <- radA_1$V4
table.radA_1$V2 <- nchar(radA_1$V10)
table.radA_1$V3 <- radA_1$V10
table.radA_1 <- table.radA_1[order(table.radA_1$V1, table.radA_1$V2),]
table.radA_1_a <- table.radA_1[V1 >= 73 & V1 <= 84170]
table.radA_1_b <- table.radA_1[V1 >= 84370 & V1 <= 110359]
table.radA_1_c <- table.radA_1[V1 >= 110434 & V1 <= 128128]
table.radA_1_NOiso <- rbind(table.radA_1_a, table.radA_1_b, table.radA_1_c)
table.radA_1_NOiso$V4 <- 1
a = 0
b = 10000
covSC_radA_1 <- data.table()
a_b_mean <- data.table()
while (a < 130000) {
  if (b < 130000) {
    a_b <- table.radA_1_NOiso[V1 >= a & V1 <= b]
    a_b_mean <- data.table(b,(sum(a_b$V4)))
    covSC_radA_1 <- rbind(covSC_radA_1, a_b_mean)
    a = b
    b = b + 10000 }
  else  {
    b = 128128
    a_b <- table.radA_1_NOiso[V1 >= a & V1 <= b]
    a_b_mean <- data.table(b,(sum(a_b$V4)))
    covSC_radA_1 <- rbind(covSC_radA_1, a_b_mean)
    a = 130000 }}

### radA_2

table.radA_2 <- data.table()
table.radA_2$V1 <- radA_2$V4
table.radA_2$V2 <- nchar(radA_2$V10)
table.radA_2$V3 <- radA_2$V10
table.radA_2 <- table.radA_2[order(table.radA_2$V1, table.radA_2$V2),]
table.radA_2_a <- table.radA_2[V1 >= 73 & V1 <= 84170]
table.radA_2_b <- table.radA_2[V1 >= 84370 & V1 <= 110359]
table.radA_2_c <- table.radA_2[V1 >= 110434 & V1 <= 128128]
table.radA_2_NOiso <- rbind(table.radA_2_a, table.radA_2_b, table.radA_2_c)
table.radA_2_NOiso$V4 <- 1
a = 0
b = 10000
covSC_radA_2 <- data.table()
a_b_mean <- data.table()
while (a < 130000) {
  if (b < 130000) {
    a_b <- table.radA_2_NOiso[V1 >= a & V1 <= b]
    a_b_mean <- data.table(b,(sum(a_b$V4)))
    covSC_radA_2 <- rbind(covSC_radA_2, a_b_mean)
    a = b
    b = b + 10000 }
  else  {
    b = 128128
    a_b <- table.radA_2_NOiso[V1 >= a & V1 <= b]
    a_b_mean <- data.table(b,(sum(a_b$V4)))
    covSC_radA_2 <- rbind(covSC_radA_2, a_b_mean)
    a = 130000 }}

### radA_3

table.radA_3 <- data.table()
table.radA_3$V1 <- radA_3$V4
table.radA_3$V2 <- nchar(radA_3$V10)
table.radA_3$V3 <- radA_3$V10
table.radA_3 <- table.radA_3[order(table.radA_3$V1, table.radA_3$V2),]
table.radA_3_a <- table.radA_3[V1 >= 73 & V1 <= 84170]
table.radA_3_b <- table.radA_3[V1 >= 84370 & V1 <= 110359]
table.radA_3_c <- table.radA_3[V1 >= 110434 & V1 <= 128128]
table.radA_3_NOiso <- rbind(table.radA_3_a, table.radA_3_b, table.radA_3_c)
table.radA_3_NOiso$V4 <- 1
a = 0
b = 10000
covSC_radA_3 <- data.table()
a_b_mean <- data.table()
while (a < 130000) {
  if (b < 130000) {
    a_b <- table.radA_3_NOiso[V1 >= a & V1 <= b]
    a_b_mean <- data.table(b,(sum(a_b$V4)))
    covSC_radA_3 <- rbind(covSC_radA_3, a_b_mean)
    a = b
    b = b + 10000 }
  else  {
    b = 128128
    a_b <- table.radA_3_NOiso[V1 >= a & V1 <= b]
    a_b_mean <- data.table(b,(sum(a_b$V4)))
    covSC_radA_3 <- rbind(covSC_radA_3, a_b_mean)
    a = 130000 }}

### radA_4

table.radA_4 <- data.table()
table.radA_4$V1 <- radA_4$V4
table.radA_4$V2 <- nchar(radA_4$V10)
table.radA_4$V3 <- radA_4$V10
table.radA_4 <- table.radA_4[order(table.radA_4$V1, table.radA_4$V2),]
table.radA_4_a <- table.radA_4[V1 >= 73 & V1 <= 84170]
table.radA_4_b <- table.radA_4[V1 >= 84370 & V1 <= 110359]
table.radA_4_c <- table.radA_4[V1 >= 110434 & V1 <= 128128]
table.radA_4_NOiso <- rbind(table.radA_4_a, table.radA_4_b, table.radA_4_c)
table.radA_4_NOiso$V4 <- 1
a = 0
b = 10000
covSC_radA_4 <- data.table()
a_b_mean <- data.table()
while (a < 130000) {
  if (b < 130000) {
    a_b <- table.radA_4_NOiso[V1 >= a & V1 <= b]
    a_b_mean <- data.table(b,(sum(a_b$V4)))
    covSC_radA_4 <- rbind(covSC_radA_4, a_b_mean)
    a = b
    b = b + 10000 }
  else  {
    b = 128128
    a_b <- table.radA_4_NOiso[V1 >= a & V1 <= b]
    a_b_mean <- data.table(b,(sum(a_b$V4)))
    covSC_radA_4 <- rbind(covSC_radA_4, a_b_mean)
    a = 130000 }}


###########################################################

#Histo

covSC_WT1$Sample <- "WT1"
covSC_WT2$Sample <- "WT2"
covSC_36$Sample <- "Mut36"
covSC_37$Sample <- "Mut37"
covSC_39$Sample <- "Mut39"
covSC_40$Sample <- "Mut40"

covSC.tot <- rbind(covSC_WT1,covSC_WT2,covSC_36,covSC_37,covSC_39,covSC_40)
covSC.tot$b <- ceiling(covSC.tot$b/1000)

names(covSC.tot)[names(covSC.tot) == "b"] = "Position"
names(covSC.tot)[names(covSC.tot) == "V2"] = "Mapped_Soft_clip"

options(scipen=5)

ggplot(covSC.tot, aes(x = Position, y = Mapped_Soft_clip))+
  geom_bar(
    aes(fill = Sample), stat = "identity", color = "white",
    position = "dodge") +
  fill_palette("jco") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 128),
                     labels = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100", "100-110", "110-120", "120-128")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))