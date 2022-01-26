#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



### PLOTTING THE RESULTS OF ANTIBODY NEUTRALIZATION ASSAY
# UPDATED: 04-Jan-2022 by Nicholas R. Minor
# ----------------------------------------------------------------- #

# This script will import the neutralization data, along with some useful info
# in the consensus FASTA, and plot the results.

# Please direct all questions about the reproducibility of this script to:
# nrminor@wisc.edu

# ----------------------------------------------------------------- #



### SETUP ####
### ----- #
list.of.packages <- c("Biostrings", "tidyverse", "plotrix")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Biostrings)
library(tidyverse)
library(plotrix)
data_filepath = args[1] ### OR INSERT YOUR FILE PATH HERE ###
setwd(data_filepath)



### FASTA PREP ####
### ---------- #
patient_fasta <- readDNAStringSet("data/alltimepoints_20211104.fasta")
seq_names <- names(patient_fasta)
fasta_df <- data.frame(seq_names) ; remove(patient_fasta)
fasta_df$seq_names <- str_replace_all(fasta_df$seq_names, fixed(" | "), ",")
fasta_df <- separate(data = fasta_df, col = 1,
                     into = c("Sample_ID", "Date"),
                     sep = ",")
fasta_df$Date <- as.Date(fasta_df$Date, "%Y %b %d")
fasta_df$day_of_infection <- c(113,124,131,159,198,297,333,388)
fasta_df$seq_platform <- c("ONT", "ONT", "ONT", "ONT", "ONT", "ONT", "Illumina", "IonTorrent")



### ANTIBODY DATA ####
### ------------- #
antibody <- read.csv("data/antibody_potency.csv")



### PLOTTING ###
### -------- #
pdf(file = "/Users/nicholasminor/Documents/informatics/E484T_paper/visuals/neutralization_assay.pdf", 
    width = 8, height = 7)
par(bty="n") # deleting the box
gap.plot(c(1:4), antibody$Chronic, 
         gap=c(1.5,9.5), 
         ylim = c(0,10), ytics = c(0.0, 0.5, 1.0, 10),
         xtics = c(1,2,3,4), xticlab = row.names(antibody),
         xlab = NA, ylab = "IC99 values µg/ml",
         cex.axis = 0.8,
         col = palette$points[2], pch = 2, lwd = 2, cex = 1.2)
grid()
# lines(c(1:4), antibody$S.614G, col = palette$lines[1], lwd = 2)
points(c(1:4), antibody$S.614G, col = palette$points[1], pch = 0, cex = 1.2)

legend(3, 2,
       legend = c("S-614G Control Virus", "Chronic Infection Virus"),
       pch = c(0, 2),
       col = c(palette$points[1], palette$points[2]),
       bty = "n", cex = 0.8, bg="transparent", xpd = T)

# abline(h=seq(1.45,1.575,.001), col="white")
axis.break(2,1.5,style="slash")
dev.off()

