#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



### PLOTTING THE RESULTS OF ANTIBODY NEUTRALIZATION ASSAY
# UPDATED: 23-Feb-2022 by Nicholas R. Minor
# ----------------------------------------------------------------- #

# This script will import the neutralization data and plot the results.

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



### ANTIBODY DATA ####
### ------------- #
antibody <- read.csv("data/antibody_potency.csv")



### PLOTTING ###
### -------- #
pdf(file = "visuals/neutralization_assay.pdf", 
    width = 8, height = 7)
palette = c("#FFD300", "#0C328A")
par(bty="n") # deleting the box
gap.plot(c(1:4), antibody$Chronic, 
         gap=c(1.5,9.5), 
         ylim = c(0,10), ytics = c(0.0, 0.5, 1.0, 10),
         xtics = c(1,2,3,4), xticlab = row.names(antibody),
         xlab = NA, ylab = "IC99 values µg/ml",
         cex.axis = 0.8,
         col = rep(palette[1], 4), pch = 2, lwd = 2, cex = 1.2)
grid()
gap.plot(c(1:4), antibody$S.614G, 
         gap=c(1.5,9.5), 
         ylim = c(0,10), ytics = c(0.0, 0.5, 1.0, 10),
         xtics = c(1,2,3,4), xticlab = row.names(antibody),
         add = TRUE,
         cex.axis = 0.8,
         col = rep(palette[2], 4), pch = 0, lwd = 2, cex = 1.2)
legend(3, 2.1,
       legend = c("S-614G Control Virus", "Chronic Infection Virus"),
       pch = c(0, 2),
       col = c(palette[2], palette[1]),
       bty = "n", cex = 0.8, bg="transparent", xpd = T)
abline(h=seq(1.45,1.575,.001), col="white")
axis.break(2,1.5,style="slash")
dev.off()

