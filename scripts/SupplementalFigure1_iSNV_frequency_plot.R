#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



## PLOTTING INTRAHOST SINGLE-NUCLEOTIDE VARIANT (iSNV) FREQUENCES
# UPDATED: 01-Mar-2022 by Nicholas R. Minor
# ----------------------------------------------------------------- #

# This script will do the following:

# 1) Import some useful information from the consensus sequence FASTA file.
# 2) Import, parse, and reformat information in annotated VCFs called from the
# raw reads from each sequencing time point.
# 3) Filter down the iSNVs to only those that are present at consensus frequency
# (≥ 0.5) by the june time point (post-diagnosis day 297).
# 4) Retrieve the raw frequencies for all those day-297 consensus variants.
# 5) Keep track of all frequencies for Spike E484A and E484T.
# 6) Loop through all these derived datasets and plot.

# It is worth noting that this script takes a very long time--2 days or so--to
# run on a personal machine, and may be worth setting up to run on a remote
# cluster. For this reason, the script also writes intermediate datasets as it
# goes, which we are including in our GitHub repository. To reproduce this
# analysis and visualization more quickly, comment out the code in between the
# read functions throughout the script, and uncomment the read functions. If you
# go this route, you will also need to reformat dates into a date format with
# lines of code such as full_vcf$DATE <- as.Date(full_vcf$DATE), which follow
# the read functions.

# Please direct all questions about the reproducibility of this script to:
# nrminor@wisc.edu

# ----------------------------------------------------------------- #



### PREPARING THE ENVIRONMENT ####
### ------------------------ #
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

list.of.packages <- c("tidyverse", "RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Biostrings)
library(tidyverse)
library(RColorBrewer)
data_filepath = args[1]  ### OR INSERT YOUR FILE PATH HERE ###
setwd(data_filepath)



### IMPORTING PATIENT DATA ####
### ---------------------- #
patient_fasta <- readDNAStringSet("data/alltimepoints_20220222.fasta")
seq_names <- names(patient_fasta)
fasta_df <- data.frame(seq_names) ; remove(patient_fasta)
fasta_df$seq_names <- str_replace_all(fasta_df$seq_names, fixed(" | "), ",")
fasta_df <- separate(data = fasta_df, col = 1,
                     into = c("Sample_ID", "Date"),
                     sep = ",")
fasta_df$Date <- as.Date(fasta_df$Date, "%Y %b %d")
fasta_df$day_of_infection <- c(113,124,131,159,198,297,333,388,415,417,482,486)
fasta_df$seq_platform <- c("ONT", "Illumina", "Illumina", "Illumina", "Illumina", "Illumina", "Illumina", "Illumina", "Illumina", "Illumina", "Illumina", "Illumina")



### IMPORTING TIME POINTS ####
### --------------------- #
# timepoint1 <- read.delim("data/individual_USA_WI-WSLH-202169_2020_ONT.ann.vcf", skip = 60)
# timepoint1$DATE <- as.Date(fasta_df$Date[1])
# timepoint1$DAY <- as.Date(fasta_df$day_of_infection[1])
# colnames(timepoint1)[10] <- "Variants"
# timepoint1 <- timepoint1[timepoint1$FILTER=="PASS",]

timepoint2 <- read.delim("data/individual_timepoint02_day124_Illumina.ann.vcf", skip = 60)
timepoint2$DATE <- as.Date(fasta_df$Date[2])
timepoint2$DAY <- fasta_df$day_of_infection[2]
colnames(timepoint2)[10] <- "Variants"
timepoint2 <- timepoint2[timepoint2$FILTER=="PASS",]

timepoint3 <- read.delim("data/individual_timepoint03_day131_Illumina.ann.vcf", skip = 60)
timepoint3$DATE <- as.Date(fasta_df$Date[3])
timepoint3$DAY <- fasta_df$day_of_infection[3]
colnames(timepoint3)[10] <- "Variants"
timepoint3 <- timepoint3[timepoint3$FILTER=="PASS",]

timepoint4 <- read.delim("data/individual_timepoint04_day159_Illumina.ann.vcf", skip = 60)
timepoint4$DATE <- as.Date(fasta_df$Date[4])
timepoint4$DAY <- fasta_df$day_of_infection[4]
colnames(timepoint4)[10] <- "Variants"
timepoint4 <- timepoint4[timepoint4$FILTER=="PASS",]

timepoint5 <- read.delim("data/individual_timepoint05_day198_Illumina.ann.vcf", skip = 60)
timepoint5$DATE <- as.Date(fasta_df$Date[5])
timepoint5$DAY <- fasta_df$day_of_infection[5]
colnames(timepoint5)[10] <- "Variants"
timepoint5 <- timepoint5[timepoint5$FILTER=="PASS",]

timepoint6 <- read.delim("data/individual_timepoint06_day297_Illumina.ann.vcf", skip = 60)
timepoint6$DATE <- as.Date(fasta_df$Date[6])
timepoint6$DAY <- fasta_df$day_of_infection[6]
colnames(timepoint6)[10] <- "Variants"
timepoint6 <- timepoint6[timepoint6$FILTER=="PASS",]

timepoint7 <- read.delim("data/individual_timepoint07_day333_Illumina.ann.vcf", skip = 60)
timepoint7$DATE <- as.Date(fasta_df$Date[7])
timepoint7$DAY <- fasta_df$day_of_infection[7]
colnames(timepoint7)[10] <- "Variants"
timepoint7 <- timepoint7[timepoint7$FILTER=="PASS",]

timepoint8 <- read.delim("data/individual_timepoint08_day388_Illumina.ann.vcf", skip = 60)
timepoint8$DATE <- as.Date(fasta_df$Date[8])
timepoint8$DAY <- fasta_df$day_of_infection[8]
colnames(timepoint8)[10] <- "Variants"
timepoint8 <- timepoint8[timepoint8$FILTER=="PASS",]

timepoint9 <- read.delim("data/individual_timepoint09_day415_Illumina.ann.vcf", skip = 60)
timepoint9$DATE <- as.Date(fasta_df$Date[9])
timepoint9$DAY <- fasta_df$day_of_infection[9]
colnames(timepoint9)[10] <- "Variants"
timepoint9 <- timepoint9[timepoint9$FILTER=="PASS",]

timepoint10 <- read.delim("data/individual_timepoint10_day417_Illumina.ann.vcf", skip = 60)
timepoint10$DATE <- as.Date(fasta_df$Date[10])
timepoint10$DAY <- fasta_df$day_of_infection[10]
colnames(timepoint10)[10] <- "Variants"
timepoint10 <- timepoint10[timepoint10$FILTER=="PASS",]

timepoint11 <- read.delim("data/individual_timepoint11_day482_Illumina.ann.vcf", skip = 60)
timepoint11$DATE <- as.Date(fasta_df$Date[11])
timepoint11$DAY <- fasta_df$day_of_infection[11]
colnames(timepoint11)[10] <- "Variants"
timepoint11 <- timepoint11[timepoint11$FILTER=="PASS",]

timepoint12 <- read.delim("data/individual_timepoint12_day488_Illumina.ann.vcf", skip = 60)
timepoint12$DATE <- as.Date(fasta_df$Date[12])
timepoint12$DAY <- fasta_df$day_of_infection[12]
colnames(timepoint12)[10] <- "Variants"
timepoint12 <- timepoint12[timepoint12$FILTER=="PASS",]

full_vcf <- rbind(timepoint2, timepoint3, timepoint4,
                  timepoint5, timepoint6, timepoint7, timepoint8,
                  timepoint9, timepoint10, timepoint11, timepoint12)


### SEPARATING OUT DEPTHS, FREQUENCIES, & EFFECTS ####
### -------------------------------------------- #
full_vcf <- full_vcf[full_vcf$QUAL>0,]
full_vcf <- full_vcf[grepl("AF",full_vcf$INFO, fixed=T),] ; rownames(full_vcf) <- NULL
full_vcf <- full_vcf[grepl("DP",full_vcf$INFO, fixed=T),] ; rownames(full_vcf) <- NULL
full_vcf <- full_vcf[grepl("ANN",full_vcf$INFO, fixed=T),] ; rownames(full_vcf) <- NULL

full_vcf$DEPTH <- lapply(full_vcf$INFO, function(i){
  info_split <- unlist(strsplit(i, split = ";"))
  depth <- info_split[grep("^DP", info_split)]
  return(str_remove_all(depth, "DP="))
})
full_vcf$DEPTH <- as.numeric(full_vcf$DEPTH)
full_vcf <- full_vcf[full_vcf$DEPTH>=20, ]

full_vcf$FREQ <- lapply(full_vcf$INFO, function(i){
  info_split <- unlist(strsplit(i, split = ";"))
  freq <- info_split[grep("^AF", info_split)]
  return(str_remove_all(freq, "AF="))
})
full_vcf$FREQ <- as.numeric(full_vcf$FREQ)

full_vcf$ANNOTATIONS <- lapply(full_vcf$INFO, function(i){
  info_split <- unlist(strsplit(i, split = ";"))
  return(info_split[grep("^ANN", info_split)])
})
full_vcf$ANNOTATIONS <- as.character(full_vcf$ANNOTATIONS)

full_vcf$VARIANT_TYPE <- lapply(full_vcf$ANNOTATIONS, function(i){
  ann <- str_replace_all(i, fixed("|"), " ; ")
  ann_split <- unlist(strsplit(ann, split = ";"))
  ann <- str_remove_all(ann_split[2], " ")
  return(str_remove_all(ann, "_variant"))
})
full_vcf$VARIANT_TYPE <- as.character(full_vcf$VARIANT_TYPE)

full_vcf$GENE <- lapply(full_vcf$ANNOTATIONS, function(i){
  ann <- str_replace_all(i, fixed("|"), " ; ")
  ann_split <- unlist(strsplit(ann, split = ";"))
  return(str_remove_all(ann_split[4], " "))
})
full_vcf$GENE <- as.character(full_vcf$GENE)

full_vcf$AA_EFFECT <- lapply(full_vcf$ANNOTATIONS, function(i){
  ann <- str_replace_all(i, fixed("|"), " ; ")
  ann_split <- unlist(strsplit(ann, split = ";"))
  ann <- str_remove_all(ann_split[11], " ")
  ann <- ifelse(ann=="",NA, ann)
  return(str_remove_all(ann, "p."))
})
full_vcf$AA_EFFECT <- as.character(full_vcf$AA_EFFECT)




### IDENTIFYING ALL UNIQUE REF_POS_ALT ####
### ---------------------------------- #
# This includes cases where there is more than one mutation at the same
# position. The code in this section will also remove some large but
# now-irrelevant objects from vector memory.
full_vcf$REF_POS_ALT <- paste(full_vcf$REF, full_vcf$POS, full_vcf$ALT, sep = "-")

write.csv(full_vcf, paste("readables/alltimepoints_variants_",
                          Sys.Date(),
                          ".csv", sep=""), quote=F, row.names=F)
# full_vcf <- read.csv("readables/alltimepoints_variants_2022-03-02.csv")
# str(full_vcf)
# full_vcf$DATE <- as.Date(full_vcf$DATE)

mutations <- full_vcf[,c("POS", "REF_POS_ALT", "GENE", "DATE", "DAY", "FREQ", 
                         "DEPTH", "QUAL", "VARIANT_TYPE", "AA_EFFECT")]
remove(full_vcf) ; remove(timepoint2) ; remove(timepoint3) ; remove(timepoint4) ; remove(timepoint5) ; remove(timepoint6) ; remove(timepoint7) ; remove(timepoint8) ; remove(timepoint9) ; remove(timepoint10) ; remove(timepoint11) ; remove(timepoint12)
mutations <- mutations[order(mutations$POS),] ; rownames(mutations) <- NULL

# mutations$SEQ_PLATFORM <- fasta_df$seq_platform[match(mutations$DATE, fasta_df$Date)]

write.csv(mutations, paste("readables/alltimepoints_variants_reduced_",
                          Sys.Date(),
                          ".csv", sep=""), quote=F, row.names=F)
# mutations <- read.csv("readables/alltimepoints_variants_reduced_2022-03-02.csv")
# str(mutations)
mutations$DATE <- as.Date(mutations$DATE, format = "%Y-%m-%d")
mutations <- mutations[order(mutations$POS),] ; rownames(mutations) <- NULL
mutations_raw <- mutations



### CREATING DATASET FOR SPIKE 484 iSNVs ####
### ------------------------------------ #
E484A <- mutations[mutations$POS==23013,] ; E484A$FREQ <- as.numeric(E484A$FREQ) ; E484A$DEPTH <- as.numeric(E484A$DEPTH)
E484A_prefilter <- E484A
E484A <- E484A[E484A$REF_POS_ALT=="A-23013-C",]
E484A <- E484A[E484A$FREQ>0,]
E484A <- E484A[order(E484A$DATE),] ; rownames(E484A) <- NULL

E484T <- mutations[mutations$POS==23012,] ; E484T$FREQ <- as.numeric(E484T$FREQ) ; E484T$DEPTH <- as.numeric(E484T$DEPTH)
E484T_prefilter <- E484T
E484T <- E484T_prefilter[E484T_prefilter$REF_POS_ALT=="G-23012-A",]
E484T <- E484T[E484T$FREQ>0,]
E484T <- E484T[order(E484T$DATE),] ; rownames(E484T) <- NULL



### CONSENSUS MUTATIONS ####
### ------------------- #
consensus_mutations <- mutations_raw
consensus_mutations <- consensus_mutations[consensus_mutations$DEPTH>=20, ]
consensus_mutations$KEEP <- NA
for (i in unique(consensus_mutations$REF_POS_ALT)){
  
  if (which(unique(consensus_mutations$REF_POS_ALT)==i) %% 1000 == 0){
    print(paste("determining if mutation", i, "reaches consensus-level frequency (0.5)",
              sep = " "))
  }
  
  mut_sub <- consensus_mutations[consensus_mutations$REF_POS_ALT==i,]
  
  if (max(mut_sub$FREQ)<0.5) {
    
    consensus_mutations[consensus_mutations$REF_POS_ALT==i,"KEEP"] <- FALSE
    
  } else {
    consensus_mutations[consensus_mutations$REF_POS_ALT==i,"KEEP"] <- TRUE
  }
}
consensus_mutations <- consensus_mutations[consensus_mutations$KEEP==T,]
rownames(consensus_mutations) <- NULL
consensus_mutations <- consensus_mutations[,-which(colnames(consensus_mutations)=="KEEP")]
write.csv(consensus_mutations,  paste("readables/consensus_mutations_",
                                      Sys.Date(),".csv", sep = ""), 
                                      quote=F, row.names = F)
# consensus_mutations <- read.csv("readables/consensus_mutations_2022-03-02.csv")
consensus_mutations$DATE <- as.Date(consensus_mutations$DATE, format = "%Y-%m-%d")



### CREATING DATASET OF MUTATIONS REACHING CONSENSUS BY DAY 297 #### 
### ----------------------------------------------------------- #
consensus_by_june <- consensus_mutations
consensus_by_june$DETECTED_BEFORE_DAY_297 <- NA
for (i in unique(consensus_by_june$REF_POS_ALT)){

    mut_sub <- consensus_by_june[consensus_by_june$REF_POS_ALT==i,]
  
  if (
    !(as.Date("2021-06-29") %in% mut_sub$DATE)
  ) {
    consensus_by_june <- consensus_by_june[!consensus_by_june$REF_POS_ALT==i,]
    
  } else {
    next
  }
  rownames(consensus_by_june) <- NULL
}
for (i in unique(consensus_by_june$REF_POS_ALT)){

  if (which(unique(consensus_by_june$REF_POS_ALT)==i) %% 100 == 0){
    print(paste("processing mutation", i, sep = " "))
  }

  mut_sub <- consensus_by_june[consensus_by_june$REF_POS_ALT==i,]

  if (min(mut_sub$DAY)==297){

    # consensus_by_june <- consensus_by_june[!consensus_by_june$REF_POS_ALT==i,]
    consensus_by_june[consensus_by_june$REF_POS_ALT==i, "DETECTED_BEFORE_DAY_297"] <- F

  } else {
    # next
    consensus_by_june[consensus_by_june$REF_POS_ALT==i, "DETECTED_BEFORE_DAY_297"] <- T
  }
  rownames(consensus_by_june) <- NULL
}
for (i in unique(consensus_by_june$REF_POS_ALT)){
  
  if (which(unique(consensus_by_june$REF_POS_ALT)==i) %% 100 == 0){
    print(paste("processing mutation", i, sep = " "))
  }
  
  mut_sub <- consensus_by_june[consensus_by_june$REF_POS_ALT==i,]
  
  if (
    mut_sub[mut_sub$DATE==as.Date("2021-06-29"), "FREQ"] < 0.5 |
    T %in% (mut_sub[mut_sub$DATE < as.Date("2021-06-29"), "FREQ"] >= 0.5)
  ) {
    consensus_by_june <- consensus_by_june[!consensus_by_june$REF_POS_ALT==i,]
    
  } else {
    next
  }
  rownames(consensus_by_june) <- NULL
}

consensus_by_june <- consensus_by_june[consensus_by_june$DEPTH>0,]
write.csv(consensus_by_june, paste("readables/consensus_by_june_variants_",
                                   Sys.Date(),".csv", sep = ""),
          quote = F, row.names = F)
# consensus_by_june <- read.csv("readables/consensus_by_june_variants_2022-03-02.csv")
# consensus_by_june$DATE <- as.Date(consensus_by_june$DATE, format = "%Y-%m-%d")

supplemental_table2 <- consensus_by_june
supplemental_table2 <- supplemental_table2[,c("POS", "REF_POS_ALT", "DAY",
                                              "GENE", "VARIANT_TYPE", "AA_EFFECT",
                                              "DETECTED_BEFORE_DAY_297")]
supplemental_table2$NO_OF_DETECTIONS <- 0
supplemental_table2$DAYS_OF_DETECTION <- 0
for (i in unique(supplemental_table2$REF_POS_ALT)){
  
  sub <- supplemental_table2[supplemental_table2$REF_POS_ALT==i,]
  supplemental_table2[supplemental_table2$REF_POS_ALT==i,
                      "NO_OF_DETECTIONS"] <- nrow(sub)
  
  supplemental_table2[supplemental_table2$REF_POS_ALT==i,
                      "DAYS_OF_DETECTION"] <- paste(sub$DAY, collapse = ";")
}
supplemental_table2$keep <- NA
supplemental_table2$keep[1] <- T
for (i in 2:nrow(supplemental_table2)){
  
  if (supplemental_table2[i,"REF_POS_ALT"]==supplemental_table2[i-1,"REF_POS_ALT"]){
    supplemental_table2$keep[i] <- F
  } else {
    supplemental_table2$keep[i] <- T
  }
  
}
supplemental_table2 <- supplemental_table2[supplemental_table2$keep==T,
                                           c("POS", "REF_POS_ALT", "GENE", 
                                             "VARIANT_TYPE", "AA_EFFECT",
                                             "NO_OF_DETECTIONS",
                                             "DAYS_OF_DETECTION",
                                             "DETECTED_BEFORE_DAY_297")]
row.names(supplemental_table2) <- NULL
write.csv(supplemental_table2, "readables/supplemental_table2_june_substitutions.csv",
          row.names = F, quote = F)



### PLOTTING LOOPS ####
# ----------------- #
pdf(file = "visuals/figsupp1_isnv_frequencies.pdf", 
    width = 8, height = 5)
# setting the plot frame and grid
plot(mutations$DAY, mutations$FREQ,
     xlim = c(100, 500),
     xlab = "Days Post-Diagnosis", 
     ylim = c(0,1), ylab = "iSNV Frequency", 
     frame.plot = F, cex.axis = 0.8, cex.lab = 0.85, las = 1,
     pch = 20, col = "darkgray", cex = 0.8, type = 'n')
grid()

# plotting mutations that increase to consensus frequency (50% or higher) at the June time point
for (i in unique(consensus_by_june$REF_POS_ALT)){
  mut_sub <- consensus_by_june[consensus_by_june$REF_POS_ALT==i,]
  
  lines(mut_sub$DAY, mut_sub$FREQ, col="gray")
  # points(mut_sub$DAY, mut_sub$FREQ, pch = 20, col = "darkgray", cex = 0.8)
  
  for (j in 1:nrow(mut_sub)){
    
    points(mut_sub$DAY[j], mut_sub$FREQ[j], pch = 20, col = "darkgray", cex = 0.8)
    
  }
  
}

E484A_plotting <- E484A
for (i in 1:nrow(E484A_plotting)){
  
  if (E484A_plotting$DATE[i] %in% E484T$DATE){
    E484A_plotting$FREQ[i] <- E484A_plotting$FREQ[i] - E484T[E484T$DATE==E484A_plotting$DATE[i], "FREQ"]
  }
  
}

# plotting E484A with confidence intervals
lines(E484A_plotting$DAY, E484A_plotting$FREQ, col="#213CAD", lwd = 2.5)
points(E484A_plotting$DAY, E484A_plotting$FREQ, pch = 20, col = "#213CAD", cex = 2)

# plotting E484T
lines(E484T$DAY, E484T$FREQ, col="#FEB815", lwd = 2.5)
points(E484T$DAY, E484T$FREQ, pch = 20, col = "#FEB815", cex = 2)

# abline(h=0.1, lty = 4)
abline(h=0.03, lty = 3)
abline(h=0.5, lty = 5)


legend(median(fasta_df$day_of_infection), 1.05, 
       legend=c("E484A", "E484T", 
                paste(length(unique(consensus_by_june$REF_POS_ALT)), 
                      "iSNVs more frequent in June", sep = " "),
                "Typical Illumina frequency cutoff",
                "consensus frequency cutoff"),
       lwd = c(NA, NA, NA, 
               1, 1, 1), 
       lty = c(NA, NA, NA,
               3, 5),
       pch = c(20, 20, 20, 
               NA, NA),
       col = c("#213CAD", "#FEB815", "gray", 
               "black", "black"),
       bty = "n", cex = 0.8, bg="transparent", ncol = 2, xpd = T,
       xjust = 0.25, yjust = 0)
dev.off()

