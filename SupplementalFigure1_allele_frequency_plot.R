### PLOTTING INTRAHOST SINGLE-NUCLEOTIDE VARIANT (iSNV) FREQUENCES
# UPDATED: 04-Jan-2022 by Nicholas R. Minor
# ----------------------------------------------------------------- #

# This script will do the following:

# 1) Import some useful information from the consensus sequence FASTA file.
# 2) Import, parse, and reformat information in VCFs called from the raw reads
# from each sequencing time point.
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
list.of.packages <- c("Biostrings", "tidyverse", "RColorBrewer")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Biostrings)
library(tidyverse)
library(RColorBrewer)
data_filepath = ""  ### INSERT YOUR FILE PATH HERE ###
setwd(data_filepath)



### IMPORTING PATIENT DATA ####
### ---------------------- #
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



### IMPORTING TIME POINTS ####
### --------------------- #
timepoint1 <- read.delim("data/individual_USA_WI-WSLH-202168_2020_ONT.vcf", skip = 55)
timepoint1$DATE <- as.Date(fasta_df$Date[1])
timepoint1$DAY <- as.Date(fasta_df$day_of_infection[1])
colnames(timepoint1)[10] <- "Variants"
timepoint1 <- timepoint1[timepoint1$FILTER=="PASS",]

timepoint2 <- read.delim("data/individual_USA_WI-UW-2731_2021_DHO_25021_ONT.vcf", skip = 55)
timepoint2$DATE <- as.Date(fasta_df$Date[2])
timepoint2$DAY <- as.Date(fasta_df$day_of_infection[2])
colnames(timepoint2)[10] <- "Variants"
timepoint2 <- timepoint2[timepoint2$FILTER=="PASS",]

timepoint3 <- read.delim("data/individual_USA_WI-UW-2731-T2_2021_DHO_25063_ONT.vcf", skip = 55)
timepoint3$DATE <- as.Date(fasta_df$Date[3])
timepoint3$DAY <- as.Date(fasta_df$day_of_infection[3])
colnames(timepoint3)[10] <- "Variants"
timepoint3 <- timepoint3[timepoint3$FILTER=="PASS",]

timepoint4 <- read.delim("data/individual_USA_WI-UW-2731-T3_2021_DHO_25216_ONT.vcf", skip = 55)
timepoint4$DATE <- as.Date(fasta_df$Date[4])
timepoint4$DAY <- as.Date(fasta_df$day_of_infection[4])
colnames(timepoint4)[10] <- "Variants"
timepoint4 <- timepoint4[timepoint4$FILTER=="PASS",]

timepoint5 <- read.delim("data/individual_USA_WI-UW-2731-T4_2021_DHO_26156_ONT.vcf", skip = 55)
timepoint5$DATE <- as.Date(fasta_df$Date[5])
timepoint5$DAY <- as.Date(fasta_df$day_of_infection[5])
colnames(timepoint5)[10] <- "Variants"
timepoint5 <- timepoint5[timepoint5$FILTER=="PASS",]

timepoint6 <- read.delim("data/june_timepoint_ARTICv3_vs_MIDNIGHT/junetimepoint_artic_midnight_cat_variants_20211129.vcf", skip = 55)
timepoint6$DATE <- as.Date(fasta_df$Date[6])
timepoint6$DAY <- as.Date(fasta_df$day_of_infection[6])
colnames(timepoint6)[10] <- "Variants"
timepoint6 <- timepoint6[timepoint6$FILTER=="PASS",]

timepoint7 <- read.delim("data/individual_WI-WSLH-217727_Illumina.vcf", skip = 55)
timepoint7$DATE <- as.Date(fasta_df$Date[7])
timepoint7$DAY <- as.Date(fasta_df$day_of_infection[7])
colnames(timepoint7)[10] <- "Variants"
timepoint7 <- timepoint7[timepoint7$FILTER=="PASS",]

timepoint8 <- read.delim("data/individual_USA_Mayo-0928_2021_IonTorrent.vcf", skip = 55)
timepoint8$DATE <- as.Date(fasta_df$Date[8])
timepoint8$DAY <- as.Date(fasta_df$day_of_infection[8])
colnames(timepoint8)[10] <- "Variants"
timepoint8 <- timepoint8[timepoint8$FILTER=="PASS",]

full_vcf <- rbind(timepoint1, timepoint2, timepoint3, timepoint4,
                  timepoint5, timepoint6, timepoint7, timepoint8)

write.csv(full_vcf, "data/alltimepoints_variants_20211214.vcf")
# full_vcf <- read.delim("data/alltimepoints_variants_20211214.vcf", skip = 55)
# str(full_vcf)
# full_vcf$DATE <- as.Date(full_vcf$DATE)



### SEPARATING OUT DEPTHS AND FREQUENCIES ####
### ------------------------------------- #
full_vcf <- full_vcf[grepl("AF",full_vcf$INFO, fixed=T),] ; rownames(full_vcf) <- NULL

for (i in 1:nrow(full_vcf)){
  print(paste("processing row", i, sep = " "))
  info <- full_vcf[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  full_vcf[i, "DEPTH"] <- info_split[grep("^DP", info_split)]
  full_vcf[i, "FREQ"] <- info_split[grep("^AF", info_split)]
}
full_vcf$DEPTH <- as.numeric(str_replace_all(full_vcf$DEPTH, "DP=", ""))
full_vcf$FREQ <- as.numeric(str_replace_all(full_vcf$FREQ, "AF=", ""))



### IDENTIFYING ALL UNIQUE REF_POS_ALT ####
### ---------------------------------- #
# This includes cases where there is more than one mutation at the same
# position. The code in this section will also remove some large but
# now-irrelevant objects from vector memory.
full_vcf$REF_POS_ALT <- paste(full_vcf$REF, full_vcf$POS, full_vcf$ALT, sep = "-")

write.csv(full_vcf, "data/alltimepoints_reduced_VCF_20211214.csv",
          quote=F, row.names=F)
# full_vcf <- read.csv("data/alltimepoints_reduced_VCF_20211214.csv")
# str(full_vcf)
# full_vcf$DATE <- as.Date(full_vcf$DATE)

mutations <- full_vcf[,c("POS", "REF_POS_ALT", "DATE", "DAY", "FREQ", "DEPTH")]
# remove(full_vcf) ; remove(timepoint1) ; remove(timepoint2) ; remove(timepoint3) ; remove(timepoint4) ; remove(timepoint5) ; remove(timepoint6) ; remove(timepoint7) ; remove(timepoint8)
mutations <- mutations[order(mutations$POS),] ; rownames(mutations) <- NULL

mutations$SEQ_PLATFORM <- fasta_df$seq_platform[match(mutations$DATE, fasta_df$Date)]
write.csv(mutations, "data/alltimepoints_mutations_20211214.csv",
          quote=F, row.names=F)
# mutations <- read.csv("/Volumes/working_ssd/e484t_manuscript/data/alltimepoints_mutations_20211214.csv")
str(mutations)
mutations$DATE <- as.Date(mutations$DATE, format = "%Y-%m-%d")
mutations <- mutations[order(mutations$POS),] ; rownames(mutations) <- NULL
mutations_raw <- mutations
E484A <- mutations[mutations$POS==23013,] ; E484A$FREQ <- as.numeric(E484A$FREQ) ; E484A$DEPTH <- as.numeric(E484A$DEPTH)
E484T <- mutations[mutations$POS==23012,] ; E484T$FREQ <- as.numeric(E484T$FREQ) ; E484T$DEPTH <- as.numeric(E484T$DEPTH)



### CONSENSUS MUTATIONS ####
### ------------------- #
consensus_mutations <- mutations_raw
consensus_mutations$KEEP <- NA
for (i in unique(consensus_mutations$REF_POS_ALT)){
  
  print(paste("processing mutation", i, sep = " "))
  
  mut_sub <- consensus_mutations[consensus_mutations$REF_POS_ALT==i,]
  
  if (max(mut_sub$FREQ)<0.5) {
    
    consensus_mutations[consensus_mutations$REF_POS_ALT==i,"KEEP"] <- FALSE
    
  } else {
    consensus_mutations[consensus_mutations$REF_POS_ALT==i,"KEEP"] <- TRUE
  }
}
consensus_mutations <- consensus_mutations[consensus_mutations$KEEP==T,]
rownames(consensus_mutations) <- NULL
consensus_mutations <- consensus_mutations[,-ncol(consensus_mutations)]
write.csv(consensus_mutations, "data/consensus_mutations_20211214.csv", quote=F, row.names = F)
# consensus_mutations <- read.csv("data/consensus_mutations_20211214.csv")
consensus_mutations$DATE <- as.Date(consensus_mutations$DATE, format = "%Y-%m-%d")



### FILTERING ####
### --------- #
ONT <- mutations[mutations$SEQ_PLATFORM=="ONT",]
ONT <- ONT[ONT$FREQ>0.1 & ONT$FREQ < 0.9, ] # 10-90% for ONT and IonTorrent Data
IonTorrent <- mutations[mutations$SEQ_PLATFORM=="IonTorrent",]
IonTorrent <- IonTorrent[IonTorrent$FREQ>0.1 & IonTorrent$FREQ < 0.9, ]
Illumina <- mutations[mutations$SEQ_PLATFORM=="Illumina",]
Illumina <- Illumina[Illumina$FREQ>0.03 & Illumina$FREQ < 0.97, ] # 3-97% for Illumina data
mutations <- rbind(ONT, Illumina, IonTorrent)
mutations <- mutations[order(mutations$POS),] ; rownames(mutations) <- NULL

mutations$KEEP <- NA # The loop below filters all indels and leaves only SNVs
for (i in as.numeric(row.names(mutations))){
  sub <- mutations[i,"REF_POS_ALT"]
  sub_split <- unlist(strsplit(sub, split = "-"))
  print(paste("processing row", i, sep = " "))
  
  if (nchar(sub_split[3])==1 & nchar(sub_split[1])==1){
    mutations$KEEP[i] <- TRUE
  } else {
    mutations$KEEP[i] <- FALSE
  }
}
mutations <- mutations[mutations$KEEP==T,]
rownames(mutations) <- NULL

mutations$KEEP <- (mutations$FREQ*mutations$DEPTH)>10 # this filters to only mutations supported by 10 or more reads
mutations <- mutations[mutations$KEEP==T,]
mutations <- mutations[,-ncol(mutations)]
rownames(mutations) <- NULL

write.csv(mutations, "data/alltimepoints_filtered_mutations_20211214.csv", quote=F, row.names = F)
# mutations <- read.csv("data/alltimepoints_filtered_mutations_20211214.csv")
# str(mutations)
# mutations$DATE <- as.Date(mutations$DATE, format = "%Y-%m-%d")
mutations <- mutations[order(mutations$POS),] ; rownames(mutations) <- NULL



### CREATING DATASET FOR SPIKE 484 iSNVs ####
### ------------------------------------ #
E484A_prefilter <- E484A
E484A <- mutations[mutations$POS==23013,]
E484A <- E484A[E484A$REF_POS_ALT=="A-23013-C",]
E484A <- rbind(E484A, E484A_prefilter[E484A_prefilter$SEQ_PLATFORM=="Illumina",])
E484A <- rbind(E484A, E484A_prefilter[E484A_prefilter$DATE==fasta_df$Date[2] & E484A_prefilter$REF_POS_ALT=="A-23013-C",])
E484A <- E484A[order(E484A$DATE),] ; rownames(E484A) <- NULL

E484T_prefilter <- E484T
E484T <- E484T_prefilter[E484T_prefilter$REF_POS_ALT=="G-23012-A",]
# E484T <- E484T[E484T$REF_POS_ALT=="G-23012-A",]

# E484T <- rbind(E484T, E484T_prefilter[E484T_prefilter$SEQ_PLATFORM=="IonTorrent" & E484T_prefilter$REF_POS_ALT=="G-23012-A",])
E484T <- E484T[order(E484T$DATE),] ; rownames(E484T) <- NULL



### CREATING DATASET OF MUTATIONS REACHING CONSENSUS BY DAY 297 #### 
### ----------------------------------------------------------- #
consensus_by_june <- consensus_mutations
for (i in unique(consensus_by_june$REF_POS_ALT)){
  print(paste("processing mutation", i, sep = " "))
  
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
  print(paste("processing mutation", i, sep = " "))
  
  mut_sub <- consensus_by_june[consensus_by_june$REF_POS_ALT==i,]
  
  if (
    !(as.Date("2021-03-22") %in% mut_sub$DATE) |
    !(as.Date("2021-02-11") %in% mut_sub$DATE) |
    !(as.Date("2021-01-14") %in% mut_sub$DATE) |
    !(as.Date("2021-01-07") %in% mut_sub$DATE) |
    !(as.Date("2020-12-27") %in% mut_sub$DATE)
  ) {
    consensus_by_june <- consensus_by_june[!consensus_by_june$REF_POS_ALT==i,]
    
  } else {
    next
  }
  rownames(consensus_by_june) <- NULL
}
consensus_by_june_step2 <- consensus_by_june
for (i in unique(consensus_by_june$REF_POS_ALT)){
  print(paste("processing mutation", i, sep = " "))
  
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

write.csv(consensus_by_june, "data/consensus_by_june_variants_20211214.csv", quote = F, row.names = F)
# consensus_by_june <- read.csv("data/consensus_by_june_variants_20211214.csv")
# consensus_by_june$DATE <- as.Date(consensus_by_june$DATE, format = "%Y-%m-%d")
E484T_allpoints <- consensus_by_june[consensus_by_june$POS==23012,]



### PLOTTING LOOPS ####
# ----------------- #
pdf(file = "visuals/allele_frequency.pdf", 
    width = 8, height = 5)
# setting the plot frame and grid
plot(mutations$DAY, mutations$FREQ,
     xlim = c(50, 450),
     xlab = "Day of Infection (Dec. 2020 - Sept. 2021)", 
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
    
    if (mut_sub$SEQ_PLATFORM[j]=="ONT"){
      points(mut_sub$DAY[j], mut_sub$FREQ[j], pch = 20, col = "darkgray", cex = 0.8)
    } else if (mut_sub$SEQ_PLATFORM[j]=="Illumina"){
      points(mut_sub$DAY[j], mut_sub$FREQ[j], pch = 17, col = "darkgray", cex = 1)
    } else {
      points(mut_sub$DAY[j], mut_sub$FREQ[j], pch = 15, col = "darkgray", cex = 1)
    }
  }
  
}

E484A_plotting <- E484A
E484A_plotting$FREQ[6:7] <- E484A_plotting$FREQ[6:7] - E484T[E484T$DATE==as.Date("2021-06-29") | E484T$DATE==as.Date("2021-08-04"),
                                                             "FREQ"]

# plotting E484A with confidence intervals
lines(E484A_plotting$DAY, E484A_plotting$FREQ, col="#213CAD", lwd = 2.5)
for (i in 1:nrow(E484A_plotting)){
  
  if (E484A_plotting$SEQ_PLATFORM[i]=="ONT"){
    points(E484A_plotting$DAY[i], E484A_plotting$FREQ[i], pch = 20, col = "#213CAD", cex = 2)
  } else if (E484A_plotting$SEQ_PLATFORM[i]=="Illumina"){
    points(E484A_plotting$DAY[i], E484A_plotting$FREQ[i], pch = 17, col = "#213CAD", cex = 1.5)
  } else {
    points(E484A_plotting$DAY[i], E484A_plotting$FREQ[i], pch = 15, col = "#213CAD", cex = 1.5)
  }
}

# plotting E484T
lines(E484T_allpoints$DAY, E484T_allpoints$FREQ, col="#FEB815", lwd = 2.5)
for (i in 1:nrow(E484A_plotting)){
  
  if (E484T_allpoints$SEQ_PLATFORM[i]=="ONT"){
    points(E484T_allpoints$DAY[i], E484T_allpoints$FREQ[i], pch = 20, col = "#FEB815", cex = 2)
  } else if (E484T_allpoints$SEQ_PLATFORM[i]=="Illumina"){
    points(E484T_allpoints$DAY[i], E484T_allpoints$FREQ[i], pch = 17, col = "#FEB815", cex = 1.5)
  } else {
    points(E484T_allpoints$DAY[i], E484T_allpoints$FREQ[i], pch = 15, col = "#FEB815", cex = 1.5)
  }
}

abline(h=0.1, lty = 4)
abline(h=0.03, lty = 3)
abline(h=0.5, lty = 5)


legend(median(fasta_df$day_of_infection), 1.05, 
       legend=c("Oxford Nanopore", "Illumina", "IonTorrent",
                "E484A", "E484T", "10 iSNVs more frequent in June",
                "ONT frequency cutoff", "Illumina frequency cutoff",
                "consensus frequency cutoff"),
       lwd = c(NA, NA, NA,
               NA, NA, NA, 
               1, 1, 1), 
       lty = c(NA, NA, NA,
               NA, NA, NA,
               4, 3, 5),
       pch = c(20, 17, 15,
               20, 20, 20, 
               NA, NA, NA),
       col = c("gray", "gray", "gray",
               "#213CAD", "#FEB815", "gray", 
               "black", "black", "black"),
       bty = "n", cex = 0.8, bg="transparent", ncol = 3, xpd = T,
       xjust = 0.25, yjust = 0)
dev.off()

