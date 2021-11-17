
### PREPARING THE ENVIRONMENT ####
library(Biostrings)
library(tidyverse)
library(RColorBrewer)
data_filepath = "/Volumes/working_ssd/e484t_manuscript/"
setwd(data_filepath)



### IMPORTING PATIENT DATA ####
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
timepoint1 <- read.delim("data/all_timepoints_raw_reads/individual_USA_WI-WSLH-202168_2020_ONT.vcf", skip = 55)
timepoint1$DATE <- as.Date(fasta_df$Date[1])
colnames(timepoint1)[10] <- "Variants"
timepoint1 <- timepoint1[timepoint1$FILTER=="PASS",]

timepoint2 <- read.delim("data/all_timepoints_raw_reads/individual_USA_WI-UW-2731_2021_DHO_25021_ONT.vcf", skip = 55)
timepoint2$DATE <- as.Date(fasta_df$Date[2])
colnames(timepoint2)[10] <- "Variants"
timepoint2 <- timepoint2[timepoint2$FILTER=="PASS",]

timepoint3 <- read.delim("data/all_timepoints_raw_reads/individual_USA_WI-UW-2731-T2_2021_DHO_25063_ONT.vcf", skip = 55)
timepoint3$DATE <- as.Date(fasta_df$Date[3])
colnames(timepoint3)[10] <- "Variants"
timepoint3 <- timepoint3[timepoint3$FILTER=="PASS",]

timepoint4 <- read.delim("data/all_timepoints_raw_reads/individual_USA_WI-UW-2731-T3_2021_DHO_25216_ONT.vcf", skip = 55)
timepoint4$DATE <- as.Date(fasta_df$Date[4])
colnames(timepoint4)[10] <- "Variants"
timepoint4 <- timepoint4[timepoint4$FILTER=="PASS",]

timepoint5 <- read.delim("data/all_timepoints_raw_reads/individual_USA_WI-UW-2731-T4_2021_DHO_26156_ONT.vcf", skip = 55)
timepoint5$DATE <- as.Date(fasta_df$Date[5])
colnames(timepoint5)[10] <- "Variants"
timepoint5 <- timepoint5[timepoint5$FILTER=="PASS",]

timepoint6 <- read.delim("data/all_timepoints_raw_reads/individual_USA_WI-UW-5350_2021_DHO_26076_ONT.vcf", skip = 55)
timepoint6$DATE <- as.Date(fasta_df$Date[6])
colnames(timepoint6)[10] <- "Variants"
timepoint6 <- timepoint6[timepoint6$FILTER=="PASS",]

timepoint7 <- read.delim("data/all_timepoints_raw_reads/individual_WI-WSLH-217727_Illumina.vcf", skip = 55)
timepoint7$DATE <- as.Date(fasta_df$Date[7])
colnames(timepoint7)[10] <- "Variants"
timepoint7 <- timepoint7[timepoint7$FILTER=="PASS",]

timepoint8 <- read.delim("data/all_timepoints_raw_reads/individual_USA_Mayo-0928_2021_IonTorrent.vcf", skip = 55)
timepoint8$DATE <- as.Date(fasta_df$Date[8])
colnames(timepoint8)[10] <- "Variants"
timepoint8 <- timepoint8[timepoint8$FILTER=="PASS",]

full_vcf <- rbind(timepoint1, timepoint2, timepoint3, timepoint4,
                  timepoint5, timepoint6, timepoint7, timepoint8)

# full_vcf <- read.delim("data/all_timepoints_raw_reads/alltimepoints_variants_20211116.vcf", skip = 55)
# full_vcf <- full_vcf[,c(1:9, 16, 14, 11, 12, 13, 15, 17, 10)]
# colnames(full_vcf)[10:17] <- c("timepoint1", "timepoint2", "timepoint3", "timepoint4", "timepoint5", "timepoint6", "timepoint7", "timepoint8")



### SEPARATING OUT DEPTHS AND FREQUENCIES ####
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
# (INCLUDING CASES WHERE THERE IS MORE THAN ONE MUTATION AT THE SAME POSITION)
# ALSO TRIMMING OUT NOW-IRRELEVANT STUFF FROM VECTOR MEMORY
full_vcf$REF_POS_ALT <- paste(full_vcf$REF, full_vcf$POS, full_vcf$ALT, sep = "-")

write_csv(full_vcf, "data/all_timepoints_raw_reads/alltimepoints_reduced_VCF_20211117.csv",
          quote=F)
# full_vcf <- read.csv("data/all_timepoints_raw_reads/alltimepoints_reduced_VCF_20211117.csv")
# str(full_vcf)

mutations <- full_vcf[,c("POS", "REF_POS_ALT", "DATE", "FREQ", "DEPTH")]
remove(full_vcf) ; remove(timepoint1) ; remove(timepoint2) ; remove(timepoint3) ; remove(timepoint4) ; remove(timepoint5) ; remove(timepoint6) ; remove(timepoint7) ; remove(timepoint8)

mutations <- mutations[order(mutations$POS),] ; rownames(mutations) <- NULL

mutations$SEQ_PLATFORM <- fasta_df$seq_platform[match(mutations$DATE, fasta_df$Date)]
E484A <- mutations[mutations$POS==23013,] ; E484A$FREQ <- as.numeric(E484A$FREQ) ; E484A$DEPTH <- as.numeric(E484A$DEPTH)
E484T <- mutations[mutations$POS==23012,] ; E484T$FREQ <- as.numeric(E484T$FREQ) ; E484T$DEPTH <- as.numeric(E484T$DEPTH)
write_csv(mutations, "data/all_timepoints_raw_reads/alltimepoints_mutations_20211117.csv",
          quote=F)
# mutations <- read.csv("/Volumes/working_ssd/e484t_manuscript/data/all_timepoints_raw_reads/alltimepoints_mutations_20211117.csv")
# mutations$DATE <- as.Date(mutations$DATE, format = "%m/%d/%y")



### FILTERING ####
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

write_csv(mutations, "data/all_timepoints_raw_reads/alltimepoints_filtered_mutations_20211117.csv",
          quote=F)


### CREATING AUXILIARY DATASETS ####
full_series <- mutations # the loop below makes a data frame of only mutations that are present for all 8 time points
for (i in unique(full_series$REF_POS_ALT)){
  mut_sub <- full_series[full_series$REF_POS_ALT==i,]
  
  if (nrow(mut_sub)<8){
    full_series <- full_series[!full_series$REF_POS_ALT==i,]
  } else {
    next
  }
  rownames(full_series) <- NULL
}

increases <- mutations # this loop creates a dataset full of mutations that increase in frequency through time
for (i in unique(increases$REF_POS_ALT)){
  mut_sub <- increases[increases$REF_POS_ALT==i,]
  
  if (mut_sub[mut_sub$DATE==max(mut_sub$DATE), "FREQ"] == 
      mut_sub[mut_sub$DATE==min(mut_sub$DATE), "FREQ"]) {
    increases <- increases[!increases$REF_POS_ALT==i,]
    
  } else if (mut_sub[mut_sub$DATE==max(mut_sub$DATE), "FREQ"] < 
             mut_sub[mut_sub$DATE==min(mut_sub$DATE), "FREQ"]){
    increases <- increases[!increases$REF_POS_ALT==i,]
    
  } else {
    next
  }
  rownames(increases) <- NULL
}

decreases <- mutations # this loop creates a dataset full of mutations that decrease in frequency through time
for (i in unique(decreases$REF_POS_ALT)){
  mut_sub <- decreases[decreases$REF_POS_ALT==i,]
  
  if (mut_sub[mut_sub$DATE==max(mut_sub$DATE), "FREQ"] == 
      mut_sub[mut_sub$DATE==min(mut_sub$DATE), "FREQ"]) {
    decreases <- decreases[!decreases$REF_POS_ALT==i,]
    
  } else if (mut_sub[mut_sub$DATE==max(mut_sub$DATE), "FREQ"] > 
             mut_sub[mut_sub$DATE==min(mut_sub$DATE), "FREQ"]){
    decreases <- decreases[!decreases$REF_POS_ALT==i,]
    
  } else {
    next
  }
  rownames(decreases) <- NULL
}

E484A_prefilter <- E484A
E484A <- mutations[mutations$POS==23013,]
E484T_prefilter <- E484T
E484T <- mutations[mutations$POS==23012,]

big_change <- mutations
for (i in unique(big_change$REF_POS_ALT)){
  mut_sub <- big_change[big_change$REF_POS_ALT==i,]
  
  if (max(mut_sub$FREQ)<0.4) {
    big_change <- big_change[!big_change$REF_POS_ALT==i,]
    
  } else {
    next
  }
  rownames(big_change) <- NULL
}



### FREQUENCY CONFIDENCE INTERVALS
mutations$UPPER <- NA
mutations$LOWER <- NA

for (i in 1:nrow(mutations)){
  p <- mutations[i, "FREQ"]
  n <- mutations[i, "DEPTH"]
  mutations[i, "UPPER"] <- p + (1.96 * sqrt( (p * (1-p))/n ))
  mutations[i, "LOWER"] <- p - (1.96 * sqrt( (p * (1-p))/n ))
}

E484A$UPPER <- NA
E484A$LOWER <- NA
for (i in 1:nrow(E484A)){
  p <- E484A[i, "FREQ"]
  n <- E484A[i, "DEPTH"]
  E484A[i, "UPPER"] <- p + (1.96 * sqrt( (p * (1-p))/n ))
  E484A[i, "LOWER"] <- p - (1.96 * sqrt( (p * (1-p))/n ))
}

E484T$UPPER <- NA
E484T$LOWER <- NA
for (i in 1:nrow(E484T)){
  p <- E484T[i, "FREQ"]
  n <- E484T[i, "DEPTH"]
  E484T[i, "UPPER"] <- p + (1.96 * sqrt( (p * (1-p))/n ))
  E484T[i, "LOWER"] <- p - (1.96 * sqrt( (p * (1-p))/n ))
}



### PLOTTING LOOPS ####
# setting the plot frame and grid
plot(mutations$DATE, mutations$FREQ,
     xlim = c(min(mutations$DATE), max(mutations$DATE)),
     xlab = "Date of Infection (Dec 2020 - Sept. 2021)", 
     ylim = c(0,1), ylab = "SNP Frequency", 
     frame.plot = F, cex.axis = 0.8, cex.lab = 0.85, las = 1, xaxt="n",
     pch = 20, col = "darkgray", cex = 0.8, type = 'n')
grid()

# adding axis for months
months_axis <- seq(min(mutations$DATE), max(mutations$DATE), by = "month")
axis(side = 1, at = months_axis,
     labels = format(months_axis, "%b"), cex.axis = 0.8)

# plotting mutation points and lines
for (i in unique(mutations$REF_POS_ALT)){
  mut_sub <- mutations[mutations$REF_POS_ALT==i,]
  
  lines(mut_sub$DATE, mut_sub$FREQ, col="gray")
  points(mut_sub$DATE, mut_sub$FREQ, pch = 20, col = "darkgray", cex = 0.8)
}

# plotting mutation confidence intervals
for (i in 1:nrow(mutations)){
  mut_sub <- mutations[i,]
  lines(x = rep(mut_sub$DATE, times = 2),
        y = c(mut_sub$UPPER, mut_sub$LOWER),
        col="darkgray")
  lines(x=c(mut_sub$DATE-1, mut_sub$DATE+1),
        y=rep(mut_sub$UPPER, times = 2), col="darkgray")
  lines(x=c(mut_sub$DATE-1, mut_sub$DATE+1),
        y=rep(mut_sub$LOWER, times = 2), col="darkgray")
}

# plotting E484A with confidence intervals
lines(E484A$DATE, E484A$FREQ, col="orange", lwd = 4)
points(E484A$DATE, E484A$FREQ, pch = 20, col = "orange", cex = 2)

for (i in 1:nrow(E484A)){
  mut_sub <- E484A[i,]
  lines(x = rep(mut_sub$DATE, times = 2),
        y = c(mut_sub$UPPER, mut_sub$LOWER),
        col="orange", lwd = 2)
  lines(x=c(mut_sub$DATE-2, mut_sub$DATE+2),
        y=rep(mut_sub$UPPER, times = 2), col="orange", lwd = 2)
  lines(x=c(mut_sub$DATE-2, mut_sub$DATE+2),
        y=rep(mut_sub$LOWER, times = 2), col="orange", lwd = 2)
}

# plotting E484T with confidence intervals
lines(E484T$DATE, E484T$FREQ, col="red", lwd = 4)
points(E484T$DATE, E484T$FREQ, pch = 20, col = "red", cex = 2.5)

for (i in 1:nrow(E484T)){
  mut_sub <- E484T[i,]
  lines(x = rep(mut_sub$DATE, times = 2),
        y = c(mut_sub$UPPER, mut_sub$LOWER),
        col="red", lwd = 2)
  lines(x=c(mut_sub$DATE-2, mut_sub$DATE+2),
        y=rep(mut_sub$UPPER, times = 2), col="red", lwd = 2)
  lines(x=c(mut_sub$DATE-2, mut_sub$DATE+2),
        y=rep(mut_sub$LOWER, times = 2), col="red", lwd = 2)
}

# legend(as.Date("2021-05-15"), 1.05, legend=c("A-23013-C", "G-23012-A", "Present at all time points"),
#        lwd = c(4, 3, 1.5),
#        col = c("orange", "red", "black"),
#        bty = "n", cex = 0.8, bg="transparent", ncol = 3, xpd = T,
#        xjust = 0.5, yjust = 0)

legend(as.Date("2021-05-15"), 1.05, legend=c("A-23013-C", "G-23012-A"),
       lwd = c(4, 3),
       col = c("orange", "red"),
       bty = "n", cex = 0.8, bg="transparent", ncol = 2, xpd = T,
       xjust = 0.5, yjust = 0)
