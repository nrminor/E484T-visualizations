#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



### PLOTTING INTRAHOST SINGLE-NUCLEOTIDE VARIANT (iSNV) FREQUENCES
# UPDATED: 04-Jan-2022 by Nicholas R. Minor
# ----------------------------------------------------------------- #

# This script will do the following:

# 1) Import VCFs called from GISAID subsample.
# 2) Import aligned FASTA sequences from GISAID subsample, plus their metadata.
# 3) Polish these three files.
# 4) Count single-nucleotide variations from Wuhan-1 (GenBank: MN908947.3)
# present at each time point, as well as in each GISAID sample.
# 5) Plot these mutations through time.

# Please direct all questions about the reproducibility of this script to:
# nrminor@wisc.edu

# ----------------------------------------------------------------- #



### PREPARING THE ENVIRONMENT ####
library(Biostrings)
library(tidyverse)
library(RColorBrewer)
data_filepath = args[1]  ### OR INSERT YOUR FILE PATH HERE ###
setwd(data_filepath)



### IMPORTING VCF & METADATA ####
gisaid_vcf <- read.delim("data/b12_enriched_global/b12_enriched_global_subsampled_aligned.vcf", skip = 13)
gisaid_vcf <- gisaid_vcf[grepl("VARSEQ",gisaid_vcf$INFO, fixed=T),] ; rownames(gisaid_vcf) <- NULL
gisaid_meta <- read.delim("data/b12_enriched_global/b12_enriched_global_subsampled_metadata.tsv")
strain_dates <- gisaid_meta[,c(1,4)]
strain_dates$date <- as.Date(strain_dates$date, format = "%Y-%m-%d")
strain_dates$day_of_infection <- as.numeric(NA)
for (i in 1:nrow(strain_dates)){
  strain_dates$day_of_infection[i] <- as.numeric(strain_dates$date[i] - as.Date("2020-09-05"))
  
}



### PREPARING THE GISAID SUBSAMPLE FASTA ####
gisaid_fasta <- readDNAStringSet("data/b12_enriched_global/b12_enriched_global.aligned.fasta")
seq_names = names(gisaid_fasta)
sequence = paste(gisaid_fasta)
gisaid_fasta_df <- data.frame(seq_names, sequence)
remove(gisaid_fasta)



### PREPARING METADATA & MATCHING WITH FASTA ####
gisaid_fasta_df <- gisaid_fasta_df[match(strain_dates$strain, gisaid_fasta_df$seq_names),] ; rownames(gisaid_fasta_df) <- NULL
gisaid_meta <- gisaid_meta[match(strain_dates$strain, gisaid_meta$strain),] ; rownames(gisaid_meta) <- NULL

# which(is.na(gisaid_fasta_df$seq_names)) # -> to_remove
# which(is.na(gisaid_meta$strain)) # -> to_remove2
gisaid_fasta_df$seq_names[which(is.na(gisaid_fasta_df$seq_names))] <- gisaid_meta$strain[which(is.na(gisaid_fasta_df$seq_names))]
which(is.na(gisaid_fasta_df$seq_names))
which(is.na(gisaid_meta$strain))

gisaid_fasta_df$date <- as.Date(strain_dates$date, format = "%Y-%m-%d")
gisaid_meta$date <- as.Date(strain_dates$date, format = "%Y-%m-%d")

gisaid_fasta_df$day_of_infection <- strain_dates$day_of_infection
gisaid_meta$day_of_infection <- strain_dates$day_of_infection



### IDENTIFYING VARIANT NAMES, TYPES, & FREQUENCIES IN VCF ####
for (i in 1:nrow(gisaid_vcf)){
  info <- gisaid_vcf[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  gisaid_vcf[i, "SAMPLES"] <- info_split[grep("^VARSEQ", info_split)]
}
gisaid_vcf$SAMPLES <- str_replace_all(gisaid_vcf$SAMPLES, "VARSEQ=", "")

for (i in 1:nrow(gisaid_vcf)){
  info <- gisaid_vcf[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  gisaid_vcf[i, "MUTATION_TYPE"] <- info_split[grep("^TYPE", info_split)]
}
gisaid_vcf$MUTATION_TYPE <- str_replace_all(gisaid_vcf$MUTATION_TYPE, "TYPE=", "")

for (i in 1:nrow(gisaid_vcf)){
  info <- gisaid_vcf[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  gisaid_vcf[i, "FREQ"] <- info_split[grep("^VF", info_split)]
}
gisaid_vcf$FREQ <- str_replace_all(gisaid_vcf$FREQ, "VF=", "")



### REMOVING ANY INDELS FROM GISAID VCF, LEAVING ONLY SNVs ####
gisaid_vcf <- gisaid_vcf[!grepl("Deletion", gisaid_vcf$MUTATION_TYPE),]



### COUNTING MUTATIONS FOR EACH SAMPLE ####
gisaid_fasta_df$distance <- rep(0, nrow(gisaid_fasta_df))
for (i in 1:length(gisaid_vcf$SAMPLES)){
  sub <- unlist(strsplit(gisaid_vcf$SAMPLES[i], split = ","))
  for (j in sub){
    add <- gisaid_fasta_df[gisaid_fasta_df$seq_names==j,"distance"] + 1
    gisaid_fasta_df[gisaid_fasta_df$seq_names==j,"distance"] <- add
  }
} 



### IDENTIFYING PANGO LINEAGES IN GISAID SUBSAMPLE ####
lineages <- as.data.frame(table(gisaid_meta$Pango.lineage))
colnames(lineages) <- c("lineage", "count")
lineages <- lineages[order(lineages$count, decreasing = T),] ; rownames(lineages) <- NULL
total_count <- sum(lineages$count) ; total_count == nrow(gisaid_fasta_df)
sum(lineages[1:7,"count"])/total_count # try to get at least 50% of the samples colored by lineage

palette <- read.delim("data/SupplementalFigure2_palette.txt", header = FALSE, sep = "\t")
palette$color <- c(rep(1, 5),
                   rep(2, 5),
                   rep(3, 5),
                   rep(4, 5))

# palette <- rev(c("#090c08","#757083","#8a95a5","#b9c6ae","#d6e5e3","#9fd8cb","#b4cded")) # palette from https://coolors.co/090c08-757083-8a95a5-b9c6ae-d6e5e3-9fd8cb-b4cded

lineages$color <- ""
lineages$color[1:4] <- palette[1:4,1] # brewer.pal(6, "Reds")[2:6]
lineages$color[5:nrow(lineages)] <- palette[5,1] # rep(brewer.pal(6, "Reds")[1], times = length(lineages$color[6:nrow(lineages)]))
lineages$label <- ""
lineages$label[1:4] <- as.character(lineages$lineage[1:4])
lineages$label[5:nrow(lineages)] <- rep("other", times = length(lineages$label[5:nrow(lineages)]))



### MOVING LINEAGE DATA INTO FASTA DF ####
gisaid_fasta_df$lineage <- gisaid_meta$Pango.lineage
gisaid_fasta_df$color <- ""
gisaid_fasta_df$label <- ""
for (i in 1:nrow(gisaid_fasta_df)){
  gisaid_fasta_df[i,"color"] <- lineages[lineages$lineage==gisaid_fasta_df$lineage[i],"color"]
  gisaid_fasta_df[i,"label"] <- lineages[lineages$lineage==gisaid_fasta_df$lineage[i],"label"]
}


### IMPORTING PATIENT DATA ####
patient_fasta <- readDNAStringSet("data/alltimepoints_20211104.fasta")
seq_names <- names(patient_fasta)
sequence <- paste(patient_fasta)
fasta_df <- data.frame(seq_names, sequence)
fasta_df$seq_names <- str_replace_all(fasta_df$seq_names, fixed(" | "), ",")
fasta_df$sequence <- str_replace_all(fasta_df$sequence, fixed("-"), "N")
fasta_df <- separate(data = fasta_df, col = 1,
                     into = c("Sample_ID", "Date"),
                     sep = ",")
fasta_df$Date <- as.Date(fasta_df$Date, "%Y %b %d")
fasta_df$day_of_infection <- c(113,124,131,159,198,297,333,388)



### COUNTING MUTATIONS #####
# reading in and reducing to usable mutations (those that contain "VARSEQ")
crude_vcf <- read.delim("data/alltimepoints_consensus_variants_20211104.vcf", skip = 20)
crude_vcf <- crude_vcf[grepl("VARSEQ",crude_vcf$INFO, fixed=T),] ; rownames(crude_vcf) <- NULL

# creating new column for dates with only a sample ID in it (for the time being)
for (i in 1:nrow(crude_vcf)){
  info <- crude_vcf[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  crude_vcf[i, "SAMPLE_ID"] <- info_split[grep("^VARSEQ", info_split)]
}
crude_vcf$SAMPLE_ID <- str_replace_all(crude_vcf$SAMPLE_ID, "VARSEQ=", "")
crude_vcf$SAMPLE_ID <- str_replace(crude_vcf$SAMPLE_ID, "-consensus", "")
crude_vcf$SAMPLE_ID <- str_replace(crude_vcf$SAMPLE_ID, "-DBaker", "")

# creating new column with only mutation type
for (i in 1:nrow(crude_vcf)){
  info <- crude_vcf[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  crude_vcf[i, "MUTATION_TYPE"] <- info_split[grep("^TYPE", info_split)]
}
crude_vcf$MUTATION_TYPE <- str_replace_all(crude_vcf$MUTATION_TYPE, "TYPE=", "")

# summing up mutations per sample
fasta_df$distance <- 0
for (i in 1:length(crude_vcf$SAMPLE_ID)){
  sub <- unlist(strsplit(crude_vcf$SAMPLE_ID[i], split = ","))
  for (j in sub){
    add <- fasta_df[fasta_df$Sample_ID==j,"distance"] + 1
    fasta_df[fasta_df$Sample_ID==j,"distance"] <- add
  }
}



### PLOTTING ####
pdf(file = "visuals/VOC_plot.pdf", 
    width = 9, height = 5)
plot(fasta_df$day_of_infection, fasta_df$distance, 
     xlim = c(min(gisaid_fasta_df$day_of_infection), max(gisaid_fasta_df$day_of_infection)), ylim = c(0,60),
     xlab = "Days Post-Diagnosis", ylab = "Genetic Distance from Wuhan-1", 
     frame.plot = F, cex.axis = 0.8, cex.lab = 0.85, las = 1, pch = 20,
     cex = 3, col = "#4B7395", type = "n")
grid()
# 
# months_axis <- seq(min(gisaid_fasta_df$date), max(gisaid_fasta_df$date), by = "month")
# axis(side = 1, at = months_axis,
#      labels = format(months_axis, "%b"), cex.axis = 0.8)

points(gisaid_fasta_df$day_of_infection, gisaid_fasta_df$distance,
       pch = 20, cex = 1, col = gisaid_fasta_df$color)

text(198, 55, labels = "Bamlanivimab\nTreatment", cex = 0.8, bty = "l")
segments(198, y0 = -2, y1 = 52, col = "red", lty = 2)
segments(198, y0 = 58, y1 = 61, col = "red", lty = 2)

# abline(VOC_lm, col = rgb(165/255,15/255,21/255,3/4), lwd = 4)

points(fasta_df$day_of_infection, fasta_df$distance,
       pch = 20, cex = 2, col = palette[11,1])
# abline(patient_lm, col = rgb(8/255,81/255,156/255,3/4), lwd = 4)
legend("topleft", legend = c(lineages$label[1:5], "patient"),
       col = c(lineages$color[1:5], palette[11,1]), bty="n",
       pch = 16, ncol = 2, xpd = T, xjust = 0.5, cex = 1)
dev.off()

