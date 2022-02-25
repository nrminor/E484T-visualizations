#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



### PLOTTING CONSENSUS MUTATIONS THROUGH TIME
# UPDATED: 24-Feb-2022 by Nicholas R. Minor
# ----------------------------------------------------------------- #

# This script will identify and plot single-nucleotide differences between each
# prolonged infection time point's consensus sequence and Wuhan-1 (GenBank:
# MN908947.3). It will also plot the genes in the SARS-CoV-2 genome, amplicon
# dropouts, and the globally unique Spike E484T substitution.

# Please direct all questions about the reproducibility of this script to:
# nrminor@wisc.edu

# ----------------------------------------------------------------- #



### SETUP ####
### ---- #
library(Biostrings)
library(tidyverse)
library(grid)
library(gridExtra)
library(gridGraphics)
data_filepath = args[1]  ### OR INSERT YOUR FILE PATH HERE ###
setwd(data_filepath)



### FASTA PREP ####
### --------- #
patient_fasta <- readDNAStringSet("data/alltimepoints_20220222.fasta")
seq_names <- names(patient_fasta)
fasta_df <- data.frame(seq_names) ; remove(patient_fasta)
fasta_df$seq_names <- str_replace_all(fasta_df$seq_names, fixed(" | "), ",")
fasta_df <- separate(data = fasta_df, col = 1,
                     into = c("Sample_ID", "Date"),
                     sep = ",")
fasta_df$Date <- as.Date(fasta_df$Date, "%Y %b %d")
fasta_df$day_of_infection <- c(113,124,131,159,198,297,333,388,415,417,482,486)
fasta_df$seq_platform <- c("ONT", "ONT", "ONT", "ONT", "ONT", "ONT", "Illumina", "ONT", "ONT", "ONT", "ONT", "ONT")



### VCF PREP ####
### ------- #
timepoint1 <- read.delim("data/USA_WI-WSLH-202168_2020_consensus_variant_table.tsv")
timepoint1 <- timepoint1[timepoint1$ALT!="N",]
timepoint1 <- timepoint1[,c("POS", "REF", "ALT", "GFF_FEATURE", "REF_CODON", 
                            "REF_AA", "ALT_CODON", "ALT_AA")]
timepoint1$DATE <- as.Date(fasta_df$Date[1])
timepoint1$DAY <- fasta_df$day_of_infection[1]

timepoint2 <- read.delim("data/USA_WI-UW-2731_2021_consensus_variant_table.tsv")
timepoint2 <- timepoint2[timepoint2$ALT!="N",]
timepoint2 <- timepoint2[,c("POS", "REF", "ALT", "GFF_FEATURE", "REF_CODON", 
                            "REF_AA", "ALT_CODON", "ALT_AA")]
timepoint2$DATE <- as.Date(fasta_df$Date[2])
timepoint2$DAY <- fasta_df$day_of_infection[2]

timepoint3 <- read.delim("data/USA_WI-UW-2731-T2_2021_consensus_variant_table.tsv")
timepoint3 <- timepoint3[timepoint3$ALT!="N",]
timepoint3 <- timepoint3[,c("POS", "REF", "ALT", "GFF_FEATURE", "REF_CODON", 
                            "REF_AA", "ALT_CODON", "ALT_AA")]
timepoint3$DATE <- as.Date(fasta_df$Date[3])
timepoint3$DAY <- fasta_df$day_of_infection[3]

timepoint4 <- read.delim("data/USA_WI-UW-2731-T3_2021_consensus_variant_table.tsv")
timepoint4 <- timepoint4[timepoint4$ALT!="N",]
timepoint4 <- timepoint4[,c("POS", "REF", "ALT", "GFF_FEATURE", "REF_CODON", 
                            "REF_AA", "ALT_CODON", "ALT_AA")]
timepoint4$DATE <- as.Date(fasta_df$Date[4])
timepoint4$DAY <- fasta_df$day_of_infection[4]

timepoint5 <- read.delim("data/USA_WI-UW-2731-T4_2021_consensus_variant_table.tsv")
timepoint5 <- timepoint5[timepoint5$ALT!="N",]
timepoint5 <- timepoint5[,c("POS", "REF", "ALT", "GFF_FEATURE", "REF_CODON", 
                            "REF_AA", "ALT_CODON", "ALT_AA")]
timepoint5$DATE <- as.Date(fasta_df$Date[5])
timepoint5$DAY <- fasta_df$day_of_infection[5]

timepoint6 <- read.delim("data/USA_WI-UW-5350_2021_consensus_variant_table.tsv")
timepoint6 <- timepoint6[timepoint6$ALT!="N",]
timepoint6 <- timepoint6[,c("POS", "REF", "ALT", "GFF_FEATURE", "REF_CODON", 
                            "REF_AA", "ALT_CODON", "ALT_AA")]
timepoint6$DATE <- as.Date(fasta_df$Date[6])
timepoint6$DAY <- fasta_df$day_of_infection[6]

timepoint7 <- read.delim("data/USA_WI-WSLH-217727_2021_consensus_variant_table.tsv")
timepoint7 <- timepoint7[timepoint7$ALT!="N",]
timepoint7 <- timepoint7[,c("POS", "REF", "ALT", "GFF_FEATURE", "REF_CODON", 
                            "REF_AA", "ALT_CODON", "ALT_AA")]
timepoint7$DATE <- as.Date(fasta_df$Date[7])
timepoint7$DAY <- fasta_df$day_of_infection[7]

timepoint8 <- read.delim("data/USA_WI-UW-PI08_2021_consensus_variant_table.tsv")
timepoint8 <- timepoint8[timepoint8$ALT!="N",]
timepoint8 <- timepoint8[,c("POS", "REF", "ALT", "GFF_FEATURE", "REF_CODON", 
                            "REF_AA", "ALT_CODON", "ALT_AA")]
timepoint8$DATE <- as.Date(fasta_df$Date[8])
timepoint8$DAY <- fasta_df$day_of_infection[8]

timepoint9 <- read.delim("data/USA_WI-UW-PI09_2021_consensus_variant_table.tsv")
timepoint9 <- timepoint9[timepoint9$ALT!="N",]
timepoint9 <- timepoint9[,c("POS", "REF", "ALT", "GFF_FEATURE", "REF_CODON", 
                            "REF_AA", "ALT_CODON", "ALT_AA")]
timepoint9$DATE <- as.Date(fasta_df$Date[9])
timepoint9$DAY <- fasta_df$day_of_infection[9]

timepoint10 <- read.delim("data/USA_WI-UW-PI10_2021_consensus_variant_table.tsv")
timepoint10 <- timepoint10[timepoint10$ALT!="N",]
timepoint10 <- timepoint10[,c("POS", "REF", "ALT", "GFF_FEATURE", "REF_CODON", 
                            "REF_AA", "ALT_CODON", "ALT_AA")]
timepoint10$DATE <- as.Date(fasta_df$Date[10])
timepoint10$DAY <- fasta_df$day_of_infection[10]

timepoint11 <- read.delim("data/USA_WI-UW-PI11_2021_consensus_variant_table.tsv")
timepoint11 <- timepoint11[timepoint11$ALT!="N",]
timepoint11 <- timepoint11[,c("POS", "REF", "ALT", "GFF_FEATURE", "REF_CODON", 
                            "REF_AA", "ALT_CODON", "ALT_AA")]
timepoint11$DATE <- as.Date(fasta_df$Date[11])
timepoint11$DAY <- fasta_df$day_of_infection[11]

timepoint12 <- read.delim("data/USA_WI-UW-PI12_2022_consensus_variant_table.tsv")
timepoint12 <- timepoint12[timepoint12$ALT!="N",]
timepoint12 <- timepoint12[,c("POS", "REF", "ALT", "GFF_FEATURE", "REF_CODON", 
                            "REF_AA", "ALT_CODON", "ALT_AA")]
timepoint12$DATE <- as.Date(fasta_df$Date[12])
timepoint12$DAY <- fasta_df$day_of_infection[12]

variants <- rbind(timepoint1, timepoint2, timepoint3, timepoint4, timepoint5, 
                  timepoint6, timepoint7, timepoint8, timepoint9, timepoint10,
                  timepoint11, timepoint12)
# sanity check before running the following line of code: 
# in the environment, do all the individual timepoint variant tables show an
# increase in the number of mutations through time?
remove(timepoint1, timepoint2, timepoint3, timepoint4, timepoint5, 
       timepoint6, timepoint7, timepoint8, timepoint9, timepoint10,
       timepoint11, timepoint12)

variants$REF_POS_ALT <- paste(variants$REF, variants$POS, variants$ALT, sep = "-")
variants <- variants[str_length(variants$ALT)==1,]



### GENE INFO AND LABELS ####
### ------------------- #
SARS_genes <- data.frame("name" = c("ORF1a", "ORF1b", "S", NA, NA, NA, NA, NA, NA, NA, "N"),
                         "start" = c(266, 13469, 21563, 25393, 26245, 26523, 27202, 27394, 27894, 28284, 28578),
                         "stop" = c(13468, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 28259, 28577, 29533),
                         "col" = c("#648FFF", "#FFB000", "#DC267F", "#FE6100", "#648FFF", "#FFB000", "#DC267F", 
                                   "#FE6100", "#648FFF", "#FFB000", "#DC267F"),
                         "left_out_names" = c(NA, NA, NA, "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "ORF9b", NA))



### LOW COVERAGE REGIONS ####
### -------------------- #
lowcov_filepath <- paste0(getwd(), "/data/low_coverage_regions", sep = "")
lowcov_docs <- list.files(lowcov_filepath)
for (i in 1:length(lowcov_docs)){
  name <- lowcov_docs[i]
  name <- str_replace(name, ".bed", "")
  name <- paste("lowcov", name, sep = "_")
  
  tmp <- read.delim(paste(lowcov_filepath, lowcov_docs[i], sep = "/"), 
                    header = F, sep = "\t")[-1]
  colnames(tmp) <- c("start", "stop", "depth")
  tmp$Date <- fasta_df$Date[i]
  tmp$day_of_infection <- fasta_df$day_of_infection[i]
  assign(name, tmp)
  
  remove(tmp) ; remove(name)
}

lowcov <- rbind(lowcov_timepoint01_day113_ONT, lowcov_timepoint02_day124_ONT,
                lowcov_timepoint03_day131_ONT, lowcov_timepoint04_day159_ONT,
                lowcov_timepoint05_day198_ONT, lowcov_timepoint06_day297_ONT,
                lowcov_timepoint07_day333_Illumina, lowcov_timepoint08_day388_ONT,
                lowcov_timepoint09_day415_ONT, lowcov_timepoint10_day417_ONT,
                lowcov_timepoint11_day482_ONT, lowcov_timepoint12_day488_ONT)
remove(lowcov_timepoint01_day113_ONT, lowcov_timepoint02_day124_ONT,
       lowcov_timepoint03_day131_ONT, lowcov_timepoint04_day159_ONT,
       lowcov_timepoint05_day198_ONT, lowcov_timepoint06_day297_ONT,
       lowcov_timepoint07_day333_Illumina, lowcov_timepoint08_day388_ONT,
       lowcov_timepoint09_day415_ONT, lowcov_timepoint10_day417_ONT,
       lowcov_timepoint11_day482_ONT, lowcov_timepoint12_day488_ONT)

lowcov$new_block <- NA
for (i in 2:nrow(lowcov)){
  if (lowcov$start[i]==lowcov$stop[i-1]){
    lowcov$new_block[i] <- FALSE
  } else {
    lowcov$new_block[i] <- TRUE
  }
}
lowcov$new_block[1] <- TRUE

new_blocks <- as.numeric(rownames(lowcov[lowcov$new_block==T,]))
for (i in 1:length(new_blocks)){
  lowcov[new_blocks[i], "stop"] <- lowcov[new_blocks[i+1]-1, "stop"]
}
lowcov <- lowcov[lowcov$new_block==TRUE,] ; row.names(lowcov) <- NULL
lowcov <- lowcov[,-3] ; lowcov <- lowcov[,-ncol(lowcov)]





### VERTICAL PLOT ####
### ------------ #
pdf("visuals/mutations_across_genome_vertical.pdf",
    width = 8, height = 5)

  plot(1,1, ylim = c(0,30000), xlim=c(-15, 450), type = "n",
      ylab = "SARS-CoV-2 Genome Position", xlab = "Days Post-Diagnosis", 
      frame.plot = F, cex.axis = 1, yaxt = "n")
  axis(side = 2, at = c(0,5000,10000,15000,20000,25000,30000), cex.axis = 0.8,
       col = "white", col.ticks = "black", las = 1)
  
  # Plotting SARS-CoV-2 Genes
  for (i in 1:nrow(SARS_genes)){
    polygon(y = c(SARS_genes[i,"start"], SARS_genes[i,"start"], 
                  SARS_genes[i,"stop"], SARS_genes[i,"stop"]),  # y-coordinates of gene 
            x = c(-30, 0, 0, -30),    # x-Coordinates of polygon
            col = SARS_genes[i,"col"], border = F, xpd = T)
    text(y = SARS_genes[i,"start"] + (SARS_genes[i,"stop"]-SARS_genes[i,"start"])/2, 
         x = -15, labels = SARS_genes[i, "name"], srt = 90, col = "black", cex = 1)
  }
  for (i in 1:nrow(SARS_genes)){
    rgb <- c(col2rgb(SARS_genes[i,"col"], alpha = T))
    rgb_faded <- round(rgb/c(1,1,1,4))
    col_tmp <- rgb(rgb_faded[1], rgb_faded[2], rgb_faded[3], rgb_faded[4], 
                   max = 255)
    polygon(y = c(SARS_genes[i,"start"], SARS_genes[i,"start"], 
                  SARS_genes[i,"stop"], SARS_genes[i,"stop"]),  # coordinates of gene 
            x = c(0, 417, 417, 0),    # Y-Coordinates of polygon
            col = col_tmp, border = F)
  }

  # Plotting regions with coverage lower than 20 reads
  for (i in 1:nrow(lowcov)){
    polygon(y = c(lowcov[i,"start"], lowcov[i,"start"],
                  lowcov[i,"stop"], lowcov[i,"stop"]),  # coordinates of low cov region
            x = c(lowcov[i, "day_of_infection"]-10, lowcov[i, "day_of_infection"]+10,
                  lowcov[i, "day_of_infection"]+10, lowcov[i, "day_of_infection"]-10),    # X-Coordinates of polygon
            col = "white", border = F)
  }
  
  # Plotting all mutations
  for (i in unique(variants$REF_POS_ALT)){
    plotting_sub <- variants[variants$REF_POS_ALT==i,]
    points(y = plotting_sub$POS, x = plotting_sub$DAY,
           type = "p", pch = 20, col = "#bcbcbc")
  }
  segments(y0 = -1200, y1 = 29800, x0=198, col = "red")
  
  # Plotting "unfixed" mutations
  fixed <- rep(NA, times = length(unique(variants$REF_POS_ALT)))
  for (i in 1:length(fixed)){
    fixed[i] <- nrow(variants[variants$REF_POS_ALT==unique(variants$REF_POS_ALT)[i],])==12
  }
  
  for (i in unique(variants$REF_POS_ALT)[!fixed]){
    plotting_sub <- variants[variants$REF_POS_ALT==i, ]
    
    points(y = plotting_sub$POS, x = plotting_sub$DAY, pch = 16, col = "black", cex = 1.2)
    
  }
  
  # Plotting S E484T specifically
  points(y = 23012, x = 297, 
         pch = 2, cex = 1.5, col = "red")
  points(y = 23012, x = 333, 
         pch = 2, cex = 1.5, col = "red")
  text(y = 23700, x = 316, labels = "E484T", col = "#7e7e7e", cex = 1)
  
  # Plotting the count of mutations at each timepoint
  for (i in unique(variants$DAY)){
    day_count <- nrow(variants[variants$DAY==i,])
    text(y = 30500, x = i, labels = day_count, 
         col = "#7e7e7e", cex = 0.5)
  }
  
  # legend(200, 28000, legend=c("ancestral", "derived during infection"), 
  #        pch = c(20, 16),
  #        col = c("#e5e5e5", "#7e7e7e", "red"),
  #        bty = "n", cex = 1, bg="transparent", ncol = 2, xpd = T, 
  #        xjust = 0.5, yjust = 0)
# }
dev.off()
