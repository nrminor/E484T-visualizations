#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



### PLOTTING CONSENSUS MUTATIONS THROUGH TIME
# UPDATED: 04-Jan-2022 by Nicholas R. Minor
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



### VCF PREP ####
### ------- #
crude_vcf <- read.delim("data/alltimepoints_consensus_variants_20211104.vcf", skip = 20)
crude_vcf <- crude_vcf[grepl("VARSEQ",crude_vcf$INFO, fixed=T),] ; rownames(crude_vcf) <- NULL

# creating a column with all dates of detection for each mutation
for (i in 1:nrow(crude_vcf)){
  info <- crude_vcf[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  crude_vcf[i, "DATES_DETECTED"] <- info_split[grep("^VARSEQ", info_split)]
}
crude_vcf$DATES_DETECTED <- str_replace_all(crude_vcf$DATES_DETECTED, "VARSEQ=", "")

# creating columns that separate out codon number, position in codon, gene, nucleotide change, amino acid change, and protein effect.
for (i in 1:nrow(crude_vcf)){
  info <- crude_vcf[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  
  if (grepl("CDSCN=", info, fixed = T)){
    crude_vcf[i, "CODON_NUMBER"] <- info_split[grep("^CDSCN", info_split)]
  } else {
    crude_vcf[i, "CODON_NUMBER"] <- NA
  }
  
  if (grepl("CDSPWC=", info, fixed = T)){
    crude_vcf[i, "CODON_POSITION"] <- info_split[grep("^CDSPWC", info_split)]
  } else {
    crude_vcf[i, "CODON_POSITION"] <- NA
  }
  
  if (grepl("CDS=", info, fixed = T)){
    crude_vcf[i, "GENE"] <- info_split[grep("^CDS=", info_split)]
  } else {
    crude_vcf[i, "GENE"] <- NA
  }
  
  if (grepl("CDNCHG=", info, fixed = T)){
    crude_vcf[i, "NUCLEOTIDE_CHANGE"] <- info_split[grep("^CDNCHG=", info_split)]
  } else {
    crude_vcf[i, "NUCLEOTIDE_CHANGE"] <- NA
  }
  
  if (grepl("AACHG=", info, fixed = T)){
    crude_vcf[i, "AA_CHANGE"] <- info_split[grep("^AACHG=", info_split)]
  } else {
    crude_vcf[i, "AA_CHANGE"] <- NA
  }
  
  if (T %in% grepl("^PE=", info_split)){
    crude_vcf[i, "PROTEIN_EFFECT"] <- info_split[grep("^PE=", info_split)]
  } else {
    crude_vcf[i, "PROTEIN_EFFECT"] <- NA
  }
  
}
crude_vcf$CODON_NUMBER <- str_replace_all(crude_vcf$CODON_NUMBER, "CDSCN=", "")
crude_vcf$CODON_POSITION <- str_replace_all(crude_vcf$CODON_POSITION, "CDSPWC=", "")
crude_vcf$GENE <- str_replace_all(crude_vcf$GENE, "CDS=", "")
crude_vcf$GENE <- str_replace_all(crude_vcf$GENE, "CDS", "")
crude_vcf$NUCLEOTIDE_CHANGE <- str_replace_all(crude_vcf$NUCLEOTIDE_CHANGE, "CDNCHG=", "")
crude_vcf$AA_CHANGE <- str_replace_all(crude_vcf$AA_CHANGE, "AACHG=", "")
crude_vcf$PROTEIN_EFFECT <- str_replace_all(crude_vcf$PROTEIN, "PE=", "")
for (i in 1:length(crude_vcf$PROTEIN_EFFECT)){
  
  if (is.na(crude_vcf$PROTEIN_EFFECT[i])){
    crude_vcf$PROTEIN_EFFECT[i] <- "noncoding"
  } else if (crude_vcf$PROTEIN_EFFECT[i]=="None"){
    crude_vcf$PROTEIN_EFFECT[i] <- "synonymous"
  } else if (crude_vcf$PROTEIN_EFFECT[i]=="Substitution"){
    crude_vcf$PROTEIN_EFFECT[i] <- "nonsynonymous"
  }
  
}

for (i in 1:nrow(crude_vcf)){
  
  if (crude_vcf$PROTEIN_EFFECT[i]=="nonsynonymous"){
    gene <- crude_vcf$GENE[i]
    codon <- crude_vcf$CODON_NUMBER[i]
    aa_sub <- crude_vcf$AA_CHANGE[i]
    aa_sub <- str_replace(aa_sub, "->", codon)
    
    crude_vcf$FULL_SUB_NAME[i] <- paste(gene, aa_sub, sep = " ")
    
  } else {
    crude_vcf$FULL_SUB_NAME[i] <- NA
  }
  
}

# parsing out and sorting the dates
for (i in 1:nrow(crude_vcf)){
  varseq <- crude_vcf[i,"DATES_DETECTED"]
  samples <- unlist(strsplit(varseq, split = ","))
  samples_ordered <- samples[match(fasta_df$Sample_ID, samples)]
  samples_ordered <- as.character(na.omit(samples_ordered))
  dates_ordered <- samples_ordered
  ref <- fasta_df[match(dates_ordered, fasta_df$Sample_ID),1:2]
  for (j in 1:length(dates_ordered)){ 
    dates_ordered[j] <- as.character(ref[j, "Date"])
  }
  if (length(dates_ordered)>1){
    dates_ordered <- paste0(dates_ordered, collapse = ",")
  }
  crude_vcf[i,"DATES_DETECTED"] <- as.character(dates_ordered)
}

# creating a column with all *days* (not dates) of detection, from the patient's first day of infection
for (i in 1:nrow(crude_vcf)){
  info <- crude_vcf[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  crude_vcf[i, "DAYS_DETECTED"] <- info_split[grep("^VARSEQ", info_split)]
}
crude_vcf$DAYS_DETECTED <- str_replace_all(crude_vcf$DAYS_DETECTED, "VARSEQ=", "")

# parsing out and sorting the days
for (i in 1:nrow(crude_vcf)){
  varseq <- crude_vcf[i,"DAYS_DETECTED"]
  samples <- unlist(strsplit(varseq, split = ","))
  samples_ordered <- samples[match(fasta_df$Sample_ID, samples)]
  samples_ordered <- as.character(na.omit(samples_ordered))
  days_ordered <- samples_ordered
  ref <- fasta_df[match(days_ordered, fasta_df$Sample_ID),c("Sample_ID","day_of_infection")]
  for (j in 1:length(days_ordered)){ 
    days_ordered[j] <- as.character(ref[j, "day_of_infection"])
  }
  if (length(days_ordered)>1){
    days_ordered <- paste0(days_ordered, collapse = ",")
  }
  crude_vcf[i,"DAYS_DETECTED"] <- as.character(days_ordered)
}

crude_vcf$EARLIEST_DATE <- NA
crude_vcf$EARLIEST_DAY <- NA

for (i in 1:nrow(crude_vcf)){
  sub <- crude_vcf$DATES_DETECTED[i]
  crude_vcf$EARLIEST_DATE[i] <- unlist(strsplit(sub, split = ","))[1]
}
crude_vcf$EARLIEST_DATE <- as.Date(crude_vcf$EARLIEST_DATE, format = "%Y-%m-%d")
for (i in 1:nrow(crude_vcf)){
  sub <- crude_vcf$DAYS_DETECTED[i]
  crude_vcf$EARLIEST_DAY[i] <- unlist(strsplit(sub, split = ","))[1]
}
crude_vcf$EARLIEST_DAY <- as.numeric(crude_vcf$EARLIEST_DAY)



crude_vcf <- crude_vcf[,c(1,2,4,5,11:ncol(crude_vcf))]
crude_vcf <- crude_vcf[,c(1:4, 6:12, 5, 14, 13, 15)]
colnames(crude_vcf)[1] <- "GENBANK_REFERENCE"
syn <- crude_vcf[crude_vcf$PROTEIN_EFFECT=="synonymous",] ; syn_count <- nrow(syn)
nonsyn <- crude_vcf[crude_vcf$PROTEIN_EFFECT=="nonsynonymous",] ; nonsyn_count <- nrow(nonsyn)

paste("There are", nonsyn_count, "nonsynonomous mutations that arose during this infection.", sep = " ")

aa_subs <- crude_vcf[,c(11,10,7,9,5,6,3,2,4,8,14,12,1)]
aa_subs$DAYS_DETECTED <- str_replace_all(aa_subs$DAYS_DETECTED, ",", ";")
write.csv(aa_subs, "data/consensus_mutations_codon_aminoacid.csv", 
          row.names = F)

pdf("visuals/mutation_spacing_by_gene.pdf")
par(mfrow=c(2,3))
for (i in unique(nonsyn$GENE)[is.na(unique(nonsyn$GENE))==F]){
  gene_sub <- nonsyn[nonsyn$GENE==i, ]
  mut_spacing <- c(0)
  
  for (j in 2:nrow(gene_sub)){
    dist <- gene_sub$POS[j] - gene_sub$POS[j-1]
    mut_spacing <- c(mut_spacing, dist)
  }
  
  hist(mut_spacing, main = i)
}
dev.off()

mutation_density_by_gene <- data.frame("GENE" = unique(crude_vcf$GENE)[is.na(unique(crude_vcf$GENE))==F], 
                                       "GENE_LENGTH" = NA, 
                                       "MUTATION_COUNT" = NA,
                                       "MUTATION_DENSITY" = NA,
                                       "NONSYN_COUNT" = NA,
                                       "NONSYN_DENSITY" = NA)
mutation_density_by_gene$GENE_LENGTH <- c((21555-266), (25384-21563), 
                                          (26220-25393), (26472-26245),
                                          (27759-27394), (28259-27894),
                                          (29533-28274))
for (i in mutation_density_by_gene$GENE){
  mutation_density_by_gene[mutation_density_by_gene$GENE==i, "MUTATION_COUNT"] <- nrow(aa_subs[aa_subs$GENE==i,])
  mutation_density_by_gene[mutation_density_by_gene$GENE==i, "NONSYN_COUNT"] <- nrow(nonsyn[nonsyn$GENE==i,])
}
for (i in 1:nrow(mutation_density_by_gene)){
  mutation_density_by_gene$MUTATION_DENSITY[i] <- mutation_density_by_gene$MUTATION_COUNT[i]/mutation_density_by_gene$GENE_LENGTH[i]
  mutation_density_by_gene$NONSYN_DENSITY[i] <- mutation_density_by_gene$NONSYN_COUNT[i]/mutation_density_by_gene$GENE_LENGTH[i]
}
write.csv(mutation_density_by_gene, "data/mutation_density_per_gene.csv", row.names=F)

density_time <- data.frame("GENE" = rep(unique(crude_vcf$GENE)[is.na(unique(crude_vcf$GENE))==F], times = 8), 
                           "DAY" = 0,
                           "GENE_LENGTH" = 0, 
                           "MUTATION_COUNT" = 0,
                           "MUTATION_DENSITY" = 0,
                           "NONSYN_COUNT" = 0,
                           "NONSYN_DENSITY" = 0)

density_time$GENE_LENGTH <- c((21555-266), (25384-21563), 
                                          (26220-25393), (26472-26245),
                                          (27759-27394), (28259-27894),
                                          (29533-28274))
density_time <- density_time[order(density_time$GENE),]
density_time$DAY <- fasta_df$day_of_infection

for (i in unique(density_time$GENE)){
  
  vcf <- crude_vcf[is.na(crude_vcf$GENE)==F, ]
  
  sub <- density_time[density_time$GENE==i,] ; rownames(sub) <- NULL
  vcf_sub <- vcf[vcf$GENE==i,] ; rownames(vcf_sub) <- NULL
  nonsyn_sub <- nonsyn[nonsyn$GENE==i,] ; rownames(nonsyn) <- NULL
  
  for (j in sub$DAY){
    
    for (k in 1:nrow(vcf_sub)){
      if (grepl(j, vcf_sub$DAYS_DETECTED[k])){
        density_time[density_time$GENE==i & density_time$DAY==j,"MUTATION_COUNT"] <- density_time[density_time$GENE==i & density_time$DAY==j,"MUTATION_COUNT"] + 1
      } else {
        next
      }
    }
    
    if (nrow(nonsyn_sub)>1){
      for (l in 1:nrow(nonsyn_sub)){
        if (grepl(j, nonsyn_sub$DAYS_DETECTED[l])){
          density_time[density_time$GENE==i & density_time$DAY==j,"NONSYN_COUNT"] <- density_time[density_time$GENE==i & density_time$DAY==j,"NONSYN_COUNT"] + 1
        } else {
          next
        }
      }
    } else {
      next
    }
  }
}
for (i in 1:nrow(density_time)){
  density_time$MUTATION_DENSITY[i] <- density_time$MUTATION_COUNT[i]/density_time$GENE_LENGTH[i]
  density_time$NONSYN_DENSITY[i] <- density_time$NONSYN_COUNT[i]/density_time$GENE_LENGTH[i]
}
write.csv(density_time, "data/mutation_density_per_gene_through_time.csv", row.names=F)

plotting_colors <- c("red", "blue", "green", "yellow", "orange", "purple", "brown")
pdf("visuals/mutations_density_through_time.pdf")
par(bty="n")
plot(1,1, xlim = c(100,400), ylim = c(0,max(density_time$MUTATION_DENSITY)),
     xlab = "Post-Diagnosis Day", ylab = "Mutation Density (Count/base pairs)",
     cex.lab=0.9)
for (i in 1:length(unique(density_time$GENE))){
  plotting_gene <- density_time[density_time$GENE==unique(density_time$GENE)[i],]
  lines(plotting_gene$DAY, plotting_gene$MUTATION_DENSITY, col = plotting_colors[i],
        lwd=3)
}
legend(160, 0.006, legend = unique(density_time$GENE), col = plotting_colors, 
       lty = rep(1, length(unique(density_time$GENE))), 
       lwd = rep(2, length(unique(density_time$GENE))), 
       cex = 0.6, bty = 'n',bg="transparent", xpd = T, ncol = 4)
dev.off()



### MUTATION LISTING ####
### ---------------- #
mut_list <- vector("list", length = nrow(crude_vcf))
names(mut_list) <- as.numeric(crude_vcf$POS)
for (i in 1:nrow(crude_vcf)){
  sub <- crude_vcf[i, "DAYS_DETECTED"]
  sub_split <- unlist(strsplit(sub, split = ","))
  days <- as.numeric(sub_split)
  mut_list[[i]] <- days
}



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
  name <- str_replace(name, "_lowcov.csv", "")
  name <- paste("lowcov", name, sep = "_")
  
  tmp <- read.csv(paste(lowcov_filepath, lowcov_docs[i], sep = "/"))
  assign(name, tmp)
  
  remove(tmp) ; remove(name)
}

lowcov_20201227_ONT$Date <- fasta_df$Date[1]
lowcov_20210107_ONT$Date <- fasta_df$Date[2]
lowcov_20210114_ONT$Date <- fasta_df$Date[3]
lowcov_20210211_ONT$Date <- fasta_df$Date[4]
lowcov_20210322_ONT$Date <- fasta_df$Date[5]
lowcov_20210629_ONT$Date <- fasta_df$Date[6]
lowcov_20210804_Illumina$Date <- fasta_df$Date[7]
lowcov_20210928_IonTorrent$Date <- fasta_df$Date[8]

lowcov_20201227_ONT$day_of_infection <- fasta_df$day_of_infection[1]
lowcov_20210107_ONT$day_of_infection <- fasta_df$day_of_infection[2]
lowcov_20210114_ONT$day_of_infection <- fasta_df$day_of_infection[3]
lowcov_20210211_ONT$day_of_infection <- fasta_df$day_of_infection[4]
lowcov_20210322_ONT$day_of_infection <- fasta_df$day_of_infection[5]
lowcov_20210629_ONT$day_of_infection <- fasta_df$day_of_infection[6]
lowcov_20210804_Illumina$day_of_infection <- fasta_df$day_of_infection[7]
lowcov_20210928_IonTorrent$day_of_infection <- fasta_df$day_of_infection[8]

lowcov <- rbind(lowcov_20201227_ONT, lowcov_20210107_ONT, lowcov_20210114_ONT, lowcov_20210211_ONT,
      lowcov_20210322_ONT, lowcov_20210629_ONT, lowcov_20210804_Illumina, lowcov_20210928_IonTorrent)
remove(lowcov_20201227_ONT, lowcov_20210107_ONT, lowcov_20210114_ONT, lowcov_20210211_ONT,
       lowcov_20210322_ONT, lowcov_20210629_ONT, lowcov_20210804_Illumina, lowcov_20210928_IonTorrent)



### VERTICAL PLOT ####
### ------------ #
pdf("visuals/mutations_across_genome_vertical.pdf",
    width = 8, height = 5)

  plot(1,1, ylim = c(0,30000), xlim=c(-15, 450), type = "n",
      ylab = "SARS-CoV-2 Genome Position", xlab = "Days Post-Diagnosis", 
      frame.plot = F, cex.axis = 1, yaxt = "n")
  axis(side = 2, at = c(0,5000,10000,15000,20000,25000,30000), cex.axis = 0.8,
       col = "white", col.ticks = "black", las = 1)
  
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

  for (i in 1:nrow(lowcov)){

    polygon(y = c(lowcov[i,"Minimum"], lowcov[i,"Minimum"],
                  lowcov[i,"Maximum"], lowcov[i,"Maximum"]),  # coordinates of low cov region
            x = c(lowcov[i, "day_of_infection"]-10, lowcov[i, "day_of_infection"]+10,
                  lowcov[i, "day_of_infection"]+10, lowcov[i, "day_of_infection"]-10),    # X-Coordinates of polygon
            col = "white", border = F)
    
  }
  
  for (i in 1:length(mut_list)){
    plotting_data <- data.frame("Day" = mut_list[[i]],
                                "Coordinate" = rep(as.numeric(names(mut_list)[i]), times = length(mut_list[[i]])))
    points(y = plotting_data$Coordinate, x = plotting_data$Day, type = "p", pch = 20, col = "#bcbcbc")
  }
  
  segments(y0 = -1200, y1 = 29800, x0=198, col = "red")
  
  fixed <- rep(NA, times = length(mut_list))
  for (i in 1:length(fixed)){
    fixed[i] <- length(mut_list[[i]])==8
  }
  
  for (i in crude_vcf[!fixed, "POS"]){
    days <- mut_list[as.numeric(names(mut_list))==i]
    
    if (length(days)==1){
      
      positions <- rep(as.numeric(names(days)), time = length(days[[1]]))
      days <- days[[1]]
      points(y = positions, x = days, pch = 16, col = "black", cex = 1.2)
      
    } else {
      
      for (j in 1:length(days)){
        
        positions <- rep(as.numeric(names(days[j])), time = length(days[j][[1]]))
        days_sub <- days[j][[1]]
        points(y = positions, x = days_sub, pch = 16, col = "black", cex = 1.2)
        
      }
      
    }
    
  }
  
  points(y = 23012, x = crude_vcf[crude_vcf$POS==23012,"EARLIEST_DAY"], 
         pch = 2, cex = 1.5, col = "red")
  points(y = 23012, x = 333, 
         pch = 2, cex = 1.5, col = "red")
  
  text(y = 23700, x = 316, labels = "E484T", col = "#7e7e7e", cex = 1)
  
  # legend(200, 28000, legend=c("ancestral", "derived during infection"), 
  #        pch = c(20, 16),
  #        col = c("#e5e5e5", "#7e7e7e", "red"),
  #        bty = "n", cex = 1, bg="transparent", ncol = 2, xpd = T, 
  #        xjust = 0.5, yjust = 0)
# }
dev.off()
