# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("Gviz", "RSamtools", "Biostrings", "GenVisR"))
# install.packages("RCircos")

library(Gviz)
library(RSamtools)
library(Biostrings)
library(tidyverse)
library(GenVisR)
library(RCircos)
library(vcfR)

# reading in and reformatting sequence data
patient_fasta <- readDNAStringSet("/Volumes/Nick Minor External HD 1/E484T/data/PT0001_alltimepoints_20210908.fasta")
seq_names = names(patient_fasta)
sequence = paste(patient_fasta)
fasta_df <- data.frame(seq_names, sequence)
fasta_df$seq_names <- str_replace_all(fasta_df$seq_names, fixed(" | "), ",")
fasta_df <- separate(data = fasta_df, col = 1,
                     into = c("Sample_ID", "Date"),
                     sep = ",")
fasta_df$Date <- as.Date(fasta_df$Date, "%Y %b %d")
fasta_df$day_of_infection <- c(113,124,131,159,198,297)

# reading in "crude" (headerless) variant data from all 6 single-sample VCFs
test1_variants = read.delim("/Volumes/Nick Minor External HD 1/E484T/data/USA_WI-WSLH-202168_2020.vcf", skip = 12)[,-10]
test2_variants = read.delim("/Volumes/Nick Minor External HD 1/E484T/data/USA_WI-UW-5350_2021.vcf", skip = 12)[,-10]
test3_variants = read.delim("/Volumes/Nick Minor External HD 1/E484T/data/USA_WI-UW-2731_2021.vcf", skip = 12)[,-10]
test4_variants = read.delim("/Volumes/Nick Minor External HD 1/E484T/data/USA_WI-UW-2731-T2_2021.vcf", skip = 12)[,-10]
test5_variants = read.delim("/Volumes/Nick Minor External HD 1/E484T/data/USA_WI-UW-2731-T3_2021.vcf", skip = 12)[,-10]
test6_variants = read.delim("/Volumes/Nick Minor External HD 1/E484T/data/USA_WI-UW-2731-T4_2021.vcf", skip = 12)[,-10]
crude_vcf_merged = rbind(test1_variants, test2_variants, test3_variants, test4_variants, test5_variants, test6_variants)

# creating new column for dates with only a sample ID in it (for the time being)
for (i in 1:nrow(crude_vcf_merged)){
  info <- crude_vcf_merged[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  crude_vcf_merged[i, "DATES_DETECTED"] <- info_split[grep("^VARSEQ", info_split)]
}
crude_vcf_merged$DATES_DETECTED <- str_replace_all(crude_vcf_merged$DATES_DETECTED, "VARSEQ=", "")

# replacing sample names with dates in date column
for (i in 1:nrow(crude_vcf_merged)){
  date <- crude_vcf_merged[i,"DATES_DETECTED"]
  ref <- fasta_df[match(date, fasta_df$Sample_ID),1:2]
  crude_vcf_merged[i,"DATES_DETECTED"] <- as.character(ref[1, "Date"])
}
crude_vcf_merged$DATES_DETECTED <- as.Date(crude_vcf_merged$DATES_DETECTED, format = "%Y-%m-%d")

# creating new column for days of infection with only a sample ID (for the time being)
for (i in 1:nrow(crude_vcf_merged)){
  info <- crude_vcf_merged[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  crude_vcf_merged[i, "DAYS_DETECTED"] <- info_split[grep("^VARSEQ", info_split)]
}
crude_vcf_merged$DAYS_DETECTED <- str_replace_all(crude_vcf_merged$DAYS_DETECTED, "VARSEQ=", "")

# replacing sample names with days of infection in day of infection column
for (i in 1:nrow(crude_vcf_merged)){
  day <- crude_vcf_merged[i,"DAYS_DETECTED"]
  ref <- fasta_df[match(day, fasta_df$Sample_ID),c(1,4)]
  crude_vcf_merged[i,"DAYS_DETECTED"] <- as.character(ref[1, "day_of_infection"])
}
crude_vcf_merged$DAYS_DETECTED <- as.numeric(crude_vcf_merged$DAYS_DETECTED)

# filtering down to only the earliest instance of every variant
variants <- crude_vcf_merged[-(1:nrow(crude_vcf_merged)),]
for (i in unique(crude_vcf_merged$POS)){
  sub <- crude_vcf_merged[crude_vcf_merged$POS==i,]
  keeper <- sub[sub$DAYS_DETECTED==min(sub$DAYS_DETECTED),]
  variants <- rbind(variants, keeper)
}

# putting together information to plot genes
SARS_genes <-
  data.frame(
    "name" = c("ORF1a", "ORF1b", "S", NA, NA, NA, NA, NA, NA, NA, "N"),
    "start" = c(
      266,
      13469,
      21563,
      25393,
      26245,
      26523,
      27202,
      27394,
      27894,
      28284,
      28578
    ),
    "stop" = c(
      13468,
      21555,
      25384,
      26220,
      26472,
      27191,
      27387,
      27759,
      28259,
      28577,
      29533
    ),
    "col" = c(
      "lightgreen",
      "orange",
      "lightblue",
      "lightgreen",
      "gold",
      "lightblue",
      "orange",
      "lightgreen",
      "lightblue",
      "gold",
      "orange"
    ),
    "left_out_names" = c(NA, NA, NA, "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "ORF9b", NA)
  )

plot(1,1, xlim = c(0,30000), ylim=c(0, 350), type = "n",
     xlab = "SARS-CoV-2 Genome Position", ylab = "Day of Infection", 
     frame.plot = F, cex.axis = 0.8, las = 1, xaxt = "n")
axis(side = 1, at = c(0,5000,10000,15000,20000,25000,30000), cex.axis = 0.8,
     col = "white", col.ticks = "black")
abline(v=SARS_genes$stop, col = "lightgrey")
abline(v=SARS_genes$start[1], col = "lightgrey")

for (i in 1:nrow(SARS_genes)){
  polygon(x = c(SARS_genes[i,"start"], SARS_genes[i,"start"], 
                SARS_genes[i,"stop"], SARS_genes[i,"stop"]),  # coordinates of gene 
          y = c(-9.5, 20, 20, -9.5),    # Y-Coordinates of polygon
          col = SARS_genes[i,"col"], border = F)
  text(x = SARS_genes[i,"start"] + (SARS_genes[i,"stop"]-SARS_genes[i,"start"])/2, 
       y = (20--9.5)/4, labels = SARS_genes[i, "name"], col = "black", cex = 0.8)
}

abline(h=198, col = "red")
# text(x = 29500, y = 210, labels = "Bamlanivumab administered", cex = 0.8, las = 3)
points(variants$POS, variants$DAYS_DETECTED, pch = 16, col = "darkgray", cex = 1.2)
points(variants[variants$POS==23012,"POS"], variants[variants$POS==23012,"DAYS_DETECTED"], 
       pch = 1, cex = 2, col = "gold")


# CIRCOS PLOTTING ####
SARS_genes_circos <- data.frame("genome" = rep(1, 11),
                                "name" = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "ORF9b", "N"),
                         "start" = c(266, 13469, 21563, 25393, 26245, 26523, 27202, 27394, 27894, 28284, 28578),
                         "stop" = c(13468, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 28259, 28577, 29533))
write.csv(SARS_genes_circos, "/Users/nicholasminor/Documents/informatics/Circos_plotting/SARS_genes_circos.csv", row.names = F, quote = F)
variants_circos <- data.frame("genome" = rep(1, nrow(variants)),
                              "position" = variants$POS,
                              "days" = variants$DAYS_DETECTED)
write.csv(variants_circos, "/Users/nicholasminor/Documents/informatics/Circos_plotting/variants.csv", row.names = F, quote = F)

