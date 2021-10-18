
### PREPARING THE ENVIRONMENT ###
library(Biostrings)
library(tidyverse)
library(RColorBrewer)
data_filepath = "/Volumes/working_ssd/e484t_manuscript/"
setwd(data_filepath)



### IMPORTING VCF & METADATA
gisaid_vcf <- read.delim("data/b112_enriched_subsampled_aligned.vcf", skip = 13)
gisaid_meta <- read.delim("data/b112_enriched_subsampled_metadata.tsv")
strain_dates <- gisaid_meta[,c(1,4)]
strain_dates$date <- as.Date(strain_dates$date, format = "%Y-%m-%d")
strain_dates$day_of_infection <- as.numeric(NA)
for (i in 1:nrow(strain_dates)){
  strain_dates$day_of_infection[i] <- as.numeric(strain_dates$date[i] - as.Date("2020-09-05"))
  
}



### PREPARING THE GISAID SUBSAMPLE FASTA ###
gisaid_fasta <- readDNAStringSet("data/b112_enriched_subsampled_aligned.fasta")
seq_names = names(gisaid_fasta)
sequence = paste(gisaid_fasta)
gisaid_fasta_df <- data.frame(seq_names, sequence)
remove(gisaid_fasta)



### PREPARING METADATA ###
gisaid_fasta_df <- gisaid_fasta_df[match(strain_dates$strain, gisaid_fasta_df$seq_names),] ; rownames(gisaid_fasta_df) <- NULL
gisaid_meta <- gisaid_meta[match(strain_dates$strain, gisaid_meta$strain),] ; rownames(gisaid_meta) <- NULL

which(is.na(gisaid_fasta_df$seq_names)) # -> to_remove
which(is.na(gisaid_meta$strain)) # -> to_remove2
{ # fixing the issue at 4940
  gisaid_meta$strain[4935:4945] == gisaid_fasta_df$seq_names[4935:4945]
  gisaid_meta$strain[4940] <- str_replace(gisaid_meta$strain[4940], "United_States", "USA")
  gisaid_fasta_df$seq_names[4940] <- gisaid_meta$strain[4940]
}
which(is.na(gisaid_fasta_df$seq_names))
which(is.na(gisaid_meta$strain))

gisaid_fasta_df$date <- as.Date(strain_dates$date, format = "%Y-%m-%d")
gisaid_meta$date <- as.Date(strain_dates$date, format = "%Y-%m-%d")

gisaid_fasta_df$day_of_infection <- strain_dates$day_of_infection
gisaid_meta$day_of_infection <- strain_dates$day_of_infection

gisaid_fasta_df$distance <- rep(0, nrow(gisaid_fasta_df))



### IDENTIFYING VARIANT NAMES & TYPES IN FASTA ###
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



### REMOVING INDELS FROM GISAID VCF, LEAVING ONLY SNVs ###
gisaid_vcf <- gisaid_vcf[!grepl("Deletion", gisaid_vcf$MUTATION_TYPE),]



### COUNTING MUTATIONS FOR EACH SAMPLE ###
for (i in 1:length(gisaid_vcf$SAMPLES)){
  sub <- unlist(strsplit(gisaid_vcf$SAMPLES[i], split = ","))
  for (j in sub){
    add <- gisaid_fasta_df[gisaid_fasta_df$seq_names==j,"distance"] + 1
    gisaid_fasta_df[gisaid_fasta_df$seq_names==j,"distance"] <- add
  }
} 



### SIMPLE LINEAR MODEL OF MUTATION INCREASE ###
# VOC_lm <- lm(gisaid_fasta_df$distance ~ gisaid_fasta_df$date) 



### IDENTIFYING PANGO LINEAGES IN GISAID SUBSAMPLE ###
lineages <- as.data.frame(table(gisaid_meta$Pango.lineage))
colnames(lineages) <- c("lineage", "count")
lineages <- lineages[order(lineages$count, decreasing = T),] ; rownames(lineages) <- NULL
total_count <- sum(lineages$count) ; total_count == nrow(gisaid_fasta_df)
sum(lineages[1:7,"count"])/total_count # try to get at least 50% of the samples colored by lineage

palette <- read.delim("/Users/nicholasminor/Documents/informatics/E484T_paper/palette.txt", header = FALSE, sep = "\t")
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



### MOVING LINEAGE DATA INTO FASTA DF ###
gisaid_fasta_df$lineage <- gisaid_meta$Pango.lineage
gisaid_fasta_df$color <- ""
gisaid_fasta_df$label <- ""
for (i in 1:nrow(gisaid_fasta_df)){
  gisaid_fasta_df[i,"color"] <- lineages[lineages$lineage==gisaid_fasta_df$lineage[i],"color"]
  gisaid_fasta_df[i,"label"] <- lineages[lineages$lineage==gisaid_fasta_df$lineage[i],"label"]
}
rgb <- col2rgb(gisaid_fasta_df$color, alpha = T)
rgb[4,] <- round(rgb[4,]/2)
rgb <- unlist(rgb)
gisaid_fasta_df$rgb <- ""
for (i in 1:ncol(rgb)){
  sub <- rgb[,i]
  names(sub) <- NULL
  gisaid_fasta_df$rgb[i] <- paste(sub, collapse = ",")
}


### IMPORTING PATIENT DATA ###
patient_fasta <- readDNAStringSet("data/PT0001_alltimepoints_20210928.fasta")
seq_names <- names(patient_fasta)
sequence <- paste(patient_fasta)
fasta_df <- data.frame(seq_names, sequence)
fasta_df$seq_names <- str_replace_all(fasta_df$seq_names, fixed(" | "), ",")
fasta_df$sequence <- str_replace_all(fasta_df$sequence, fixed("-"), "N")
fasta_df <- separate(data = fasta_df, col = 1,
                     into = c("Sample_ID", "Date"),
                     sep = ",")
fasta_df$Date <- as.Date(fasta_df$Date, "%Y %b %d")
fasta_df$day_of_infection <- c(113,124,131,159,198,297,388)



# COUNTING MUTATIONS ####
# Counting the number of variants at each sampling event
test1_variants = read.delim("data/USA_WI-WSLH-202168_2020.vcf", skip = 12)[,-10]
test2_variants = read.delim("data/USA_WI-UW-5350_2021.vcf", skip = 12)[,-10]
test3_variants = read.delim("data/USA_WI-UW-2731_2021.vcf", skip = 12)[,-10]
test4_variants = read.delim("data/USA_WI-UW-2731-T2_2021.vcf", skip = 12)[,-10]
test5_variants = read.delim("data/USA_WI-UW-2731-T3_2021.vcf", skip = 12)[,-10]
test6_variants = read.delim("data/USA_WI-UW-2731-T4_2021.vcf", skip = 12)[,-10]
test7_variants = read.delim("data/USA_CA-Spietz-09282021.vcf", skip = 12)[,-10]
crude_vcf_merged = rbind(test1_variants, test2_variants, test3_variants, test4_variants, test5_variants, test6_variants, test7_variants)
crude_vcf_merged <- crude_vcf_merged[grepl("VARSEQ",crude_vcf_merged$INFO, fixed=T),] ; rownames(crude_vcf_merged) <- NULL

# creating new column for dates with only a sample ID in it (for the time being)
for (i in 1:nrow(crude_vcf_merged)){
  info <- crude_vcf_merged[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  crude_vcf_merged[i, "SAMPLE_ID"] <- info_split[grep("^VARSEQ", info_split)]
}
crude_vcf_merged$SAMPLE_ID <- str_replace_all(crude_vcf_merged$SAMPLE_ID, "VARSEQ=", "")
crude_vcf_merged$SAMPLE_ID <- str_replace(crude_vcf_merged$SAMPLE_ID, "-consensus", "")
crude_vcf_merged$SAMPLE_ID <- str_replace(crude_vcf_merged$SAMPLE_ID, "-DBaker", "")

# creating new column with only mutation type
for (i in 1:nrow(crude_vcf_merged)){
  info <- crude_vcf_merged[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  crude_vcf_merged[i, "MUTATION_TYPE"] <- info_split[grep("^TYPE", info_split)]
}
crude_vcf_merged$MUTATION_TYPE <- str_replace_all(crude_vcf_merged$MUTATION_TYPE, "TYPE=", "")

# summing up mutations per sample
fasta_df$distance <- NA
for (i in unique(crude_vcf_merged$SAMPLE_ID)){
  fasta_df[fasta_df$Sample_ID==i, "distance"] <- nrow(crude_vcf_merged[crude_vcf_merged$SAMPLE_ID==i,])
}



### SIMPLE LINEAR MODEL OF PATIENT MUTATIONS THROUGH TIME
# patient_lm <- lm(fasta_df$distance ~ fasta_df$Date)



### PLOTTING ###
plot(fasta_df$Date, fasta_df$distance, 
     xlim = c(min(gisaid_fasta_df$date), max(gisaid_fasta_df$date)+31), ylim = c(0,60),
     xlab = "Date of Infection (Sept. 2020 - Sept. 2021)", ylab = "Genetic Distance from Wuhan-1", 
     frame.plot = F, cex.axis = 0.8, cex.lab = 0.85, las = 1, pch = 20, xaxt="n",
     cex = 3, col = "#4B7395", type = "n")
grid()

months_axis <- seq(min(gisaid_fasta_df$date), max(gisaid_fasta_df$date)+31, by = "month")
axis(side = 1, at = months_axis,
     labels = format(months_axis, "%b"), cex.axis = 0.8)

points(gisaid_fasta_df$date, gisaid_fasta_df$distance,
       pch = 20, cex = 1, col = gisaid_fasta_df$color)

points(gisaid_fasta_df$date, gisaid_fasta_df$distance,
       pch = 20, cex = 1, col = gisaid_fasta_df$color)

text(as.Date("2021-03-22"), 45, labels = "Bamlanivimab\nTreatment", cex = 0.8, bty = "l")
segments(as.Date("2021-03-22"), y0 = -2, y1 = 42, col = "red", lty = 2)
segments(as.Date("2021-03-22"), y0 = 48, y1 = 61, col = "red", lty = 2)

# abline(VOC_lm, col = rgb(165/255,15/255,21/255,3/4), lwd = 4)

points(fasta_df$Date, fasta_df$distance,
       pch = 20, cex = 2, col = palette[11,1])
# abline(patient_lm, col = rgb(8/255,81/255,156/255,3/4), lwd = 4)
legend("topleft", legend = c(lineages$label[1:5], "patient"),
       col = c(lineages$color[1:5], palette[11,1]), bty="n",
       pch = 16, ncol = 2, xpd = T, xjust = 0.5, cex = 0.9)

