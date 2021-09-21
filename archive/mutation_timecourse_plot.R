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

# HOMEMADE FIGURE ####
mutations = read.vcfR("/Volumes/Nick Minor External HD 1/E484T/data/alltimepoints_variants.vcf")

crude_vcf = read.delim("/Volumes/Nick Minor External HD 1/E484T/data/alltimepoints_variants.vcf", skip = 13)
for (i in 1:nrow(crude_vcf)){
  info <- crude_vcf[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  crude_vcf[i, "DATES_DETECTED"] <- info_split[grep("^VARSEQ", info_split)]
}
crude_vcf$DATES_DETECTED <- str_replace_all(crude_vcf$DATES_DETECTED, "VARSEQ=", "")

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

for (i in 1:nrow(crude_vcf)){
  info <- crude_vcf[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  crude_vcf[i, "DAYS_DETECTED"] <- info_split[grep("^VARSEQ", info_split)]
}
crude_vcf$DAYS_DETECTED <- str_replace_all(crude_vcf$DAYS_DETECTED, "VARSEQ=", "")
for (i in 1:nrow(crude_vcf)){
  varseq <- crude_vcf[i,"DAYS_DETECTED"]
  samples <- unlist(strsplit(varseq, split = ","))
  samples_ordered <- samples[match(fasta_df$Sample_ID, samples)]
  samples_ordered <- as.character(na.omit(samples_ordered))
  days_ordered <- samples_ordered
  ref <- fasta_df[match(days_ordered, fasta_df$Sample_ID),c(1,4)]
  for (j in 1:length(days_ordered)){ 
    days_ordered[j] <- as.character(ref[j, "day_of_infection"])
  }
  if (length(days_ordered)>1){
    days_ordered <- paste0(days_ordered, collapse = ",")
  }
  crude_vcf[i,"DAYS_DETECTED"] <- as.character(days_ordered)
}

# plotting based on date:
mutation_dates <- matrix(ncol = length(crude_vcf), nrow = 30000)
mutation_dates <- as.data.frame(mutation_dates)
mutation_dates$coordinate <- 1:30000
mutation_dates$mut_type <- NA

mut_list <- vector("list", length = nrow(crude_vcf))
names(mut_list) <- as.numeric(crude_vcf$POS)
for (i in 1:nrow(crude_vcf)){
  sub <- crude_vcf[i, "DATES_DETECTED"]
  sub_split <- unlist(strsplit(sub, split = ","))
  dates <- as.Date(sub_split, format = "%Y-%m-%d")
  mut_list[[i]] <- dates
}

plot(1,fasta_df$Date[3], xlim = c(0,30000), ylim=c(fasta_df$Date[1], fasta_df$Date[6]), type = "n",
     xlab = "SARS-CoV-2 Genome Position", ylab = "Mutation Date", frame.plot = F, cex.axis = 0.8, las = 1)
for (i in 1:length(mut_list)){
  plotting_data <- data.frame("Date" = mut_list[[i]],
                              "Coordinate" = rep(as.numeric(names(mut_list)[i]), times = length(mut_list[[i]])))
  points(x = plotting_data$Coordinate, plotting_data$Date, type = "p", pch = 20)
}

# plotting based on day of infection:
mut_list <- vector("list", length = nrow(crude_vcf))
names(mut_list) <- as.numeric(crude_vcf$POS)
for (i in 1:nrow(crude_vcf)){
  sub <- crude_vcf[i, "DAYS_DETECTED"]
  sub_split <- unlist(strsplit(sub, split = ","))
  days <- as.numeric(sub_split)
  mut_list[[i]] <- days
}

plot(1,1, xlim = c(0,30000), ylim=c(0, 350), type = "n",
     xlab = "SARS-CoV-2 Genome Position", ylab = "Day of Infection", frame.plot = F, cex.axis = 0.8, las = 1)
for (i in 1:length(mut_list)){
  plotting_data <- data.frame("Date" = mut_list[[i]],
                              "Coordinate" = rep(as.numeric(names(mut_list)[i]), times = length(mut_list[[i]])))
  points(x = plotting_data$Coordinate, plotting_data$Date, type = "p", pch = 20)
}

# plotting first day of infection when particular mutations were detected:
mut_list <- vector("list", length = nrow(crude_vcf))
names(mut_list) <- as.numeric(crude_vcf$POS)
for (i in 1:nrow(crude_vcf)){
  sub <- crude_vcf[i, "DAYS_DETECTED"]
  sub_split <- unlist(strsplit(sub, split = ","))
  days <- as.numeric(sub_split)
  mut_list[[i]] <- days
}

plot(1,1, xlim = c(0,30000), ylim=c(0, 350), type = "n",
     xlab = "SARS-CoV-2 Genome Position", ylab = "Day of Infection", frame.plot = F, cex.axis = 0.8, las = 1)
for (i in 1:length(mut_list)){
  plotting_data <- data.frame("Date" = mut_list[[i]][1],
                              "Coordinate" = as.numeric(names(mut_list)[i]))
  points(x = plotting_data$Coordinate, plotting_data$Date, type = "p", pch = 20)
}




# end of section

# GViz ####
gtrack <- GenomeAxisTrack(range = IRanges(start = 21563,
                                          end = 25384,
                                          names = "spike"))
plotTracks(gtrack, from = 0, to = 29533, showId = TRUE, labelPos = "below")


# GenVisR Lollipop Plots
data <- brcaMAF[brcaMAF$Hugo_Symbol == 'TP53',c('Hugo_Symbol', 'amino_acid_change_WU')]
data <- as.data.frame(cbind(data, 'ENST00000269305'))
colnames(data) <- c('gene', 'amino_acid_change', 'transcript_name')

# Call lolliplot
lolliplot(data)

fall <- Waterfall(vepObject, recurrence = 0.4)


# CIRCOS ####
data("UCSC.HG19.Human.CytoBandIdeogram")
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
RCircos.Set.Core.Components(cyto.info, chr.exclude=NULL, tracks.inside = 10,
                            tracks.outside = 0)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

data("RCircos.Gene.Label.Data")
RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data, track.num = 1, side = "in")
RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data, name.col = 4, track.num = 2, side = "in")















# CODE GRAVEYARD ####
# SARS2_gene_track <- GenomeAxisTrack(range = IRanges(start = c(266, 13469, 21563, 25393, 26245, 26523, 27202, 27394, 27760, 27894, 28284, 28285),
#                                           end = c(13468, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 27887, 28259, 28577, 29533),
#                                           names = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "ORF9b", "N")))
# plotTracks(SARS2_gene_track, from = 0, to = 29533, showId = TRUE)
