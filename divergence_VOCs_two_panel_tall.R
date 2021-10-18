library(Biostrings)
library(tidyverse)
data_filepath = "/Volumes/working_ssd/e484t_manuscript"
setwd(data_filepath)

# READING IN FASTA ####
# reading in/re-formatting fasta
patient_fasta <- readDNAStringSet("data/PT0001_alltimepoints_20210928.fasta")
seq_names = names(patient_fasta)
sequence = paste(patient_fasta)
fasta_df <- data.frame(seq_names, sequence)
fasta_df$seq_names <- str_replace_all(fasta_df$seq_names, fixed(" | "), ",")
fasta_df <- separate(data = fasta_df, col = 1,
                     into = c("Sample_ID", "Date"),
                     sep = ",")
fasta_df$Date <- as.Date(fasta_df$Date, "%Y %b %d")
fasta_df$day_of_infection <- c(113,124,131,159,198,297)

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

# summing up mutations per sample
fasta_df$distance <- NA
for (i in unique(crude_vcf_merged$SAMPLE_ID)){
  fasta_df[fasta_df$Sample_ID==i, "distance"] <- nrow(crude_vcf_merged[crude_vcf_merged$SAMPLE_ID==i,])
}

# PREPARING MUTATION DATA ####
crude_vcf = read.delim("data/alltimepoints_variants.vcf", skip = 13)

# creating a column with all dates of detection for each mutation
for (i in 1:nrow(crude_vcf)){
  info <- crude_vcf[i, "INFO"]
  info_split <- unlist(strsplit(info, split = ";"))
  crude_vcf[i, "DATES_DETECTED"] <- info_split[grep("^VARSEQ", info_split)]
}
crude_vcf$DATES_DETECTED <- str_replace_all(crude_vcf$DATES_DETECTED, "VARSEQ=", "")

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

# creating a column with all *days* (note dates) of detection, from the patient's first day of infection
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
  ref <- fasta_df[match(days_ordered, fasta_df$Sample_ID),c(1,4)]
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

# first for dates
mut_list1 <- vector("list", length = nrow(crude_vcf))
names(mut_list1) <- as.numeric(crude_vcf$POS)
for (i in 1:nrow(crude_vcf)){
  sub <- crude_vcf[i, "DATES_DETECTED"]
  sub_split <- unlist(strsplit(sub, split = ","))
  dates <- as.Date(sub_split, format = "%Y-%m-%d")
  mut_list1[[i]] <- dates
}

# second for days
mut_list2 <- vector("list", length = nrow(crude_vcf))
names(mut_list2) <- as.numeric(crude_vcf$POS)
for (i in 1:nrow(crude_vcf)){
  sub <- crude_vcf[i, "DAYS_DETECTED"]
  sub_split <- unlist(strsplit(sub, split = ","))
  days <- as.numeric(sub_split)
  mut_list2[[i]] <- days
}

SARS_genes <- data.frame("name" = c("ORF1a", "ORF1b", "S", NA, NA, NA, NA, NA, NA, NA, "N"),
                         "start" = c(266, 13469, 21563, 25393, 26245, 26523, 27202, 27394, 27894, 28284, 28578),
                         "stop" = c(13468, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 28259, 28577, 29533),
                         "col" = c("lightgreen", "orange", "lightblue", "lightgreen", "gold", "lightblue", "orange", 
                                   "lightgreen", "lightblue", "gold", "orange"),
                         "left_out_names" = c(NA, NA, NA, "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "ORF9b", NA))

# PLOTTING MUTATIONS ####
matrix.layout<-matrix(c(1,1,1,1,2,2),nrow=3,byrow = TRUE)
matrix.layout
layout(matrix.layout)
# layout.show(2)
par(mar=c(5.1, 4.1, 4.1, 2.1))

{plot(1,1, ylim = c(0,30000), xlim=c(0, 350), type = "n",
      ylab = "SARS-CoV-2 Genome Position", xlab = "Day of Infection", 
      frame.plot = F, cex.axis = 0.8, yaxt = "n")
  axis(side = 2, at = c(0,5000,10000,15000,20000,25000,30000), cex.axis = 0.8,
       col = "white", col.ticks = "black", las = 1)
  segments(y0=SARS_genes$stop, 
           x0 = rep(-10, times = length(SARS_genes$stop)),
           x1 = rep(350, times = length(SARS_genes$stop)),
           col = "lightgrey")
  segments(y0 = SARS_genes$start[1], x0 = -10, x1 = 350, col = "lightgrey")
  
  for (i in 1:nrow(SARS_genes)){
    polygon(y = c(SARS_genes[i,"start"], SARS_genes[i,"start"], 
                  SARS_genes[i,"stop"], SARS_genes[i,"stop"]),  # coordinates of gene 
            x = c(-9.5, 20, 20, -9.5),    # Y-Coordinates of polygon
            col = SARS_genes[i,"col"], border = F)
    text(y = SARS_genes[i,"start"] + (SARS_genes[i,"stop"]-SARS_genes[i,"start"])/2, 
         x = (20--9.5)/4, labels = SARS_genes[i, "name"], col = "black", cex = 0.8)
  }
  
  for (i in 1:length(mut_list2)){
    plotting_data <- data.frame("Date" = mut_list2[[i]],
                                "Coordinate" = rep(as.numeric(names(mut_list2)[i]), times = length(mut_list2[[i]])))
    points(y = plotting_data$Coordinate, x = plotting_data$Date, type = "p", pch = 20, col = "#e5e5e5")
  }
  
  fixed <- rep(NA, times = length(mut_list2))
  for (i in 1:length(fixed)){
    fixed[i] <- length(mut_list2[[i]])==6
  }
  
  segments(y0 = -1200, y1 = 29800, x0=198, col = "red")
  # text(x = 29500, y = 210, labels = "Bamlanivumab administered", cex = 0.8, las = 3)
  points(y = crude_vcf[!fixed, "POS"], x = crude_vcf[!fixed, "EARLIEST_DAY"], pch = 16, col = "#7e7e7e", cex = 1.2)
  points(y = 23012, x = crude_vcf[crude_vcf$POS==23012,"EARLIEST_DAY"], 
         pch = 1, cex = 2, col = "gold")
  text(y = 23000, x = 315, labels = "E484T", col = "#7e7e7e", cex = 0.8)
  
  legend(175, 34000, legend=c("fixed mutation", "novel mutation", "antibody treatment"), 
         pch = c(20, 16, NA), lty = c(NA, NA, 1),
         col = c("#e5e5e5", "#7e7e7e", "red"),
         bty = "o", cex = 0.8, box.col = "white", ncol = 3, xpd = T, xjust = 0.5)
}

# PLOTTING DIVERGENCE ####
{par(mar=c(5.1,4.1,2,2.1))
  set.seed(14)
  simulation1 <- data.frame("distance" = rnorm(1000, 35, 5),
                            "day_of_infection" = rnorm(1000, 100, 100))
  simulation2 <- data.frame("distance" = rnorm(1000, 40, 5),
                            "day_of_infection" = rnorm(1000, 200, 100))
  simulation <- rbind(simulation1, simulation2)
  simul_lm <- lm(simulation$distance ~ simulation$day_of_infection)
  patient_lm <- lm(fasta_df$distance ~ fasta_df$day_of_infection)
  
  lineages <- as.data.frame(table(gisaid_meta$Pango.lineage))
  
  plot(fasta_df$day_of_infection, fasta_df$distance, ylim = c(25, 50), xlim = c(0, 350),
       xlab = NA, ylab = "Genetic Distance from Wuhan Variant", frame.plot = F,
       cex.axis = 0.8, cex.lab = 0.85, las = 1, pch = 20, cex = 3, col = "#4B7395", type = "n", xaxt = "n")
  axis(side = 3, at = c(0,50,100,150,200,250,300,350), cex.axis = 0.8,
       col.ticks = "black", las = 1)
  grid()
  text(198, 45, labels = "Antibody\nTreatment", cex = 0.8, bty = "l")
  segments(x0 = 198, y0 = 22, y1 = 43, col = "red", lty = 2)
  segments(x0 = 198, y0 = 47, y1 = 51, col = "red", lty = 2)
  
  points(simulation$day_of_infection, simulation$distance,
         pch = 20, cex = 1, col = rgb(75/255,115/255,149/255,0.5))
  abline(simul_lm, col = rgb(75/255,115/255,149/255,1), lwd = 4)
  
  points(fasta_df$day_of_infection, fasta_df$distance,
         pch = 20, cex = 2, col = rgb(227/255,152/255,58/255, 0.8))
  abline(patient_lm, col = rgb(227/255,152/255,58/255,1), lwd = 4)
  legend(175, 24, legend = c("US VOCs", "Patient Mutations"),
         col = c(rgb(75/255,115/255,149/255,1), rgb(227/255,152/255,58/255,1)),
         lwd = c(4,4), ncol = 2, xpd = T, xjust = 0.5, bty = "n", cex = 0.8)
}

