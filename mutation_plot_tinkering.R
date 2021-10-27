# SETUP ####
library(Biostrings)
library(tidyverse)
library(grid)
library(gridExtra)
library(gridGraphics)
data_filepath = "/Volumes/working_ssd/e484t_manuscript/"
setwd(data_filepath)

# FASTA PREP ####
patient_fasta <- readDNAStringSet("data/alltimepoints_20211025.fasta")
seq_names = names(patient_fasta)
sequence = paste(patient_fasta)
fasta_df <- data.frame(seq_names, sequence)
fasta_df$seq_names <- str_replace_all(fasta_df$seq_names, fixed(" | "), ",")
fasta_df <- separate(data = fasta_df, col = 1,
                     into = c("Sample_ID", "Date"),
                     sep = ",")
fasta_df$Date <- as.Date(fasta_df$Date, "%Y %b %d")
fasta_df$day_of_infection <- c(113,124,131,159,198,297,333,388)

# VCF PREP ####
crude_vcf = read.delim("data/alltimepoints_variants_20211021.vcf", skip = 13)
crude_vcf <- crude_vcf[grepl("VARSEQ",crude_vcf$INFO, fixed=T),] ; rownames(crude_vcf) <- NULL

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

# MUTATION LISTING ####
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

# GENE INFO AND LABELS ####
SARS_genes <- data.frame("name" = c("ORF1a", "ORF1b", "S", NA, NA, NA, NA, NA, NA, NA, "N"),
                         "start" = c(266, 13469, 21563, 25393, 26245, 26523, 27202, 27394, 27894, 28284, 28578),
                         "stop" = c(13468, 21555, 25384, 26220, 26472, 27191, 27387, 27759, 28259, 28577, 29533),
                         "col" = c("#648FFF", "#FFB000", "#DC267F", "#FE6100", "#648FFF", "#FFB000", "#DC267F", 
                                   "#FE6100", "#648FFF", "#FFB000", "#DC267F"),
                         "left_out_names" = c(NA, NA, NA, "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "ORF9b", NA))

# HORIZONTAL PLOT ####
pdf("/Users/nicholasminor/Documents/informatics/E484T_paper/visuals/mutations_across_genome.pdf",
    width = 8, height = 8)
{plot(1,1, xlim = c(0,30000), ylim=c(0, 400), type = "n",
      xlab = "SARS-CoV-2 Genome Position", ylab = "Day of Infection", 
      frame.plot = F, cex.axis = 1, las = 1, xaxt = "n")
  axis(side = 1, at = c(0,5000,10000,15000,20000,25000,30000), cex.axis = 1,
       col = "white", col.ticks = "black")
  segments(x0=SARS_genes$stop, 
           y0 = rep(-10, times = length(SARS_genes$stop)),
           y1 = rep(400, times = length(SARS_genes$stop)),
           col = "lightgrey")
  segments(x0 = SARS_genes$start[1], y0 = -10, y1 = 400, col = "lightgrey")
  
  for (i in 1:nrow(SARS_genes)){
    polygon(x = c(SARS_genes[i,"start"], SARS_genes[i,"start"], 
                  SARS_genes[i,"stop"], SARS_genes[i,"stop"]),  # coordinates of gene 
            y = c(-9.5, 20, 20, -9.5),    # Y-Coordinates of polygon
            col = SARS_genes[i,"col"], border = F)
    text(x = SARS_genes[i,"start"] + (SARS_genes[i,"stop"]-SARS_genes[i,"start"])/2, 
         y = (20--9.5)/4, labels = SARS_genes[i, "name"], col = "black", cex = 1)
  }
  
  for (i in 1:length(mut_list2)){
    plotting_data <- data.frame("Date" = mut_list2[[i]],
                                "Coordinate" = rep(as.numeric(names(mut_list2)[i]), times = length(mut_list2[[i]])))
    points(x = plotting_data$Coordinate, plotting_data$Date, type = "p", pch = 20, col = "#e5e5e5")
  }
  
  fixed <- rep(NA, times = length(mut_list2))
  for (i in 1:length(fixed)){
    fixed[i] <- length(mut_list2[[i]])==8
  }
  
  segments(x0 = -1200, x1 = 29800, y0=198, col = "red")
  # text(x = 29500, y = 210, labels = "Bamlanivumab administered", cex = 1, las = 3)
  points(crude_vcf[!fixed, "POS"], crude_vcf[!fixed, "EARLIEST_DAY"], pch = 16, col = "#7e7e7e", cex = 1.2)
  points(23012, crude_vcf[crude_vcf$POS==23012,"EARLIEST_DAY"], 
         pch = 1, cex = 2, col = "gold")
  text(23000, 320, labels = "E484T", col = "#7e7e7e", cex = 1)
  
  legend(15000,400, legend=c("fixed mutation", "novel mutation", "antibody treatment"), 
         pch = c(20, 16, NA), lty = c(NA, NA, 1),
         col = c("#e5e5e5", "#7e7e7e", "red"),
         bty = "o", cex = 1, box.col = "white", ncol = 3, xpd = T, xjust = 0.5)
}
dev.off() 

# FLIPPED VERTICAL PLOT WITH BIGGER FONTS ###
pdf("/Users/nicholasminor/Documents/informatics/E484T_paper/visuals/mutations_across_genome.pdf",
    width = 8, height = 8)
{plot(1,1, xlim = c(0,30000), ylim=c(400, 0), type = "n",
      ylab = "Day of Infection", xlab='',
      frame.plot = F, cex.axis = 1, las = 1, xaxt = "n")
  axis(side = 3, at = c(0,5000,10000,15000,20000,25000,30000), cex.axis = 1,
       col = "white", col.ticks = "black")
  segments(x0=SARS_genes$stop, 
           y0 = rep(-10, times = length(SARS_genes$stop)),
           y1 = rep(400, times = length(SARS_genes$stop)),
           col = "lightgrey")
  segments(x0 = SARS_genes$start[1], y0 = -10, y1 = 400, col = "lightgrey")
  
  for (i in 1:nrow(SARS_genes)){
    polygon(x = c(SARS_genes[i,"start"], SARS_genes[i,"start"], 
                  SARS_genes[i,"stop"], SARS_genes[i,"stop"]),  # coordinates of gene 
            y = c(-9.5, 20, 20, -9.5),    # Y-Coordinates of polygon
            col = SARS_genes[i,"col"], border = F)
    text(x = SARS_genes[i,"start"] + (SARS_genes[i,"stop"]-SARS_genes[i,"start"])/2, 
         y = (20--9.5)/4, labels = SARS_genes[i, "name"], col = "black", cex = 1)
  }
  
  for (i in 1:length(mut_list2)){
    plotting_data <- data.frame("Date" = mut_list2[[i]],
                                "Coordinate" = rep(as.numeric(names(mut_list2)[i]), times = length(mut_list2[[i]])))
    points(x = plotting_data$Coordinate, plotting_data$Date, type = "p", pch = 20, col = "#e5e5e5")
  }
  
  fixed <- rep(NA, times = length(mut_list2))
  for (i in 1:length(fixed)){
    fixed[i] <- length(mut_list2[[i]])==8
  }
  
  segments(x0 = -1200, x1 = 29800, y0=198, col = "red")
  # text(x = 29500, y = 210, labels = "Bamlanivumab administered", cex = 1, las = 3)
  points(crude_vcf[!fixed, "POS"], crude_vcf[!fixed, "EARLIEST_DAY"], pch = 16, col = "#7e7e7e", cex = 1.2)
  points(23012, crude_vcf[crude_vcf$POS==23012,"EARLIEST_DAY"], 
         pch = 1, cex = 2, col = "gold")
  text(23000, 320, labels = "E484T", col = "#7e7e7e", cex = 1)
  
  legend(15000,370, legend=c("fixed mutation", "novel mutation", "antibody treatment"), 
         pch = c(20, 16, NA), lty = c(NA, NA, 1),
         col = c("#e5e5e5", "#7e7e7e", "red"),
         bty = "o", cex = 1, box.col = "white", ncol = 3, xpd = T, xjust = 0.5)
}
# dev.off() 

# VERTICAL PLOT ####
pdf("/Users/nicholasminor/Documents/informatics/E484T_paper/visuals/mutations_across_genome_vertical.pdf",
    width = 8, height = 7)
mut_plot <- function(){plot(1,1, ylim = c(0,30000), xlim=c(0, 400), type = "n",
      ylab = "SARS-CoV-2 Genome Position", xlab = NA, 
      frame.plot = F, cex.axis = 1, yaxt = "n")
  axis(side = 2, at = c(0,5000,10000,15000,20000,25000,30000), cex.axis = 0.8,
       col = "white", col.ticks = "black", las = 1)
  segments(y0=SARS_genes$stop, 
           x0 = rep(-10, times = length(SARS_genes$stop)),
           x1 = rep(400, times = length(SARS_genes$stop)),
           col = "lightgrey")
  segments(y0 = SARS_genes$start[1], x0 = -10, x1 = 400, col = "lightgrey")
  
  for (i in 1:nrow(SARS_genes)){
    polygon(y = c(SARS_genes[i,"start"], SARS_genes[i,"start"], 
                  SARS_genes[i,"stop"], SARS_genes[i,"stop"]),  # coordinates of gene 
            x = c(-9.5, 20, 20, -9.5),    # Y-Coordinates of polygon
            col = SARS_genes[i,"col"], border = F)
    text(y = SARS_genes[i,"start"] + (SARS_genes[i,"stop"]-SARS_genes[i,"start"])/2, 
         x = (20-9.5)/4, labels = SARS_genes[i, "name"], srt = 90, col = "black", cex = 1)
  }
  
  for (i in 1:length(mut_list2)){
    plotting_data <- data.frame("Day" = mut_list2[[i]],
                                "Coordinate" = rep(as.numeric(names(mut_list2)[i]), times = length(mut_list2[[i]])))
    points(y = plotting_data$Coordinate, x = plotting_data$Day, type = "p", pch = 20, col = "#e5e5e5")
  }
  
  fixed <- rep(NA, times = length(mut_list2))
  for (i in 1:length(fixed)){
    fixed[i] <- length(mut_list2[[i]])==8
  }
  
  segments(y0 = -1200, y1 = 29800, x0=198, col = "red")
  # text(x = 29500, y = 210, labels = "Bamlanivumab administered", cex = 1, las = 3)
  points(y = crude_vcf[!fixed, "POS"], x = crude_vcf[!fixed, "EARLIEST_DAY"], pch = 16, col = "#7e7e7e", cex = 1.2)
  
  for (i in crude_vcf[!fixed, "POS"]){
    days <- unlist(mut_list2[as.numeric(names(mut_list2))==i])
    days <- as.data.frame(days)
    coordinates <- rep(as.numeric(names(mut_list2[as.numeric(names(mut_list2))==i])),
                       times = length(days$days))
    points(y = coordinates, x = days$days, pch = 16, col = "#7e7e7e", cex = 1.2)
  }
  
  points(y = 23012, x = crude_vcf[crude_vcf$POS==23012,"EARLIEST_DAY"], 
         pch = 1, cex = 2, col = "gold")
  points(y = 23012, x = 333, 
         pch = 1, cex = 2, col = "gold")
  
  text(y = 23700, x = 316, labels = "E484T", col = "#7e7e7e", cex = 1)
  
  # legend(200, 28000, legend=c("ancestral", "derived during infection"), 
  #        pch = c(20, 16),
  #        col = c("#e5e5e5", "#7e7e7e", "red"),
  #        bty = "n", cex = 1, bg="transparent", ncol = 2, xpd = T, 
  #        xjust = 0.5, yjust = 0)
}
mut_plot()
# grid.echo()
# mut_plot <- grid.grab()
# mut_plot <- recordPlot()
dev.off()

