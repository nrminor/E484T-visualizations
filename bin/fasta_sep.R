#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("The working directory and path to data must be specified.", call.=FALSE)
} else if (length(args)>1) {
  stop("This R script only accepts one argument: The path to the merged FASTA file of consensus sequences.", call.=FALSE)
} 

### PREPARING THE ENVIRONMENT ####
library(Biostrings)
library(tidyverse)

data_filepath = args[1]

fasta <- readDNAStringSet(data_filepath)
seq_names <- names(fasta)

for (i in 1:length(seq_names)){
  if (grepl(" ", seq_names[i])){
    seq_names[i] <- strsplit(seq_names[i], " | ")[[1]][1]
  } else {next}
}

fasta <- read.delim(data_filepath, header = F)
sample_locs <- data.frame("sample" = seq_names, "start" = which(grepl(">", fasta[,1])), "stop" = NA)
sample_locs$stop <- c((sample_locs$start - 1)[2:nrow(sample_locs)], nrow(fasta))

for (i in 1:nrow(sample_locs)){
    seq_name <- paste(">", seq_names[i], sep = "")
    sub <- fasta[sample_locs$start[i]:sample_locs$stop[i],1]
    sub <- sub[-1]
    sub <- str_replace_all(sub, "-", "N")
    sub <- c(seq_name[1], sub)
    filename <- str_replace_all(string = sample_locs$sample[i], pattern = "/", replacement = "_")
    write.table(sub, paste(filename, ".fasta", sep = ""), 
                quote = F, row.names = F, col.names = F)
    print(paste("finished with sample", i, "of", 
                length(seq_names), sep = " "))
    print(paste(round(i/length(seq_names) * 100), 
                "% of samples separated.", sep = ""))
  }



