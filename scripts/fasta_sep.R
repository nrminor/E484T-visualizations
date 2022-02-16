#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("The data filepath and temp filepath must be specified.", call.=FALSE)
} else if (length(args)>2) {
  stop("This R script only accepts two arguments.", call.=FALSE)
} 

### PREPARING THE ENVIRONMENT ####
library(Biostrings)
library(tidyverse)
data_filepath = args[1]  ### OR INSERT YOUR FILE PATH HERE ###
setwd(data_filepath)

fasta <- readDNAStringSet("b12_enriched_global_subsampled_sequences.fasta")
seq_names <- names(fasta)
fasta <- read.delim("b12_enriched_global_subsampled_sequences.fasta", header = F)
sample_locs <- data.frame("sample" = seq_names, "start" = which(grepl(">", fasta[,1])), "stop" = NA)
sample_locs$stop <- c((sample_locs$start - 1)[2:nrow(sample_locs)], nrow(fasta))

for (i in 1:nrow(sample_locs)){
  sub <- fasta[sample_locs$start[i]:sample_locs$stop[i],1]
  name <- str_replace_all(string = sample_locs$sample[i], pattern = "/", replacement = "_")
  write.table(sub, paste(getwd(), "/tmp/", name, ".fasta", sep = ""), 
              quote = F, row.names = F, col.names = F)
  print(paste("finished with sample", i, "of 4992;", sep = " "))
  print(paste(round(i/4992 * 100), "% of samples separated.", sep = ""))
}


