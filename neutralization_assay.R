# SETUP ####
library(Biostrings)
library(tidyverse)
library(grid)
library(gridExtra)
library(gridGraphics)
library(ggplot2)
library(plotrix)
data_filepath = "/Volumes/working_ssd/e484t_manuscript/"
setwd(data_filepath)



# FASTA PREP ####
patient_fasta <- readDNAStringSet("data/alltimepoints_20211104.fasta")
seq_names = names(patient_fasta)
sequence = paste(patient_fasta)
fasta_df <- data.frame(seq_names, sequence)
fasta_df$seq_names <- str_replace_all(fasta_df$seq_names, fixed(" | "), ",")
fasta_df <- separate(data = fasta_df, col = 1,
                     into = c("Sample_ID", "Date"),
                     sep = ",")
fasta_df$Date <- as.Date(fasta_df$Date, "%Y %b %d")
fasta_df$day_of_infection <- c(113,124,131,159,198,297,333,388)



# ANTIBODY DATA ####
antibody <- read.csv("data/antibody_potency.csv")
antibody$Chronic[1] <- 2.5


# PLOTTING ###
pdf(file = "/Users/nicholasminor/Documents/informatics/E484T_paper/visuals/neutralization_assay.pdf", 
    width = 8, height = 8)
plot(c(1:4), antibody$S.614G, type = "n",
     ylim = c(0, 2.5), frame.plot = F,
     xlab = NA, xaxt="n", ylab = "IC99 values µg/ml")
grid()
axis(side = 1, at = c(1,2,3,4), labels = row.names(antibody), 
     cex.axis = 0.8)

palette <- data.frame("points" = c("#2D339F", "#E5B326"),
                      "lines" = c("#484DAA", "#F4C950"))
lines(c(1:4), antibody$Chronic, col = palette$lines[2], lwd = 2)
points(c(1:4), antibody$S.614G, col = palette$points[1], pch = 0, cex = 1.2)
points(c(1:4), antibody$Chronic, col = palette$points[2], pch = 2, cex = 1.2)


dev.off()

pdf(file = "/Users/nicholasminor/Documents/informatics/E484T_paper/visuals/neutralization_assay.pdf", 
    width = 8, height = 7)
par(bty="n") # deleting the box
gap.plot(c(1:4), antibody$Chronic, 
         gap=c(1.5,9.5), 
         ylim = c(0,10), ytics = c(0.0, 0.5, 1.0, 10),
         xtics = c(1,2,3,4), xticlab = row.names(antibody),
         xlab = NA, ylab = "IC99 values µg/ml",
         cex.axis = 0.8,
         col = palette$points[2], pch = 2, lwd = 2, cex = 1.2)
grid()
# lines(c(1:4), antibody$S.614G, col = palette$lines[1], lwd = 2)
points(c(1:4), antibody$S.614G, col = palette$points[1], pch = 0, cex = 1.2)

legend(3, 2,
       legend = c("S-614G Control Virus", "Chronic Infection Virus"),
       pch = c(0, 2),
       col = c(palette$points[1], palette$points[2]),
       bty = "n", cex = 0.8, bg="transparent", xpd = T)

# abline(h=seq(1.45,1.575,.001), col="white")
axis.break(2,1.5,style="slash") 

# par(new=T)
# plot(c(1,2),c(10,6.25),
#       ylim = c(0,10), xlim = c(1,4), xaxt="n", yaxt = "n", xlab = NA, ylab = NA,
#       type="l", col = palette$points[2], lwd = 2)
# abline(h=seq(7.4, 7.65,.001), col="white")
axis.break(2,7.5,style="slash")
dev.off()

