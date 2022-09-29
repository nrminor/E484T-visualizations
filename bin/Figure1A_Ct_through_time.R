#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)



### PLOTTING QRT-PCR CT VALUES OVER THE COURSE OF A 486-DAY PROLONGED INFECTION
# UPDATED: 29-Mar-2022 by Nicholas R. Minor
# ----------------------------------------------------------------- #

# This script will simply import and plot the available Ct values. 

# Please direct all questions about the reproducibility of this script to:
# nrminor@wisc.edu

# ----------------------------------------------------------------- #



# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("The working directory must be specified here.", call.=FALSE)
} else if (length(args)>1) {
  stop("This R script only accepts one argument.", call.=FALSE)
} 



### PREPARING THE ENVIRONMENT AND READING IN DATA ####
### --------------------------------------------- #
list.of.packages <- "ggplot2"
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
invisible(lapply(list.of.packages, library, character.only = TRUE))

ct_values = args[1]

Ct <- read.csv(ct_values)
Ct$SARS.CoV.2.PCR.Ct <- round(Ct$SARS.CoV.2.PCR.Ct, digits = 1)
Ct <- Ct[is.na(Ct$SARS.CoV.2.PCR.Ct)==F,]



### PLOTTING DATA ####
### ------------- #
pdf(file = "fig1a_ct_values_thru_time.pdf", 
    width = 8, height = 6)

Ct_plot <- ggplot(data = Ct, aes(x = Encounter.days.post.first.SARS.CoV.2.test, 
                                  y = SARS.CoV.2.PCR.Ct)) +
  coord_cartesian(ylim = c(30,10)) + 
  scale_x_continuous(n.breaks = 5) +
  geom_segment(aes(x = Encounter.days.post.first.SARS.CoV.2.test[6],
                   xend=Encounter.days.post.first.SARS.CoV.2.test[6],
                   y=9,yend=31), 
               colour = "red", size=0.5) +
  geom_segment(aes(x = 432,xend=432,
                   y=9,yend=31), 
               colour = "black", size=0.5, linetype = "dashed") +
  geom_point(size = 4) + 
  # geom_label(aes(x=Encounter.days.post.first.SARS.CoV.2.test, 
  #                y = SARS.CoV.2.PCR.Ct, label = SARS.CoV.2.PCR.Ct), 
  #            nudge_y = c(1.2, 1.2, -1.2, 1.2, -1.2, 1.2, -1.2, 1.2, 1.2, 1.2, -1.2, 1.2)) + 
  theme(panel.background = element_rect(fill = "#FFFFFF"), 
        panel.grid = element_line(color = "#D3D3D3"),
        axis.line = element_line(color = "black"), 
        legend.key.height = unit(1, 'cm')) + 
  labs(x = "Post-Diagnosis Day",
       y = "PCR Cycle Threshold (Ct)",
       col="") + 
  theme(legend.position = "none") + 
  theme(plot.margin = unit(par("mar"), "line")) +
  theme(text = element_text(size = 14))
Ct_plot

dev.off()

