### PLOTTING QRT-PCR CT VALUES OVER THE COURSE OF A 433-DAY PROLONGED INFECTION
# UPDATED: 04-Jan-2022 by Nicholas R. Minor
# ----------------------------------------------------------------- #

# This script will simply import and plot the available Ct values. 

# Please direct all questions about the reproducibility of this script to:
# nrminor@wisc.edu

# ----------------------------------------------------------------- #



### PREPARING THE ENVIRONMENT AND READING IN DATA ####
### --------------------------------------------- #
library(ggplot2)

data_filepath = "/Volumes/working_ssd/e484t_manuscript/"
setwd(data_filepath)
Ct <- read.csv("data/Ct_timeline.csv")
Ct$Encounter.date <- as.Date(Ct$Encounter.date, format = "%m/%d/%y")
Ct$SARS.CoV.2.PCR.Ct <- round(Ct$SARS.CoV.2.PCR.Ct, digits = 1)



### PLOTTING DATA ####
### ------------- #
pdf(file = "/Users/nicholasminor/Documents/informatics/E484T_paper/visuals/Ct_values.pdf", 
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
  geom_label(aes(x=Encounter.days.post.first.SARS.CoV.2.test, 
                 y = SARS.CoV.2.PCR.Ct, label = SARS.CoV.2.PCR.Ct), 
             nudge_y = c(1.2, 1.2, -1.2, 1.2, -1.2, 1.2, -1.2, 1.2, 1.2, 1.2, -1.2, 1.2)) + 
  theme(panel.background = element_rect(fill = "#FFFFFF"), 
        panel.grid = element_line(color = "#D3D3D3"),
        axis.line = element_line(color = "black"), 
        legend.key.height = unit(1, 'cm')) + 
  labs(x = "Day of Infection",
       y = "PCR Cycle Threshold (Ct)",
       col="") + 
  theme(legend.position = "none") + 
  theme(plot.margin = unit(par("mar"), "line")) +
  theme(text = element_text(size = 14))
Ct_plot

dev.off()

