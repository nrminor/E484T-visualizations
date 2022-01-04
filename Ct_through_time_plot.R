library(ggplot2)
library(cowplot)

data_filepath = "/Volumes/working_ssd/e484t_manuscript/"
setwd(data_filepath)
Ct <- read.csv("data/Ct_timeline.csv")
Ct$Encounter.date <- as.Date(Ct$Encounter.date, format = "%m/%d/%y")
Ct$SARS.CoV.2.PCR.Ct <- round(Ct$SARS.CoV.2.PCR.Ct, digits = 1)

# ppi <- 300
# png("visuals/Ct_values_treatments.png", width = 8*ppi, height = 8*ppi, res = ppi)
# svg("visuals/Ct_values_treatments.svg")
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
  # geom_text(aes(x=Encounter.days.post.first.SARS.CoV.2.test[6], 
  #               y= 9, label = "Bamlanivumab administered"), 
  #           color = "black", vjust = 1, size = 4) +
  # geom_line(linetype = 2) + 
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
  # scale_color_gradient(name = "Relative\nViral Load",
  #                      breaks = c(0.04, 0.06, 0.08),
  #                      labels = c("Low","Medium","High"))
Ct_plot

# treatments <- ggplot(data = Ct, aes(x = Encounter.days.post.first.SARS.CoV.2.test, 
#                                     y = SARS.CoV.2.PCR.Ct)) +
#   scale_x_continuous(n.breaks = 10, position = "top") + 
#   theme(axis.title.x=element_blank(),
#         axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
#         axis.line.y = element_blank(), axis.line.x = element_line(color = "black"),
#         panel.background = element_rect(fill = "#FFFFFF"))
# 
# plot_grid(Ct_plot, treatments, align = "v", ncol = 1, rel_heights = c(.7, .3), axis = "tblr")
# # ggsave("visuals/Ct_values_treatments.svg", device = svg, dpi = "print") # a "print" dpi is 300, the same as what we set ppi to above 
# ggsave("/Users/nicholasminor/Documents/informatics/E484T_paper/visuals/Ct_values_treatments.pdf", 
#        device = pdf, dpi = "print",
#        height = 6,
#        width = 8) # a "print" dpi is 300, the same as what we set ppi to above 
dev.off()

