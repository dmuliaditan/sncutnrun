library(ggplot2)
setwd("D:/snCUT_RUN/metadata")

data  <- read.csv("20200114 UMRs FRiP rank.csv", header = T)

a <- ggplot(data = data, mapping = aes(x = factor(antibody), y = UMR, fill = factor(antibody))) +
  geom_violin(lwd=1) +
  geom_boxplot(width=0.3, lwd=1) +
  geom_hline(yintercept = 500, linetype = 'dashed', size=0.5 ) +
  geom_hline(yintercept = 10000, linetype = 'dashed', size=0.75 ) +
  geom_hline(yintercept = 50000, linetype = 'dashed', size=0.5 ) +
  scale_y_continuous(trans = 'log10', labels = scales::comma,position = "right",
                     breaks = c(1,100,1000,10000,100000)) +
  scale_fill_brewer(palette = "Accent") +
  theme_classic() +
  #coord_flip() +
  ylab("Uniquely Mapped Reads\n") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=36),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5),
        axis.line = element_blank())
  
  
a

c <- ggplot(data = data, mapping = aes(x = factor(antibody), y = Fraction_tA_iPeaks*100, fill = factor(antibody))) +
  geom_violin(lwd=1) +
  geom_boxplot(width=0.3, lwd=1) +
  geom_hline(yintercept = 20, linetype = 'dashed', size=0.5 ) +
  geom_hline(yintercept = 40, linetype = 'dashed', size=0.75 ) +
  geom_hline(yintercept = 60, linetype = 'dashed', size=0.5 ) +
  scale_fill_brewer(palette = "Accent") +
  ylim(c(0,75))+
  scale_y_continuous(position = "right") +
  theme_classic() +
  #coord_flip() +
  ylab("Percentage of reads in peaks\n") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=36),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5),
        axis.line = element_blank()
        )


c
