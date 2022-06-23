#Check that UMR and FRiP are not confounding factors in preMET analysis
### HN120preMET UMR count
my_comparisons_HN120 <- list( c("HN120PRI", "HN120MET"), c("HN120PRI", "HN120preMET"), c("HN120MET", "HN120preMET") )
ggplot(subset(HN@meta.data, subclust %in% c("HN120PRI", "HN120MET", "HN120preMET")), 
       mapping = aes(x = subclust, y = UMRs, fill = subclust))+
  geom_violin(lwd=1) +
  geom_boxplot(width=0.3, lwd=1, fill = "white") +
  scale_fill_brewer(palette = "Accent") +
  theme_bw() +
  ylim(c(0,130000)) +
  stat_compare_means(comparisons = my_comparisons_HN120, size = 8, label = "p.signif", vjust = 0.1, label.y = c(100000, 110000, 122000)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=36),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5  ),
        axis.line = element_blank())

### HN120preMET FRiP count
ggplot(subset(HN@meta.data, subclust %in% c("HN120PRI", "HN120MET", "HN120preMET")), 
       mapping = aes(x = subclust, y = Pct_reads_in_peaks, fill = subclust))+
  geom_violin(lwd=1) +
  geom_boxplot(width=0.3, lwd=1,fill = "white") +
  scale_fill_brewer(palette = "Accent") +
  theme_bw() +
  ylim(c(0,100)) +
  ylab("Percentage of reads in peaks\n") +
  stat_compare_means(comparisons = my_comparisons_HN120, size = 8, label = "p.signif", label.y = c(73,80,87)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5),
        axis.line = element_blank())

#Linear regression to show that FRIP are not correlated with module scores
ggplot(subset(HN@meta.data, subclust %in% c("HN120preMET")), 
       mapping = aes(x = Pct_reads_in_peaks, y = HN120MET))+
  geom_point() +
  geom_smooth(method = "lm") +
  stat_regline_equation(label.y = 9, aes(label = ..eq.label..), size = 7) +
  stat_regline_equation(label.y = 8.5, aes(label = ..rr.label..), size = 7) +
  theme_bw() +
  ylab("HN120PRI Module Score") +
  xlab("HN120preMET - Pct reads in Peaks") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size=28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5),
        axis.line = element_blank()) 

#Linear regression to show that UMRs are not correlated with module scores
ggplot(subset(HN@meta.data, subclust %in% c("HN120preMET")), 
       mapping = aes(x = UMRs, y = HN120MET))+
  geom_point() +
  geom_smooth(method = "lm") +
  stat_regline_equation(label.y = 9, aes(label = ..eq.label..), size = 7) +
  stat_regline_equation(label.y = 8.5, aes(label = ..rr.label..), size = 7) +
  theme_bw() +
  ylab("HN120PRI Module Score") +
  xlab("HN120preMET - UMRs") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size=28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5),
        axis.line = element_blank()) 
