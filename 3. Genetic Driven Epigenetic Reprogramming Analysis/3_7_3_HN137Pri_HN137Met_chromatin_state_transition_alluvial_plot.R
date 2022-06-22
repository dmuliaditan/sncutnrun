#3.7.3 HN137Pri > HN137Met chromatin state alluvial plot

#After correlating RNAseq Z-score, CNV, and chromatin-state statically in HN137Met, we look at the transition between HN137Pri and HN137Met
#First, we look at the chromatin state changes using an alluvial plot

#Version: 22/06/2022
#Author: Daniel Muliaditan

#Load existing dataset and required packages
setwd("D:/snCUT_RUN/scripts")

library(networkD3)
library(dplyr)
library(tidyr)
library(ggalluvial)
library(rGREAT)

#Plot alluvial plot
#With the results from the table, make an excel sheet with the transitions
#Repeat for HN120Pri-HN120PCR, HN137Pri-HN137Met, HN137Pri-HN137PCR

table <- as.data.frame(table)
links <- separate(data = table, col = Var1, into = c("source", "target"), sep = "\\-")
colnames(links) <- c("source", "target", "value")
links$source <- as.factor(links$source)
links$target <- as.factor(links$target)
links$value <- as.numeric(links$value)
links <- links[-19,] #Remove E4-E4 transition
links

is_alluvia_form(as.data.frame(links), axes = 1:3, silent = TRUE)
ggplot(links, aes(y = value, axis1 = source, axis2 = target)) +
  geom_alluvium(aes(fill = source),width = 1/12) +
  geom_stratum(width = 1/12, color = "black", size = 1) +
  #geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("source", "target"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("HN137Pri - HN137Met Chromatin State Transition") +
  ylab("") +
  theme(title = element_text(size = 22),
        legend.title = element_blank(),
        legend.text = element_text(size = 24),
        axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()) #Plot was further refined in powerpoint

#Looking at pathways, for example the regions gaining H3K27ac from H3K4me3 alone
chrom_subset <- which(segment_combined$state == "E1-E2")
subseg <- segment_combined[chrom_subset,]

#rGREAT analysis
job <- submitGreatJob(gr = subseg, species = 'hg38', rule = 'twoClosest',
                      adv_twoDistance = 2000)
tb <- getEnrichmentTables(job, category = "GO")
names(tb)
head(tb$`GO Biological Process`, n=10) #Look at affected pathways
availableOntologies(job)
res <- plotRegionGeneAssociationGraphs(job) #Most regions gaining H3K27ac from H3K4me3 alone seems to be concentrated near the promoter, which is expected
