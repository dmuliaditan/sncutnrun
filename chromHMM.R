#R script to perform ChromHMM alluvial plots and subsetting genomic regions based on transitions
library(networkD3)
library(dplyr)
library(xlsx)
library(ggalluvial)
library(rGREAT)

#Import Excel sheet - HN137Pri > HN137Met
links <- read.xlsx2(file = "E:/snCUT_RUN/results/ChromHMM/ChromHMM_transitions.xlsx", sheetIndex = 3, header = F)[c(1:24),]
colnames(links) <- c("source", "target", "value")
links$source <- as.factor(links$source)
links$target <- as.factor(links$target)
links$value <- as.numeric(links$value)

is_alluvia_form(as.data.frame(links), axes = 1:3, silent = TRUE)
links <- as.data.frame(links)
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
        axis.ticks = element_blank())

head(links)

#Subset genomic regions based on transitions
setwd("E:/snCUT_RUN/results/ChromHMM/outputdir/5_states")
bed1 <- read.table(file = "HN137Pri_5_segments_binned_sorted.bed", header = F)
bed2 <- read.table(file = "HN137Met_5_segments_binned_sorted.bed", header = F)

comparison <- data.frame(bin=paste(bed1$V1,bed1$V2,bed1$V3, sep = "-"),transition = paste0(bed1$V4,"-",bed2$V4))

table <- table(comparison$transition)
table

segment_combined <- data.frame(chr=bed1$V1,
                               start=bed1$V2,
                               end=bed1$V3,
                               state=paste0(bed1$V4,"-",bed2$V4))

setwd("E:/snCUT_RUN/results/ChromHMM")

#H3K4me3 gaining H3K27ac during HN137Pri-HN137Met transition
E1_E2 <- which(segment_combined$state == "E1-E2")
subseg <- segment_combined[E1_E2,]
write.table(x = subseg,
            file = "HN137Pri_HN137Met_E1_E2_regions.bed", 
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

#H3K4me3 losing H3K4me3 but gaining H3K27ac during HN137Pri-HN137Met transition
E1_E3 <- which(segment_combined$state == "E1-E3")
subseg <- segment_combined[E1_E3,]
write.table(x = subseg,
            file = "HN137Pri_HN137Met_E1_E3_regions.bed", 
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

#H3K4me3 becoming unmodified during HN137Pri-HN137Met transition
E1_E4 <- which(segment_combined$state == "E1-E4")
subseg <- segment_combined[E1_E4,]
write.table(x = subseg,
            file = "HN137Pri_HN137Met_E1_E4_regions.bed", 
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

#H3K4me3 becoming heterochromatin during HN137Pri-HN137Met transition
E1_E5 <- which(segment_combined$state == "E1-E5")
subseg <- segment_combined[E1_E5,]
write.table(x = subseg,
            file = "HN137Pri_HN137Met_E1_E5_regions.bed", 
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")
            
#H3K27ac becoming unmodified during HN137Pri-HN137Met transition
E3_E4 <- which(segment_combined$state == "E3-E4")
subseg <- segment_combined[E3_E4,]
write.table(x = subseg,
            file = "HN137Pri_HN137Met_E3_E4_regions.bed", 
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

#Unmodified gaining H3K27me3 during HN137Pri-HN137Met transition
E4_E5 <- which(segment_combined$state == "E4-E5")
subseg <- segment_combined[E4_E5,]
write.table(x = subseg,
            file = "HN137Pri_HN137Met_E4_E5_regions.bed", 
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")

#H3K27me3 becoming unmodified during HN137Pri-HN137Met transition
E5_E4 <- which(segment_combined$state == "E5-E4")
subseg <- segment_combined[E5_E4,]
write.table(x = subseg,
            file = "HN137Pri_HN137Met_E5_E4_regions.bed", 
            quote = F,
            col.names = F,
            row.names = F,
            sep = "\t")
