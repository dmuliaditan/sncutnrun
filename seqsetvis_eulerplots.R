#This code is to extract Signac peak files and analyze through seqsetvis and to create euler plots

#Load necessary packages
library(eulerr)
library(seqsetvis)
library(GenomicRanges)
library(data.table)
library(cowplot)
library(bedr)
theme_set(cowplot::theme_cowplot())
library(enrichTF)
require(rtracklayer)
library(RColorBrewer)
library(VennDiagram)
library(BSgenome.Hsapiens.UCSC.hg38)


#Import data
#Extract H3K4me3 peaks from H3K4me3 object
load("K4me3_final.RData")
sample_names <- c("HN120PRI", "HN120MET", "HN120PCR", "HN137PRI", "HN137MET", "HN137PCR")
for (b in seq_along(sample_names)) {
  print(paste(sample_names[b]))
  peaks_index <- which(grepl(peaks$peak_called_in,pattern = sample_names[b]) == T)
  data <- data.frame(chr=peaks@seqnames[peaks_index],peaks@ranges[peaks_index])
  write.table(x = data,file = paste0(sample_names[b],"_K4me3_peaks_Signac.txt"),
              quote = F, sep = "\t", col.names = F, row.names = F)
}

#Extract H3K27ac peaks from H3K27ac object
load("K27ac_final.RData")
sample_names <- c("HN120PRI", "HN120MET", "HN120PCR", "HN137PRI", "HN137MET", "HN137PCR")
for (b in seq_along(sample_names)) {
  print(paste(sample_names[b]))
  peaks_index <- which(grepl(peaks$peak_called_in,pattern = sample_names[b]) == T)
  data <- data.frame(chr=peaks@seqnames[peaks_index],peaks@ranges[peaks_index])
  write.table(x = data,file = paste0(sample_names[b],"_K27ac_peaks_Signac.txt"),
              quote = F, sep = "\t", col.names = F, row.names = F)
}

sample_list <- c('HN120PRI_K4me3',
                 'HN120PRI_K27ac',
                 'HN120MET_K4me3',
                 'HN120MET_K27ac',
                 'HN120PCR_K4me3',
                 'HN120PCR_K27ac',
                 'HN137PRI_K4me3',
                 'HN137PRI_K27ac',
                 'HN137MET_K4me3',
                 'HN137MET_K27ac',
                 'HN137PCR_K4me3',
                 'HN137PCR_K27ac')

bedlist <- NULL
for (a in seq_along(sample_list)) {
  print(paste(sample_list[a]))
  bed_granges <- bed_to_granges(file = paste0(sample_list[a],"_peaks_Signac.txt"))
  bedlist[[a]] <- bed_granges
  bed_granges <- as.data.frame(bed_granges)
  print(nrow(bed_granges))
}
names(bedlist) <- sub(pattern = "_", replacement = " ", x = sample_list, fixed = T)

olaps <- ssvOverlapIntervalSets(bedlist, maxgap = 10)
ssvFeatureBars(olaps) +
  ylab("") +
  scale_x_discrete(labels = names(bedlist)) +
  ggtitle("Number of Signac Peaks") +
  scale_fill_brewer(labels = names(bedlist), palette = "Paired") +
  theme(legend.title = element_text(size = 0),
        axis.text.x = element_text(size = 18, angle = 45, vjust = 1),
        legend.text = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        title = element_text(size = 24))
