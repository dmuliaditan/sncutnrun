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

#Plot number of Signac peaks
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

ssvFeatureBinaryHeatmap(olaps)
ssvFeatureVenn(object = olaps[,c(1,3)], line_width = 1, line_alpha = 0.5) +
  theme(legend.text = element_text(size=14))
ssvFeatureUpset(olaps)

#Venn diagrams with venneuler
euler <- cbind(
  "HN120PRI K4me3" = olaps$`HN120PRI K4me3`,
  "HN120MET K4me3"= olaps$`HN120MET K4me3`
  )

plot(euler(euler), 
     quantities = list(cex = 4),
     labels = F,
     fill = c("#EA15D3", "#D3EA15"),
     edges =  F,
     alpha = 1,
     ) # "#D3EA15" for Met, "#15D3EA" for PCR
