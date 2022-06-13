#In this section, we will call peaks from the dataset and plot coverage of the genomic regions of interest

#Version: 13/06/2022
#Author: Daniel Muliaditan

#Load necessary packages and set seed
library(ggplot2)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(tidyverse)
library(chromVAR)
library(SummarizedExperiment)
library(dplyr)
library(ggpubr)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(SeuratWrappers)
library(Matrix)
library(monocle3)
library(rGREAT)
set.seed(1234)

#Set the working directory and load the dataset
setwd("/mnt/raid5/cutnrun/userFiles/daniel/DM112021/K27ac")
load(paste0(histone,"_postclustering.RData")) 

#Call peaks with MACS2
peaks <- CallPeaks(
  object = HN,
  assay = "bin",
  group.by = "Experiment",
  broad = F,
  macs2.path = "/home/cutnrun/miniconda3/bin/macs2",
  outdir = getwd(),
  combine.peaks = T,
  additional.args = "--nomodel -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits -B"
)

#Plot coverage gene(s) of interest, below genes are example genes
coverage_genes <- c("ARID5B", "BASP1", "EGFR", "FGFR3", "CDR2L", "RAI14",
                    "NEK2", "AURKA", "TGFBR2", "HDAC9", "BRD4", "TGFBR1", 
                    "BRD2", "GMNN", "ASXL1", "HDAC9", "HDAC1", 
                    "ATM", "ANKRD10", "PDLIM1", "RIF1")

for (j in c(coverage_genes)) {
  print(paste(j))
  png(filename = paste0("coverageplot_",j,".png"), width = 1440, height = 720)
  g <- CoveragePlot(
    object = HN,
    assay = "bin",
    region = j,
    ranges = peaks,
    ranges.title = "MACS2",
    annotation = T,
    peaks = T,
    extend.upstream = 20000,
    extend.downstream = 20000
    #, idents = c("HN120PRI", "HN120MET", "HN120PCR") #Unhash to look at specific cell lines
  )
  print(g)
  dev.off()
}

save.image(file = paste0(histone,"_postPeakCalling.RData")) #Checkpoint
