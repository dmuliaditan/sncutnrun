#In this section, we perform differential peak calling and differential motif analysis

#Version: 12/07/2022
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
setwd("D:/snCUT_RUN/signac/DM_052022")
load(paste0(histone,"_tfmotifs.RData")) 

#Motif analysis: everything at the same time
#Make a list of the comparisons to be made
comparisons <- list()
comparisons[[1]] <- c("HN120PRI", "HN120MET")
comparisons[[2]] <- c("HN120PRI", "HN120PCR")
comparisons[[3]] <- c("HN120MET", "HN120PRI")
comparisons[[4]] <- c("HN120PCR", "HN120PRI")
comparisons[[5]] <- c("HN137PRI", "HN137MET")
comparisons[[6]] <- c("HN137PRI", "HN137PCR")
comparisons[[7]] <- c("HN137MET", "HN137PRI")
comparisons[[8]] <- c("HN137PCR", "HN137PRI")

Idents(HN) <- HN@meta.data$Experiment
Idents(HN) <- factor(Idents(HN), 
                     levels = c("HN120PRI", "HN120MET", "HN120PCR", "HN137PRI", "HN137MET", "HN137PCR"))

for (o in 1:8) {
  print(paste(o))
  da_peaks <- FindMarkers(
    object = HN,
    ident.1 = comparisons[[o]][[1]],
    ident.2 = comparisons[[o]][[2]],
    only.pos = TRUE,
    test.use = 'LR',
    latent.vars = 'nCount_peak',
    min.pct = 0.05,
    assay = 'peak'
  )
  
  #Get top differentially accessible peaks
  top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
  
  #Get all peaks of HN137Pri and HN137Met
  cells <- HN$Experiment == comparisons[[o]][[1]] | HN$Experiment == comparisons[[o]][[2]]  # choose cells for the cluster you want
  counts <- GetAssayData(HN, assay = 'peak', slot = 'counts')
  peaks.use <- which(Matrix::rowSums(counts[, cells]) > 0)
  
  peaks.matched <- MatchRegionStats(
    meta.feature=meta.feature[names(peaks.use),],
    query.feature=meta.feature[top.da.peak,],
    n=nrow(meta.feature[names(peaks.use),]),
    verbose = T
  )
  
  #Test enrichment
  enriched.motifs <- FindMotifs(
    object = HN,
    features = top.da.peak,
    assay = 'peak',
    background=peaks.matched
  )
  
  write.table(x = enriched.motifs[,c(1,6:8)], 
              file = paste0(histone,"_enriched_motifs_",comparisons[[o]][[1]], "_vs_",comparisons[[o]][[2]],".txt"),
              sep = "\t", 
              quote = F, 
              row.names = F, 
              col.names = T)
  
  #Save peaks as bed file to do GREAT analysis
  startend_peaks <- sub("^[^-]*-", "", top.da.peak)
  great_peaks <- data.frame(chr=sub("-.*", "", top.da.peak),
                            start=as.numeric(sub("-.*", "", startend_peaks)),
                            end= as.numeric(sub(".*-", "", top.da.peak)))
  write.table(x = great_peaks, 
              file = paste0(histone,"_peaks_",comparisons[[o]][[1]], "_vs_",comparisons[[o]][[2]],".bed"),
              col.names = F, 
              row.names = F, 
              sep = "\t", 
              quote = F)
  
}

#Specific comparisons
da_peaks <- FindMarkers(
  object = HN,
  ident.1 = 'HN137MET',
  ident.2 = 'HN137PRI',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peak',
  min.pct = 0.05,
  assay = 'peak'
)

#Get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

#Get all peaks of HN137Pri and HN137Met
cells <- HN$Experiment == "HN137PRI" | HN$Experiment == "HN137MET"  # choose cells for the cluster you want
counts <- GetAssayData(HN, assay = 'peak', slot = 'counts')
peaks.use <- which(Matrix::rowSums(counts[, cells]) > 0)

peaks.matched <- MatchRegionStats(
  meta.feature=meta.feature[names(peaks.use),],
  query.feature=meta.feature[top.da.peak,],
  n=nrow(meta.feature[names(peaks.use),]),
  verbose = T
)

#Test enrichment
enriched.motifs <- FindMotifs(
  object = HN,
  features = top.da.peak,
  assay = 'peak',
  background=peaks.matched
)
enriched.motifs[1:100,]

HN120PRI_motifs <- enriched.motifs[1:100,]
HN120MET_motifs <- enriched.motifs[1:100,]
HN120preMET_motifs <- enriched.motifs[1:100,]

#Save motif images and save ChromVAR images of motif score per cell
enriched.motifs[enriched.motifs$motif.name == "FOXC1",]
pattern <- ":|-|\\(|\\."

for (k in which(grepl(pattern,setdiff(bedlist[[1]],bedlist[[2]])) == F)){
  print(paste(k))
  png(filename = paste0(enriched.motifs$motif.name[k],"_motif.png"), width = 1000, height = 150, units = "px")
  l = MotifPlot(
    object = HN,
    motifs = rownames(enriched.motifs)[enriched.motifs$motif.name == "KLF11"],
    assay = "peak", facet = "wrap"
  ) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          title = element_blank(),
          strip.text = element_blank()
    )
  print(l)
  dev.off()
  
  png(filename = paste0(enriched.motifs$motif.name[k],"_motif_ChromVAR.png"), width = 1000, height = 1000, units = "px")
  motif_ind <- which(enriched.motifs$motif.name == enriched.motifs$motif.name[k])
  p = FeaturePlot(
    object = HN,
    features = enriched.motifs$motif[motif_ind],
    min.cutoff = 'q10',
    max.cutoff = 'q90',
    pt.size = 3
  ) +
    ggtitle(enriched.motifs$motif.name[motif_ind]) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(axis.title.x = element_text(size = 24),
          axis.title.y = element_text(size = 24),
          axis.text = element_text(size = 18),
          legend.text = element_text(size = 18),
          title = element_text(size=30),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  print(p)
  dev.off()
}

#Find specific motif

for (t in c("KLF5", "TEAD4", "FOSL2")) {
  print(t)
  TF <- t
  motif_ind <- which(enriched.motifs$motif.name == TF)
  motifplot <- FeaturePlot(
    object = HN,
    features = enriched.motifs$motif[motif_ind],
    min.cutoff = 'q10',
    max.cutoff = 'q90',
    pt.size = 1,
    cells = rownames(HN@meta.data)[which(grepl(pattern = "HN137", x = HN$Experiment) == T)]
  ) +
    ggtitle(enriched.motifs$motif.name[motif_ind]) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text = element_text(size = 9),
          legend.text = element_text(size = 9),
          title = element_text(size=13),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5))
  ggsave(filename = paste0("D:/snCUT_RUN/figures/figures/supplementary_figs/", TF,"_HN137_motifs.pdf"),
         plot = motifplot, width = 1024, height = 1024, units = "px", dpi = 300)
}

