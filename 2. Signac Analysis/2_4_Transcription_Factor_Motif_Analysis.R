#In this section, we perform transcription factor motif analysis on the peak calls from step 2_3

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
load(paste0(histone,"_postPeakCalling.RData")) 

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

#Create peak assay, this assay will be used to call the TF motifs
#Create peak feature matrix
peak_matrix <- FeatureMatrix(
  fragments = Fragments(HN),
  features = peaks
)

#Add the peak matrix to the peak assay
HN[["peak"]] <- CreateChromatinAssay(
  counts = peak_matrix[,colnames(HN)],
  sep = c(":", "-"),
  fragments = fragment.path,
)

# Add motif information to the peak assay, this step will take a while!
HN <- AddMotifs(
  object = HN,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm,
  assay = 'peak'
)

#Get meta features, these are features such as GC content that will be used for the background signal later on
meta.feature <- GetAssayData(HN, 
                             assay="peak", 
                             slot="meta.features")

#Run ChromVAR to get the TF motif signal for all cell lines across the dataset
HN <- RunChromVAR(
  object = HN,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = 'peak'
)

#Save ChromVAR results
write.table (HN$chromvar@data, paste0(histone, "_deviation.txt"), 
             row.names=TRUE, 
             col.names=TRUE)

save.image(file = paste0(histone,"_tfmotifs.RData"))
