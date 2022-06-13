################################################
############ snCUT&RUN data analysis ###########
############# v. 13 June 2022 ##################

######## All code were processed with ##########
############## R version 4.1 ###################
########## Author: Daniel Muliaditan ###########

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

#For Signac analysis, two files are needed: 
#1) A fragment file, which is basically a bed file with the following columns:
#   chr, start, end, barcode, number of PCR duplicates of that fragment
#2) Barcode file with single-cell metadata
#See section 1 (Data Preprocessing for the code to create these files)

#Fetch the fragment file and create the bin matrix
#Define the path of the fragment file
setwd("/mnt/raid5/cutnrun/userFiles/daniel/DM112021/K27ac")

#For H3K27ac or H3K4me3
#Also import single-cell metadata  and genome
fragment.path <- "/mnt/raid5/cutnrun/signac/K27ac_combined_sorted.bed.gz"
barcodes <- read.table(file = "/mnt/raid5/cutnrun/signac/K27ac_barcodes_new.txt", header = T)
colnames(barcodes) <- c("Barcode", "UMRs", "Date", "Experiment", "Reads_in_Peaks", 
                        "Reads_in_Blacklist", "Pct_reads_in_peaks", "Pct_reads_in_blacklist")
histone <- "H3K27ac" #Set which histone modification is being analysed
genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38) #Specify genome build (here hg38 was used)

#Make the fragment object from the fragment file
fragments <- CreateFragmentObject(fragment.path)

#Create a sparse count matrix (bin matrix). Rows will be single-cells, columns will be 10kbins
#It takes a long time: >30min. Parallelization is necessary!
#Once done, save it as .rds file to prevent having to create bin matrix again 
bin_matrix <- GenomeBinMatrix( fragments = fragments, genome, binsize = 10000, process_n = 13)

#Save the above results after the first generation 
saveRDS(bin_matrix, file = paste0(histone,"_bin_matrix.rds"))

#In all subsequent runs we can just load it from saved data
bin_matrix <- readRDS(paste0(histone,"_bin_matrix.rds"))

#Make chromatin assay
chrom_assay <- CreateChromatinAssay(
  counts = bin_matrix,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = fragment.path,
  min.cells = 0,
  min.features = 200
)

#Process the metadata for the cells, filter out cells with <200 fragments
low_umr_ind <- which(barcodes$UMRs <= 201)
barcodes <- barcodes[-low_umr_ind,]
barcodes <- droplevels.data.frame(x = barcodes)

#Create Seurat object
HN <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "bin",
  meta.data = barcodes
)

dim(HN@meta.data)
HN@meta.data$Barcode <- row.names(HN@meta.data)
HN <- subset(
  x = HN,
  subset = Barcode %in% barcodes$Barcode
)
HN
HN@meta.data <- HN@meta.data[,-c(4:11)] #This step is necessary due to some errors when making the metadata of the assay
order <- barcodes[match(rownames(HN@meta.data), barcodes$Barcode),]

#Check if barcode order of "order" is equal to "HN@meta.data"
all.equal(target = rownames(HN@meta.data), current = order$Barcode)

#Append the barcodes data
HN@meta.data <- cbind(HN@meta.data, order[,c(2:8)])
HN
HN[['bin']]
granges(HN)
HN@meta.data$Experiment <- factor(HN@meta.data$Experiment, levels = c("HN120PRI", "HN120MET", "HN120PCR",
                                                                      "HN137PRI", "HN137MET", "HN137PCR"))
HN@meta.data$Type <- ifelse(test = grepl(pattern = "PRI", x = HN@meta.data$Experiment),
                            yes = "Primary", no = ifelse(test = grepl(pattern = "MET", x = HN@meta.data$Experiment),
                                                         yes = "Metastatic", no = "CisplatinResistant"))
HN@meta.data$Type <- factor(HN@meta.data$Type, levels = c("Primary", "Metastatic", "CisplatinResistant"))

#Extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#Change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
Annotation(HN) <- annotations #Add the gene information to the object

#QC filtering:
#QC1: Compute nucleosome signal score per cell
HN <- NucleosomeSignal(object = HN)
HN$nucleosome_group <- ifelse(HN$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = HN, group.by = 'nucleosome_group')

#QC2: Compute TSS enrichment score per cell
HN <- TSSEnrichment(object = HN, fast = FALSE)
plot(HN$TSS.enrichment)
HN$high.tss <- ifelse(HN$TSS.enrichment > 0.3, 'High', 'Low')
TSSPlot(HN, group.by = 'high.tss') + NoLegend() + ylim(c(0,5))

#QC3: UMRs <1000 and >100000
#QC4: Percent reads in peaks
#QC5: Percent reads in blacklist

#Reads in peaks and blacklist should be in percentages
HN@meta.data$Pct_reads_in_peaks <- HN@meta.data$Pct_reads_in_peaks*100
HN@meta.data$Pct_reads_in_blacklist <- HN@meta.data$Pct_reads_in_blacklist*100

#Plot QC figures
VlnPlot(object = HN, features = c('Pct_reads_in_peaks', 
                                  'UMRs', 
                                  'TSS.enrichment', 
                                  'Pct_reads_in_blacklist', 
                                  'nucleosome_signal'),
        pt.size = 0.1, ncol = 5)

#At this point save the preprocessed Signac object, before any subsetting
saveRDS(HN, file=paste0(histone,"_orig.rds"))
