#In this section, we will subset the dataset by filtering low quality cells, cluster the cells by performing dimensional reduction with UMAP
#and look at batch effects and QC measures per cell

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

#Load the dataset
HN <- readRDS(paste0(histone,"_orig.rds")) #In all subsequent runs we can just load it from saved data
HN

#Filter cells:
HN <- subset(
  x = HN,
  subset = 
    UMRs > 1000 &
    UMRs < 100000 &
    Pct_reads_in_peaks > 25 &
    Pct_reads_in_blacklist < 1 &
    nucleosome_signal < 5 &
    TSS.enrichment > 0.5
)
HN

#Perform TFIDF normalization and singular value decomposition
HN <- RunTFIDF(HN, verbose = T)
HN <- FindTopFeatures(HN, min.cutoff = 'q5')
HN <- RunSVD(HN, verbose = T)
DepthCor(HN)#Remove dim 1 if component has close to -1 correlation

#Non-linear dimension reduction and visualization
HN <- RunUMAP(object = HN, reduction.key = "UMAP_", 
              umap.method = "uwot",
              n.neighbors = 10, 
              metric = "manhattan", 
              n.epochs = 500,
              min.dist = 0.1,
              spread = 1,
              set.op.mix.ratio = 1,
              reduction = 'lsi', 
              dims = 1:29)
HN <- FindNeighbors(object = HN, 
                    distance.matrix = F,
                    k =17,
                    annoy.metric = "manhattan", 
                    reduction = 'lsi', dims = 2:30,
                    verbose = T)
HN <- FindClusters(object = HN, 
                   modularity.fxn = 1,
                   resolution = 0.5,
                   algorithm = 2,
                   n.start = 100,
                   n.iter = 100,
                   verbose = T)

DimPlot(object = HN, label = TRUE, pt.size = 2) + NoLegend()

#UMAP by Experiment
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
DimPlot(object = HN, 
        group.by="subclust", 
        cols = cbPalette, 
        shape.by = "Type",
        label = F, 
        pt.size = 2
)  +
  ggtitle(histone) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(title = element_text(size = 30),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 24),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  guides(color = guide_legend(override.aes = list(size=10)))


#UMAP by batch
for (p in unique(HN@meta.data$Experiment)) {
  print(paste(p))
  subset <- HN@meta.data[which(HN@meta.data$Experiment == p),]
  
  for (q in seq_along(unique(subset$Date))) {
    print(paste(q))
    index <- which(subset$Date == unique(subset$Date)[q])
    HN@meta.data$Date[which(HN@meta.data$Date == unique(subset$Date)[q])] <- paste0(p, " ", "Rep", " ", q)
  }
}

DimPlot(object = HN, 
        group.by="Date", 
        label = F, 
        pt.size = 2) +
  ggtitle(histone) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(title = element_text(size = 30),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) #Figure 2

#Change default identities from Seurat cluster number to cell type
Idents(HN) <- HN@meta.data$Experiment
DefaultAssay(HN) <- 'bin'

#Plot QC measurements per cell
FeaturePlot(
  object = HN,
  features = c('Pct_reads_in_peaks'),
  pt.size = 2,
  cols = c("blue", "red")
  
) +
  ggtitle("Percent reads in Peaks") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        title = element_text(size = 24),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))


#QC per cluster
VlnPlot(object = HN, 
        features = c('Pct_reads_in_peaks', 'UMRs', 'TSS.enrichment', 
                     'Pct_reads_in_blacklist', 'nucleosome_signal'),
        pt.size = 0.1, ncol = 5)

save.image(file = paste0(histone,"_postclustering.RData"))
