#4.1. Module Score analysis and HN120preMET / HN137prePCR analysis

#Version: 23/06/2022
#Author: Daniel Muliaditan

#In this section, we will discover epigenetic heterogeneity, and find HN120preMET and HN137prePCR cells

#Load existing dataset and required packages
setwd("D:/snCUT_RUN/scripts")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#AA4499", "#332288",  
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

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
library(rGREAT)
library(ggtern)
set.seed(1234)

#Set the working directory and load the dataset
setwd("D:/snCUT_RUN/signac/DM_052022")
load(paste0(histone,"_tfmotifs.RData")) 

#HN120preMET analysis
#Calculate cell line specific modules
head(HN@meta.data)
#HN@meta.data <- HN@meta.data[,-c(21:27)] #Optional, remove previous module scores

#Rename the identities correctly
Idents(HN) <- HN@meta.data$Experiment
Idents(HN) <- factor(Idents(HN), 
                     levels = c("HN120PRI", "HN120MET", "HN120PCR", "HN137PRI", "HN137MET", "HN137PCR"))

#Find top 50 differential peaks to create a cell line specific 'Module'
HN.markers <- FindAllMarkers(object = HN, 
                             assay = "peak", 
                             only.pos = TRUE, 
                             min.pct = 0.2, 
                             test.use = 'LR',
                             latent.vars = 'nCount_peak')
top50 <- HN.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)

#Plot marker genes in a heatmap
DoHeatmap(object = HN, 
          assay = "peak", 
          features = top50$gene, 
          label = TRUE, 
          slot="data", 
          size = 5, 
          angle = 90) +
  ggplot2::theme(axis.text.y = element_blank(),
                 legend.text = element_text(size = 16),
                 legend.title = element_blank()) 

closest_genes <- ClosestFeature(HN, regions = top50$gene)
top50 %>% print(n = 300)  

#Create modules
HN120PRI_gene<-top50$gene[1:50]  
HN120MET_gene<-top50$gene[51:100]  
HN120PCR_gene<-top50$gene[101:150]
HN137PRI_gene<-top50$gene[151:200]  
HN137MET_gene<-top50$gene[201:250]  
HN137PCR_gene<-top50$gene[251:300]  
modules <- list(HN120PRI_gene, HN120MET_gene, HN120PCR_gene,
                HN137PRI_gene, HN137MET_gene, HN137PCR_gene)
names(modules) <- c("HN120PRI", "HN120MET", "HN120PCR",
                    "HN137PRI", "HN137MET", "HN137PCR")

modules_HN120 <- list(HN120PRI_gene, HN120MET_gene)
names(modules_HN120) <- c("HN120PRI", "HN120MET")

#Save module peaks as bed files
for (k in names(modules)) {
  print(paste(k))
  startend_peaks <- sub("^[^-]*-", "", modules[[k]])
  great_peaks <- data.frame(chr=sub("-.*", "", modules[[k]]),
                            start=as.numeric(sub("-.*", "", startend_peaks)),
                            end= as.numeric(sub(".*-", "", modules[[k]])))
  write.table(x = great_peaks, 
              file = paste0("D:/snCUT_RUN/data/peaks/",k,"_marker_peaks_",histone,".bed"),
              col.names = F, 
              row.names = F, 
              sep = "\t", 
              quote = F)
}

#Calculate module score for each single-cell
HN <- AddChromatinModule(object = HN, 
                         features = modules, 
                         genome = BSgenome.Hsapiens.UCSC.hg38, 
                         assay = 'peak', 
                         verbose = TRUE)

#Plot the cell line specific module scores across all cell lines
VlnPlot(object=HN, features='HN120PCR',
        idents = c("HN120PRI", "HN120MET", "HN120PCR"),
        cols = cbPalette[1:3]) +
  geom_boxplot(width=0.3, lwd=1, fill = "white") +
  coord_cartesian(ylim = c(-5,12), expand = T) +
  xlab("")  +
  theme(title = element_text(size = 36),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 36),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 24),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "none")

#HN120preMET analysis
#Put module scores as a data frame
module.scores <- HN@meta.data[,21:26]
module.scores[is.na(module.scores)] <- 0

#Plot ternary plots
HN120_metadata <- HN@meta.data[which(grepl(pattern = "HN120", HN@meta.data$Experiment) == T),]
HN120_metadata2 <- data.frame(Experiment=HN120_metadata$Experiment,
                              HN120PRI=HN120_metadata$HN120PRI,
                              HN120MET=HN120_metadata$HN120MET,
                              HN120PCR=HN120_metadata$HN120PCR)
rownames(HN120_metadata2) <- rownames(HN120_metadata)
HN120_metadata2 <- na.omit(HN120_metadata2)
range(HN120_metadata2$HN120PRI)
min(HN120_metadata2$HN120PRI)

#Normalize module scores to take into account module scores below zero
pos_normalised_HN120 <- NULL
for (i in seq_along(HN120_metadata2[,1])) {
  print(paste(i))
  row <- HN120_metadata2[i,-1]
  if (min(row) < 0) {
    row=row+abs(min(row))
  }
  pos_normalised_HN120 <- rbind(pos_normalised_HN120,row)
  
}

HN120_metadata2 <- cbind(HN120_metadata2,pos_normalised_HN120)
colnames(HN120_metadata2) <- c("Experiment", "HN120PRI_prenorm", "HN120MET_prenorm", "HN120PCR_prenorm",
                               "HN120PRI_postnorm", "HN120MET_postnorm", "HN120PCR_postnorm")

HN120_metadata3 <- HN120_metadata2[,c(5:7)]
head(HN120_metadata3)

#Process the data so that the total score of each cell is 1
pos_normalised_HN120 <- NULL
for (i in seq_along(HN120_metadata3[,1])) {
  print(paste(i))
  row <- HN120_metadata3[i,]
  rowsum <- sum(row)
  row <- row/rowsum
  pos_normalised_HN120 <- rbind(pos_normalised_HN120,row)
  
}
head(pos_normalised_HN120)
pos_normalised_HN120$Experiment <- HN120_metadata2$Experiment
pos_normalised_HN120$Experiment <- droplevels(pos_normalised_HN120$Experiment)

#Plot ternary plot
ggtern(data = pos_normalised_HN120, aes(HN120MET_postnorm,HN120PRI_postnorm, HN120PCR_postnorm, fill=Experiment, shape=Experiment))  +
  geom_mask() +
  geom_point(size = 3, shape = 21, color="black") + 
  Tlab("HN120Pri") +
  Tarrowlab("HN120PRI module") +
  Llab("HN120Met") +
  Larrowlab("HN120MET module") +
  Rlab("HN120PCR") +
  Rarrowlab("HN120PCR module") +
  scale_fill_manual(values = cbPalette[1:3]) +
  theme_custom(
    base_size = 30,
    base_family = "",
    tern.plot.background = "white",
    tern.panel.background = "white",
    col.T = cbPalette[1],
    col.L = cbPalette[2],
    col.R = cbPalette[3],
    col.grid.minor = "grey"
  ) +
  theme_arrowlarge() +
  theme(legend.position = c(0.3,1),
        legend.justification = c(1, 1),
        tern.axis.arrow = element_line(size=2),
        axis.text.x = element_text(vjust = 1))

#Subset HN120Pri cells with higher HN120Met scores and relabel them as HN120preMET
pos_normalised_HN120$subclust <- 0
pos_normalised_HN120$subclust <- ifelse(test = pos_normalised_HN120$HN120PRI > 0.4 &
                                    pos_normalised_HN120$HN120MET > 0.4 &
                                    pos_normalised_HN120$HN120PCR < 0.2 &
                                    pos_normalised_HN120$Experiment == "HN120PRI", yes = "HN120preMET", no = pos_normalised$Experiment)
table(pos_normalised_HN120$Experiment)
table(pos_normalised_HN120$subclust)
HN120preMET_index <- rownames(pos_normalised)[which(pos_normalised$subclust == "HN120preMET")] #Extract HN120preMET cell barcodes

#Append the meta data
HN@meta.data$subclust <- ifelse(rownames(HN@meta.data) %in% HN120preMET_index, "HN120preMET", levels(HN@meta.data$Experiment)[HN@meta.data$Experiment])
HN@meta.data$subclust <- factor(HN@meta.data$subclust, levels = c("HN120PRI", "HN120MET", "HN120PCR",
                                                                  "HN120preMET", "HN137PRI", 
                                                                  "HN137MET", "HN137PCR"))

#Change default identity to subclust identity
Idents(HN) <- factor(HN@meta.data$subclust, levels = c("HN120PRI", "HN120MET", "HN120PCR",
                                                       "HN120preMET", "HN137PRI", 
                                                       "HN137MET", "HN137PCR"))
pos_normalised[rownames(pos_normalised) %in% rownames(HN@meta.data)[which(HN@meta.data$subclust == "HN120preMET")],]


#Dimplot HN120 with HN120preMET
DimPlot(object = HN, group.by="subclust", label = F, pt.size = 3, 
        cells = rownames(HN@meta.data)[which(HN@meta.data$subclust == "HN120PRI" |
                                               HN@meta.data$subclust == "HN120MET" |
                                               HN@meta.data$subclust == "HN120PCR" |
                                               HN@meta.data$subclust == "HN120preMET")]) +
  scale_colour_manual(values = c("#7F20DF", "#CD3238", "#32CDC7", "#80DF20")) +
  ggtitle(histone) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 24),
        title = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))

#Plot HN120Met module score with HN120preMET cells
VlnPlot(object=HN, features='HN120MET',
        idents = c("HN120PRI", "HN120MET", "HN120preMET")) +
  xlab("")  +
  coord_cartesian(ylim = c(-5, 12), expand = T) +
  ggtitle("HN120MET Module Score") +
  geom_boxplot(width=0.3, lwd=1, fill = "white") +
  scale_fill_manual(values = cbPalette[c(1:2,7)]) +
  theme(title = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))

save.image(paste0("230622_",histone,"_HN120preMET.RData"))


##################################### HN137prePCR analysis #########################################
#Repeat the same thing, but for HN137prePCR
head(HN@meta.data)
HN@meta.data <- HN@meta.data[,-c(21:27)]
Idents(HN) <- HN@meta.data$Experiment
Idents(HN) <- factor(Idents(HN), 
                     levels = c("HN120PRI", "HN120MET", "HN120PCR", "HN137PRI", "HN137MET", "HN137PCR"))

#Find top 50 differential peaks to create a cell line specific 'Module'
HN.markers <- FindAllMarkers(object = HN, 
                             assay = "peak", 
                             only.pos = TRUE, 
                             min.pct = 0.2, 
                             test.use = 'LR',
                             latent.vars = 'nCount_peak')
top50 <- HN.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)

#Plot marker genes
DoHeatmap(object = HN, 
          assay = "peak", 
          features = top50$gene, 
          label = TRUE, 
          slot="data", 
          size = 5, 
          angle = 90) +
  ggplot2::theme(axis.text.y = element_blank(),
                 legend.text = element_text(size = 16),
                 legend.title = element_blank()) 

closest_genes <- ClosestFeature(HN, regions = top50$gene)
top50 %>% print(n = 300)  

#Create modules
HN120PRI_gene<-top50$gene[1:50]  
HN120MET_gene<-top50$gene[51:100]  
HN120PCR_gene<-top50$gene[101:150]
HN137PRI_gene<-top50$gene[151:200]  
HN137MET_gene<-top50$gene[201:250]  
HN137PCR_gene<-top50$gene[251:300]  
modules <- list(HN120PRI_gene, HN120MET_gene, HN120PCR_gene,
                HN137PRI_gene, HN137MET_gene, HN137PCR_gene)
names(modules) <- c("HN120PRI", "HN120MET", "HN120PCR",
                    "HN137PRI", "HN137MET", "HN137PCR")

#Save module peaks as bed files
for (k in names(modules)) {
  print(paste(k))
  startend_peaks <- sub("^[^-]*-", "", modules[[k]])
  great_peaks <- data.frame(chr=sub("-.*", "", modules[[k]]),
                            start=as.numeric(sub("-.*", "", startend_peaks)),
                            end= as.numeric(sub(".*-", "", modules[[k]])))
  write.table(x = great_peaks, 
              file = paste0("E:/snCUT_RUN/data/peaks/",k,"_marker_peaks_",histone,".bed"),
              col.names = F, 
              row.names = F, 
              sep = "\t", 
              quote = F)
}

#Calculate module score for each single-cell
HN <- AddChromatinModule(object = HN, 
                         features = modules, 
                         genome = BSgenome.Hsapiens.UCSC.hg38, 
                         assay = 'peak', 
                         verbose = TRUE)
module.scores <- HN@meta.data[,21:26]
module.scores[is.na(module.scores)] <- 0

#Extract HN137 cells
HN137_metadata <- HN@meta.data[which(grepl(pattern = "HN137", HN@meta.data$Experiment) == T),]
HN137_metadata2 <- data.frame(Experiment=HN137_metadata$Experiment,
                              HN137PRI=HN137_metadata$HN137PRI,
                              HN137MET=HN137_metadata$HN137MET,
                              HN137PCR=HN137_metadata$HN137PCR,
                              row.names = rownames(HN137_metadata))
HN137_metadata2 <- na.omit(HN137_metadata2)
range(HN137_metadata2$HN137PRI)
min(HN137_metadata2$HN137PRI)
head(HN137_metadata2)

pos_normalised_HN137 <- NULL
for (i in seq_along(HN137_metadata2[,1])) {
  print(paste(i))
  row <- HN137_metadata2[i,-1]
  if (min(row) < 0) {
    row=row+abs(min(row))
  }
  pos_normalised_HN137 <- rbind(pos_normalised_HN137,row)
  
}

HN137_metadata2 <- cbind(HN137_metadata2,pos_normalised_HN137)
colnames(HN137_metadata2) <- c("Experiment", "HN137PRI_prenorm", "HN137MET_prenorm", "HN137PCR_prenorm",
                               "HN137PRI_postnorm", "HN137MET_postnorm", "HN137PCR_postnorm")

HN137_metadata3 <- HN137_metadata2[,c(5:7)]
head(HN137_metadata3)

#Normalise so the scores add to 1
pos_normalised_HN137 <- NULL
for (i in seq_along(HN137_metadata3[,1])) {
  print(paste(i))
  row <- HN137_metadata3[i,]
  rowsum <- sum(row)
  row <- row/rowsum
  pos_normalised_HN137 <- rbind(pos_normalised_HN137,row)
  
}
head(pos_normalised_HN137)
pos_normalised_HN137$Experiment <- HN137_metadata2$Experiment
pos_normalised_HN137$Experiment <- droplevels(pos_normalised_HN137$Experiment)

#Plot ternary plot
library(ggtern)
ggtern(data = pos_normalised_HN137, aes(HN137MET_postnorm,HN137PRI_postnorm, HN137PCR_postnorm, fill=Experiment, shape=Experiment))  +
  geom_mask() +
  geom_point(size = 3, shape = 21, color="black") + 
  Tlab("HN137Pri") +
  Tarrowlab("HN137PRI module") +
  Llab("HN137MET") +
  Larrowlab("HN137MET module") +
  Rlab("HN137PCR") +
  Rarrowlab("HN137PCR module") +
  scale_fill_manual(values = cbPalette[4:6]) +
  theme_custom(
    base_size = 30,
    base_family = "",
    tern.plot.background = "white",
    tern.panel.background = "white",
    col.T = cbPalette[4],
    col.L = cbPalette[5],
    col.R = cbPalette[6],
    col.grid.minor = "grey"
  ) +
  theme_arrowlarge() +
  theme(legend.position = c(0.3,1),
        legend.justification = c(1, 1),
        tern.axis.arrow = element_line(size=2),
        axis.text.x = element_text(vjust = 1))

pos_normalised_HN137$subclust <- 0
pos_normalised_HN137$subclust <- ifelse(test = pos_normalised_HN137$HN137PRI > 0.4 &
                                          pos_normalised_HN137$HN137MET < 0.2 &
                                          pos_normalised_HN137$HN137PCR > 0.4 &
                                          pos_normalised_HN137$Experiment == "HN137PRI", yes = "HN137prePCR", no = pos_normalised_HN137$Experiment)
table(pos_normalised_HN137$Experiment)
table(pos_normalised_HN137$subclust)
HN137prePCR_index <- rownames(pos_normalised_HN137)[which(pos_normalised_HN137$subclust == "HN137prePCR")]

HN@meta.data$subclust <- ifelse(rownames(HN@meta.data) %in% HN137prePCR_index, "HN137prePCR", levels(HN@meta.data$Experiment)[HN@meta.data$Experiment])
table(HN@meta.data$subclust)
HN@meta.data$subclust <- factor(HN@meta.data$subclust, levels = c("HN120PRI", "HN120MET", "HN120PCR",
                                                                  "HN137PRI", "HN137MET", "HN137PCR",
                                                                  "HN137prePCR"))
table(HN@meta.data$subclust)

#Change default identity to subclust identity
Idents(HN) <- factor(HN@meta.data$subclust, levels = c("HN120PRI", "HN120MET", "HN120PCR",
                                                       "HN137PRI", "HN137MET", "HN137PCR",
                                                       "HN137prePCR"))
pos_normalised_HN137[rownames(pos_normalised_HN137) %in% rownames(HN@meta.data)[which(HN@meta.data$subclust == "HN137prePCR")],]

#Dimplot HN137 with HN137prePCR
DimPlot(object = HN, group.by="subclust", label = F, pt.size = 3, 
        cells = rownames(HN@meta.data)[which(HN@meta.data$subclust == "HN137PRI" |
                                               HN@meta.data$subclust == "HN137MET" |
                                               HN@meta.data$subclust == "HN137PCR" |
                                               HN@meta.data$subclust == "HN137prePCR")]) +
  scale_colour_manual(values = cbPalette[c(4:6,8)]) +
  ggtitle(histone) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 24),
        title = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))

#Redo marker analysis with HN137prePCR subpopulation
HN.markers <- FindAllMarkers(object = HN, 
                             assay = "peak", 
                             only.pos = TRUE, 
                             min.pct = 0.2, 
                             test.use = 'LR',
                             latent.vars = 'nCount_peak')
top50 <- HN.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)
closest_genes <- ClosestFeature(HN, regions = top50$gene)
top50$gene2<-closest_genes$gene_name
top50 %>% print(n = 350) 
DoHeatmap(object = HN, 
          assay = "peak", 
          features = top50$gene[151:350], 
          label = F, 
          slot="counts", 
          size = 5, 
          angle = 90,
          cells = rownames(HN@meta.data)[which(HN@meta.data$subclust == "HN137PRI" |
                                                 HN@meta.data$subclust == "HN137MET" |
                                                 HN@meta.data$subclust == "HN137PCR" |
                                                 HN@meta.data$subclust == "HN137prePCR")]) +
  ggplot2::theme(axis.text = element_blank(),
                 legend.position = "none"
  ) 

DoHeatmap(object = HN, 
          assay = "peak", 
          features = top50$gene[c(151:200, 251:300, 301:350)], 
          label = F, 
          slot="counts", 
          size = 5, 
          angle = 90,
          cells = rownames(HN@meta.data)[which(HN@meta.data$subclust == "HN137PRI" |
                                                 HN@meta.data$subclust == "HN137PCR" |
                                                 HN@meta.data$subclust == "HN137prePCR" )]) +
  ggplot2::theme(axis.text = element_blank(),
                 legend.position = "none"
  ) 

top50 %>% print(n = 350)  
HN120PRI_gene<-top50$gene[1:50]  
HN120MET_gene<-top50$gene[51:100]  
HN120PCR_gene<-top50$gene[101:150]
HN137PRI_gene<-top50$gene[151:200]  
HN137MET_gene<-top50$gene[201:250]  
HN137PCR_gene<-top50$gene[251:300]
HN137prePCR_gene<-top50$gene[301:350]
modules <- list(HN120PRI_gene, HN120MET_gene, HN120PCR_gene, 
                HN137PRI_gene, HN137MET_gene, HN137PCR_gene, HN137prePCR_gene)
names(modules) <- c("HN120PRI", "HN120MET", "HN120PCR", 
                    "HN137PRI", "HN137MET", "HN137PCR", "HN137prePCR")

#Re-add chromatin modules
head(HN@meta.data)
HN@meta.data <- HN@meta.data[,-c(21:26)]
HN <- AddChromatinModule(object = HN, modules, BSgenome.Hsapiens.UCSC.hg38, assay = 'peak', verbose = TRUE)

#Plot new module scores: HN137
my_comparisons <- list(c("HN137PRI", "HN137PCR"), c("HN137PCR", "HN137prePCR"), c("HN137PRI", "HN137prePCR"))
VlnPlot(object=HN, features='HN137PCR',
        idents = c("HN137PRI", "HN137MET", "HN137PCR","HN137prePCR")) +
  xlab("")  +
  coord_cartesian(ylim = c(-5, 12), expand = T) +
  ggtitle("HN137PCR Module Score") +
  geom_boxplot(width=0.3, lwd=1, fill = "white") +
  scale_fill_manual(values = cbPalette[c(4:6,8)]) +
  #stat_compare_means(comparisons = my_comparisons, size = 8, label.y = c(11,13,15), label = "p.signif", na.rm = T) +
  theme(title = element_text(size = 24),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_blank(),
        legend.text = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) 


#Check that UMR and FRiP are not confounding factors in HN137prePCR analysis
### HN137prePCR UMR count
my_comparisons_HN137 <- list( c("HN137PRI", "HN137PCR"), c("HN137PRI", "HN137prePCR"), c("HN137PCR", "HN137prePCR") )
ggplot(subset(HN@meta.data, subclust %in% c("HN137PRI", "HN137PCR", "HN137prePCR")), 
       mapping = aes(x = subclust, y = UMRs, fill = subclust))+
  geom_violin(lwd=1) +
  geom_boxplot(width=0.3, lwd=1, fill = "white") +
  scale_fill_manual(values = cbPalette[c(4,6,8)]) +
  theme_bw() +
  coord_cartesian(ylim = (c(0,130000)), expand = T) +
  stat_compare_means(comparisons = my_comparisons_HN137, size = 8, label = "p.signif", vjust = 0.1, label.y = c(100000, 110000, 122000)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=36),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5  ),
        axis.line = element_blank())

### HN137prePCR FRiP count
ggplot(subset(HN@meta.data, subclust %in% c("HN137PRI", "HN137PCR", "HN137prePCR")), 
       mapping = aes(x = subclust, y = Pct_reads_in_peaks, fill = subclust))+
  geom_violin(lwd=1) +
  geom_boxplot(width=0.3, lwd=1, fill = "white") +
  scale_fill_manual(values = cbPalette[c(4,6,8)]) +
  theme_bw() +
  ylim(c(0,90)) +
  ylab("Percentage of reads in peaks\n") +
  stat_compare_means(comparisons = my_comparisons_HN137, size = 8, label = "p.signif") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5),
        axis.line = element_blank())

#Linear regression to show that UMRs are not correlated with module scores
ggplot(subset(HN@meta.data, subclust %in% c("HN137prePCR")), 
       mapping = aes(x = UMRs, y = HN137PCR))+
  geom_point() +
  geom_smooth(method = "lm") +
  stat_regline_equation(label.y = 9, aes(label = ..eq.label..), size = 7) +
  stat_regline_equation(label.y = 8.5, aes(label = ..rr.label..), size = 7) +
  theme_bw() +
  ylab("HN137PCR Module Score") +
  xlab("HN137prePCR - UMRs") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size=28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5),
        axis.line = element_blank()) 

ggplot(subset(HN@meta.data, subclust %in% c("HN137prePCR")), 
       mapping = aes(x = Pct_reads_in_peaks, y = HN137PCR))+
  geom_point() +
  geom_smooth(method = "lm") +
  stat_regline_equation(label.y = 9, aes(label = ..eq.label..), size = 7) +
  stat_regline_equation(label.y = 8.5, aes(label = ..rr.label..), size = 7) +
  theme_bw() +
  ylab("HN137PCR Module Score") +
  xlab("HN137prePCR - Fraction Reads in Peaks") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size=28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5),
        axis.line = element_blank()) 
