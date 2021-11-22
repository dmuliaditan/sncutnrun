
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


################################################
############ snCUT&RUN data analysis ###########
################################################
############# v. 11 October 2021 ###############

#For Signac analysis, two files are needed: 
#1) A fragment file, which is basically a bed file with the following columns:
#   chr, start, end, barcode, number of PCR duplicates of that fragment
#2) Barcode file with single-cell metadata


#Fetch the fragment file and create the bin matrix
#Define the path of the fragment file
setwd("/mnt/raid5/cutnrun/userFiles/daniel/DM112021/K27ac")

#For H3K27ac
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
                                                                      "HN120PRI", "HN120MET", "HN137PCR"))


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
HN <- readRDS(paste0(histone,"_orig.rds")) #In all subsequent runs we can just load it from saved data

#Filter cells:
HN <- subset(
  x = HN,
  subset = 
    UMRs > 1000 &
    UMRs < 100000 &
    Pct_reads_in_peaks > 25 &
    Pct_reads_in_blacklist < 0.5 &
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
                 dims = 2:30)
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
DimPlot(object = HN, 
        group.by="Experiment", 
        #cols = "Set3", 
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
  theme(axis.title.x = element_text(size = 24),
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
  max.cutoff = 'q95'
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


#Trajectory analysis with Monocle3, HN120 and HN137 separately
#Here HN120 cell lines were used as en example
HN120 <- HN[, HN$Experiment %in% c("HN120PRI", "HN120MET", "HN120PCR")]
HN.cds <- as.cell_data_set(x = HN120)
HN.cds <- cluster_cells(cds = HN.cds, 
                        reduction_method = "UMAP",
                        cluster_method = 'leiden', 
                        num_iter = 100, 
                        partition_qval = 0.01, 
                        verbose = T)
HN.cds <- learn_graph(HN.cds, 
                      use_partition = TRUE, 
                      verbose = T)

#Order cells, use primary cells as root
HN.cds <- order_cells(HN.cds, 
                      reduction_method = "UMAP",
                      root_cells = rownames(HN120@meta.data)[HN120@meta.data$Experiment == "HN120PRI"])

#Plot pseudotime
plot_cells(
  cds = HN.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = F,
  cell_size = 2
) +
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 24),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
remove(HN120)# Make space

#Call peaks with MACS2
peaks <- CallPeaks(
  object = HN,
  assay = "bin",
  group.by = "Experiment",
  broad = T,
  macs2.path = "/home/cutnrun/miniconda3/bin/macs2",
  outdir = getwd(),
  combine.peaks = T
)

save.image(file = paste0(histone,".RData")) #Checkpoint, optional

#Plot coverage gene(s) of interest
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


############################ Motif analysis ######################################
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

#Create peak assay
#Create peak feature matrix
peak_matrix <- FeatureMatrix(
  fragments = Fragments(HN),
  features = peaks
)
HN[["peak"]] <- CreateChromatinAssay(
  counts = peak_matrix[,colnames(HN)],
  sep = c(":", "-"),
  fragments = fragment.path,
)

# Add motif information, this step will take a while!
HN <- AddMotifs(
  object = HN,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm,
  assay = 'peak'
  
)

#Get meta features
meta.feature <- GetAssayData(HN, 
                             assay="peak", 
                             slot="meta.features")

#Run ChromVAR
HN <- RunChromVAR(
  object = HN,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = 'peak'
  
)

#Save ChromVAR results
write.table (HN$chromvar@data, paste0(histone, "_deviation.txt"), 
             row.names=TRUE, 
             col.names=TRUE)
#It's good to save image here

#Motif analysis: everything at the same time
comparisons <- list()
comparisons[[1]] <- c("HN120PRI", "HN120MET")
comparisons[[2]] <- c("HN120PRI", "HN120PCR")
comparisons[[3]] <- c("HN120MET", "HN120PRI")
comparisons[[4]] <- c("HN120PCR", "HN120PRI")
comparisons[[5]] <- c("HN120PRI", "HN120MET")
comparisons[[6]] <- c("HN120PRI", "HN137PCR")
comparisons[[7]] <- c("HN120MET", "HN120PRI")
comparisons[[8]] <- c("HN137PCR", "HN120PRI")

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
  peaks.matched <- MatchRegionStats(
    meta.feature=meta.feature,
    query.feature=meta.feature[top.da.peak,],
    n=nrow(meta.feature),
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
              file = paste0(histone,"_enriched_motifs",comparisons[[o]][[1]], "_vs_",comparisons[[o]][[2]],".txt"),
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
              file = paste0(histone,"_peaks",comparisons[[o]][[1]], "_vs_",comparisons[[o]][[2]],".bed"),
              col.names = F, 
              row.names = F, 
              sep = "\t", 
              quote = F)
  
}

#Specific comparisons
da_peaks <- FindMarkers(
  object = HN,
  ident.1 = 'HN120PCR',
  ident.2 = 'HN120PRI',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peak',
  min.pct = 0.05,
  assay = 'peak'
)

#Get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
peaks.matched <- MatchRegionStats(
  meta.feature=meta.feature,
  query.feature=meta.feature[top.da.peak,],
  n=nrow(meta.feature),
  verbose = T
)

#Test enrichment
enriched.motifs <- FindMotifs(
  object = HN,
  features = top.da.peak,
  assay = 'peak',
  background=peaks.matched
)
enriched.motifs[1:20,]

#Save motif images and save ChromVAR images of motif score per cell
for (k in c(1:10)){
  print(paste(k))
  png(filename = paste0(enriched.motifs$motif.name[k],"_motif.png"), width = 1000, height = 150, units = "px")
  l = MotifPlot(
    object = HN,
    motifs = rownames(enriched.motifs)[k],
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
motif_ind <- which(enriched.motifs$motif.name == "SP1")
FeaturePlot(
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


################### HN120PRIMet and HN120PRIPCR analysis ###########################
Idents(HN) <- factor(Idents(HN), 
                     levels = c("HN120PRI", "HN120MET", "HN120PCR", "HN120PRI", "HN120MET", "HN137PCR"))

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
HN120PRI_gene<-top50$gene[151:200]  
HN120MET_gene<-top50$gene[201:250]  
HN137PCR_gene<-top50$gene[251:300]  
modules <- list(HN120PRI_gene, HN120MET_gene, HN120PCR_gene,
                HN120PRI_gene, HN120MET_gene, HN137PCR_gene)
names(modules) <- c("HN120PRI", "HN120MET", "HN120PCR",
                    "HN120PRI", "HN120MET", "HN137PCR")

#Save module peaks as bed files
#for (k in names(modules)) {
# print(paste(k))
# startend_peaks <- sub("^[^-]*-", "", modules[[k]])
# great_peaks <- data.frame(chr=sub("-.*", "", modules[[k]]),
#                           start=as.numeric(sub("-.*", "", startend_peaks)),
#                           end= as.numeric(sub(".*-", "", modules[[k]])))
# write.table(x = great_peaks, 
#             file = paste0("E:/snCUT_RUN/data/peaks/",k,"_marker_peaks_",histone,".bed"),
#             col.names = F, 
#             row.names = F, 
#             sep = "\t", 
#             quote = F)
#}

#Calculate module score for each single-cell
HN <- AddChromatinModule(object = HN, 
                          features = modules, 
                          genome = BSgenome.Hsapiens.UCSC.hg38, 
                          assay = 'peak', 
                          verbose = TRUE)


my_comparisons <- list( c("HN120PRI", "HN120MET"), c("HN120MET", "HN120PCR"), c("HN120PRI", "HN120PCR"))
VlnPlot(object=HN, features='HN120PCR', 
        idents = c("HN120PRI", "HN120MET", "HN120PCR"), cols = c("#F8766D", "#D39200", "#00BA38")) +
  geom_boxplot(width=0.3, lwd=1, fill = "white") +
  xlab("")  +
  ylim(c(-5, 15)) +
  theme(title = element_text(size = 28),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 24),
        legend.text = element_text(size = 24),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  stat_compare_means(comparisons = my_comparisons, size = 8, label.y = c(9,11,13), label = "p.signif")

#Put module score as a data frame
module.scores <- HN@meta.data[,21:26]
module.scores[is.na(module.scores)] <- 0

#HN120PRIMet analysis
HN120_metadata <- HN@meta.data[which(grepl(pattern = "HN120", HN@meta.data$Experiment) == T),]
HN120_index <- which(rownames(module.scores) %in% row.names(HN120_metadata))
module.scores.HN120 <- module.scores[HN120_index,]
module.scores.HN120$orig_ident <- HN120_metadata$Experiment

module.scores.HN120$Ident <- 0
module.scores.HN120$Ident <- ifelse(test = module.scores.HN120$HN120PRI > quantile(module.scores.HN120$HN120PRI)[4]
                                    & module.scores.HN120$HN120MET > quantile(module.scores.HN120$HN120MET)[4]
                                    & module.scores.HN120$orig_ident == "HN120PRI", 
                                    yes = "HN120PRIMet", module.scores.HN120$orig_ident)

HN120PRIMet_index <- rownames(module.scores.HN120)[which(module.scores.HN120$Ident == "HN120PRIMet")]
table(module.scores.HN120$orig_ident) #Old ident
table(module.scores.HN120$Ident) #New ident

HN@meta.data$subclust <- ifelse(rownames(HN@meta.data) %in% HN120PRIMet_index, "HN120PRIMet", HN@meta.data$Experiment)
HN@meta.data$subclust <- factor(HN@meta.data$subclust, levels = c("HN120PRI", "HN120MET", "HN120PCR",
                                                                  "HN120PRIMet", "HN120PRI", 
                                                                  "HN120MET", "HN137PCR"))

#Change default identity to subclust identity
Idents(HN) <- factor(HN@meta.data$subclust, levels = c("HN120PRI", "HN120MET", "HN120PCR",
                                                       "HN120PRIMet", "HN120PRI", 
                                                       "HN120MET", "HN137PCR"))

#Dimplot HN120 with HN120PRIMet
DimPlot(object = HN, group.by="subclust", label = F, pt.size = 3, 
        cells = rownames(HN@meta.data)[which(HN@meta.data$subclust == "HN120PRI" |
                                               HN@meta.data$subclust == "HN120MET" |
                                               HN@meta.data$subclust == "HN120PCR" |
                                               HN@meta.data$subclust == "HN120PRIMet")]) +
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

#Redo marker analysis with HN120PRIMet subpopulation
HN.markers <- FindAllMarkers(object = HN, 
                             assay = "peak", 
                             only.pos = TRUE, 
                             min.pct = 0.2, 
                             test.use = 'LR',
                             latent.vars = 'nCount_peak')
top50 <- HN.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)
top50 <- top50[which(grepl(pattern = "HN120PCR", x = top50$cluster) == T),] 
closest_genes <- ClosestFeature(HN, regions = top50$gene)
top50$gene2<-closest_genes$gene_name
DoHeatmap(object = HN, 
          assay = "peak", 
          features = top50$gene, 
          label = F, 
          slot="counts", 
          size = 5, 
          angle = 90,
          cells = rownames(HN@meta.data)[which(HN@meta.data$subclust == "HN120PRI" |
                                                 HN@meta.data$subclust == "HN120MET" |
                                                 HN@meta.data$subclust == "HN120PCR" |
                                                 HN@meta.data$subclust == "HN120PRIMet")]) +
  ggplot2::theme(axis.text = element_blank(),
                 legend.position = "none"
                 ) 

top50 %>% print(n = 350)  
HN120PRI_gene<-top50$gene[1:15]  
HN120MET_gene<-top50$gene[16:65]  
HN120PCR_gene<-top50$gene[66:115]
HN120PRIMet_gene<-top50$gene[116:165]
HN120PRI_gene<-top50$gene[166:215]  
HN120MET_gene<-top50$gene[216:265]  
HN137PCR_gene<-top50$gene[266:315]  
modules <- list(HN120PRI_gene, HN120MET_gene, HN120PCR_gene, HN120PRIMet_gene,
                HN120PRI_gene, HN120MET_gene, HN137PCR_gene)
names(modules) <- c("HN120PRI", "HN120MET", "HN120PCR", "HN120PRIMet",
                    "HN120PRI", "HN120MET", "HN137PCR")

#Re-add chromatin modules
HN@meta.data <- HN@meta.data[,-c(21:26)]
HN <- AddChromatinModule(object = HN, modules, BSgenome.Hsapiens.UCSC.hg38, assay = 'peak', verbose = TRUE)

#Plot new module scores: HN120
my_comparisons <- list( c("HN120PRI", "HN120MET"), c("HN120MET", "HN120PRIMet"), c("HN120PRI", "HN120PRIMet"))
VlnPlot(object=HN, features='HN120MET',
        idents = c("HN120PRI", "HN120MET", "HN120PRIMet","HN120PCR")) +
  xlab("")  +
  ylim(c(-5, 15)) +
  ggtitle("HN120MET Module Score") +
  geom_boxplot(width=0.3, lwd=1, fill = "white") +
  scale_fill_manual(values = c("#F8766D", "#D39200", "#00BA38", "#DB72FB")) +
  theme(title = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  stat_compare_means(comparisons = my_comparisons, size = 8, label.y = c(8,10,12), label = "p.signif")

#Check that UMR and FRiP are not confounding factors in PriMet analysis
### HN120PRIMet UMR count
my_comparisons_HN120 <- list( c("HN120PRI", "HN120MET"), c("HN120PRI", "HN120PRIMet"), c("HN120MET", "HN120PRIMet") )
ggplot(subset(HN@meta.data, subclust %in% c("HN120PRI", "HN120MET", "HN120PRIMet")), 
       mapping = aes(x = subclust, y = UMRs, fill = subclust))+
  geom_violin(lwd=1) +
  geom_boxplot(width=0.3, lwd=1) +
  scale_fill_brewer(palette = "Accent") +
  theme_bw() +
  ylim(c(0,130000)) +
  stat_compare_means(comparisons = my_comparisons_HN120, size = 8, label = "p.signif", vjust = 0.1, label.y = c(100000, 110000, 122000)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=36),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5  ),
        axis.line = element_blank())

### HN120PRIMet FRiP count
ggplot(subset(HN@meta.data, subclust %in% c("HN120PRI", "HN120MET", "HN120PRIMet")), 
       mapping = aes(x = subclust, y = Pct_reads_in_peaks, fill = subclust))+
  geom_violin(lwd=1) +
  geom_boxplot(width=0.3, lwd=1) +
  scale_fill_brewer(palette = "Accent") +
  theme_bw() +
  ylim(c(0,70)) +
  ylab("Percentage of reads in peaks\n") +
  stat_compare_means(comparisons = my_comparisons_HN120, size = 8, label = "p.signif") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5),
        axis.line = element_blank())

#Linear regression to show that UMRs are not correlated with module scores
ggplot(subset(HN@meta.data, subclust %in% c("HN120PRIMet")), 
       mapping = aes(x = Pct_reads_in_peaks, y = HN120MET))+
  geom_point() +
  geom_smooth(method = "lm") +
  stat_regline_equation(label.y = 9, aes(label = ..eq.label..), size = 7) +
  stat_regline_equation(label.y = 8.5, aes(label = ..rr.label..), size = 7) +
  theme_bw() +
  ylab("HN120MET Module Score") +
  xlab("HN120PRIMet - Pct reads in Peaks") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size=28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5),
        axis.line = element_blank()) 


#HN120PRIPCR analysis
HN@meta.data <- HN@meta.data[,-c(21:27)]
HN137_metadata <- HN@meta.data[which(grepl(pattern = "HN137", HN@meta.data$Experiment) == T),]
HN137_index <- which(rownames(module.scores) %in% row.names(HN137_metadata))
module.scores.HN137 <- module.scores[HN137_index,]
module.scores.HN137$orig_ident <- HN137_metadata$Experiment

module.scores.HN137$Ident <- 0
module.scores.HN137$Ident <- ifelse(test = module.scores.HN137$HN120PRI > quantile(module.scores.HN137$HN120PRI)[4]
                                    & module.scores.HN137$HN137PCR > quantile(module.scores.HN137$HN137PCR)[4]
                                    & module.scores.HN137$orig_ident == "HN120PRI", 
                                    yes = "HN120PRIPCR", module.scores.HN137$orig_ident)

HN120PRIPCR_index <- rownames(module.scores.HN137)[which(module.scores.HN137$Ident == "HN120PRIPCR")]
table(module.scores.HN137$orig_ident) #Old ident
table(module.scores.HN137$Ident) #New ident

HN@meta.data$subclust <- ifelse(rownames(HN@meta.data) %in% HN120PRIPCR_index, "HN120PRIPCR", HN@meta.data$Experiment)
HN@meta.data$subclust <- factor(HN@meta.data$subclust, levels = c("HN120PRI", "HN120MET", "HN120PCR",
                                                                 "HN120PRI", "HN120MET", 
                                                                 "HN137PCR", "HN120PRIPCR"))

#Change default identity to subclust identity
Idents(HN) <- factor(HN@meta.data$subclust, levels = c("HN120PRI", "HN120MET", "HN120PCR",
                                                       "HN120PRI", "HN120MET", "HN137PCR", "HN120PRIPCR"))

#Dimplot HN137 with HN120PRIPCR
DimPlot(object = HN, group.by="subclust", label = F, pt.size = 3, 
        cells = rownames(HN@meta.data)[which(HN@meta.data$subclust == "HN120PRI" |
                                               HN@meta.data$subclust == "HN120MET" |
                                               HN@meta.data$subclust == "HN137PCR" |
                                               HN@meta.data$subclust == "HN120PRIPCR")]) +
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

#Redo marker analysis with HN120PRIPCR subpopulation
HN.markers <- FindAllMarkers(object = HN, 
                             assay = "peak", 
                             only.pos = TRUE, 
                             min.pct = 0.2, 
                             test.use = 'LR',
                             latent.vars = 'nCount_peak')
top50 <- HN.markers %>% group_by(cluster) %>% top_n(50, avg_log2FC)
DoHeatmap(object = HN, assay = "peak", features = top50$gene, label = TRUE, slot="data", size = 5, angle = 90) +
  ggplot2::theme(axis.text.y = element_blank(),
                 legend.text = element_text(size = 16),
                 legend.title = element_blank()) 

top50 %>% print(n = 350)  
HN120PRI_gene<-top50$gene[1:50]  
HN120MET_gene<-top50$gene[51:100]  
HN120PCR_gene<-top50$gene[101:150]
HN120PRI_gene<-top50$gene[151:183]  
HN120MET_gene<-top50$gene[184:233]  
HN137PCR_gene<-top50$gene[234:283]  
HN120PRIPCR_gene<-top50$gene[284:333]
modules <- list(HN120PRI_gene, HN120MET_gene, HN137PCR_gene,
                HN120PRI_gene, HN120MET_gene, HN137PCR_gene, HN120PRIPCR_gene)
names(modules) <- c("HN120PRI", "HN120MET", "HN137PCR",
                    "HN120PRI", "HN120MET", "HN137PCR", "HN120PRIPCR")

#Re-add chromatin modules
HN@meta.data <- HN@meta.data[,-c(21:26)]
HN <- AddChromatinModule(object = HN, modules, BSgenome.Hsapiens.UCSC.hg38, assay = 'peak', verbose = TRUE)

#Plot new module scores: HN137
my_comparisons <- list( c("HN120PRI", "HN137PCR"), c("HN137PCR", "HN120PRIPCR"), c("HN120PRI", "HN120PRIPCR"))
VlnPlot(object=HN, features='HN137PCR') +
  xlab("")  +
  ylim(c(-6, 25)) +
  ggtitle("HN137PCR Module Score") +
  scale_fill_manual(values = c("#F8766D", "#D39200", "#00BA38", "#DB72FB", "#00C19F", "#619CFF", "#FF61C3")) +
  theme(title = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  stat_compare_means(comparisons = my_comparisons, size = 8, label.y = c(14,18,22), label = "p.signif")

#Check that UMR and FRiP are not confounding factors in PriPCR analysis
### HN120PRIPCR UMR count
my_comparisons_HN137 <- list( c("HN120PRI", "HN137PCR"), c("HN120PRI", "HN120PRIPCR"), c("HN137PCR", "HN120PRIPCR") )
ggplot(subset(HN@meta.data, subclust %in% c("HN120PRI", "HN137PCR", "HN120PRIPCR")), 
       mapping = aes(x = subclust, y = UMRs, fill = subclust))+
  geom_violin(lwd=1) +
  geom_boxplot(width=0.3, lwd=1) +
  scale_fill_brewer(palette = "Accent") +
  theme_bw() +
  ylim(c(0,130000)) +
  stat_compare_means(comparisons = my_comparisons_HN137, size = 8, label = "p.signif", vjust = 0.1, label.y = c(100000, 110000, 122000)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=36),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5  ),
        axis.line = element_blank())

### HN120PRIPCR FRiP count
ggplot(subset(HN@meta.data, subclust %in% c("HN120PRI", "HN137PCR", "HN120PRIPCR")), 
       mapping = aes(x = subclust, y = Pct_reads_in_peaks, fill = subclust))+
  geom_violin(lwd=1) +
  geom_boxplot(width=0.3, lwd=1) +
  scale_fill_brewer(palette = "Accent") +
  theme_bw() +
  ylim(c(0,70)) +
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
ggplot(subset(HN@meta.data, subclust %in% c("HN120PRIPCR")), 
       mapping = aes(x = UMRs, y = HN137PCR))+
  geom_point() +
  geom_smooth(method = "lm") +
  stat_regline_equation(label.y = 9, aes(label = ..eq.label..), size = 7) +
  stat_regline_equation(label.y = 8.5, aes(label = ..rr.label..), size = 7) +
  theme_bw() +
  ylab("HN137PCR Module Score") +
  xlab("HN120PRIPCR - UMRs") +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size=28),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(colour = "black", fill=NA,size=0.5),
        axis.line = element_blank()) 

#HN120PRIMet motif analysis
da_peaks_primet <- FindMarkers(
  object = HN,
  ident.1 = 'HN120PRIPCR',
  ident.2 = 'HN120PRI',
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peak',
  min.pct = 0.05,
  assay = 'peak'
)

top.da.peak_primet <- rownames(da_peaks_primet[da_peaks_primet$p_val < 0.005, ])
peaks.matched <- MatchRegionStats(
  meta.feature=meta.feature,
  query.feature=meta.feature[top.da.peak_primet,],
  n=nrow(meta.feature),
  verbose = T
)


#Test enrichment
enriched.motifs <- FindMotifs(
  object = HN,
  features = top.da.peak_primet,
  assay = 'peak',
  background=peaks.matched
)
enriched.motifs[1:50,]
sig_ind <- which(enriched.motifs$pvalue < 0.005)
motifs <- data.frame(motifs=enriched.motifs$motif.name[sig_ind])
write.table(x = motifs,
            file = "HN120PCR_HN120PRI_TF_enriched.bed",
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

MotifPlot(
  object = HN,
  motifs = rownames(enriched.motifs)[1:4],
  assay = "peak"
) +
  theme(axis.text = element_text(size = 18),
        title = element_text(size=22),
        strip.text = element_text(size=28)
  )

enriched.motifs[1:50,]

enriched.motifs[order(enriched.motifs$fold.enrichment, decreasing = T),]

motif_ind <- which(enriched.motifs$motif.name == "HSF2")
FeaturePlot(
  object = HN,
  features = enriched.motifs$motif[motif_ind],
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 2
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

module.scores[hn137_pripcr.index,]

#rGREAT analysis
da_peaks <- FindMarkers(
  object = HN,
  ident.1 = c('HN120PCR'),
  ident.2 = c('HN120PRI'),
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peak',
  min.pct = 0.05,
  assay = 'peak'
)

#Get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
startend_peaks <- sub("^[^-]*-", "", top.da.peak)
great_peaks <- data.frame(chr=sub("-.*", "", top.da.peak),
                          start=as.numeric(sub("-.*", "", startend_peaks)),
                          end= as.numeric(sub(".*-", "", top.da.peak)))

job <- submitGreatJob(gr = great_peaks, species = 'hg38', rule = 'twoClosest',
                      adv_twoDistance = 2000)
tb <- getEnrichmentTables(job, category = "GO")

names(tb)
head(tb$`GO Molecular Function`, n=30)
availableOntologies(job)
res <- plotRegionGeneAssociationGraphs(job)
plotRegionGeneAssociationGraphs(job, type = 2)
availableCategories(job)
availableOntologies(job, category = "Phenotype")

GO_table <- tb$`GO Biological Process`[1:10,c(2,9)]
write.table(x = GO_table, file = "GO_K4me3_HN120PRI_enriched_vs_HN137PCR.txt", quote = F, sep = "\t", row.names = F, col.names = F)


#Only use in windows
DefaultAssay(HN) <- 'bin'
fragment<-Fragments(HN)
fragment<-UpdatePath(fragment[[1]],"E:/snCUT_RUN/data/signac_fragments/K4me3_combined_sorted.bed.gz")
Fragments(HN)<-NULL
HN <- SetAssayData(HN, slot = "fragments", new.data = fragment)

#Coverageplot
Idents(HN) <- factor(HN@meta.data$Experiment, levels = c("HN120PRI", "HN120MET", "HN120PCR", "HN120PRI", "HN120MET", "HN137PCR"))

cov_plot <- CoveragePlot(
  object = HN,
  region = "TP63",
  annotation = T,
  peaks = F,
  idents = c("HN120PRI", "HN120MET", "HN137PCR"),
  extend.upstream = 10000,
  extend.downstream = 10000
)
cov_plot +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))
