#RNA-seq analysis with R after salmon quantification
#Tutorial RNA-seq workflow: gene-level exploratory analysis and differential expression, Michael I. Love, Simon Anders, Vladislav Kim and Wolfgang Huber, 16 October, 2019
#https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

#Set working directory
setwd("/mnt/raid5/cutnrun/rnaseq")

#Load required libraries
library(tximeta)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)
library(apeglm)
library(genefilter)
library(mygene)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(fgsea)
library(tidyverse)

#Locate salmon output (.quant files)
dir <- "/mnt/raid5/cutnrun/rnaseq/salmon_output"

#Import salmon output to R with tximeta
list.files(dir) #Check if the data are complete
list.files(file.path(dir, "quants"))

#Create an appropriate metadata file detailing library and sample names, as well as necessary annotations
csvfile <- file.path(dir, "sample_meta.csv")
coldata <- read.csv(csvfile, stringsAsFactors=FALSE, header=T) #Import into R
coldata

#Refactorize data so that the data are of correct type and order
coldata$Sample <- factor(coldata$Sample, levels = c("HN120Pri", "HN120Met", "HN120PCR", "HN137Pri", "HN137Met", "HN137PCR"))
coldata$Type <-  factor(coldata$Type, levels = c("Primary", "Metastatic", "CisplatinResistant"))
coldata$Patient <- factor(coldata$Patient)
coldata$names <- coldata$Library
coldata$files <- file.path(dir, coldata$names, "quant.sf")
file.exists(coldata$files) #Check if data are in the correct directory

#Use tximeta to import all the RNAseq .quant files
se <- tximeta(coldata)
dim(se)
head(rownames(se))

#Process to gene level
gse <- summarizeToGene(se) 
dim(gse)
head(rownames(gse))
data(gse)
gse

#Use the SummarizedExperiment package to explore the data
assayNames(gse)
head(assay(gse), 3)
colSums(assay(gse))
rowRanges(gse)
seqinfo(rowRanges(gse))
colData(gse)

#Final assay is a SummarizedExperiment object with the following components
#coldata: the samples metadata
#ranges: ranges of the features (transcripts)
#count: the actual data itself

#Now we are ready for the actual DESeq2 analysis
#Analysis with DESeq2
#Check how many mapped fragments are per sample
round(colSums(assay(gse)) / 1e6, 1 )

#Construct DESeq2 object
#Make sure to have the correct design!
dds <- DESeqDataSet(gse, design = ~ Sample)

#Exploratory analysis and visualization
nrow(dds)

#Filtering: keep only features with at least 10 counts or more in 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
nrow(dds)

#Exploration with transformations to mitigate variance differences due to low or very high counts
#With VST:
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

#With rlog:
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

#Check transformation effects
dds <- estimateSizeFactors(dds)
df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
         mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
colnames(df)[1:2] <- c("x", "y")  
lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  
  

####### Z-score calculation ######
#Use the vst transformed data and the scale function to create the raw Z-score object
vsd <- assay(vst(dds))
Z <- t(scale(t(vsd)))
#Remove rows with NaN
Z <- Z[complete.cases(Z),] #Final Z-score matrix
Z[1:5,1:5]

#Change ENSG gene annotation to more readable gene SYMBOL
gene.list <- rownames(Z)
gene.list <- substr(x=gene.list,start=1,stop = 15)
symbol.list <- getGenes(geneids=gene.list, fields='symbol')
head(symbol.list)
rownames(Z) <- symbol.list$symbol #Replace with the new geneIDs
Z[1:5,1:5]

#Take mean Z-score for the triplicates as representative Z-score of that sample
Z <- as.data.frame(Z)
Z$HN120PRI <- rowMeans(Z[,1:3])
Z$HN120MET <- rowMeans(Z[,4:6])
Z$HN120PCR <- rowMeans(Z[,7:9])
Z$HN137PRI <- rowMeans(Z[,10:12])
Z$HN137MET <- rowMeans(Z[,13:15])
Z$HN137PCR <- rowMeans(Z[,16:18])
write.table(x = Z[,19:24], file = "/home/daniel/daniel_new/DM_SNCUTRUN_RNAseq_HN120_HN137_Zscore_allgenes.txt", quote = F, sep = "\t",
            row.names = T, col.names = T)
            
####### Clustering and similarity analysis ########
#There are multiple ways to check which samples are similar, for example:
#1. Through calculating the euclidean distances between the samples
#2. Through calculating the Poisson distances between the samples
#3. Through PCA analysis

#Calculate distances between samples
#1. Calculate Euclidean distances between samples:
sampleDists <- dist(t(assay(rld)))
sampleDists
#Visualize
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Type, rld$Patient, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#2. Calculate Poisson distances
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$Type, dds$Patient, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

#3. Visualize with PCA
#From DESeq2 package
plotPCA(rld, intgroup = c("Type", "Patient"))

#From scratch
pcaData <- plotPCA(rld, intgroup = c( "Type", "Patient"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = Type, shape = Patient)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with rlog data")

#With glmpca
#…we propose the use of GLM-PCA, a generalization of PCA to exponential family likelihoods. 
# GLM-PCA operates on raw counts, avoiding the pitfalls of normalization. 
# We also demonstrate that applying PCA to deviance or Pearson residuals provides a useful and fast approximation to GLM-PCA.
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$Type <- dds$Type
gpca.dat$Patient <- dds$Patient
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = Type, shape = Patient)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA") +
  xlab("Dim1") +
  ylab("Dim2") +
  theme_bw() +
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 28),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        title = element_text(size = 24))

###### Differential analysis ######
#Running the differential expression pipeline
dds <- DESeq(dds)

#Building the results table
res <- results(dds)
res
resultsNames(dds)
res <- results(dds, contrast=c("Sample","HN137Met", "HN137Pri"))
mcols(res, use.names = TRUE)
summary(res)

#Lower the threshold
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)

#Count genes with >=2 fold change
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

#Multiple testing
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
sum(res$padj < 0.1, na.rm=TRUE)
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

#Bland-Altman plot
resultsNames(dds)
res <- lfcShrink(dds, coef="Sample_HN137Met_vs_HN137Pri", type="apeglm")
plotMA(res, ylim = c(-10, 10))

#Gene clustering
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)
mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("Sample","Type")])
pheatmap(mat, annotation_col = anno)

#Annotating results
#Change ENSG geneIDs to more readable gene symbols
columns(org.Hs.eg.db)
ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
res$enstrans <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENSEMBL",
                     keytype="ENSEMBL",
                     multiVals="first")
res

#Check top differential genes
resOrdered <- res[order(res$pvalue),]
head(resOrdered)
write.table(x = resOrdered[,c(7,6,2)],
            file = "HN137Met_vs_HN137Pri_top_differentially regulated gene.txt", 
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T
)

#Check top upregulated genes
resOrdered <- res[order(res$log2FoldChange, decreasing = TRUE ),]
resOrdered
write.table(x = resOrdered[,c(7,6,2)],
            file = "HN137Met_vs_HN137Pri_top_upregulated gene.txt", 
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T
)

plotCounts(dds, gene = "ENSG00000137693.14" , intgroup=c("Sample"))

#Gene set enrichment analysis using fgsea package
res$row <- rownames(res)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
ens2symbol <- as_tibble(data.frame(ENSEMBL=rownames(res), SYMBOL=res$symbol))
ens2symbol

res2 <- inner_join(as_tibble(res), ens2symbol, by=c("row"="ENSEMBL"))
res2

res3 <- res2 %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))
res3

ranks <- deframe(res3)
head(ranks, 20)

# Load the pathways into a named list
pathways.hallmark <- gmtPathways("h.all.v7.4.symbols.gmt")
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
#Tidy the results
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

save.image("RNAseq_DESeq_analysis_02122021.RData")
