#RNA-seq analysis with R
#Tutorial RNA-seq workflow: gene-level exploratory analysis and differential expression, Michael I. Love, Simon Anders, Vladislav Kim and Wolfgang Huber, 16 October, 2019
#https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

#Import salmon output to R with tximeta
setwd("/mnt/raid5/cutnrun/rnaseq")
dir <- "/mnt/raid5/cutnrun/rnaseq/salmon_output"

list.files(dir)

list.files(file.path(dir, "quants"))

csvfile <- file.path(dir, "sample_meta.csv")
coldata <- read.csv(csvfile, stringsAsFactors=FALSE, header=T)
coldata

coldata$Sample <- factor(coldata$Sample, levels = c("HN120Pri", "HN120Met", "HN120PCR", "HN137Pri", "HN137Met", "HN137PCR"))
coldata$Type <-  factor(coldata$Type, levels = c("Primary", "Metastatic", "CisplatinResistant"))
coldata$Patient <- factor(coldata$Patient)
coldata$names <- coldata$Library
coldata$files <- file.path(dir, coldata$names, "quant.sf")

coldata <- coldata[coldata$Patient == "HN137",]

file.exists(coldata$files)

library("tximeta")
se <- tximeta(coldata)
dim(se)
head(rownames(se))

gse <- summarizeToGene(se) 

dim(gse)
head(rownames(gse))

data(gse)
gse

library(SummarizedExperiment)
assayNames(gse)
head(assay(gse), 3)
colSums(assay(gse))
rowRanges(gse)
seqinfo(rowRanges(gse))
colData(gse)

#Final assay is a SummarizedExperiment with the following components
#coldata: the samples metadata
#ranges: ranges of the features (transcripts)
#count: the actual data itself

#Analysis with DESeq2
#Check how many mapped fragments
round( colSums(assay(gse)) / 1e6, 1 )

#Construct DESeq2 object
library("DESeq2")
dds <- DESeqDataSet(gse, design = ~ Sample)

#Exploratory analysis and visualization
nrow(dds)

#Keep only features with at least 10 counts or more in 3 samples
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
library("dplyr")
library("ggplot2")

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

#Calculate distances between samples
#Calculate Euclidean distances between samples:
sampleDists <- dist(t(assay(rld)))
sampleDists
#Visualize
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Type, rld$Patient, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#Calculate Poisson distances
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$Type, dds$Patient, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

#Visualize with PCA
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
library("glmpca")
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

#Differential analysis
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

#Other comparisons
#Check between HN137Pri and HN137Met
#results(dds, contrast = c("Sample", "Metastatic", "Primary"))

#Multiple testing
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
sum(res$padj < 0.1, na.rm=TRUE)
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

#Plotting results
#Count Plot
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("Sample"))

#Beeswarmplot
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("Type","Sample"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Type, y = count, color = Sample)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

#Bland-Altman plot
library("apeglm")
resultsNames(dds)
res <- lfcShrink(dds, coef="Type_Metastatic_vs_Primary", type="apeglm")
plotMA(res, ylim = c(-5, 5))

#Gene clustering
library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)
mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("Sample","Type")])
pheatmap(mat, annotation_col = anno)

#Annotating results
library("AnnotationDbi")
library("org.Hs.eg.db")
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
#Check top differential genes
resOrdered <- res[order(res$pvalue),]
head(resOrdered)

#Check top upregulated genes
resOrdered <- res[order(res$log2FoldChange, decreasing = TRUE ),]
resOrdered
plotCounts(dds, gene = "ENSG00000137693.14" , intgroup=c("Sample"))

#Check gene expression of ChromHMM transition states
#Import data

#H3K4me3 gaining H3K27ac
primet_nearest_genes <- read.table(file = "HN137Pri_HN137Met_E1_E2_regions_nearest_gene.bed", header = F)
head(primet_nearest_genes)
length(unique(primet_nearest_genes$V8))
genes <- unique(primet_nearest_genes$V8)
gene_index <- which(genes %in% resOrdered$symbol)
length(gene_index)
genes <- genes[gene_index]
rnaseq_gene_index <- which(resOrdered$symbol %in% genes)
length(rnaseq_gene_index)
res_subset <- resOrdered[rnaseq_gene_index,]
summary(res_subset$log2FoldChange)
E2_E3 <- data.frame(transition="E2_E3", log2FoldChange=res_subset$log2FoldChange)

#Unmodified gaining H3K27ac
primet_nearest_genes <- read.table(file = "HN137Pri_HN137Met_E4_E3_regions_nearest_gene.bed", header = F)
head(primet_nearest_genes)
length(unique(primet_nearest_genes$V8))
genes <- unique(primet_nearest_genes$V8)
gene_index <- which(genes %in% resOrdered$symbol)
length(gene_index)
genes <- genes[gene_index]
rnaseq_gene_index <- which(resOrdered$symbol %in% genes)
length(rnaseq_gene_index)
res_subset <- resOrdered[rnaseq_gene_index,]
summary(res_subset$log2FoldChange)
E3_E4 <- data.frame(transition="E3_E4", log2FoldChange=res_subset$log2FoldChange)

#Combine data
trans_change <- rbind(E2_E3, E3_E4)
ggplot(trans_change, aes(x=factor(transition), y=log2FoldChange, fill = factor(transition))) +
       geom_boxplot()

save.image("RNAseq_DESeq_analysis_16112021.RData")
