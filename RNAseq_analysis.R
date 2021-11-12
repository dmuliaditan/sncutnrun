#RNA-seq analysis with R
#Tutorial RNA-seq workflow: gene-level exploratory analysis and differential expression, Michael I. Love, Simon Anders, Vladislav Kim and Wolfgang Huber, 16 October, 2019
#https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

setwd(/mnt/raid5/cutnrun/rnaseq)
dir <- /mnt/raid5/cutnrun/rnaseq/salmon_output

list.files(dir)

list.files(file.path(dir, "quants"))

csvfile <- file.path(dir, "sample_meta.csv")
coldata <- read.csv(csvfile, stringsAsFactors=FALSE, header=T)
coldata

coldata <- coldata[c(1:2,4),]
coldata$names <- coldata$Library
coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf.gz")
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

assayNames(gse)
head(assay(gse), 3)
colSums(assay(gse))
rowRanges(gse)
seqinfo(rowRanges(gse))
colData(gse)



