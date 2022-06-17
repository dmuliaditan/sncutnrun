#These series of code (3.4. - 3.9.) describes the analysis to correlate gene expression, gene copy number and gene chromatin state
#Previously gene expression data was retrieved using RNAseq analysis (see section on RNAseq analysis)
#CNV data was retrieved with CNVkit (see section on CNV calling)
#Chromatin state data was achieved with ChromHMM (see section on ChromHMM)

#3.4. This section describes the static CNV-RNAseq-Zscore correlation in the HN137Met cell line
#Version 10/06/2022
#Daniel Muliaditan

#First, set some colourblind palettes for plotting and set the working directory
setwd("D:/snCUT_RUN/scripts")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#AA4499", "#332288",  
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
#Load RNAseq z-scores
rna <- read.table(file = "D:/snCUT_RUN/results/rnaseq/DM_SNCUTRUN_RNAseq_HN120_HN137_Zscore_allgenes.txt", header = T)

#Load HN137MET CNV data
cnv <- read.table("D:/snCUT_RUN/wes/CNV/cnvkit/HN137MET.call.cns", header = T)

#Make a dataframe with the RNA Z-score and gene symbols
rna_cnv <- data.frame(GENE = rna$GENE,
                      RNAseq_Zscore=rna$HN137MET)

#Use the CNV data to find the copy number of the genes
rna_cnv$gene_CN = "NA"
for (j in seq_along(rna_cnv$GENE)) {
  print(paste(rna_cnv$GENE[j]))
  segind <- which(grepl(pattern = paste0(",", rna_cnv$GENE[j], ","), x = cnv$gene) == T)
  rna_cnv$gene_CN[j] <- ifelse(length(segind) == 1, yes = cnv$cn[segind], no = NA)
}

dim(rna_cnv)

#Remove rows with incomplete CNV data
rna_cnv_complete <- rna_cnv[complete.cases(rna_cnv),]
dim(rna_cnv_complete)

#Remove low prevalent copy numbers
CNind <- c("0", names(which(table(rna_cnv_complete$gene_CN) < 10)))
zeroCNind <- which(rna_cnv_complete$gene_CN %in% CNind)
rna_cnv_complete <- rna_cnv_complete[-zeroCNind,]
rna_cnv_complete$gene_CN <- as.factor(rna_cnv_complete$gene_CN)

#Plot CN RNA correlation
library(ggplot2)
ggplot2::ggplot(data = rna_cnv_complete, mapping = aes(x=reorder(factor(gene_CN),RNAseq_Zscore,FUN=median), y=RNAseq_Zscore,  fill = gene_CN)) +
  geom_boxplot() +
  theme_bw() +
  ylab("RNAseq Z score") +
  xlab("Gene Copy Number") +
  scale_fill_manual(values = safe_colorblind_palette[7:11]) +
  theme(legend.position = "none",
        axis.text = element_text(size = 30),
        axis.title = element_text(size = 40))

save.image(file = "25052022_RNAseq_Zscore_CNV_correlation.RData")
