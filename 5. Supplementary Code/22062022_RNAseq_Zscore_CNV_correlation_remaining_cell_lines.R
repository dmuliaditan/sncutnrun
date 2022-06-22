#This code describes the correlation between CN-profile and RNA-seq Z-score in the remaining cell lines other than HN137Met

#Version 22/06/2022
#Daniel Muliaditan

#First, set some colourblind palettes for plotting and set the working directory
setwd("D:/snCUT_RUN/scripts")

library(ggplot2)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#AA4499", "#332288",  
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

rna <- read.table(file = "D:/snCUT_RUN/results/rnaseq/DM_SNCUTRUN_RNAseq_HN120_HN137_Zscore_allgenes.txt", header = T)

for (a in c("HN120PRI", "HN120MET", "HN120PCR", "HN137PRI", "HN137PCR")) {
  print(a)
  #Load RNAseq z-scores
  
  #Load HN137MET CNV data
  cnv <- read.table(paste0("D:/snCUT_RUN/wes/CNV/cnvkit/",a,"_subsampled_sorted.call.cns"), header = T)
  
  #Make a dataframe with the RNA Z-score and gene symbols
  sample_ind <- which(colnames(rna) == a)
  rna_cnv <- data.frame(GENE = rna$GENE,
                        RNAseq_Zscore=rna[,sample_ind])
  
  #Use the CNV data to find the copy number of the genes
  rna_cnv$gene_CN = "NA"
  for (b in seq_along(rna_cnv$GENE)) {
    print(paste(rna_cnv$GENE[b]))
    segind <- which(grepl(pattern = paste0(",", rna_cnv$GENE[b], ","), x = cnv$gene) == T)
    rna_cnv$gene_CN[b] <- ifelse(length(segind) == 1, yes = cnv$cn[segind], no = NA)
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
  rna_cnv_complete <- droplevels.data.frame(rna_cnv_complete)
  
  #Plot CN RNA correlation
  
  p <- ggplot2::ggplot(data = rna_cnv_complete, mapping = aes(x=reorder(factor(gene_CN),RNAseq_Zscore,FUN=median), y=RNAseq_Zscore,  fill = gene_CN)) +
    geom_boxplot() +
    theme_bw() +
    ylab("RNAseq Z score") +
    xlab("Gene Copy Number") +
    scale_fill_manual(values = safe_colorblind_palette) +
    theme(legend.position = "none",
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 15))
  ggsave(filename = paste0("D:/snCUT_RUN/figures/figures/supplementary_figs/",a,"_geneCN_RNAseq_Zscore_correlation.png"),
         plot = p, width = 1024, height = 1024, units = "px", dpi = 300)
  
   
}

save.image(file = "22062022_RNAseq_Zscore_CNV_correlation_remaining_cell_lines.RData")
