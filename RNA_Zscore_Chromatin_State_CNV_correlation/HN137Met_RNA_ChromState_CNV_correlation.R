#Static Chromatin State-CNV-RNAseq-Zscore correlation in the HN137Met cell line

#Load existing dataset and required packages
setwd("D:/snCUT_RUN/scripts")
load(file = "25052022_RNAseq_Zscore_chromatin_state_static_correlation.RData")

library(ggplot2)

#Gene CN, chromatin state, RNA Z-score correlation
rna_epi_cnv <- rna_epi
rna_epi_cnv$gene_CN = "NA"
for (j in seq_along(rna_epi_cnv$GENE)) {
  print(rna_epi_cnv$GENE[j])
  segind <- which(grepl(pattern = paste0(",", rna_epi_cnv$GENE[j], ","), x = cnv$gene) == T)
  rna_epi_cnv$gene_CN[j] <- ifelse(length(segind) == 1, yes = cnv$cn[segind], no = NA)
}

dim(rna_epi_cnv)
rna_epi_cnv <- rna_epi_cnv[complete.cases(rna_epi_cnv),]
dim(rna_epi_cnv)
zeroCNind <- which(rna_epi_cnv$gene_CN == 0)
rna_epi_cnv <- rna_epi_cnv[-zeroCNind,]
rna_epi_cnv$gene_CN <- as.numeric(rna_epi_cnv$gene_CN)
rna_epi_cnv$gene_CN <- ifelse(test = rna_epi_cnv$gene_CN >= 4, 
                                       yes = "4+", no = rna_epi_cnv$gene_CN)

ggplot2::ggplot(data = rna_epi_cnv) +
  geom_boxplot(mapping = aes(x=STATE, y=RNAseq_Zscore,  fill = factor(gene_CN))) +
  theme_bw() +
  ylab("RNAseq Z score") +
  xlab("Chromatin State") +
  labs(fill="Gene Copy Number") +
  scale_fill_manual(values = safe_colorblind_palette) +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 28),
        axis.text.y = element_text(size = 24),
        axis.text.x = element_text(size = 18),
        axis.title = element_text(size = 28))

save.image(file = "01062022_RNAseq_Zscore_chromatin_state_CNV_static_correlation.RData")
