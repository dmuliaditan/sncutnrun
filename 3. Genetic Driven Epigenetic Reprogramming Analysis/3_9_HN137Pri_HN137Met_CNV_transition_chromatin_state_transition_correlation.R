#3.9. HN137Pri, HN137Met State transition analysis, look at copy number affected regions

#Version: 10/06/2022
#Author: Daniel Muliaditan

#Here, we look at genes with a change in CN and investigate what are the accompanying chromatin state changes in these genes

#Load existing dataset and required packages
setwd("D:/snCUT_RUN/scripts")
load(file = "01062022_DownUpregulated_genes_HN137Pri_HN137Met_chromatin_state_transition_correlation.RData")

library(ggplot2)
library(stringr)

#Load HN137Pri and HN137Met CNV data
cnv_pri <- read.table("D:/snCUT_RUN/wes/CNV/cnvkit/HN137PRI_subsampled_sorted.call.cns", header = T)
cnv_met <- read.table("D:/snCUT_RUN/wes/CNV/cnvkit/HN137MET.call.cns", header = T)

#Create a gene list (in this case taken from the RNA-seq annotation)
genes <- data.frame(GENE = rna$GENE)

#Extract CNV data for the genes
genes$HN137Pri_geneCN = "NA"
genes$HN137Met_geneCN = "NA"
for (j in seq_along(genes$GENE)) {
  print(paste(genes$GENE[j]))
  segind <- which(grepl(pattern = paste0(",", genes$GENE[j], ","), x = cnv_pri$gene) == T)
  genes$HN137Pri_geneCN[j] <- ifelse(length(segind) == 1, yes = cnv_pri$cn[segind], no = NA)
  
  segind <- which(grepl(pattern = paste0(",", genes$GENE[j], ","), x = cnv_met$gene) == T)
  genes$HN137Met_geneCN[j] <- ifelse(length(segind) == 1, yes = cnv_met$cn[segind], no = NA)
}

genes_cnv <- genes #backup data frame
dim(genes_cnv)
genes_cnv_complete <- genes_cnv[complete.cases(genes_cnv),] #Remove genes without CN annotation
dim(genes_cnv_complete)

zeroCNind <- which(genes_cnv_complete$HN137Pri_geneCN == 0 | 
                     genes_cnv_complete$HN137Met_geneCN == 0 ) #Remove genes with 0 CN
genes_cnv_complete <- genes_cnv_complete[-zeroCNind,]

genes_cnv_complete$HN137Pri_geneCN <- as.numeric(genes_cnv_complete$HN137Pri_geneCN)
genes_cnv_complete$HN137Met_geneCN <- as.numeric(genes_cnv_complete$HN137Met_geneCN)

#Recategorize high copy number genes
genes_cnv_complete$HN137Pri_geneCN <- ifelse(test = genes_cnv_complete$HN137Pri_geneCN >= 4, 
                                             yes = "4+", no = genes_cnv_complete$HN137Pri_geneCN)

genes_cnv_complete$HN137Met_geneCN <- ifelse(test = genes_cnv_complete$HN137Met_geneCN >= 4, 
                                             yes = "4+", no = genes_cnv_complete$HN137Met_geneCN)


table(genes_cnv_complete$HN137Pri_geneCN)
table(genes_cnv_complete$HN137Met_geneCN)

#Annotate CN transition
genes_cnv_complete$CN_transition <- paste(genes_cnv_complete$HN137Pri_geneCN, 
                                          genes_cnv_complete$HN137Met_geneCN, sep = " > ")

cnv_epi_transition <- genes_cnv_complete #2nd backup data frame

#Add in chromatin state data
cnv_epi_transition <- cnv_epi_transition[which(cnv_epi_transition$GENE %in% genechrom$GENE),] #intersect genes occurring in the CNV data and epigenetic data
reorder_index <- match(x = cnv_epi_transition$GENE ,table = genechrom$GENE) #Reorder to match CNV gene order
chromstate_genes2 <- genechrom[reorder_index,]
cnv_epi_transition$EPISTATE <- factor(chromstate_genes2$STATE, levels = c("+H3K4me3", "+H3K27ac", "+H3K4me3/+H3K27ac","-H3K4me3", "-H3K27ac", "-H3K4me3/-H3K27ac",
                                                                          "No chromatin state change", "Ambiguous"))

table(cnv_epi_transition$EPISTATE)

#Remove genes with ambiguous chromatin state calling
cnv_epi_transition <- cnv_epi_transition[-which(cnv_epi_transition$EPISTATE == "Ambiguous"),]

cnv_epi_transition <- droplevels.data.frame(cnv_epi_transition)

#Remove the infrequent occurring transitions
low_freq_transition <- names(which(table(cnv_epi_transition$CN_transition) < 10))
cnv_epi_transition <- cnv_epi_transition[-which(cnv_epi_transition$CN_transition %in% low_freq_transition),] 
cnv_epi_transition <- droplevels.data.frame(cnv_epi_transition)

#Plot the CNV transitions stratified by chromatin state change
ggplot(data = cnv_epi_transition, mapping = aes(x = CN_transition, y= ..count../sum(..count..), fill=str_wrap(EPISTATE,20))) +
  geom_bar(position = "fill") +
  theme_bw() +
  ylab("Proportion of Genes") +
  xlab("Gene Copy Number Transition") +
  scale_fill_manual(values = cbPalette) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.height=unit(1.7, "cm"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 0.95, size = 22),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 15, b = 0, l = 0)))

cnv_epi_transition[cnv_epi_transition$EPISTATE == "+H3K4me3/+H3K27ac",]

#Reduce number of categories by merging some transitions
no_cn_change <- c("1 > 1", "2 > 2", "3 > 3", "4+ > 4+")
gain_1cn <- c("1 > 2", "2 > 3", "3 > 4+")
gain_2cn <- c("1 > 3", "1 > 4+")
lose_1cn <- c("2 > 1", "3 > 2")
lose_2cn <- c("3 > 1", "4+ > 1")

cnv_epi_transition$group <- 0
cnv_epi_transition$group <- ifelse(test = cnv_epi_transition$CN_transition %in% no_cn_change, yes = "No change in CN",
                                   no = ifelse(cnv_epi_transition$CN_transition %in% gain_1cn, yes = "Gain 1 CN", 
                                               no = ifelse(cnv_epi_transition$CN_transition %in% gain_2cn, yes = "Gain 2 CN",
                                                           no = ifelse(cnv_epi_transition$CN_transition %in% lose_1cn, yes = "Lose 1 CN",
                                                                       no = "Lose 2 CN"))))

#Plot again
ggplot(data = cnv_epi_transition, mapping = aes(x = group, y= ..count../sum(..count..), fill=str_wrap(EPISTATE,20))) +
  geom_bar(position = "fill") +
  theme_bw() +
  ylab("Proportion of Genes") +
  xlab("Gene Copy Number Transition") +
  scale_fill_manual(values = cbPalette) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.height=unit(1.45, "cm"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 0.95, size = 22),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 15, b = 0, l = 0)))

#Plot chromatin state, stratified by CN transition category
ggplot(data = cnv_epi_transition, mapping = aes(x = EPISTATE, y= ..count../sum(..count..), fill=group)) +
  geom_bar(position = "fill") +
  theme_bw() +
  ylab("Proportion of Genes") +
  xlab("") +
  scale_fill_manual(values = c(safe_colorblind_palette, cbPalette)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.height=unit(0.65, "cm"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 0.95, size = 18),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 15, b = 0, l = 0)))

save.image("01062022_HN137Pri_HN137Met_CNV_transition_chromatin_state_transition_correlation.RData")
