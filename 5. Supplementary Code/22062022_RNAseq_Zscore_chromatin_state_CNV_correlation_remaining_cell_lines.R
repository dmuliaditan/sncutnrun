#3.5./3.6. Static Chromatin State-RNAseq-Zscore correlation in the remaining cell lines

#Version: 22/06/2022
#Author: Daniel Muliaditan

#After correlating gene expression and copy number, now correlate gene expression and chromatin state,
#more specifically in the presence and absence of activating marks
#After correlating gene-expression/CNV and gene-expression/chromatin state separately, 
#we also correlate gene-expression, CNV and chromatin state together

#Load existing dataset and required packages
setwd("D:/snCUT_RUN/scripts")

library(ggpubr)
library(stringr)
library(ggplot2)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#AA4499", "#332288",  
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

rna <- read.table(file = "D:/snCUT_RUN/results/rnaseq/DM_SNCUTRUN_RNAseq_HN120_HN137_Zscore_allgenes.txt", header = T)

for (a in c("HN120Pri", "HN120Met", "HN120PCR", "HN137Pri", "HN137PCR")) {
  print(a)
  
  #Load chromatin states
  chromstate <- read.table(paste0("D:/snCUT_RUN/results/ChromHMM/",a,"_binned_promoters_nearest_gene.bed"))
  
  #Combine into a big data.frame
  chromstate <- chromstate[,c(8,4)]
  colnames(chromstate) <- c("GENE", "STATE")
  table(chromstate$STATE)
  length(unique(chromstate$GENE))
  
  #Remove genes with NA
  chromstate <- chromstate[which(is.na(chromstate$GENE) == F),]
  head(chromstate, n = 50)
  #Remove H3K27me3 genes
  E5 <- which(chromstate$STATE == "E5")
  chromstate <- chromstate[-E5,]
  chromstate <- droplevels.data.frame(chromstate)
  
  #Reannotate chromatin state
  chromstate$STATE <- ifelse(test = chromstate$STATE == "E1", yes = "H3K4me3+/H3K27ac-", 
                             no = ifelse(test = chromstate$STATE == "E2", yes = "H3K4me3+/H3K27ac+",
                                         no = ifelse(test = chromstate$STATE == "E3", yes = "H3K4me3-/H3K27ac+",
                                                     no = "H3K4me3-/H3K27ac-")))
  head(chromstate, n = 50)
  table(chromstate$STATE)
  
  #Annotate each gene with a chromatin state
  #Some genes have multiple gene annotations due to the high resolution of ChromHMM annotation (200bp)
  #Due to the data being concentrated at the promoter regions (5000bp), a gene would typically have a maximum of 5000/200 = 25 bins
  #For a gene to be classified as being positive for a specific mark, it needs to have at least 2 bins (400bp) with an activating mark
  
  #Create an empty list
  genechrom <- NULL
  
  #Loop through all genes
  for (j in unique(chromstate$GENE)) {
    print(j)
    gene_index <- chromstate[chromstate$GENE == j,]
    
    #In case of one unique index
    if (length(unique(gene_index$STATE)) == 1) {
      genechrom <- rbind(genechrom, gene_index[1,])
      next
    }
    
    #In case of multiple chromatin states
    statetable <- as.data.frame(table(gene_index$STATE))
    activesum <- sum(statetable$Freq[statetable$Var1 != "H3K4me3-/H3K27ac-"])
    
    #If there are more than 2 bins with at least one activating mark, assume the gene to be active, if not, categorize it as non-active
    if (activesum < 2) {
      print("This gene has no activating marks!")
      finalannot <- data.frame(GENE=j, STATE= "H3K4me3-/H3K27ac-")
      genechrom <- rbind(genechrom, finalannot)
      next
    }
    
    if (activesum >= 2) {
      print("This gene has activating marks!")
      activetable <- statetable[which(statetable$Var1 != "H3K4me3-/H3K27ac-"),]
      
      
      #If from the activating states there is only one kind, annotate the gene as having that chromatin state
      if (nrow(activetable) == 1) {
        finalannot <- data.frame(GENE=j, STATE= activetable$Var1)
        genechrom <- rbind(genechrom, finalannot)
        next
      }
      
      #If a gene has multiple annotation states, annotate accordingly
      if (nrow(activetable) > 1) {
        
        #If more than 20% of activating bins is H3K4me3+/H3K27ac+, annotate the gene as having both marks
        if ("H3K4me3+/H3K27ac+" %in% unique(activetable$Var1) == T) {
          
          if ((activetable$Freq[activetable$Var1 == "H3K4me3+/H3K27ac+"]) / sum(activetable$Freq) >= 0.2) {
            finalannot <- data.frame(GENE=j, STATE= "H3K4me3+/H3K27ac+")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          
          #If less than 20% of the bins is H3K4me3+/H3K27ac+, compare the ratio between H3K4me3 and H3K27ac
          if ((activetable$Freq[activetable$Var1 == "H3K4me3+/H3K27ac+"]) / sum(activetable$Freq) < 0.2) {
            resttable <- activetable[which(activetable$Var1 != "H3K4me3+/H3K27ac+"),]
            
            #If all bins are the same annotation
            if (nrow(resttable) == 1) {
              finalannot <- data.frame(GENE=j, STATE= resttable$Var1)
              genechrom <- rbind(genechrom, finalannot)
              next
              
            }
            
            
            #If more than 60% of bins is annotated with H3K4me3+/H3K27ac-, categorize as H3K4me3+/H3K27ac-
            if (resttable$Freq[resttable$Var1 == "H3K4me3+/H3K27ac-"] / sum(resttable$Freq) > 0.6 &
                resttable$Freq[resttable$Var1 == "H3K4me3-/H3K27ac+"] / sum(resttable$Freq) < 0.4) {
              
              finalannot <- data.frame(GENE=j, STATE= "H3K4me3+/H3K27ac-")
              genechrom <- rbind(genechrom, finalannot)
              next
              
            }
            
            #If more than 60% of bins is annotated with H3K4me3-/H3K27ac+, categorize as H3K4me3-/H3K27ac+
            if (resttable$Freq[resttable$Var1 == "H3K4me3-/H3K27ac+"] / sum(resttable$Freq) > 0.6 &
                resttable$Freq[resttable$Var1 == "H3K4me3+/H3K27ac-"] / sum(resttable$Freq) < 0.4) {
              
              finalannot <- data.frame(GENE=j, STATE= "H3K4me3-/H3K27ac+")
              genechrom <- rbind(genechrom, finalannot)
              next
              
            }
            
            
            #If ratio of H3K4me3+ bins and H3K27ac+ bins are similar, categorize as H3K4me3+/H3K27ac+
            if (resttable$Freq[resttable$Var1 == "H3K4me3+/H3K27ac-"] / sum(resttable$Freq) > 0.4 &
                resttable$Freq[resttable$Var1 == "H3K4me3-/H3K27ac+"] / sum(resttable$Freq) > 0.4) {
              
              finalannot <- data.frame(GENE=j, STATE= "H3K4me3+/H3K27ac+")
              genechrom <- rbind(genechrom, finalannot)
              next
              
            }
            
            
          }
        }
        
        if ("H3K4me3+/H3K27ac+" %in% unique(activetable$Var1) == F) {
          resttable <- activetable[which(activetable$Var1 != "H3K4me3+/H3K27ac+"),]
          
          #If all bins are the same annotation
          if (nrow(resttable) == 1) {
            finalannot <- data.frame(GENE=j, STATE= resttable$Var1)
            genechrom <- rbind(genechrom, finalannot)
            next
            
          }
          
          #If more than 60% of bins is annotated with H3K4me3+/H3K27ac-, categorize as H3K4me3+/H3K27ac-
          if (resttable$Freq[resttable$Var1 == "H3K4me3+/H3K27ac-"] / sum(resttable$Freq) > 0.6 &
              resttable$Freq[resttable$Var1 == "H3K4me3-/H3K27ac+"] / sum(resttable$Freq) < 0.4) {
            
            finalannot <- data.frame(GENE=j, STATE= "H3K4me3+/H3K27ac-")
            genechrom <- rbind(genechrom, finalannot)
            next
            
          }
          
          #If more than 60% of bins is annotated with H3K4me3-/H3K27ac+, categorize as H3K4me3-/H3K27ac+
          if (resttable$Freq[resttable$Var1 == "H3K4me3-/H3K27ac+"] / sum(resttable$Freq) > 0.6 &
              resttable$Freq[resttable$Var1 == "H3K4me3+/H3K27ac-"] / sum(resttable$Freq) < 0.4) {
            
            finalannot <- data.frame(GENE=j, STATE= "H3K4me3-/H3K27ac+")
            genechrom <- rbind(genechrom, finalannot)
            next
            
          }
          
          #If ratio of H3K4me3+ bins and H3K27ac+ bins are similar, categorize as H3K4me3+/H3K27ac+
          if (resttable$Freq[resttable$Var1 == "H3K4me3+/H3K27ac-"] / sum(resttable$Freq) > 0.4 &
              resttable$Freq[resttable$Var1 == "H3K4me3-/H3K27ac+"] / sum(resttable$Freq) > 0.4) {
            
            finalannot <- data.frame(GENE=j, STATE= "H3K4me3+/H3K27ac+")
            genechrom <- rbind(genechrom, finalannot)
            next
            
          }
          
          
        }
        
      }
      next
    }
    
    prevalent_state <- "Ambiguous"
    prevalent_index <- data.frame(GENE=j, STATE= prevalent_state)
    genechrom <- rbind(genechrom, prevalent_index)
    next
    
    
  }
  
  table(genechrom$STATE)
  
  b <- toupper(a)
  
  #Retrieve RNA data
  sample_ind <- which(colnames(rna) == b)
  rna2 <- data.frame(GENE = rna$GENE,
                        RNAseq_Zscore=rna[,sample_ind])
  
  rna_epi <- rna2[which(rna2$GENE %in% genechrom$GENE),]
  
  
  #Reorder the RNA en EPI gene names so they match
  head(rna_epi)
  head(genechrom)
  reorder_index <- match(x = rna_epi$GENE,table = genechrom$GENE)
  genechrom2 <- genechrom[reorder_index,]
  head(genechrom2)
  
  #Add chromatin state data
  rna_epi$STATE <- factor(genechrom2$STATE, levels = c("H3K4me3+/H3K27ac+",
                                                       "H3K4me3-/H3K27ac+",
                                                       "H3K4me3+/H3K27ac-", 
                                                       "H3K4me3-/H3K27ac-"))
  
  
  #Make a list of comparisons for statistical calculation
  my_comparisons <- list( c("H3K4me3+/H3K27ac+", "H3K4me3-/H3K27ac+"), 
                          c("H3K4me3-/H3K27ac+", "H3K4me3+/H3K27ac-"), 
                          c("H3K4me3+/H3K27ac+", "H3K4me3-/H3K27ac-"),
                          c("H3K4me3+/H3K27ac-", "H3K4me3-/H3K27ac-"))
  
  #Plot the image
  p <- ggplot2::ggplot(data = rna_epi, aes(x=STATE, y=RNAseq_Zscore, fill = str_wrap(STATE,20))) +
    geom_boxplot() +
    theme_bw() +
    ylab("RNAseq Z score") +
    xlab("Epigenetic State") +
    scale_fill_manual(values = safe_colorblind_palette[8:12]) +
    coord_cartesian(ylim = c(-2.5, 4.2), expand = T) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 10),
          legend.key.height=unit(1.8, "cm"),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 14)) +
    stat_compare_means(comparisons = my_comparisons, 
                       method = "wilcox.test", size = 6, label.y = c(2.2,2.7,3.2,3.7), label = "p.signif")
  ggsave(filename = paste0("D:/snCUT_RUN/figures/figures/supplementary_figs/",a,"_chromatin_state_RNAseq_Zscore_correlation.png"),
         plot = p, width = 2048, height = 1024, units = "px", dpi = 300)
  
  
  #RNAseq Z-score, chromatin state, CNV correlation
  rna_epi_cnv <- rna_epi
  
  cnv <- read.table(paste0("D:/snCUT_RUN/wes/CNV/cnvkit/",a,"_subsampled_sorted.call.cns"), header = T)
  rna_epi_cnv$gene_CN = "NA"
  
  for (j in seq_along(rna_epi_cnv$GENE)) {
    print(rna_epi_cnv$GENE[j])
    segind <- which(grepl(pattern = paste0(",", rna_epi_cnv$GENE[j], ","), x = cnv$gene) == T)
    rna_epi_cnv$gene_CN[j] <- ifelse(length(segind) == 1, yes = cnv$cn[segind], no = NA)
  }
  
  dim(rna_epi_cnv)
  rna_epi_cnv <- rna_epi_cnv[complete.cases(rna_epi_cnv),] #Remove genes with NA CN
  dim(rna_epi_cnv)
  
  #Remove genes with 0 CN
  zeroCNind <- which(rna_epi_cnv$gene_CN == 0)
  rna_epi_cnv <- rna_epi_cnv[-zeroCNind,]
  rna_epi_cnv$gene_CN <- as.numeric(rna_epi_cnv$gene_CN)
  
  #Recategorize genes with higher CN as 4+ CN
  rna_epi_cnv$gene_CN <- ifelse(test = rna_epi_cnv$gene_CN >= 4, 
                                yes = "4+", no = rna_epi_cnv$gene_CN)
  
  #Plot the figure
  q <- ggplot2::ggplot(data = rna_epi_cnv) +
    geom_boxplot(mapping = aes(x=STATE, y=RNAseq_Zscore,  fill = factor(gene_CN)), lwd=0.3, outlier.stroke = 0.3) +
    theme_bw() +
    ylab("RNAseq Z score") +
    xlab("Chromatin State") +
    labs(fill="Gene Copy Number") +
    scale_fill_manual(values = safe_colorblind_palette) +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 7.5),
          axis.title = element_text(size = 12))
  
  ggsave(filename = paste0("D:/snCUT_RUN/figures/figures/supplementary_figs/",a,"_chromatin_state_CNV_RNAseq_Zscore_correlation.png"),
         plot = q, width = 2048, height = 1024, units = "px", dpi = 300)
  
  
}


#save.image(file = "22062022_RNAseq_Zscore_chromatin_state_static_correlation_other_cell_lines.RData")
