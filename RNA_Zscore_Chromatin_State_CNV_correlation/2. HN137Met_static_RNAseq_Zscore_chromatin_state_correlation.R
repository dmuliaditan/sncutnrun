#2. Static Chromatin State-RNAseq-Zscore correlation in the HN137Met cell line

#Version: 10/06/2022
#Author: Daniel Muliaditan

#After correlating gene expression and copy number, now correlate gene expression and chromatin state,
#more specifically in the presence and absence of activating marks

#Load existing dataset and required packages
setwd("D:/snCUT_RUN/scripts")
load(file = "25052022_RNAseq_Zscore_CNV_correlation.RData")

library(ggpubr)
library(stringr)
library(ggplot2)

#Extract chromatin state data from ChromHMM results
#HN137Pri --> HN137Met transition data was used to extract HN137Met annotation states
#The second state (e.g. in E1_E2 it will be E2), is the state in HN137Met
states <- c("E1_E1", "E1_E2", "E1_E3", "E1_E4", "E1_E5",
            "E2_E1", "E2_E2", "E2_E3", "E2_E4", "E2_E5",
            "E3_E1", "E3_E2", "E3_E3", "E3_E4", "E3_E5",
            "E4_E1", "E4_E2", "E4_E3", "E4_E5",
            "E5_E1", "E5_E2", "E5_E3", "E5_E4", "E5_E5")

#Import all the chromatin state data
epistate <- list()
for (k in seq_along(states)) {
  print(paste(states[k]))
  epistate[[k]] <- read.table(paste0("D:/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_",states[k],"_regions_nearest_gene.bed"), header = F)
}
names(epistate) <- states
epistate[[1]]

#HN137Met RNAseq Zscore - chromatin state correlation
#Annotate HN137Met epigenetic states
#H3K27me3+ positive genes were removed due to low numbers

E1_genes <- data.frame(GENE=rbind(epistate[[1]],
                                  epistate[[6]],
                                  epistate[[11]],
                                  epistate[[16]],
                                  epistate[[20]]),
                       STATE="H3K4me3+/H3K27ac-")
E2_genes <- data.frame(GENE=rbind(epistate[[2]],
                                  epistate[[7]],
                                  epistate[[12]],
                                  epistate[[17]],
                                  epistate[[21]]),
                       STATE="H3K4me3+/H3K27ac+")
E3_genes <- data.frame(GENE=rbind(epistate[[3]],
                                  epistate[[8]],
                                  epistate[[13]],
                                  epistate[[18]],
                                  epistate[[22]]),
                       STATE="H3K4me3-/H3K27ac+")
E4_genes <- data.frame(GENE=rbind(epistate[[4]],
                                  epistate[[9]],
                                  epistate[[14]],
                                  epistate[[23]]),
                       STATE="H3K4me3-/H3K27ac-")

#Combine into a big data.frame
chromstate <- rbind(E1_genes, E2_genes, E3_genes, E4_genes)
chromstate <- chromstate[,c(8,10)]
colnames(chromstate) <- c("GENE", "STATE")
table(chromstate$STATE)
length(unique(chromstate$GENE))

#Remove genes with NA
chromstate <- chromstate[which(is.na(chromstate$GENE) == F),]
head(chromstate, n = 50)

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

#Retrieve RNA data
rna2 <- data.frame(GENE = rna$GENE,
                           RNAseq_Zscore=rna$HN137MET)

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
ggplot2::ggplot(data = rna_epi, aes(x=STATE, y=RNAseq_Zscore, fill = str_wrap(STATE,20))) +
  geom_boxplot() +
  theme_bw() +
  ylab("RNAseq Z score") +
  xlab("Epigenetic State") +
  scale_fill_manual(values = safe_colorblind_palette[8:12]) +
  coord_cartesian(ylim = c(-2.5, 4.2), expand = T) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.height=unit(1.8, "cm"),
        axis.text.y = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 28)) +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test", size = 8, label.y = c(2.2,2.7,3.2,3.7), label = "p.signif")

save.image(file = "25052022_RNAseq_Zscore_chromatin_state_static_correlation.RData")
