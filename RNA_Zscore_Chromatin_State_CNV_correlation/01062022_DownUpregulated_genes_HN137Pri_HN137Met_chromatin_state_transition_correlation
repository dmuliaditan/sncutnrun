##### HN137Pri, HN137Met State transition analysis

#Load existing dataset and required packages
setwd("D:/snCUT_RUN/scripts")
load(file = "01062022_RNAseq_Zscore_chromatin_state_CNV_static_correlation.RData")

library(ggplot2)

#Define genes gaining or losing H3K4me3 and H3K27ac
gain_k4me3 <- data.frame(GENE=rbind(epistate[[11]], epistate[[12]],epistate[[16]],epistate[[20]]), STATE="+H3K4me3")
gain_k27ac <- data.frame(GENE=rbind(epistate[[2]], epistate[[3]],epistate[[18]],epistate[[22]]), STATE="+H3K27ac")
gain_both <- data.frame(GENE=rbind(epistate[[17]],epistate[[21]]), STATE="+H3K4me3/+H3K27ac")
lose_k4me3 <- data.frame(GENE=rbind(epistate[[4]],epistate[[8]]), STATE="-H3K4me3")
lose_k27ac <- data.frame(GENE=rbind(epistate[[6]],epistate[[14]]), STATE="-H3K27ac")
lose_both <- data.frame(GENE=rbind(epistate[[9]],epistate[[10]]), STATE="-H3K4me3/-H3K27ac")
nochange <- data.frame(GENE=rbind(epistate[[1]],epistate[[7]],epistate[[13]],epistate[[24]]), STATE="No chromatin state change") 

#Combine all annotations
chromstate <- rbind(gain_k4me3, gain_k27ac, gain_both, lose_k4me3, lose_k27ac, lose_both, nochange)
chromstate <- chromstate[,c(8,10)]
colnames(chromstate) <- c("GENE", "STATE")
table(chromstate$STATE)
length(unique(chromstate$GENE))
chromstate <- chromstate[which(is.na(chromstate$GENE) == F),]

#An example of the annotation
chromstate[chromstate$GENE == "STAT1",]

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
  activesum <- sum(statetable$Freq[statetable$Var1 != "No chromatin state change"])
  
  #If there are more than 2 bins with at least one transition, assume the gene to be in transition, if not, categorize it as no-chromatin state change
  if (activesum < 2) {
    print("This gene is not transitioning")
    finalannot <- data.frame(GENE=j, STATE= "No chromatin state change")
    genechrom <- rbind(genechrom, finalannot)
    next
  }
  
  if (activesum >= 2) {
    print("This gene is transitioning")
    activetable <- statetable[which(statetable$Var1 != "No chromatin state change"),]
    
    
    #If from the transitioning states there is only one kind, annotate the gene as having that transition
    if (nrow(activetable) == 1) {
      finalannot <- data.frame(GENE=j, STATE= activetable$Var1)
      genechrom <- rbind(genechrom, finalannot)
      next
    }
    
    #If a gene has multiple transitions, annotate accordingly
    if (nrow(activetable) > 1) {
      
      #If more than 20% of transitioning bins are +H3K4me3/+H3K27ac, annotate the gene as gaining both marks
      if ("+H3K4me3/+H3K27ac" %in% unique(activetable$Var1) == T) {
        
        if ((activetable$Freq[activetable$Var1 == "+H3K4me3/+H3K27ac"]) / sum(activetable$Freq) >= 0.2) {
          finalannot <- data.frame(GENE=j, STATE= "+H3K4me3/+H3K27ac")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        
        #If less than 20% of the bins is +H3K4me3/+H3K27ac, compare the ratio between H3K4me3 and H3K27ac
        if ((activetable$Freq[activetable$Var1 == "+H3K4me3/+H3K27ac"]) / sum(activetable$Freq) < 0.2) {
          resttable <- activetable[which(activetable$Var1 != "+H3K4me3/+H3K27ac"),]
          
          #If all bins are the same annotation
          if (nrow(resttable) == 1) {
            finalannot <- data.frame(GENE=j, STATE= resttable$Var1)
            genechrom <- rbind(genechrom, finalannot)
            next
            
          }
          
          #If there are more than 2 states
          if (nrow(resttable) > 2) {
            finalannot <- data.frame(GENE=j, STATE= resttable$Var1[which.max(resttable$Freq)])
            genechrom <- rbind(genechrom, finalannot)
            next
            
          }
          
          #If only two states contested
          #In case of -H3K27ac and +H3K4me3
          if ("-H3K27ac" %in% unique(resttable$Var1) && "+H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K27ac"] / sum(resttable$Freq) > 0.4 &
                                                                                                   resttable$Freq[resttable$Var1 == "+H3K4me3"] / sum(resttable$Freq) > 0.4)) {
            
            finalannot <- data.frame(GENE=j, STATE= "Ambiguous")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          if ("-H3K27ac" %in% unique(resttable$Var1) && "+H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K27ac"] / sum(resttable$Freq) > 0.6)) {
            
            finalannot <- data.frame(GENE=j, STATE= "-H3K27ac")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          if ("-H3K27ac" %in% unique(resttable$Var1) && "+H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K4me3"] / sum(resttable$Freq) > 0.6)) {
            
            finalannot <- data.frame(GENE=j, STATE= "+H3K4me3")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          
          #In case of +H3K27ac and -H3K4me3
          if ("+H3K27ac" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K27ac"] / sum(resttable$Freq) > 0.4 &
                                                                                                   resttable$Freq[resttable$Var1 == "-H3K4me3"] / sum(resttable$Freq) > 0.4)) {
            
            finalannot <- data.frame(GENE=j, STATE= "Ambiguous")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          if ("+H3K27ac" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K27ac"] / sum(resttable$Freq) > 0.6)) {
            
            finalannot <- data.frame(GENE=j, STATE= "+H3K27ac")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          if ("+H3K27ac" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K4me3"] / sum(resttable$Freq) > 0.6)) {
            
            finalannot <- data.frame(GENE=j, STATE= "-H3K4me3")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          
          #In case of +H3K27ac and -H3K27ac
          if ("+H3K27ac" %in% unique(resttable$Var1) && "-H3K27ac" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K27ac"] / sum(resttable$Freq) > 0.4 &
                                                                                                   resttable$Freq[resttable$Var1 == "-H3K27ac"] / sum(resttable$Freq) > 0.4)) {
            
            finalannot <- data.frame(GENE=j, STATE= "Ambiguous")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          if ("+H3K27ac" %in% unique(resttable$Var1) && "-H3K27ac" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K27ac"] / sum(resttable$Freq) > 0.6)) {
            
            finalannot <- data.frame(GENE=j, STATE= "+H3K27ac")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          if ("+H3K27ac" %in% unique(resttable$Var1) && "-H3K27ac" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K27ac"] / sum(resttable$Freq) > 0.6)) {
            
            finalannot <- data.frame(GENE=j, STATE= "-H3K27ac")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          
          #In case of +H3K4me3 and -H3K4me3
          if ("+H3K4me3" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K4me3"] / sum(resttable$Freq) > 0.4 &
                                                                                                   resttable$Freq[resttable$Var1 == "-H3K4me3"] / sum(resttable$Freq) > 0.4)) {
            
            finalannot <- data.frame(GENE=j, STATE= "Ambiguous")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          if ("+H3K4me3" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K4me3"] / sum(resttable$Freq) > 0.6)) {
            
            finalannot <- data.frame(GENE=j, STATE= "+H3K4me3")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          if ("+H3K4me3" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K4me3"] / sum(resttable$Freq) > 0.6)) {
            
            finalannot <- data.frame(GENE=j, STATE= "-H3K4me3")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          
          #In case of +H3K27ac and +H3K4me3
          if ("+H3K27ac" %in% unique(resttable$Var1) && "+H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K27ac"] / sum(resttable$Freq) > 0.4 &
                                                                                                   resttable$Freq[resttable$Var1 == "+H3K4me3"] / sum(resttable$Freq) > 0.4)) {
            
            finalannot <- data.frame(GENE=j, STATE= "+H3K4me3/+H3K27ac")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          if ("+H3K27ac" %in% unique(resttable$Var1) && "+H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K27ac"] / sum(resttable$Freq) > 0.6)) {
            
            finalannot <- data.frame(GENE=j, STATE= "+H3K27ac")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          if ("+H3K27ac" %in% unique(resttable$Var1) && "+H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K4me3"] / sum(resttable$Freq) > 0.6)) {
            
            finalannot <- data.frame(GENE=j, STATE= "+H3K4me3")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          
          #In case of -H3K27ac and -H3K4me3
          if ("-H3K27ac" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K27ac"] / sum(resttable$Freq) > 0.4 &
                                                                                                   resttable$Freq[resttable$Var1 == "-H3K4me3"] / sum(resttable$Freq) > 0.4)) {
            
            finalannot <- data.frame(GENE=j, STATE= "-H3K4me3/-H3K27ac")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          if ("-H3K27ac" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K27ac"] / sum(resttable$Freq) > 0.6)) {
            
            finalannot <- data.frame(GENE=j, STATE= "-H3K27ac")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          if ("-H3K27ac" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K4me3"] / sum(resttable$Freq) > 0.6)) {
            
            finalannot <- data.frame(GENE=j, STATE= "-H3K4me3")
            genechrom <- rbind(genechrom, finalannot)
            next
          }
          
        }
      }
      
      if ("+H3K4me3/+H3K27ac" %in% unique(activetable$Var1) == F) {
        resttable <- activetable[which(activetable$Var1 != "+H3K4me3/+H3K27ac"),]
        
        #If all bins are the same annotation
        if (nrow(resttable) == 1) {
          finalannot <- data.frame(GENE=j, STATE= resttable$Var1)
          genechrom <- rbind(genechrom, finalannot)
          next
          
        }
        
        #If there are more than 2 states
        if (nrow(resttable) > 2) {
          finalannot <- data.frame(GENE=j, STATE= resttable$Var1[which.max(resttable$Freq)])
          genechrom <- rbind(genechrom, finalannot)
          next
          
        }
        
        #If only two states contested
        #In case of -H3K27ac and +H3K4me3
        if ("-H3K27ac" %in% unique(resttable$Var1) && "+H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K27ac"] / sum(resttable$Freq) > 0.4 &
                                                                                                 resttable$Freq[resttable$Var1 == "+H3K4me3"] / sum(resttable$Freq) > 0.4)) {
          
          finalannot <- data.frame(GENE=j, STATE= "Ambiguous")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        if ("-H3K27ac" %in% unique(resttable$Var1) && "+H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K27ac"] / sum(resttable$Freq) > 0.6)) {
          
          finalannot <- data.frame(GENE=j, STATE= "-H3K27ac")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        if ("-H3K27ac" %in% unique(resttable$Var1) && "+H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K4me3"] / sum(resttable$Freq) > 0.6)) {
          
          finalannot <- data.frame(GENE=j, STATE= "+H3K4me3")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        
        #In case of +H3K27ac and -H3K4me3
        if ("+H3K27ac" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K27ac"] / sum(resttable$Freq) > 0.4 &
                                                                                                 resttable$Freq[resttable$Var1 == "-H3K4me3"] / sum(resttable$Freq) > 0.4)) {
          
          finalannot <- data.frame(GENE=j, STATE= "Ambiguous")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        if ("+H3K27ac" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K27ac"] / sum(resttable$Freq) > 0.6)) {
          
          finalannot <- data.frame(GENE=j, STATE= "+H3K27ac")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        if ("+H3K27ac" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K4me3"] / sum(resttable$Freq) > 0.6)) {
          
          finalannot <- data.frame(GENE=j, STATE= "-H3K4me3")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        
        #In case of +H3K27ac and -H3K27ac
        if ("+H3K27ac" %in% unique(resttable$Var1) && "-H3K27ac" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K27ac"] / sum(resttable$Freq) > 0.4 &
                                                                                                 resttable$Freq[resttable$Var1 == "-H3K27ac"] / sum(resttable$Freq) > 0.4)) {
          
          finalannot <- data.frame(GENE=j, STATE= "Ambiguous")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        if ("+H3K27ac" %in% unique(resttable$Var1) && "-H3K27ac" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K27ac"] / sum(resttable$Freq) > 0.6)) {
          
          finalannot <- data.frame(GENE=j, STATE= "+H3K27ac")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        if ("+H3K27ac" %in% unique(resttable$Var1) && "-H3K27ac" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K27ac"] / sum(resttable$Freq) > 0.6)) {
          
          finalannot <- data.frame(GENE=j, STATE= "-H3K27ac")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        
        #In case of +H3K4me3 and -H3K4me3
        if ("+H3K4me3" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K4me3"] / sum(resttable$Freq) > 0.4 &
                                                                                                 resttable$Freq[resttable$Var1 == "-H3K4me3"] / sum(resttable$Freq) > 0.4)) {
          
          finalannot <- data.frame(GENE=j, STATE= "Ambiguous")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        if ("+H3K4me3" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K4me3"] / sum(resttable$Freq) > 0.6)) {
          
          finalannot <- data.frame(GENE=j, STATE= "+H3K4me3")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        if ("+H3K4me3" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K4me3"] / sum(resttable$Freq) > 0.6)) {
          
          finalannot <- data.frame(GENE=j, STATE= "-H3K4me3")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        
        #In case of +H3K27ac and +H3K4me3
        if ("+H3K27ac" %in% unique(resttable$Var1) && "+H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K27ac"] / sum(resttable$Freq) > 0.4 &
                                                                                                 resttable$Freq[resttable$Var1 == "+H3K4me3"] / sum(resttable$Freq) > 0.4)) {
          
          finalannot <- data.frame(GENE=j, STATE= "+H3K4me3/+H3K27ac")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        if ("+H3K27ac" %in% unique(resttable$Var1) && "+H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K27ac"] / sum(resttable$Freq) > 0.6)) {
          
          finalannot <- data.frame(GENE=j, STATE= "+H3K27ac")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        if ("+H3K27ac" %in% unique(resttable$Var1) && "+H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "+H3K4me3"] / sum(resttable$Freq) > 0.6)) {
          
          finalannot <- data.frame(GENE=j, STATE= "+H3K4me3")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        
        #In case of -H3K27ac and -H3K4me3
        if ("-H3K27ac" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K27ac"] / sum(resttable$Freq) > 0.4 &
                                                                                                 resttable$Freq[resttable$Var1 == "-H3K4me3"] / sum(resttable$Freq) > 0.4)) {
          
          finalannot <- data.frame(GENE=j, STATE= "-H3K4me3/-H3K27ac")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        if ("-H3K27ac" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K27ac"] / sum(resttable$Freq) > 0.6)) {
          
          finalannot <- data.frame(GENE=j, STATE= "-H3K27ac")
          genechrom <- rbind(genechrom, finalannot)
          next
        }
        if ("-H3K27ac" %in% unique(resttable$Var1) && "-H3K4me3" %in% unique(resttable$Var1) && (resttable$Freq[resttable$Var1 == "-H3K4me3"] / sum(resttable$Freq) > 0.6)) {
          
          finalannot <- data.frame(GENE=j, STATE= "-H3K4me3")
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

#Plot chromatin state - RNA-expression transition correlation
#Read in RNA-seq Top Differentially Regulated genes between HN137Pri and HN137Met
met_rna <- read.table("D:/snCUT_RUN/results/rnaseq/HN137Met_vs_HN137Pri_top_differentially regulated gene.txt", header = T)
head(met_rna)

#Subset only the top differentially expressed genes
met_rna <- met_rna[which(met_rna$padj < 0.001),] 
met_rna <- met_rna[which(met_rna$log2FoldChange > 1 | met_rna$log2FoldChange < -0.5),]
met_rna_noNA <- met_rna[which(is.na(met_rna$symbol) == F),]
met_rna_fin <- met_rna_noNA[which(duplicated(met_rna_noNA$symbol) == F),]

#Reorder chromatin transition data so it matches the RNAseq gene order
met_rna_epi <- met_rna_fin[which(met_rna_fin$symbol %in% genechrom$GENE),]
reorder_index <- match(x = met_rna_epi$symbol,table = genechrom$GENE)
chromstate_genes2 <- genechrom[reorder_index,]
head(met_rna_epi)
head(chromstate_genes2)

#Add the chromstate data to the RNA-seq Z-score data
met_rna_epi$EPISTATE <- factor(chromstate_genes2$STATE, levels = c("+H3K4me3", "+H3K27ac", "+H3K4me3/+H3K27ac", 
                                                                   "-H3K4me3", "-H3K27ac", "-H3K4me3/-H3K27ac",
                                                                   "No chromatin state change", "Ambiguous"))
#Categorize down- and upregulated genes
met_rna_epi$group <- 0
met_rna_epi$group <- ifelse(met_rna_epi$log2FoldChange > 1, yes = "Upregulated", no = "Downregulated")
met_rna_epi$group <- factor(met_rna_epi$group)
table(met_rna_epi$EPISTATE)
table(met_rna_epi$group)

#Drop genes with ambiguous chromatin states
met_rna_epi <- met_rna_epi[-which(met_rna_epi$EPISTATE == "Ambiguous"),]
met_rna_epi <- droplevels.data.frame(met_rna_epi)

#Plot down and upregulated genes stratified by chromatin state transition
ggplot(data = met_rna_epi, mapping = aes(x = group, y= ..count../sum(..count..), fill=str_wrap(EPISTATE,20))) +
  geom_bar(position = "fill") +
  theme_bw() +
  ylab("Proportion of Genes") +
  xlab("") +
  scale_fill_manual(values = cbPalette) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.height=unit(1.48, "cm"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 0.95, size = 22),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 15, b = 0, l = 0)))

#Plot chromatin state transition stratified by down and upregulated genes
ggplot(data = met_rna_epi, mapping = aes(x = EPISTATE, y= ..count../sum(..count..), fill=group)) +
  geom_bar(position = "fill") +
  theme_bw() +
  ylab("Proportion of Genes") +
  xlab("") +
  scale_fill_manual(values = cbPalette) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.key.height=unit(1.8, "cm"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 0.95, size = 18),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24, margin = margin(t = 0, r = 15, b = 0, l = 0)))

save.image(file = "01062022_DownUpregulated_genes_HN137Pri_HN137Met_chromatin_state_transition_correlation.RData")
