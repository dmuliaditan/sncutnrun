#4. Pri > Progressed chromatin state alluvial plots for remaining transitions

#Version: 19/07/2022
#Author: Daniel Muliaditan

#Load existing dataset and required packages
setwd("D:/snCUT_RUN/scripts")

library(networkD3)
library(dplyr)
library(tidyr)
library(ggalluvial)
library(rGREAT)


for (a in "HN120Pri") {
  print(a)

  for (b in c("HN120Met", "HN120PCR")) {
    print(paste0(a," - ",b))
    
    #Load ChromHMM final .bed files. See section on ChromHMM for the code to make this .bed files
    #It is good to note that the genome has been segmented to 200bp to equalize the bin size across the different samples
    bed1 <- read.table(file = paste0("D:/snCUT_RUN/results/ChromHMM/outputdir/5_states/",a,"_5_segments_binned_sorted.bed"), header = F)
    bed2 <- read.table(file = paste0("D:/snCUT_RUN/results/ChromHMM/outputdir/5_states/",b,"_5_segments_binned_sorted.bed"), header = F)
    
    #Make the data frames to combine the chromatin state annotation between HN137Pri and HN137Met
    #This data frame has the unique bin identifier and chromatin state change
    comparison <- data.frame(bin=paste(bed1$V1,bed1$V2,bed1$V3, sep = "-"),transition = paste0(bed1$V4,"-",bed2$V4))
    table <- table(comparison$transition)
    
    #Most regions do not transition, but rather unmodified chromatin that remained unmodified (E4-E4)
    
    #This is a more conventional .bed file with chr, start, end and transition columns
    #segment_combined <- data.frame(chr=bed1$V1,
     #                              start=bed1$V2,
      #                             end=bed1$V3,
       #                            state=paste0(bed1$V4,"-",bed2$V4))
    
    #Now, save all the various transition data per transition type as .bed file, except E4-E4 transition
    #transition <- c("E1-E1", "E1-E2", "E1-E3", "E1-E4", "E1-E5",
     #               "E2-E1", "E2-E2", "E2-E3", "E2-E4", "E2-E5",
      #              "E3-E1", "E3-E2", "E3-E3", "E3-E4", "E3-E5",
       #             "E4-E1", "E4-E2", "E4-E3", "E4-E5",
        #            "E5-E1", "E5-E2", "E5-E3", "E5-E4", "E5-E5")
    
   # for (j in transition) {
    #  print(j)
     # chrom_subset <- which(segment_combined$state == j)
      #subseg <- segment_combined[chrom_subset,]
      #write.table(x = subseg,
       #           file = paste0("D:/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_",j,"_regions.bed"), 
        #          quote = F,
         #         col.names = F,
          #        row.names = F,
           #       sep = "\t")
    #}
    
    #Plot alluvial plot
    #With the results from the table, make an excel sheet with the transitions
    #Repeat for HN120Pri-HN120PCR, HN137Pri-HN137Met, HN137Pri-HN137PCR
    
    table <- as.data.frame(table)
    links <- separate(data = table, col = Var1, into = c("source", "target"), sep = "\\-")
    colnames(links) <- c("source", "target", "value")
    links$source <- as.factor(links$source)
    links$target <- as.factor(links$target)
    links$value <- as.numeric(links$value)
    links <- links[-19,] #Remove E4-E4 transition
    links
    
    is_alluvia_form(as.data.frame(links), axes = 1:3, silent = TRUE)
    c <- ggplot(links, aes(y = value, axis1 = source, axis2 = target)) +
      geom_alluvium(aes(fill = source),width = 1/12) +
      geom_stratum(width = 1/12, color = "black", size = 1) +
      #geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
      scale_x_discrete(limits = c("source", "target"), expand = c(.05, .05)) +
      scale_fill_brewer(type = "qual", palette = "Set1") +
      ggtitle(paste0(a, " - ", b, " Chromatin State Transition")) +
      ylab("") +
      theme(title = element_text(size = 22),
            legend.title = element_blank(),
            legend.text = element_text(size = 24),
            axis.text = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank()) 
      
      ggsave(filename = paste0("D:/snCUT_RUN/figures/figures/supplementary_figs/",a, " - ", b, " Chromatin State Transition.svg"),
             device = "svg", plot = c, width = 10, height = 10, units = "cm", dpi = 300)
   
  }
}

for (a in c("HN137Pri")) {
  print(a)
  
  for (b in c("HN137PCR")) {
    print(b)
    
    #Load ChromHMM final .bed files. See section on ChromHMM for the code to make this .bed files
    #It is good to note that the genome has been segmented to 200bp to equalize the bin size across the different samples
    bed1 <- read.table(file = paste0("D:/snCUT_RUN/results/ChromHMM/outputdir/5_states/",a,"_5_segments_binned_sorted.bed"), header = F)
    bed2 <- read.table(file = paste0("D:/snCUT_RUN/results/ChromHMM/outputdir/5_states/",b,"_5_segments_binned_sorted.bed"), header = F)
    
    #Make the data frames to combine the chromatin state annotation between HN137Pri and HN137Met
    #This data frame has the unique bin identifier and chromatin state change
    comparison <- data.frame(bin=paste(bed1$V1,bed1$V2,bed1$V3, sep = "-"),transition = paste0(bed1$V4,"-",bed2$V4))
    table <- table(comparison$transition)
    #Most regions do not transition, but rather unmodified chromatin that remained unmodified (E4-E4)
    
    #This is a more conventional .bed file with chr, start, end and transition columns
    #segment_combined <- data.frame(chr=bed1$V1,
    #                              start=bed1$V2,
    #                             end=bed1$V3,
    #                            state=paste0(bed1$V4,"-",bed2$V4))
    
    #Now, save all the various transition data per transition type as .bed file, except E4-E4 transition
    #transition <- c("E1-E1", "E1-E2", "E1-E3", "E1-E4", "E1-E5",
    #               "E2-E1", "E2-E2", "E2-E3", "E2-E4", "E2-E5",
    #              "E3-E1", "E3-E2", "E3-E3", "E3-E4", "E3-E5",
    #             "E4-E1", "E4-E2", "E4-E3", "E4-E5",
    #            "E5-E1", "E5-E2", "E5-E3", "E5-E4", "E5-E5")
    
    # for (j in transition) {
    #  print(j)
    # chrom_subset <- which(segment_combined$state == j)
    #subseg <- segment_combined[chrom_subset,]
    #write.table(x = subseg,
    #           file = paste0("D:/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_",j,"_regions.bed"), 
    #          quote = F,
    #         col.names = F,
    #        row.names = F,
    #       sep = "\t")
    #}
    
    #Plot alluvial plot
    #With the results from the table, make an excel sheet with the transitions
    #Repeat for HN120Pri-HN120PCR, HN137Pri-HN137Met, HN137Pri-HN137PCR
    
    table <- as.data.frame(table)
    links <- separate(data = table, col = Var1, into = c("source", "target"), sep = "\\-")
    colnames(links) <- c("source", "target", "value")
    links$source <- as.factor(links$source)
    links$target <- as.factor(links$target)
    links$value <- as.numeric(links$value)
    links <- links[-19,] #Remove E4-E4 transition
    links
    
    is_alluvia_form(as.data.frame(links), axes = 1:3, silent = TRUE)
    c <- ggplot(links, aes(y = value, axis1 = source, axis2 = target)) +
      geom_alluvium(aes(fill = source),width = 1/12) +
      geom_stratum(width = 1/12, color = "black", size = 1) +
      #geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
      scale_x_discrete(limits = c("source", "target"), expand = c(.05, .05)) +
      scale_fill_brewer(type = "qual", palette = "Set1") +
      ggtitle(paste0(a, " - ", b, " Chromatin State Transition")) +
      ylab("") +
      theme(title = element_text(size = 22),
                      legend.title = element_blank(),
                      legend.text = element_text(size = 24),
                      axis.text = element_blank(),
                      panel.border = element_blank(),
                      panel.grid = element_blank(),
                      axis.line = element_blank(),
                      axis.ticks = element_blank()) 
      ggsave(filename = paste0("D:/snCUT_RUN/figures/figures/supplementary_figs/",a, " - ", b, " Chromatin State Transition.svg"),
                     device = "svg", plot = c, width = 10, height = 10, units = "cm", dpi = 300)
              
              
  }
  
}








#Looking at pathways, for example the regions gaining H3K27ac from H3K4me3 alone
chrom_subset <- which(segment_combined$state == "E1-E2")
subseg <- segment_combined[chrom_subset,]

#rGREAT analysis
job <- submitGreatJob(gr = subseg, species = 'hg38', rule = 'twoClosest',
                      adv_twoDistance = 2000)
tb <- getEnrichmentTables(job, category = "GO")
names(tb)
head(tb$`GO Biological Process`, n=10) #Look at affected pathways
availableOntologies(job)
res <- plotRegionGeneAssociationGraphs(job) #Most regions gaining H3K27ac from H3K4me3 alone seems to be concentrated near the promoter, which is expected
