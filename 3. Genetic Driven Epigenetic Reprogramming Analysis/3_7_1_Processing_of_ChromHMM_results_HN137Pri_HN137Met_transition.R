#3.7.1. This section describes the script used to further process ChromHMM data in R
#Version 22/06/2022
#Daniel Muliaditan

#Required input:
#-ChromHMM output, binned to 200bp segments

#Load ChromHMM final .bed files. See section on ChromHMM for the code to make this .bed files
#In this case, HN137Pri and HN137Met data shall be combined as we will be looking into the chromatin state transition 
#between HN137Pri and HN137Met later on.
#It is good to note that the genome has been segmented to 200bp to equalize the bin size across the different samples
#Read in .bed files
bed1 <- read.table(file = "D:/snCUT_RUN/results/ChromHMM/outputdir/5_states/HN137Pri_5_segments_binned_sorted.bed", header = F)
bed2 <- read.table(file = "D:/snCUT_RUN/results/ChromHMM/outputdir/5_states/HN137Met_5_segments_binned_sorted.bed", header = F)

#Check the chromosome positions across the bed files are equal
nrow(bed1)
nrow(bed2)
head(bed1, n = 200)
head(bed2, n = 200)

#Make the data frames to combine the chromatin state annotation between HN137Pri and HN137Met
#This data frame has the unique bin identifier and chromatin state change
comparison <- data.frame(bin=paste(bed1$V1,bed1$V2,bed1$V3, sep = "-"),transition = paste0(bed1$V4,"-",bed2$V4))
table <- table(comparison$transition)
table 
#Most regions do not transition, but rather unmodified chromatin that remained unmodified (E4-E4)

#This is a more conventional .bed file with chr, start, end and transition columns
segment_combined <- data.frame(chr=bed1$V1,
                               start=bed1$V2,
                               end=bed1$V3,
                               state=paste0(bed1$V4,"-",bed2$V4))

#As an example look at regions of H3K4me3 that gain H3K27ac during HN137Pri-HN137Met transition
chrom_subset <- which(segment_combined$state == "E1-E2")
subseg <- segment_combined[chrom_subset,]
head(subseg)
#13060 bins with the transition H3K4me3+ > H3K4me3+/H3K27ac+

#Now, save all the various transition data per transition type as .bed file, except E4-E4 transition
transition <- c("E1_E1", "E1_E2", "E1_E3", "E1_E4", "E1_E5",
                "E2_E1", "E2_E2", "E2_E3", "E2_E4", "E2_E5",
                "E3_E1", "E3_E2", "E3_E3", "E3_E4", "E3_E5",
                "E4_E1", "E4_E2", "E4_E3", "E4_E5",
                "E5_E1", "E5_E2", "E5_E3", "E5_E4", "E5_E5")


for (j in transition) {
  print(j)
  chrom_subset <- which(segment_combined$state == j)
  subseg <- segment_combined[chrom_subset,]
  write.table(x = subseg,
              file = paste0("D:/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_",j,"_regions.bed"), 
                            quote = F,
                            col.names = F,
                            row.names = F,
                            sep = "\t")
}
