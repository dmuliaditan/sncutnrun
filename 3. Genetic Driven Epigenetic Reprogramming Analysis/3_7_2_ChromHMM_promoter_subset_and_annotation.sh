#!/usr/bin/env bash

#3.7.2. This section describes the script used to subset binned data to include only promoter bins and to annotate the bins with the nearest genes in the HN137Pri > HN137Met transition
#Version 22/06/2022
#Daniel Muliaditan

#Required input:
#-Output .bed files for each transition as created in section 3.1.2.
#-Reference .bed files with promoter regions. This file was attained using the promoters() function in R.

#Required tools: bedtools

#Output: .bed files of each chromatin state transition, annotated with gene names and subset to promoter only regions.

for j in E1 E2 E3 E4 E5
do

  for k in E1 E2 E3 E4 E5
  
  do
  
  #Subset regions to include only promoter regions
  echo "Intersect regions that are in promoters:" "$j" ">" "$k"
  bedtools intersect -a '/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_'"$j"'_'"$k"'_regions.bed' \
  -b /mnt/d/genome_references/hg38_reference/promoter_regions_hg38_sorted.bed -sorted -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai >\
  '/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_'"$j"'_'"$k"'_promoters.bed'
  
  #Annotate with the closest genes in promoters
  bedtools closest -a '/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_'"$j"'_'"$k"'_promoters.bed' \
  -b /mnt/d/genome_references/hg38_reference/gencode.v32.basic.annotation.bed -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai -d \
  > '/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_'"$j"'_'"$k"'_regions_nearest_gene.bed'
  
  done
 
done





